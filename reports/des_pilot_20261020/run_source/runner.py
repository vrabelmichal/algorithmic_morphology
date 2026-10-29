"""Config-based per-object CAS pipeline. Read-only inputs, checkpointed outputs."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import logging
import os
from pathlib import Path
import re
import subprocess
import importlib.metadata
import uuid

import numpy as np
import pandas as pd
from .metrics import measure_cas, PROFILES
from .report import render_report, sha256_file

LOG = logging.getLogger('algorithmic_morphology')


def _jsonable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(k):_jsonable(v) for k,v in value.items()}
    if isinstance(value, (list,tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value,(str,int,float,bool)) or value is None:
        return value
    return repr(value)


def _sync(path):
    with Path(path).open('rb') as stream:
        os.fsync(stream.fileno())
    fd = os.open(Path(path).parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _publish_json(path, obj):
    tmp = path.with_name('.'+path.name+'.'+uuid.uuid4().hex+'.tmp')
    with tmp.open('x',encoding='utf-8') as f:
        json.dump(obj,f,indent=2,sort_keys=True,default=str)
        f.flush();os.fsync(f.fileno())
    os.replace(tmp,path)
    _sync(path)


def _config(path):
    path = Path(path).expanduser().resolve()
    spec = importlib.util.spec_from_file_location('morphology_run_config',path)
    if spec is None or spec.loader is None:
        raise ValueError('invalid_config_module')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # As in nucleated-lsbg, configs are trusted Python.
    settings = {k:v for k,v in vars(module).items() if not k.startswith('_')
                and isinstance(v,(str,int,float,bool,list,tuple,dict,Path,type(None)))}
    return settings, path


def _catalog(settings):
    path = Path(settings['lsb_catalog_fits_path'])
    if path.suffix.lower() == '.parquet':
        df = pd.read_parquet(path)
    else:
        from astropy.table import Table
        df = Table.read(path).to_pandas()
    if 'coadd_object_id' not in df and 'COADD_OBJECT_ID' in df:
        df = df.rename(columns={'COADD_OBJECT_ID':'coadd_object_id'})
    if 'coadd_object_id' not in df:
        raise ValueError('catalog_missing_coadd_object_id')
    df['coadd_object_id'] = df['coadd_object_id'].astype(str)
    if df['coadd_object_id'].duplicated().any():
        raise ValueError('catalog_duplicate_object_ids')
    return df


def _read_fits(path):
    from astropy.io import fits
    with fits.open(path,memmap=False) as hdus:
        hdu = next((h for h in hdus if h.data is not None and np.ndim(h.data)==2),None)
        if hdu is None:
            raise ValueError('fits_missing_2d_image')
        return np.asarray(hdu.data,dtype=float),hdu.header.copy()


def _resolve_file(directory, ident, pattern):
    if directory is None:
        return None
    # No catalog string may inject path traversal.
    ident = str(ident)
    if not re.fullmatch(r'[A-Za-z0-9._-]+',ident) or ident.startswith('.'):
        raise ValueError('unsafe_object_id')
    filename = pattern.format(coadd_object_id=ident)
    if Path(filename).name != filename:
        raise ValueError('unsafe_file_pattern')
    return Path(directory)/filename


def _center(settings, row, header, shape):
    # In WCS images use object RA/Dec to avoid assuming centered stamps.
    ra = row.get('ra',row.get('RA'))
    dec = row.get('dec',row.get('DEC'))
    if ra is not None and dec is not None and np.isfinite(float(ra)) and np.isfinite(float(dec)):
        try:
            from astropy.wcs import WCS
            x,y = WCS(header).celestial.world_to_pixel_values(float(ra),float(dec))
            if np.isfinite(x) and np.isfinite(y):
                return float(x),float(y)
        except (ValueError,TypeError,KeyError,AttributeError):
            pass
    if settings.get('allow_geometric_center',False):
        return (shape[1]-1)/2,(shape[0]-1)/2
    raise ValueError('center_missing_wcs_or_opt_in')


def _one(settings, row, band):
    ident = row['coadd_object_id']
    img_dirs = settings.get('band_sky_subtracted_dir_paths')
    if img_dirs is None:
        if not settings.get('input_is_sky_subtracted',False):
            raise ValueError('input_not_declared_sky_subtracted')
        img_dirs = settings.get('band_cutouts_dir_paths') or {band:settings.get('cutouts_dir_path')}
    primary = img_dirs.get(band)
    fallback = settings.get('fallback_band_cutouts_dir_paths',{}).get(band)
    fmt = settings.get('cutout_file_name_format','LSB_ID{coadd_object_id}.fits')
    selected = _resolve_file(primary,ident,fmt)
    if selected is None or not selected.exists():
        selected = _resolve_file(fallback,ident,fmt)
    if selected is None or not selected.exists():
        raise FileNotFoundError('missing_cutout')
    image,header = _read_fits(selected)
    masks = settings.get('band_mask_dir_paths',{})
    mask_dir = masks.get(band)
    mask_file = _resolve_file(mask_dir,ident,settings.get('mask_file_name_format','mask{coadd_object_id}.fits'))
    if mask_file is None or not mask_file.exists():
        if settings.get('require_mask',True):
            raise FileNotFoundError('missing_contamination_mask')
        contamination = None
    else:
        mask_im,mask_header = _read_fits(mask_file)
        contamination = (mask_im != 0) | ~np.isfinite(mask_im)
        if mask_im.shape != image.shape:
            raise ValueError('mask_shape_mismatch')
        from astropy.wcs import WCS
        iw, mw = WCS(header).celestial, WCS(mask_header).celestial
        if iw.has_celestial and mw.has_celestial:
            y, x = np.indices((3, 3), dtype=float)
            x *= (image.shape[1]-1)/2
            y *= (image.shape[0]-1)/2
            mx, my = mw.world_to_pixel_values(*iw.pixel_to_world_values(x, y))
            if not np.allclose(mx, x, atol=1e-3) or not np.allclose(my, y, atol=1e-3):
                raise ValueError('mask_wcs_mismatch')
        elif settings.get('require_mask_wcs', False):
            raise ValueError('mask_wcs_missing')
    center = _center(settings,row,header,image.shape)
    measurements = measure_cas(image,mask=contamination,center=center,
                               skybox_size=int(settings.get('skybox_size',32)),
                               center_search=float(settings.get('center_search',2.0)))
    measurements.update({'image_file':str(selected.resolve()),'mask_file':str(mask_file.resolve()) if mask_file and mask_file.exists() else None,
                         'image_sha256':sha256_file(selected),
                         'mask_sha256':sha256_file(mask_file) if mask_file and mask_file.exists() else None,
                         'image_unit':header.get('BUNIT'),
                         'image_height':image.shape[0], 'image_width':image.shape[1],
                         'center_source':'catalog_wcs' if not settings.get('allow_geometric_center', False) else 'wcs_or_geometric_opt_in'})
    return measurements


def _checkpoint(run_dir, rows, n, state):
    if not rows:
        return
    part = run_dir/'parts'/f'part-{n:06d}.parquet'
    tmp = part.with_name('.'+part.name+'.'+uuid.uuid4().hex+'.tmp')
    pd.DataFrame(rows).to_parquet(tmp,engine='pyarrow',index=False)
    os.replace(tmp,part)
    _sync(part)
    state['parts'].append({'file':str(part.relative_to(run_dir)),'rows':len(rows),'sha256':sha256_file(part)})
    state['checkpointed'] += len(rows)
    _publish_json(run_dir/'progress.json',state)
    LOG.info('Checkpoint %s; %d rows durable; failures=%d',part.name,state['checkpointed'],state['failed'])


def _latest(root,run_id):
    temp = root/('.latest.'+uuid.uuid4().hex)
    os.symlink(run_id,temp)
    os.replace(temp,root/'latest')


def run(config_file, *, limit=None, chunk_size=50):
    if chunk_size<1:
        raise ValueError('chunk_size_must_be_positive')
    cfg,source = _config(config_file)
    catalog = _catalog(cfg)
    if limit is not None:
        if limit<0: raise ValueError('limit_must_be_nonnegative')
        catalog = catalog.head(limit)
    bands = cfg.get('band_order') or ['g']
    if len(set(bands))!=len(bands): raise ValueError('duplicate_band_names')
    basedir = Path(cfg.get('morphology_output_dir_path',Path('.outputs')/source.stem))
    basedir.mkdir(parents=True,exist_ok=True)
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')+'_'+uuid.uuid4().hex[:8]
    run_dir = basedir/run_id
    run_dir.mkdir(exist_ok=False)
    (run_dir/'parts').mkdir(exist_ok=False)
    cfg_bytes = source.read_bytes()
    chash = hashlib.sha256(cfg_bytes).hexdigest()
    effective={k:_jsonable(v) for k,v in cfg.items()}
    _publish_json(run_dir/'effective_config.json',{'settings':effective,'source':str(source),'sha256':chash})
    code_files = list(Path(__file__).parent.glob('*.py'))
    provenance = {
        'git_commit': subprocess.check_output(['git','rev-parse','HEAD'], text=True).strip(),
        'git_dirty': bool(subprocess.check_output(['git','status','--porcelain'], text=True).strip()),
        'code_sha256': {p.name:sha256_file(p) for p in code_files},
        'versions': {name:importlib.metadata.version(name) for name in
                     ['numpy','scipy','pandas','pyarrow','astropy','matplotlib']},
        'catalog_sha256':sha256_file(Path(cfg['lsb_catalog_fits_path'])),
        'method_version':'cas_discrete_v2',
        'parameters':{'skybox_size':cfg.get('skybox_size',32),
                      'center_search':cfg.get('center_search',2.0),
                      'petrosian_eta':0.2,'aperture_factor':1.5},
    }
    state={'run_id':run_id,'status':'running','count':len(catalog)*len(bands),
           'checkpointed':0,'failed':0,'invalid':0,'provenance':provenance,'parts':[],'methods':[p.__dict__ for p in PROFILES],
           'config_sha256':chash,'catalog':str(cfg['lsb_catalog_fits_path'])}
    _publish_json(run_dir/'progress.json',state)
    LOG.info('Starting run=%s; galaxies=%d; bands=%s',run_id,len(catalog),bands)
    rows=[]; part_no=0
    try:
        for record in catalog.to_dict('records'):
            for band in bands:
                item={'run_id':run_id,'coadd_object_id':str(record['coadd_object_id']),'band':band}
                try:
                    item.update(_one(cfg,record,band))
                    item['status']='ok' if all(item.get(k) is not None for k in ['A_raw','C']) else 'invalid'
                    item['error_type']=None
                    state['invalid'] += int(item['status']=='invalid')
                except Exception as exc:  # one bad galaxy/band must not abort catalog
                    item.update({'status':'error','error_type':type(exc).__name__,'error_message':str(exc)})
                    state['failed'] += 1
                    LOG.warning('object=%s band=%s: %s: %s',item['coadd_object_id'],band,type(exc).__name__,exc)
                rows.append(item)
                if len(rows)>=chunk_size:
                    part_no += 1
                    _checkpoint(run_dir,rows,part_no,state)
                    rows=[]
        if rows:
            part_no += 1
            _checkpoint(run_dir,rows,part_no,state)
        table = run_dir/f'measurements__{run_id}.parquet'
        tmp = table.with_name('.'+table.name+'.'+uuid.uuid4().hex+'.tmp')
        # Stream partitions through PyArrow; avoid collecting all galaxies in RAM.
        import pyarrow as pa
        import pyarrow.parquet as pq
        paths = [run_dir/entry['file'] for entry in state['parts']]
        if not paths:
            raise ValueError('empty_catalog_or_no_partitions')
        schema = pa.unify_schemas([pq.read_schema(path) for path in paths],
                                  promote_options='permissive')
        writer = pq.ParquetWriter(tmp, schema)
        seen = set()
        recorded = 0
        try:
            for path in paths:
                tab = pq.read_table(path)
                keys = list(zip(tab['coadd_object_id'].to_pylist(),tab['band'].to_pylist()))
                if len(keys) != len(set(keys)) or any(k in seen for k in keys):
                    raise RuntimeError('duplicate_object_band_keys')
                seen.update(keys)
                recorded += len(keys)
                aligned = pa.Table.from_arrays(
                    [tab[field.name].cast(field.type) if field.name in tab.column_names
                     else pa.nulls(tab.num_rows,type=field.type)
                     for field in schema],schema=schema)
                writer.write_table(aligned)
        finally:
            writer.close()
        if recorded != state['count']:
            raise RuntimeError('partitions_failed_integrity_check')
        os.replace(tmp,table)
        _sync(table)
        report=run_dir/f'report__{run_id}.html'
        html=render_report(table,run_id,chash,inputs={'catalog':str(cfg['lsb_catalog_fits_path']),'provenance':provenance})
        report_temp=report.with_name('.'+report.name+'.'+uuid.uuid4().hex+'.tmp')
        with report_temp.open('x',encoding='utf-8') as f:
            f.write(html); f.flush(); os.fsync(f.fileno())
        os.replace(report_temp,report)
        _sync(report)
        state.update({'status':'complete_with_errors' if state['failed'] else 'complete',
                      'table':table.name,'table_sha256':sha256_file(table),'report':report.name})
        _publish_json(run_dir/'manifest.json',state)
        _publish_json(run_dir/'progress.json',state)
        _latest(basedir,run_id)
        LOG.info('Run complete: %s; %s; %s',run_id,table,report)
        return run_dir
    except BaseException:
        state['status']='incomplete'
        _publish_json(run_dir/'progress.json',state)
        LOG.exception('Run incomplete: %s',run_id)
        raise
