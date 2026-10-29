import hashlib
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip('pyarrow')
fits = pytest.importorskip('astropy.io.fits')

from algorithmic_morphology.runner import run


def test_pipeline_checkpoint_report_and_nonoverwrite(tmp_path, monkeypatch):
    data=tmp_path/'images'; masks=tmp_path/'masks'
    data.mkdir(); masks.mkdir()
    y,x=np.indices((96,96))
    img=40*np.exp(-((x-47.5)**2+(y-47.5)**2)/(2*9**2))
    fits.writeto(data/'LSB_ID1.fits',img,overwrite=False)
    fits.writeto(masks/'mask1.fits',np.zeros((96,96),dtype=np.uint8),overwrite=False)
    fits.BinTableHDU.from_columns([fits.Column(name='COADD_OBJECT_ID',format='K',array=[1,2])]).writeto(tmp_path/'catalog.fits')
    c=tmp_path/'example.py'
    c.write_text(f'''from pathlib import Path
lsb_catalog_fits_path=Path({str(tmp_path/'catalog.fits')!r})
band_order=['g']
band_sky_subtracted_dir_paths={{'g':Path({str(data)!r})}}
band_mask_dir_paths={{'g':Path({str(masks)!r})}}
allow_geometric_center=True
workers=2
morphology_output_dir_path=Path({str(tmp_path/'output')!r})
''')
    first=run(c,chunk_size=1)
    assert len(list((first/'parts').glob('*.parquet')))==2
    table=next(first.glob('measurements__*.parquet'))
    report=next(first.glob('report__*.html'))
    import pandas as pd
    df=pd.read_parquet(table)
    assert list(df['status'])==['ok','error']
    text=report.read_text()
    assert hashlib.sha256(table.read_bytes()).hexdigest() in text
    assert 'lotz2004_statmorph_like_v1' in text
    assert 'conselice2003_like_v1' in text
    assert (tmp_path/'output'/'latest').resolve()==first.resolve()
    second=run(c,chunk_size=1)
    assert second!=first
    assert table.exists() and report.exists()
    assert (tmp_path/'output'/'latest').resolve()==second.resolve()
    # A report failure preserves the prior completed latest and durable parts.
    def fail_report(*args,**kwargs):
        raise RuntimeError('report failed')
    monkeypatch.setattr('algorithmic_morphology.runner.render_report',fail_report)
    with pytest.raises(RuntimeError,match='report failed'):
        run(c,chunk_size=1)
    assert (tmp_path/'output'/'latest').resolve()==second.resolve()
    import json
    incomplete=[p for p in (tmp_path/'output').iterdir() if p.is_dir() and p not in [first,second] and not p.is_symlink()]
    assert len(incomplete)==1
    assert json.loads((incomplete[0]/'progress.json').read_text())['status']=='incomplete'
    assert len(list((incomplete[0]/'parts').glob('*.parquet')))==2



def test_mask_wcs_mismatch_rejected(tmp_path):
    from astropy.wcs import WCS
    from algorithmic_morphology.runner import _one
    data = tmp_path/'data'; masks=tmp_path/'masks'
    data.mkdir(); masks.mkdir()
    w = WCS(naxis=2)
    w.wcs.crpix=[48,48]; w.wcs.cdelt=[-.000073,.000073]
    w.wcs.crval=[10,-2]; w.wcs.ctype=['RA---TAN','DEC--TAN']
    header=w.to_header()
    fits.writeto(data/'LSB_ID1.fits',np.ones((96,96)),header)
    header['CRPIX1'] += 2
    fits.writeto(masks/'mask1.fits',np.zeros((96,96)),header)
    settings={'band_sky_subtracted_dir_paths':{'g':data},'band_mask_dir_paths':{'g':masks}}
    with pytest.raises(ValueError,match='mask_wcs_mismatch'):
        _one(settings,{'coadd_object_id':'1','RA':10,'DEC':-2},'g')


def test_git_metadata_is_optional(monkeypatch):
    from algorithmic_morphology.runner import _git_provenance
    def unavailable(*args, **kwargs):
        raise FileNotFoundError('git unavailable')
    monkeypatch.setattr('algorithmic_morphology.runner.subprocess.check_output',unavailable)
    assert _git_provenance()=={'git_commit':None,'git_dirty':None}
