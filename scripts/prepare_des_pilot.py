"""Inventory existing DES products and write an availability-selected pilot config.

No source products are modified. Selection uses file availability, never CAS values.
"""
import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd
from astropy.table import Table
from astropy.wcs import WCS
from algorithmic_morphology.runner import _read_fits
from algorithmic_morphology.report import sha256_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-project', type=Path, default=Path('../nucleated_lsbg'))
    parser.add_argument('--catalog', type=Path, default=Path('../clustering_lsbg/.outputs/rejected_lsb_datasets/LSBG_DETR_and_ViT_candidates/v3/concatenated-detr-vit-all.fits'))
    parser.add_argument('--output', type=Path, default=Path('.outputs/des_inputs'))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    roots = list((args.source_project/'.outputs/sia/sky_subtracted_multiband_1024_wide_legacy_fliplr_object_centered').glob('**/r/final_sky_subtracted'))
    if len(roots) != 1:
        raise ValueError(f'Expected one product root, found {len(roots)}')
    root = roots[0].parent.parent.resolve()
    bands = ['g','r','i']
    ids = sorted({p.name.removeprefix('LSB_ID').removesuffix('.fits.gz')
                  for band in bands for p in (root/band/'final_sky_subtracted').glob('LSB_ID*.fits.gz')}, key=int)
    catalog = Table.read(args.catalog).to_pandas()
    sample = catalog[catalog.COADD_OBJECT_ID.astype(str).isin(ids)].copy()
    if len(sample) != len(ids):
        raise ValueError('Product IDs absent from catalog')
    sample = sample.sort_values('COADD_OBJECT_ID')
    selected = args.output/'catalog.parquet'
    sample.to_parquet(selected, index=False)
    rows = []
    for row in sample.to_dict('records'):
        ident = str(row['COADD_OBJECT_ID'])
        for band in bands:
            image = root/band/'final_sky_subtracted'/f'LSB_ID{ident}.fits.gz'
            mask = root/band/'final_combined_obj_removed_binary_mask'/f'mask{ident}.fits.gz'
            item = dict(coadd_object_id=ident, band=band, image=str(image), mask=str(mask),
                        image_exists=image.exists(), mask_exists=mask.exists())
            if image.exists() and mask.exists():
                data, header = _read_fits(image)
                contamination, mh = _read_fits(mask)
                x,y = WCS(header).celestial.world_to_pixel_values(row['RA'], row['DEC'])
                mx,my = WCS(mh).celestial.world_to_pixel_values(row['RA'], row['DEC'])
                item.update(shape=str(data.shape), center_x=float(x), center_y=float(y),
                            wcs_mask_offset=float(np.hypot(x-mx,y-my)),
                            mask_fraction=float((contamination!=0).mean()),
                            mask_binary=bool(np.isin(contamination,[0,1]).all()),
                            image_sha256=sha256_file(image),mask_sha256=sha256_file(mask))
                if not (0 <= x < data.shape[1] and 0 <= y < data.shape[0]):
                    raise ValueError(f'Target outside cutout: {ident} {band}')
                item['target_masked'] = bool(contamination[int(round(float(y))), int(round(float(x)))])
            rows.append(item)
    pd.DataFrame(rows).to_csv(args.output/'inventory.csv',index=False)
    meta = dict(source_catalog=str(args.catalog.resolve()),source_catalog_sha256=sha256_file(args.catalog),
                selected_catalog_sha256=sha256_file(selected),objects=len(sample),rows=len(rows),
                selection='Union of available precomputed sky-subtracted g/r/i product IDs; sorted by ID; no morphology selection',
                product_root=str(root),orientation='Upstream legacy fliplr object-centered products; no additional flip',
                upstream_git=__import__('subprocess').check_output(['git','-C',str(args.source_project),'rev-parse','HEAD'],text=True).strip())
    (args.output/'selection.json').write_text(json.dumps(meta,indent=2))
    cfg = f'''from pathlib import Path
lsb_catalog_fits_path = Path({str(selected.resolve())!r})
band_order = {bands!r}
band_sky_subtracted_dir_paths = {{b: Path({str(root)!r})/b/'final_sky_subtracted' for b in band_order}}
band_mask_dir_paths = {{b: Path({str(root)!r})/b/'final_combined_obj_removed_binary_mask' for b in band_order}}
cutout_file_name_format = 'LSB_ID{{coadd_object_id}}.fits.gz'
mask_file_name_format = 'mask{{coadd_object_id}}.fits.gz'
require_mask = True
require_mask_wcs = True
allow_geometric_center = False
skybox_size = 32
center_search = 2.0
morphology_output_dir_path = Path({str(Path('.outputs/des_pilot').resolve())!r})
'''
    (args.output/'config.py').write_text(cfg)
    print(json.dumps(meta,indent=2))


if __name__ == '__main__':
    main()
