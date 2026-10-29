"""Match the complete upstream catalog; inventory products without selecting survivors."""
import argparse
import importlib
import json
from pathlib import Path
import sys
import pandas as pd
from astropy.table import Table
from algorithmic_morphology.report import sha256_file

DEFAULT='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_1024_wide_legacy_fliplr_object_centered__galfit_1024_sanitized__per_band_mag_reff__xy_2_5arcsec'
PRODUCTS='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_256__per_band_magnitude_optimization__r_constant_offset'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-project',type=Path,default=Path('../nucleated_lsbg'))
    parser.add_argument('--catalog-config',default=DEFAULT)
    parser.add_argument('--product-config',default=PRODUCTS)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=6)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    sys.path.insert(0,str(args.source_project.resolve()))
    catalog_config=importlib.import_module('configs.'+args.catalog_config)
    product_config=importlib.import_module('configs.'+args.product_config)
    catalog=Table.read(catalog_config.lsb_catalog_fits_path).to_pandas()
    if Path(catalog_config.lsb_catalog_fits_path).resolve()!=Path(product_config.lsb_catalog_fits_path).resolve():
        raise ValueError('Catalog configs differ; explicit reconciliation required')
    root=product_config.base_output_dir_path__for_masks
    bands=list(catalog_config.band_order)
    ids=catalog.COADD_OBJECT_ID.astype(str).tolist()
    if len(set(ids))!=len(ids):raise ValueError('Duplicate source IDs')
    rows=[];coverage=[]
    for b in bands:
        images={p.name.removeprefix('LSB_ID').removesuffix('.fits.gz') for p in (root/b/'final_sky_subtracted').iterdir() if p.name.endswith('.fits.gz')}
        masks={p.name.removeprefix('mask').removesuffix('.fits.gz') for p in (root/b/'final_combined_obj_removed_binary_mask').iterdir() if p.name.endswith('.fits.gz')}
        for ident in ids:
            rows.append(dict(coadd_object_id=ident,band=b,image_exists=ident in images,mask_exists=ident in masks))
        coverage.append(dict(band=b,requested=len(ids),images=len(set(ids)&images),masks=len(set(ids)&masks),paired=len(set(ids)&images&masks)))
    pd.DataFrame(rows).to_parquet(args.output/'inventory.parquet',index=False)
    pd.DataFrame(coverage).to_csv(args.output/'coverage.csv',index=False)
    metadata=dict(catalog=str(catalog_config.lsb_catalog_fits_path),catalog_sha256=sha256_file(catalog_config.lsb_catalog_fits_path),
                  objects=len(ids),requested_rows=len(rows),selection='Complete upstream catalog in original order; no availability or metric filter',
                  catalog_config=args.catalog_config,product_config=args.product_config,
                  product_root=str(root.resolve()),coverage=coverage,
                  source_configs={str(Path(c.__file__).resolve()):sha256_file(c.__file__) for c in [catalog_config,product_config]})
    (args.output/'selection.json').write_text(json.dumps(metadata,indent=2))
    config=f'''from pathlib import Path
lsb_catalog_fits_path=Path({str(Path(catalog_config.lsb_catalog_fits_path).resolve())!r})
band_order={bands!r}
_product_root=Path({str(root.resolve())!r})
band_sky_subtracted_dir_paths={{b:_product_root/b/'final_sky_subtracted' for b in band_order}}
band_mask_dir_paths={{b:_product_root/b/'final_combined_obj_removed_binary_mask' for b in band_order}}
cutout_file_name_format='LSB_ID{{coadd_object_id}}.fits.gz'
mask_file_name_format='mask{{coadd_object_id}}.fits.gz'
require_mask=True
require_mask_wcs=True
allow_geometric_center=False
skybox_size=32
center_search=2.0
workers={args.workers}
morphology_output_dir_path=Path({str(Path('.outputs/full_catalog').resolve())!r})
upstream_catalog_config={args.catalog_config!r}
upstream_product_config={args.product_config!r}
'''
    (args.output/'config.py').write_text(config)
    print(json.dumps(metadata,indent=2))


if __name__=='__main__':main()
