from pathlib import Path
lsb_catalog_fits_path=Path('/mnt/data_sm870qvo/lsbg/analysis/clustering_lsbg/.outputs/rejected_lsb_datasets/LSBG_DETR_and_ViT_candidates/v3/concatenated-detr-vit-all.fits')
band_order=['r', 'g', 'i']
_product_root=Path('/media/vrabel/DATA_KIN4TB/lsbg/analysis/nucleated_lsbg/.outputs/sia/sky_subtracted_multiband_256/reff_mul_2_0/circ_upper_0_8/reff_mul_for_source_1_5/circ_upper_for_source_0_8/init_sep_2_0/final_sep_high_2_0/final_sep_low_2_0/initial_dilation_9/outside_dilation_9/object_non_masked_radius_2_0')
band_sky_subtracted_dir_paths={b:_product_root/b/'final_sky_subtracted' for b in band_order}
band_mask_dir_paths={b:_product_root/b/'final_combined_obj_removed_binary_mask' for b in band_order}
cutout_file_name_format='LSB_ID{coadd_object_id}.fits.gz'
mask_file_name_format='mask{coadd_object_id}.fits.gz'
require_mask=True
require_mask_wcs=True
allow_geometric_center=False
skybox_size=32
center_search=2.0
workers=10
morphology_output_dir_path=Path('/mnt/data_sm870qvo/lsbg/analysis/algorithmic_morphology/.outputs/full_catalog')
upstream_catalog_config='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_1024_wide_legacy_fliplr_object_centered__galfit_1024_sanitized__per_band_mag_reff__xy_2_5arcsec'
upstream_product_config='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_256__per_band_magnitude_optimization__r_constant_offset'
