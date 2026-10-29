from pathlib import Path
lsb_catalog_fits_path = Path('/mnt/data_sm870qvo/lsbg/analysis/algorithmic_morphology/.outputs/des_inputs_v2/catalog.parquet')
band_order = ['g', 'r', 'i']
band_sky_subtracted_dir_paths = {b: Path('/mnt/data_sm870qvo/lsbg/analysis/nucleated_lsbg/.outputs/sia/sky_subtracted_multiband_1024_wide_legacy_fliplr_object_centered/reff_mul_2_0/circ_upper_0_8/reff_mul_for_source_1_5/circ_upper_for_source_0_8/init_sep_2_0/final_sep_high_2_0/final_sep_low_2_0/initial_dilation_9/outside_dilation_9/object_non_masked_radius_2_0')/b/'final_sky_subtracted' for b in band_order}
band_mask_dir_paths = {b: Path('/mnt/data_sm870qvo/lsbg/analysis/nucleated_lsbg/.outputs/sia/sky_subtracted_multiband_1024_wide_legacy_fliplr_object_centered/reff_mul_2_0/circ_upper_0_8/reff_mul_for_source_1_5/circ_upper_for_source_0_8/init_sep_2_0/final_sep_high_2_0/final_sep_low_2_0/initial_dilation_9/outside_dilation_9/object_non_masked_radius_2_0')/b/'final_combined_obj_removed_binary_mask' for b in band_order}
cutout_file_name_format = 'LSB_ID{coadd_object_id}.fits.gz'
mask_file_name_format = 'mask{coadd_object_id}.fits.gz'
require_mask = True
require_mask_wcs = True
allow_geometric_center = False
skybox_size = 32
center_search = 2.0
morphology_output_dir_path = Path('/mnt/data_sm870qvo/lsbg/analysis/algorithmic_morphology/.outputs/des_pilot')
