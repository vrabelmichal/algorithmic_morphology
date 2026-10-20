"""Example configuration: fill in local paths and use pre-sky-subtracted images."""
from pathlib import Path

lsb_catalog_fits_path = Path('/path/to/galaxies.fits')
band_order = ['g','r','i']
band_sky_subtracted_dir_paths = {
    'g':Path('/path/to/sky_subtracted/g'),
    'r':Path('/path/to/sky_subtracted/r'),
    'i':Path('/path/to/sky_subtracted/i'),
}
# Convention: 1 = contaminating source (masked), 0 = sky/target retained.
band_mask_dir_paths = {
    'g':Path('/path/to/combined_contaminant_masks/g'),
    'r':Path('/path/to/combined_contaminant_masks/r'),
    'i':Path('/path/to/combined_contaminant_masks/i'),
}
cutout_file_name_format = 'LSB_ID{coadd_object_id}.fits'
mask_file_name_format = 'mask{coadd_object_id}.fits'
require_mask = True
require_mask_wcs = False  # Set True for products expected to carry aligned WCS
allow_geometric_center = False  # Prefer RA/Dec + cutout WCS
skybox_size = 32
center_search = 2.0
morphology_output_dir_path = Path('.outputs/pilot')
