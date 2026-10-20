# Additional aperture diagnostics, v1

The full-catalog run retains the existing CAS method and both S profiles, and
adds separately named aperture diagnostics. They use the same unmasked finite
pixels inside the measured circular 1.5 Rpetro aperture at the optimized A
center. These are not statmorph's segmentation-based Gini/M20 values. Published
classification boundaries should not be applied to these aperture adaptations.

The input sky-subtracted image stays signed. No alteration of CAS photometry
occurs. A missing Petrosian aperture gives null diagnostics. At least ten valid
pixels and positive absolute flux are required. Output names explicitly encode
`aperture_v1`; `aperture_diagnostics_version` records the shared definition.

- **Gini_aperture_v1**: sort absolute pixel fluxes X_i in ascending order and
  calculate sum[(2i−n−1)X_i] / [(n−1) sum X_i], with i starting at one.
- **M20_aperture_v1**: calculate the centroid using signed flux, total signed
  second central moment, and the corresponding moment of the brightest pixels
  through the first 20% signed-total-flux crossing. Include all pixels tied at
  the crossing brightness. Return log10(bright moment / total moment), or null
  with a flag for nonpositive total flux or moments. Negative noise can make
  signed moment estimates unstable; this is exposed, not silently clipped.
- **axis_ratio_positive_aperture_v1**: calculate the flux-weighted covariance
  matrix using max(I,0) weights; q=sqrt(smaller eigenvalue / larger eigenvalue).
- **ellipticity_positive_aperture_v1**: 1−q using those positive weights.
- **position_angle_positive_aperture_v1**: major covariance eigenvector angle
  in [0,180) degrees, measured from increasing image x toward increasing image
  y. This is an image-coordinate angle, not a celestial position angle. Its
  physical interpretation becomes unstable for nearly round objects.

Positive weights are confined to the explicitly named shape diagnostics. They
can bias faint/noisy shapes, and no PSF deconvolution is applied. The Gini and
M20 formulas are motivated by Lotz et al. (2004), but use our declared aperture
rather than their source segmentation. The reference code also illustrates
absolute flux for Gini and signed moments for M20:

- https://arxiv.org/abs/astro-ph/0311352
- https://statmorph.readthedocs.io/en/latest/_modules/statmorph/statmorph.html

Tests cover uniform and single-pixel Gini limits, analytic resolved Gaussian
M20, known Gaussian axis ratio, reflection invariance, masked contaminants and
nonpositive-flux handling. They validate these definitions, not a physical
classification or exact historical/segmentation-code equivalence.

# Full upstream selection

The recent wide/sanitized configurations and the latest 256-pixel
`per_band_magnitude_optimization__r_constant_offset` configuration all point to
`concatenated-detr-vit-all.fits` (84,200 IDs). The full run uses that exact catalog
in its original order and r/g/i band order, without an availability selection.
It reads the existing two-pass sky-subtracted 256-pixel products and masks shared
by the upstream 256 configurations. It does not use GALFIT models or radii, run
new GALFIT fits, or substitute raw sky-containing cutouts for these products.
Wide-image experiments have much smaller product coverage and remain distinct.

The input tree is a symlink into `/media/vrabel/DATA_KIN4TB`; directory discovery
must resolve that link explicitly. Availability is 84,188 r, 84,172 g and 84,144 i
image-mask pairs; all 252,600 requested rows must remain in the final table,
including the 96 missing combinations and any measurement/validation failures.
