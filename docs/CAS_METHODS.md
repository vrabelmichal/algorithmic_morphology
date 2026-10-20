# Implemented CAS pilot methods (shared v2, S profiles v1)

**Status:** both S profiles and both A outputs are implemented on fixed sky-subtracted
pixel arrays. This is *not* a bitwise clone of the original papers or
`statmorph`. The profile names deliberately end in `_like_v1`.
Synthetic tests cover core invariants and a 75-candidate DES pilot has run;
**DES-specific calibration and independent reference matching remain outstanding**.

## Papers and reference code

- Conselice (2003), ApJS 147, 1:
  https://arxiv.org/abs/astro-ph/0303065
- Lotz, Primack & Madau (2004), AJ 128, 163:
  https://arxiv.org/abs/astro-ph/0311352
- Rodriguez-Gomez et al. (2019), MNRAS 483, 4140:
  https://academic.oup.com/mnras/article/483/3/4140/5237717
- statmorph **0.7.2** reference implementation:
  https://statmorph.readthedocs.io/en/latest/_modules/statmorph/statmorph.html
- Conselice method description and alternative inner exclusion:
  https://academic.oup.com/mnras/article/418/2/801/1068863
- Background-bias caveat: Thorp et al. (2021),
  https://doi.org/10.1093/mnras/stab2201

## Shared conventions: `cas_discrete_v2`

The first real-data run uses shared method v2. Compared with v1, bilinear
rotation now requires all contributing pixels to be unmasked (nearest-neighbor
mask sampling could admit zero-filled contaminated neighbors), A_sky uses the
same symmetrically valid aperture flux denominator as A_raw, and center
optimization never accepts a result worse than the initial center. Aperture
sampling is restricted to a bounding box for speed without moving coordinates.
S profile parameters remain v1; centers can change under the corrected A method.
The exact source used for the run is archived with matching SHA-256 hashes in
`reports/des_pilot_20261020/run_source/`.


The input image must be **sky-subtracted**. Mask True = contaminating source
excluded, False = preserve target galaxy and sky. NaN/inf pixels are
excluded. Pixels are *not clipped positive* for photometry or the original
C growth curve. Pixels with nonpositive total galaxy flux or invalid
Petrosian radius are explicitly flagged.

1. **Radius**: compute discrete one-pixel circular annular sums about the
   input center using **pixel centers**, with eta(r) = mean[annulus] /
   mean[enclosed] and first descending crossing of eta=0.2, starting at
   three pixels. Linear interpolate the crossing. The aperture is
   R=1.5*Rpetro and must fit inside the image with sufficient margin for
   center optimization. A missing crossing **does not** trigger a fallback
   to a GALFIT radius.
2. **A center**: bounded Powell minimization (initial center ± up to two
   pixels, further limited to 0.25*Rpetro) of *raw A*. The optimized center
   is used for raw and corrected variants as well as C and S.
3. **C**: 5*log10(r80/r20), for 20% and 80% of the total signed flux
   within R. Distances are sorted by pixel-center radius. Piecewise-linear
   interpolation in cumulative flux is used. Strongly nonmonotonic or
   nonpositive growth curves are flagged. This approximation differs from
   statmorph's fractional-pixel `photutils` aperture root-finding.
4. **Rotation**: bilinear sample of sky-subtracted data at (2xc-x,2yc-y).
   All bilinear contributors to a rotated pixel must be unmasked, and the
   original pixel must be unmasked. Any out-of-image rotated position is invalid.
5. **A_raw**: sum |I-I180| / sum |I| over valid symmetrically unmasked
   pixels within R. **No additional factor of 2** in the denominator.
6. **Sky tile**: choose the cleanest (smallest unmasked standard deviation)
   of four corner square tiles. Each must have at least 95% valid pixels,
   lie entirely outside R+2 pixels, and be at least 8 pixels wide.
   Requested size defaults to 32 and may shrink to N/4 for small stamps.
   If none qualifies, the sky correction is **null**, not zero.
7. **A_sky**: average |B-B180| on valid tile pixels, multiplied by the
   galaxy's number of valid A pixels / sum|I| on the identical valid A pixels. **A_corrected=A_raw-A_sky**.
   These corrections are estimates and may be biased at low S/N; retain
   both and use noise/PSF experiments before interpreting them physically.

For reproducibility, the method version and exact profile parameters are
stored in the result rows and run manifest. Clumpiness numbers from different
profiles must *not* be compared to published thresholds without matching
all conventions.

## Parallel clumpiness profiles

For **both** profiles we compute the positive residual
max(I - boxcar(I), 0), integrated in a circular annulus from Rinner to
1.5Rpetro, divided by the signed galaxy flux in the **same annulus**.
Convolutions use `scipy.ndimage.uniform_filter` with `reflect` boundary
conditions and a **mask-normalized numerator and denominator**, so
contaminating pixels do not create spurious smoothing halos. Smooth
kernel size is `int(width_fraction * Rpetro)` pixels; sizes below two
pixels produce a null result and a flag.

| Profile output prefix | Boxcar width | Exclude center | Normalization |
| --- | --- | --- | --- |
| `S_lotz2004_statmorph_like_v1` | 0.25 Rpetro | inner 0.25 Rpetro | ×1 |
| `S_conselice2003_like_v1` | 0.30 Rpetro | inner 0.075 Rpetro | ×10 |

The original Conselice (2003) paper establishes 0.30 Rpetro smoothing and
×10 normalization, but does not settle every implementation choice in the
same way as modern code. The retained pilot's 0.075 value equals **0.05×the outer aperture radius
(1.5 Rpetro)**. This is a pilot adaptation, not the literal prescription in
Appendix A of Hambleton et al. (2011), which states **0.05×Rpetro** and describes
Rpetro/6 smoothing. The earlier wording incorrectly attributed our 0.075 value
to that appendix. Citation checking during report revision 2 corrected this
description; the user-approved profile parameters and measured values remain
unchanged. Other implementations also differ in inner exclusion. Any numerical
change would require a separately versioned profile and a matched comparison.

The statmorph 0.7.2 reference uses 0.25 for both scale and inner exclusion,
and ×1 normalization. Unlike our mask-normalized kernels and discrete
apertures, statmorph zero-fills masked pixels and uses fractional photometry.
The labels **like** are intentional.

For both profiles, store `_raw`, `_sky`, `_corrected` and diagnostic
kernel and inner-radius columns. The sky term uses the mean positive
high-frequency residual from the same selected blank tile (smoothed with
the same kernel), times the number of valid pixels in the galaxy annulus,
divided by galaxy annular flux, times the profile factor. When the
sky kernel is larger than the selected tile, record null corrections and
a quality flag. **Corrected S may be negative** and is retained.

## Data output

Each object-band has C, r20/r80/Rpetro, A_raw, A_sky, A_corrected,
both S_raw/sky/corrected, initial/optimized centers, mask fraction,
quality_flags, row status and input paths. A failed galaxy gets an explicit
error row. Each run emits incrementally committed Parquet parts, then a
canonical Parquet table, an HTML report showing all variants and their
direct comparisons (including count of valid measurements), and a
manifest. The report embeds the precise table SHA-256. The mutable
`latest` symlink changes only after complete publication.

## Important restrictions and next validation tasks

- A [real DES candidate pilot](../reports/des_pilot_20261020/README.md) has
  measured 162 available image-mask pairs for 75 IDs. All 225 requested rows
  remain in the table. The 129 `no_clean_sky_tile` flags are a major limitation:
  only two g, two r and 29 i measurements have sky corrections. No preferred
  A/S definition or physical classification follows from this selected sample.
- 144 deterministic synthetic noise/PSF/masked-contaminant cases were run;
  these stress tests are not a DES-specific uncertainty calibration.
- Check cutout FITS WCS, image units, band flips and mask alignment against
  the source project. `configs/example.py` uses placeholders and requires
  pointing at existing sky-subtracted images/masks.
- The corner-only sky estimator may miss valid sky regions or accidentally
  choose contaminated sky at low surface brightness. Compare alternative
  sky estimation methods and the Thorp et al. low-S/N correction.
- Quantify the impact of using pixel-center binning rather than
  fractional-pixel photometry, including C and Rpetro comparisons to a
  pinned statmorph version. It is possible one method is scientifically
  inappropriate for small, faint galaxies.
- Test size/PSF sensitivity and central nucleation with injected galaxies.
  Do not infer that either profile is better from agreement or histogram
  shape alone; compare failure rates, stability vs S/N, and known truth in
  simulations.
- Current runner does not resume a failed run and intentionally creates a
  **new** timestamped output directory on each invocation. Completed
  checkpoint partitions remain inspectable. It supports FITS/Parquet
  catalogs, but only FITS image cutouts in this first implementation.
- End-to-end pipeline tests require astropy and pyarrow; dependencies are
  declared in pyproject.toml and checked by GitHub Actions.


## Shared method v3: full-catalog optimizer robustness

The complete 84,200-ID run uses `cas_discrete_v3`. Its formulas, sky selection,
mask interpolation and S profile parameters match v2. SciPy Powell's rare
empty-line-search-direction ValueError now triggers the documented initial-center
fallback with `center_optimization_failed`, rather than losing the entire row.
The first full-catalog attempt exposed this issue and was stopped with its
checkpoints preserved. A regression test reproduces the exception. The new full
run starts in a fresh directory; no scientific rows are mixed across attempts.

See [extended diagnostics](EXTENDED_METHODS.md) for the added aperture Gini/M20,
shape outputs, exact upstream catalog identity and symlink-aware input inventory.
