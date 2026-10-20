# CAS morphology: literature review (2026-10-20)

**Status: literature background and open scientific questions.** The first pilot now implements both raw and corrected asymmetry, plus two explicit, versioned clumpiness profiles (see [CAS_METHODS.md](CAS_METHODS.md)). They are *-like* adaptations requiring real-data and matched-code validation, not exact historical reimplementations.

CAS measures non-parametric **concentration (C)**, **180-degree asymmetry (A)**, and **clumpiness / smoothness (S)**. A fixed algorithm is repeatable on fixed input, but the results depend on image noise, PSF, band, sky, masks, radii, centering and processing details. These choices matter especially for faint, diffuse LSBGs.

## Literature

1. Bershady, Jangren & Conselice (2000), *Structural and Photometric Classification of Galaxies I* — concentration and morphological calibration: https://arxiv.org/abs/astro-ph/0002262
2. Conselice, Bershady & Jangren (2000), *The Asymmetry of Galaxies*, ApJ 529, 886 — center minimization and noise: https://ned.ipac.caltech.edu/NED::Refcode/2000ApJ...529..886C
3. Conselice (2003), *The Relationship between Stellar Light Distributions of Galaxies and Their Formation Histories*, ApJS 147, 1 — original CAS framework: https://arxiv.org/abs/astro-ph/0303065
4. Lotz, Primack & Madau (2004), *A New Non-Parametric Approach to Galaxy Morphological Classification*, AJ 128, 163 — widely used implementation conventions: https://arxiv.org/abs/astro-ph/0311352
5. Rodriguez-Gomez et al. (2019), *The optical morphologies of galaxies in the IllustrisTNG simulation*, MNRAS 483, 4140 — **statmorph** library: https://academic.oup.com/mnras/article/483/3/4140/5237717 . Its documented reference source is https://statmorph.readthedocs.io/en/latest/_modules/statmorph/statmorph.html ; authors explicitly do **not** call this an official implementation.
6. Conselice's review (2014), CAS section: https://ned.ipac.caltech.edu/level5/March14/Conselice/Conselice2.html
7. *Advanced morphological galaxy classification* (2011), Appendix A, alternative clumpiness implementations: https://academic.oup.com/mnras/article/418/2/801/1068863
8. *The Lyman alpha reference sample IV* (2015), central smoothness exclusion: https://doi.org/10.1051/0004-6361/201425053
9. *Towards robust determination of non-parametric morphologies in marginal astronomical data* (2021), low-S/N bias: https://academic.oup.com/mnras/article/507/1/886/6332310

## Equations and disagreements

### Concentration

The familiar definition is

    C = 5 log10(r80 / r20)

Here r20 and r80 enclose 20% and 80% of a *specified total galaxy flux*. A widely used CAS choice defines total flux inside a **circular** aperture of 1.5 times the Petrosian radius R_p, with R_p where local-to-enclosed mean brightness equals 0.2. Some implementations instead use 2 R_p, an elliptical aperture, or a different fractional pair. Explicitly specify how the Petrosian annulus is measured (some use 0.8–1.25 R), whether masked pixels contribute to aperture area, and what happens when a faint, noisy growth curve is nonmonotonic. **Do not silently substitute GALFIT r_eff for R_p.** If total flux is nonpositive, or r20/r80 cannot be determined, flag invalid instead of forcing a plausible number.

### Asymmetry

One common convention is

    A_raw = sum |I - rotate180(I, center)| / sum |I|

where sums are within a declared measurement aperture, commonly 1.5 R_p. Find the center by minimizing A in a constrained subpixel search, rather than assuming the cutout center, catalogue center or GALFIT center. The sky contribution must be estimated from a genuinely blank region and normalized appropriately; report both raw and corrected values. The *rotated contaminant mask* must be combined with the original mask to avoid artificial differences. Rotation interpolation and edge behavior affect A.

**Important inconsistency:** published definitions also use **2 sum|I|** in the denominator; such values are half the above raw values. Their absolute thresholds are not interchangeable. At low S/N, classical sky subtraction can itself be biased (source 9). Keep raw A, sky correction, corrected A, optimized x/y and center/sky quality flags.

### Clumpiness / smoothness

Often defined from *positive* high-frequency residuals

    S_base = [sum max(I - smooth(I), 0) - scaled sky residual] / galaxy flux

over a declared aperture, typically with an unresolved central region excluded. Precisely which region is excluded from the numerator and denominator must be specified.

**There is no one accepted S setting.** Conselice's original form multiplies by **10**; Lotz/statmorph do **not**. Conselice (2003) and the later review describe boxcar smoothing of about **0.3 R_p**; statmorph 0.7.2 defaults to **0.25 R_p** and excludes the inner **0.25 R_p**. Another paper (source 7) discusses a central **0.05 R_p** exclusion with smoothing **R_p/6**; source 8 describes **0.3 R_p** inner exclusion. These produce materially different values. Published numerical S values cannot be compared without matching all rules.

In particular, zero-filling masked stars *before* boxcar convolution may produce spurious residuals at the mask edge. We must agree whether to use mask-normalized convolution or a documented inpainting procedure; validate both if useful. A negative corrected S may be a valid noise-corrected result, not necessarily a broken measurement.

## Decisions to settle before algorithm implementation

| ID | Question | Proposed baseline for discussion | Status |
| --- | --- | --- | --- |
| D1 | Which named reference? | A **Lotz (2004) / statmorph-compatible** variant, identified by exact code/settings; optionally a distinct Conselice (2003) variant later | OPEN |
| D2 | R_p and outer aperture? | eta=0.2, measured circular R_p, outer radius 1.5 R_p; explicitly define annulus, interpolation and edge failures | OPEN |
| D3 | Center? | Minimize A around the target centroid with bounded subpixel search, retaining optimized x/y | OPEN |
| D4 | A denominator/sky? | Unscaled ratio, store A_raw, A_sky, A_corrected, test optional low-S/N correction separately | OPEN |
| D5 | S definition? | Select exactly one smoothing scale, filter, kernel edges, central exclusion, sky correction, and x1/x10 normalization **after a comparison** | OPEN |
| D6 | Aperture geometry? | Circular CAS main values; future elliptical variants under distinct names (q and PA must be explicit) | OPEN |
| D7 | Mask semantics? | Target-preserving contaminant masks, symmetrized for A, controlled smoothing near masks | OPEN |
| D8 | Invalid and negative pixels? | Preserve signed sky-subtracted values; null+flag unmeasurable values; never treat -99 as science | OPEN |
| D9 | Bands, PSF, units? | Independent g/r/i results with input PSF, pixel scale and sky provenance | OPEN |
| D10 | Validation? | Matched statmorph method/options + analytic synthetic images + faint DES galaxy sample | OPEN |

**Priority discussion:** D1, D2, D4, D5 and D7. Do not implement measurements until these have explicit choices.

## Validation plan (not yet run)

Test noise-free symmetric profiles and offset-clump variants, noise/sky gradients, bad background, no sky pixels, central bright point sources, elongated profiles, mask holes, cutout edges, missing FITS, NaNs, invalid total flux and mask orientation. Compare with statmorph only when the **entire** definition is matched, including pre-processing. For a small real DES LSBG pilot, plot A/C/S across bands versus S/N, R_p/PSF, masking fraction and quality flags; quantify failures rather than discarding them. A reference disagreement can indicate a mistaken interpretation in *either* implementation, not necessarily a regression.


## Follow-up scientific decision, 2026-10-20

The user requested **both A_raw and A_corrected**, and **both S profiles** to be measured on the same galaxies, preserving exact method names in Parquet and comparing distributions/scatter in the report before choosing any default. The question is now whether the v1 discrete/approximate methods recover the intended literature-based behavior on real cutouts, not which profile to silently choose. Detailed implemented values (0.25/0.25/×1 and 0.30/0.075/×10) and their limitations are recorded in [CAS_METHODS.md](CAS_METHODS.md).
