# Algorithmic Morphology

Explainable, literature-motivated measurements of galaxy-image morphology.
The first implementation measures **concentration (C)**, **raw and
sky-corrected 180° asymmetry (A)**, and **two named clumpiness/smoothness
(S) profiles** for each galaxy and photometric band.

**Status: pilot implementation exercised on 75 locally available DES candidates; scientific calibration remains open.**
The profiles are documented adaptations, not exact numerical reproductions
of every implementation in the literature.

## Metrics produced

- `C = 5 log10(r80/r20)` in a Petrosian aperture.
- `A_raw`, `A_sky` and `A_corrected`, using the *same* optimized center.
- `S_lotz2004_statmorph_like_v1_{raw,sky,corrected}`: boxcar width
  0.25 Rpetro, central exclusion 0.25 Rpetro, ×1 normalization.
- `S_conselice2003_like_v1_{raw,sky,corrected}`: boxcar width
  0.30 Rpetro, central exclusion 0.075 Rpetro, ×10 normalization.

The raw and corrected A values and the two S definitions are
**saved together and compared** in the HTML report. No winner is selected.

See [implemented methods and their limitations](docs/CAS_METHODS.md)
and the [scientific literature review](docs/CAS_LITERATURE_REVIEW.md).

## Full upstream catalog processing

The full run selects the same **84,200 IDs** as the recent `nucleated_lsbg`
configurations, with no availability filter: **252,600 object-band rows**.
Existing 256-pixel sky-subtracted products provide 252,504 image-mask pairs.
The earlier 75-object report missed this symlinked product tree and is a pilot,
not the full catalog result. [Processing status](reports/full_catalog_20261020_progress.html)
is explicitly provisional until the full analysis is published.

The run also evaluates [aperture Gini/M20 and shape diagnostics](docs/EXTENDED_METHODS.md),
with metric-specific validity counts. These are separately versioned aperture
adaptations rather than statmorph segmentation-based outputs.

```bash
python scripts/prepare_full_catalog.py --output .outputs/full_inputs_new
python scripts/run_full_catalog.py --config .outputs/full_inputs_new/config.py --inputs .outputs/full_inputs_new --output reports/full_catalog_new
```

The second command processes every upstream ID, then automatically publishes
its self-contained HTML findings, plots, per-object metrics and provenance.

## Installation and pilot execution

Python 3.10+:

```bash
python -m pip install -e '.[test]'
python -m pytest -q
cp configs/example.py configs/my_run.py
# Edit the paths to existing sky-subtracted FITS cutouts, masks and catalog
python -m algorithmic_morphology run --config-file-name configs/my_run.py --limit 50
```

The configuration accepts `lsb_catalog_fits_path`, `band_order`, and
related Python/Path conventions from
[nucleated-lsbg](https://github.com/vrabelmichal/nucleated-lsbg).
It requires sky-subtracted input; original raw cutouts must not be
mistaken for already-background-subtracted images. Mask polarity is:
**nonzero = contaminating object excluded**; the sky and target galaxy
must remain unmasked. Unless explicitly opted out of, masks are required.
WCS + catalog coordinates are used to find target centers by default.

The command writes **Parquet checkpoints periodically** without
overwriting prior runs. Each unique UTC-timestamped directory contains
immutable output data, provenance, a SHA-256-bound HTML report with
distributions and direct comparisons, progress and final manifests.
`latest` is a symlink to the most recent *completed* run only.
Log messages with UTC timestamps go to stdout. A missing galaxy/band is
recorded as an error row, not silently skipped.

The [2026-10-20 analysis report](reports/des_pilot_20261020/index.html) includes
225 requested object-band rows: 162 available images measured and 63 missing
inputs recorded. Raw C/A and both S profiles are retained. Only 33 images pass
the current corner-sky criterion, so corrected results are strongly selected.
The exported directory includes the checksum-bound HTML, Parquet measurements,
input inventory and 144 synthetic noise/PSF/mask experiments. Scientific
reference matching and DES-specific calibration remain necessary.

To inventory local upstream products and run a new pilot:

```bash
python scripts/prepare_des_pilot.py --output .outputs/des_inputs_new
python -m algorithmic_morphology run --config-file-name .outputs/des_inputs_new/config.py --chunk-size 15
python scripts/validate_sensitivity.py --output .outputs/sensitivity_new
python scripts/analyze_pilot.py --run-dir <completed-run-directory> --inputs .outputs/des_inputs_new --sensitivity .outputs/sensitivity_new --output <new-report-directory>
```

The current shared method is `cas_discrete_v2`: bilinear rotation excludes any
masked contributor, and raw/sky A use the same valid-pixel flux denominator.
Both original S profile parameter names remain unchanged.

For architecture and remaining work see [project design](docs/PROJECT_DESIGN.md).

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) and [AGENTS.md](AGENTS.md).
