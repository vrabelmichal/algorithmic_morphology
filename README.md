# Algorithmic Morphology

Explainable, literature-motivated measurements of galaxy-image morphology.
The first implementation measures **concentration (C)**, **raw and
sky-corrected 180° asymmetry (A)**, and **two named clumpiness/smoothness
(S) profiles** for each galaxy and photometric band.

**Status: complete 84,200-object upstream catalog processed; scientific calibration remains open.**
The profiles are documented adaptations, not exact numerical reproductions
of every implementation in the literature.

## Metrics produced

- `C = 5 log10(r80/r20)` in a Petrosian aperture.
- Aperture Gini/M20 and positive-flux axis ratio, ellipticity and image-coordinate position angle (see [definitions](docs/EXTENDED_METHODS.md)).
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
not the full catalog result. The [completed HTML analysis](reports/full_catalog_20261020/index.html)
contains all 252,600 requested rows: 247,800 with finite C and raw A, 4,704 with
invalid CAS results, and 96 explicit missing-cutout rows. Individual metrics have
separate validity counts; corrected A is available for 242,335 rows. The report
now reads as a scientific paper: executive summary, methods and citations,
19 grouped figures with robust shared ranges and full-range companions, a
searchable column dictionary, GALFIT comparisons and reproducibility appendices.
Per-object measurements and their hashes remain unchanged.

The GALFIT comparison uses the 84,061-object single-Sersic export from
`per_band_mag_reff__xy_2_5arcsec`, matched on exact ID and band, with declared
warning and constraint-boundary sensitivity cohorts. It relates morphology to
mean effective surface brightness, fitted radius, magnitude, Sersic index,
axis ratio and stored reduced chi-square. Every data source and filtering rule
is recorded in the HTML and `paper_v2_data/galfit_sources.json`.

To build this report revision from existing full-catalog measurements:

```bash
python scripts/build_paper_report.py --report-dir reports/full_catalog_20261020 --source-project ../nucleated_lsbg
```

The revision builder refuses to overwrite an existing revision. Revision 1 and
its validation record remain preserved as `report_v1.html` and `validation.json`;
revision 2 has `report_v2_manifest.json` and `report_v2_validation.json`.

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

The pilot used shared method `cas_discrete_v2`: bilinear rotation excludes any
masked contributor, and raw/sky A use the same valid-pixel flux denominator.
Both original S profile parameter names remain unchanged. The full run uses
`cas_discrete_v3`, adding a flagged initial-center fallback for a rare SciPy
Powell empty-direction exception; the formulas and sky criterion are unchanged.

For architecture and remaining work see [project design](docs/PROJECT_DESIGN.md).

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) and [AGENTS.md](AGENTS.md).

## Human and DECaLS10 context and embedding atlas

The [current HTML report](reports/full_catalog_20261020/index.html) is revision 3.
It preserves the full scientific report and adds 25 classification figures and
an interactive selector for 63 morphology panels on 21 saved PCA/UMAP spaces.
The report directory is portable: keep `paper_v3_assets/` beside the HTML.
Earlier v1/v2 HTML and validation records remain available.

The latest completed human products use the 14 October export: 33,164 labeled
objects, 32,288 with at least two morphology responses. Continuous comparisons
use raw equal-weight answer fractions; selected classes use the producer's
published consensus queries unchanged. The primary machine input is the newer
84,104-object, background-subtracted/masked, compatible Nano DECaLS10 run,
separate from the older 26,209-object table and human-label fine-tuned models.
Sources, vocabularies, denominators, counts and SHA-256 fingerprints are explicit.

All ten DECaLS10 probabilities and selected human morphology, LSB, visible-nucleus
and fit-quality answers are compared with C, both A outputs, both S profiles,
aperture Gini/M20 and positive-flux shape. The analysis includes matched cohorts,
≥3-vote sensitivity, class effect sizes and rank adjustment for fitted mean
surface brightness, radius and axis ratio. These are descriptive associations,
not calibrated morphology classes or classifier accuracy measurements.

Embedding figures are implemented and produced in `clustering_lsbg`; this project
copies its products and provenance. Saved coordinates remain unchanged. Exact
ID joins, all finite raw feature values, missing measurements and color tails
are independently audited in `report_v3_validation.json`.

Reproduction is staged so figures do not require rerunning image measurements:

```bash
python scripts/prepare_classification_analysis.py --sources /path/to/sources.json --output /new/frozen/input
python scripts/plot_classification_analysis.py --input /new/frozen/input
# In clustering_lsbg, run its render_algorithmic_morphology_overlays module.
python scripts/build_classification_report.py --input /new/frozen/input --embedding-manifest /path/to/clustering/manifest.json
python scripts/validate_classification_report.py --report-dir reports/full_catalog_20261020
```

The source JSON keys and current paths are recorded in
`paper_v3_data/manifest.json`. Use new output directories/revisions; builders
refuse to overwrite existing revisions. See [the analysis protocol](docs/CLASSIFICATION_CONTEXT.md).

Report revision 4 embeds the interactive correlation explorer directly in
**section 6.4**, with an **Interactive correlations** navigation link. Its class,
band and cohort selectors use the same frozen 2,040 association records as the
static figures. It is embedded through inline `srcdoc`, so it also works when
opening the report from disk. A separate `association_explorer.html` is provided.
Reproduce this presentation-only addition with
`python scripts/integrate_association_explorer.py` (new revision required if v4 exists).
