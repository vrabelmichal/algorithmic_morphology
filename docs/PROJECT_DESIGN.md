# Proposed project architecture (2026-10-20)

**Historical architecture proposal.** A first working pilot now exists (see [implemented methods](CAS_METHODS.md) and README.md). Some sections below remain aspirational, particularly resume and independent reference validation. The supported CLI is `python -m algorithmic_morphology run --config-file-name <config.py>`.

## Goals and inputs

Run documented algorithmic morphology features per galaxy **and per band**, avoiding trained black-box classifications. Interpretability and exact reproducibility for fixed inputs and pinned software are goals; low-S/N systematic errors are still possible.

I inspected `vrabelmichal/nucleated-lsbg` on `main`, particularly:
- `configs/sky_subtracted_2_onmr2__detr_and_vit_all_multiband_256.py`
- `configs/sky_subtracted_2__detr_and_vit_all.py`
- `functions/load_data.py`
- `docs/galaxy_masking_procedure.md`
- `common/output.py`

Compatible configuration names worth retaining when semantics match: `lsb_catalog_fits_path`, `band_order`, `band_cutouts_dir_paths`, `fallback_band_cutouts_dir_paths`, `cutouts_dir_path`, `Npixels`, `coadd_object_id` / `COADD_OBJECT_ID`, and Path-based `base_output_dir_path`. Accept FITS or Parquet catalog input where practical. Reuse Python config modules and typed/path validation, but DO NOT force unrelated GALFIT options into a morphology schema.

The referenced pipeline has multi-band cutouts, a two-pass SEP background subtraction, background RMS images, segmentation maps and combined contaminant masks. Its documented final object-removed mask convention is **1 = excluded contaminant**, **0 = retained**, with the central galaxy preserved and the background remaining unmasked. The morphology processor should first **read these precomputed products** without re-running GALFIT or reconstructing sky. Resolve band-specific WCS/flip, FITS units, mask shape and orientation, source identity, and target center explicitly. Configuration radii based on fitted r_eff are *not* CAS Petrosian radii.

## Proposed CLI and config design

Proposed future commands, **not implemented**:

    python -m algorithmic_morphology validate --config-file-name configs/example.py
    python -m algorithmic_morphology run --config-file-name configs/example.py --limit 100
    python -m algorithmic_morphology report --run-dir .outputs/.../<run-id>/

Config categories: catalog, source identifiers, bands and image locations, precomputed sky/RMS and mask paths, science profile+version, parallelism/chunk size, output root, report options, logging. Preserve the *effective serialized configuration* and its source hash in each run. Validate all selected files/columns before the costly work, allow per-band missing inputs to become failure rows. No implicit calls to an external repository.

## Durable incremental Parquet and non-overwriting runs

Proposed immutable run layout (illustrative):

    .outputs/<dataset>/<config>/<UTCtimestamp>_<uniqueid>/
      effective_config.json
      manifest.json
      progress.json
      parts/part-000001.parquet
      parts/part-000002.parquet
      measurements__<run-id>.parquet
      report__<run-id>.html
    .outputs/<dataset>/<config>/latest -> <UTCtimestamp>_<uniqueid>/

Use exclusive directory creation and collision-safe UTC IDs. **No overwrite by default**, even after failure. Write independent temporary Parquet chunks at regular bounded intervals (for example, 50-200 object-band rows), then atomically rename each finished chunk. Parquet must *not* be naively appended in place. Periodically persist a small atomic checkpoint manifest with partition filenames, completed object-band keys and error counts. An interrupted run is inspectable; only explicit resume of the same incomplete run may reuse its directory, after verifying inputs/config and avoiding duplicate keys.

On completion, assemble a canonical Parquet via streaming row groups (without loading the entire sample into RAM), verify row counts and unique keys, atomically publish it, and create an HTML report from that **exact table**. Store table SHA-256 digest and filename in the report metadata and manifest. Then atomically mark the run complete and update the `latest` symlink **only** after both data and report succeed. A failed run never replaces the previous successful latest. Report re-generation never silently replaces the existing HTML and always preserves the data digest; `latest` must not be consulted to resolve a report for a specific historical run.

Proposed table row fields: run_id, coadd_object_id, band, method_id/version, image/mask/sky source info, processing status, error type, quality flags, C, r20/r80/R_p, aperture radius, A raw/sky/corrected, S raw/sky/corrected, asymmetry center x/y, relevant mask fraction, sky and S/N diagnostics. Use **null** for invalid measurements, not -99. Keep even valid negative sky-corrected measurements and mark confidence/quality separately. Each metric should publish its exact method parameters in the manifest.

## Progress and robust processing

Python logging to **stdout** by default with ISO-8601 UTC timestamp and log level, run ID and stage. Include input validation, start/finish of each major stage, counts, progress/throughput, chunk writes, failures, analysis and final completion. Store per-object failures without stopping other bands. Clearly distinguish *completed with failed objects* from *all successful* and from *interrupted*; do not silently hide failures by filtering those rows. Optional mirror file may be added; stdout is canonical.

Atomic operations should guard against process crashes: no partial Parquet file should be considered complete, no partial HTML should be published, and the manifest and `latest` should never point to an uncommitted product. Record tool versions, code git SHA, effective settings and input file identifiers/fingerprints for reproducibility.

## Automatic minimal HTML report

A self-contained HTML document for each completed run. Show run ID, exact table filename and SHA-256, input/config/method provenance, object-band counts, completion/flag/error summaries and **simple distribution histograms** of C, A, S per band. Add pairwise scatter plots C–A, A–S and C–S where useful. Label the number of valid points on each plot. Make failed and invalid fractions visible, and avoid presenting a survivor-only distribution as the full population. Never mix data from different runs or unlabelled method variants. Data outputs and reports should be named so they are visually associated, but checksum-based association is the robust guarantee.

## Implementation phases after scientific agreement

1. Validate sample real cutouts and masks in g/r/i: units, WCS, origin, sky, orientation, target retention.
2. Implement pure documented CAS numerical functions with synthetic fixtures and matched statmorph cross-checks. Version the scientific method, including **all** choices that change output.
3. Build catalog/band runner with checkpointed Parquet, failure handling and crash/restart/non-overwrite tests.
4. Build self-contained HTML report, table-hash validation and simple visualizations.
5. Run a representative faint-galaxy pilot, document systematic effects/failures; only then scale out.

Later extensions can include Gini/M20, shape/outer asymmetry, MID and other interpretable features, with separate literature-backed method profiles. No trained classifier, GALFIT refitting, overwriting or silently changing scientific definitions is in scope for the initial implementation.

## Pilot execution update (2026-10-20)

The [availability-selected DES report](../reports/des_pilot_20261020/README.md)
now contains 75 candidates, 162 measured image-mask pairs and 63 explicit
missing-band rows. It includes exact input hashes, an archived source snapshot,
144 synthetic noise/PSF/mask stress cases, finite-only summaries and comparison
plots. Shared method v2 fixes bilinear mask support and the A sky denominator;
both named S profiles remain v1. Report failures preserve the previous latest,
and input WCS disagreement is rejected. Independent reference-code agreement,
DES-specific uncertainty calibration and resume remain future work. Corner-sky
coverage (33/162 images) is the principal observed pilot limitation.
