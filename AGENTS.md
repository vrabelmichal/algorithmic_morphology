# Agent guidance

Read README.md, docs/CAS_LITERATURE_REVIEW.md and docs/PROJECT_DESIGN.md.

Scientific methods are **not yet agreed**. Do not silently implement one of several inconsistent CAS variants. Record the equation, source, exact aperture and Petrosian rules, center optimization, sky-noise correction, smoothing scale, central exclusion, mask behavior and normalization for each version. Require synthetic and matched-reference tests before adopting results. Fixed-algorithm results can still be biased by noise or masks.

Use precomputed target-preserving contaminant masks where valid (background unmasked; other bright sources excluded). Symmetrize for asymmetry and avoid artifacts from masked pixels during filtering. Confirm shape, WCS/orientation and pixel scale. Do not substitute GALFIT r_eff for Petrosian R_p without an explicitly different named statistic.

Never overwrite completed run data by default. Persist atomic Parquet checkpoints during processing. Include quality flags and per-object error rows. Distinguish failure/NaN from measured zero; avoid -99 sentinel values and preserve valid negative corrected indices. Timestamp major-stage logs to stdout. Promote a latest symlink only for a fully completed table+HTML report with exact table checksum/provenance; do not mix runs.

Keep compatible input/config names from nucleated-lsbg when semantics match. Do not import its entire GALFIT pipeline or assume local input paths exist for every developer. Do not state tests have passed unless run.
