# Agent guidance

Read README.md, docs/CAS_METHODS.md, docs/CAS_LITERATURE_REVIEW.md
and docs/PROJECT_DESIGN.md before editing scientific code.

User-approved comparison: **keep A raw and sky-corrected** together,
plus **both** named S profiles (Lotz/statmorph-like and Conselice-like).
Do not drop an output, pick a winner, or silently change any definition.
The v1 methods are explicitly *literature-motivated approximations*, not
a claim of exact statmorph or historical Conselice compatibility.

Any method change requires a new method/profile version; log all numerical
choices, formulas, input fingerprints and quality flags. Preserve signed
sky-subtracted pixels and valid negative corrected indices. Mask
**contaminants only**, not the target or the sky. Prevent zeros in masked
star regions from creating smoothing halos. Flag invalid/no-sky results
rather than substituting zero or -99.

Respect the verified input/config conventions of nucleated-lsbg, but do
not silently reuse GALFIT r_eff as the Petrosian radius. Do not assume
all cutouts are sky-subtracted or properly oriented.

Make output partitions durably recoverable. Do not overwrite old runs by
default. Promote 'latest' only after the exact Parquet table and matching
SHA-256-bound HTML report have been committed. Timestamp stage/progress
logs to stdout and include per-object errors.

Before declaring scientific accuracy, test real DES LSBGs, masking and
orientation, noise and PSF effects, and matched reference code settings.
A discrepancy might be an error in the reference interpretation rather
than in the new implementation. Tests must be run, not merely written.
