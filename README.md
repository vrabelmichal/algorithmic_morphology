# Algorithmic Morphology

Explainable, literature-motivated measurements of galaxy-image morphology.
The first implementation measures **concentration (C)**, **raw and
sky-corrected 180° asymmetry (A)**, and **two named clumpiness/smoothness
(S) profiles** for each galaxy and photometric band.

**Status: pilot implementation; needs validation on real DES LSBG images.**
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

**No DES results have been calculated in this PR.** Scientific validation
against known simulated galaxies and reference implementations remains
necessary before drawing astrophysical conclusions.

For architecture and remaining work see [project design](docs/PROJECT_DESIGN.md).

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) and [AGENTS.md](AGENTS.md).
