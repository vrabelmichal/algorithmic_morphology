# Complete upstream catalog morphology

Open [the HTML analysis](index.html).

84,200 objects; 252,600 requested object-band rows.

Full canonical table: `/mnt/data_sm870qvo/lsbg/analysis/algorithmic_morphology/.outputs/full_catalog/20261020T181011149128Z_48a36fe9/measurements__20261020T181011149128Z_48a36fe9.parquet`.
SHA-256: `4505526fbb05f410dcba1a3ac0063e684191eea19a40b376a77b21f774c3d9c7`.

Compact table: `measurements_metrics.parquet` (separate digest in `table_provenance.json`).

All 84,200 upstream IDs appear once in each of g/r/i. Status totals: 247,800 ok,
4,704 invalid CAS rows and 96 missing cutouts. An invalid CAS row can retain
usable extended diagnostics. The HTML report gives each metric's finite sample,
measured findings, comparison plots and scientific interpretation limits.

Verified canonical/compact table hashes, exact source snapshot hashes, numerical
identity of the compact export, exact IDs and unique keys, finite metric counts
and all 48 embedded plot images of revision 1. See `validation.json`.


## Scientific report revision 2

The default `index.html` is now the paper-style [revision 2](report_v2.html):
executive summary; methods with typeset equations and references; grouped,
comparably scaled distributions; full-range companions and tail counts;
GALFIT relationships; discussion; searchable column definitions and appendices.
It has 19 main grouped figures plus four full-range distribution companions.

The primary GALFIT table contains 84,061 unique IDs from the per-band-radius,
single-Sersic configuration. The smaller 5,406-ID radius-order-1 variant is a
separate configuration-sensitivity comparison. Exact join/quality cohorts,
source fingerprints, correlations, binned trends and plot limits are in
`paper_v2_data/`. No CAS/S definition or morphology table value changed.

`validation.json` and its HTML digest describe preserved `report_v1.html`.
The current revision has its own `report_v2_manifest.json` and
`report_v2_validation.json`; `report_history.json` maps revision files and hashes.
