# Classification context: report revision 3

This extension consumes completed producer artifacts, not new classifier runs.
It leaves all 252,600 morphology rows and all scientific method definitions
unchanged. It adds classification association protocol v1 and the clustering
project's independently versioned external overlays.

## Frozen source selection, 22 October 2026

- Human vote mixtures: `classification_analysis/outputs/label-mixtures/20261022T100624Z`.
  Human counts/fractions only. Its older Zoobot columns are deliberately ignored.
- Published consensus: `classification_analysis/outputs/classification-analysis/20261014T170523`.
  Exact object-level membership and response-count products are reused.
- DECaLS10: `zoobot_classification/.outputs/classification/convnext_nano_aplpy_irg_linear_bgsub_masked_obj302317964_20260617`.
  84,104 valid ten-category vectors; matching Nano checkpoint/configuration.
  Tiny/base directory alternatives record a checkpoint/config-encoder mismatch;
  curated subsets are not substituted for full predictions.

These products join exactly on string object IDs, never row positions or a
numeric-ID fallback. All 84,200 morphology IDs remain, including missing labels,
failed images and metric-specific nulls. Human tables contain 33,164 voted
objects; at least two morphology responses leave 32,288. All of these have
matching current DECaLS10 predictions. Raw human morphology counts independently
match between the latest mixture table and main-report records.

An adjacent human run-metadata JSON was updated during this work. Its current
copy was frozen, with both the initial inventory digest and snapshot digest
recorded. All actual feature/membership tables retain their original hashes.
No personal reviewer identifiers, email addresses or comment text are exported.

## Semantics and estimands

Human answer 0 means **Not sure (Irr/other)**, not purely irregular. The other
answers are Featureless, LTG/spiral and ETG/elliptical. Raw fractions are equal-
weight responses with question-specific denominators; fewer than two answers
becomes missing, not zero. Class distribution queries preserve the producer's
seniority-weighted morphology/LSB plurality (senior weight 1.5) and explicit
splits. The two representations are labeled distinctly.

Machine probabilities retain their original numerical values. Validation
requires all ten values finite and within [0,1], with sum within 1e-4 of one.
No renormalization, replacement of invalid values or human/machine class merger
is performed. The model class names/order come from the local producer metadata;
the public checkpoint card does not itself document the label mapping.

Each correlation is per-band Spearman correlation on all finite eligible pairs.
The four cohorts are available pairs, matched ≥2-human/valid-machine pairs,
≥3 human morphology votes, and GALFIT-unflagged objects with finite controls.
The last also reports correlation of linearly residualized ranks after
adjustment for fitted mean effective surface brightness, effective radius and
axis ratio. This is a descriptive sensitivity check, not causal adjustment.

Selected-class effects are rank-biserial: `2U/(n_class*n_other)-1`. Human
nonmembers are only objects with sufficient comparable answers, including
other classes and splits. Machine top-class comparisons also export a maximum-
probability ≥0.5 sensitivity. Every effect includes class/other finite counts
and the selected class's median and 16–84% range. Rare classes are highlighted.

Probability trends use ten fixed bins, at least 20 finite pairs, with median
curves and 16–84% object spreads. No significance, calibration, held-out accuracy
or physical class frequency is inferred. Probabilities are compositional;
human fractions are discrete and commonly based on only two or three votes.

Both A values and both raw/corrected S profiles remain. Ellipticity equals
`1-q_positive`, so its Spearman correlations reverse the q signs exactly.
Image-coordinate position angle is axial and orientation-dependent; ordinary
scalar correlations would be misleading and are not included.

## Embedding ownership

`clustering_lsbg` owns rendering on the supplied GALFIT and external feature
output trees. There are seven feature families, with PCA and UMAP at 15/120
neighbors, each with 59,780 loose-selection objects. All 21 current 2D spaces
receive g/r/i panels and raw-valued browser exports. Saved coordinates are never
refit or changed. N-dimensional BIRCH inputs and archived histories are excluded
from this 2D replay and explicitly documented.

The report copies 63 upstream PNGs byte-for-byte. Each has shared robust feature
color ranges, grey missing points and finite/tail counts. DECaLS10 and human-
trained v3/v4 feature spaces are labeled separately. Supervised-space agreement
with training labels is not an independent validation. Projection geometry,
selection, shared images, noise, PSF and overplotting limit physical conclusions.

## Validation and preservation

The v3 validator verifies the unchanged morphology digest, previous HTML digest,
all source/product hashes, exact human response counts and consensus membership,
unaltered model vectors, independently recomputed available-pair correlations,
histogram tails, all 63 exported coordinates/features and PNG validity, the full
new column dictionary, HTML anchors and gallery assets. Float32 embedding CSVs
are verified by exact round-trip to their original float32 dtype; float64
coordinates use 1e-13 absolute/relative tolerance. The image arrays used for
rendering remain the original saved coordinate arrays.

The full morphology suite passes 26 tests. Clustering adds strict identity and
browser-export tests and retains the existing saved-alignment/export regressions.
`report_v3_manifest.json`, `report_v3_validation.json` and `report_history.json`
bind report revision, input table and reproducibility evidence. v1/v2 reports
and their validation products remain preserved.
