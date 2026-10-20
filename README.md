# Algorithmic Morphology

Reproducible, explainable features of galaxy-image morphology. The first
scientific target is non-parametric **concentration, asymmetry, and
clumpiness/smoothness (CAS)** measured for each galaxy and image band.

**Status: literature and architecture review.** There is **no CAS
implementation or executable processing pipeline yet**. The literature
contains inconsistent definitions, especially for asymmetry normalization
and clumpiness smoothing and scaling. Agreeing on the exact scientific
method before implementing it is intentional.

- [CAS literature, equations, conflicting methods and open decisions](docs/CAS_LITERATURE_REVIEW.md)
- [Proposed configuration, inputs, incremental Parquet, reporting and logging](docs/PROJECT_DESIGN.md)
- [Agent development guidance](AGENTS.md)

The intended processor will reuse DES galaxy cutout/config/mask conventions
from [nucleated-lsbg](https://github.com/vrabelmichal/nucleated-lsbg)
where appropriate. Results will be checkpointed as Parquet during processing;
outputs will be non-overwriting, associated with one immutable timestamped
run, and accompanied by a self-contained HTML distribution/QA report tied
to exact output file checksums. Standard output will carry timestamped
progress and errors.

Next steps: agree on the scientific CAS profile, build and validate
measurements, then implement the catalog runner and report. No
trained morphology classification is planned for this initial scope.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Discuss substantial choices before
implementing; cite literature and preserve method/input provenance.
