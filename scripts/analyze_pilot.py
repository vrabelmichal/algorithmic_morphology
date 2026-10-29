"""Export a portable, checksum-bound pilot analysis from one completed run."""
import argparse
from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
from algorithmic_morphology.report import summarize, sha256_file
from algorithmic_morphology.metrics import PROFILES


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--sensitivity',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--run-source',type=Path,help='Historical source snapshot; defaults to installed package')
    args=parser.parse_args()
    manifest=json.loads((args.run_dir/'manifest.json').read_text())
    table=args.run_dir/manifest['table']
    if sha256_file(table)!=manifest['table_sha256']:
        raise ValueError('Table digest mismatch')
    import algorithmic_morphology
    source_root=args.run_source or Path(algorithmic_morphology.__file__).parent
    for name,digest in manifest['provenance']['code_sha256'].items():
        if sha256_file(source_root/name)!=digest:
            raise ValueError(f'Code differs from run: {name}; pass --run-source with its archived source')
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'run_source').mkdir()
    for name in manifest['provenance']['code_sha256']:
        shutil.copy2(source_root/name,args.output/'run_source'/name)
    for source in [table,args.run_dir/manifest['report'],args.run_dir/'manifest.json',args.run_dir/'effective_config.json',
                   args.inputs/'inventory.csv',args.inputs/'selection.json',args.inputs/'config.py',args.sensitivity/'sensitivity.csv']:
        shutil.copy2(source,args.output/source.name)
    frame=pd.read_parquet(table)
    summary=summarize(frame)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2))
    stats=pd.DataFrame(summary['statistics'])
    stats.to_csv(args.output/'statistics.csv',index=False)
    comparisons=[]
    for band,part in frame.groupby('band'):
        for suffix in ['raw','corrected']:
            lotz=f'S_{PROFILES[0].name}_{suffix}'; cons=f'S_{PROFILES[1].name}_{suffix}'
            both=part[[lotz,cons]].dropna()
            comparisons.append(dict(band=band,comparison=f'S {suffix}',matched=len(both),
                                    spearman=both[lotz].corr(both[cons],method='spearman'),
                                    median_conselice_div10_minus_lotz=(both[cons]/10-both[lotz]).median()))
    pd.DataFrame(comparisons).to_csv(args.output/'profile_comparisons.csv',index=False)
    band_rows=[]
    for band,part in frame.groupby('band'):
        band_rows.append(dict(band=band,requested=len(part),available=int(part.image_file.notna().sum()),
                              C=int(part.C.notna().sum()),A_raw=int(part.A_raw.notna().sum()),
                              A_corrected=int(part.A_corrected.notna().sum()),
                              S_lotz_corrected=int(part[f'S_{PROFILES[0].name}_corrected'].notna().sum()),
                              S_conselice_corrected=int(part[f'S_{PROFILES[1].name}_corrected'].notna().sum())))
    counts=pd.DataFrame(band_rows)
    counts.to_csv(args.output/'counts.csv',index=False)
    sensitivity=pd.read_csv(args.sensitivity/'sensitivity.csv')
    smetrics=['C','A_raw','A_corrected']+[f'S_{p.name}_{v}' for p in PROFILES for v in ['raw','corrected']]
    synth=sensitivity.groupby(['morphology','psf_sigma','noise_sigma','mask_mode'])[smetrics].median()
    synth.to_csv(args.output/'sensitivity_summary.csv')
    cross=[]
    for metric in smetrics:
        pivot=frame.pivot(index='coadd_object_id',columns='band',values=metric)
        for b1,b2 in [('g','r'),('r','i'),('g','i')]:
            both=pivot[[b1,b2]].dropna()
            cross.append(dict(metric=metric,bands=f'{b1}-{b2}',matched=len(both),
                              median_difference=(both[b1]-both[b2]).median(),
                              spearman=both[b1].corr(both[b2],method='spearman')))
    pd.DataFrame(cross).to_csv(args.output/'band_comparisons.csv',index=False)
    (args.output/'README.md').write_text(f'''# DES CAS availability-selected pilot — 2026-10-20

Run `{manifest['run_id']}`, method `cas_discrete_v2`; both original S profile names retained.
Table SHA-256: `{manifest['table_sha256']}`.

The sample is the union of 75 candidate IDs with precomputed products in the local
nucleated_lsbg legacy-fliplr object-centered 1024-pixel output tree. It is an
availability sample from the 84,200-row candidate catalog, not a representative
LSBG population or a confirmed-galaxy sample. No CAS-based selection was made.
All 225 requested object-band rows are retained, including 63 absent images.
Source files were read without modification or an additional image flip.

## Coverage

```
{counts.to_string(index=False)}
```

`ok` means finite C and raw A, not a guarantee of a usable sky correction or
scientific accuracy. The HTML report includes all variants, finite-only sample
quantiles, missing counts and flag frequencies. `statistics.csv` gives p16,
median and p84; these describe this selected sample, not measurement errors.

## Profile comparisons

```
{pd.DataFrame(comparisons).to_string(index=False,float_format=lambda v:f'{v:.5g}')}
```

The division by ten above removes only the Conselice normalization. Remaining
differences also reflect the kernel size and excluded central radius; they do
not establish which profile is better. Correlations use paired finite rows.
`band_comparisons.csv` likewise uses only common object IDs in each band pair;
independent g/r/i distributions have different sample membership.

## Input and quality checks

All 162 existing image-mask pairs have identical target WCS positions, target
positions within 0.5 pixel of the stamp geometric center, and an unmasked target
pixel. The runner also checks WCS agreement on a 3×3 grid across each stamp.
The masks contain 0/1/2 values in g/r (overlapping upstream mask contributions)
and 0/1 in i; all nonzero values are contaminants. WCS agreement establishes
alignment between stored arrays, not an independent astrometric calibration.
The upstream legacy horizontal flip has not been reinterpreted or reapplied.
FITS DATE-OBS normalization warnings occurred; no data arrays were altered.

Flags: `{json.dumps(summary['quality_flags'],sort_keys=True)}`.
Errors: `{json.dumps(summary['errors'],sort_keys=True)}`.

## Controlled sensitivity experiment

`sensitivity.csv` contains 144 deterministic synthetic measurements: symmetric
Gaussian and offset-clump scenes, PSF sigma 0/1/2 pixels, noise sigma 0/0.5/2/5,
three random seeds and an injected masked contaminant variant. The target
Gaussian has peak 40 and sigma 12 pixels; the offset clump has peak 12, sigma 2.
All metrics, including Petrosian radius and center, are remeasured in every case.
`sensitivity_summary.csv` gives medians across the three seeds. These are stress
experiments with known scenes, not a DES noise or PSF calibration. Three seeds
are insufficient for quantitative uncertainty estimates. A masked scene also
removes some genuine galaxy light under the contaminant footprint.

## Interpretation limits and next validation

The corner-only sky estimator often cannot find a ≥95% valid 32-pixel tile in
these masked fields. Missing corrections remain null. The corrected subset
therefore has additional selection bias. Negative corrected values are retained.
The cleanest-corner choice, integer-aligned sky rotation versus subpixel galaxy
rotation, noise, PSF, target-mask extent and discrete apertures may bias values.
No morphology classes or published CAS thresholds are assigned. Both A outputs
and both S profiles remain necessary for subsequent calibration.

Matched statmorph/historical implementation agreement and DES-specific injection
recovery, PSF matching and uncertainty calibration remain unperformed. These are
required before physical classification or population conclusions. The present
report establishes end-to-end execution and exposes limitations of this pilot.

## Reproduction

Install this repository and its declared test dependencies, then run:

```bash
python scripts/prepare_des_pilot.py --output .outputs/des_inputs_new
python -m algorithmic_morphology run --config-file-name .outputs/des_inputs_new/config.py --chunk-size 15
python scripts/validate_sensitivity.py --output .outputs/sensitivity_new
python scripts/analyze_pilot.py --run-dir <completed-run-directory> --inputs .outputs/des_inputs_new --sensitivity .outputs/sensitivity_new --output <new-report-directory>
```

The preparation script accepts `--source-project` and `--catalog` for other local
locations. The exported config records the actual local paths. Run manifests
record package versions, code hashes, catalog/config fingerprints and numerical
settings; table rows contain image/mask SHA-256. The run began with uncommitted
fixes, so code hashes, rather than its base git commit alone, identify the code.
''')
    print(counts.to_string(index=False))
    print(json.dumps(summary['quality_flags']))


if __name__=='__main__':
    main()
