"""Publish the complete-catalog HTML findings and compact review artifacts."""
import argparse
from html import escape
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from algorithmic_morphology.metrics import PROFILES
from algorithmic_morphology.extended import EXTRA_METRICS
from algorithmic_morphology.report import sha256_file, summarize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    parser.add_argument('--inputs',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    manifest=json.loads((args.run_dir/'manifest.json').read_text())
    selection=json.loads((args.inputs/'selection.json').read_text())
    table=args.run_dir/manifest['table']
    if sha256_file(table)!=manifest['table_sha256']:raise ValueError('Table checksum differs')
    frame=pd.read_parquet(table)
    expected=selection['objects']*len(selection['coverage'])
    if len(frame)!=expected or frame.duplicated(['coadd_object_id','band']).any():raise ValueError('Catalog coverage invariant failed')
    from astropy.table import Table
    ids=Table.read(selection['catalog'])['COADD_OBJECT_ID']
    if set(frame.coadd_object_id)!=set(map(str,ids)):raise ValueError('Upstream IDs differ')
    args.output.mkdir(parents=True,exist_ok=False)
    for p in [args.run_dir/'manifest.json',args.run_dir/'effective_config.json',args.inputs/'selection.json',args.inputs/'coverage.csv',args.inputs/'config.py']:
        shutil.copy2(p,args.output/p.name)
    shutil.copytree(args.run_dir/'run_source',args.output/'run_source')
    # Keep source hashes and local full table in provenance; export compact metrics
    # for reviewing without duplicating unique image paths and file fingerprints.
    numeric=[c for c in frame if pd.api.types.is_numeric_dtype(frame[c])]
    columns=list(dict.fromkeys(['coadd_object_id','band','status','quality_flags','error_type','error_message','method_version','aperture_diagnostics_version','extended_quality_flags']+numeric))
    columns=[c for c in columns if c in frame]
    compact=args.output/'measurements_metrics.parquet'
    frame[columns].to_parquet(compact,index=False)
    summary=summarize(frame)
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2))
    stats=pd.DataFrame(summary['statistics']);stats.to_csv(args.output/'statistics.csv',index=False)
    metrics=['C','A_raw','A_corrected']+[f'S_{p.name}_{v}' for p in PROFILES for v in ['raw','corrected']]+EXTRA_METRICS
    counts=[]
    for band,part in frame.groupby('band'):
        for metric in metrics:
            vals=pd.to_numeric(part[metric],errors='coerce').replace([np.inf,-np.inf],np.nan)
            counts.append(dict(band=band,metric=metric,requested=len(part),valid=int(vals.notna().sum()),missing=int(vals.isna().sum())))
    counts=pd.DataFrame(counts);counts.to_csv(args.output/'metric_coverage.csv',index=False)
    status=frame.groupby(['band','status']).size().unstack(fill_value=0);status.to_csv(args.output/'statuses.csv')
    paired=[]
    for band,part in frame.groupby('band'):
        for suffix in ['raw','corrected']:
            x=f'S_{PROFILES[0].name}_{suffix}';y=f'S_{PROFILES[1].name}_{suffix}'
            both=part[[x,y]].dropna()
            paired.append(dict(band=band,comparison=f'S {suffix}',paired=len(both),spearman=both[x].corr(both[y],method='spearman'),median_conselice_div10_minus_lotz=(both[y]/10-both[x]).median()))
        for x,y in [('C','Gini_aperture_v1'),('Gini_aperture_v1','M20_aperture_v1')]:
            both=part[[x,y]].dropna()
            paired.append(dict(band=band,comparison=f'{x} vs {y}',paired=len(both),spearman=both[x].corr(both[y],method='spearman'),median_conselice_div10_minus_lotz=None))
    pd.DataFrame(paired).to_csv(args.output/'paired_comparisons.csv',index=False)
    bands=[]
    for metric in metrics:
        pivot=frame.pivot(index='coadd_object_id',columns='band',values=metric)
        for b1,b2 in [('g','r'),('g','i'),('r','i')]:
            both=pivot[[b1,b2]].dropna()
            bands.append(dict(metric=metric,bands=f'{b1}-{b2}',paired=len(both),median_difference=(both[b1]-both[b2]).median(),spearman=both[b1].corr(both[b2],method='spearman')))
    pd.DataFrame(bands).to_csv(args.output/'band_comparisons.csv',index=False)
    extended_flags=frame.extended_quality_flags.dropna().str.split('|').explode()
    extended_flags=extended_flags[extended_flags!=''].value_counts().to_dict()
    table_info=dict(full_table=str(table.resolve()),full_table_sha256=manifest['table_sha256'],
                    compact_table=compact.name,compact_table_sha256=sha256_file(compact),
                    requested_rows=expected,unique_objects=frame.coadd_object_id.nunique(),
                    matched_upstream_ids=True,source_sha256=selection['catalog_sha256'])
    (args.output/'table_provenance.json').write_text(json.dumps(table_info,indent=2))
    stats_html=stats[stats.metric.isin(metrics)].to_html(index=False,float_format=lambda x:f'{x:.5g}',na_rep='missing')
    latest=selection['catalog_config'];products=selection['product_config']
    intro=f'''<h1>Complete upstream catalog: morphology analysis</h1>
<p><strong>{selection['objects']:,} galaxies requested; {expected:,} object–band rows retained.</strong>
The object IDs exactly match the upstream catalog. There is no availability or morphology selection.
The earlier 75-ID report was an initial pilot and is superseded for coverage by this run.</p>
<h2>Dataset and configuration alignment</h2><p>The recent upstream configuration
<code>{escape(latest)}</code> and the 256-pixel configuration
<code>{escape(products)}</code> both select the same <code>concatenated-detr-vit-all.fits</code>
catalog. This run uses its complete original ordering and r/g/i bands. It reads the existing
sky-subtracted 256-pixel images and target-preserving masks shared by the upstream 256-pixel
configurations. The wide-image experimental products have different coverage and preprocessing;
those are not substituted for the baseline images. The symlink to the external data drive is
resolved explicitly. GALFIT models and fitted radii are not inputs to these measurements.</p>
<p>Input availability: {sum(r['paired'] for r in selection['coverage']):,} image–mask pairs;
{expected-sum(r['paired'] for r in selection['coverage']):,} missing combinations before measurement.
The statuses below also retain WCS, aperture, photometry and other failures.</p>
<h2>Processing outcomes</h2>{status.to_html()}
<h2>Every evaluated metric and its usable sample</h2>{counts.to_html(index=False)}
<p>C, raw/sky/corrected A, and both raw/sky/corrected S profiles are retained together.
Gini, M20, axis ratio, ellipticity and position angle are additional aperture diagnostics.
A null value is a failed or unavailable measurement, never an implicit zero.</p>
<h2>Additional metric definitions</h2><ul>
<li>Gini_aperture_v1 measures absolute flux inequality inside the unmasked circular CAS aperture.</li>
<li>M20_aperture_v1 measures the signed second moment of the brightest 20% of signed total flux,
relative to the signed total second moment at the flux centroid. Nonpositive moments are flagged.</li>
<li>The shape diagnostics use explicitly positive flux weights in that same aperture.
Axis ratio is the square root of the covariance eigenvalue ratio; ellipticity is 1−q.
The angle is in image coordinates, not celestial position angle, and is unstable for round sources.</li></ul>
<p>These aperture adaptations are not statmorph's segmentation-based Gini/M20. They have no
validated classification thresholds here. The formula background is documented in
<a href="https://arxiv.org/abs/astro-ph/0311352">Lotz et al. (2004)</a> and the
<a href="https://statmorph.readthedocs.io/en/latest/_modules/statmorph/statmorph.html">statmorph reference source</a>.
Exact formulas and limitations are in the repository's <code>docs/EXTENDED_METHODS.md</code>.</p>
<h2>Sample summaries</h2><p>p16–p84 describes the measured sample, not errors on an individual galaxy.
Use paired comparisons for cross-band statements; valid samples can differ by band and metric.</p>{stats_html}
<h2>Paired comparisons</h2>{pd.DataFrame(paired).to_html(index=False,float_format=lambda x:f'{x:.5g}',na_rep='—')}
<h2>Interpretation and remaining work</h2><p>No user decision prevented this processing.
Reference matching, realistic noise/PSF injection recovery and segmentation sensitivity remain
scientific validation tasks, not permission gates. The 144-case pilot stress test showed raw-A
noise sensitivity and negative corrected-A bias; it is not a DES uncertainty calibration. The
corner-only sky method remains unchanged under shared method v3, (v3 adds a flagged fallback for the SciPy empty-direction optimizer error), so the larger run is directly traceable without
silently replacing its scientific definition. Both S profiles remain v1.</p>
<p>Signed sky-subtracted values and valid negative corrections are retained. Gini is computed
from absolute pixel values; the shape diagnostics' positive weights can bias noisy images.
No morphology classes, merger labels, preferred S profile or population-causal conclusions
are assigned from these uncalibrated values.</p>
<p>Extended diagnostic flags: <code>{escape(json.dumps({str(k):int(v) for k,v in extended_flags.items()}))}</code>.</p>
<h2>Download and reproduce</h2><p><a href="measurements_metrics.parquet">Compact per-object measurements</a> ·
<a href="metric_coverage.csv">Metric coverage</a> · <a href="statistics.csv">Summary statistics</a> ·
<a href="band_comparisons.csv">Matched band comparisons</a> · <a href="table_provenance.json">Exact table provenance</a>.</p>
<p>The full canonical Parquet remains at <code>{escape(str(table.resolve()))}</code>, SHA-256
<code>{manifest['table_sha256']}</code>. The compact export omits file paths and fingerprints;
it has its own SHA-256 and is not presented as the same file. The archived source snapshot,
manifest and effective config identify the exact measurement code and inputs.</p>
<hr><h2>Run-bound distributions, quality flags and provenance</h2>'''
    original=(args.run_dir/manifest['report']).read_text()
    output=original.replace('<body>','<body>'+intro,1)
    output=output.replace('<title>Algorithmic morphology — profile comparison</title>','<title>Complete DES catalog — morphology analysis</title>')
    (args.output/'index.html').write_text(output)
    (args.output/'README.md').write_text(f'# Complete upstream catalog morphology\n\nOpen [the HTML analysis](index.html).\n\n{selection["objects"]:,} objects; {expected:,} requested object-band rows.\n\nFull canonical table: `{table.resolve()}`.\nSHA-256: `{manifest["table_sha256"]}`.\n\nCompact table: `measurements_metrics.parquet` (separate digest in `table_provenance.json`).\n')
    print(args.output.resolve())
    print(status.to_string())


if __name__=='__main__':main()
