"""Preserve report v2 and publish a portable v3 with class and embedding context.

Consumes clustering-owned artifacts; never renders or refits an embedding here.
"""
import argparse
from datetime import datetime,timezone
from html import escape
import json
from pathlib import Path
import re
import shutil

import numpy as np
import pandas as pd
from algorithmic_morphology.classification_data import HUMAN_LABELS,CLASS_LABELS
from algorithmic_morphology.report import sha256_file


def table(frame):
    return frame.to_html(index=False,border=0,classes='data',float_format=lambda x:f'{x:.3f}',escape=True)


def gallery(name,items,default=0):
    choices=''.join(f'<option value="{i}"'+(' selected' if i==default else '')+f'>{escape(x["title"])}</option>' for i,x in enumerate(items))
    item=items[default]
    return f'''<div class="analysis-gallery" id="{name}"><label for="{name}-select">Choose a figure</label>
<select id="{name}-select" data-gallery="{name}">{choices}</select>
<figure><a class="full-figure" href="{item['file']}"><img loading="lazy" src="{item['file']}" alt="{escape(item['title'])}"></a>
<figcaption>{escape(item['caption'])}</figcaption><p><a class="full-figure" href="{item['file']}">Open full-resolution figure</a></p></figure></div>'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--input',type=Path,required=True)
    p.add_argument('--embedding-manifest',type=Path,required=True);p.add_argument('--report-dir',type=Path,default=Path('reports/full_catalog_20261020'))
    a=p.parse_args();root=a.report_dir.resolve();source=a.input.resolve()
    if (root/'report_v3.html').exists():raise FileExistsError('Report v3 exists; preserve it and choose a new revision')
    data=root/'paper_v3_data';assets=root/'paper_v3_assets';data.mkdir(exist_ok=False);assets.mkdir(exist_ok=False)
    frozen=json.loads((source/'manifest.json').read_text());contract=json.loads((source/'feature_contract.json').read_text())
    if sha256_file(source/'joined.parquet')!=contract['table_sha256']:raise ValueError('Joined input hash mismatch')
    df=pd.read_parquet(source/'joined.parquet');assoc=pd.read_csv(source/'associations.csv')
    embedded=json.loads(a.embedding_manifest.read_text());items=json.loads((source/'figures.json').read_text())
    for item in items:
        src=Path(item['file']);dst=assets/src.name;shutil.copyfile(src,dst);item['file']=dst.relative_to(root).as_posix()
    embedding_items=[]
    for i,record in enumerate(embedded['embeddings']):
        src=Path(record['png'])
        if sha256_file(src)!=record['png_sha256']:raise ValueError('Embedding PNG digest mismatch')
        dst=assets/f'embedding_{i:02d}_{record["band"]}.png';shutil.copyfile(src,dst)
        embedding_items.append(dict(title=record['title']+' · '+record['projection']+' · '+record['band']+' band',
            file=dst.relative_to(root).as_posix(),caption=f'{record["rows"]:,} loose-selection galaxies; {record["matched"]:,} exact morphology ID matches. '+
            'All saved points are retained; grey means that metric is unavailable. Each panel has its own finite count and color-tail audit. '+
            'UMAP and PCA coordinates are dimensionless and are not physical distances. These are overlays on previously fitted spaces, not new morphology embeddings. '+
            'Human-trained v3/v4 spaces are separate models from the DECaLS10 classifier; label agreement in a supervised space is not independent validation.'))
    for filename in ['associations.csv','coverage.csv','class_effects.csv','probability_trends.csv','class_distribution_range_audit.csv','feature_contract.json','manifest.json']:
        shutil.copyfile(source/filename,data/filename)
    shutil.copyfile(a.embedding_manifest,data/'clustering_overlay_manifest.json')
    # Save only external labels, once per object. Morphology and GALFIT tables already exist in the bundle.
    morphology_cols=set(pd.read_parquet(root/'measurements_metrics.parquet').columns)
    labelcols=['coadd_object_id',*[c for c in df if c not in morphology_cols and not c.startswith('gf_')]]
    labels=df.loc[df.band=='r',labelcols].reset_index(drop=True);labels.to_parquet(data/'classification_labels.parquet',index=False)
    mapping=json.loads(Path(frozen['sources']['human_mapping']['path']).read_text());(data/'human_query_mapping.json').write_text(json.dumps(mapping,indent=2))
    (data/'figure_index.json').write_text(json.dumps({'classification':items,'embeddings':embedding_items},indent=2))
    defs=[]
    for key,label in HUMAN_LABELS.items():
        defs.append(dict(column=key,meaning=label,definition='Raw equal-weight answer fraction; question-specific denominator ≥2. Missing means insufficient/no responses; not zero.'))
    for i,label in enumerate(CLASS_LABELS):defs.append(dict(column=f'class_{i}',meaning='DECaLS10: '+label,definition='Supplied ten-way model probability; finite, within [0,1], vector sum within 1e-4 of one. Invalid vectors are null; no renormalization.'))
    defs.extend([dict(column='human_vote_count',meaning='Human morphology responses',definition='Latest named reviewer vote per object; anonymous rows separately retained by producer.'),
                 dict(column='top_morphology',meaning='Raw vote plurality with explicit ties',definition='Label-mixture producer output; not the seniority-weighted consensus used for class distributions.'),
                 dict(column='zoobot_top',meaning='Largest supplied class probability',definition='0–9 class index; maximum probability is not a human vote or proof of class membership.'),
                 dict(column='zoobot_valid',meaning='Probability-vector validation',definition='All ten entries finite and in [0,1], sum≈1; missing for absent predictions.'),
                 dict(column='normalized_entropy',meaning='Human vote-mixture entropy',definition='−Σp ln p / ln4. A discrete vote-disagreement summary, not calibration error.'),
                 dict(column='zoobot_top_probability',meaning='Maximum supplied model probability',definition='Confidence sensitivity at ≥0.5 is exported separately; no such cut in the main analysis.')])
    for c in labels:
        if c.startswith('analysis-'):defs.append(dict(column=c,meaning='Published human query membership',definition='Boolean copied unchanged from classification_analysis. Exact expression and min-response requirement in human_query_mapping.json.'))
        if c.startswith('count_'):defs.append(dict(column=c,meaning='Raw per-category human morphology count',definition='Numerator for the matching p_ fraction; not seniority-weighted.'))
        if c.endswith('_count') and c.startswith('p_'):defs.append(dict(column=c,meaning='Question-specific answer denominator',definition='Number of comparable human responses; fraction is null when fewer than two.'))
    pd.DataFrame(defs).to_csv(data/'classification_column_dictionary.csv',index=False)
    cov=pd.read_csv(source/'coverage.csv').set_index('band').loc[['g','r','i']].reset_index()
    r=df[df.band=='r'];counts=r.groupby('zoobot_top').size().reindex(range(10),fill_value=0)
    top=pd.DataFrame({'Class':CLASS_LABELS,'Top-class objects':counts.to_numpy(),'Share of valid predictions (%)':100*counts.to_numpy()/r.zoobot_valid.eq(True).sum()})
    def rho(metric,target,cohort='available',stat='spearman'):
        return float(assoc[(assoc.band=='r')&(assoc.metric==metric)&(assoc.target==target)&(assoc.cohort==cohort)][stat].iloc[0])
    findings=f'''<p><strong>Visible nucleus votes provide a clear example of structural association.</strong> In r, the visible-nucleus answer fraction correlates with aperture Gini (ρ={rho('Gini_aperture_v1','p_visible_nucleus'):.3f}), M₂₀ (ρ={rho('M20_aperture_v1','p_visible_nucleus'):.3f}) and concentration (ρ={rho('C','p_visible_nucleus'):.3f}). After rank adjustment for GALFIT mean effective surface brightness, radius and axis ratio, the first two become {rho('Gini_aperture_v1','p_visible_nucleus','galfit_adjusted','partial_rank'):.3f} and {rho('M20_aperture_v1','p_visible_nucleus','galfit_adjusted','partial_rank'):.3f}. Shared brightness/size/shape information explains part of the unadjusted association; the residual is still descriptive. Visible nucleus means the human answer to that question, not a physically confirmed nuclear star cluster.</p>
<p><strong>Class names are not interchangeable with individual indices.</strong> The r-band raw-asymmetry correlation with DECaLS10 disturbed probability is {rho('A_raw','class_0'):.3f}, versus {rho('A_corrected','class_0'):.3f} for corrected asymmetry. The ETG vote fraction has only ρ={rho('C','p_etg'):.3f} with concentration. Neither result supports replacing human or learned classes with a single fixed CAS threshold; corrected indices remain noise-sensitive and can be negative.</p>
<p><strong>Probability and top-class results answer different questions.</strong> A class can correlate with morphology while rarely being the largest predicted probability. The cigar-shaped smooth class is top for only {counts.loc[4]:,} objects and the unbarred tight-spiral class for {counts.loc[6]:,}. Their hard-class effects are fragile. Disturbed and round-smooth predictions dominate this supplied sample, but their proportions are model outputs under domain shift, not measured physical class frequencies. Probabilities are compositional: increasing one necessarily changes the available mass for the others.</p>'''
    methods='''<p>Associations use each object once within each band and all finite eligible pairs. Spearman ρ is the Pearson correlation of midranks, so monotonic relationships and ties are handled without assuming Gaussian errors. No pooled-band correlation treats repeated observations of one object as independent. The four views are: all available pairs; the matched ≥2-human-vote/valid-DECaLS10 sample; ≥3 morphology votes; and a GALFIT-unflagged control cohort with linear residual correlation of ranks after adjusting for mean effective surface brightness, effective radius and axis ratio. The export includes the unadjusted correlation on that same control cohort. Partial rank adjustment removes only the fitted rank-linear component, not all confounding, and is not causal.</p>
<p>For selected classes, rank-biserial separation is 2U/(n₁n₀)−1: positive means that a randomly selected class member tends to have a larger metric than an eligible nonmember; ties count one half in U. Human nonmembers include other answers and split outcomes, but never objects without ≥2 comparable responses to the relevant question. The model view uses the largest probability and separately exports a ≥0.5-maximum-probability sensitivity. No p-value, accuracy estimate or independent predictive test is inferred from this exploratory survey of many associations.</p>
<p>Probability trends use fixed bins of width 0.1, minimum 20 finite pairs, median lines and the 16–84% object spread. Quantiles use all finite measurements before display cropping. Human fractions are discrete and dominated by two or three responses; the ≥3-vote view is a sample-selection check, not a monotonic certainty guarantee. All images use the retained v3 CAS measurements and both v1 S profiles without redefining them. Positive-flux ellipticity is exactly 1−q₊, so its rank associations are the sign reversal of q₊ and are not repeated. Image position angle is axial and depends on image orientation; ordinary linear correlations with that angle are intentionally not interpreted as class associations.</p>'''
    source_table=pd.DataFrame([
        {'Product':'Human vote mixtures','Run / input':'22 Oct 2026 analysis; 14 Oct export','Objects':33164},
        {'Product':'Published human query membership','Run / input':'14 Oct 2026 17:05:23 analysis','Objects':84200},
        {'Product':'DECaLS10 Nano predictions','Run / input':'bgsub + masked, obj302317964, 17 Jun 2026','Objects':frozen['zoobot_input_rows']},
        {'Product':'Saved embedding membership','Run / input':'Loose GALFIT / external-feature trees','Objects':embedded['embeddings'][0]['rows']},
    ])
    classified=f'''<section id="classification"><h2>6. Human labels and DECaLS10 classification</h2>
<h3>6.1 Sources, vocabulary and coverage</h3><p>We use the newest completed human-label products located at this analysis freeze and the newest compatible Nano DECaLS10 background-subtracted/masked run. The ten-category checkpoint is <code>mwalmsley/zoobot-finetuned-decals10</code>, distinct from the human-label fine-tuned v3/v4 models. Its class names and class order are those recorded by the local producer, whose mapping is preserved; the public model card itself has no descriptive label mapping. See <a href="#ref8">the model source [8]</a> and <a href="#ref9">Zoobot documentation [9]</a>.</p>{table(source_table)}
<p>The latest label-mixture report used an older 26,209-object DECaLS10 table. We reuse only its human counts/fractions and ignore its class_ columns, joining the newer 84,104-vector table directly. Nano/tiny/base directory names alone are insufficient provenance: the alternative tiny/base metadata record a checkpoint/config-encoder mismatch. We choose the matching Nano run, retain all predictions, and do not silently substitute curated subsets or human-trained models. Every available candidate metadata file and the source hashes are recorded.</p>{table(cov)}
<p>Human morphology answers are Featureless (−1), Not sure (Irr/other) (0), LTG/spiral (1), and ETG/elliptical (2). The 0 answer combines uncertainty and irregularity and is not a pure irregular class. Continuous views use equal-weight raw fractions; the class-distribution views copy the producer's published seniority-weighted plurality selections (senior weight 1.5 for morphology/LSB) unchanged, with explicit split outcomes. Visible-nucleus and failed-fitting consensus have their producer-specific definitions. At least two comparable responses are required. These four human categories are not aliases for the ten machine categories.</p>
<details><summary>DECaLS10 top-class coverage and rare classes</summary>{table(top)}</details>
<h3>6.2 Analysis definitions and sensitivity checks</h3>{methods}
<h3>6.3 Measured associations and interpretation</h3>{findings}
<h3>6.4 Correlation and class-separation atlas</h3>{gallery('class-matrices',items[:5])}
<h3>6.5 Human consensus distributions</h3>{gallery('class-distributions',items[5:8],default=1)}
<h3>6.6 Probability and answer-fraction trends</h3>{gallery('class-trends',items[8:])}
<details><summary>Column meanings, source lineage and full analysis tables</summary>{table(pd.DataFrame(defs))}
<p>Raw reviewer names, emails and comments are not included. Human labels are once-per-object; the separate morphology table supplies g/r/i rows.</p>
<ul>{''.join(f'<li><a href="paper_v3_data/{name}">{escape(label)}</a></li>' for name,label in [('classification_labels.parquet','Object-level labels, raw fractions and published memberships'),('associations.csv','All 2,040 finite-pair correlation / sensitivity records'),('class_effects.csv','Class contrasts, exact sample counts and quantile summaries'),('probability_trends.csv','Probability-bin counts, medians and percentile spreads'),('class_distribution_range_audit.csv','Every class-histogram display range and tail count'),('classification_column_dictionary.csv','New column dictionary'),('human_query_mapping.json','Exact published human query definitions'),('manifest.json','Input paths, file hashes and source selection inventory')])}</ul></details></section>'''
    embedding_section=f'''<section id="embedding-morphology"><h2>7. Morphology on saved clustering spaces</h2>
<p>The plots in this section are created by <code>clustering_lsbg</code>, using its existing saved-embedding replay pattern. Every current PCA-2D and UMAP-2D artifact in the two requested output trees is included: 7 feature families × 3 projections (PCA, UMAP with 15 neighbors, UMAP with 120 neighbors) × 3 photometric bands = 63 panels, each containing all 10 morphology features. N-dimensional BIRCH inputs and archived projection-history copies are not treated as new two-dimensional spaces. No embedding is recomputed, no row-order match is assumed, and no morphology value is used to change the coordinates.</p>
<p>Each embedding has {embedded['embeddings'][0]['rows']:,} loose-selection galaxies, a subset of the full 84,200-object morphology catalog. Comparisons across maps therefore concern that saved selection. The maps retain all saved coordinates and use exact index-to-object joins. Point colors use pooled g/r/i robust limits; missing measurements remain grey. Corresponding Conselice-like S color limits are exactly ten times Lotz-like limits to expose normalization. Color-tail counts are annotated and exported; no data value is changed. A fixed random drawing order (seed 42) limits systematic overplot order, but dense scatter remains subject to occlusion.</p>
<p><strong>Read the geometry cautiously.</strong> GALFIT-only spaces reflect fitted size, brightness and shape, while external spaces include DECaLS10 probabilities/representations or human-trained v3/v4 features. Feature-family agreement can reflect shared images and training information. In particular, a human-trained embedding is not an independent test against the labels used in training. UMAP distances, empty regions and orientation are projection-dependent; morphology gradients are visual context, not cluster significance or physical class boundaries.</p>
{gallery('embedding-atlas',embedding_items,default=next(i for i,x in enumerate(embedded['embeddings']) if x['family']=='galfit_g_i' and x['band']=='r' and '120' in x['projection']))}
<details><summary>Clustering output locations, hashes and export conventions</summary><p>Original PNGs and per-band audit manifests remain beside each <code>embedding.parquet</code> under <code>embedding_overview_external/algorithmic_morphology/full_catalog_20261022_v2/&lt;band&gt;/</code>. Raw coordinate/metric exports and JSON continuous-column definitions use the established <code>embedding_exports/embedding_overview_external/algorithmic_morphology/</code> hierarchy. The report copies the clustering-produced PNGs byte-for-byte into its portable asset folder.</p><p><a href="paper_v3_data/clustering_overlay_manifest.json">Complete upstream overlay manifest, source/coordinate/image/export SHA-256 and all finite/tail counts</a></p><pre>{escape(str(a.embedding_manifest.resolve()))}</pre></details></section>'''
    base=(root/'report_v2.html').read_text()
    html=base.replace('scientific analysis · revision 2','scientific analysis · revision 3')
    html=html.replace('<h2>6. Discussion','<h2>8. Discussion').replace('6 · Discussion','8 · Discussion')
    html=html.replace('<a href="#discussion">','<a href="#classification">6 · Human / DECaLS10</a><a href="#embedding-morphology">7 · Embedding atlas</a><a href="#discussion">',1)
    html=html.replace('<section id="discussion">',classified+embedding_section+'<section id="discussion">',1)
    extra=f'''<p><strong>Classification extension, 22 October 2026.</strong> The analysis now relates all retained morphology features to 33,164 human-labeled galaxies (32,288 with ≥2 votes), 84,104 DECaLS10 prediction vectors, and all 21 requested saved PCA/UMAP spaces across g/r/i. Nucleus visibility tracks light concentration, while single-index class mappings remain weak or confounded. <a href="#classification">Class results</a> · <a href="#embedding-morphology">Embedding atlas</a>.</p>'''
    html=html.replace('<section id="summary"><h2>Executive summary</h2>','<section id="summary"><h2>Executive summary</h2>'+extra,1)
    html=html.replace('</body>',f'''<script>const scientificGalleries={json.dumps({'class-matrices':items[:5],'class-distributions':items[5:8],'class-trends':items[8:],'embedding-atlas':embedding_items})};
document.querySelectorAll('select[data-gallery]').forEach(select=>select.addEventListener('change',()=>{{const box=document.getElementById(select.dataset.gallery),item=scientificGalleries[select.dataset.gallery][Number(select.value)];box.querySelector('img').src=item.file;box.querySelector('img').alt=item.title;box.querySelector('figcaption').textContent=item.caption;box.querySelectorAll('.full-figure').forEach(a=>a.href=item.file);}}));</script></body>''')
    html=html.replace('</style>','select{max-width:100%;padding:.7em;font:inherit;background:var(--card,#fff);color:inherit;border:1px solid var(--border,#aaa);border-radius:6px}.analysis-gallery label{display:block;margin-bottom:.5em}.analysis-gallery img{width:100%;height:auto}.analysis-gallery{margin-bottom:2em}</style>',1)
    # Append references without changing the preserved measurement-method bibliography.
    html=html.replace('</ol></section>', '<li id="ref8"><a href="https://huggingface.co/mwalmsley/zoobot-finetuned-decals10">Walmsley: zoobot-finetuned-decals10 checkpoint</a>; class order/names from the saved local producer metadata.</li><li id="ref9"><a href="https://zoobot.readthedocs.io/en/latest/">Zoobot documentation</a>: classifier, representations and training conventions.</li></ol></section>',1) if '<li id="ref7"' in html else html
    # Ensure external reference links resolve even if base bibliography markup changes.
    if 'id="ref8"' not in html:html=html.replace('</main>','<section><h2>Classification references</h2><p id="ref8"><a href="https://huggingface.co/mwalmsley/zoobot-finetuned-decals10">[8] DECaLS10 checkpoint</a></p><p id="ref9"><a href="https://zoobot.readthedocs.io/en/latest/">[9] Zoobot documentation</a></p></section></main>')
    (root/'report_v3.html').write_text(html)
    manifest=dict(report_version='scientific_report_v3',built_utc=datetime.now(timezone.utc).isoformat(),
        previous_report_sha256=sha256_file(root/'report_v2.html'),morphology_sha256=sha256_file(root/'measurements_metrics.parquet'),
        classification_source_manifest=frozen,clustering_manifest=str(a.embedding_manifest.resolve()),
        classification_figures=len(items),embedding_figures=len(embedding_items),
        analysis_protocol=dict(probability_sum_tolerance=1e-4,human_min_answers=2,vote_sensitivity_min_answers=3,
            correlation='Spearman midranks; per-band pairwise finite values',partial_rank_controls=['gf_mu_mean','gf_re_arcsec','gf_q'],
            partial_rank_cohort='gf_unflagged with all controls finite',class_effect='2*U/(n_class*n_other)-1',
            probability_bins=10,minimum_bin_pairs=20,spread_percentiles=[16,84],significance_testing=False,
            position_angle='Axial/image-coordinate quantity excluded from ordinary scalar correlations',ellipticity='Exactly 1-q_positive; rank correlations reverse sign'),
        files_sha256={str(path.relative_to(root)):sha256_file(path) for directory in [data,assets] for path in sorted(directory.rglob('*')) if path.is_file()},
        source_code={str(Path(f).resolve()):sha256_file(f) for f in [__file__,'algorithmic_morphology/classification_data.py','algorithmic_morphology/classification_plots.py','scripts/prepare_classification_analysis.py','scripts/plot_classification_analysis.py']})
    manifest['html_sha256']=sha256_file(root/'report_v3.html');(root/'report_v3_manifest.json').write_text(json.dumps(manifest,indent=2))
    temporary=root/'index.html.tmp';temporary.write_text(html);temporary.replace(root/'index.html')
    print(f'Published {root / "index.html"}; {len(items)} classification figures; {len(embedding_items)} clustering panels',flush=True)


if __name__=='__main__':main()
