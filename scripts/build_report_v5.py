"""Build report revision 5: the same frozen results, reorganised by C, A and S.

Presentation-only revision. Reads the immutable morphology table and the hashed
v2/v3 analysis tables; never re-measures, refits or redefines any quantity.
Every figure is a single-axes PNG; multi-panel layouts are HTML grids.
"""
import argparse
import base64
from datetime import datetime,timezone
from html import escape
from io import BytesIO
import json
import logging
from multiprocessing import Pool
from pathlib import Path
import re
import sys
import time

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from algorithmic_morphology.report import sha256_file
from algorithmic_morphology.report_data import BANDS,robust_limits
from algorithmic_morphology.paper_plots import LABELS,L,K
from algorithmic_morphology.classification_data import HUMAN_LABELS,CLASS_LABELS
from algorithmic_morphology.report_v5_plots import (Panels,embedding_panel,binned_trends,FAMILIES,PLAIN,PROPERTIES,TARGETS,COHORTS,
    HUMAN_SELECTIONS,EFFECT_SELECTIONS,BAND_COLORS)

LOG=logging.getLogger('report_v5')
VERSION='scientific_report_v5'
KEY_TARGETS=['p_visible_nucleus','p_etg','class_0']
EMBED_METRICS=[m for f in FAMILIES.values() for m in f['metrics']]
FAMILY_STYLE={'C':('#2f5aa8','C','Concentration'),'A':('#b5591a','A','Asymmetry'),'S':('#1d7d74','S','Clumpiness'),
              'D':('#5b6470','G·M20·q','Supplementary diagnostics'),'X':('#4b4f8a','All','All metrics together'),'P':('#6b7280','Rₚ','Shared aperture')}


# ---------------------------------------------------------------- references
REFS={
 'petrosian1976':'Petrosian, V. (1976), <em>Surface brightness and evolution of galaxies</em>, ApJ 209, L1. <a href="https://ui.adsabs.harvard.edu/abs/1976ApJ...209L...1P/abstract">ADS</a>',
 'graham2005p':'Graham, A. W., Driver, S. P., Petrosian, V., et al. (2005), <em>Total galaxy magnitudes and effective radii from Petrosian magnitudes and radii</em>, AJ 130, 1535. <a href="https://ui.adsabs.harvard.edu/abs/2005AJ....130.1535G/abstract">ADS</a>',
 'binney1998':'Binney, J. &amp; Merrifield, M. (1998), <em>Galactic Astronomy</em>, Princeton University Press (book; chapters on galaxy photometry, surface brightness and light profiles).',
 'kent1985':'Kent, S. M. (1985), <em>CCD surface photometry of field galaxies. II. Bulge/disk decompositions</em>, ApJS 59, 115. <a href="https://ui.adsabs.harvard.edu/abs/1985ApJS...59..115K/abstract">ADS</a>',
 'bershady2000':'Bershady, M. A., Jangren, A. &amp; Conselice, C. J. (2000), <em>Structural and photometric classification of galaxies. I. Calibration based on a nearby galaxy sample</em>, AJ 119, 2645. <a href="https://ui.adsabs.harvard.edu/abs/2000AJ....119.2645B/abstract">ADS</a>',
 'abraham1996':'Abraham, R. G., Tanvir, N. R., Santiago, B. X., et al. (1996), <em>Galaxy morphology to I=25 mag in the Hubble Deep Field</em>, MNRAS 279, L47. <a href="https://ui.adsabs.harvard.edu/abs/1996MNRAS.279L..47A/abstract">ADS</a>',
 'conselice2000':'Conselice, C. J., Bershady, M. A. &amp; Jangren, A. (2000), <em>The asymmetry of galaxies: physical morphology for nearby and high-redshift galaxies</em>, ApJ 529, 886. <a href="https://ui.adsabs.harvard.edu/abs/2000ApJ...529..886C/abstract">ADS</a>',
 'conselice2003':'Conselice, C. J. (2003), <em>The relationship between stellar light distributions of galaxies and their formation histories</em>, ApJS 147, 1. <a href="https://arxiv.org/abs/astro-ph/0303065">arXiv</a>',
 'lotz2004':'Lotz, J. M., Primack, J. &amp; Madau, P. (2004), <em>A new nonparametric approach to galaxy morphological classification</em>, AJ 128, 163. <a href="https://arxiv.org/abs/astro-ph/0311352">arXiv</a>',
 'abraham2003':'Abraham, R. G., van den Bergh, S. &amp; Nair, P. (2003), <em>A new approach to galaxy morphology. I. Analysis of the Sloan Digital Sky Survey Early Data Release</em>, ApJ 588, 218. <a href="https://ui.adsabs.harvard.edu/abs/2003ApJ...588..218A/abstract">ADS</a>',
 'rodriguezgomez2019':'Rodriguez-Gomez, V., Snyder, G. F., Lotz, J. M., et al. (2019), <em>The optical morphologies of galaxies in the IllustrisTNG simulation</em>, MNRAS 483, 4140 (statmorph). <a href="https://doi.org/10.1093/mnras/sty3341">Paper</a>; <a href="https://statmorph.readthedocs.io/en/latest/">statmorph documentation</a>',
 'hambleton2011':'Hambleton, K. M., Gibson, B. K., Brook, C. B., et al. (2011), <em>Advanced morphological galaxy classification: a comparison of observed and simulated galaxies</em>, MNRAS 418, 801. <a href="https://arxiv.org/abs/1107.6045">arXiv</a>',
 'thorp2021':'Thorp, M. D., Bottrell, C., Ellison, S. L., et al. (2021), <em>Towards robust determination of non-parametric morphologies in marginal astronomical data</em>, MNRAS 507, 886. <a href="https://arxiv.org/abs/2107.13578">arXiv</a>',
 'sersic1963':'Sérsic, J. L. (1963), <em>Influence of the atmospheric and instrumental dispersion on the brightness distribution in a galaxy</em>, Boletín de la Asociación Argentina de Astronomía 6, 41. <a href="https://ui.adsabs.harvard.edu/abs/1963BAAA....6...41S/abstract">ADS</a>',
 'graham2005':'Graham, A. W. &amp; Driver, S. P. (2005), <em>A concise reference to (projected) Sérsic R^1/n quantities</em>, PASA 22, 118. <a href="https://ui.adsabs.harvard.edu/abs/2005PASA...22..118G/abstract">ADS</a>',
 'peng2002':'Peng, C. Y., Ho, L. C., Impey, C. D. &amp; Rix, H.-W. (2002), <em>Detailed structural decomposition of galaxy images</em>, AJ 124, 266. <a href="https://arxiv.org/abs/astro-ph/0204182">arXiv</a>',
 'peng2010':'Peng, C. Y., Ho, L. C., Impey, C. D. &amp; Rix, H.-W. (2010), <em>Detailed decomposition of galaxy images. II. Beyond axisymmetric models</em>, AJ 139, 2097. <a href="https://arxiv.org/abs/0912.0731">arXiv</a>',
 'haussler2013':'Häußler, B., Bamford, S. P., Vika, M., et al. (2013), <em>MegaMorph – multiwavelength measurement of galaxy structure: complete Sérsic profile information from modern surveys</em>, MNRAS 430, 330 (GALFITM). <a href="https://ui.adsabs.harvard.edu/abs/2013MNRAS.430..330H/abstract">ADS</a>',
 'abbott2018':'DES Collaboration, Abbott, T. M. C., et al. (2018), <em>The Dark Energy Survey Data Release 1</em>, ApJS 239, 18. <a href="https://ui.adsabs.harvard.edu/abs/2018ApJS..239...18A/abstract">ADS</a>',
 'tanoglidis2021':'Tanoglidis, D., Drlica-Wagner, A., Wei, K., et al. (2021), <em>Shadows in the dark: low-surface-brightness galaxies discovered in the Dark Energy Survey</em>, ApJS 252, 18. <a href="https://ui.adsabs.harvard.edu/abs/2021ApJS..252...18T/abstract">ADS</a>',
 'spearman1904':'Spearman, C. (1904), <em>The proof and measurement of association between two things</em>, American Journal of Psychology 15, 72. <a href="https://doi.org/10.2307/1412159">doi</a>',
 'cureton1956':'Cureton, E. E. (1956), <em>Rank-biserial correlation</em>, Psychometrika 21, 287. <a href="https://doi.org/10.1007/BF02289138">doi</a>',
 'kerby2014':'Kerby, D. S. (2014), <em>The simple difference formula: an approach to teaching nonparametric correlation</em>, Comprehensive Psychology 3, 11.IT.3.1. <a href="https://doi.org/10.2466/11.IT.3.1">doi</a>',
 'mannwhitney1947':'Mann, H. B. &amp; Whitney, D. R. (1947), <em>On a test of whether one of two random variables is stochastically larger than the other</em>, Annals of Mathematical Statistics 18, 50. <a href="https://doi.org/10.1214/aoms/1177730491">doi</a>',
 'jolliffe2016':'Jolliffe, I. T. &amp; Cadima, J. (2016), <em>Principal component analysis: a review and recent developments</em>, Phil. Trans. R. Soc. A 374, 20150202. <a href="https://doi.org/10.1098/rsta.2015.0202">doi</a>',
 'mcinnes2018':'McInnes, L., Healy, J. &amp; Melville, J. (2018), <em>UMAP: Uniform Manifold Approximation and Projection for dimension reduction</em>. <a href="https://arxiv.org/abs/1802.03426">arXiv:1802.03426</a>',
 'walmsley2022':'Walmsley, M., Lintott, C., Géron, T., et al. (2022), <em>Galaxy Zoo DECaLS: detailed visual morphology measurements from volunteers and deep learning for 314 000 galaxies</em>, MNRAS 509, 3966. <a href="https://ui.adsabs.harvard.edu/abs/2022MNRAS.509.3966W/abstract">ADS</a>',
 'walmsley2023':'Walmsley, M., Allen, C., Aussel, B., et al. (2023), <em>Zoobot: adaptable deep learning models for galaxy morphology</em>, Journal of Open Source Software 8, 5312. <a href="https://doi.org/10.21105/joss.05312">doi</a>; <a href="https://zoobot.readthedocs.io/en/latest/">documentation</a>',
 'decals10':'Walmsley, M., <em>zoobot-finetuned-decals10</em> model checkpoint. <a href="https://huggingface.co/mwalmsley/zoobot-finetuned-decals10">Hugging Face</a>; class order and names from the saved local producer metadata.',
}


def cite(*keys):
    return '['+', '.join(f'[[cite:{k}]]' for k in keys)+']'


# ---------------------------------------------------------------- small helpers
def table_html(frame,digits=4):
    return '<div class="table-scroll">'+frame.to_html(index=False,float_format=lambda x:f'{x:.{digits}g}',na_rep='—',border=0,classes='data-table')+'</div>'


def equation(text):
    fig=plt.figure(figsize=(11,.8));fig.text(.01,.45,text,fontsize=17,va='center');buffer=BytesIO()
    fig.savefig(buffer,format='svg',bbox_inches='tight',transparent=True);plt.close(fig)
    return '<img class="equation" alt="'+escape(text)+'" src="data:image/svg+xml;base64,'+base64.b64encode(buffer.getvalue()).decode()+'">'


def badge(fam):
    color,letter,name=FAMILY_STYLE[fam]
    return f'<span class="badge" style="--fam:{color}"><b>{escape(letter)}</b> {escape(name)}</span>'


def fig(n,panel=None):return f'[[fig:{n}:{panel}]]' if panel else f'[[fig:{n}]]'


class Registry:
    """Figure specifications; numbers are assigned by order of appearance in the final HTML."""

    def __init__(self,root,panels):
        self.root=root;self.panels=panels;self.items={}

    def add(self,key,family,title,why,look,data,panels,cols=3,extra=None):
        assert key not in self.items,key
        self.items[key]=dict(key=key,family=family,title=title,why=why,look=look,data=data,panels=panels,cols=cols,extra=extra or [])
        return f'[[figure:{key}]]'

    def render(self,key,number):
        f=self.items[key]
        def grid(panels,offset,cols):
            out=[]
            for j,(path,cap) in enumerate(panels):
                letter=chr(97+offset+j);rel=Path(path).relative_to(self.root).as_posix()
                out.append(f'<div class="panel" id="fig-{key}-{letter}"><a href="{rel}" target="_blank" rel="noopener"><img loading="lazy" src="{rel}" alt="Figure {number}{letter}: {escape(cap)}"></a><p class="panel-cap"><b>({letter})</b> {escape(cap)}</p></div>')
            return f'<div class="panel-grid cols-{cols}">'+''.join(out)+'</div>'
        body=grid(f['panels'],0,f['cols']);offset=len(f['panels'])
        for summary,panels,cols in f['extra']:
            body+=f'<details class="more-panels"><summary>{escape(summary)} ({len(panels)} panels)</summary>{grid(panels,offset,cols)}</details>';offset+=len(panels)
        names={Path(p).stem for p,_ in f['panels']}|{Path(p).stem for _,ps,_ in f['extra'] for p,_ in ps}
        rows=[r for r in self.panels.records if r['panel'] in names];audit=''
        if rows:
            t=pd.DataFrame(rows);t['series']=t.series.map(lambda s:PLAIN.get(s.split(' | ')[0].removesuffix('/10'),s.split(' | ')[0])+(' / 10' if s.split(' | ')[0].endswith('/10') else '')+(' | '+s.split(' | ')[1] if ' | ' in s else ''))
            audit='<details class="audit"><summary>Display windows and off-screen counts for every panel</summary>'+table_html(t,5)+'</details>'
        return (f'<figure class="fig" id="fig-{key}" style="--fam:{FAMILY_STYLE[f["family"]][0]}"><div class="fig-head">{badge(f["family"])}<h4>Figure {number}. {escape(f["title"])}</h4></div>{body}'
                f'<figcaption><p><b>Why this figure.</b> {f["why"]}</p><p><b>What to look for.</b> {f["look"]}</p><p><b>Data used.</b> {f["data"]}</p></figcaption>{audit}</figure>')


def resolve(html,registry):
    order=[]
    for key in re.findall(r'\[\[figure:([\w-]+)\]\]',html):
        if key in order:raise ValueError('Figure placed twice: '+key)
        order.append(key)
    missing=set(registry.items)-set(order)
    if missing:raise ValueError(f'Registered but not placed: {missing}')
    number={k:i+1 for i,k in enumerate(order)}
    for key in order:html=html.replace(f'[[figure:{key}]]',registry.render(key,number[key]))
    def figref(m):
        key,panel=m.group(1),m.group(2)
        if key not in number:raise ValueError('Unknown figure reference '+key)
        target=f'fig-{key}-{panel}' if panel else f'fig-{key}'
        return f'<a class="figref" href="#{target}">Figure {number[key]}{panel or ""}</a>'
    html=re.sub(r'\[\[fig:([\w-]+)(?::([a-z]))?\]\]',figref,html)
    cites=[]
    for k in re.findall(r'\[\[cite:(\w+)\]\]',html):
        if k not in cites:cites.append(k)
    for k in cites:
        if k not in REFS:raise ValueError('Unknown reference '+k)
        html=html.replace(f'[[cite:{k}]]',f'<a href="#ref-{k}">{cites.index(k)+1}</a>')
    refs='<ol class="refs">'+''.join(f'<li id="ref-{k}">{REFS[k]}</li>' for k in cites)+'</ol>'
    return html.replace('[[references]]',refs),order,cites


# ---------------------------------------------------------------- main
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,default=Path('reports/full_catalog_20261020'))
    p.add_argument('--workers',type=int,default=10);a=p.parse_args()
    logging.basicConfig(level=logging.INFO,stream=sys.stdout,format='%(asctime)sZ %(levelname)s %(message)s',datefmt='%Y-%m-%dT%H:%M:%S');logging.Formatter.converter=time.gmtime
    root=a.report_dir.resolve();assets=root/'paper_v5_assets';data=root/'paper_v5_data'
    if (root/'report_v5.html').exists():raise FileExistsError('report_v5.html exists; preserve it and choose a new revision')
    assets.mkdir(exist_ok=False);data.mkdir(exist_ok=False)

    # ---- inputs, all hash-checked against earlier manifests
    LOG.info('Verifying frozen inputs')
    prov=json.loads((root/'table_provenance.json').read_text());m2=json.loads((root/'report_v2_manifest.json').read_text());m3=json.loads((root/'report_v3_manifest.json').read_text())
    morph_path=root/'measurements_metrics.parquet'
    if sha256_file(morph_path)!=prov['compact_table_sha256']:raise ValueError('Morphology table changed')
    inputs={'measurements_metrics.parquet':prov['compact_table_sha256']}
    for name,digest in {**m2['files_sha256'],**m3['files_sha256']}.items():
        if name.startswith(('paper_v2_data/','paper_v3_data/')) and name.endswith(('.csv','.parquet','.json')):
            if sha256_file(root/name)!=digest:raise ValueError('Input changed: '+name)
            inputs[name]=digest
    frame=pd.read_parquet(morph_path);frame['coadd_object_id']=frame.coadd_object_id.astype(str)
    assert len(frame)==252600 and frame.coadd_object_id.nunique()==84200
    gf=pd.read_parquet(root/'paper_v2_data/galfit_join.parquet')
    joined=frame.merge(gf,on=['coadd_object_id','band'],how='left',validate='one_to_one',sort=False)
    for c in ['gf_present','gf_unflagged','gf_interior']:joined[c]=joined[c].fillna(False).astype(bool)
    labels=pd.read_parquet(root/'paper_v3_data/classification_labels.parquet')
    classed=frame.merge(labels,on='coadd_object_id',how='left',validate='many_to_one',sort=False);assert len(classed)==len(frame)
    corr=pd.read_csv(root/'paper_v2_data/galfit_correlations.csv');gtrends=pd.read_csv(root/'paper_v2_data/galfit_binned_trends.csv')
    assoc=pd.read_csv(root/'paper_v3_data/associations.csv');effects=pd.read_csv(root/'paper_v3_data/class_effects.csv');ptrends=pd.read_csv(root/'paper_v3_data/probability_trends.csv')
    contract=json.loads((root/'paper_v3_data/feature_contract.json').read_text());lim={k:tuple(v) for k,v in contract['limits'].items()}
    emb=json.loads((root/'paper_v3_data/clustering_overlay_manifest.json').read_text())['embeddings']
    stats=pd.read_csv(root/'statistics.csv');bandcmp=pd.read_csv(root/'band_comparisons.csv');paircmp=pd.read_csv(root/'paired_comparisons.csv')
    gcov=pd.read_csv(root/'paper_v2_data/galfit_coverage.csv');ccov=pd.read_csv(root/'paper_v3_data/coverage.csv');summary=json.loads((root/'summary.json').read_text())
    fitsens=pd.read_csv(root/'paper_v2_data/fit_configuration_sensitivity.csv')
    # v2 stored no binned q+ trends. Recompute them with the v2 algorithm, after proving it reproduces stored v2 rows.
    for m in ['C','A_corrected']:
        key=['property','band','left'];pd.testing.assert_frame_equal(binned_trends(joined,m).sort_values(key).reset_index(drop=True),gtrends[gtrends.metric==m].sort_values(key).reset_index(drop=True),check_exact=False,rtol=1e-12)
    qtrend=binned_trends(joined,'axis_ratio_positive_aperture_v1');qtrend.to_csv(data/'galfit_binned_trends_axis_ratio.csv',index=False)
    gtrends=pd.concat([gtrends,qtrend],ignore_index=True)

    P=Panels(assets);R=Registry(root,P)
    def st(metric,band,col):return float(stats[(stats.metric==metric)&(stats.band==band)][col].iloc[0])
    def neg(metric,band):return 100*st(metric,band,'negative')/st(metric,band,'valid')
    def gr(metric,prop,band,cohort='gf_unflagged'):return float(corr[(corr.cohort==cohort)&(corr.metric==metric)&(corr.property==prop)&(corr.band==band)].spearman.iloc[0])
    def ar(metric,target,band='r',cohort='available',stat='spearman'):return float(assoc[(assoc.band==band)&(assoc.metric==metric)&(assoc.target==target)&(assoc.cohort==cohort)][stat].iloc[0])
    def ef(metric,sel,band='r'):return float(effects[(effects.band==band)&(effects.metric==metric)&(effects.selection==sel)&(effects.confidence=='all')].rank_biserial.iloc[0])
    def bc(metric,bands,col='median_difference'):return float(bandcmp[(bandcmp.metric==metric)&(bandcmp.bands==bands)][col].iloc[0])
    def strongest(metric,band='r',k=2):
        s=assoc[(assoc.band==band)&(assoc.metric==metric)&(assoc.cohort=='available')].copy();s['abs']=s.spearman.abs()
        return ', '.join(f'{TARGETS[r.target]} (ρ={r.spearman:+.2f})' for r in s.nlargest(k,'abs').itertuples())
    def trio(fn):return '/'.join(f'{fn(b):.2f}' for b in BANDS)
    def rng(metric,band):return f'{st(metric,band,"median"):.3g} (16th–84th percentile {st(metric,band,"p16"):.3g} to {st(metric,band,"p84"):.3g})'
    h2=int(ccov.set_index('band').loc['r','human_two_plus']);zv=int(ccov.set_index('band').loc['r','zoobot_valid'])
    gmain={r.band:int(r.unflagged) for r in gcov.itertuples()}

    # ---- shared data-description sentences
    D_ALL=f'All 84,200 catalog galaxies per band; each curve uses only that band’s finite values (counts in the legend). Missing cutouts and flagged measurements are excluded, not set to zero.'
    D_GF=(f'The main “unflagged” GALFIT cohort (successful single-Sérsic fit, finite physical parameters, no parameter warning): {gmain["g"]:,}/{gmain["r"]:,}/{gmain["i"]:,} galaxies in g/r/i. '
          'Lines are medians in 16 common x bins, shading the 16th–84th percentile spread of galaxies (not an error bar); bins with fewer than 30 pairs are blank. Numbers come unchanged from <code>paper_v2_data/galfit_binned_trends.csv</code>.')
    D_HH=('Galaxies with at least two human morphology answers, split by the published seniority-weighted consensus (ETG, LTG, Featureless, Not sure) or no consensus (“split”). '
          'Each curve is normalised to its own class size, so shapes compare even though class sizes differ. Display window: pooled g/r/i 0.5–99.5 percentiles with 3.5% padding.')
    D_AS=(f'All rank correlations from <code>paper_v3_data/associations.csv</code>. Human values are raw answer fractions (only galaxies with ≥2 answers to that question; {h2:,} for morphology); DECaLS10 values are model probabilities ({zv:,} valid vectors). '
          'Four cohorts: all available pairs; galaxies with both ≥2 human votes and DECaLS10; ≥3 human votes; and the GALFIT-unflagged cohort after removing the rank-linear part explained by mean surface brightness, effective radius and axis ratio (partial rank correlation).')
    D_EF=('Rank-biserial separations from <code>paper_v3_data/class_effects.csv</code>. Human rows use published consensus memberships among galaxies with ≥2 comparable answers; DECaLS10 rows use the class with the largest model probability. '
          'The right-hand block repeats DECaLS10 rows only for galaxies whose largest probability is at least 0.5 (human rows are blank there by design).')
    D_PT=('Galaxies with a valid human fraction (≥2 answers) or DECaLS10 probability, binned in steps of 0.1. Points are bin medians, shading the 16th–84th percentile spread; bins with fewer than 20 galaxies are omitted. '
          'Numbers come unchanged from <code>paper_v3_data/probability_trends.csv</code>.')
    default_space=next(i for i,e in enumerate(emb) if e['title']=='GALFIT g+i fitted properties' and e['projection']=='UMAP (120 neighbors)' and e['band']=='r')
    D_EM=(f'The {emb[default_space]["rows"]:,} galaxies of the saved “loose selection” used by the clustering project (a subset of the catalog), placed at their saved coordinates in the <em>GALFIT g+i fitted properties, UMAP (120 neighbours)</em> map. '
          'Each hexagonal cell is coloured by the median value of its galaxies (at least 3 finite values); grey cells contain galaxies but too few finite values. Colour limits are pooled g/r/i robust limits; Conselice-like limits are ten times the Lotz-like limits.')

    # ---- embedding maps (one PNG per metric per saved map), rendered in parallel
    LOG.info('Rendering %d single-metric embedding maps',len(emb)*len(EMBED_METRICS))
    tasks=[];emb_index=[]
    for i,e in enumerate(emb):
        if sha256_file(e['export'])!=e['export_sha256']:raise ValueError('Embedding export changed: '+e['export'])
        d=pd.read_csv(e['export'],usecols=['embedding_x','embedding_y',*EMBED_METRICS])
        if len(d)!=e['rows']:raise ValueError('Embedding row count changed')
        counts={}
        for m in EMBED_METRICS:
            v=pd.to_numeric(d[m],errors='coerce').to_numpy(float);counts[m]=int(np.isfinite(v).sum())
            tasks.append((assets/f'emb_{i:02d}_{m}.png',(d.embedding_x.to_numpy(float),d.embedding_y.to_numpy(float)),v,lim[m],PLAIN[m],f'{PLAIN[m]} · {e["band"]} band',110))
        emb_index.append(dict(index=i,title=e['title'],projection=e['projection'],band=e['band'],rows=e['rows'],matched=e['matched'],finite=counts,export_sha256=e['export_sha256']))
    with Pool(a.workers) as pool:
        for j,_ in enumerate(pool.imap_unordered(embedding_panel,tasks,chunksize=8)):
            if (j+1)%100==0:LOG.info('Embedding maps: %d/%d',j+1,len(tasks))
    def emb_panels(metric,bands=('r',)):
        out=[]
        for b in bands:
            i=next(k for k,e in enumerate(emb) if e['title']==emb[default_space]['title'] and e['projection']==emb[default_space]['projection'] and e['band']==b)
            out.append((assets/f'emb_{i:02d}_{metric}.png',f'{PLAIN[metric]}, {b} band ({emb_index[i]["finite"][metric]:,} finite of {emb[i]["rows"]:,})'))
        return out

    LOG.info('Rendering single-panel figures')
    def hist_set(metric,title_prefix):
        zoom=P.band_histogram(frame,metric,robust_limits([frame[metric]]),f'{metric}_dist',f'{title_prefix}: all three bands')
        full=P.band_histogram(frame,metric,robust_limits([frame[metric]],quantiles=(0,1)),f'{metric}_dist_full',f'{title_prefix}: full finite range')
        return zoom,full
    def galfit_panels(metric):
        yl=robust_limits([joined[metric]])
        return [(P.galfit_trend(gtrends,prop,metric,yl,f'{metric}__{prop}',f'{PLAIN[metric]} vs {PLAIN[prop]}'),f'{PLAIN[metric]} against {PLAIN[prop]}') for prop in PROPERTIES]
    def human_hist(metric,band):
        return (P.human_class_histogram(classed,metric,band,lim[metric],f'{metric}__human_{band}',f'{PLAIN[metric]} by human consensus · {band} band'),f'{PLAIN[metric]}, {band} band')
    def assoc_panel(metric):
        return (P.association_heatmap(assoc,metric,f'{metric}__assoc',f'{PLAIN[metric]}: rank correlation with human answers and DECaLS10 probabilities'),f'{PLAIN[metric]}: all targets, bands and cohorts')
    def effect_panel(metric):
        return (P.effect_heatmap(effects,metric,f'{metric}__effects',f'{PLAIN[metric]}: separation of each class from the rest'),f'{PLAIN[metric]}: class separation')
    def trend(metric,target):
        return (P.probability_trend(ptrends,target,metric,lim[metric],f'{metric}__trend_{target}',f'{PLAIN[metric]} vs {TARGETS[target]}'),f'{PLAIN[metric]} against {TARGETS[target]}')
    def all_trends(metrics):
        return [trend(m,t) for t in TARGETS if t not in KEY_TARGETS for m in metrics]

    # ======================================================== shared aperture section
    rp_zoom=P.band_histogram(frame,'rpetro',robust_limits([frame.rpetro]),'rpetro_dist','Petrosian radius: all three bands')
    unfl=joined[joined.gf_unflagged]
    rp_vs=[P.density(unfl,b,'gf_re_pixels','rpetro',robust_limits([unfl.gf_re_pixels]),robust_limits([unfl.rpetro]),f'rpetro_vs_re_{b}',f'{b} band',equality=True)[0] for b in BANDS]
    sec_aperture=f'''<section id="aperture"><h2>3. The shared measurement aperture {badge('P')}</h2>
<p>Every index in this report (C, A, S and the supplementary diagnostics) is measured inside the same circle: 1.5 times the galaxy’s <a href="#g-petrosian">Petrosian radius</a> R<sub>P</sub>, centred on the asymmetry-minimising centre. The aperture therefore sets the “ruler” for all three CAS indices. {fig('rpetro_dist')} shows how large that circle is, and {fig('rpetro_vs_re')} shows that the image-based Petrosian radius and the model-based <a href="#g-reff">GALFIT effective radius</a> are related but not interchangeable. This is why the analysis never substitutes one radius for the other.</p>
{R.add('rpetro_dist','P','Petrosian radius distribution in g, r and i',
 'R<sub>P</sub> is the size scale on which C, A and S are measured. Knowing its typical value (in 0.263″ DES pixels) tells you how many pixels each measurement averages over, and therefore how much noise and blurring by the <a href="#g-psf">PSF</a> can influence it.',
 'The median and spread per band; whether bands agree (the same galaxy should have a similar R<sub>P</sub> in all bands) and how many galaxies have a very small R<sub>P</sub> (few pixels, so measurements are more noise- and PSF-limited).',
 D_ALL+' The full finite range is in panel (b).',[(rp_zoom,'Petrosian radius R_P in pixels, robust display window')],cols=1,
 extra=[('Full finite range',[(P.band_histogram(frame,'rpetro',robust_limits([frame.rpetro],quantiles=(0,1)),'rpetro_dist_full','Petrosian radius: full finite range'),'R_P, full finite range')],1)])}
{R.add('rpetro_vs_re','P','Petrosian radius versus GALFIT effective radius, per band',
 'A reader might expect the two “sizes” to be the same thing. They are not: R<sub>P</sub> is defined from the observed light profile, R<sub>e</sub> from a fitted model that encloses half of the model light. For a pure exponential (Sérsic n=1) profile, R<sub>P</sub> at η=0.2 is roughly twice R<sub>e</sub> '+cite('graham2005p')+'.',
 'Whether points follow a tight line (a fixed ratio) or scatter widely, and the dashed equality line, which shows that R<sub>P</sub> is systematically larger than R<sub>e</sub>. The Spearman ρ in each title measures how well the ranking of galaxies by one radius agrees with the other. Note also the sparse population to the right of the main ridge, with large fitted R<sub>e</sub> but ordinary R<sub>P</sub>; these fits deserve visual inspection before their R<sub>e</sub> is trusted.',
 'Main unflagged GALFIT cohort; both radii in pixels. Each panel is one band.',[(p,f'{b} band') for p,b in zip(rp_vs,BANDS)],cols=3)}
</section>'''

    # ======================================================== C
    c_zoom,c_full=hist_set('C','Concentration C')
    sec_c=f'''<section id="concentration" class="family" style="--fam:{FAMILY_STYLE['C'][0]}"><h2>4. Concentration (C) {badge('C')}</h2>
<p class="lead"><a href="#g-concentration">Concentration</a> asks: <em>how strongly is the light piled up toward the centre?</em> C=5 log<sub>10</sub>(r<sub>80</sub>/r<sub>20</sub>), the ratio of the radii that enclose 80% and 20% of the light inside the aperture {cite('kent1985','bershady2000','conselice2003')}. A flat, diffuse disc has C≈2.5–3; a bright nucleus or bulge pushes C upward. This section shows the C distribution ({fig('C_dist')}), how C relates to fitted galaxy structure ({fig('C_galfit')}), and how it compares with human and machine classifications ({fig('C_human')}–{fig('C_embed')}).</p>
<h3 id="c-dist">4.1 How concentrated are these galaxies?</h3>
<p>The median C is {rng('C','g')} in g, {rng('C','r')} in r and {rng('C','i')} in i ({fig('C_dist')}). These values sit at the diffuse end of the range found for bright nearby galaxies {cite('bershady2000','conselice2003')}, as expected for a catalog dominated by <a href="#g-sb">low-surface-brightness</a> candidates. For the same galaxies, C is slightly lower in g than in r and i: the matched median difference is {bc('C','g-r'):.3f} (g−r) and {bc('C','g-i'):.3f} (g−i). Redder light being more centrally concentrated is one possible explanation; band-dependent noise, PSF size and aperture also contribute and are not separated here.</p>
{R.add('C_dist','C','Distribution of concentration C in g, r and i',
 'This is the basic picture of the C values in the catalog and the reference against which all later C plots should be read.',
 'Where the peak lies, how wide the distribution is, and whether the three band curves are shifted relative to each other. The curves show the percentage of each band’s finite sample per bin; tails outside the window are listed in the audit table and shown in panel (b).',
 D_ALL,[(c_zoom,'Concentration C, robust display window')],cols=1,extra=[('Full finite range including outliers',[(c_full,'Concentration C, full finite range')],1)])}
<h3 id="c-galfit">4.2 Concentration and fitted structure (GALFIT)</h3>
<p>An independent description of each galaxy’s light profile comes from the <a href="#g-galfit">GALFIT</a> single-<a href="#g-sersic">Sérsic</a> fit. The most direct check of C is the Sérsic index n, which also measures how peaked a profile is. C rises with n in every band (Spearman ρ = {trio(lambda b:gr('C','gf_n',b))} for g/r/i; {fig('C_galfit','a')}), so the non-parametric C and the model-based n tell a consistent story. The relations with effective radius (ρ = {trio(lambda b:gr('C','gf_re_arcsec',b))}; {fig('C_galfit','c')}) and mean surface brightness (ρ = {trio(lambda b:gr('C','gf_mu_mean',b))}; {fig('C_galfit','b')}) are weaker and change with band. C is therefore much less sensitive to how faint a galaxy is than raw asymmetry ({fig('A_galfit')}) or Gini ({fig('D_galfit')}), whose surface-brightness correlations reach |ρ|≈0.7–0.9.</p>
{R.add('C_galfit','C','Concentration C against six GALFIT properties',
 'C and the GALFIT parameters describe the same images in different ways. Agreement with the Sérsic index supports C as a measure of light concentration, and trends with brightness or size reveal measurement sensitivities.',
 'Panel (a) is the key comparison: a rising line means galaxies with more peaked fitted profiles also have higher C. In (b) and (d) a flat line means C does not depend on how faint the galaxy is. In (c) a trend with size can reveal the influence of resolution, because small galaxies are more affected by blurring. Panels (e) and (f) are checks on shape and fit quality. Line style also distinguishes the bands: g solid, r dashed, i dotted.',
 D_GF,galfit_panels('C'),cols=2)}
<h3 id="c-classes">4.3 Concentration and human / DECaLS10 classifications</h3>
<p>Human volunteers answered morphology, nucleus and fitting questions; the DECaLS10 <a href="#g-zoobot">Zoobot</a> model assigns ten class probabilities {cite('walmsley2022','walmsley2023','decals10')}. A clear C association is with the human visible-nucleus answer: ρ={ar('C','p_visible_nucleus'):.3f} in r ({fig('C_assoc')}, {fig('C_trends','a')}), and galaxies with a consensus visible nucleus have higher C than the rest (rank-biserial {ef('C','analysis-visible-nucleus-visible-min2'):+.3f} in r; {fig('C_effects')}). By contrast, the human ETG fraction has only ρ={ar('C','p_etg'):.3f} with C. The four human morphology consensus classes overlap strongly in C ({fig('C_human')}), so C alone cannot reproduce these labels. The strongest r-band associations overall are {strongest('C')}.</p>
{R.add('C_human','C','Concentration C for each human morphology consensus class',
 'If human classes such as ETG (elliptical-like) or Featureless correspond to different light concentrations, their C distributions should be shifted relative to each other.',
 'Horizontal shifts between curves (different typical C) and how much they overlap. Strong overlap means C cannot separate these classes by itself. Panel (a) shows r band; g and i are in the expandable block.',
 D_HH,[human_hist('C','r')],cols=1,extra=[('g and i bands',[human_hist('C','g'),human_hist('C','i')],2)])}
{R.add('C_assoc','C','Concentration C: rank correlations with every human answer and DECaLS10 class',
 'This single table-like heatmap summarises every association between C and the 17 external labels, in all three bands and four comparison cohorts, so that a pattern can be checked for consistency at a glance.',
 'Red cells mean higher C goes with a larger answer fraction/probability; blue the opposite; white means no monotonic association. Read along a row: if the colour stays the same across bands and cohorts, the association is robust to band choice and sample selection. The rightmost block shows what remains after controlling for brightness, size and shape.',
 D_AS,[assoc_panel('C')],cols=1)}
{R.add('C_effects','C','Concentration C: how well each class separates from all other galaxies',
 'Correlations use continuous fractions; this figure asks the complementary question for hard class labels: does a randomly chosen class member tend to have a higher C than a randomly chosen non-member?',
 'Values near +1 or −1 mean almost complete separation; near 0 means none. Compare the human nucleus rows with the morphology rows, and the DECaLS10 rows with and without the ≥0.5 confidence cut.',
 D_EF,[effect_panel('C')],cols=1)}
{R.add('C_trends','C','Concentration C along human answer fractions and DECaLS10 probabilities',
 'Heatmaps reduce each relation to one number. These trend lines show the <em>shape</em> of the three most informative relations: visible nucleus, ETG answer and DECaLS10 Disturbed probability.',
 'Whether the median C climbs steadily, saturates or does nothing as the fraction or probability rises, and how wide the spread remains. Human fractions take only a few discrete values (e.g. 0, ½, 1 for two votes). The other 14 targets are in the expandable block.',
 D_PT,[trend('C',t) for t in KEY_TARGETS],cols=3,extra=[('All other human answers and DECaLS10 classes',all_trends(['C']),3)])}
{R.add('C_embed','C','Concentration C across a saved clustering map',
 'The clustering project placed galaxies on 2-D maps (here a <a href="#g-umap">UMAP</a> of GALFIT g+i properties). Colouring the map by C shows whether concentrated galaxies occupy a particular region, i.e. whether the clustering features already carry concentration information.',
 'Smooth colour gradients (C varies systematically across the map) versus a speckled pattern (unrelated). Axes have no physical units. Other maps and bands are in <a href="#embedding-atlas">Section 8</a>.',
 D_EM,emb_panels('C'),cols=1,extra=[('g and i bands',emb_panels('C',('g','i')),2)])}
</section>'''

    # ======================================================== A
    a_raw_z,a_raw_f=hist_set('A_raw','Raw asymmetry A_raw');a_cor_z,a_cor_f=hist_set('A_corrected','Sky-corrected asymmetry A_corr');a_sky_z,a_sky_f=hist_set('A_sky','Sky asymmetry term A_sky')
    alim=robust_limits([frame.A_raw,frame.A_corrected])
    a_over=[P.pair_histogram(frame,b,['A_raw','A_corrected'],['Raw','Sky-corrected'],alim,f'A_overlay_{b}',f'Raw and corrected asymmetry · {b} band',r'Asymmetry $A$ (dimensionless)') for b in BANDS]
    a_pair=[];a_rho={}
    for b in BANDS:
        path,info=P.density(frame,b,'A_raw','A_corrected',robust_limits([frame.A_raw]),robust_limits([frame.A_corrected]),f'A_paired_{b}',f'{b} band');a_pair.append(path);a_rho[b]=info['rho']
    def gal_pairs(metrics):
        cols=[galfit_panels(m) for m in metrics]
        return [cols[j][i] for i in range(len(PROPERTIES)) for j in range(len(metrics))]
    sec_a=f'''<section id="asymmetry" class="family" style="--fam:{FAMILY_STYLE['A'][0]}"><h2>5. Asymmetry (A) {badge('A')}</h2>
<p class="lead"><a href="#g-asymmetry">Asymmetry</a> asks: <em>does the galaxy look the same after rotating it by 180°?</em> A<sub>raw</sub> is the summed absolute difference between the image and its rotated copy divided by the summed absolute light {cite('abraham1996','conselice2000')}. Image noise is never symmetric, so even a perfectly symmetric galaxy has A<sub>raw</sub>&gt;0. We therefore also estimate the same quantity on a blank sky patch, A<sub>sky</sub>, and keep A<sub>corr</sub>=A<sub>raw</sub>−A<sub>sky</sub> {cite('conselice2003')}. <strong>Both raw and corrected values are always shown side by side</strong>: raw is the direct measurement, corrected is the noise-adjusted estimate, and neither is preferred here.</p>
<h3 id="a-dist">5.1 How large is the asymmetry, and how large is the noise term?</h3>
<p>Median A<sub>raw</sub> is {rng('A_raw','r')} in r ({fig('A_dist','a')}), but the median sky term is almost the same size, {st('A_sky','r','median'):.3g} ({fig('A_dist','c')}). As a result A<sub>corr</sub> is negative for {neg('A_corrected','g'):.0f}%/{neg('A_corrected','r'):.0f}%/{neg('A_corrected','i'):.0f}% of finite g/r/i measurements ({fig('A_dist','b')}). This is the central fact about asymmetry in this faint sample: <em>most of the measured raw asymmetry is comparable to what noise alone produces.</em> Negative corrected values are valid estimates (noise in the sky patch can exceed that in the galaxy) and are kept. {fig('A_overlay')} overlays the two versions per band and {fig('A_paired')} shows them galaxy by galaxy: the rank agreement between raw and corrected is only ρ={a_rho['g']:.2f}/{a_rho['r']:.2f}/{a_rho['i']:.2f} (g/r/i), so the correction reorders galaxies substantially.</p>
{R.add('A_dist','A','Distributions of raw asymmetry, sky-corrected asymmetry and the sky term',
 'Raw and corrected asymmetry can only be interpreted together with the size of the noise correction. These three panels show the three ingredients of A<sub>corr</sub>=A<sub>raw</sub>−A<sub>sky</sub>.',
 'Compare the location of (a) and (c): when they are similar, the corrected values in (b) cluster near zero and many become negative (left of the grey zero line). Also compare the bands: the i band has a larger sky term.',
 D_ALL+' Full finite ranges are in panels (d)–(f).',[(a_raw_z,'Raw asymmetry A_raw'),(a_cor_z,'Sky-corrected asymmetry A_corr'),(a_sky_z,'Sky term A_sky')],cols=3,
 extra=[('Full finite ranges including outliers',[(a_raw_f,'A_raw, full range'),(a_cor_f,'A_corr, full range'),(a_sky_f,'A_sky, full range')],3)])}
{R.add('A_overlay','A','Raw and sky-corrected asymmetry on a common axis, one band per panel',
 'Placing raw and corrected asymmetry on the same axis shows directly how far the noise correction moves the whole distribution.',
 'The horizontal shift between the solid blue (raw) and dashed orange (corrected) curves, and the fraction of the orange curve left of zero. These are distribution comparisons; the galaxy-by-galaxy pairing is shown in '+fig('A_paired')+'.',
 D_ALL+' Each curve uses its own finite sample (counts can differ, because the corrected value also needs a clean sky patch).',[(p,f'{b} band') for p,b in zip(a_over,BANDS)],cols=3)}
{R.add('A_paired','A','Raw versus sky-corrected asymmetry for the same galaxies',
 'If the sky correction were a constant offset, galaxies would lie on a narrow band parallel to the equality line, and their ranking would be unchanged. This figure tests that.',
 'The width of the cloud and the Spearman ρ in each title. A wide cloud and a moderate ρ mean the correction differs strongly from galaxy to galaxy, so conclusions can depend on which version is used.',
 'Galaxies with both values finite in that band. Colour shows the number of galaxies per hexagonal cell on a logarithmic scale.',[(p,f'{b} band') for p,b in zip(a_pair,BANDS)],cols=3)}
<h3 id="a-galfit">5.2 Asymmetry and fitted structure (GALFIT)</h3>
<p>The most important relation is with <a href="#g-sb">mean effective surface brightness</a> ({fig('A_galfit','c')}). A<sub>raw</sub> rises steeply toward fainter galaxies (ρ = {trio(lambda b:gr('A_raw','gf_mu_mean',b))} for g/r/i), exactly what is expected if noise dominates the raw value: fainter galaxies have less signal relative to the same sky noise. After the sky correction the trend reverses (ρ = {trio(lambda b:gr('A_corrected','gf_mu_mean',b))}; {fig('A_galfit','d')}), so the correction more than removes the faint-galaxy rise. This is consistent with substantial noise sensitivity of both versions and does <em>not</em> establish a physical link between asymmetry and surface brightness {cite('thorp2021')}.</p>
{R.add('A_galfit','A','Raw (left) and sky-corrected (right) asymmetry against six GALFIT properties',
 'Shows whether asymmetry tracks galaxy structure (Sérsic index, size) or mainly brightness, which would indicate noise sensitivity. Raw and corrected are placed side by side for every property.',
 'Rows are GALFIT properties; the left column is A<sub>raw</sub>, the right column A<sub>corr</sub>. Pay most attention to the surface-brightness row (c, d) and the magnitude row (g, h): strong opposite trends there indicate that noise, not intrinsic asymmetry, drives much of the signal.',
 D_GF,gal_pairs(['A_raw','A_corrected']),cols=2)}
<h3 id="a-classes">5.3 Asymmetry and human / DECaLS10 classifications</h3>
<p>Asymmetry is the CAS index meant to trace disturbed or merging systems. In r, A<sub>raw</sub> correlates with the DECaLS10 Disturbed probability with ρ={ar('A_raw','class_0'):.3f}, whereas A<sub>corr</sub> gives ρ={ar('A_corrected','class_0'):.3f} ({fig('A_assoc')}). The strongest r-band associations are {strongest('A_raw')} for raw and {strongest('A_corrected')} for corrected asymmetry. The human class distributions ({fig('A_human')}), class separations ({fig('A_effects')}) and trend lines ({fig('A_trends')}) should be read with the noise dominance of {fig('A_dist')} in mind: a class that is on average fainter will look more asymmetric in A<sub>raw</sub> for that reason alone. The GALFIT-adjusted cohort in {fig('A_assoc')} partially removes this effect.</p>
{R.add('A_human','A','Raw and sky-corrected asymmetry for each human morphology consensus class',
 'Tests whether galaxies that humans call irregular/unsure or LTG (spiral) appear more asymmetric than ETG or Featureless galaxies, as would be expected if A traced visible irregularity.',
 'Shifts between class curves in (a) raw and (b) corrected asymmetry. A shift present in raw but absent in corrected suggests a brightness/noise origin. g and i bands are in the expandable block.',
 D_HH,[human_hist('A_raw','r'),human_hist('A_corrected','r')],cols=2,extra=[('g and i bands',[human_hist(m,b) for b in ('g','i') for m in ('A_raw','A_corrected')],2)])}
{R.add('A_assoc','A','Raw and sky-corrected asymmetry: rank correlations with every human answer and DECaLS10 class',
 'One heatmap per asymmetry version, with identical layout, so the effect of the sky correction on every association can be compared cell by cell.',
 'Compare the same cell in (a) raw and (b) corrected. Associations that weaken or change sign after correction, or after the GALFIT adjustment (rightmost block), are likely to be driven by brightness rather than by shape.',
 D_AS,[assoc_panel('A_raw'),assoc_panel('A_corrected')],cols=1)}
{R.add('A_effects','A','Raw and sky-corrected asymmetry: separation of each class from all other galaxies',
 'Shows, for hard class labels, whether members tend to be more (positive) or less (negative) asymmetric than non-members.',
 'Rows with strong colour in (a) but pale colour in (b) indicate a separation that disappears after the noise correction.',
 D_EF,[effect_panel('A_raw'),effect_panel('A_corrected')],cols=1)}
{R.add('A_trends','A','Raw (left) and sky-corrected (right) asymmetry along key human fractions and DECaLS10 probabilities',
 'Shows the shape of the asymmetry relation with visible nucleus, ETG answer and DECaLS10 Disturbed probability for both versions.',
 'Each row is one target; compare left (raw) and right (corrected). The Disturbed row (e, f) is the most directly relevant to the classic use of asymmetry. The remaining 14 targets are in the expandable block.',
 D_PT,[trend(m,t) for t in KEY_TARGETS for m in ('A_raw','A_corrected')],cols=2,extra=[('All other human answers and DECaLS10 classes',all_trends(['A_raw','A_corrected']),2)])}
{R.add('A_embed','A','Raw and sky-corrected asymmetry across a saved clustering map',
 'Shows whether asymmetric galaxies occupy a particular region of the clustering map, and whether raw and corrected asymmetry highlight the same regions.',
 'Compare the colour patterns of (a) and (b). Because A<sub>raw</sub> follows surface brightness, a gradient in (a) may simply reflect a brightness gradient across the map.',
 D_EM,emb_panels('A_raw')+emb_panels('A_corrected'),cols=2,extra=[('g and i bands',[p for b in ('g','i') for m in ('A_raw','A_corrected') for p in emb_panels(m,(b,))],2)])}
</section>'''

    # ======================================================== S
    sv=[L+'_raw',L+'_corrected',K+'_raw',K+'_corrected']
    s_hist=[];s_full=[]
    for m in sv:
        role='raw' if m.endswith('raw') else 'corrected'
        lo,hi=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10]);sc=10 if m.startswith(K) else 1
        s_hist.append((P.band_histogram(frame,m,(lo*sc,hi*sc),f'{m}_dist',PLAIN[m]),PLAIN[m]))
        lo,hi=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10],quantiles=(0,1))
        s_full.append((P.band_histogram(frame,m,(lo*sc,hi*sc),f'{m}_dist_full',PLAIN[m]+': full range'),PLAIN[m]+', full range'))
    s_over=[];s_pair=[];s_rho={}
    for role in ['raw','corrected']:
        lims=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10])
        for b in BANDS:
            s_over.append((P.pair_histogram(frame,b,[L+'_'+role,K+'_'+role],['Lotz-like','Conselice-like ÷ 10'],lims,f'S_overlay_{role}_{b}',f'{role.capitalize()} clumpiness · {b} band',r'Clumpiness $S_{\rm L}$ or $S_{\rm C}/10$',scales=(1,10)),f'{role.capitalize()}, {b} band'))
        for b in BANDS:
            path,info=P.density(frame,b,L+'_'+role,K+'_'+role,lims,lims,f'S_paired_{role}_{b}',f'{role.capitalize()} · {b} band',div_y=10,equality=True);s_pair.append((path,f'{role.capitalize()}, {b} band'));s_rho[(role,b)]=info['rho']
    pc=paircmp.set_index(['band','comparison'])
    sec_s=f'''<section id="clumpiness" class="family" style="--fam:{FAMILY_STYLE['S'][0]}"><h2>6. Clumpiness (S) {badge('S')}</h2>
<p class="lead"><a href="#g-clumpiness">Clumpiness</a> (also called smoothness) asks: <em>how much of the light sits in small-scale knots rather than in a smooth distribution?</em> The image is blurred with a box filter, the blurred version is subtracted, and the positive leftover light is divided by the total light {cite('conselice2003')}. Different authors choose different blurring widths, central exclusion zones and scale factors, so <strong>two named versions are kept throughout</strong>: a <em>Lotz/statmorph-like</em> profile (blur 0.25 R<sub>P</sub>, central 0.25 R<sub>P</sub> excluded, factor 1) {cite('lotz2004','rodriguezgomez2019')} and a <em>Conselice-like</em> profile (blur 0.30 R<sub>P</sub>, central 0.075 R<sub>P</sub> excluded, factor 10) {cite('conselice2003','hambleton2011')}. As for asymmetry, each has a raw and a sky-corrected value: four S columns in total. Both are literature-motivated approximations, not exact reproductions.</p>
<h3 id="s-dist">6.1 Distributions of the four clumpiness outputs</h3>
<p>{fig('S_dist')} shows all four outputs. Raw values are large (median Lotz-like {st(L+'_raw','r','median'):.3g}, Conselice-like {st(K+'_raw','r','median'):.3g} in r) but, as for asymmetry, the sky term is nearly as large, so the corrected medians are small ({st(L+'_corrected','r','median'):.3g} and {st(K+'_corrected','r','median'):.3g}) and {neg(L+'_corrected','r'):.0f}% (Lotz-like) / {neg(K+'_corrected','r'):.0f}% (Conselice-like) of corrected r-band values are negative. Because the Conselice-like profile multiplies by ten, its axes are drawn at exactly ten times the Lotz-like axes. After dividing by ten, the two profiles almost coincide ({fig('S_overlay')}); galaxy by galaxy their rank agreement is ρ={pc.loc[('r','S raw'),'spearman']:.3f} raw and ρ={pc.loc[('r','S corrected'),'spearman']:.3f} corrected in r ({fig('S_paired')}). The different kernel and central exclusion therefore change individual values but hardly the ordering of galaxies. This does not show that either profile is preferable.</p>
{R.add('S_dist','S','Distributions of the four clumpiness outputs (two profiles × raw/corrected)',
 'All four stored S outputs side by side, with the bands overlaid, so that profile and correction effects can be separated visually.',
 'Rows are profiles (Lotz-like top, Conselice-like bottom); columns are raw (left) and sky-corrected (right). Compare left with right to see the size of the noise correction; compare top with bottom (allowing for the ×10 axis) to see the profile difference.',
 D_ALL+' Within each correction state, the Conselice-like display window is exactly ten times the Lotz-like window. Full finite ranges are in panels (e)–(h).',s_hist,cols=2,extra=[('Full finite ranges including outliers',s_full,2)])}
{R.add('S_overlay','S','Lotz-like and Conselice-like ÷ 10 clumpiness on a common axis',
 'Removing only the factor-of-ten normalisation shows how much of the difference between the two profiles is a scale factor and how much comes from their different blur width and central exclusion.',
 'Top row raw, bottom row corrected; one band per column. Near-identical curves mean the remaining definitional differences matter little for the distribution as a whole.',
 D_ALL+' The Conselice-like values are divided by ten <em>for display only</em>; stored values are unchanged.',s_over,cols=3)}
{R.add('S_paired','S','Lotz-like versus Conselice-like ÷ 10 clumpiness for the same galaxies',
 'Tests whether the two profiles rank individual galaxies in the same order.',
 'Closeness of the cloud to the dashed equality line and the Spearman ρ in each title. Top row raw, bottom row corrected. The corrected panels show more scatter because both sky estimates add noise.',
 'Galaxies with both profiles finite in that band; colour is the number of galaxies per cell on a logarithmic scale.',s_pair,cols=3)}
<h3 id="s-galfit">6.2 Clumpiness and fitted structure (GALFIT)</h3>
<p>Raw clumpiness behaves like raw asymmetry: it rises toward fainter mean surface brightness (Lotz-like ρ = {trio(lambda b:gr(L+'_raw','gf_mu_mean',b))}, Conselice-like ρ = {trio(lambda b:gr(K+'_raw','gf_mu_mean',b))}; {fig('S_galfit_L','c')}, {fig('S_galfit_K','c')}). After correction the relation becomes ρ = {trio(lambda b:gr(L+'_corrected','gf_mu_mean',b))} (Lotz-like) and {trio(lambda b:gr(K+'_corrected','gf_mu_mean',b))} (Conselice-like). As for asymmetry, this strongly suggests that pixel noise contributes much of the raw high-frequency residual in these faint galaxies.</p>
{R.add('S_galfit_L','S','Lotz-like clumpiness, raw (left) and corrected (right), against six GALFIT properties',
 'Shows whether the Lotz-like clumpiness follows structural parameters or mainly brightness.',
 'Rows are GALFIT properties; left raw, right corrected. The surface-brightness row (c, d) is the main noise diagnostic; the size row (e, f) shows whether small galaxies, where the blur kernel is only a few pixels wide, behave differently.',
 D_GF,gal_pairs([L+'_raw',L+'_corrected']),cols=2)}
{R.add('S_galfit_K','S','Conselice-like clumpiness, raw (left) and corrected (right), against six GALFIT properties',
 'The same comparison as '+fig('S_galfit_L')+' for the second profile, laid out identically so the two can be compared panel by panel.',
 'Whether the trends match those of the Lotz-like profile after allowing for the ×10 scale on the y axis.',
 D_GF,gal_pairs([K+'_raw',K+'_corrected']),cols=2)}
<h3 id="s-classes">6.3 Clumpiness and human / DECaLS10 classifications</h3>
<p>Clumpy light is expected in star-forming spirals and irregulars. The strongest r-band associations are: Lotz-like raw – {strongest(L+'_raw')}; Lotz-like corrected – {strongest(L+'_corrected')}; Conselice-like raw – {strongest(K+'_raw')}; Conselice-like corrected – {strongest(K+'_corrected')} ({fig('S_assoc')}). {fig('S_human')}, {fig('S_effects')} and {fig('S_trends')} show distributions, class separations and trends for all four outputs; {fig('S_embed')} shows them on the clustering map.</p>
{R.add('S_human','S','Clumpiness for each human morphology consensus class (r band)',
 'Tests whether human-identified LTG/spiral or irregular galaxies are clumpier than ETG or Featureless ones.',
 'Panel order: (a) Lotz-like raw, (b) Lotz-like corrected, (c) Conselice-like raw, (d) Conselice-like corrected. Look for class shifts present in both profiles; shifts present only in raw values may reflect brightness differences between classes. g and i bands are in the expandable block.',
 D_HH,[human_hist(m,'r') for m in sv],cols=2,extra=[('g and i bands',[human_hist(m,b) for b in ('g','i') for m in sv],2)])}
{R.add('S_assoc','S','Clumpiness: rank correlations with every human answer and DECaLS10 class',
 'One heatmap per stored clumpiness output, identically laid out, so profile and correction effects on each association can be compared cell by cell.',
 'Panel order as in '+fig('S_human')+'. Compare (a) with (c) and (b) with (d) for profile agreement, and raw with corrected for the effect of the sky term; check the GALFIT-adjusted block for associations that survive brightness/size/shape control.',
 D_AS,[assoc_panel(m) for m in sv],cols=1)}
{R.add('S_effects','S','Clumpiness: separation of each class from all other galaxies',
 'Shows, for hard class labels, whether members tend to be clumpier (positive) or smoother (negative) than non-members, for each stored output.',
 'Panel order as in '+fig('S_human')+'. Consistent colour across all four panels indicates a separation that does not depend on the profile or correction choice.',
 D_EF,[effect_panel(m) for m in sv],cols=1)}
{R.add('S_trends','S','Clumpiness along key human fractions and DECaLS10 probabilities',
 'Shows the shape of the clumpiness relation with visible nucleus, ETG answer and DECaLS10 Disturbed probability, for all four outputs.',
 'Each group of four panels is one target, ordered Lotz-like raw, Lotz-like corrected, Conselice-like raw, Conselice-like corrected. The remaining 14 targets are in the expandable block.',
 D_PT,[trend(m,t) for t in KEY_TARGETS for m in sv],cols=2,extra=[('All other human answers and DECaLS10 classes',all_trends(sv),2)])}
{R.add('S_embed','S','Clumpiness across a saved clustering map',
 'Shows where clumpy galaxies lie on the clustering map, for each stored output.',
 'Panel order as in '+fig('S_human')+'. Compare raw with corrected: regions that are bright in raw clumpiness but not in corrected clumpiness are likely low-surface-brightness regions of the map.',
 D_EM,[p for m in sv for p in emb_panels(m)],cols=2,extra=[('g and i bands',[p for b in ('g','i') for m in sv for p in emb_panels(m,(b,))],2)])}
</section>'''

    # ======================================================== supplementary diagnostics
    dm=['Gini_aperture_v1','M20_aperture_v1','axis_ratio_positive_aperture_v1']
    d_hist=[];d_full=[]
    for m,lims in [('Gini_aperture_v1',None),('M20_aperture_v1',None),('axis_ratio_positive_aperture_v1',(0,1)),('ellipticity_positive_aperture_v1',(0,1)),('position_angle_positive_aperture_v1',(0,180))]:
        d_hist.append((P.band_histogram(frame,m,lims or robust_limits([frame[m]]),f'{m}_dist',PLAIN[m]),PLAIN[m]))
        if lims is None:d_full.append((P.band_histogram(frame,m,robust_limits([frame[m]],quantiles=(0,1)),f'{m}_dist_full',PLAIN[m]+': full range'),PLAIN[m]+', full range'))
    gm=[];gm_rho={}
    for b in BANDS:
        path,info=P.density(frame,b,'Gini_aperture_v1','M20_aperture_v1',robust_limits([frame.Gini_aperture_v1]),robust_limits([frame.M20_aperture_v1]),f'gini_m20_{b}',f'{b} band');gm.append((path,f'{b} band'));gm_rho[b]=info['rho']
    qq=[(P.density(unfl,b,'gf_q','axis_ratio_positive_aperture_v1',(.3,1),(0,1),f'q_compare_{b}',f'{b} band',equality=True)[0],f'{b} band') for b in BANDS]
    gc=[(P.density(frame,b,'C','Gini_aperture_v1',robust_limits([frame.C]),robust_limits([frame.Gini_aperture_v1]),f'C_vs_gini_{b}',f'{b} band')[0],f'{b} band') for b in BANDS]
    sec_d=f'''<section id="diagnostics" class="family" style="--fam:{FAMILY_STYLE['D'][0]}"><h2>7. Supplementary diagnostics: Gini, M<sub>20</sub> and shape {badge('D')}</h2>
<p class="lead">These quantities are <strong>not</strong> part of the original three CAS indices. They are extra descriptions of the same aperture light, kept in a separate section so they are not confused with C, A or S. <a href="#g-gini">Gini</a> measures how unequally the light is shared among pixels {cite('abraham2003','lotz2004')}; <a href="#g-m20">M<sub>20</sub></a> measures how far from the centre the brightest 20% of the light lies {cite('lotz2004')}; the <a href="#g-axisratio">positive-flux axis ratio</a> q<sub>+</sub> measures how elongated the light is. All are aperture adaptations and differ from the segmentation-map versions in the literature.</p>
<h3>7.1 Distributions and relations to CAS</h3>
<p>{fig('D_dist')} shows the distributions. Gini and M<sub>20</sub> are anti-correlated (ρ={gm_rho['g']:.2f}/{gm_rho['r']:.2f}/{gm_rho['i']:.2f}; {fig('D_gm20')}), as in bright galaxies {cite('lotz2004')}. Surprisingly, Gini is almost unrelated to concentration (ρ={pc.loc[('g','C vs Gini_aperture_v1'),'spearman']:.2f} in g; {fig('D_cgini')}), although both describe how centrally the light is distributed. The explanation is in {fig('D_galfit')}: Gini is dominated by surface brightness (ρ = {trio(lambda b:gr('Gini_aperture_v1','gf_mu_mean',b))} with mean surface brightness), because in faint galaxies the absolute-flux ranking of pixels is largely a ranking of noise. The aperture axis ratio agrees only loosely with the GALFIT model axis ratio ({fig('D_qcompare')}).</p>
{R.add('D_dist','D','Distributions of Gini, M20 and positive-flux shape',
 'Reference distributions for the supplementary diagnostics, which appear in the GALFIT, classification and embedding figures of this section.',
 'Gini (a): 0 means all pixels equally bright, 1 means all light in one pixel. M<sub>20</sub> (b): more negative means the brightest light is more central. Shape (c–e): q<sub>+</sub> near 1 means round. The image angle (e) should be flat for randomly oriented galaxies; peaks would indicate an image-orientation artefact.',
 D_ALL+' Natural ranges are used for q<sub>+</sub>, ellipticity (0–1) and angle (0–180°). Full ranges for Gini and M<sub>20</sub> are in panels (f)–(g).',d_hist,cols=3,extra=[('Full finite ranges for Gini and M20',d_full,2)])}
{R.add('D_gm20','D','Gini versus M20 for the same galaxies',
 'Gini–M<sub>20</sub> diagrams are widely used to separate normal from merging galaxies in bright samples '+cite('lotz2004')+'. This figure shows where this faint catalog lies in that plane.',
 'The overall anti-correlation and whether there is a separate tail of high-Gini/high-M<sub>20</sub> galaxies. No published classification lines are drawn, because our aperture definitions differ from the ones they were calibrated on.',
 'Galaxies with both values finite in that band; colour shows galaxies per cell on a logarithmic scale.',gm,cols=3)}
{R.add('D_cgini','D','Concentration versus Gini for the same galaxies',
 'Both C and Gini are often called “concentration-like”. This figure checks whether they agree in this sample.',
 'A shapeless cloud (ρ close to zero) means the two indices respond to different properties of these faint images.',
 'Galaxies with both values finite in that band; colour shows galaxies per cell on a logarithmic scale.',gc,cols=3)}
{R.add('D_qcompare','D','Positive-flux aperture axis ratio versus GALFIT axis ratio',
 'q<sub>+</sub> is computed from image moments with negative (noise) pixels set to zero, which biases faint shapes toward roundness. Comparing it with the fitted model axis ratio checks how much shape information it retains.',
 'Distance from the dashed equality line. Note the GALFIT axis ratio is constrained to 0.3–1 and shared across bands in the fit.',
 'Main unflagged GALFIT cohort.',qq,cols=3)}
<h3>7.2 Relations to GALFIT and to classifications</h3>
<p>{fig('D_galfit')} compares the three diagnostics with all GALFIT properties; {fig('D_human')}, {fig('D_assoc')}, {fig('D_effects')} and {fig('D_trends')} relate them to the human and DECaLS10 labels, and {fig('D_embed')} shows them on the clustering map. In r, Gini associates most strongly with {strongest('Gini_aperture_v1')}; after GALFIT adjustment the visible-nucleus correlation of Gini is {ar('Gini_aperture_v1','p_visible_nucleus',cohort='galfit_adjusted',stat='partial_rank'):.3f} (unadjusted {ar('Gini_aperture_v1','p_visible_nucleus'):.3f}).</p>
{R.add('D_galfit','D','Gini (left), M20 (middle) and axis ratio q+ (right) against six GALFIT properties',
 'Shows which fitted properties each supplementary diagnostic follows, especially their sensitivity to surface brightness.',
 'Rows are GALFIT properties, columns the three diagnostics. The surface-brightness row (d–f) shows the strong brightness dependence of Gini.',
 D_GF,gal_pairs(dm),cols=3)}
{R.add('D_human','D','Gini, M20 and q+ for each human morphology consensus class (r band)',
 'Tests whether the human morphology classes differ in these diagnostics.',
 'Class shifts in (a) Gini, (b) M<sub>20</sub> and (c) q<sub>+</sub>. g and i are in the expandable block.',
 D_HH,[human_hist(m,'r') for m in dm],cols=3,extra=[('g and i bands',[human_hist(m,b) for b in ('g','i') for m in dm],3)])}
{R.add('D_assoc','D','Gini, M20 and q+: rank correlations with every human answer and DECaLS10 class',
 'The same overview as for C, A and S, for the supplementary diagnostics.',
 'Panels (a) Gini, (b) M<sub>20</sub>, (c) q<sub>+</sub>. Note that the visible-nucleus association of Gini weakens after GALFIT adjustment (rightmost block). Ellipticity is exactly 1−q<sub>+</sub> and image angle is not a linear quantity, so neither has its own panel.',
 D_AS,[assoc_panel(m) for m in dm],cols=1)}
{R.add('D_effects','D','Gini, M20 and q+: separation of each class from all other galaxies',
 'Class-label separations for the supplementary diagnostics.',
 'Panels (a) Gini, (b) M<sub>20</sub>, (c) q<sub>+</sub>. Edge-on DECaLS10 classes should separate in q<sub>+</sub> if it carries real shape information.',
 D_EF,[effect_panel(m) for m in dm],cols=1)}
{R.add('D_trends','D','Gini, M20 and q+ along key human fractions and DECaLS10 probabilities',
 'Shape of the relation with visible nucleus, ETG answer and DECaLS10 Disturbed probability.',
 'Each row is one target; columns are Gini, M<sub>20</sub> and q<sub>+</sub>. Other targets are in the expandable block.',
 D_PT,[trend(m,t) for t in KEY_TARGETS for m in dm],cols=3,extra=[('All other human answers and DECaLS10 classes',all_trends(dm),3)])}
{R.add('D_embed','D','Gini, M20 and q+ across a saved clustering map',
 'Shows where these diagnostics vary across the clustering map.',
 'Compare (a) Gini with the A<sub>raw</sub> map in '+fig('A_embed')+': similar patterns would point to a shared brightness dependence.',
 D_EM,[p for m in dm for p in emb_panels(m)],cols=3,extra=[('g and i bands',[p for b in ('g','i') for m in dm for p in emb_panels(m,(b,))],3)])}
</section>'''

    # ======================================================== overview + explorer
    allm=['C','A_raw','A_corrected',*sv,'Gini_aperture_v1','M20_aperture_v1']
    mats=[(P.galfit_matrix(corr,b,allm,f'galfit_matrix_{b}',f'{b} band: morphology vs GALFIT'),f'{b} band') for b in BANDS]
    main_c=corr[corr.cohort=='gf_unflagged'].set_index(['band','metric','property']);int_c=corr[corr.cohort=='gf_interior'].set_index(['band','metric','property'])
    delta=(main_c.spearman-int_c.spearman).abs().dropna()
    explorer=(root/'association_explorer.html').read_text()
    explorer=explorer.replace('</style>','body{padding:16px;box-sizing:border-box;background:#fff;color:#222}</style>')
    explorer=explorer.replace('</body>',"<script>new ResizeObserver(()=>parent.postMessage({type:'morphology-explorer-height',height:Math.ceil(document.body.getBoundingClientRect().height)},'*')).observe(document.body);</script></body>")
    sens=fitsens.rename(columns={'band':'Band','property':'GALFIT property','paired':'Galaxies','median_primary_minus_constant':'Median primary − alternate','p16':'16th pct','p84':'84th pct'})
    sens['GALFIT property']=sens['GALFIT property'].map(PLAIN)
    sec_x=f'''<section id="overview"><h2>8. All metrics together {badge('X')}</h2>
<p>Sections 4–7 deliberately treat one metric family at a time. This section puts them next to each other. {fig('galfit_matrix')} condenses every GALFIT trend of {fig('C_galfit')}, {fig('A_galfit')}, {fig('S_galfit_L')}, {fig('S_galfit_K')} and {fig('D_galfit')} into one number per pair. Removing galaxies whose fitted n or q lies at a fit boundary changes these correlations by a median of {delta.median():.3f} (maximum {delta.max():.3f}), so they are not driven by boundary fits.</p>
{R.add('galfit_matrix','X','Rank correlations of every morphology output with every GALFIT property',
 'A compact summary for comparing families: which metrics follow brightness, which follow profile shape, and which follow size.',
 'Rows are grouped by family (white separators): C; A raw/corrected; four S outputs; Gini and M<sub>20</sub>. Columns are GALFIT properties. A strong mean-surface-brightness column for raw A and raw S, and for Gini, versus a strong Sérsic-n cell for C, is the key contrast.',
 'Main unflagged GALFIT cohort; every finite pair, no display cropping. Exact values and pair counts: <code>paper_v2_data/galfit_correlations.csv</code>.',mats,cols=1)}
<h3 id="explorer">8.1 Interactive correlation explorer</h3>
<p>The explorer below re-draws any single column of {fig('C_assoc')}, {fig('A_assoc')}, {fig('S_assoc')} and {fig('D_assoc')} as bars, for every morphology output at once. Choose a human answer or DECaLS10 class, a band and a cohort. <a href="association_explorer.html">Open it separately</a>.</p>
<iframe id="association-explorer-frame" title="Interactive morphology and classification rank correlations" style="display:block;width:100%;height:1200px;border:0;border-radius:8px;background:white" srcdoc="{escape(explorer,quote=True)}"></iframe>
<h3>8.2 Sensitivity to the GALFIT configuration</h3><p>A smaller alternate GALFIT export (radius-order-1 / constant offset) is compared on common galaxies with unflagged fits in both. Differences are primary minus alternate. This checks how much the fitted properties used throughout Sections 4–8 depend on the fit configuration.</p>{table_html(sens)}
</section>'''

    # ======================================================== embedding atlas
    spaces=[];key_of={}
    for e in emb_index:
        k=(e['title'],e['projection'])
        if k not in key_of:key_of[k]=len(spaces);spaces.append(dict(title=e['title'],projection=e['projection'],bands={}))
        spaces[key_of[k]]['bands'][e['band']]=e['index']
    atlas_blocks=''
    for fam in ['C','A','S','D']:
        cells=''.join(f'<div class="panel"><a data-metric="{m}" href="paper_v5_assets/emb_{default_space:02d}_{m}.png" target="_blank" rel="noopener"><img loading="lazy" data-metric="{m}" src="paper_v5_assets/emb_{default_space:02d}_{m}.png" alt="{escape(PLAIN[m])} on the selected embedding"></a><p class="panel-cap" data-metric-cap="{m}">{escape(PLAIN[m])}</p></div>' for m in FAMILIES[fam]['metrics'])
        atlas_blocks+=f'<div class="atlas-family" style="--fam:{FAMILY_STYLE[fam][0]}">{badge(fam)}<div class="panel-grid cols-{3 if fam=="D" else 2}">{cells}</div></div>'
    space_opts=''.join(f'<option value="{i}"'+(' selected' if s['title']==emb[default_space]['title'] and s['projection']==emb[default_space]['projection'] else '')+f'>{escape(s["title"])} · {escape(s["projection"])}</option>' for i,s in enumerate(spaces))
    sec_e=f'''<section id="embedding-atlas"><h2>9. Embedding atlas: every saved map, every metric</h2>
<p>The <a href="#g-umap">PCA and UMAP</a> maps were computed by the clustering project from GALFIT properties alone or combined with DECaLS10 / human-trained Zoobot features; they were <em>not</em> computed from these morphology values {cite('jolliffe2016','mcinnes2018')}. There are 7 feature families × 3 projections × 3 bands = {len(emb)} saved maps. Choose a map and a band; the panels below always show the same ten morphology outputs, grouped by family. Each map is drawn as single-metric panels (hexagonal cells coloured by the median value), replacing the earlier ten-panel composite images. The <em>default</em> map is the one used in {fig('C_embed')}, {fig('A_embed')}, {fig('S_embed')} and {fig('D_embed')}.</p>
<p><strong>How to read the maps.</strong> Position and distance on a UMAP map are not physical; only neighbourhoods are meaningful. A smooth colour gradient means a metric varies systematically with the clustering features. A map built partly from human-trained features is not an independent test of human labels. All {emb[0]['rows']:,} saved points are drawn (grey where too few finite values).</p>
<div class="atlas-controls"><label>Saved map<select id="atlas-space">{space_opts}</select></label><label>Band<select id="atlas-band"><option>g</option><option selected>r</option><option>i</option></select></label></div>
<p id="atlas-caption" class="meta"></p>{atlas_blocks}</section>'''

    # ======================================================== front matter, glossary, methods, discussion
    neg_a=stats[stats.metric=='A_corrected'].set_index('band');neg_a=100*neg_a.negative/neg_a.valid
    executive=f'''<section id="summary"><h2>Summary</h2>
<p>This report describes the morphology of all <strong>84,200</strong> galaxy candidates in the DES low-surface-brightness catalog, measured in the g, r and i images (252,600 object–band rows; 96 cutouts missing). The core of the analysis is the classical <strong>CAS</strong> system — Concentration, Asymmetry and clumpiness (S) {cite('conselice2003')} — which has its own section each. Supplementary diagnostics (Gini, M<sub>20</sub>, shape) are reported separately. Every measured quantity is compared with fitted structure (GALFIT), human classifications, the DECaLS10 machine classifier and saved clustering maps.</p>
<div class="summary-grid">
<div class="card" style="--fam:{FAMILY_STYLE['C'][0]}">{badge('C')}<p><b>Concentration behaves as expected.</b> It follows the fitted Sérsic index (ρ≈{gr('C','gf_n','r'):.2f} in r), depends little on surface brightness and tracks the human visible-nucleus answer (ρ={ar('C','p_visible_nucleus'):.2f}). See {fig('C_galfit')} and {fig('C_assoc')}.</p><a href="#concentration">Go to Section 4 →</a></div>
<div class="card" style="--fam:{FAMILY_STYLE['A'][0]}">{badge('A')}<p><b>Asymmetry is dominated by noise.</b> The sky term is about as large as the raw value; {neg_a.min():.0f}–{neg_a.max():.0f}% of corrected values are negative, and raw A rises strongly toward fainter galaxies (ρ≈{gr('A_raw','gf_mu_mean','r'):.2f}). See {fig('A_dist')} and {fig('A_galfit')}.</p><a href="#asymmetry">Go to Section 5 →</a></div>
<div class="card" style="--fam:{FAMILY_STYLE['S'][0]}">{badge('S')}<p><b>The two clumpiness profiles agree in ranking.</b> After removing the ×10 factor they rank galaxies almost identically (ρ={pc.loc[('r','S raw'),'spearman']:.3f} raw in r), but like A they are strongly affected by the sky term. See {fig('S_overlay')} and {fig('S_paired')}.</p><a href="#clumpiness">Go to Section 6 →</a></div>
<div class="card" style="--fam:{FAMILY_STYLE['D'][0]}">{badge('D')}<p><b>Gini mostly measures brightness here.</b> It correlates with mean surface brightness at ρ≈{gr('Gini_aperture_v1','gf_mu_mean','r'):.2f} and is nearly unrelated to C. See {fig('D_galfit')} and {fig('D_cgini')}.</p><a href="#diagnostics">Go to Section 7 →</a></div>
</div>
<p>All values are descriptive. The catalog is a supplied candidate list, not a complete or representative galaxy sample, and none of the indices has been calibrated for DES noise and PSF; physical interpretation therefore needs the tests listed in <a href="#discussion">Section 10</a>.</p></section>'''

    howto=f'''<section id="howto"><h2>1. How to read this report</h2>
<p><strong>Colour tags.</strong> Every section heading and every figure carries a tag naming the metric family it is about: {badge('C')} {badge('A')} {badge('S')} {badge('D')}. Shared and cross-metric material is tagged {badge('P')} or {badge('X')}. The coloured stripe on the left edge of each figure matches its tag.</p>
<p><strong>Figures and panels.</strong> Each figure is a set of single plots (panels) labelled (a), (b), (c)… and arranged in a grid that collapses to one column on narrow screens. Click any panel to open it at full resolution. Some figures have an expandable block (“g and i bands”, “Full finite range”, “All other …”) with additional panels that continue the lettering. Below each figure, three short paragraphs explain <em>why</em> it is shown, <em>what to look for</em>, and <em>which data</em> were used; an expandable table lists the display window and the number of values outside it for every panel. Figures are always referenced from the text.</p>
<p><strong>Bands.</strong> In all plots the photometric bands use the same colours and line styles: <span style="color:{BAND_COLORS['g']}">g solid</span>, <span style="color:{BAND_COLORS['r']}">r dashed</span>, <span style="color:{BAND_COLORS['i']}">i dotted</span>. g is the bluest band, i the reddest.</p>
<p><strong>Display windows.</strong> Most axes show the pooled 0.5th–99.5th percentiles plus 3.5% padding, so the bulk of the data is visible. Nothing is deleted: histogram percentages are of the <em>whole</em> finite sample (the visible bars need not add to 100%), full ranges are available, and every off-screen count is tabulated. Medians, percentiles and correlations are always computed before any display cropping.</p>
<p><strong>Samples.</strong> Different comparisons need different galaxies, so each figure states its sample. The table summarises them.</p>
{table_html(pd.DataFrame([
 {'Sample':'Morphology catalog','Galaxies':'84,200 per band','Used for':'All distributions (Figures with “Distribution” in their title)'},
 {'Sample':'GALFIT main (unflagged) cohort','Galaxies':f'{gmain["g"]:,} / {gmain["r"]:,} / {gmain["i"]:,} (g/r/i)','Used for':'All GALFIT trends and correlation matrices'},
 {'Sample':'Human labels, ≥2 morphology answers','Galaxies':f'{h2:,}','Used for':'Human consensus distributions, correlations, class separations'},
 {'Sample':'DECaLS10 valid probability vectors','Galaxies':f'{zv:,}','Used for':'DECaLS10 correlations and class separations'},
 {'Sample':'Saved clustering selection','Galaxies':f'{emb[0]["rows"]:,}','Used for':'Embedding maps (Section 9)'},
]))}
</section>'''

    glossary=f'''<section id="glossary"><h2>2. Key concepts in plain language</h2>
<p>This section explains every technical term used in the report in everyday language, with references for readers who want the details. Terms are linked from the text where they first matter.</p>
<h3 id="g-image">Galaxy images, bands and sky subtraction</h3><p>The images come from the Dark Energy Survey (DES), which photographed a large part of the southern sky through several colour filters {cite('abbott2018')}. We use three: <b>g</b> (blue-green light), <b>r</b> (red) and <b>i</b> (near-infrared). The galaxies are candidates from the DES low-surface-brightness catalog {cite('tanoglidis2021')}. Every image also contains the faint glow of the night sky; the image providers estimated and subtracted it, so a pixel with no galaxy light should be near zero, scattering randomly above and below. That random scatter is <b>noise</b>. Because these galaxies are very faint, their light per pixel is often comparable to the noise, which is the single most important fact for interpreting this report. Light from unrelated stars and galaxies is excluded with a <b>mask</b>; masked pixels are ignored, never set to zero.</p>
<h3 id="g-psf">Seeing and the point-spread function (PSF)</h3><p>The atmosphere and telescope blur every image slightly, so a point-like star appears as a small fuzzy disc. The shape of that blur is called the point-spread function. Blurring smooths out small details and makes galaxies look more centrally concentrated or smoother than they are, especially galaxies that are only a few pixels across {cite('binney1998','conselice2000')}. The g, r and i images have slightly different blurring.</p>
<h3 id="g-petrosian">Petrosian radius</h3><p>To measure a galaxy you first need to know how big it is, but galaxies have no sharp edge — they fade gradually into the sky. Petrosian {cite('petrosian1976')} proposed a clever way to define a size from the image itself: walk outwards from the centre, and at each radius compare the brightness of a thin ring with the average brightness of everything inside that ring. Near the centre the ring is about as bright as the average; far out it becomes much fainter. The <b>Petrosian radius</b> R<sub>P</sub> is where the ring has dropped to a chosen fraction (here η=0.2, i.e. 20%) of the inside average. Because it uses a <em>ratio</em> of brightnesses, the same galaxy placed farther away or observed with a different exposure gets the same R<sub>P</sub> in angular units, which is why it is widely used for galaxy surveys {cite('graham2005p')}. All our indices are measured inside a circle of radius 1.5 R<sub>P</sub>, see <a href="#aperture">Section 3</a>.</p>
<h3 id="g-concentration">Concentration (C)</h3><p>Imagine collecting a galaxy’s light in growing circles. r<sub>20</sub> is the circle that contains 20% of the light and r<sub>80</sub> the one that contains 80%. If most light is in a bright centre, r<sub>20</sub> is tiny compared to r<sub>80</sub>; if the light is spread evenly, they are closer. C=5 log<sub>10</sub>(r<sub>80</sub>/r<sub>20</sub>) turns this ratio into a convenient number: roughly 3 for disc galaxies and above 4 for ellipticals in bright, nearby galaxy samples {cite('kent1985','bershady2000','conselice2003')}.</p>
<h3 id="g-asymmetry">Asymmetry (A) and the sky correction</h3><p>Rotate the image of a galaxy by half a turn and subtract it from the original. A perfectly symmetric galaxy leaves nothing; a lopsided or disturbed one leaves residual light. A<sub>raw</sub> is the leftover divided by the total {cite('abraham1996','conselice2000')}. But noise is never symmetric, so even a perfect galaxy leaves some residual. To estimate how much, the same calculation is done on an empty patch of sky (A<sub>sky</sub>), scaled to the galaxy’s size and brightness, and subtracted: A<sub>corr</sub>=A<sub>raw</sub>−A<sub>sky</sub> {cite('conselice2003')}. If the galaxy is barely brighter than the noise, A<sub>corr</sub> can be negative. That is a legitimate result of an uncertain estimate, not an error {cite('thorp2021')}.</p>
<h3 id="g-clumpiness">Clumpiness (S)</h3><p>Blur the galaxy image a little, then subtract the blurred version from the original. Smooth light disappears; small bright knots (for example star-forming regions) remain. Clumpiness is the remaining positive light divided by the total light, after excluding the centre (which is always peaked) {cite('conselice2003','lotz2004')}. Results depend on how much you blur and how much of the centre you exclude, which is why two named versions are kept (<a href="#clumpiness">Section 6</a>). As for asymmetry, an empty-sky estimate is subtracted to give the corrected value, which can be negative.</p>
<h3 id="g-gini">Gini coefficient (G)</h3><p>Borrowed from economics, where it measures how unequally wealth is shared in a population. Here the “population” is the pixels and the “wealth” is their brightness: G=0 if all pixels are equally bright, G=1 if one pixel has all the light {cite('abraham2003','lotz2004')}. We compute it on absolute pixel values inside our circular aperture, which differs from the usual definition on a segmentation map.</p>
<h3 id="g-m20">M<sub>20</sub></h3><p>Take the brightest pixels that together contain 20% of the light and measure how spread out they are around the centre, relative to all the light. A value near −2 or below means the brightest light is concentrated in the middle; values closer to 0 mean bright spots are off-centre, as in a merger or a galaxy with bright clumps {cite('lotz2004')}.</p>
<h3 id="g-axisratio">Axis ratio and ellipticity</h3><p>The axis ratio q is the short axis of the galaxy’s outline divided by its long axis: 1 for round, small for elongated. Ellipticity is 1−q. Our aperture version q<sub>+</sub> ignores negative pixels, so in noisy images it tends toward roundness; the GALFIT q is from a model fit.</p>
<h3 id="g-sersic">Sérsic profile and Sérsic index n</h3><p>Galaxy brightness typically falls off smoothly from the centre. Sérsic {cite('sersic1963')} described this with one formula, I(R)∝exp[−b<sub>n</sub>(R/R<sub>e</sub>)<sup>1/n</sup>], where the index n controls the shape: n≈1 is an exponential disc (gently falling), n≈4 a classical elliptical (sharply peaked centre with extended wings), n&lt;1 an even flatter core {cite('graham2005')}. n is thus a model-based counterpart of concentration.</p>
<h3 id="g-reff">Effective radius R<sub>e</sub></h3><p>The radius containing half of the light of the fitted Sérsic model, measured along the long axis. It is a property of the model, whereas R<sub>P</sub> is measured directly from the image; for an exponential disc R<sub>P</sub>≈2R<sub>e</sub> {cite('graham2005p')}.</p>
<h3 id="g-sb">Magnitude and surface brightness</h3><p>Astronomers measure brightness in <b>magnitudes</b>, a logarithmic scale where <em>larger numbers mean fainter</em>: 5 magnitudes is a factor of 100 in brightness. <b>Surface brightness</b> (mag arcsec<sup>−2</sup>) is brightness per patch of sky — how bright the galaxy looks per unit area, independent of its distance in the nearby universe. Low-surface-brightness galaxies have most of their light spread thinly, close to the brightness of the night sky {cite('binney1998')}. We use the GALFIT <em>mean effective surface brightness</em> ⟨μ⟩<sub>e</sub>, the average within R<sub>e</sub> {cite('graham2005')}.</p>
<h3 id="g-galfit">GALFIT, GALFITM and reduced χ²</h3><p>GALFIT fits a mathematical model (here a single Sérsic profile) to a galaxy image by adjusting position, brightness, size, n and shape until the model matches the pixels as closely as possible {cite('peng2002','peng2010')}; GALFITM fits several bands together {cite('haussler2013')}. The fits were done by the companion <code>nucleated_lsbg</code> project; we only read them. The reduced χ² measures how well the model fits relative to the noise (≈1 is a good fit). The “unflagged” cohort keeps fits that succeeded and raised no parameter warnings.</p>
<h3 id="g-spearman">Spearman rank correlation ρ</h3><p>Sort the galaxies by one quantity and separately by another; ρ measures how similar the two orderings are {cite('spearman1904')}. ρ=+1: perfectly the same order; −1: exactly reversed; 0: unrelated. It does not assume a straight-line relation and is not upset by a few extreme values. As a rough guide, |ρ|&lt;0.1 is negligible, 0.1–0.3 weak, 0.3–0.5 moderate and &gt;0.5 strong. With tens of thousands of galaxies even weak values are statistically precise, so the question is whether they are <em>meaningful</em>, not whether they are “significant”.</p>
<h3 id="g-partial">Partial rank correlation (GALFIT-adjusted)</h3><p>Two quantities can be correlated only because both depend on a third (e.g. both depend on brightness). The partial rank correlation first removes, from both rankings, the part that can be predicted linearly from the ranks of mean surface brightness, effective radius and axis ratio, and then correlates what remains. It removes only that simple linear part, not all shared influences.</p>
<h3 id="g-rankbiserial">Rank-biserial separation</h3><p>For a yes/no class, pick one class member and one non-member at random. The rank-biserial value is the probability that the member has the larger metric value minus the probability that it has the smaller one {cite('cureton1956','kerby2014')}. +1 means every member is above every non-member, 0 no separation. It is computed from the Mann–Whitney statistic {cite('mannwhitney1947')}.</p>
<h3 id="g-percentile">Medians, percentiles and spreads</h3><p>The median is the middle value: half the galaxies are above it. The 16th and 84th percentiles enclose the middle 68% of galaxies (for a bell-shaped distribution, roughly one standard deviation either side). Shaded bands in trend plots show this galaxy-to-galaxy spread, which is much wider than the uncertainty of the median itself.</p>
<h3 id="g-zoobot">Human labels, Zoobot and DECaLS10</h3><p>Human reviewers answered questions for {h2:,} galaxies (morphology type, visible nucleus, low surface brightness, fit quality). An answer fraction is the share of reviewers who chose that answer; a consensus class follows the published weighted-vote rule. Zoobot is a deep-learning model trained on Galaxy Zoo volunteer classifications {cite('walmsley2022','walmsley2023')}; the DECaLS10 version outputs probabilities for ten classes such as “Disturbed” or “Round smooth” {cite('decals10')}. It was trained on a different survey (DECaLS), so its outputs for DES galaxies are model judgements, not ground truth.</p>
<h3 id="g-umap">PCA and UMAP maps</h3><p>When each galaxy is described by many numbers, it helps to draw them on a 2-D map in which similar galaxies lie close together. PCA does this by finding the directions of greatest variation {cite('jolliffe2016')}; UMAP builds a map that preserves local neighbourhoods {cite('mcinnes2018')}. The axes of these maps have no physical meaning, and in UMAP the distances between distant groups are not reliable.</p>
</section>'''

    eq=equation(r'$C=5\log_{10}(r_{80}/r_{20}),\quad A_{\rm raw}=\frac{\sum|I-I_{180}|}{\sum|I|},\quad A_{\rm corr}=A_{\rm raw}-A_{\rm sky},\quad S=f\,\frac{\sum\max(I-\widetilde I,0)}{\sum I}-S_{\rm sky}$')
    methods=f'''<section id="methods"><h2>Measurement method in brief</h2>
<p>All measurements use the shared CAS implementation <code>cas_discrete_v3</code> with both S profiles at version v1 and the aperture diagnostics at <code>aperture_diagnostics_v1</code>. Inputs are signed, sky-subtracted 256-pixel DES cutouts; contaminants are masked, the target and usable sky are kept. R<sub>P</sub> is the first descending η=0.2 crossing; the aperture is 1.5 R<sub>P</sub>; the asymmetry centre is found by bounded minimisation; sky terms come from the cleanest of four corner tiles that is ≥95% valid and outside the aperture. When no such tile exists the corrected values are missing, not zero. Full definitions, flags and literature differences are in <code>docs/CAS_METHODS.md</code> and Appendix A.</p>{eq}
<p>This revision (v5) changes only presentation: it re-draws the frozen v2/v3 analysis tables and the morphology table as single-panel figures. No measurement, profile, cohort or statistic is redefined; the input SHA-256 values are recorded in <code>report_v5_manifest.json</code>.</p></section>'''

    discussion=f'''<section id="discussion"><h2>10. Interpretation, limitations and next tests</h2>
<h3>{badge('C')} Concentration</h3><p>C is the most robust of the three indices in this sample. It agrees with the fitted Sérsic index ({fig('C_galfit')}) and depends little on surface brightness. Its absolute values are still affected by the PSF for small galaxies and by the choice of aperture, so comparisons with published C values need matched settings.</p>
<h3>{badge('A')} Asymmetry</h3><p>The sky term is of the same size as the raw signal ({fig('A_dist')}), and raw A follows surface brightness ({fig('A_galfit')}). Raw A is therefore largely a noise indicator in these faint galaxies, and corrected A is a small difference of two noisy numbers. The choice of the minimum-scatter sky corner, and the different interpolation used for galaxy and sky, are possible sources of bias. Injection tests with realistic DES noise and PSF are needed before either version is used to identify disturbed galaxies {cite('thorp2021')}.</p>
<h3>{badge('S')} Clumpiness</h3><p>The two profiles rank galaxies almost identically ({fig('S_paired')}) but, like asymmetry, are strongly affected by the sky term ({fig('S_galfit_L')}, {fig('S_galfit_K')}). Agreement between profiles does not make either correct; the same noise affects both. Neither should be compared with published thresholds without matching every convention.</p>
<h3>{badge('D')} Supplementary diagnostics</h3><p>Aperture Gini mostly follows surface brightness ({fig('D_galfit')}), and q<sub>+</sub> is biased toward round shapes. No published Gini–M<sub>20</sub> classification boundaries apply to these aperture versions.</p>
<h3>{badge('X')} Across all metrics</h3><p>Magnitude, effective radius, axis ratio and mean surface brightness are mathematically linked in a Sérsic model, so their correlations with morphology are not independent effects. Human and machine labels are themselves influenced by brightness and image quality. All results are descriptive; no p-values or causal claims are made. Next tests: injected galaxies with realistic sky and PSF, alternative blank-sky estimates, mask perturbations, and matched comparisons with reference code (e.g. statmorph {cite('rodriguezgomez2019')}).</p></section>'''

    flags=pd.DataFrame([{'Flag':k,'Rows':v} for k,v in summary['quality_flags'].items()])
    dictionary=pd.read_csv(root/'paper_v2_data/column_dictionary.csv')
    appendix=f'''<section id="appendix"><h2>Appendices</h2>
<details><summary>A. Column definitions and units ({len(dictionary)} entries)</summary><input class="search" data-search="dictionary" aria-label="Filter column definitions" placeholder="Filter column names or meanings…"><div id="dictionary">{table_html(dictionary)}</div></details>
<details><summary>B. Summary statistics of every metric</summary>{table_html(stats)}</details>
<details><summary>C. Band-to-band comparisons of the same galaxies</summary>{table_html(bandcmp)}{table_html(paircmp)}</details>
<details><summary>D. Quality flags</summary>{table_html(flags)}<p>Flag totals overlap. All 96 error rows are missing cutouts. An invalid CAS status can still leave other diagnostics usable.</p></details>
<details><summary>E. Data coverage</summary>{table_html(gcov)}{table_html(ccov)}</details>
<details><summary>F. Reproduction, inputs and earlier revisions</summary><pre>python scripts/build_report_v5.py --report-dir reports/full_catalog_20261020
python scripts/validate_report_v5.py --report-dir reports/full_catalog_20261020</pre><p>The builder refuses to overwrite <code>report_v5.html</code>. Earlier revisions remain unchanged: <a href="report_v1.html">v1</a>, <a href="report_v2.html">v2</a>, <a href="report_v3.html">v3</a>, <a href="report_v4.html">v4</a>. Downloads: <a href="measurements_metrics.parquet">morphology table</a>, <a href="paper_v5_data/figure_index.json">figure index</a>, <a href="paper_v5_data/plot_range_audit.csv">all display windows and tail counts</a>, <a href="paper_v2_data/galfit_correlations.csv">GALFIT correlations</a>, <a href="paper_v3_data/associations.csv">classification associations</a>.</p></details></section>'''

    nav=[('summary','Summary'),('howto','1 Reading guide'),('glossary','2 Concepts'),('aperture','3 Aperture'),('concentration','4 C'),('asymmetry','5 A'),('clumpiness','6 S'),('diagnostics','7 Gini/M20/shape'),('overview','8 All metrics'),('embedding-atlas','9 Embedding atlas'),('discussion','10 Discussion'),('appendix','Appendices'),('references','References')]
    css='''html{scroll-behavior:smooth;scroll-padding-top:4rem}body{font:16px/1.65 system-ui,sans-serif;color:#25303a;background:#fff;margin:0}header,main{max-width:1280px;margin:auto;padding:1.5rem 2rem}
h1{font-family:Georgia,serif;font-size:2.2rem;line-height:1.2;max-width:950px}h2{font-family:Georgia,serif;font-size:1.8rem;border-bottom:1px solid #d6dce0;padding-bottom:.4rem;margin-top:3.2rem;display:flex;flex-wrap:wrap;align-items:center;gap:.6rem}h3{margin-top:2.2rem;display:flex;flex-wrap:wrap;align-items:center;gap:.5rem}p,li{max-width:105ch}
nav{position:sticky;top:0;z-index:5;background:#fff;border-block:1px solid #d6dce0;padding:.55rem 1rem;display:flex;gap:1.1rem;overflow-x:auto;white-space:nowrap;font-size:.85rem}a{color:inherit;text-decoration:underline;text-decoration-color:#4a96ba;text-underline-offset:3px}.figref{font-weight:600}
.subtitle,.meta{color:#566675}.badge{display:inline-flex;align-items:center;gap:.35em;font:600 .78rem/1 system-ui,sans-serif;letter-spacing:.02em;padding:.32em .6em;border-radius:999px;background:color-mix(in srgb,var(--fam) 14%,white);color:color-mix(in srgb,var(--fam) 75%,black);border:1px solid color-mix(in srgb,var(--fam) 45%,white);white-space:nowrap}.badge b{background:var(--fam);color:#fff;border-radius:999px;padding:.12em .45em}
section.family{border-left:5px solid var(--fam);padding-left:1.2rem;margin-left:-1.45rem}.lead{font-size:1.04rem}
.summary-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:1rem;margin:1.2rem 0}.card{border:1px solid #d6dce0;border-top:4px solid var(--fam);border-radius:8px;padding:1rem;background:#fafbfc}.card p{margin:.6rem 0}
figure.fig{margin:2rem 0 2.8rem;border-left:4px solid var(--fam);padding:.2rem 0 .2rem 1rem}.fig-head{display:flex;flex-wrap:wrap;align-items:center;gap:.6rem}.fig-head h4{margin:.4rem 0;font-size:1.05rem}
.panel-grid{display:grid;gap:14px;margin:.8rem 0}.cols-1{grid-template-columns:minmax(0,1fr);max-width:900px}.cols-2{grid-template-columns:repeat(2,minmax(0,1fr))}.cols-3{grid-template-columns:repeat(3,minmax(0,1fr))}
.panel{margin:0;scroll-margin-top:4rem}.panel img{display:block;width:100%;height:auto;border:1px solid #e1e5e8;border-radius:4px;background:#fff}.panel-cap{font-size:.84rem;color:#4c5863;margin:.3rem 0 0}
figcaption{font-size:.92rem;color:#3d4a55;background:#f5f7f9;border-radius:6px;padding:.7rem 1rem;margin-top:.6rem}figcaption p{margin:.35rem 0;max-width:none}
.equation{max-width:100%;height:auto;margin:.7rem 0}.table-scroll{overflow:auto;margin:1rem 0}table{border-collapse:collapse;font-size:.85rem;width:100%}th,td{text-align:left;vertical-align:top;padding:.45rem .65rem;border-bottom:1px solid #d6dce0}th{background:#f1f5f7}
pre{white-space:pre-wrap;overflow-wrap:anywhere;background:#f1f5f7;padding:1rem;font-size:.82rem;border-radius:6px}details{border:1px solid #d6dce0;padding:.6rem 1rem;border-radius:5px;margin:.8rem 0}summary{cursor:pointer;font-weight:600}details.audit,details.more-panels{font-size:.9rem}
.search{font:inherit;padding:.5rem;width:min(90%,650px);margin:.7rem 0}select{max-width:100%;padding:.55em;font:inherit;border:1px solid #aab3ba;border-radius:6px;background:#fff}.atlas-controls{display:flex;flex-wrap:wrap;gap:1rem;margin:1rem 0}.atlas-controls label{display:flex;flex-direction:column;font-size:.85rem;gap:.25rem}
.atlas-family{border-left:4px solid var(--fam);padding-left:1rem;margin:1.2rem 0}.refs li{margin:.35rem 0}
@media(max-width:1000px){.cols-3{grid-template-columns:repeat(2,minmax(0,1fr))}}
@media(max-width:700px){header,main{padding:1rem}.summary-grid,.cols-2,.cols-3{grid-template-columns:minmax(0,1fr)}h1{font-size:1.65rem}h2{font-size:1.4rem}section.family{margin-left:0;padding-left:.7rem;border-left-width:3px}figure.fig{padding-left:.6rem}}
@media print{nav{position:static}figure{break-inside:avoid}details{display:block}}'''
    atlas_js='const ATLAS='+json.dumps(dict(spaces=spaces,entries=emb_index,metrics=EMBED_METRICS,plain={m:PLAIN[m] for m in EMBED_METRICS}))+''';
(()=>{const s=document.getElementById('atlas-space'),b=document.getElementById('atlas-band'),cap=document.getElementById('atlas-caption');
function draw(){const space=ATLAS.spaces[Number(s.value)],idx=space.bands[b.value],e=ATLAS.entries[idx],pad=String(idx).padStart(2,'0');
cap.textContent=space.title+' · '+space.projection+' · '+b.value+' band · '+e.rows.toLocaleString()+' saved galaxies, '+e.matched.toLocaleString()+' matched to the morphology table.';
for(const m of ATLAS.metrics){const src='paper_v5_assets/emb_'+pad+'_'+m+'.png';document.querySelectorAll('[data-metric="'+m+'"]').forEach(el=>{if(el.tagName==='IMG')el.src=src;else el.href=src;});
const c=document.querySelector('[data-metric-cap="'+m+'"]');if(c)c.textContent=ATLAS.plain[m]+' — '+e.finite[m].toLocaleString()+' finite values';}}
s.addEventListener('change',draw);b.addEventListener('change',draw);draw();})();'''
    other_js='''document.querySelectorAll('[data-search]').forEach(input=>input.addEventListener('input',()=>{const q=input.value.toLowerCase();document.querySelectorAll('#'+input.dataset.search+' tbody tr').forEach(row=>row.hidden=!row.textContent.toLowerCase().includes(q));}));
window.addEventListener('message',event=>{const f=document.getElementById('association-explorer-frame');if(f&&event.source===f.contentWindow&&event.data?.type==='morphology-explorer-height'&&Number.isFinite(event.data.height)){f.style.height=Math.max(500,Math.min(4000,event.data.height+16))+'px';}});
function openTarget(){const h=location.hash.slice(1);const el=h&&document.getElementById(h);if(el){const d=el.closest('details');if(d)d.open=true;}}window.addEventListener('hashchange',openTarget);openTarget();
document.addEventListener('click',e=>{const a=e.target.closest('a.figref');if(a){const el=document.getElementById(a.getAttribute('href').slice(1));const d=el&&el.closest('details');if(d)d.open=true;}});'''
    body=executive+howto+glossary+methods+sec_aperture+sec_c+sec_a+sec_s+sec_d+sec_x+sec_e+discussion+methods+appendix+'<section id="references"><h2>References</h2>[[references]]</section>'
    html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DES LSB galaxy morphology: concentration, asymmetry and clumpiness</title><style>{css}</style></head><body>
<header><p class="meta">DES low-surface-brightness candidate catalog · scientific report · revision 5</p><h1>Concentration, asymmetry and clumpiness of 84,200 DES low-surface-brightness galaxy candidates</h1>
<p class="subtitle">One section per CAS index, each compared with fitted structure, human and machine classifications and clustering maps; supplementary diagnostics kept separate.</p><p class="meta">Measurements frozen 20 October 2026 · report revision 5 built {datetime.now(timezone.utc):%d %B %Y} · g, r and i bands</p></header>
<nav aria-label="Report sections">{''.join(f'<a href="#{k}">{t}</a>' for k,t in nav)}</nav><main>{body}<p><a href="#summary">Back to summary ↑</a></p></main>
<script>{atlas_js}
{other_js}</script></body></html>'''
    html,order,cites=resolve(html,R)
    # Every figure must be referenced from running text outside its own figure block.
    text_only=re.sub(r'<figure class="fig".*?</figure>','',html,flags=re.S)
    unref=[k for k in order if f'href="#fig-{k}"' not in text_only and f'href="#fig-{k}-' not in text_only]
    if unref:raise ValueError(f'Figures not referenced from the text: {unref}')
    out=root/'report_v5.html';out.write_text(html)
    pd.DataFrame(P.records).to_csv(data/'plot_range_audit.csv',index=False)
    figure_index=[dict(number=i+1,key=k,family=R.items[k]['family'],title=R.items[k]['title'],
                       panels=[Path(p).relative_to(root).as_posix() for p,_ in R.items[k]['panels']]+[Path(p).relative_to(root).as_posix() for _,ps,_ in R.items[k]['extra'] for p,_ in ps]) for i,k in enumerate(order)]
    (data/'figure_index.json').write_text(json.dumps(dict(figures=figure_index,embedding_maps=emb_index,default_embedding=default_space),indent=2))
    manifest=dict(report_version=VERSION,built_utc=datetime.now(timezone.utc).isoformat(),html_sha256=sha256_file(out),
                  presentation_only=True,inputs_sha256=inputs,embedding_exports_sha256={e['export']:e['export_sha256'] for e in emb},
                  figures=len(order),panels=sum(len(f['panels']) for f in figure_index),references=cites,
                  display_policy=dict(png='256-colour median-cut palette, no dither',quantiles=[.005,.995],padding=.035,histogram_bins=60,embedding_hexbin_gridsize=70,embedding_min_cell_count=3,outliers_removed=False),
                  band_palette=BAND_COLORS,previous_reports={f'report_v{i}.html':sha256_file(root/f'report_v{i}.html') for i in range(1,5)},
                  files_sha256={p.relative_to(root).as_posix():sha256_file(p) for d in [assets,data] for p in sorted(d.iterdir()) if p.is_file()},
                  source_code_sha256={str(Path(f).resolve()):sha256_file(f) for f in [__file__,'algorithmic_morphology/report_v5_plots.py','algorithmic_morphology/paper_plots.py','algorithmic_morphology/report_data.py']})
    (root/'report_v5_manifest.json').write_text(json.dumps(manifest,indent=2))
    LOG.info('Wrote %s: %d figures, %d panels (+%d atlas maps), %d references',out,len(order),manifest['panels'],len(tasks),len(cites))


if __name__=='__main__':main()
