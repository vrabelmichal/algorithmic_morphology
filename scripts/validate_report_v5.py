"""Audit report revision 5: hashes, images, links, figure references and key numbers."""
import argparse
from datetime import datetime,timezone
from html.parser import HTMLParser
import json
from pathlib import Path
import re

import pandas as pd
from PIL import Image

from algorithmic_morphology.report import sha256_file


class Audit(HTMLParser):
    def __init__(self):
        super().__init__();self.ids=set();self.anchors=[];self.images=[];self.links=[]

    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'id' in a:self.ids.add(a['id'])
        href=a.get('href') or ''
        if tag=='a' and href.startswith('#'):self.anchors.append(href[1:])
        elif tag=='a' and href and not href.startswith(('http://','https://')):self.links.append(href)
        if tag=='img' and not a['src'].startswith('data:'):self.images.append(a['src'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,default=Path('reports/full_catalog_20261020'));a=p.parse_args()
    root=a.report_dir.resolve();meta=json.loads((root/'report_v5_manifest.json').read_text());html=(root/'report_v5.html').read_text()
    checks={}
    assert sha256_file(root/'report_v5.html')==meta['html_sha256'];checks['html_sha256']=meta['html_sha256']
    for name,digest in meta['inputs_sha256'].items():assert sha256_file(root/name)==digest,name
    checks['frozen_inputs_unchanged']=len(meta['inputs_sha256'])
    for name,digest in meta['files_sha256'].items():assert sha256_file(root/name)==digest,name
    checks['asset_and_data_files']=len(meta['files_sha256'])
    for name,digest in meta['previous_reports'].items():assert sha256_file(root/name)==digest,name
    history=json.loads((root/'report_history.json').read_text())
    for v in ['v1','v2','v3','v4']:assert sha256_file(root/history[v]['html'])==history[v]['sha256'],v
    checks['previous_revisions_unchanged']=sorted(meta['previous_reports'])
    assert '[[' not in html and ']]' not in html,'Unresolved placeholder'
    parser=Audit();parser.feed(html)
    missing=[x for x in parser.anchors if x not in parser.ids];assert not missing,missing
    checks['internal_anchors']=len(parser.anchors)
    for src in set(parser.images):
        with Image.open(root/src) as im:im.verify()
    for href in set(parser.links):assert (root/href.split('#')[0]).exists(),href
    checks['images_verified']=len(set(parser.images));checks['local_links']=len(set(parser.links))
    index=json.loads((root/'paper_v5_data/figure_index.json').read_text())
    figures=index['figures'];assert len(figures)==meta['figures']
    # Each figure: rendered once, numbered in order, and referenced from running text outside itself.
    text=re.sub(r'<figure class="fig".*?</figure>','',html,flags=re.S)
    for f in figures:
        assert html.count(f'id="fig-{f["key"]}"')==1,f['key']
        assert f'<h4>Figure {f["number"]}. ' in html
        assert re.search(rf'href="#fig-{re.escape(f["key"])}(-[a-z])?"',text),'unreferenced '+f['key']
        for panel in f['panels']:
            with Image.open(root/panel) as im:im.verify()
    checks['figures_referenced_from_text']=len(figures);checks['figure_panels']=sum(len(f['panels']) for f in figures)
    for e in index['embedding_maps']:
        for m in e['finite']:
            with Image.open(root/f'paper_v5_assets/emb_{e["index"]:02d}_{m}.png') as im:im.verify()
    checks['embedding_maps']=sum(len(e['finite']) for e in index['embedding_maps'])
    # Every cited reference has a numbered entry, and every glossary target exists.
    refs=re.findall(r'<li id="ref-(\w+)">',html);assert refs==meta['references']
    for k in refs:assert f'href="#ref-{k}"' in text,k
    checks['references']=len(refs)
    # Spot-check numbers in prose against independent source tables.
    corr=pd.read_csv(root/'paper_v2_data/galfit_correlations.csv');corr=corr[corr.cohort=='gf_unflagged'].set_index(['metric','property','band']).spearman
    stats=pd.read_csv(root/'statistics.csv').set_index(['metric','band'])
    expected=['/'.join(f'{corr[(m,"gf_n" if m=="C" else "gf_mu_mean",b)]:.2f}' for b in 'gri') for m in ['C','A_raw','A_corrected']]
    neg=stats.loc['A_corrected'];neg=(100*neg.negative/neg.valid)
    expected+=[f'{neg.min():.0f}–{neg.max():.0f}%',f'{stats.loc[("C","r"),"median"]:.3g}']
    for s in expected:assert s in html,s
    checks['prose_numbers_checked']=expected
    out=dict(validated_utc=datetime.now(timezone.utc).isoformat(),report='report_v5.html',checks=checks,status='passed')
    (root/'report_v5_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))


if __name__=='__main__':main()
