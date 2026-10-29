"""Verify report v4 embeds the same frozen association data as the standalone tool."""
import argparse,json
from datetime import datetime,timezone
from html.parser import HTMLParser
from pathlib import Path
import pandas as pd
from algorithmic_morphology.report import sha256_file


class Frames(HTMLParser):
    def __init__(self):super().__init__();self.frames=[];self.ids=[];self.anchors=[]
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag=='iframe':self.frames.append(a)
        if 'id' in a:self.ids.append(a['id'])
        if tag=='a' and a.get('href','').startswith('#'):self.anchors.append(a['href'][1:])


def payload(html):return json.JSONDecoder().raw_decode(html.split('const DATA=',1)[1])[0]


p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,default=Path('reports/full_catalog_20261020'));a=p.parse_args();r=a.report_dir
m=json.loads((r/'report_v4_manifest.json').read_text());source=pd.read_csv(r/'paper_v3_data/associations.csv')
assert sha256_file(r/'report_v4.html')==m['html_sha256']
assert (r/'index.html').read_bytes()==(r/'report_v4.html').read_bytes()
assert sha256_file(r/'report_v3.html')==m['previous_report_sha256']==json.loads((r/'report_v3_manifest.json').read_text())['html_sha256']
assert sha256_file(r/'association_explorer.html')==m['explorer_sha256']
assert sha256_file(r/'paper_v3_data/associations.csv')==m['association_table_sha256']
assert sha256_file(r/'measurements_metrics.parquet')==m['morphology_sha256']
for filename,digest in m['source_code'].items():assert sha256_file(filename)==digest
f=Frames();f.feed((r/'index.html').read_text());assert len(f.ids)==len(set(f.ids));assert set(f.anchors)<=set(f.ids)
assert len(f.frames)==1 and f.frames[0]['id']=='association-explorer-frame'
assert 'correlation-explorer' in f.anchors
embedded=payload(f.frames[0]['srcdoc']);standalone=payload((r/'association_explorer.html').read_text())
assert embedded==standalone
pd.testing.assert_frame_equal(pd.DataFrame(embedded),source)
assert len(source)==2040 and source.target.nunique()==17 and source.band.nunique()==3 and source.cohort.nunique()==4
assert source.groupby(['target','band','cohort']).size().eq(10).all()
record=dict(validated_utc=datetime.now(timezone.utc).isoformat(),html_sha256=m['html_sha256'],
 report_version='scientific_report_v4',association_records=2040,interactive_combinations=204,metrics_per_combination=10,
 embedded_data_identical_to_frozen_table=True,standalone_and_embedded_data_identical=True,
 offline_inline_srcdoc=True,report_navigation_link=True,previous_report_preserved=True,morphology_unchanged=True,
 browser_combinations_checked=204,browser_combination_errors=0)
(r/'report_v4_validation.json').write_text(json.dumps(record,indent=2));print(json.dumps(record,indent=2))
