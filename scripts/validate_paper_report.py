"""Validate report revision 2 without remeasuring galaxies or changing any table."""
import argparse
import base64
from datetime import datetime,timezone
from html.parser import HTMLParser
from io import BytesIO
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from algorithmic_morphology.report import sha256_file
from algorithmic_morphology.report_data import GF_COLUMNS


class HtmlAudit(HTMLParser):
    def __init__(self):
        super().__init__();self.ids=[];self.anchors=[];self.images=0
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if 'id' in attrs:self.ids.append(attrs['id'])
        if tag=='a' and attrs.get('href','').startswith('#'):self.anchors.append(attrs['href'][1:])
        if tag=='img':
            src=attrs['src'];assert src.startswith('data:image/'),'Non-self-contained image'
            payload=base64.b64decode(src.split(',',1)[1]);self.images+=1
            if src.startswith('data:image/png'):Image.open(BytesIO(payload)).verify()
            else:assert b'<svg' in payload


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,required=True)
    root=p.parse_args().report_dir.resolve();m=json.loads((root/'report_v2_manifest.json').read_text())
    assert sha256_file(root/'report_v2.html')==m['html_sha256']
    assert (root/'index.html').read_bytes()==(root/'report_v2.html').read_bytes()
    assert sha256_file(root/'measurements_metrics.parquet')==m['morphology']['compact_table_sha256']
    for f,digest in m['files_sha256'].items():assert sha256_file(root/f)==digest,f
    for f,digest in m['source_code_sha256'].items():assert sha256_file(f)==digest,f
    audit=HtmlAudit();audit.feed((root/'index.html').read_text())
    assert len(audit.ids)==len(set(audit.ids));assert set(audit.anchors)<=set(audit.ids)
    assert audit.images==27;assert len(m['figures'])==19
    frame=pd.read_parquet(root/'measurements_metrics.parquet')
    joined=pd.read_parquet(root/'paper_v2_data/galfit_join.parquet')
    assert len(frame)==len(joined)==252600
    assert not joined.duplicated(['coadd_object_id','band']).any()
    pd.testing.assert_frame_equal(frame[['coadd_object_id','band']],joined[['coadd_object_id','band']])
    raw=pd.read_parquet(m['galfit']['primary']['path']);raw['coadd_object_id']=raw.coadd_object_id.astype(str)
    assert sha256_file(m['galfit']['primary']['path'])==m['galfit']['primary']['sha256']
    assert not raw.coadd_object_id.duplicated().any();raw=raw.set_index('coadd_object_id')
    for band,part in joined.groupby('band'):
        assert len(part)==84200 and part.coadd_object_id.nunique()==84200
        assert part.gf_present.sum()==84061
        assert set(part.loc[part.gf_present,'coadd_object_id'])==set(raw.index)
        suffix='' if band=='g' else '_'+band
        present=part[part.gf_present]
        for dest,source in GF_COLUMNS.items():
            expected=pd.to_numeric(raw.loc[present.coadd_object_id,source+suffix],errors='coerce').replace([np.inf,-np.inf],np.nan).to_numpy()
            np.testing.assert_allclose(present[dest].to_numpy(),expected,equal_nan=True)
    ranges=pd.read_csv(root/'paper_v2_data/plot_range_audit.csv')
    assert (ranges.finite==ranges.shown+ranges.below+ranges.above).all();assert (ranges.upper>ranges.lower).all()
    hist=ranges[ranges.figure.str.endswith(('_zoom','_full'))]
    for r in hist.itertuples():
        col=r.series.removesuffix('/10');v=pd.to_numeric(frame.loc[frame.band==r.band,col],errors='coerce').to_numpy(dtype=float)
        if r.series.endswith('/10'):v=v/10
        v=v[np.isfinite(v)]
        assert len(v)==r.finite and (v<r.lower).sum()==r.below and (v>r.upper).sum()==r.above
    for key in ['A_corrected','S_lotz2004_statmorph_like_v1_raw','S_conselice2003_like_v1_corrected']:
        rows=ranges[(ranges.series==key)&ranges.figure.str.endswith('_zoom')]
        assert len(rows)==3 and rows.upper.max()-rows.lower.min()<100
    dictionary=pd.read_csv(root/'paper_v2_data/column_dictionary.csv')
    assert set(frame.columns)<=set(dictionary.Column)
    old_digest=sha256_file(root/'report_v1.html')
    if (root/'validation.json').exists():assert old_digest==json.loads((root/'validation.json').read_text())['html_sha256']
    record=dict(validated_utc=datetime.now(timezone.utc).isoformat(),report_version='scientific_report_v2',
                morphology_rows=len(frame),morphology_objects=frame.coadd_object_id.nunique(),morphology_hash_unchanged=True,
                source_digests_match=True,galfit_rows=len(raw),join_unique_and_complete=True,
                galfit_values_equal_to_source_with_verified_band_suffixes=True,range_audit_rows=len(ranges),
                histogram_tails_independently_verified=True,main_figures=19,embedded_images=27,
                complete_column_dictionary=True,internal_anchors_resolve=True,previous_report_preserved=True,
                html_sha256=m['html_sha256'])
    (root/'report_v2_validation.json').write_text(json.dumps(record,indent=2))
    (root/'report_history.json').write_text(json.dumps({'v1':{'html':'report_v1.html','sha256':old_digest,'validation':'validation.json'},
                                                     'v2':{'html':'report_v2.html','sha256':m['html_sha256'],'validation':'report_v2_validation.json'},
                                                     'default':'index.html','measurement_table_unchanged':True},indent=2))
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
