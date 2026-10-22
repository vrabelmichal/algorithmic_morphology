"""Audit v3 data, plots and embedding exports against independent source tables."""
import argparse,base64
from datetime import datetime,timezone
from html.parser import HTMLParser
from io import BytesIO
import json
from pathlib import Path
from PIL import Image
import numpy as np
import pandas as pd
from algorithmic_morphology.report import sha256_file
from algorithmic_morphology.classification_data import HUMAN_LABELS


class Audit(HTMLParser):
    def __init__(self,root):
        super().__init__();self.root=root;self.ids=[];self.anchors=[];self.images=0
    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if 'id' in a:self.ids.append(a['id'])
        if tag=='a' and a.get('href','').startswith('#'):self.anchors.append(a['href'][1:])
        if tag=='img':
            src=a['src'];self.images+=1
            if src.startswith('data:image/png'):Image.open(BytesIO(base64.b64decode(src.split(',',1)[1]))).verify()
            elif src.startswith('data:image/svg'):assert b'<svg' in base64.b64decode(src.split(',',1)[1])
            else:Image.open(self.root/src).verify()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,required=True);a=p.parse_args();root=a.report_dir.resolve()
    meta=json.loads((root/'report_v3_manifest.json').read_text());source=meta['classification_source_manifest'];data=root/'paper_v3_data'
    assert sha256_file(root/'report_v3.html')==meta['html_sha256'];assert (root/'index.html').read_bytes()==(root/'report_v3.html').read_bytes()
    assert sha256_file(root/'report_v2.html')==meta['previous_report_sha256']==json.loads((root/'report_v2_manifest.json').read_text())['html_sha256']
    assert sha256_file(root/'measurements_metrics.parquet')==meta['morphology_sha256']==source['morphology_sha256']
    for name,digest in meta['files_sha256'].items():assert sha256_file(root/name)==digest,name
    for name,digest in meta['source_code'].items():assert sha256_file(name)==digest,name
    for record in source['sources'].values():assert sha256_file(record['path'])==record['sha256']
    labels=pd.read_parquet(data/'classification_labels.parquet');assert len(labels)==84200 and labels.coadd_object_id.nunique()==84200
    raw=pd.read_parquet(source['sources']['zoobot_predictions']['path']);raw.id_str=raw.id_str.astype(str);raw=raw.set_index('id_str')
    labeled=labels.set_index('coadd_object_id');assert int(labeled.zoobot_valid.eq(True).sum())==len(raw)==84104
    np.testing.assert_array_equal(labeled.loc[raw.index,[f'class_{i}' for i in range(10)]].to_numpy(),raw.to_numpy())
    sums=raw.to_numpy(float).sum(axis=1);assert np.isfinite(raw).all().all() and (raw>=0).all().all() and np.allclose(sums,1,atol=1e-4,rtol=0)
    records=pd.read_parquet(source['sources']['human_records']['path']);records.galaxyId=records.galaxyId.astype(str);records=records.set_index('galaxyId').loc[labeled.index]
    mixes=pd.read_parquet(source['sources']['human_mixtures']['path']);mixes.galaxyId=mixes.galaxyId.astype(str);mixes=mixes.set_index('galaxyId')
    for field in ['count_featureless','count_irregular_or_unsure','count_ltg','count_etg','human_vote_count']:
        np.testing.assert_array_equal(labeled.loc[mixes.index,field].to_numpy(),mixes[field].to_numpy())
    for source_col,dest in [('featurelessVotes','count_featureless'),('irregularVotes','count_irregular_or_unsure'),('ltgVotes','count_ltg'),('etgVotes','count_etg')]:
        np.testing.assert_array_equal(records.loc[mixes.index,source_col],mixes[dest])
    for col in list(HUMAN_LABELS)[:4]:
        expected=mixes[col].where(mixes.human_vote_count>=2)
        np.testing.assert_allclose(labeled.loc[mixes.index,col],expected,equal_nan=True)
    for field,num,den in [('p_lsb','lsbVotes','lsbComparableVotes'),('p_visible_nucleus','visibleNucleusVotes','visibleNucleusComparableVotes'),('p_failed_fitting','failedFittingVotes','failedFittingComparableVotes')]:
        expected=(records[num]/records[den]).where(records[den]>=2)
        np.testing.assert_allclose(labeled[field],expected,equal_nan=True)
    memberships=pd.read_parquet(source['sources']['human_membership']['path']);memberships.galaxyId=memberships.galaxyId.astype(str);memberships=memberships.set_index('galaxyId').loc[labeled.index]
    for col in memberships:np.testing.assert_array_equal(labeled[col],memberships[col])
    morph=pd.read_parquet(root/'measurements_metrics.parquet');morph.coadd_object_id=morph.coadd_object_id.astype(str)
    joined=morph.merge(labels,on='coadd_object_id',validate='many_to_one',how='left',sort=False)
    assert len(joined)==252600 and not joined.duplicated(['coadd_object_id','band']).any()
    assoc=pd.read_csv(data/'associations.csv');metrics=source['feature_contract']['metrics']
    for row in assoc[assoc.cohort.eq('available')].itertuples():
        pair=joined.loc[joined.band.eq(row.band),[row.metric,row.target]].replace([np.inf,-np.inf],np.nan).dropna()
        assert len(pair)==row.n
        expected=pair[row.metric].corr(pair[row.target],method='spearman')
        np.testing.assert_allclose(row.spearman,expected,equal_nan=True,rtol=1e-12,atol=1e-12)
    ranges=pd.read_csv(data/'class_distribution_range_audit.csv')
    for row in ranges.itertuples():
        vals=joined.loc[joined.band.eq(row.figure[-1])&joined[row.selection].eq(True),row.metric].to_numpy(float);vals=vals[np.isfinite(vals)]
        assert len(vals)==row.finite and (vals<row.lower).sum()==row.below and (vals>row.upper).sum()==row.above
    embedding=json.loads((data/'clustering_overlay_manifest.json').read_text());assert len(embedding['embeddings'])==63
    for key in ['source_code','export_source_code']:
        for filename,digest in embedding.get(key,{}).items():assert sha256_file(filename)==digest
    unique=set();rows=[]
    for item in embedding['embeddings']:
        path=Path(item['embedding']);assert sha256_file(path)==item['embedding_sha256'];unique.add(str(path))
        saved=pd.read_parquet(path);assert saved.index.name=='coadd_object_id' and not saved.index.has_duplicates
        exported=pd.read_csv(item['export'],dtype={'coadd_object_id':str});assert sha256_file(item['export'])==item['export_sha256']
        if 'export_config' in item:
            assert sha256_file(item['export_config'])==item['export_config_sha256']
            config=json.loads(Path(item['export_config']).read_text())
            assert config['colorMapping']['column']=='C' and set(config['continuousColumns'])==set(metrics)
            assert exported.embedding_shape.eq('circle').all() and exported.color.str.match(r'^#[0-9a-f]{6}$').all()
            assert exported.loc[~np.isfinite(exported.C),'color'].eq('#c4c4c4').all()
        assert exported.coadd_object_id.tolist()==saved.index.astype(str).tolist()
        actual=exported[['embedding_x','embedding_y']].to_numpy();expected=saved[['embedding_0','embedding_1']].to_numpy()
        if expected.dtype==np.float32:
            # pandas emits shortest float32 round-trip decimal strings in CSV.
            np.testing.assert_array_equal(actual.astype(np.float32),expected)
        else:np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=1e-13)
        original=morph[morph.band.eq(item['band'])].set_index('coadd_object_id').reindex(exported.coadd_object_id)
        assert original.status.notna().sum()==item['matched']==item['rows']==len(saved)
        for metric in metrics:
            np.testing.assert_allclose(exported[metric],original[metric],equal_nan=True,rtol=1e-13,atol=1e-13)
        assert sha256_file(item['png'])==item['png_sha256'];Image.open(item['png']).verify()
        for audit in item['metrics']:
            vals=exported[audit['metric']].to_numpy(float);finite=np.isfinite(vals)
            assert finite.sum()==audit['finite'] and (~finite).sum()==audit['missing']
            assert (vals[finite]<audit['vmin']).sum()==audit['below'] and (vals[finite]>audit['vmax']).sum()==audit['above']
        rows.append({k:item[k] for k in ['title','projection','band','rows','matched']})
    assert len(unique)==21
    pd.DataFrame(rows).to_csv(data/'embedding_coverage_validation.csv',index=False)
    # This validation product is written after manifest creation and has its own digest below.
    audit=Audit(root);audit.feed((root/'index.html').read_text());assert len(audit.ids)==len(set(audit.ids));assert set(audit.anchors)<=set(audit.ids)
    figures=json.loads((data/'figure_index.json').read_text());assert len(figures['classification'])==25 and len(figures['embeddings'])==63
    for item in [*figures['classification'],*figures['embeddings']]:Image.open(root/item['file']).verify()
    definitions=pd.read_csv(data/'classification_column_dictionary.csv');assert set(labels)-{'coadd_object_id'}<=set(definitions.column)
    history=json.loads((root/'report_history.json').read_text());history['v3']={'html':'report_v3.html','sha256':meta['html_sha256'],'validation':'report_v3_validation.json'};history['default']='index.html'
    (root/'report_history.json').write_text(json.dumps(history,indent=2))
    report=dict(validated_utc=datetime.now(timezone.utc).isoformat(),report_version='scientific_report_v3',html_sha256=meta['html_sha256'],
        morphology_rows=len(morph),unique_objects=len(labels),human_objects=len(mixes),human_two_plus=int(labeled.human_vote_count.ge(2).sum()),
        zoobot_vectors=len(raw),probabilities_identical_to_source=True,human_counts_equal_between_latest_producers=True,
        human_membership_unchanged=True,all_available_pair_correlations_independently_recomputed=True,
        classification_figures=25,embedding_spaces=len(unique),band_panels=63,embedding_values_and_coordinates_verified=True,
        plot_tails_independently_verified=True,all_pngs_valid=True,internal_navigation_resolves=True,
        morphology_unchanged=True,previous_report_preserved=True,
        coverage_validation_sha256=sha256_file(data/'embedding_coverage_validation.csv'))
    (root/'report_v3_validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))


if __name__=='__main__':main()
