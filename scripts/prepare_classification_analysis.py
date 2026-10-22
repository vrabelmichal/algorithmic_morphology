"""Freeze a joined, non-remeasured morphology/human/DECaLS10 analysis input."""
import argparse,json,logging
from datetime import datetime,timezone
from pathlib import Path
import pandas as pd
import numpy as np
from algorithmic_morphology.classification_data import human_features,probabilities,join_labels,association_table
from algorithmic_morphology.paper_plots import METRICS,LABELS,L,K
from algorithmic_morphology.report_data import robust_limits
from algorithmic_morphology.report import sha256_file


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--sources',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    root=Path('reports/full_catalog_20261020');paths=json.loads(a.sources.read_text())
    sources={k:dict(path=str(Path(v).resolve()),sha256=sha256_file(v)) for k,v in paths.items()}
    human=human_features(*[pd.read_parquet(paths[k]) for k in ['human_mixtures','human_records','human_membership']])
    z=probabilities(pd.read_parquet(paths['zoobot_predictions']))
    zm=json.loads(Path(paths['zoobot_metadata']).read_text())
    assert zm['checkpoint_spec']=='mwalmsley/zoobot-finetuned-decals10' and not zm['checkpoint_config_encoder_mismatch']
    morph=pd.read_parquet(root/'measurements_metrics.parquet');gf=pd.read_parquet(root/'paper_v2_data/galfit_join.parquet')
    gfcols=['coadd_object_id','band',*[c for c in gf if c.startswith('gf_')]]
    morph=morph.merge(gf[gfcols],on=['coadd_object_id','band'],validate='one_to_one',how='left',sort=False)
    joined=join_labels(morph,human,z);joined.to_parquet(a.output/'joined.parquet',index=False)
    metrics=[*METRICS,'axis_ratio_positive_aperture_v1']
    limits={m:list(robust_limits([joined[m]])) for m in metrics};limits['axis_ratio_positive_aperture_v1']=[0,1]
    for suffix in ['raw','corrected']:
        lo,hi=robust_limits([joined[L+'_'+suffix],joined[K+'_'+suffix]/10]);limits[L+'_'+suffix]=[lo,hi];limits[K+'_'+suffix]=[lo*10,hi*10]
    contract=dict(version='algorithmic_morphology_overlay_v1',metrics=metrics,labels={m:LABELS[m] for m in metrics},limits=limits,
                  table=str((a.output/'joined.parquet').resolve()),table_sha256=sha256_file(a.output/'joined.parquet'),
                  range_policy='Pooled g/r/i 0.5–99.5% with 3.5% padding; corresponding Conselice S axes x10; q [0,1]. Values unchanged.')
    (a.output/'feature_contract.json').write_text(json.dumps(contract,indent=2))
    association_table(joined,metrics).to_csv(a.output/'associations.csv',index=False)
    coverage=[]
    for band,d in joined.groupby('band'):
        coverage.append(dict(band=band,requested=len(d),human_any=int(d.human_vote_count.gt(0).sum()),human_two_plus=int(d.human_vote_count.ge(2).sum()),
                             human_three_plus=int(d.human_vote_count.ge(3).sum()),zoobot_valid=int(d.zoobot_valid.eq(True).sum()),
                             matched_two_plus=int((d.human_vote_count.ge(2)&d.zoobot_valid.eq(True)).sum())))
    pd.DataFrame(coverage).to_csv(a.output/'coverage.csv',index=False)
    # Record all available DECaLS10 candidate metadata rather than infer dates from directory names alone.
    candidates=[]
    for meta in sorted(Path(paths['zoobot_predictions']).parents[1].glob('*/classification_run_metadata.json')):
        md=json.loads(meta.read_text());pred=meta.with_name('predictions_zoobot_finetuned_decals10.parquet')
        candidates.append(dict(run=meta.parent.name,metadata=str(meta),prediction_exists=pred.exists(),objects=md.get('n_objects'),
                               checkpoint=md.get('checkpoint_spec'),encoder_mismatch=md.get('checkpoint_config_encoder_mismatch'),
                               selected=pred.resolve()==Path(paths['zoobot_predictions']).resolve(),metadata_sha256=sha256_file(meta)))
    manifest=dict(built_utc=datetime.now(timezone.utc).isoformat(),sources=sources,zoobot_candidates=candidates,zoobot_metadata=zm,
                  rows=len(joined),objects=joined.coadd_object_id.nunique(),zoobot_input_rows=len(z),zoobot_invalid=int((~z.zoobot_valid).sum()),
                  morphology_sha256=sha256_file(root/'measurements_metrics.parquet'),feature_contract=contract,
                  selection='Latest completed human runs; newest compatible Nano DECaLS10 bgsub/masked obj302317964 run, full predictions not curated subset.',
                  older_mixture_zoobot_ignored=True,source_code={str(Path(f).resolve()):sha256_file(f) for f in [__file__,'algorithmic_morphology/classification_data.py']})
    (a.output/'manifest.json').write_text(json.dumps(manifest,indent=2));print(json.dumps(coverage,indent=2),flush=True)


if __name__=='__main__':main()
