"""Exact joins and descriptive morphology associations with external labels."""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, rankdata
from .report_data import id_strings

HUMAN_LABELS={'p_featureless':'Featureless','p_irregular_or_unsure':'Not sure (Irr/other)',
              'p_ltg':'LTG / spiral','p_etg':'ETG / elliptical',
              'p_lsb':'Is LSB: yes','p_visible_nucleus':'Visible nucleus: yes','p_failed_fitting':'Failed fitting: yes'}
CLASS_LABELS=['Disturbed','Merging','Round smooth','Intermediate round smooth','Cigar-shaped smooth',
              'Barred spiral','Unbarred tight spiral','Unbarred loose spiral','Edge-on without bulge','Edge-on with bulge']


def keyed(frame,key):
    frame=frame.copy();frame['coadd_object_id']=id_strings(frame[key])
    if frame.coadd_object_id.duplicated().any():raise ValueError(f'Duplicate IDs in {key}')
    return frame


def probabilities(frame):
    frame=keyed(frame,'id_str');columns=[f'class_{i}' for i in range(10)]
    values=frame[columns].to_numpy(dtype=float)
    valid=np.isfinite(values).all(axis=1)&(values>=0).all(axis=1)&(values<=1).all(axis=1)&np.isclose(values.sum(axis=1),1,atol=1e-4,rtol=0)
    frame['zoobot_valid']=valid
    # Preserve supplied float32 probabilities; never renormalize or fabricate.
    frame.loc[~valid,columns]=np.nan
    frame['zoobot_top']=pd.Series(pd.NA,index=frame.index,dtype='Int64')
    frame.loc[valid,'zoobot_top']=np.argmax(values[valid],axis=1)
    frame['zoobot_top_probability']=frame[columns].max(axis=1)
    return frame[['coadd_object_id',*columns,'zoobot_valid','zoobot_top','zoobot_top_probability']]


def human_features(mixtures,records,membership):
    mix=keyed(mixtures,'galaxyId');rec=keyed(records,'galaxyId');mem=keyed(membership,'galaxyId')
    columns=['coadd_object_id','human_vote_count','top_morphology','normalized_entropy',*list(HUMAN_LABELS)[:4],
             'count_featureless','count_irregular_or_unsure','count_ltg','count_etg']
    result=rec[['coadd_object_id']].merge(mix[columns],how='left',on='coadd_object_id',validate='one_to_one')
    total=result[[f'count_{k}' for k in ['featureless','irregular_or_unsure','ltg','etg']]].sum(axis=1,min_count=1)
    if not np.allclose(total.dropna(),result.loc[total.notna(),'human_vote_count']):raise ValueError('Human vote denominators disagree')
    for key,num,den in [('p_lsb','lsbVotes','lsbComparableVotes'),('p_visible_nucleus','visibleNucleusVotes','visibleNucleusComparableVotes'),
                        ('p_failed_fitting','failedFittingVotes','failedFittingComparableVotes')]:
        values=rec.set_index('coadd_object_id').reindex(result.coadd_object_id)
        n=pd.to_numeric(values[den]).to_numpy(float);v=pd.to_numeric(values[num]).to_numpy(float)
        if np.any((v<0)|(v>n)):raise ValueError('Invalid human answer counts')
        result[key]=np.divide(v,n,out=np.full(len(n),np.nan),where=n>=2)
        result[key+'_count']=n
    for col in list(HUMAN_LABELS)[:4]:result.loc[result.human_vote_count<2,col]=np.nan
    # Export published query membership unchanged; do not redefine weighted consensus.
    result=result.merge(mem.drop(columns='galaxyId'),on='coadd_object_id',how='left',validate='one_to_one')
    return result


def join_labels(morphology,human,zoobot):
    result=morphology.copy();result['coadd_object_id']=id_strings(result.coadd_object_id)
    if result.duplicated(['coadd_object_id','band']).any():raise ValueError('Duplicate morphology keys')
    for external in [human,zoobot]:
        result=result.merge(external,on='coadd_object_id',how='left',validate='many_to_one',sort=False)
    if len(result)!=len(morphology):raise ValueError('Join changed morphology row count')
    return result


def partial_rank(x,y,controls):
    """Linear residual correlation of ranks, descriptive rather than causal."""
    x=rankdata(x);y=rankdata(y)
    design=np.column_stack([np.ones(len(x)),*[rankdata(controls[:,j]) for j in range(controls.shape[1])]])
    a=x-design@np.linalg.lstsq(design,x,rcond=None)[0]
    b=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
    if np.std(a)<1e-10 or np.std(b)<1e-10:return np.nan
    return float(np.corrcoef(a,b)[0,1])


def association_table(frame,metrics):
    rows=[]
    for band,part in frame.groupby('band'):
        for target in [*HUMAN_LABELS,*[f'class_{i}' for i in range(10)]]:
            for metric in metrics:
                for cohort in ['available','human_zoobot_matched','three_plus_human_votes','galfit_adjusted']:
                    mask=part[target].notna()&np.isfinite(part[metric])
                    if cohort=='human_zoobot_matched':mask &= part.human_vote_count.ge(2)&part.zoobot_valid.eq(True)
                    if cohort=='three_plus_human_votes':mask &= part.human_vote_count.ge(3)
                    cols=[metric,target]
                    if cohort=='galfit_adjusted':
                        mask &= part.gf_unflagged.eq(True)&part.gf_re_arcsec.gt(0)
                        cols+=['gf_mu_mean','gf_re_arcsec','gf_q']
                    pair=part.loc[mask,cols].replace([np.inf,-np.inf],np.nan).dropna()
                    rho=pair[metric].corr(pair[target],method='spearman') if len(pair)>2 and pair[metric].nunique()>1 and pair[target].nunique()>1 else np.nan
                    adjusted=partial_rank(pair[metric],pair[target],pair[cols[2:]].to_numpy()) if cohort=='galfit_adjusted' and len(pair)>20 else np.nan
                    rows.append(dict(band=band,target=target,metric=metric,cohort=cohort,n=len(pair),spearman=rho,partial_rank=adjusted))
    return pd.DataFrame(rows)


def class_effect(values,selected):
    a=np.asarray(values)[np.asarray(selected,dtype=bool)];b=np.asarray(values)[~np.asarray(selected,dtype=bool)]
    a=a[np.isfinite(a)];b=b[np.isfinite(b)]
    effect=2*mannwhitneyu(a,b,alternative='two-sided').statistic/(len(a)*len(b))-1 if len(a) and len(b) else np.nan
    return len(a),len(b),effect
