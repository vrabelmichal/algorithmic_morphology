"""Validated GALFIT joins and transparent plotting limits for report revision 2."""
from pathlib import Path
import numpy as np
import pandas as pd

BANDS=('g','r','i')
GF_COLUMNS={
 'gf_mag':'sersic_mag','gf_re_arcsec':'sersic_reff_arcsec',
 'gf_re_pixels':'sersic_reff_pixels','gf_mu_mean':'sersic_mean_mue',
 'gf_mu_e':'sersic_mue','gf_n':'sersic_n','gf_q':'sersic_q',
 'gf_chi2red':'sersic_chi2red',
}
WARNING_PARAMETERS=('x','y','mag','reff_pixels','n','q')


def finite(values):
    return pd.to_numeric(pd.Series(values),errors='coerce').replace([np.inf,-np.inf],np.nan).to_numpy(dtype=float)


def robust_limits(arrays,quantiles=(.005,.995),positive=False,padding=.035):
    values=np.concatenate([finite(a) for a in arrays])
    values=values[np.isfinite(values)]
    if positive:values=values[values>0]
    if not len(values):return (1.,10.) if positive else (0.,1.)
    lo,hi=np.quantile(values,quantiles)
    if lo==hi:
        spread=max(abs(lo)*.05,.05);lo-=spread;hi+=spread
    if positive:
        lo=max(lo,np.nextafter(0.,1.));ratio=(hi/lo)**padding
        return float(lo/ratio),float(hi*ratio)
    spread=(hi-lo)*padding
    return float(lo-spread),float(hi+spread)


def id_strings(series):
    if series.isna().any():raise ValueError('Missing GALFIT object ID')
    # Preserve large integer identifiers; never round floating IDs silently.
    if pd.api.types.is_float_dtype(series):
        if not np.equal(series,np.floor(series)).all():raise ValueError('Fractional object ID')
        series=series.astype('int64')
    return series.astype(str)


def galfit_long(table):
    table=table.copy();table['coadd_object_id']=id_strings(table.coadd_object_id)
    if table.coadd_object_id.duplicated().any():raise ValueError('Duplicate GALFIT object IDs')
    parts=[]
    for band in BANDS:
        suffix='' if band=='g' else '_'+band
        part=pd.DataFrame({'coadd_object_id':table.coadd_object_id,'band':band})
        for dest,source in GF_COLUMNS.items():part[dest]=finite(table[source+suffix])
        part['gf_success']=table.success.eq(True).to_numpy()
        warning_cols=['sersic_'+k+'_warning'+suffix for k in WARNING_PARAMETERS]
        part['gf_warning_count']=table[warning_cols].ne(0).sum(axis=1).to_numpy()
        core=part[['gf_mag','gf_re_arcsec','gf_n','gf_q','gf_mu_mean','gf_mu_e','gf_chi2red']]
        part['gf_finite_physical']=(np.isfinite(core).all(axis=1)&(part.gf_re_arcsec>0)&(part.gf_n>0)&(part.gf_q>0)&(part.gf_q<=1)&(part.gf_chi2red>0))
        part['gf_unflagged']=part.gf_success&part.gf_finite_physical&part.gf_warning_count.eq(0)
        part['gf_at_n_bound']=np.isclose(part.gf_n,.2,atol=.001,rtol=0)|np.isclose(part.gf_n,5.,atol=.001,rtol=0)
        part['gf_at_q_bound']=np.isclose(part.gf_q,.3,atol=.001,rtol=0)|np.isclose(part.gf_q,1.,atol=.001,rtol=0)
        part['gf_interior']=part.gf_unflagged&~part.gf_at_n_bound&~part.gf_at_q_bound
        part['gf_present']=True
        parts.append(part)
    return pd.concat(parts,ignore_index=True)


def join_galfit(morphology,galfit):
    if morphology.duplicated(['coadd_object_id','band']).any():raise ValueError('Duplicate morphology keys')
    left=morphology.copy();left['coadd_object_id']=id_strings(left.coadd_object_id)
    result=left.merge(galfit_long(galfit),how='left',on=['coadd_object_id','band'],validate='one_to_one',sort=False)
    assert len(result)==len(left)
    for key in ['gf_present','gf_success','gf_finite_physical','gf_unflagged','gf_interior','gf_at_n_bound','gf_at_q_bound']:
        result[key]=result[key].fillna(False).astype(bool)
    return result


def relation_stats(frame,metrics,properties):
    rows=[]
    for cohort in ['gf_unflagged','gf_interior']:
        for band,part in frame[frame[cohort]].groupby('band',sort=True):
            for metric in metrics:
                for prop in properties:
                    data=part[[metric,prop]].replace([np.inf,-np.inf],np.nan).dropna()
                    rho=data[metric].corr(data[prop],method='spearman') if len(data)>2 else np.nan
                    rows.append(dict(cohort=cohort,band=band,metric=metric,property=prop,n=len(data),spearman=rho))
    return pd.DataFrame(rows)
