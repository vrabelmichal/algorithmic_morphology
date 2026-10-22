import numpy as np
import pandas as pd
import pytest
from algorithmic_morphology.report_data import robust_limits,join_galfit,GF_COLUMNS,BANDS


def test_limits_resist_outlier_without_deleting_values():
    values=np.r_[np.linspace(-.2,.3,10000),1e9,np.nan,np.inf]
    before=values.copy();lo,hi=robust_limits([values])
    assert -.3<lo<0 and .2<hi<.4
    np.testing.assert_equal(values,before)


def test_constant_positive_and_empty_limits():
    lo,hi=robust_limits([[2]*100]);assert lo<2<hi
    lo,hi=robust_limits([[0,1,10,np.inf]],positive=True);assert 0<lo<hi
    assert robust_limits([[np.nan]])==(0.,1.)


def fixture():
    f=pd.DataFrame({'coadd_object_id':[123,456],'success':[True,True]})
    for band in BANDS:
        suffix='' if band=='g' else '_'+band
        for dest,source in GF_COLUMNS.items():f[source+suffix]=[2.,3.]
        f['sersic_q'+suffix]=[.7,.8]
        f['sersic_mag'+suffix]=[dict(g=20,r=19,i=18)[band],21]
        for key in ['x','y','mag','reff_pixels','n','q']:f['sersic_'+key+'_warning'+suffix]=0
    return f


def test_join_suffixes_coverage_and_warnings():
    morph=pd.DataFrame({'coadd_object_id':['123','123','123','999'],'band':['g','r','i','g']})
    table=fixture();table.loc[0,'sersic_mag_warning_r']=1
    out=join_galfit(morph,table)
    assert list(out.gf_mag[:3])==[20,19,18]
    assert len(out)==4 and not out.gf_present.iloc[3]
    assert list(out.gf_unflagged[:3])==[True,False,True]
    assert table.sersic_mag_r.iloc[0]==19


def test_duplicate_ids_fail_instead_of_cartesian_join():
    table=fixture();table.loc[1,'coadd_object_id']=123
    with pytest.raises(ValueError,match='Duplicate GALFIT'):
        join_galfit(pd.DataFrame({'coadd_object_id':['123'],'band':['g']}),table)
