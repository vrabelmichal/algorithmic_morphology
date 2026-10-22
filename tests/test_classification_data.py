import numpy as np
import pandas as pd
import pytest
from algorithmic_morphology.classification_data import probabilities,keyed,join_labels,partial_rank,class_effect


def test_probability_validation_preserves_and_rejects():
    d=pd.DataFrame(np.eye(10)[:3],columns=[f'class_{i}' for i in range(10)])
    d['id_str']=['1','2','3'];d.loc[1,'class_0']=-.1;d.loc[2,'class_2']=.9
    out=probabilities(d)
    assert out.zoobot_valid.tolist()==[True,False,False]
    assert out.loc[0,'class_0']==1 and out.loc[1,'class_0']!=out.loc[1,'class_0']


def test_exact_many_to_one_join_preserves_band_rows():
    m=pd.DataFrame({'coadd_object_id':['2','1','2'],'band':['g','g','r']})
    h=pd.DataFrame({'coadd_object_id':['1'],'p':[.25]})
    z=pd.DataFrame({'coadd_object_id':['2'],'p_z':[.9]})
    out=join_labels(m,h,z)
    assert out.coadd_object_id.tolist()==['2','1','2'] and len(out)==3
    assert out.p_z.tolist()[::2]==[.9,.9] and np.isnan(out.p.iloc[0])
    with pytest.raises(ValueError):keyed(pd.DataFrame({'id':['1','1']}),'id')


def test_rank_effect_and_adjustment():
    assert class_effect([1,2,3,4],[False,False,True,True])==(2,2,1.)
    rng=np.random.default_rng(42);c=rng.normal(size=4000);x=c+rng.normal(size=4000);y=c+rng.normal(size=4000)
    assert abs(partial_rank(x,y,c[:,None]))<.06
    assert partial_rank(x,x,c[:,None])>.9999
