import numpy as np
import pytest
from algorithmic_morphology.extended import aperture_diagnostics


def test_gini_uniform_and_single_bright_pixel():
    image=np.ones((32,32));center=(15.5,15.5)
    assert aperture_diagnostics(image,None,center,10)['Gini_aperture_v1']==pytest.approx(0)
    image[:]=0;image[15,15]=1
    assert aperture_diagnostics(image,None,center,10)['Gini_aperture_v1']==pytest.approx(1)


def test_resolved_gaussian_m20_and_ellipticity():
    y,x=np.indices((192,192));image=np.exp(-((x-95.5)**2+(y-95.5)**2)/(2*12**2))
    r=aperture_diagnostics(image,None,(95.5,95.5),60)
    expected=np.log10(1-(1-np.log(.8))*.8)
    assert r['M20_aperture_v1']==pytest.approx(expected,abs=.025)
    assert r['ellipticity_positive_aperture_v1']==pytest.approx(0,abs=1e-12)
    elongated=np.exp(-((x-95.5)**2/(2*12**2)+(y-95.5)**2/(2*6**2)))
    r=aperture_diagnostics(elongated,None,(95.5,95.5),60)
    assert r['axis_ratio_positive_aperture_v1']==pytest.approx(.5,abs=.001)


def test_masked_star_and_reflection_invariance():
    y,x=np.indices((64,64));image=np.exp(-((x-31.5)**2+(y-31.5)**2)/100)
    mask=np.zeros_like(image,dtype=bool);mask[20:25,35:40]=True
    reference=aperture_diagnostics(image,mask,(31.5,31.5),25)
    image[mask]=1e8
    measured=aperture_diagnostics(image,mask,(31.5,31.5),25)
    reflected=aperture_diagnostics(image[:,::-1],mask[:,::-1],(31.5,31.5),25)
    for key in ['Gini_aperture_v1','M20_aperture_v1','axis_ratio_positive_aperture_v1']:
        assert reference[key]==pytest.approx(measured[key])
        assert reflected[key]==pytest.approx(measured[key])


def test_nonpositive_signed_flux_keeps_gini_but_flags_m20():
    r=aperture_diagnostics(-np.ones((32,32)),None,(15.5,15.5),10)
    assert r['Gini_aperture_v1']==0
    assert r['M20_aperture_v1'] is None
    assert 'm20_nonpositive_flux' in r['extended_quality_flags']
