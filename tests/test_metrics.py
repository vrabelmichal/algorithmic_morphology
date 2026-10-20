import numpy as np
import pytest
from algorithmic_morphology.metrics import measure_cas, _boxcar


def galaxy(shape=(96, 96), shift=0):
    y,x = np.indices(shape)
    return 40*np.exp(-((x-47.5)**2+(y-47.5)**2)/(2*9**2)) + shift


def test_symmetric_galaxy_two_s_profiles():
    r = measure_cas(galaxy(), rpetro=17, center_search=0)
    assert r['A_raw'] < 1.e-9
    assert r['C'] > 0
    assert r['S_lotz2004_statmorph_like_v1_raw'] is not None
    assert r['S_conselice2003_like_v1_raw'] is not None
    assert r['S_conselice2003_like_v1_raw'] != r['S_lotz2004_statmorph_like_v1_raw']


def test_asymmetric_clump_increases_raw():
    base=galaxy()
    changed=base.copy(); changed[40:43, 62:65] += 20
    a=measure_cas(base,rpetro=17,center_search=0)
    b=measure_cas(changed,rpetro=17,center_search=0)
    assert b['A_raw']>a['A_raw']
    assert b['S_lotz2004_statmorph_like_v1_raw']>a['S_lotz2004_statmorph_like_v1_raw']


def test_mask_normalized_smooth_prevents_zero_halo():
    a=np.ones((40,40)); a[10:13,10:13] = 10000
    mask=np.zeros_like(a,dtype=bool); mask[10:13,10:13]=True
    s=_boxcar(a,~mask,7)
    assert np.allclose(s,1)


def test_background_correction_and_flag():
    r=measure_cas(galaxy()+np.random.default_rng(1).normal(0,.01,(96,96)),rpetro=17,center_search=0)
    assert r['A_sky'] is not None
    assert r['A_corrected']==pytest.approx(r['A_raw']-r['A_sky'])
    assert r['S_conselice2003_like_v1_corrected']==pytest.approx(r['S_conselice2003_like_v1_raw']-r['S_conselice2003_like_v1_sky'])


def test_no_sky_returns_raw_only():
    im=galaxy((64,64))
    r=measure_cas(im,center=(31.5,31.5),rpetro=16,center_search=0)
    assert r['A_raw'] is not None
    assert r['A_corrected'] is None
    assert 'no_clean_sky_tile' in r['quality_flags']


def test_nonfinite_and_bad_mask():
    im=galaxy(); im[1,1]=np.nan
    with pytest.raises(ValueError, match='mask_shape_mismatch'):
        measure_cas(im,np.zeros((30,30)))
    r=measure_cas(im,rpetro=17,center_search=0)
    assert r['C']>0


def test_petrosian_and_truncation():
    r=measure_cas(galaxy(),center_search=0)
    assert r['rpetro']>0
    out=measure_cas(galaxy(),rpetro=32,center_search=0)
    assert out['A_raw'] is None
    assert 'aperture_truncated' in out['quality_flags']
