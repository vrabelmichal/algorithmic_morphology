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


def test_bilinear_rotation_excludes_masked_contributors():
    from algorithmic_morphology.metrics import _asymmetry_samples
    image = np.ones((32,32))
    good = np.ones_like(image,dtype=bool)
    good[12,12] = False
    image[12,12] = 0
    original, rotated = _asymmetry_samples(image,good,(15.2,15.3),10)
    assert np.all(original == 1)
    assert np.allclose(rotated,1)


def test_sky_asymmetry_uses_same_masked_aperture_denominator():
    from algorithmic_morphology.metrics import _asymmetry_samples, _sky_tile, _rotate
    image = galaxy() + np.random.default_rng(55).normal(0,.5,(96,96))
    mask = np.zeros_like(image,dtype=bool)
    mask[41:47,58:64] = True
    result = measure_cas(image,mask,center=(47.2,47.3),rpetro=17,center_search=0)
    good = ~mask
    cleaned = np.where(good,image,0)
    original, _ = _asymmetry_samples(cleaned,good,(47.2,47.3),25.5)
    tile, valid = _sky_tile(cleaned,good,(47.2,47.3),25.5)
    center = ((tile.shape[1]-1)/2,(tile.shape[0]-1)/2)
    skyok = valid & (_rotate(valid.astype(float),center,order=0)>.5)
    mean = np.abs(tile-_rotate(tile,center))[skyok].mean()
    assert result['method_version'] == 'cas_discrete_v3'
    assert result['A_sky'] == pytest.approx(mean*len(original)/np.abs(original).sum())


def test_center_search_cannot_increase_asymmetry():
    image = galaxy() + np.random.default_rng(9).normal(0,2,(96,96))
    fixed = measure_cas(image,rpetro=17,center_search=0)
    optimized = measure_cas(image,rpetro=17,center_search=2)
    assert optimized['A_raw'] <= fixed['A_raw']


def test_concentration_matches_resolved_gaussian_growth_curve():
    y,x=np.indices((192,192))
    image=40*np.exp(-((x-95.5)**2+(y-95.5)**2)/(2*12**2))
    result=measure_cas(image,center_search=0)
    # Exact infinite circular Gaussian has r_f=sigma*sqrt(-2*log(1-f)).
    expected=2.5*np.log10(np.log(5)/np.log(1.25))
    assert result['C'] == pytest.approx(expected,abs=.03)


def test_powell_empty_direction_falls_back_to_initial_center(monkeypatch):
    def fail(*args,**kwargs):
        raise ValueError('zero-size array to reduction operation maximum which has no identity')
    monkeypatch.setattr('algorithmic_morphology.metrics.minimize',fail)
    r=measure_cas(galaxy(),rpetro=17)
    assert r['A_raw'] == pytest.approx(0)
    assert r['center_asym_x']==47.5
    assert 'center_optimization_failed' in r['quality_flags']
