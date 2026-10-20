"""Explicit CAS variants; see docs/CAS_METHODS.md for limits and references.

These implementations are literature-motivated, NOT claimed bitwise equivalents
of statmorph or the original authors' unreleased measurement programs.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.ndimage import map_coordinates, uniform_filter
from scipy.optimize import minimize


@dataclass(frozen=True)
class SmoothnessProfile:
    name: str
    width_fraction: float
    excluded_fraction: float
    scale: float


PROFILES = (
    SmoothnessProfile("lotz2004_statmorph_like_v1", 0.25, 0.25, 1.0),
    SmoothnessProfile("conselice2003_like_v1", 0.30, 0.075, 10.0),
)


def _radius(shape, center):
    y, x = np.indices(shape)
    return np.hypot(x - center[0], y - center[1])


def _rotate(image, center, order=1, cval=0):
    y, x = np.indices(image.shape, dtype=float)
    return map_coordinates(image, [2 * center[1] - y, 2 * center[0] - x],
                           order=order, mode="constant", cval=cval,
                           prefilter=False)


def petrosian_radius(image, good, center, eta=0.2, minimum=3.0):
    """First annular eta crossing using one-pixel circular pixel-center bins.

    This discrete-annulus method is *not* statmorph's fractional photutils
    aperture implementation. No fallback to a fitted effective radius.
    """
    edge = min(center[0], center[1], image.shape[1]-1-center[0],
               image.shape[0]-1-center[1])
    limit = int(edge - 1)
    if limit < max(5, minimum + 2):
        raise ValueError("insufficient_aperture_space")
    rad = _radius(image.shape, center)
    bins = np.floor(rad[good]).astype(int)
    count = np.bincount(bins, minlength=limit+1)
    flux = np.bincount(bins, weights=image[good], minlength=limit+1)
    count_cum, flux_cum = np.cumsum(count), np.cumsum(flux)
    ratios = np.full(limit+1, np.nan)
    for r in range(max(1, int(math.ceil(minimum))), limit):
        if count[r] < 3 or count_cum[r] == 0 or flux_cum[r] <= 0:
            continue
        local = flux[r] / count[r]
        enclosed = flux_cum[r] / count_cum[r]
        ratios[r] = local / enclosed
    indices = np.where(np.isfinite(ratios))[0]
    for prev, cur in zip(indices[:-1], indices[1:]):
        if cur != prev + 1:
            continue
        if ratios[prev] > eta >= ratios[cur]:
            fraction = (ratios[prev] - eta) / (ratios[prev] - ratios[cur])
            return float(prev + fraction)
    raise ValueError("petrosian_crossing_not_found")


def _quantile_radius(image, good, center, outer_radius, fraction):
    rad = _radius(image.shape, center)
    inside = good & (rad <= outer_radius)
    sorted_indices = np.argsort(rad[inside])
    distances = rad[inside][sorted_indices]
    pixel_flux = image[inside][sorted_indices]
    if pixel_flux.size < 10 or float(np.sum(pixel_flux)) <= 0:
        raise ValueError("nonpositive_or_empty_total_flux")
    cumsum = np.cumsum(pixel_flux)
    total = cumsum[-1]
    # In low-S/N data signed pixels make the growth curve nonmonotonic.
    if np.any(np.diff(cumsum) < -0.01 * total):
        raise ValueError("nonmonotone_growth_curve")
    above = np.flatnonzero(cumsum >= fraction * total)
    if above.size == 0:
        raise ValueError("flux_fraction_not_found")
    k = int(above[0])
    if k == 0:
        return float(distances[0])
    delta = cumsum[k] - cumsum[k-1]
    t = (fraction*total - cumsum[k-1]) / delta if delta > 0 else 0
    return float(distances[k-1] + t*(distances[k]-distances[k-1]))


def _asymmetry_samples(image, good, center, radius):
    # Sample only the aperture bounding box; coordinates remain in the full image.
    x0 = max(0, int(np.floor(center[0] - radius)))
    x1 = min(image.shape[1], int(np.ceil(center[0] + radius)) + 1)
    y0 = max(0, int(np.floor(center[1] - radius)))
    y1 = min(image.shape[0], int(np.ceil(center[1] + radius)) + 1)
    y, x = np.mgrid[y0:y1, x0:x1]
    coords = [2 * center[1] - y, 2 * center[0] - x]
    rotated = map_coordinates(image, coords, order=1, mode="constant", cval=0,
                              prefilter=False)
    # All bilinear contributors must be valid, including at mask boundaries.
    support = map_coordinates(good.astype(float), coords, order=1,
                              mode="constant", cval=0, prefilter=False)
    inside = (np.hypot(x-center[0], y-center[1]) <= radius)
    inside &= good[y0:y1, x0:x1] & (support >= 1 - 1e-12)
    return image[y0:y1, x0:x1][inside], rotated[inside]


def _asymmetry_raw(image, good, center, radius):
    original, rotated = _asymmetry_samples(image, good, center, radius)
    denom = float(np.abs(original).sum())
    if original.size < 10 or denom <= 0:
        return float("inf"), 0
    return float(np.abs(original - rotated).sum() / denom), int(original.size)


def _boxcar(image, good, width):
    # Mask-normalized smoothing: a masked star never contributes zero flux
    # to the kernel and never produces a fake residual ring.
    vals = uniform_filter(np.where(good, image, 0.0), size=width, mode="reflect")
    weights = uniform_filter(good.astype(float), size=width, mode="reflect")
    return np.divide(vals, weights, out=np.full(image.shape, np.nan), where=weights > 0)


def _sky_tile(image, good, center, forbidden_radius, width=32):
    """Choose clean full-size corner tile outside CAS aperture, or return None."""
    ny, nx = image.shape
    k = min(width, ny // 4, nx // 4)
    if k < 8:
        return None
    candidates = []
    for ys in (slice(0, k), slice(ny-k, ny)):
        for xs in (slice(0, k), slice(nx-k, nx)):
            sub = image[ys, xs]
            mask = good[ys, xs]
            rr = _radius(image.shape, center)[ys, xs]
            if mask.mean() < 0.95 or np.any(rr <= forbidden_radius + 2):
                continue
            candidates.append((np.std(sub[mask]), sub, mask))
    if not candidates:
        return None
    _, sub, mask = min(candidates, key=lambda t:t[0])
    return sub.copy(), mask.copy()


def measure_cas(image, mask=None, center=None, *, rpetro=None,
                skybox_size=32, center_search=2.0):
    """Measure C, raw/corrected A and both S variants on one 2D image.

    Mask is True for contaminating sources, NOT the sky or target galaxy.
    Invalid individual metrics are None with human-readable quality flags.
    A correction is sky residual *valid aperture pixel count / galaxy |flux|;
    S uses the analogous mean positive sky residual scaled by annular count.
    """
    arr = np.asarray(image, dtype=float)
    if arr.ndim != 2 or min(arr.shape) < 16:
        raise ValueError("image_must_be_2d_at_least_16_pixels")
    if mask is None:
        contamination = np.zeros(arr.shape, dtype=bool)
    else:
        contamination = np.asarray(mask, dtype=bool)
        if contamination.shape != arr.shape:
            raise ValueError("mask_shape_mismatch")
    if not np.isfinite(center_search) or center_search < 0:
        raise ValueError("invalid_center_search")
    if skybox_size < 8:
        raise ValueError("invalid_skybox_size")
    good = (~contamination) & np.isfinite(arr)
    cleaned = np.where(good, arr, 0.0)
    if center is None:
        center = ((arr.shape[1]-1)/2, (arr.shape[0]-1)/2)
    center = np.asarray(center, dtype=float)
    if center.shape != (2,) or not np.isfinite(center).all():
        raise ValueError("invalid_center")
    if not (0 <= center[0] < arr.shape[1] and 0 <= center[1] < arr.shape[0]):
        raise ValueError("center_outside_cutout")

    out = {"method_version": "cas_discrete_v2", "center_initial_x": float(center[0]),
           "center_initial_y": float(center[1]), "mask_fraction": float(1-good.mean())}
    flags = []
    try:
        rp = float(rpetro) if rpetro is not None else petrosian_radius(cleaned, good, center)
        if not np.isfinite(rp) or rp <= 0:
            raise ValueError("invalid_rpetro")
        aperture = 1.5 * rp
        margin = min(center[0], center[1], arr.shape[1]-1-center[0], arr.shape[0]-1-center[1])
        if aperture >= margin - center_search - 0.5:
            raise ValueError("aperture_truncated")
    except ValueError as exc:
        out.update({"rpetro": None, "C": None, "r20": None, "r80": None,
                    "A_raw": None, "A_sky": None, "A_corrected": None,
                    "center_asym_x": None, "center_asym_y": None})
        flags.append(str(exc))
        for profile in PROFILES:
            for suffix in ("raw", "sky", "corrected"):
                out[f"S_{profile.name}_{suffix}"] = None
        out["quality_flags"] = "|".join(sorted(set(flags)))
        return out

    out["rpetro"] = rp
    out["aperture_radius"] = aperture
    # Minimize the *raw* A once; both variants use the identical center.
    search = max(0.0, min(float(center_search), 0.25*rp))
    bounds = [(center[0]-search, center[0]+search),
              (center[1]-search, center[1]+search)]
    if search > 0:
        opt = minimize(lambda pos: _asymmetry_raw(cleaned, good, pos, aperture)[0],
                       center, method="Powell", bounds=bounds,
                       options={"xtol":1e-3,"ftol":1e-5,"maxiter":120})
        best = opt.x if (opt.success and np.isfinite(opt.fun)
                             and opt.fun <= _asymmetry_raw(cleaned, good, center, aperture)[0]) else center
        if not opt.success:
            flags.append("center_optimization_failed")
        if np.any(np.isclose(best, np.array(bounds)[:,0], atol=1e-2)) or np.any(np.isclose(best, np.array(bounds)[:,1], atol=1e-2)):
            flags.append("center_search_boundary")
    else:
        best = center
    out["center_asym_x"],out["center_asym_y"] = map(float,best)
    a_raw, a_pixels = _asymmetry_raw(cleaned, good, best, aperture)
    if not np.isfinite(a_raw):
        a_raw = None
        flags.append("asymmetry_no_flux")
    out["A_raw"] = a_raw
    a_original, _ = _asymmetry_samples(cleaned, good, best, aperture)
    denom = float(np.abs(a_original).sum())
    out["asymmetry_valid_pixels"] = a_pixels
    out["asymmetry_absolute_flux"] = denom
    sky = _sky_tile(cleaned, good, best, aperture, width=skybox_size)
    if sky is None:
        flags.append("no_clean_sky_tile")
        asym_sky = None
    else:
        skyimg, skyvalid = sky
        skycenter = ((skyimg.shape[1]-1)/2, (skyimg.shape[0]-1)/2)
        skyrot = _rotate(skyimg, skycenter)
        skyvalidrot = _rotate(skyvalid.astype(float), skycenter, order=0) > 0.5
        skyok = skyvalid & skyvalidrot
        sky_mean = float(np.abs(skyimg[skyok]-skyrot[skyok]).mean()) if skyok.any() else float("nan")
        asym_sky = float(sky_mean*a_pixels/denom) if np.isfinite(sky_mean) and denom > 0 else None
    out["A_sky"] = asym_sky
    out["A_corrected"] = a_raw-asym_sky if a_raw is not None and asym_sky is not None else None

    try:
        r20 = _quantile_radius(cleaned, good, best, aperture, 0.2)
        r80 = _quantile_radius(cleaned, good, best, aperture, 0.8)
        if r20 <= 0 or r80 <= r20:
            raise ValueError("invalid_flux_quantiles")
        out.update({"r20":r20,"r80":r80,"C":float(5*np.log10(r80/r20))})
    except ValueError as exc:
        flags.append(str(exc))
        out.update({"r20":None,"r80":None,"C":None})

    radial = _radius(arr.shape, best)
    for profile in PROFILES:
        prefix = f"S_{profile.name}"
        kernel = int(profile.width_fraction*rp)
        out[prefix+"_boxcar_width_pixels"] = kernel
        out[prefix+"_inner_radius"] = profile.excluded_fraction*rp
        region = good & (radial > profile.excluded_fraction*rp) & (radial <= aperture)
        flux = float(cleaned[region].sum())
        if kernel < 2 or region.sum() < 10 or flux <= 0:
            flags.append(prefix+"_insufficient_scale_or_flux")
            out.update({prefix+"_raw":None,prefix+"_sky":None,prefix+"_corrected":None})
            continue
        smooth = _boxcar(cleaned, good, kernel)
        resid = np.maximum(cleaned-smooth, 0.0)
        # Valid where the normalized kernel has at least one good pixel.
        region = region & np.isfinite(smooth)
        numerator = float(resid[region].sum())
        raw = profile.scale*numerator/flux
        if sky is None:
            sky_term = None
        else:
            skyimg, skyvalid = sky
            if kernel > min(skyimg.shape):
                flags.append(prefix+"_sky_kernel_exceeds_tile")
                sky_term = None
            else:
                sky_smooth = _boxcar(skyimg, skyvalid, kernel)
                sky_resid = np.maximum(skyimg-sky_smooth, 0.0)
                sky_term = (profile.scale * float(np.nanmean(sky_resid[skyvalid]))
                            * int(region.sum()) / flux)
        out.update({prefix+"_raw":float(raw), prefix+"_sky":sky_term,
                    prefix+"_corrected":raw-sky_term if sky_term is not None else None})
    out["quality_flags"] = "|".join(sorted(set(flags)))
    return out
