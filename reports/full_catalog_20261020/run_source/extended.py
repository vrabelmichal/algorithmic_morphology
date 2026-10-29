"""Versioned aperture diagnostics; not statmorph segmentation-based outputs.

Gini uses absolute signed pixels. M20 uses signed moments at the signed-flux
centroid and includes pixels through the first crossing of 20% total flux.
Shape uses positive flux weights explicitly; it never changes CAS photometry.
"""
import numpy as np


EXTRA_METRICS = ['Gini_aperture_v1','M20_aperture_v1','ellipticity_positive_aperture_v1',
                 'axis_ratio_positive_aperture_v1','position_angle_positive_aperture_v1']


def aperture_diagnostics(image, mask, center, radius):
    result = {name:None for name in EXTRA_METRICS}
    result['aperture_diagnostics_version']='aperture_diagnostics_v1'
    if radius is None or center[0] is None:
        return result
    image=np.asarray(image,dtype=float)
    y,x=np.indices(image.shape,dtype=float)
    region=np.isfinite(image)&(np.hypot(x-center[0],y-center[1])<=radius)
    if mask is not None:
        region &= ~np.asarray(mask,dtype=bool)
    result['aperture_valid_pixels']=int(region.sum())
    geometric=np.hypot(x-center[0],y-center[1])<=radius
    result['aperture_mask_fraction']=float(1-region.sum()/geometric.sum())
    flux=image[region]; xx=x[region]; yy=y[region]; n=len(flux)
    flags=[]
    if n<10 or np.abs(flux).sum()<=0:
        result['extended_quality_flags']='insufficient_aperture_pixels_or_flux'
        return result
    ordered=np.sort(np.abs(flux))
    result['Gini_aperture_v1']=float(np.sum((2*np.arange(1,n+1)-n-1)*ordered)/((n-1)*ordered.sum()))
    total=flux.sum()
    if total>0:
        cx=np.sum(flux*xx)/total;cy=np.sum(flux*yy)/total
        moments=flux*((xx-cx)**2+(yy-cy)**2)
        order=np.argsort(-flux,kind='stable'); cumulative=np.cumsum(flux[order])
        crossing=np.flatnonzero(cumulative>=.2*total)
        if len(crossing):
            # Include all ties at the crossing brightness to avoid orientation bias.
            threshold=flux[order[crossing[0]]]
            bright=flux>=threshold
            numerator=moments[bright].sum();denominator=moments.sum()
            if numerator>0 and denominator>0:
                result['M20_aperture_v1']=float(np.log10(numerator/denominator))
            else:flags.append('m20_nonpositive_moment')
        result['centroid_signed_x']=float(cx);result['centroid_signed_y']=float(cy)
    else:flags.append('m20_nonpositive_flux')
    positive=np.maximum(flux,0); weight=positive.sum()
    if weight>0:
        cx=np.sum(positive*xx)/weight;cy=np.sum(positive*yy)/weight
        dx=xx-cx;dy=yy-cy
        cov=np.array([[np.sum(positive*dx*dx),np.sum(positive*dx*dy)],
                      [np.sum(positive*dx*dy),np.sum(positive*dy*dy)]])/weight
        values,vectors=np.linalg.eigh(cov)
        if values[1]>0 and values[0]>=0:
            ratio=float(np.sqrt(values[0]/values[1]));v=vectors[:,1]
            result['axis_ratio_positive_aperture_v1']=ratio
            result['ellipticity_positive_aperture_v1']=1-ratio
            result['position_angle_positive_aperture_v1']=float(np.degrees(np.arctan2(v[1],v[0]))%180)
    result['extended_quality_flags']='|'.join(flags)
    return result
