"""Deterministic synthetic noise/PSF/mask experiments; no physical calibration claim."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from algorithmic_morphology.metrics import measure_cas


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    y,x=np.indices((192,192)); center=(95.5,95.5)
    base=40*np.exp(-((x-center[0])**2+(y-center[1])**2)/(2*12**2))
    clump=12*np.exp(-((x-116)**2+(y-89)**2)/(2*2**2))
    rows=[]
    for morphology,truth in [('symmetric',base),('clump',base+clump)]:
        for psf_sigma in [0,1,2]:
            convolved=gaussian_filter(truth,psf_sigma) if psf_sigma else truth
            for noise_sigma in [0,.5,2,5]:
                for mask_mode in ['none','contaminant']:
                    for seed in range(3):
                        image=convolved+np.random.default_rng(seed).normal(0,noise_sigma,base.shape)
                        mask=np.zeros_like(image,dtype=bool)
                        if mask_mode=='contaminant':
                            star=(x-112)**2+(y-109)**2<4**2
                            image[star]+=1000
                            mask[star]=True
                        result=measure_cas(image,mask,center=center)
                        rows.append(dict(morphology=morphology,psf_sigma=psf_sigma,
                                         noise_sigma=noise_sigma,mask_mode=mask_mode,seed=seed,**result))
    pd.DataFrame(rows).to_csv(args.output/'sensitivity.csv',index=False)
    print(f'Wrote {len(rows)} synthetic measurements to {args.output}')


if __name__=='__main__':
    main()
