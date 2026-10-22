"""Write standalone diagnostic plots and a target/mask inspection sheet."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from algorithmic_morphology.metrics import PROFILES
from algorithmic_morphology.runner import _read_fits


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report-dir',type=Path,required=True)
    args=parser.parse_args()
    frame=pd.read_parquet(next(args.report_dir.glob('measurements__*.parquet')))
    diagnostics=[]
    for row in frame[frame.image_file.notna()].itertuples():
        mask,_=_read_fits(row.mask_file)
        good=mask==0; y,x=np.indices(mask.shape)
        aperture=np.hypot(x-row.center_asym_x,y-row.center_asym_y)<=row.aperture_radius
        corners=[good[:32,:32],good[:32,-32:],good[-32:,:32],good[-32:,-32:]]
        diagnostics.append(dict(coadd_object_id=row.coadd_object_id,band=row.band,
                                aperture_mask_fraction=float(1-good[aperture].mean()),
                                best_corner_valid_fraction=max(float(c.mean()) for c in corners)))
    pd.DataFrame(diagnostics).to_csv(args.report_dir/'mask_diagnostics.csv',index=False)
    colors={'g':'#187a5c','r':'#a64338','i':'#6856a5'}
    fig,axes=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    for band,part in frame.groupby('band'):
        axes[0,0].hist(part.C.dropna(),bins=np.linspace(1.8,4.2,22),histtype='step',label=f'{band} n={part.C.notna().sum()}',color=colors[band])
        axes[0,1].scatter(part.A_raw,part.A_corrected,s=16,label=f'{band} n={part.A_corrected.notna().sum()}',color=colors[band],alpha=.7)
        axes[1,0].scatter(part[f'S_{PROFILES[0].name}_raw'],part[f'S_{PROFILES[1].name}_raw']/10,s=12,label=band,color=colors[band],alpha=.6)
    axes[0,0].set(xlabel='C',ylabel='Available candidates',title='Concentration; availability-selected sample')
    axes[0,1].axhline(0,color='grey',lw=1)
    axes[0,1].set(xlabel='A raw',ylabel='A corrected',title='Only 33 images have a corner sky correction')
    axes[1,0].plot([0,1.5],[0,1.5],'--',color='grey',lw=1)
    axes[1,0].set(xlabel='S Lotz-like raw',ylabel='S Conselice-like raw / 10',title='Removing ×10 still leaves profile differences')
    s=pd.read_csv(args.report_dir/'sensitivity_summary.csv')
    sub=s[(s.morphology=='symmetric')&(s.mask_mode=='none')&(s.psf_sigma==0)]
    axes[1,1].plot(sub.noise_sigma,sub.A_raw,'o-',label='A raw')
    axes[1,1].plot(sub.noise_sigma,sub.A_corrected,'s-',label='A corrected')
    axes[1,1].axhline(0,color='grey',lw=1)
    axes[1,1].set(xlabel='Injected Gaussian noise sigma (peak=40)',ylabel='Median A, three seeds',title='Symmetric synthetic galaxy: noise sensitivity')
    for ax in axes.flat:
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    fig.savefig(args.report_dir/'overview.png',dpi=160)
    plt.close(fig)
    # Select three g-band rows spanning stored full-image mask fraction.
    candidates=frame[(frame.band=='g')&frame.image_file.notna()].sort_values('mask_fraction')
    chosen=candidates.iloc[[0,len(candidates)//2,len(candidates)-1]]
    fig,axes=plt.subplots(3,2,figsize=(9,12),constrained_layout=True)
    for row_index,(_,row) in enumerate(chosen.iterrows()):
        data,_=_read_fits(row.image_file);mask,_=_read_fits(row.mask_file)
        cx,cy=row.center_initial_x,row.center_initial_y
        x0,y0=int(cx)-90,int(cy)-90
        stamp=data[y0:y0+180,x0:x0+180]
        ax=axes[row_index,0];v=np.nanpercentile(stamp,[5,99])
        ax.imshow(stamp,origin='lower',cmap='gray',vmin=v[0],vmax=v[1])
        ax.scatter([cx-x0],[cy-y0],marker='+',color='cyan',s=50)
        ax.add_patch(plt.Circle((row.center_asym_x-x0,row.center_asym_y-y0),row.aperture_radius,fill=False,color='orange',lw=1))
        ax.set_title(f'ID {row.coadd_object_id}, g; C={row.C:.2f}, A raw={row.A_raw:.2f}')
        axes[row_index,1].imshow(mask[y0:y0+180,x0:x0+180]!=0,origin='lower',cmap='gray',vmin=0,vmax=1)
        axes[row_index,1].scatter([cx-x0],[cy-y0],marker='+',color='cyan',s=50)
        axes[row_index,1].set_title(f'Contaminants white; full-image mask={row.mask_fraction:.1%}')
        for ax in axes[row_index]: ax.set(xlabel='Cutout x pixel',ylabel='Cutout y pixel')
    fig.savefig(args.report_dir/'target_mask_checks.png',dpi=130)
    plt.close(fig)


if __name__=='__main__':
    main()
