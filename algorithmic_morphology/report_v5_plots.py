"""Single-panel figures for report revision 5, grouped by morphology family.

Every function draws exactly one matplotlib axes and saves one PNG. Multi-panel
layouts are composed in HTML, so desktop grids can collapse to a single column
on small screens. Plot crops never alter values; off-screen counts are audited.
"""
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator,FuncFormatter,NullFormatter,NullLocator
from PIL import Image

from .report_data import BANDS,finite,robust_limits
from .paper_plots import LABELS,L,K
from .classification_data import HUMAN_LABELS,CLASS_LABELS

# Validated with the dataviz palette checker (lightness, chroma, CVD, contrast).
BAND_COLORS={'g':'#00916b','r':'#cc4a30','i':'#7a4fc4'}
# Line style is a secondary, colour-independent band encoding.
BAND_STYLES={'g':'-','r':'--','i':':'}
PAIR_COLORS=('#2a62b8','#d9761a')
CLASS_COLORS=['#2a78d6','#eb6834','#1baf7a','#eda100','#e87ba4']
CLASS_STYLES=['-','--','-.',':',(0,(5,1,1,1))]

FAMILIES={
 'C':dict(name='Concentration',metrics=['C']),
 'A':dict(name='Asymmetry',metrics=['A_raw','A_corrected']),
 'S':dict(name='Clumpiness',metrics=[L+'_raw',L+'_corrected',K+'_raw',K+'_corrected']),
 'D':dict(name='Supplementary aperture diagnostics',metrics=['Gini_aperture_v1','M20_aperture_v1','axis_ratio_positive_aperture_v1']),
}
PLAIN={
 'C':'Concentration C','A_raw':'Raw asymmetry A_raw','A_sky':'Sky asymmetry term A_sky','A_corrected':'Sky-corrected asymmetry A_corr',
 L+'_raw':'Lotz-like clumpiness S_L, raw',L+'_corrected':'Lotz-like clumpiness S_L, sky-corrected',
 K+'_raw':'Conselice-like clumpiness S_C, raw',K+'_corrected':'Conselice-like clumpiness S_C, sky-corrected',
 'Gini_aperture_v1':'Aperture Gini G','M20_aperture_v1':'Aperture M20',
 'axis_ratio_positive_aperture_v1':'Positive-flux axis ratio q+','ellipticity_positive_aperture_v1':'Positive-flux ellipticity 1-q+',
 'position_angle_positive_aperture_v1':'Image major-axis angle','rpetro':'Petrosian radius R_P',
 'gf_mu_mean':'Mean effective surface brightness','gf_re_arcsec':'GALFIT effective radius R_e','gf_re_pixels':'GALFIT effective radius R_e (pixels)',
 'gf_mag':'GALFIT total magnitude','gf_n':'GALFIT Sersic index n','gf_q':'GALFIT axis ratio q','gf_chi2red':'GALFIT reduced chi-square',
}
FAMILY_OF={m:f for f,v in FAMILIES.items() for m in v['metrics']}
PROPERTIES=['gf_n','gf_mu_mean','gf_re_arcsec','gf_mag','gf_q','gf_chi2red']
LOG_PROPERTIES={'gf_re_arcsec','gf_chi2red'}
HUMAN_SELECTIONS={
 'analysis-morphology-etg-min2':'ETG / elliptical','analysis-morphology-ltg-min2':'LTG / spiral',
 'analysis-morphology-featureless-min2':'Featureless','analysis-morphology-irregular-min2':'Not sure (Irr/other)',
 'analysis-undecidable-morphology':'Split (no consensus)',
}
EFFECT_SELECTIONS={
 'analysis-morphology-etg-min2':'Human: ETG / elliptical','analysis-morphology-ltg-min2':'Human: LTG / spiral',
 'analysis-morphology-featureless-min2':'Human: Featureless','analysis-morphology-irregular-min2':'Human: Not sure (Irr/other)',
 'analysis-visible-nucleus-visible-min2':'Human: visible nucleus','analysis-visible-nucleus-not-visible-min2':'Human: nucleus not visible',
 'analysis-is-lsb-yes-min2':'Human: is LSB','analysis-is-lsb-no-min2':'Human: not LSB',
 'analysis-failed-fitting-flagged-min2':'Human: failed fitting','analysis-failed-fitting-good-fit-min2':'Human: good fitting',
 **{f'zoobot_top_{i}':'DECaLS10 top: '+name for i,name in enumerate(CLASS_LABELS)},
}
TARGETS={**{k:'Human: '+v for k,v in HUMAN_LABELS.items()},**{f'class_{i}':'DECaLS10: '+name for i,name in enumerate(CLASS_LABELS)}}
COHORTS=[('available','All available'),('human_zoobot_matched','Human+DECaLS10 matched'),('three_plus_human_votes','≥3 human votes'),('galfit_adjusted','GALFIT-adjusted')]


def style():
    plt.rcParams.update({'font.size':10.5,'axes.titlesize':11.5,'axes.labelsize':10.5,'font.family':'DejaVu Sans',
                         'mathtext.fontset':'stix','axes.spines.top':False,'axes.spines.right':False,'legend.fontsize':9,
                         'axes.edgecolor':'#8a949c','axes.labelcolor':'#25303a','xtick.color':'#4c5863','ytick.color':'#4c5863',
                         'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.alpha':.16,'grid.linewidth':.6})


def write_png(fig,path,dpi):
    """Save as a 256-colour palette PNG: visually lossless for flat plot graphics, ~3x smaller."""
    buffer=BytesIO();fig.savefig(buffer,format='png',dpi=dpi,bbox_inches='tight');plt.close(fig)
    image=Image.open(buffer).convert('RGB').quantize(256,method=Image.Quantize.MEDIANCUT,dither=Image.Dither.NONE)
    image.save(path,optimize=True);return path


def plain_log_axis(axis,lo,hi):
    ticks=[t for t in [.2,.3,.5,.7,1,1.2,1.5,2,3,5,7,10,15,20,30,50,70,100] if lo<=t<=hi]
    axis.set_major_locator(FixedLocator(ticks));axis.set_major_formatter(FuncFormatter(lambda v,_:f'{v:g}'))
    axis.set_minor_locator(NullLocator());axis.set_minor_formatter(NullFormatter())


class Panels:
    """Writes single-axes PNGs and records the display window of every series."""

    def __init__(self,root,dpi=140):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True);self.dpi=dpi;self.records=[];style()

    def save(self,fig,name):
        return write_png(fig,self.root/(name+'.png'),self.dpi)

    def audit(self,panel,series,values,limits,band=None,axis='x'):
        v=finite(values);v=v[np.isfinite(v)];lo,hi=limits
        self.records.append(dict(panel=panel,series=series,band=band,axis=axis,finite=len(v),lower=lo,upper=hi,
                                 below=int((v<lo).sum()),above=int((v>hi).sum()),shown=int(((v>=lo)&(v<=hi)).sum())))

    # ---------- distributions ----------
    def band_histogram(self,frame,column,limits,name,title,scale=1.,xlabel=None):
        fig,ax=plt.subplots(figsize=(6.2,4.2));bins=np.linspace(*limits,61)
        for band in BANDS:
            vals=finite(frame.loc[frame.band==band,column])/scale;vals=vals[np.isfinite(vals)]
            if len(vals):ax.hist(vals,bins=bins,weights=np.full(len(vals),100/len(vals)),histtype='step',lw=1.8,color=BAND_COLORS[band],ls=BAND_STYLES[band],label=f'{band} band (n={len(vals):,})')
            self.audit(name,column,vals,limits,band)
        if limits[0]<0<limits[1]:ax.axvline(0,color='#6b7280',lw=.8)
        ax.set(xlim=limits,xlabel=xlabel or LABELS[column],ylabel='Share of finite sample per bin (%)',title=title);ax.legend(frameon=False)
        return self.save(fig,name)

    def pair_histogram(self,frame,band,columns,labels,limits,name,title,xlabel,scales=(1,1)):
        fig,ax=plt.subplots(figsize=(6.2,4.2));bins=np.linspace(*limits,71);part=frame[frame.band==band]
        for col,label,color,ls,scale in zip(columns,labels,PAIR_COLORS,['-','--'],scales):
            vals=finite(part[col])/scale;vals=vals[np.isfinite(vals)]
            ax.hist(vals,bins=bins,weights=np.full(len(vals),100/max(1,len(vals))),histtype='step',lw=1.9,ls=ls,color=color,label=f'{label} (n={len(vals):,})')
            self.audit(name,col+('/10' if scale==10 else ''),vals,limits,band)
        if limits[0]<0<limits[1]:ax.axvline(0,color='#6b7280',lw=.8)
        ax.set(xlim=limits,xlabel=xlabel,ylabel='Share of finite sample per bin (%)',title=title);ax.legend(frameon=False)
        return self.save(fig,name)

    def density(self,frame,band,x,y,xl,yl,name,title,div_y=1,equality=False,xlabel=None,ylabel=None):
        fig,ax=plt.subplots(figsize=(6.2,4.8))
        part=frame[frame.band==band][[x,y]].replace([np.inf,-np.inf],np.nan).dropna()
        a=part[x].to_numpy(float);b=part[y].to_numpy(float)/div_y
        inside=(a>=xl[0])&(a<=xl[1])&(b>=yl[0])&(b<=yl[1])
        if inside.any():
            h=ax.hexbin(a[inside],b[inside],gridsize=48,extent=(*xl,*yl),mincnt=1,bins='log',cmap='cividis',rasterized=True,linewidths=0)
            fig.colorbar(h,ax=ax,label='Galaxies per cell (log scale)',shrink=.85)
        rho=part[x].corr(part[y],method='spearman') if len(part)>2 else np.nan
        if equality:
            lo=max(xl[0],yl[0]);hi=min(xl[1],yl[1]);ax.plot([lo,hi],[lo,hi],'--',color='#cc4a30',lw=1.1,label='equality line');ax.legend(frameon=False,loc='upper left')
        ax.set(xlim=xl,ylim=yl,xlabel=xlabel or LABELS[x],ylabel=ylabel or (LABELS[y]+(' / 10' if div_y==10 else '')),
               title=f'{title}\npaired n={len(part):,}; shown {inside.sum():,}; Spearman ρ={rho:.3f}');ax.grid(alpha=.1)
        self.audit(name,x,a,xl,band,'x');self.audit(name,y+('/10' if div_y==10 else ''),b,yl,band,'y')
        return self.save(fig,name),dict(n=len(part),rho=rho)

    # ---------- GALFIT ----------
    def galfit_trend(self,trends,prop,metric,ylimits,name,title):
        fig,ax=plt.subplots(figsize=(6.2,4.3));log=prop in LOG_PROPERTIES
        part=trends[(trends.property==prop)&(trends.metric==metric)]
        for band in BANDS:
            d=part[part.band==band].sort_values('left')
            centers=np.sqrt(d.left*d.right) if log else (d.left+d.right)/2
            ax.plot(centers,d['median'],color=BAND_COLORS[band],ls=BAND_STYLES[band],lw=1.9,label=f'{band} band (pairs={int(d.n.sum()):,})')
            ax.fill_between(centers,d.p16,d.p84,color=BAND_COLORS[band],alpha=.09,lw=0)
        ax.set(xlim=(part.left.min(),part.right.max()),ylim=ylimits,xlabel=LABELS[prop],ylabel=LABELS[metric],title=title)
        if log:ax.set_xscale('log');plain_log_axis(ax.xaxis,part.left.min(),part.right.max())
        if prop in ('gf_mu_mean','gf_mag'):ax.annotate('fainter →',xy=(.98,.02),xycoords='axes fraction',ha='right',fontsize=8.5,color='#4c5863')
        if ylimits[0]<0<ylimits[1]:ax.axhline(0,color='#6b7280',lw=.7)
        ax.legend(frameon=False,fontsize=8.5);return self.save(fig,name)

    def galfit_matrix(self,correlations,band,metrics,name,title):
        sub=correlations[(correlations.cohort=='gf_unflagged')&(correlations.band==band)]
        mat=sub.pivot(index='metric',columns='property',values='spearman').reindex(index=metrics,columns=PROPERTIES)
        fig,ax=plt.subplots(figsize=(7.4,6.4))
        im=ax.imshow(mat.to_numpy(float),vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        for yy in range(len(metrics)):
            for xx in range(len(PROPERTIES)):
                v=mat.iloc[yy,xx];ax.text(xx,yy,f'{v:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(v)>.55 else '#1f2933')
        ax.set_xticks(range(len(PROPERTIES)),['Sersic n','Mean SB','R_e','Magnitude','Axis ratio q','Red. χ²'],rotation=35,ha='right')
        ax.set_yticks(range(len(metrics)),[PLAIN[m] for m in metrics]);ax.grid(False)
        for yy in [.5,2.5,6.5]:ax.axhline(yy,color='white',lw=3)
        fig.colorbar(im,ax=ax,label='Spearman rank correlation ρ',shrink=.8);ax.set_title(title)
        return self.save(fig,name)

    # ---------- human / DECaLS10 ----------
    def association_heatmap(self,assoc,metric,name,title):
        targets=list(TARGETS);columns=[(c,b) for c,_ in COHORTS for b in BANDS]
        mat=np.full((len(targets),len(columns)),np.nan)
        for j,(cohort,band) in enumerate(columns):
            stat='partial_rank' if cohort=='galfit_adjusted' else 'spearman'
            s=assoc[(assoc.metric==metric)&(assoc.cohort==cohort)&(assoc.band==band)].set_index('target')[stat]
            mat[:,j]=s.reindex(targets).to_numpy(float)
        fig,ax=plt.subplots(figsize=(9.6,7.6))
        im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        for yy in range(len(targets)):
            for xx in range(len(columns)):
                v=mat[yy,xx]
                if np.isfinite(v):ax.text(xx,yy,f'{v:.2f}',ha='center',va='center',fontsize=7.6,color='white' if abs(v)>.55 else '#1f2933')
        ax.set_xticks(range(len(columns)),[b for _,b in columns]);ax.set_yticks(range(len(targets)),[TARGETS[t] for t in targets]);ax.grid(False)
        for xx in [2.5,5.5,8.5]:ax.axvline(xx,color='white',lw=3)
        ax.axhline(len(HUMAN_LABELS)-.5,color='white',lw=3)
        sec=ax.secondary_xaxis('top');sec.set_xticks([1,4,7,10],[label for _,label in COHORTS],fontsize=8.5);sec.tick_params(length=0)
        ax.set_xlabel('Photometric band within each comparison cohort')
        fig.colorbar(im,ax=ax,label='Rank correlation (partial for GALFIT-adjusted)',shrink=.75);ax.set_title(title,pad=26)
        return self.save(fig,name)

    def effect_heatmap(self,effects,metric,name,title):
        selections=list(EFFECT_SELECTIONS);columns=[(c,b) for c in ['all','probability_ge_0.5'] for b in BANDS]
        mat=np.full((len(selections),len(columns)),np.nan)
        for j,(conf,band) in enumerate(columns):
            s=effects[(effects.metric==metric)&(effects.confidence==conf)&(effects.band==band)].set_index('selection').rank_biserial
            mat[:,j]=s.reindex(selections).to_numpy(float)
        fig,ax=plt.subplots(figsize=(7.4,8.2))
        im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        for yy in range(len(selections)):
            for xx in range(len(columns)):
                v=mat[yy,xx]
                if np.isfinite(v):ax.text(xx,yy,f'{v:.2f}',ha='center',va='center',fontsize=8,color='white' if abs(v)>.55 else '#1f2933')
        ax.set_xticks(range(len(columns)),[b for _,b in columns]);ax.set_yticks(range(len(selections)),list(EFFECT_SELECTIONS.values()));ax.grid(False)
        ax.axvline(2.5,color='white',lw=3);ax.axhline(9.5,color='white',lw=3)
        sec=ax.secondary_xaxis('top');sec.set_xticks([1,4],['All valid','Top probability ≥ 0.5\n(DECaLS10 rows only)'],fontsize=8.5);sec.tick_params(length=0)
        ax.set_xlabel('Photometric band')
        fig.colorbar(im,ax=ax,label='Rank-biserial separation (class vs others)',shrink=.75);ax.set_title(title,pad=30)
        return self.save(fig,name)

    def human_class_histogram(self,frame,metric,band,limits,name,title):
        fig,ax=plt.subplots(figsize=(6.2,4.2));bins=np.linspace(*limits,45);d=frame[frame.band==band]
        for j,(sel,label) in enumerate(HUMAN_SELECTIONS.items()):
            v=d.loc[d[sel].eq(True),metric].to_numpy(float);v=v[np.isfinite(v)]
            if len(v):ax.hist(v,bins=bins,weights=np.full(len(v),100/len(v)),histtype='step',lw=1.7,color=CLASS_COLORS[j],ls=CLASS_STYLES[j],label=f'{label} (n={len(v):,})')
            self.audit(name,metric+' | '+label,v,limits,band)
        if limits[0]<0<limits[1]:ax.axvline(0,color='#6b7280',lw=.7)
        ax.set(xlim=limits,xlabel=LABELS[metric],ylabel='Share of class sample per bin (%)',title=title);ax.legend(frameon=False,fontsize=8)
        return self.save(fig,name)

    def probability_trend(self,trends,target,metric,ylimits,name,title):
        fig,ax=plt.subplots(figsize=(6.2,4.2));part=trends[(trends.target==target)&(trends.metric==metric)]
        for band in BANDS:
            d=part[part.band==band].sort_values('lower')
            ax.plot(d.probability_median,d['median'],'o',ls=BAND_STYLES[band],ms=4,lw=1.7,color=BAND_COLORS[band],label=f'{band} band (n={int(d.n.sum()):,})')
            ax.fill_between(d.probability_median,d.p16,d.p84,color=BAND_COLORS[band],alpha=.09,lw=0)
        if ylimits[0]<0<ylimits[1]:ax.axhline(0,color='#6b7280',lw=.7)
        ax.set(xlim=(-.02,1.02),ylim=ylimits,xlabel=('Human answer fraction' if target in HUMAN_LABELS else 'DECaLS10 model probability')+' — '+TARGETS[target].split(': ',1)[1],
               ylabel=LABELS[metric],title=title);ax.legend(frameon=False,fontsize=8.5)
        return self.save(fig,name)


def binned_trends(joined,metric):
    """Same binning as paper_plots.PlotBook.relations (v2): 16 common x bins, >=30 pairs."""
    selected=joined[joined.gf_unflagged];rows=[]
    for prop in PROPERTIES:
        log=prop in LOG_PROPERTIES;xl=robust_limits([selected[prop]],positive=log)
        edges=np.geomspace(*xl,17) if log else np.linspace(*xl,17)
        for band in BANDS:
            part=selected[selected.band==band][[prop,metric]].replace([np.inf,-np.inf],np.nan).dropna()
            bins=np.digitize(part[prop],edges)-1
            for i in range(16):
                vals=part.loc[bins==i,metric].to_numpy()
                q=np.quantile(vals,[.16,.5,.84]) if len(vals)>=30 else [np.nan]*3
                rows.append(dict(property=prop,metric=metric,band=band,left=edges[i],right=edges[i+1],n=len(vals),p16=q[0],median=q[1],p84=q[2]))
    return pd.DataFrame(rows)


def embedding_cell_medians(task):
    """Worker (revision 6): the finite per-cell medians that embedding_panel colours, without drawing."""
    coords,values=task;x,y=coords;ok=np.isfinite(values);extent=(x.min(),x.max(),y.min(),y.max())
    fig,ax=plt.subplots()
    h=ax.hexbin(x[ok],y[ok],C=values[ok],reduce_C_function=np.median,gridsize=70,mincnt=3,extent=extent)
    cells=np.asarray(h.get_array(),dtype=float);plt.close(fig);return cells[np.isfinite(cells)]


def embedding_panel(task):
    """Worker: one median-colour hexbin map of one metric on one saved embedding."""
    path,coords,values,limits,label,title,dpi=task
    style();fig,ax=plt.subplots(figsize=(5.4,4.6))
    x,y=coords;ok=np.isfinite(values);extent=(x.min(),x.max(),y.min(),y.max())
    # Grey layer: every saved point, so cells without enough finite values stay visible.
    ax.hexbin(x,y,gridsize=70,mincnt=1,extent=extent,cmap=matplotlib.colors.ListedColormap(['#d9dde1']),linewidths=0,rasterized=True)
    # Clipping is monotone, so it only saturates colour at the display-window edges.
    h=ax.hexbin(x[ok],y[ok],C=np.clip(values[ok],*limits),reduce_C_function=np.median,gridsize=70,mincnt=3,extent=extent,
                vmin=limits[0],vmax=limits[1],cmap='cividis',linewidths=0,rasterized=True)
    fig.colorbar(h,ax=ax,label='Median '+label+' per cell',shrink=.85)
    ax.set(xticks=[],yticks=[],xlabel='Embedding axis 1 (dimensionless)',ylabel='Embedding axis 2 (dimensionless)',title=title);ax.grid(False)
    return str(write_png(fig,path,dpi))
