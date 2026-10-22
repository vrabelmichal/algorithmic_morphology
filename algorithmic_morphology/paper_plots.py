"""Scientific presentation of fixed measurements; plotting crops never alter data."""
from pathlib import Path
import base64
from io import BytesIO
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .report_data import BANDS,finite,robust_limits
from .metrics import PROFILES

COLORS={'g':'#14806c','r':'#be4c3b','i':'#6a58a3'}
L='S_'+PROFILES[0].name
K='S_'+PROFILES[1].name
METRICS=['C','A_raw','A_corrected',L+'_raw',L+'_corrected',K+'_raw',K+'_corrected','Gini_aperture_v1','M20_aperture_v1']
LABELS={
 'C':r'Concentration $C=5\log_{10}(r_{80}/r_{20})$',
 'A_raw':r'Raw asymmetry $A_{\rm raw}$',
 'A_sky':r'Sky asymmetry correction $A_{\rm sky}$',
 'A_corrected':r'Sky-corrected asymmetry $A_{\rm corr}$',
 L+'_raw':r'Lotz-like clumpiness $S_{\rm L,raw}$',
 L+'_corrected':r'Lotz-like clumpiness $S_{\rm L,corr}$',
 K+'_raw':r'Conselice-like clumpiness $S_{\rm C,raw}$',
 K+'_corrected':r'Conselice-like clumpiness $S_{\rm C,corr}$',
 'Gini_aperture_v1':r'Aperture Gini coefficient $G$',
 'M20_aperture_v1':r'Brightest-20% moment $M_{20}$',
 'axis_ratio_positive_aperture_v1':r'Positive-flux axis ratio $q_{+}$',
 'ellipticity_positive_aperture_v1':r'Positive-flux ellipticity $1-q_{+}$',
 'position_angle_positive_aperture_v1':r'Image major-axis angle $\theta_{+}$ (degrees)',
 'gf_mu_mean':'Mean effective surface brightness\n'+r'$\langle\mu\rangle_e$ (mag arcsec$^{-2}$)',
 'gf_mu_e':r'Surface brightness at $R_e$, $\mu_e$ (mag arcsec$^{-2}$)',
 'gf_re_arcsec':r'GALFIT effective radius $R_e$ (arcsec)',
 'gf_re_pixels':r'GALFIT effective radius $R_e$ (pixels)',
 'gf_mag':r'GALFIT total apparent magnitude $m$ (mag)',
 'gf_n':r'GALFIT S\'ersic index $n$',
 'gf_q':r'GALFIT minor/major axis ratio $q$',
 'gf_chi2red':r'Stored GALFIT reduced $\chi^2$',
 'rpetro':r'Measured Petrosian radius $R_{\rm P}$ (pixels)',
}
# Avoid Python raw-string apostrophe ambiguity in the display name.
LABELS['gf_n']=r'GALFIT Sersic index $n$'
PROPERTIES=['gf_mu_mean','gf_re_arcsec','gf_mag','gf_n','gf_q','gf_chi2red']
SHORT={key:label.split('$')[0].strip() for key,label in LABELS.items()}
SHORT['gf_chi2red']='Stored GALFIT reduced chi-square'


class PlotBook:
    def __init__(self,root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.records=[];self.figures=[]
        plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':11,
                             'font.family':'DejaVu Sans','mathtext.fontset':'stix',
                             'axes.spines.top':False,'axes.spines.right':False,
                             'legend.fontsize':9,'figure.facecolor':'white','savefig.facecolor':'white'})

    def save(self,fig,name,title,caption,full=None):
        # Files also support export into documents independently of the HTML.
        path=self.root/(name+'.png');fig.savefig(path,dpi=150,bbox_inches='tight');plt.close(fig)
        item=dict(name=name,title=title,caption=caption,file=str(path),full_file=full,range_rows=[r for r in self.records if r['figure']==name])
        self.figures.append(item);return item

    def audit(self,figure,series,values,limits,band=None,axis='x'):
        v=finite(values);valid=v[np.isfinite(v)];lo,hi=limits
        self.records.append(dict(figure=figure,series=series,band=band,axis=axis,finite=len(valid),
                                 lower=lo,upper=hi,below=int((valid<lo).sum()),above=int((valid>hi).sum()),
                                 shown=int(((valid>=lo)&(valid<=hi)).sum())))

    def histogram(self,ax,frame,column,limits,figure,*,scale=1.,label=None,band_colors=True):
        bins=np.linspace(*limits,61)
        for band in BANDS:
            vals=finite(frame.loc[frame.band==band,column])/scale;vals=vals[np.isfinite(vals)]
            if len(vals):
                ax.hist(vals,bins=bins,weights=np.full(len(vals),100/len(vals)),histtype='step',lw=1.6,
                        color=COLORS[band],label=f'{band} band (n={len(vals):,})')
            self.audit(figure,column+('/10' if scale==10 else ''),vals,limits,band)
        ax.set_xlim(*limits);ax.set_xlabel(label or LABELS[column]);ax.set_ylabel('Finite sample per bin (%)')
        ax.grid(alpha=.18);ax.legend(loc='best',frameon=False)

    def grouped_histograms(self,frame,columns,name,title,limits=None):
        figures=[]
        for mode in ['zoom','full']:
            fig,axes=plt.subplots(1,len(columns),figsize=(5*len(columns),3.9),squeeze=False,constrained_layout=True)
            for ax,column in zip(axes.flat,columns):
                lim=(limits or {}).get(column) or robust_limits([frame[column]],quantiles=(0,1) if mode=='full' else (.005,.995))
                self.histogram(ax,frame,column,lim,name+'_'+mode)
            caption='Bands share each property’s bin edges and limits. Percentages use all finite values, including off-screen tails; the cropped histogram need not sum to 100%. '
            caption+= 'Main view uses pooled 0.5th–99.5th percentiles plus 3.5% padding. Full-range view includes all finite extrema. Exact per-band tail counts are in the plot-range audit.'
            item=self.save(fig,name+'_'+mode,title,caption);figures.append(item)
        figures[0]['full_file']=figures[1]['file'];return figures[0]

    def asymmetry_overlay(self,frame):
        lim=robust_limits([frame.A_raw,frame.A_corrected]);fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
        for ax,band in zip(axes,BANDS):
            part=frame[frame.band==band]
            for col,style,color in [('A_raw','-','#315a86'),('A_corrected','--','#d37b21')]:
                vals=finite(part[col]);vals=vals[np.isfinite(vals)];bins=np.linspace(*lim,71)
                ax.hist(vals,bins=bins,weights=np.full(len(vals),100/max(1,len(vals))),histtype='step',lw=1.8,linestyle=style,color=color,label=f'{"Raw" if col=="A_raw" else "Corrected"} (n={len(vals):,})')
                self.audit('asymmetry_overlay',col,vals,lim,band)
            ax.axvline(0,color='#777',lw=.8);ax.set(xlim=lim,title=f'{band} band',xlabel=r'Asymmetry $A$ (dimensionless)',ylabel='Finite sample per bin (%)');ax.legend(frameon=False);ax.grid(alpha=.18)
        return self.save(fig,'asymmetry_overlay','Raw and sky-corrected asymmetry on a common scale',
                         'Identical bins and limits across all three panels. Each line uses its own finite sample; raw/corrected counts can differ. These are distribution comparisons, not paired difference estimates.')

    def clumpiness(self,frame):
        results=[]
        for mode in ['zoom','full']:
            fig,axes=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
            for row,role in enumerate(['raw','corrected']):
                lim=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10],quantiles=(0,1) if mode=='full' else (.005,.995))
                self.histogram(axes[row,0],frame,L+'_'+role,lim,'clumpiness_'+mode)
                self.histogram(axes[row,1],frame,K+'_'+role,tuple(v*10 for v in lim),'clumpiness_'+mode)
                axes[row,0].set_title(f'{"Raw" if role=="raw" else "Sky-corrected"}: Lotz-like, factor 1')
                axes[row,1].set_title(f'{"Raw" if role=="raw" else "Sky-corrected"}: Conselice-like, factor 10')
            item=self.save(fig,'clumpiness_'+mode,'Both clumpiness definitions and both correction states',
                'Columns preserve each profile’s stored units. Conselice-axis limits and bin widths are exactly ten times the Lotz limits within each row, using their pooled normalization-adjusted quantiles. Raw and corrected rows use separate ranges. All three bands share limits within a panel; omitted tails are counted in the audit.');results.append(item)
        results[0]['full_file']=results[1]['file']
        fig,axes=plt.subplots(2,3,figsize=(14,7.5),constrained_layout=True)
        for row,role in enumerate(['raw','corrected']):
            lim=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10]);bins=np.linspace(*lim,61)
            for ax,band in zip(axes[row],BANDS):
                part=frame[frame.band==band]
                for prefix,scale,color,text in [(L,1,'#315a86','Lotz-like'),(K,10,'#d37b21','Conselice-like / 10')]:
                    vals=finite(part[prefix+'_'+role])/scale;vals=vals[np.isfinite(vals)]
                    ax.hist(vals,bins=bins,weights=np.full(len(vals),100/max(1,len(vals))),histtype='step',lw=1.5,color=color,label=text)
                    self.audit('clumpiness_overlay',prefix+'_'+role+('/10' if scale==10 else ''),vals,lim,band)
                ax.set(xlim=lim,title=f'{band} band · {"raw" if role=="raw" else "sky-corrected"}',xlabel=r'Clumpiness $S_{\rm L}$ or $S_{\rm C}/10$',ylabel='Finite sample per bin (%)');ax.axvline(0,color='#777',lw=.7);ax.grid(alpha=.18);ax.legend(frameon=False)
        overlay=self.save(fig,'clumpiness_overlay','Clumpiness profiles after removing only the normalization factor',
            'The orange curve divides the stored Conselice-like values by ten for display only. Smoothing scale and central exclusion still differ. Both stored profiles and their negative corrected values are unchanged.')
        return results[0],overlay

    def density(self,frame,x,y,name,title,xlimits=None,ylimits=None,div_y=1):
        xl=xlimits or robust_limits([frame[x]]);yl=ylimits or robust_limits([frame[y]/div_y])
        fig,axes=plt.subplots(1,3,figsize=(14,4),constrained_layout=True)
        for ax,band in zip(axes,BANDS):
            part=frame[frame.band==band][[x,y]].replace([np.inf,-np.inf],np.nan).dropna()
            a=part[x].to_numpy();b=part[y].to_numpy()/div_y
            inside=(a>=xl[0])&(a<=xl[1])&(b>=yl[0])&(b<=yl[1])
            if inside.any():
                h=ax.hexbin(a[inside],b[inside],gridsize=48,extent=(*xl,*yl),mincnt=1,bins='log',cmap='cividis',rasterized=True)
                fig.colorbar(h,ax=ax,label='Objects per cell (log color scale)',shrink=.8)
            rho=part[x].corr(part[y],method='spearman')
            ax.set(xlim=xl,ylim=yl,xlabel=LABELS[x],ylabel=LABELS[y]+(' / 10' if div_y==10 else ''),title=f'{band}: paired n={len(part):,}\nshown={inside.sum():,}; rank correlation={rho:.3f}')
            if x=='A_raw' or (x.startswith('S_') and y.startswith('S_')) or x=='gf_q':
                lo=max(xl[0],yl[0]);hi=min(xl[1],yl[1]);ax.plot([lo,hi],[lo,hi],'--',color='#be4c3b',lw=1)
            self.audit(name,x,a,xl,band,'x');self.audit(name,y+('/10' if div_y==10 else ''),b,yl,band,'y')
            ax.grid(alpha=.12)
        return self.save(fig,name,title,'Only paired finite values are used. Panel titles give the full paired count, displayed count and Spearman rank correlation calculated before display cropping. Limits are shared across bands. Dashed equality lines appear only where quantities are comparable.')

    def relations(self,joined):
        selected=joined[joined.gf_unflagged];result=[];trend_rows=[]
        for prop in PROPERTIES:
            log=prop in ['gf_re_arcsec','gf_chi2red']
            xl=robust_limits([selected[prop]],positive=log)
            edges=np.geomspace(*xl,17) if log else np.linspace(*xl,17)
            centers=np.sqrt(edges[:-1]*edges[1:]) if log else (edges[:-1]+edges[1:])/2
            fig,axes=plt.subplots(3,3,figsize=(14,11),constrained_layout=True)
            for ax,metric in zip(axes.flat,METRICS):
                yl=robust_limits([joined[metric]])
                for band in BANDS:
                    part=selected[selected.band==band][[prop,metric]].replace([np.inf,-np.inf],np.nan).dropna()
                    bins=np.digitize(part[prop],edges)-1
                    med=[];low=[];high=[]
                    for i in range(16):
                        vals=part.loc[bins==i,metric].to_numpy()
                        q=np.quantile(vals,[.16,.5,.84]) if len(vals)>=30 else [np.nan]*3
                        low.append(q[0]);med.append(q[1]);high.append(q[2])
                        trend_rows.append(dict(property=prop,metric=metric,band=band,left=edges[i],right=edges[i+1],n=len(vals),p16=q[0],median=q[1],p84=q[2]))
                    ax.plot(centers,med,color=COLORS[band],lw=1.7,label=f'{band} (paired n={len(part):,})')
                    ax.fill_between(centers,low,high,color=COLORS[band],alpha=.09)
                    self.audit('relations_'+prop,prop,part[prop],xl,band,'x')
                    self.audit('relations_'+prop,metric,part[metric],yl,band,'y')
                ax.set(xlim=xl,ylim=yl,xlabel=LABELS[prop],ylabel=LABELS[metric]);ax.grid(alpha=.18)
                if log:ax.set_xscale('log')
                ax.legend(frameon=False,fontsize=8,loc='best')
            result.append(self.save(fig,'relations_'+prop,'Morphology versus '+SHORT[prop],
                'Main GALFIT cohort: successful exports with finite positive radius/index/axis ratio/reduced chi-square and no x, y, magnitude, radius, index or axis-ratio warning. Lines are medians in 16 common x bins; shading is the 16th–84th percentile spread, not an uncertainty on the median. Bins with fewer than 30 pairs have gaps. All finite y values contribute to summaries before display cropping. Brighter magnitude/surface brightness is toward smaller x. Radius and reduced chi-square use logarithmic x axes. Band-specific paired counts are shown in every legend.'))
        return result,pd.DataFrame(trend_rows)


def figure_html(item,number):
    def img(path):return 'data:image/png;base64,'+base64.b64encode(Path(path).read_bytes()).decode()
    full=''
    if item.get('full_file'):
        full=f'<details><summary>Show full finite range, including outliers</summary><img loading="lazy" src="{img(item["full_file"])}" alt="Full-range companion for figure {number}"></details>'
    from html import escape
    rows=item.get('range_rows',[])
    audit=''
    if rows:
        table=pd.DataFrame(rows).drop(columns=['figure'])
        table['series']=table['series'].map(lambda s:LABELS.get(s.removesuffix('/10'),s)+(' / 10' if s.endswith('/10') else ''))
        # Plain descriptions in HTML; mathematical typesetting remains on axes.
        table['series']=table['series'].str.replace(r'\$[^$]*\$','',regex=True).str.strip()
        audit='<details><summary>Display limits and off-screen tail counts</summary>'+table.to_html(index=False,float_format=lambda x:f'{x:.5g}')+'</details>'
    return f'<figure id="fig-{item["name"]}"><h4>Figure {number}. {escape(item["title"])}</h4><img loading="lazy" src="{img(item["file"])}" alt="{escape(item["title"])}"><figcaption>{escape(item["caption"])}</figcaption>{audit}{full}<p class="plot-export"><a href="paper_v2_assets/{item["name"]}.png">Open full-resolution PNG</a></p></figure>'
