"""Readable class distributions and probability trends on exact finite pairs."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .classification_data import HUMAN_LABELS,CLASS_LABELS,class_effect
from .paper_plots import COLORS

HUMAN_QUERIES={
 'analysis-morphology-etg-min2':'Human ETG / elliptical',
 'analysis-morphology-ltg-min2':'Human LTG / spiral',
 'analysis-morphology-featureless-min2':'Human Featureless',
 'analysis-morphology-irregular-min2':'Human Not sure (Irr/other)',
 'analysis-visible-nucleus-visible-min2':'Human visible nucleus',
 'analysis-visible-nucleus-not-visible-min2':'Human nucleus not visible',
 'analysis-is-lsb-yes-min2':'Human LSB',
 'analysis-is-lsb-no-min2':'Human Non-LSB',
 'analysis-failed-fitting-flagged-min2':'Human failed fitting',
 'analysis-failed-fitting-good-fit-min2':'Human good fitting',
}


def target_label(target):
    return HUMAN_LABELS.get(target,'DECaLS10: '+CLASS_LABELS[int(target.split('_')[1])] if target.startswith('class_') else target)


def query_denominator(frame,query):
    if 'morphology' in query:return frame.human_vote_count.ge(2)
    if 'nucleus' in query:return frame.p_visible_nucleus_count.ge(2)
    if 'failed-fitting' in query:return frame.p_failed_fitting_count.ge(2)
    return frame.p_lsb_count.ge(2)


def effects(frame,metrics):
    rows=[]
    for band,d in frame.groupby('band'):
        for query in [*HUMAN_QUERIES,*[f'zoobot_top_{i}' for i in range(10)]]:
            if query.startswith('zoobot'):
                mask=d.zoobot_valid.eq(True);selected=d.zoobot_top.eq(int(query.split('_')[-1])).fillna(False)
            else:mask=query_denominator(d,query);selected=d[query].eq(True)
            for confidence in ['all','probability_ge_0.5'] if query.startswith('zoobot') else ['all']:
                eligible=mask & (d.zoobot_top_probability.ge(.5) if confidence!='all' else True)
                for metric in metrics:
                    vals=d.loc[eligible,metric].to_numpy(float);sel=selected[eligible].to_numpy(bool)
                    n1,n0,e=class_effect(vals,sel)
                    a=vals[sel & np.isfinite(vals)]
                    q=np.quantile(a,[.16,.5,.84]) if len(a) else [np.nan]*3
                    rows.append(dict(band=band,selection=query,confidence=confidence,metric=metric,n_class=n1,n_other=n0,
                                     rank_biserial=e,p16=q[0],median=q[1],p84=q[2]))
    return pd.DataFrame(rows)


def create_plots(frame,assoc,contract,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    metrics=contract['metrics'];labels=contract['labels'];limits=contract['limits'];items=[];trends=[];audit=[]
    plt.rcParams.update({'font.size':10,'mathtext.fontset':'stix','axes.spines.top':False,'axes.spines.right':False})
    def save(fig,name,title,caption):
        path=output/(name+'.png');fig.savefig(path,dpi=130,bbox_inches='tight');plt.close(fig)
        items.append(dict(name=name,title=title,file=str(path),caption=caption))
    for cohort,stat in [('available','spearman'),('human_zoobot_matched','spearman'),('galfit_adjusted','partial_rank'),('three_plus_human_votes','spearman')]:
        fig,axes=plt.subplots(3,2,figsize=(22,19),constrained_layout=True)
        for row,band in enumerate(['g','r','i']):
            for col,targets in enumerate([list(HUMAN_LABELS),[f'class_{i}' for i in range(10)]]):
                ax=axes[row,col];part=assoc[(assoc.band==band)&(assoc.cohort==cohort)]
                matrix=part.pivot(index='metric',columns='target',values=stat).reindex(index=metrics,columns=targets)
                n=part.pivot(index='metric',columns='target',values='n').reindex(index=metrics,columns=targets)
                art=ax.imshow(matrix.to_numpy(),vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
                for iy in range(len(metrics)):
                    for ix in range(len(targets)):
                        val=matrix.iloc[iy,ix]
                        if np.isfinite(val):ax.text(ix,iy,f'{val:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(val)>.65 else 'black')
                ax.set_xticks(range(len(targets)),[target_label(t).replace('DECaLS10: ','') for t in targets],rotation=40,ha='right')
                ax.set_yticks(range(len(metrics)),[labels[m] for m in metrics],fontsize=10)
                ax.set_title(f'{band} band · '+('Human answer fractions' if col==0 else 'DECaLS10 probabilities')+f'\nfinite-pair n: {int(n.min().min()):,}–{int(n.max().max()):,}')
                fig.colorbar(art,ax=ax,shrink=.7,label='Partial rank correlation' if stat=='partial_rank' else r'Spearman rank correlation $\rho$')
        save(fig,'correlations_'+cohort,cohort.replace('_',' ').capitalize()+' probability associations',
             'Numbers are finite-pair descriptive rank correlations, not accuracies. Human fractions require ≥2 responses to that question. '+
             ('Linear residuals of ranks adjusted for GALFIT mean effective surface brightness, effective radius and axis ratio; same finite control cohort.' if stat=='partial_rank' else 'Each cell uses all eligible finite pairs; no display-window cut.'))
    effect=effects(frame,metrics);effect.to_csv(output.parent/'class_effects.csv',index=False)
    fig,axes=plt.subplots(3,2,figsize=(24,20),constrained_layout=True)
    for row,band in enumerate(['g','r','i']):
        for col,selections in enumerate([list(HUMAN_QUERIES),[f'zoobot_top_{i}' for i in range(10)]]):
            ax=axes[row,col];part=effect[(effect.band==band)&(effect.confidence=='all')]
            vals=part.pivot(index='metric',columns='selection',values='rank_biserial').reindex(index=metrics,columns=selections)
            art=ax.imshow(vals,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
            for iy in range(len(metrics)):
                for ix in range(len(selections)):
                    v=vals.iloc[iy,ix]
                    if np.isfinite(v):ax.text(ix,iy,f'{v:.2f}',ha='center',va='center',fontsize=9,color='white' if abs(v)>.65 else 'black')
            names=[HUMAN_QUERIES[s] if col==0 else CLASS_LABELS[int(s.split('_')[-1])] for s in selections]
            ax.set_xticks(range(len(selections)),names,rotation=40,ha='right');ax.set_yticks(range(len(metrics)),[labels[m] for m in metrics])
            ax.set_title(f'{band} band · '+('Published human consensus selections' if col==0 else 'DECaLS10 largest-probability class'))
            fig.colorbar(art,ax=ax,shrink=.7,label='Rank-biserial effect (class versus eligible others)')
    save(fig,'class_effects','Selected classes: distributional separation','Effect = P(class value > other) − P(class value < other), including ties through the Mann–Whitney statistic. Positive means larger morphology values in the named class. Human contrasts use only galaxies with ≥2 answers for that question; other includes split outcomes. Zoobot top-class contrasts use valid prediction vectors. Some top classes are rare; exact finite sample sizes and ≥0.5 confidence sensitivity are in the companion table.')
    for band in ['g','r','i']:
        d=frame[frame.band==band];fig,axes=plt.subplots(3,4,figsize=(19,13),constrained_layout=True)
        selections=list(HUMAN_QUERIES)[:4]+['analysis-undecidable-morphology'];names=[HUMAN_QUERIES[s].replace('Human ','') for s in selections[:4]]+['Split consensus']
        for ax,metric in zip(axes.flat,metrics):
            bins=np.linspace(*limits[metric],45)
            for j,(sel,name) in enumerate(zip(selections,names)):
                v=d.loc[d[sel].eq(True),metric].to_numpy(float);v=v[np.isfinite(v)]
                if len(v):ax.hist(v,bins=bins,weights=np.full(len(v),100/len(v)),histtype='step',lw=1.5,color=plt.get_cmap('tab10')(j),label=f'{name} (n={len(v):,})')
                audit.append(dict(figure='human_distributions_'+band,metric=metric,selection=sel,finite=len(v),below=int((v<bins[0]).sum()),above=int((v>bins[-1]).sum()),lower=bins[0],upper=bins[-1]))
            ax.set(xlabel=labels[metric],ylabel='Finite class sample per bin (%)',xlim=limits[metric]);ax.grid(alpha=.15);ax.legend(fontsize=7)
        for ax in axes.flat[len(metrics):]:ax.axis('off')
        fig.suptitle(f'{band} band · Human morphology consensus (≥2 responses) · All catalog candidates',fontsize=15)
        save(fig,'human_distributions_'+band,f'{band} band: human morphology distributions','Published weighted-consensus memberships are used unchanged. Four classes and split outcomes are distinct; unclassified objects are not coded as any class. Histograms use shared robust display windows and normalize by the full finite class sample; off-screen counts are downloadable.')
    for target in [*HUMAN_LABELS,*[f'class_{i}' for i in range(10)]]:
        fig,axes=plt.subplots(3,4,figsize=(19,13),constrained_layout=True)
        edges=np.linspace(0,1,11)
        for ax,metric in zip(axes.flat,metrics):
            for band in ['g','r','i']:
                d=frame[(frame.band==band)&frame[target].notna()&np.isfinite(frame[metric])];xs=[];med=[];low=[];high=[]
                for j in range(10):
                    mask=d[target].ge(edges[j])&(d[target].le(edges[j+1]) if j==9 else d[target].lt(edges[j+1]))
                    vals=d.loc[mask,metric].to_numpy(float)
                    if len(vals)<20:continue
                    qs=np.quantile(vals,[.16,.5,.84]);xs.append(d.loc[mask,target].median());low.append(qs[0]);med.append(qs[1]);high.append(qs[2])
                    trends.append(dict(target=target,band=band,metric=metric,lower=edges[j],upper=edges[j+1],n=len(vals),probability_median=xs[-1],p16=qs[0],median=qs[1],p84=qs[2]))
                ax.plot(xs,med,'o-',ms=3,color=COLORS[band],label=f'{band} (n={len(d):,})');ax.fill_between(xs,low,high,color=COLORS[band],alpha=.10)
            ax.set(xlabel='Human answer fraction' if target in HUMAN_LABELS else 'DECaLS10 model probability',ylabel=labels[metric],xlim=(-.02,1.02),ylim=limits[metric]);ax.grid(alpha=.15);ax.legend(fontsize=8)
        for ax in axes.flat[len(metrics):]:ax.axis('off')
        fig.suptitle(target_label(target)+' · All eligible catalog candidates',fontsize=15)
        save(fig,'trend_'+target,target_label(target)+': morphology trends','Lines show within-bin medians; shading is the 16–84% object spread, not uncertainty on the median. Common probability bins of width 0.1, minimum 20 finite pairs. All finite morphology values contribute before the robust y-axis display window. Human vote fractions are discrete, often based on only two or three votes; they are not calibrated probabilities. Missing bins do not imply zero values.')
    pd.DataFrame(trends).to_csv(output.parent/'probability_trends.csv',index=False)
    pd.DataFrame(audit).to_csv(output.parent/'class_distribution_range_audit.csv',index=False)
    return items
