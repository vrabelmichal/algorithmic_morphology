"""Build revision 2: a navigable scientific report from immutable morphology results.

Reuses fixed measurements; never runs CAS, changes profiles or overwrites source data.
"""
import argparse
import base64
from datetime import datetime,timezone
from html import escape
from io import BytesIO
import json
import logging
import os
from pathlib import Path
import shutil
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from algorithmic_morphology.report import sha256_file
from algorithmic_morphology.report_data import BANDS,GF_COLUMNS,galfit_long,join_galfit,relation_stats,finite,robust_limits
from algorithmic_morphology.paper_plots import PlotBook,METRICS,PROPERTIES,LABELS,L,K,figure_html,COLORS

LOG=logging.getLogger('paper_report')
PRIMARY='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_256__per_band_mag_reff__xy_2_5arcsec'
ALTERNATE='sky_subtracted_2_onmr2__detr_and_vit_all_multiband_256__per_band_magnitude_optimization__r_constant_offset'


def table_html(frame):
    return '<div class="table-scroll">'+frame.to_html(index=False,float_format=lambda x:f'{x:.4g}',na_rep='unavailable',border=0,classes='data-table')+'</div>'


def short_label(key):
    names={'C':'Concentration','A_raw':'Raw asymmetry','A_corrected':'Corrected asymmetry',L+'_raw':'Lotz-like S: raw',L+'_corrected':'Lotz-like S: corrected',K+'_raw':'Conselice-like S: raw',K+'_corrected':'Conselice-like S: corrected','Gini_aperture_v1':'Aperture Gini','M20_aperture_v1':'Aperture M20',
           'gf_mu_mean':'Mean effective surface brightness','gf_re_arcsec':'Effective radius','gf_mag':'Total apparent magnitude','gf_n':'Sersic index','gf_q':'GALFIT axis ratio','gf_chi2red':'Reduced chi-square'}
    return names.get(key,key)


def equation(text):
    fig=plt.figure(figsize=(11,.85));fig.text(.02,.45,text,fontsize=17,va='center');buffer=BytesIO();fig.savefig(buffer,format='svg',bbox_inches='tight',transparent=True);plt.close(fig)
    return '<img class="equation" alt="'+escape(text)+'" src="data:image/svg+xml;base64,'+base64.b64encode(buffer.getvalue()).decode()+'">'


def dictionary(columns):
    exact={
      'coadd_object_id':('—','Stable DES catalog identifier; join key together with band.'),
      'band':('g, r, i','Photometric band; measurements and GALFIT columns are matched within the same band.'),
      'status':('category','ok: finite C and raw A; invalid: one or both unavailable; error: input or processing failure. Other metrics can survive an invalid CAS status.'),
      'quality_flags':('text','Pipe-separated CAS diagnostics; the list of flag meanings follows below.'),
      'error_type':('text','Python exception category for error rows; absent for measured rows.'),
      'error_message':('text','Input/processing reason; the full run contains only missing_cutout error rows.'),
      'method_version':('text','Shared CAS definition identifier: cas_discrete_v3 for this run.'),
      'aperture_diagnostics_version':('text','Extended diagnostic definition: aperture_diagnostics_v1.'),
      'extended_quality_flags':('text','Reasons an extended diagnostic is unavailable, including nonpositive flux or second moment.'),
      'rpetro':('pixels','Measured circular Petrosian radius, eta=0.2, about the catalog input center. It is not GALFIT effective radius.'),
      'aperture_radius':('pixels','Outer measurement radius, 1.5 times the measured Petrosian radius.'),
      'r20':('pixels','Radius enclosing 20% of the signed aperture flux about the optimized A center.'),
      'r80':('pixels','Radius enclosing 80% of the signed aperture flux about the optimized A center.'),
      'C':('dimensionless','5 log10(r80/r20); larger values describe more centrally concentrated aperture light.'),
      'A_raw':('dimensionless','Sum of absolute 180-degree residuals / absolute galaxy flux on the identical valid symmetric aperture pixels.'),
      'A_sky':('dimensionless','Mean absolute rotation residual in the chosen blank corner, scaled by valid galaxy pixel count / galaxy absolute flux.'),
      'A_corrected':('dimensionless','A_raw − A_sky; valid negative values are retained.'),
      'mask_fraction':('fraction','Fraction of the entire stamp excluded by the contaminant mask or nonfinite pixels.'),
      'aperture_mask_fraction':('fraction','Fraction of the circular aperture excluded by contaminant/nonfinite pixels.'),
      'asymmetry_valid_pixels':('pixels','Number of original and bilinearly rotated valid pixels used in raw A.'),
      'asymmetry_absolute_flux':('input flux units','Sum |I| on precisely the pixels used by A; the common raw/sky correction denominator.'),
      'aperture_valid_pixels':('pixels','Finite, unmasked pixels inside the aperture used for Gini, M20 and shape.'),
      'Gini_aperture_v1':('dimensionless','Absolute-flux inequality of pixels in our circular aperture; 0 is uniform, 1 is a single bright pixel. Not segmentation Gini.'),
      'M20_aperture_v1':('log10 ratio','Signed second moment of brightest 20% of aperture flux / total signed second moment, about the signed centroid. Nonpositive moments are invalid; signed noise can produce atypical positive ratios.'),
      'axis_ratio_positive_aperture_v1':('dimensionless','Square root of the minor/major covariance-eigenvalue ratio using max(I,0) weights.'),
      'ellipticity_positive_aperture_v1':('dimensionless','1 minus the positive-flux aperture axis ratio.'),
      'position_angle_positive_aperture_v1':('degrees','Major-axis angle from increasing image x toward increasing image y, modulo 180 degrees. Not a celestial PA; unstable for round objects.'),
      'image_unit':('FITS BUNIT','Stored input image unit, where present; no unit is invented for missing BUNIT.'),
      'image_height':('pixels','Input stamp height.'),'image_width':('pixels','Input stamp width.'),
    }
    for axis in ['x','y']:
        exact['center_initial_'+axis]=('zero-based pixels','Catalog coordinate transformed through input FITS WCS; '+axis+' coordinate.')
        exact['center_asym_'+axis]=('zero-based pixels','Bounded raw-A minimizing center (or flagged initial-center fallback); '+axis+' coordinate.')
        exact['centroid_signed_'+axis]=('zero-based pixels','Signed-flux centroid inside the aperture for the extended moment diagnostics; '+axis+' coordinate.')
    for prefix,title,width,inner,scale in [(L,'Lotz/statmorph-like',.25,.25,1),(K,'Conselice-like',.30,.075,10)]:
        for role in ['raw','sky','corrected']:
            desc={'raw':'Positive boxcar residual summed over the declared annulus / signed annular flux.', 'sky':'Mean positive boxcar residual in the selected sky tile, scaled to annular galaxy pixel count / signed annular flux.', 'corrected':'Raw minus sky term; valid negative values are retained.'}[role]
            exact[prefix+'_'+role]=('dimensionless',f'{title}, factor {scale}; smoothing {width} R_P and inner exclusion {inner} R_P. '+desc)
        exact[prefix+'_boxcar_width_pixels']=('pixels',f'Actual integer boxcar width int({width} R_P); widths below two pixels are invalid.')
        exact[prefix+'_inner_radius']=('pixels',f'Central exclusion radius {inner} R_P for the {title} profile.')
    gf={
      'gf_mag':('mag','sersic_mag','Stored total single-Sersic apparent magnitude; the source configuration uses a photometric zeropoint of 30. No additional extinction, redshift or K correction is applied here.'),
      'gf_re_arcsec':('arcsec','sersic_reff_arcsec','Stored semi-major-axis effective (half-model-light) radius. The producer converts fitted pixels to arcsec; this analysis does not substitute it for R_P.'),
      'gf_re_pixels':('pixels','sersic_reff_pixels','Stored fitted semi-major-axis effective radius in pixels.'),
      'gf_mu_mean':('mag arcsec^-2','sersic_mean_mue','Stored mean surface brightness inside the effective ellipse; fainter corresponds to larger numbers. Equals m+2.5 log10(2 pi q R_e^2) for the model.'),
      'gf_mu_e':('mag arcsec^-2','sersic_mue','Stored Sersic surface brightness at R_e, distinct from the mean inside R_e.'),
      'gf_n':('dimensionless','sersic_n','Sersic shape index; shared across bands in this GALFITM configuration, constrained to 0.2–5.'),
      'gf_q':('dimensionless','sersic_q','Model minor/major axis ratio; shared across bands, constrained to 0.3–1.'),
      'gf_chi2red':('dimensionless','sersic_chi2red','Stored CHI2NU/reduced chi-square from the fitter; not recomputed using exported NDOF.'),
    }
    rows=[]
    for col in columns:
        if col not in exact:raise ValueError('Column lacks dictionary entry: '+col)
        unit,desc=exact[col];rows.append({'Column':col,'Units':unit,'Meaning':desc})
    for col,(unit,source,desc) in gf.items():rows.append({'Column':col,'Units':unit,'Meaning':f'{source} + band suffix → {col}. '+desc})
    extras={
      'gf_present':'Object ID is present in the primary GALFIT export.',
      'gf_success':'Producer success flag is True; it means a returned fit, not scientific validation.',
      'gf_warning_count':'Number of nonzero x, y, magnitude, effective-radius, n and q parameter warning flags for this band.',
      'gf_finite_physical':'Finite magnitude, mean/at-radius surface brightness, radius, n, q and reduced chi-square, with radius>0, n>0, 0<q<=1 and chi-square>0.',
      'gf_unflagged':'Main relation-plot cohort: present, success, finite_physical, and warning_count=0.',
      'gf_at_n_bound':'n within absolute 0.001 of configured 0.2 or 5 bounds.',
      'gf_at_q_bound':'q within absolute 0.001 of configured 0.3 or 1 bounds; a round object is not necessarily a failed fit.',
      'gf_interior':'Sensitivity cohort: main unflagged cohort minus n/q boundary cases. This is a selection check, not a claim those excluded fits are wrong.',
    }
    for col,desc in extras.items():rows.append({'Column':col,'Units':'flag/count','Meaning':desc})
    return pd.DataFrame(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report-dir',type=Path,required=True)
    p.add_argument('--source-project',type=Path,default=Path('../nucleated_lsbg'));a=p.parse_args()
    logging.basicConfig(level=logging.INFO,stream=sys.stdout,format='%(asctime)sZ %(levelname)s %(message)s',datefmt='%Y-%m-%dT%H:%M:%S');logging.Formatter.converter=time.gmtime
    root=a.report_dir.resolve();data=root/'paper_v2_data';assets=root/'paper_v2_assets'
    if (root/'report_v2.html').exists():raise FileExistsError('Report v2 already exists; preserve this revision and choose a new version')
    data.mkdir(exist_ok=False);assets.mkdir(exist_ok=False)
    provenance=json.loads((root/'table_provenance.json').read_text());morph_path=root/'measurements_metrics.parquet'
    if sha256_file(morph_path)!=provenance['compact_table_sha256']:raise ValueError('Morphology table changed')
    frame=pd.read_parquet(morph_path);assert len(frame)==252600 and frame.coadd_object_id.nunique()==84200
    sys.path.insert(0,str(a.source_project.resolve()))
    pc=__import__('configs.'+PRIMARY,fromlist=['x']);ac=__import__('configs.'+ALTERNATE,fromlist=['x'])
    gfpath=Path(pc.single_component_fits_table_path).resolve();altpath=Path(ac.single_component_fits_table_path).resolve()
    LOG.info('Loading immutable morphology and primary GALFIT source %s',gfpath)
    gf=pd.read_parquet(gfpath);alt=pd.read_parquet(altpath)
    joined=join_galfit(frame,gf);alternative=galfit_long(alt)
    joined_cols=['coadd_object_id','band']+[c for c in joined if c.startswith('gf_')]
    joined[joined_cols].to_parquet(data/'galfit_join.parquet',index=False)
    coverage=[]
    for band,part in joined.groupby('band'):
        coverage.append(dict(band=band,requested=len(part),matched=int(part.gf_present.sum()),unmatched=int((~part.gf_present).sum()),
                             unflagged=int(part.gf_unflagged.sum()),interior=int(part.gf_interior.sum()),
                             parameter_warning=int((part.gf_present&part.gf_warning_count.gt(0)).sum()),
                             invalid_core=int((part.gf_present&~part.gf_finite_physical).sum()),
                             near_n_boundary=int((part.gf_unflagged&part.gf_at_n_bound).sum()),near_q_boundary=int((part.gf_unflagged&part.gf_at_q_bound).sum())))
    coverage=pd.DataFrame(coverage);coverage.to_csv(data/'galfit_coverage.csv',index=False)
    LOG.info('GALFIT coverage: %s',coverage.to_dict('records'))
    correlations=relation_stats(joined,METRICS+['axis_ratio_positive_aperture_v1','ellipticity_positive_aperture_v1'],PROPERTIES)
    correlations.to_csv(data/'galfit_correlations.csv',index=False)
    sensitivity=joined[['coadd_object_id','band','gf_unflagged']+list(GF_COLUMNS)].merge(alternative,on=['coadd_object_id','band'],suffixes=('_primary','_constant'),validate='one_to_one')
    sensitivity=sensitivity[sensitivity.gf_unflagged_primary&sensitivity.gf_unflagged_constant]
    sensrows=[]
    for band,part in sensitivity.groupby('band'):
        for col in ['gf_mag','gf_re_arcsec','gf_mu_mean','gf_n','gf_q']:
            delta=part[col+'_primary']-part[col+'_constant'];sensrows.append(dict(band=band,property=col,paired=len(part),median_primary_minus_constant=float(delta.median()),p16=float(delta.quantile(.16)),p84=float(delta.quantile(.84))))
    sensitivity_table=pd.DataFrame(sensrows);sensitivity_table.to_csv(data/'fit_configuration_sensitivity.csv',index=False)
    unflagged=joined[joined.gf_unflagged]
    residual=unflagged.gf_mag+2.5*np.log10(2*np.pi*unflagged.gf_q*unflagged.gf_re_arcsec**2)-unflagged.gf_mu_mean
    source_meta=dict(primary=dict(config=PRIMARY,path=str(gfpath),sha256=sha256_file(gfpath),modified_utc=datetime.fromtimestamp(gfpath.stat().st_mtime,timezone.utc).isoformat(),rows=len(gf),unique_ids=gf.coadd_object_id.nunique()),
                     alternate=dict(config=ALTERNATE,path=str(altpath),sha256=sha256_file(altpath),modified_utc=datetime.fromtimestamp(altpath.stat().st_mtime,timezone.utc).isoformat(),rows=len(alt),unique_ids=alt.coadd_object_id.nunique()),
                     band_suffixes={'g':'','r':'_r','i':'_i'},source_columns=GF_COLUMNS,
                     mean_surface_brightness_formula_max_abs_residual=float(residual.abs().max()),
                     config_sha256={str(Path(c.__file__).resolve()):sha256_file(c.__file__) for c in [pc,ac]},
                     upstream_code_sha256={str(p):sha256_file(p) for p in [a.source_project/'galfit_functions/galfit_output_reading.py',a.source_project/'functions/statistics.py']},
                     main_cohort='success; finite core; positive radius/n/q/reduced-chi-square; q<=1; no x/y/mag/Re/n/q warning',
                     interior_cohort='main cohort excluding n or q within 0.001 of configured bounds')
    (data/'galfit_sources.json').write_text(json.dumps(source_meta,indent=2))
    dictionary_table=dictionary(frame.columns);dictionary_table.to_csv(data/'column_dictionary.csv',index=False)
    LOG.info('Rendering grouped robust distributions and full-range companions')
    book=PlotBook(assets)
    structure=book.grouped_histograms(frame,['C','Gini_aperture_v1','M20_aperture_v1'],'structure','Concentration and aperture light-distribution diagnostics')
    asym=book.grouped_histograms(frame,['A_raw','A_corrected'],'asymmetry','Raw and sky-corrected asymmetry by band')
    asym_overlay=book.asymmetry_overlay(frame)
    clump,clump_overlay=book.clumpiness(frame)
    shape=book.grouped_histograms(frame,['axis_ratio_positive_aperture_v1','ellipticity_positive_aperture_v1','position_angle_positive_aperture_v1'],'shape','Positive-flux aperture shape',limits={'axis_ratio_positive_aperture_v1':(0,1),'ellipticity_positive_aperture_v1':(0,1),'position_angle_positive_aperture_v1':(0,180)})
    joint=book.density(frame,'Gini_aperture_v1','M20_aperture_v1','gini_m20','Joint aperture Gini and brightest-light moment')
    paired_a=book.density(frame,'A_raw','A_corrected','asymmetry_paired','Paired raw and corrected asymmetry')
    paired_s=[]
    for role in ['raw','corrected']:
        lim=robust_limits([frame[L+'_'+role],frame[K+'_'+role]/10])
        paired_s.append(book.density(frame,L+'_'+role,K+'_'+role,'clumpiness_paired_'+role,'Paired '+role+' clumpiness with normalization removed',xlimits=lim,ylimits=lim,div_y=10))
    LOG.info('Rendering relations to six GALFIT properties')
    relation_figs,trends=book.relations(joined);trends.to_csv(data/'galfit_binned_trends.csv',index=False)
    qcomparison=book.density(unflagged,'gf_q','axis_ratio_positive_aperture_v1','axis_ratio_comparison','Image-aperture and GALFIT model axis ratios',xlimits=(.3,1),ylimits=(0,1))
    sizecomparison=book.density(unflagged,'gf_re_pixels','rpetro','radius_comparison','Petrosian radius and GALFIT effective radius are distinct')
    # Correlation overview uses all finite pairs in the declared cohort, not zoomed bins.
    fig,axes=plt.subplots(1,3,figsize=(17,7),constrained_layout=True)
    labels=['Mean surface\nbrightness','Effective\nradius','Magnitude','Sersic n','Axis ratio','Reduced\nchi-square']
    for ax,band in zip(axes,BANDS):
        sub=correlations[(correlations.cohort=='gf_unflagged')&(correlations.band==band)]
        mat=sub.pivot(index='metric',columns='property',values='spearman').reindex(index=METRICS,columns=PROPERTIES)
        im=ax.imshow(mat,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        ax.set_xticks(range(len(PROPERTIES)),labels,rotation=60,ha='right');ax.set_yticks(range(len(METRICS)),[short_label(k) for k in METRICS]);ax.set_title(f'{band} band')
        for y in range(len(METRICS)):
            for x in range(len(PROPERTIES)):
                value=mat.iloc[y,x];ax.text(x,y,f'{value:.2f}',ha='center',va='center',color='white' if abs(value)>.55 else '#222',fontsize=8)
    fig.colorbar(im,ax=axes,label=r'Spearman rank correlation $\rho_S$',shrink=.7)
    matrix=book.save(fig,'galfit_correlation_matrix','Descriptive morphology–GALFIT rank correlations','Computed on all finite matched pairs in the main unflagged cohort. The complete pair counts and the interior-cohort sensitivity results are downloadable. Many tests share quantities and galaxies; these are not independent discovery claims or causal effects.')
    pd.DataFrame(book.records).to_csv(data/'plot_range_audit.csv',index=False)
    plot_index=[];counter=0
    def figure(item):
        nonlocal counter
        counter+=1;plot_index.append(dict(number=counter,name=item['name'],title=item['title']));return figure_html(item,counter)
    stats=pd.read_csv(root/'statistics.csv');metric_coverage=pd.read_csv(root/'metric_coverage.csv');pair_comparisons=pd.read_csv(root/'paired_comparisons.csv');band_comparisons=pd.read_csv(root/'band_comparisons.csv');summary=json.loads((root/'summary.json').read_text())
    # Dynamic executive findings avoid pinning prose to a previous table revision.
    a_stats=stats[stats.metric=='A_corrected'].set_index('band');negative=100*a_stats.negative/a_stats.valid
    c_pairs=band_comparisons[band_comparisons.metric=='C'].set_index('bands')
    primary_corr=correlations[(correlations.cohort=='gf_unflagged')&correlations.metric.isin(METRICS)]
    mu_rows=primary_corr[(primary_corr.metric=='A_raw')&(primary_corr.property=='gf_mu_mean')].set_index('band')
    cn=primary_corr[(primary_corr.metric=='C')&(primary_corr.property=='gf_n')].set_index('band')
    gm=primary_corr[(primary_corr.metric=='Gini_aperture_v1')&(primary_corr.property=='gf_mu_mean')].set_index('band')
    executive=f'''<div class="summary-grid"><div><strong>84,200</strong><span>upstream galaxy IDs, all retained</span></div><div><strong>252,600</strong><span>object–band rows; 96 missing images</span></div><div><strong>{len(gf):,}</strong><span>unique primary GALFIT exports</span></div></div>
<p>We compare concentration, rotation asymmetry, two clumpiness definitions, aperture Gini/M20 and positive-flux shape in the DES g, r and i images. The sample is the complete supplied candidate catalog, not an independently confirmed or population-representative galaxy sample.</p>
<ul><li><strong>Noise corrections dominate the asymmetry interpretation.</strong> Corrected A is negative in {negative.min():.1f}–{negative.max():.1f}% of its finite per-band measurements. Both raw and corrected values are retained; neither is calibrated here as a physical disturbance classifier.</li>
<li><strong>Band differences are measurable but not uniquely physical.</strong> The matched-object median concentration differences are {c_pairs.loc['g-r','median_difference']:.3f} for g−r and {c_pairs.loc['g-i','median_difference']:.3f} for g−i. Noise, PSF and measurement apertures remain potential explanations alongside stellar populations.</li>
<li><strong>Model-based comparisons now use identified GALFIT data.</strong> In the unflagged matched cohort, raw asymmetry versus mean effective surface brightness has rank correlations g={mu_rows.loc['g','spearman']:.3f}, r={mu_rows.loc['r','spearman']:.3f}, i={mu_rows.loc['i','spearman']:.3f}. These are descriptive associations, not causal tests; surface brightness also depends mathematically on magnitude, size and axis ratio.</li>
<li><strong>Concentration tracks fitted profile shape.</strong> Its correlation with Sersic n is g={cn.loc['g','spearman']:.3f}, r={cn.loc['r','spearman']:.3f}, i={cn.loc['i','spearman']:.3f}. Aperture Gini instead has a strong negative association with mean effective surface brightness (correlations {gm.spearman.min():.3f} to {gm.spearman.max():.3f}), making surface-brightness and noise sensitivity central to its interpretation.</li>
<li><strong>Presentation zooms expose the bulk without deleting tails.</strong> Distribution panels combine bands, compare definitions on matched scales, and offer full-range companions. Every view reports finite samples and off-screen counts.</li></ul>'''
    coverage_display=coverage.rename(columns={'band':'Band','requested':'Requested','matched':'GALFIT matched','unmatched':'No GALFIT row','unflagged':'Main cohort','interior':'Interior cohort','parameter_warning':'Parameter warning','invalid_core':'Invalid core','near_n_boundary':'n near bound','near_q_boundary':'q near bound'})
    sources_table=pd.DataFrame([{'Role':'Primary model comparison','Model':'Single Sersic; per-band magnitude and effective radius; shared n/q/position; ±2.5 arcsec center constraints','Objects':len(gf),'Configuration':PRIMARY},
                                {'Role':'Configuration sensitivity only','Model':'Single Sersic; magnitude per band, radius order 1 / constant-offset configuration','Objects':len(alt),'Configuration':ALTERNATE}])
    eq1=equation(r'$C=5\log_{10}(r_{80}/r_{20}),\quad A_{\rm raw}=\frac{\sum|I-I_{180}|}{\sum|I|},\quad A_{\rm corr}=A_{\rm raw}-A_{\rm sky}$')
    eq2=equation(r'$S=f\,\frac{\sum_{R_{\rm in}<r\leq 1.5R_{\rm P}}\max(I-\widetilde I,0)}{\sum_{R_{\rm in}<r\leq 1.5R_{\rm P}} I}-S_{\rm sky}$')
    eq3=equation(r'$G=\frac{\sum_{i=1}^{N}(2i-N-1)|I|_{(i)}}{(N-1)\sum_i |I_i|},\qquad M_{20}=\log_{10}\!\left(\frac{\sum_{\rm brightest\ 20\%} I_i r_i^2}{\sum_{\rm aperture} I_i r_i^2}\right)$')
    eq4=equation(r'$\langle\mu\rangle_e=m+2.5\log_{10}(2\pi q R_e^2),\qquad R_e\ {\rm in\ arcsec}$')
    sources_table=sources_table.drop(columns='Configuration')
    methods=f'''<p>The CAS implementation is <code>cas_discrete_v3</code>; both S parameter profiles retain their v1 names. Inputs are signed, sky-subtracted 256-pixel images. Contaminants and nonfinite pixels are excluded, while the target and usable sky remain. Catalog RA/Dec and FITS WCS set initial centers; image/mask WCS alignment was checked. The first descending one-pixel-annulus Petrosian crossing at η=0.2 defines R<sub>P</sub>. The circular measurement aperture is 1.5 R<sub>P</sub>, with no GALFIT-radius substitution.</p>
{eq1}
<p>Raw A is minimized in a bounded subpixel Powell search around the initial center (at most ±2 pixels and 0.25 R<sub>P</sub>). Bilinear rotation requires every contributor to be valid; raw and sky A share the same galaxy-pixel denominator. A rare empty-direction optimizer error falls back to the initial center with a flag. The cleanest of four corner tiles must be at least 95% valid and outside the aperture plus two pixels. The default sky tile is 32 pixels. Missing sky corrections remain missing, and corrected values can be negative. This is a documented approximation to CAS practice, not numerical parity with historical codes [<a href="#ref1">1</a>–<a href="#ref3">3</a>].</p>
{eq2}
{table_html(pd.DataFrame([{'Profile':'Lotz/statmorph-like v1','Boxcar width':'0.25 R_P','Inner exclusion':'0.25 R_P','Factor f':1},{'Profile':'Conselice-like v1','Boxcar width':'0.30 R_P','Inner exclusion':'0.075 R_P','Factor f':10}]))}
<p>The boxcar uses mask-normalized smoothing and reflective edges, with actual width int(width fraction × R<sub>P</sub>). This avoids zero-filled contaminant holes creating artificial smoothing halos. The positive residual and signed denominator use the same annulus. These profiles deliberately differ in more than normalization; raw, sky and corrected results for both remain in the table [<a href="#ref1">1</a>–<a href="#ref3">3</a>]. The retained 0.075 R<sub>P</sub> exclusion equals 0.05 of our outer aperture radius. It is a pilot adaptation: Appendix A of <a href="#ref6">Hambleton et al. (2011)</a> states 0.05 R<sub>P</sub> and describes R<sub>P</sub>/6 smoothing. These numerical definitions have been retained for comparison, without claiming exact historical equivalence. The earlier citation wording is corrected in <code>docs/CAS_METHODS.md</code>.</p>
<p>Aperture Gini sorts absolute fluxes, while aperture M20 uses signed flux and signed second moments about the signed centroid. Bright pixels are included through the first 20%-flux crossing, with ties retained. Nonpositive moments or total flux yield missing M20 with flags. Positive-flux shape uses max(I,0) covariance weights: q<sub>+</sub> is the square root of the minor/major eigenvalue ratio, ellipticity is 1−q<sub>+</sub>, and the angle is modulo 180° in image coordinates. These aperture adaptations differ from segmentation-based Gini/M20 [<a href="#ref2">2</a>, <a href="#ref3">3</a>].</p>
{eq3}
<p>The GALFIT relation columns are read from the existing single-Sersic exports [<a href="#ref4">4</a>]; they are not refits, and no two-component nuclear parameters or BIC-derived class labels are used. The g-band columns have no suffix, while r and i append <code>_r</code> and <code>_i</code>. The radius is a semi-major-axis half-model-light radius in arcseconds; the Petrosian radius remains a separate image-based quantity. The stored mean effective surface brightness is consistent to {source_meta['mean_surface_brightness_formula_max_abs_residual']:.2g} mag on the main cohort with <a href="#ref5">Graham &amp; Driver (2005)</a>:</p>
{eq4}
<p>Surface brightness at R<sub>e</sub>, μ<sub>e</sub>, is a different stored column. Larger magnitudes or surface-brightness values mean fainter light. No extra extinction, K correction, cosmological dimming or physical-distance conversion is introduced in this analysis. The producer uses magnitude zeropoint 30; we do not assume additional calibration beyond its stored values.</p>'''
    plotting='''<p>For each quantity, bands share bin edges and axis limits. Zoomed limits use the pooled 0.5th–99.5th percentiles plus 3.5% padding; bounded shape parameters use their natural ranges. Histograms show percentage of the <em>entire finite sample per bin</em>, not a renormalized cropped sample. Full-range views are expandable. Each figure provides its limits and below/above-window counts; the complete machine-readable audit is downloadable.</p>
<p>For clumpiness distributions, Conselice axes are ten times the corresponding Lotz axes, making normalization transparent. Overlay and paired panels explicitly divide Conselice S by ten for display only. Raw and corrected states have separate comparison rows. No stored value, outlier or negative correction is altered.</p>
<p>GALFIT trend plots use 16 common bins in x: linear for magnitude, surface brightness, n and q, logarithmic for radius and stored reduced χ². Each band’s line is the median, with 16th–84th percentile spread; bins with fewer than 30 pairs are left empty. The spread is not an error bar on the median. All finite y values in the x bin contribute before display cropping. Rank correlations use all finite pairs in the declared cohort, without percentile clipping. No significance thresholds, multiplicity-adjusted discoveries or causal conclusions are claimed.</p>'''
    morphology_sections=f'''<h3 id="structure">4.1 Concentration and light distribution</h3><p>These indices describe different aspects of the aperture light distribution. Gini/M20 here share the circular aperture rather than the segmentation used by common reference implementations; no published classification boundaries are applied.</p>{figure(structure)}{figure(joint)}
<h3 id="asymmetry">4.2 Rotation asymmetry and sky correction</h3><p>Band overlays isolate distribution shifts. Paired density panels then show the actual raw/corrected association on the same objects; negative corrected values are visible rather than truncated at zero.</p>{figure(asym)}{figure(asym_overlay)}{figure(paired_a)}
<h3 id="clumpiness">4.3 Both clumpiness profiles</h3><p>Read each profile’s axis factor explicitly. The factor-of-ten display adjustment does not remove differences in the smoothing kernel and central exclusion.</p>{figure(clump)}{figure(clump_overlay)}{''.join(figure(item) for item in paired_s)}
<h3 id="shape">4.4 Aperture shape</h3><p>Positive-flux weighting can bias faint-source shapes toward the background distribution. Image position angles are not celestial angles, and their interpretation is unstable for round sources. We therefore do not compute ordinary linear position-angle correlations with GALFIT angles.</p>{figure(shape)}'''
    main_cohort=correlations[correlations.cohort=='gf_unflagged'].set_index(['band','metric','property'])
    interior=correlations[correlations.cohort=='gf_interior'].set_index(['band','metric','property'])
    delta=(main_cohort.spearman-interior.spearman).abs().dropna()
    cohort_note=f'''<p>The main relation cohort requires a returned successful fit, finite core quantities, positive R<sub>e</sub>, n, q and reduced χ², q≤1, and no x/y/magnitude/R<sub>e</sub>/n/q parameter warning. This is called <em>unflagged</em>, not “validated”: extreme or biased fits can still pass. No arbitrary magnitude, size, surface-brightness or χ² threshold is imposed. Near-boundary cases are counted, not silently removed.</p>
<p>A separate sensitivity cohort removes n within 0.001 of 0.2 or 5 and q within 0.001 of 0.3 or 1, the configured bounds. For the reported property/metric combinations, the median absolute change in Spearman correlation is {delta.median():.3f}; the largest is {delta.max():.3f}. These cohorts differ in selection, so stability is diagnostic rather than proof that either is physically unbiased. All per-band pair counts and correlations for both cohorts are available in the CSV.</p>'''
    relation_findings=f'''<h3>5.1 Empirical associations</h3><p><strong>Surface brightness:</strong> raw A increases strongly toward fainter model mean surface brightness, with g/r/i rank correlations {mu_rows.loc['g','spearman']:.3f}/{mu_rows.loc['r','spearman']:.3f}/{mu_rows.loc['i','spearman']:.3f}. Corrected A slopes downward in the binned trends. This is consistent with substantial signal-to-noise and sky-correction sensitivity; it does not identify an intrinsic asymmetry–surface-brightness law.</p>
<p><strong>Profile shape and size:</strong> C correlates positively with fitted Sersic n (g/r/i: {cn.loc['g','spearman']:.3f}/{cn.loc['r','spearman']:.3f}/{cn.loc['i','spearman']:.3f}). Its effective-radius associations are weaker and band dependent; the size panels and the separate R_P-versus-R_e figure show why these radii cannot be interchanged. The fitted profile and the non-parametric light distribution provide related, but different, descriptions of the same image.</p>
<p><strong>Aperture Gini:</strong> its strong negative surface-brightness correlations ({gm.spearman.min():.3f} to {gm.spearman.max():.3f}) are much larger in magnitude than its fitted-n correlations. A noise contribution to the absolute-flux distribution is one possible explanation for the flattening at low surface brightness; this interpretation requires injection/recovery testing. The pattern motivates targeted noise/aperture experiments, not a surface-brightness-independent physical classification.</p>'''
    relation_html=figure(matrix)+''.join(figure(item) for item in relation_figs)+figure(qcomparison)+figure(sizecomparison)
    # Human-readable correlation table for the requested principal photometric properties.
    corrtab=correlations[(correlations.cohort=='gf_unflagged')&correlations.property.isin(PROPERTIES[:3])&correlations.metric.isin(['C','A_raw','A_corrected','Gini_aperture_v1','M20_aperture_v1'])].copy()
    corrtab['metric']=corrtab.metric.map(short_label);corrtab['property']=corrtab.property.map(short_label)
    corrtab=corrtab.drop(columns='cohort').rename(columns={'band':'Band','metric':'Morphology','property':'GALFIT property','n':'Finite pairs','spearman':'Spearman correlation'})
    flags=pd.DataFrame([{'Flag':k,'Rows':v} for k,v in summary['quality_flags'].items()])
    flag_explanations='''<ul><li><code>no_clean_sky_tile</code>: no corner satisfies the valid-pixel and target-separation tests; corrected A/S are missing.</li><li><code>nonmonotone_growth_curve</code>: signed cumulative flux has a sufficiently large downward step; C is withheld.</li><li><code>insufficient_scale_or_flux</code>: S kernel below two pixels, too few annular pixels or nonpositive annular flux.</li><li><code>asymmetry_no_flux</code>: fewer than ten symmetrically valid pixels or nonpositive absolute-flux denominator.</li><li><code>center_optimization_failed</code>: bounded optimization did not succeed or hit the rare empty-direction error; the initial-center fallback is flagged.</li><li><code>aperture_truncated</code>: the full aperture and center-search margin do not fit inside the stamp.</li><li><code>petrosian_crossing_not_found</code>: no valid descending eta=0.2 crossing was found; no fitted-radius fallback is used.</li><li><code>center_search_boundary</code>: the optimized center lies near a search bound.</li><li><code>nonpositive_or_empty_total_flux</code>: insufficient pixels or nonpositive signed aperture total for concentration.</li></ul>'''
    provenance_meta=dict(report_version='scientific_report_v2',built_utc=datetime.now(timezone.utc).isoformat(),morphology=provenance,
                         galfit=source_meta,display_policy={'quantiles':[.005,.995],'padding':.035,'histogram_bins':60,'outliers_removed_from_measurements':False,'trend_bins':16,'trend_minimum_pairs':30},
                         figures=plot_index,source_code_sha256={str(p):sha256_file(p) for p in [Path(__file__),Path('algorithmic_morphology/report_data.py'),Path('algorithmic_morphology/paper_plots.py')]})
    # Detailed information is retained in expandable appendices rather than interrupting the paper.
    appendix=f'''<details><summary>A. All column definitions and units ({len(dictionary_table)} entries)</summary><p>Exact machine column names are searchable here and in the downloadable dictionary. S raw/sky/corrected values have separate entries. GALFIT source names use the documented band suffix.</p><input class="search" data-search="dictionary" aria-label="Filter column definitions" placeholder="Filter column names or meanings…"><div id="dictionary">{table_html(dictionary_table)}</div></details>
<details><summary>B. Every morphology metric’s valid sample and summary statistics</summary>{table_html(metric_coverage)}{table_html(stats)}</details>
<details><summary>C. Quality flags and failure meanings</summary>{table_html(flags)}{flag_explanations}<p>All 96 error rows are missing cutouts; 4,704 rows have invalid C or raw A. An invalid CAS status does not imply that every extended diagnostic is missing. Flag totals overlap.</p></details>
<details><summary>D. Paired cross-band and clumpiness comparisons</summary>{table_html(pair_comparisons)}{table_html(band_comparisons)}</details>
<details><summary>E. Exact data paths, hashes and fit/configuration provenance</summary><p>The morphology file and its checksum are unchanged. The report uses the compact export of the same canonical table; the full table retains input-image fingerprints. GALFIT values are read from exported tables, not recomputed fits. The previous HTML revision remains available as <a href="report_v1.html">report v1</a>.</p><pre>{escape(json.dumps(provenance_meta,indent=2,default=str))}</pre></details>
<details><summary>F. Reproduction and downloadable analysis tables</summary><pre>python scripts/build_paper_report.py --report-dir reports/full_catalog_20261020 --source-project ../nucleated_lsbg</pre><p>This revision builder refuses to replace an existing report_v2.html. For a new revision, preserve the current report and choose a new version. Measurements are never recalculated by this command.</p><ul><li><a href="measurements_metrics.parquet">Immutable morphology measurements</a></li><li><a href="paper_v2_data/galfit_join.parquet">GALFIT-only join keyed by object and band</a></li><li><a href="paper_v2_data/galfit_sources.json">Source files, hashes and suffix mapping</a></li><li><a href="paper_v2_data/galfit_coverage.csv">Join and quality-cohort coverage</a></li><li><a href="paper_v2_data/galfit_correlations.csv">All finite-pair correlations and interior sensitivity</a></li><li><a href="paper_v2_data/galfit_binned_trends.csv">Binned medians and percentile spreads</a></li><li><a href="paper_v2_data/plot_range_audit.csv">Every plot’s display limits and tail counts</a></li><li><a href="paper_v2_data/column_dictionary.csv">Complete column dictionary</a></li><li><a href="paper_v2_data/fit_configuration_sensitivity.csv">Primary versus radius-order-1 / constant-offset fit comparison</a></li></ul></details>'''
    refs='''<ol><li id="ref1">Conselice (2003), <em>The Relationship between Stellar Light Distributions of Galaxies and Their Formation Histories</em>, ApJS 147, 1. <a href="https://arxiv.org/abs/astro-ph/0303065">Paper</a>. CAS motivation; our exact conventions and named adaptations are specified above.</li>
<li id="ref2">Lotz, Primack &amp; Madau (2004), <em>A New Non-Parametric Approach to Galaxy Morphological Classification</em>, AJ 128, 163. <a href="https://arxiv.org/abs/astro-ph/0311352">Paper</a>. Gini/M20 and morphology robustness; our aperture selection differs.</li>
<li id="ref3">Rodriguez-Gomez et al. (2019), <em>The optical morphologies of galaxies in the IllustrisTNG simulation</em>, MNRAS 483, 4140. <a href="https://doi.org/10.1093/mnras/sty3341">Paper</a>; <a href="https://statmorph.readthedocs.io/en/latest/api.html">statmorph API</a>; <a href="https://statmorph.readthedocs.io/en/latest/_modules/statmorph/statmorph.html">reference source</a>. Reference conventions, not an equivalence claim.</li>
<li id="ref4">Peng et al. (2010), <em>Detailed Decomposition of Galaxy Images. II. Beyond Axisymmetric Models</em>, AJ 139, 2097. <a href="https://arxiv.org/abs/0912.0731">Paper</a>. GALFIT structural modeling; this analysis reads the sibling project's stored multiband fits.</li>
<li id="ref5">Graham &amp; Driver (2005), <em>A concise reference to projected Sersic quantities</em>, PASA 22, 118. <a href="https://ned.ipac.caltech.edu/level5/March05/Graham/Graham_contents.html">Paper</a>. Effective and mean effective surface-brightness conventions.</li>
<li id="ref6">Hambleton et al. (2011), <em>Advanced morphological galaxy classification: a comparison of observed and simulated galaxies</em>, MNRAS 418, 801. <a href="https://arxiv.org/abs/1107.6045">Paper</a>; <a href="https://arxiv.org/html/1107.6045v2#A1">Appendix A</a>. Central-exclusion and smoothing conventions differ from our retained pilot profile.</li>
<li id="ref7">Thorp et al. (2021), <em>Towards robust determination of non-parametric morphologies in marginal astronomical data</em>, MNRAS 507, 886. <a href="https://arxiv.org/abs/2107.13578">Paper</a>. Noise and resolution robustness context.</li></ol>'''
    figure_navigation='<details id="figure-index"><summary>Jump to a specific figure</summary><ol>'+''.join(f'<li><a href="#fig-{v["name"]}">Figure {v["number"]}. {escape(v["title"])}</a></li>' for v in plot_index)+'</ol></details>'
    nav=[('summary','Executive summary'),('data','1 · Data and matching'),('methods','2 · Methods'),('plotting','3 · Reading the figures'),('distributions','4 · Distributions'),('galfit','5 · GALFIT relations'),('discussion','6 · Discussion'),('appendix','Appendices / columns'),('references','References')]
    css='''html{scroll-behavior:smooth;scroll-padding-top:4rem}body{font:16px/1.65 system-ui,sans-serif;color:var(--foreground,#25303a);background:var(--background,#fff);margin:0}header,main{max-width:1180px;margin:auto;padding:1.5rem 2rem}header{padding-bottom:.5rem}h1{font-family:Georgia,serif;font-size:2.25rem;line-height:1.2;max-width:950px}h2{font-family:Georgia,serif;font-size:1.8rem;border-bottom:1px solid var(--border,#d6dce0);padding-bottom:.4rem;margin-top:3rem}h3{margin-top:2.3rem}h4{font-size:1.05rem}p,li{max-width:100ch}nav{position:sticky;top:0;z-index:5;background:var(--background,#fff);border-block:1px solid var(--border,#d6dce0);padding:.6rem 1rem;display:flex;gap:1.2rem;overflow-x:auto;white-space:nowrap;font-size:.85rem}a{color:inherit;text-decoration:underline;text-decoration-color:#4a96ba;text-underline-offset:3px}.subtitle,.meta{color:var(--muted-foreground,#566675)}.summary-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:1rem;margin:1rem 0}.summary-grid div{padding:1rem;background:var(--card,#f1f5f7);border:1px solid var(--border,#d6dce0);border-radius:8px}.summary-grid strong{display:block;font-size:1.7rem}.summary-grid span{font-size:.9rem}figure{margin:2rem 0 2.6rem}figure>img{display:block;width:100%;height:auto;border:1px solid var(--border,#ddd);border-radius:4px;background:#fff}figcaption{font-size:.9rem;color:var(--muted-foreground,#52616c);margin-top:.6rem;max-width:110ch}.equation{max-width:100%;height:auto;background:#fff;border-radius:5px;margin:.7rem 0}.table-scroll{overflow:auto;margin:1rem 0}table{border-collapse:collapse;font-size:.87rem;width:100%}th,td{text-align:left;vertical-align:top;padding:.55rem .7rem;border-bottom:1px solid var(--border,#d6dce0)}th{background:var(--card,#f1f5f7)}td code,code{overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;background:var(--code-background,#f1f5f7);color:var(--code-foreground,#25303a);padding:1rem;font-size:.82rem;border-radius:6px}details{border:1px solid var(--border,#d6dce0);padding:.8rem 1rem;border-radius:5px;margin:1rem 0}summary{cursor:pointer;font-weight:600}details img{max-width:100%}.search{font:inherit;padding:.5rem;width:min(90%,650px);margin:.7rem 0}.section-links{display:flex;gap:1rem;flex-wrap:wrap}.badge{font-size:.85rem;letter-spacing:.05em;text-transform:uppercase}.note{border-left:3px solid #d7963b;padding-left:1rem}@media(max-width:700px){header,main{padding:1rem}.summary-grid{grid-template-columns:1fr}h1{font-size:1.7rem}h2{font-size:1.4rem}nav{gap:1rem}}@media print{nav{position:static}details{break-inside:avoid}figure{break-inside:avoid}body{font-size:11pt}header,main{padding:0;max-width:none}}'''
    html=f'''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DES galaxy morphology: distributions and structural context</title><style>{css}</style></head><body><header><p class="badge">DES candidate catalog · scientific analysis · revision 2</p><h1>Galaxy morphology across bands and its relation to fitted structure</h1><p class="subtitle">Concentration, asymmetry, clumpiness, aperture light distribution and GALFIT properties for the complete 84,200-object input catalog.</p><p class="meta">Measured 20 October 2026 · g, r and i bands · Detailed methods and reproducibility below</p></header><nav aria-label="Report sections">{''.join(f'<a href="#{key}">{text}</a>' for key,text in nav)}</nav><main>
<section id="summary"><h2>Executive summary</h2>{executive}</section>
<section id="data"><h2>1. Data, models and object matching</h2><p>The morphology table contains the same 84,200 upstream catalog IDs once in each of g, r and i: 252,600 requested rows. It uses existing two-pass sky-subtracted 256-pixel images and contaminant masks in <code>nucleated_lsbg</code>. Ninety-six input image combinations are missing; all errors and invalid individual measurements remain in the data. The earlier 75-object wide-image pilot is not used for these full-catalog distributions.</p>
<p>For structural context, the primary data are the <strong>single-Sersic, per-band-radius GALFITM export</strong> below. We choose it explicitly because it covers 84,061 objects, whereas the radius-order-1 / constant-offset export contains 5,406. Both use the same parent catalog and shared 256-pixel image/mask product tree. The primary configuration fits magnitude and radius separately by band while sharing n, q and position. These are stored model fits, with different measurement assumptions from our non-parametric circular apertures.</p>{table_html(sources_table)}
<p>Joins are exact on stable object ID and same-band label, validated one-to-one. Unmatched objects remain in the morphology sample and disappear only from the declared paired-fit analysis. No fit properties are taken from the source-catalog proxy columns in place of actual GALFIT results.</p>{table_html(coverage_display)}<details><summary>Exact GALFIT file locations and fingerprints</summary><pre>{escape(json.dumps(source_meta,indent=2))}</pre></details></section>
<section id="methods"><h2>2. Measurement and model conventions</h2>{methods}</section>
<section id="plotting"><h2>3. Reading the figures and analysis cohorts</h2>{plotting}{cohort_note}{figure_navigation}</section>
<section id="distributions"><h2>4. Distributions and direct comparisons</h2><div class="section-links"><a href="#structure">Structure</a><a href="#asymmetry">Asymmetry</a><a href="#clumpiness">Clumpiness</a><a href="#shape">Shape</a></div>{morphology_sections}</section>
<section id="galfit"><h2>5. Morphology in GALFIT structural context</h2><p>These panels compare non-parametric image measurements with fitted mean surface brightness, effective radius, total apparent magnitude, Sersic index, model axis ratio and stored reduced χ². Relationships are computed in the main unflagged cohort; paired counts vary by metric and band. GALFIT n and q are shared across bands in the source model, so the three band panels are not independent object samples.</p>{relation_findings}{relation_html}
<details><summary>Numerical correlations for surface brightness, radius and magnitude</summary>{table_html(corrtab)}</details>
<h3>5.2 Sensitivity to fit configuration</h3><p>The smaller radius-order-1 / constant-offset export is used only on common object IDs with unflagged fits in both configurations. Differences below are per-band primary minus radius-order-1 / constant-offset values. The alternate config sets radius optimization order 1 (rather than 1,1,1); the stored band radii need not be numerically identical. This is a limited common-cohort sensitivity check; its 5,406-row source is not substituted into the broader primary comparison.</p>{table_html(sensitivity_table)}</section>
<section id="discussion"><h2>6. Interpretation, limitations and next tests</h2><p>The strongest descriptive patterns combine intrinsic structure, surface brightness, resolution and measurement noise. Apparent magnitude, effective radius, axis ratio and mean effective surface brightness are mathematically coupled in a Sersic model; correlations with them must not be presented as independent physical drivers. Band-specific noise/PSF differences, conditional GALFIT availability and fit constraints can contribute to the observed associations.</p>
<p>Large raw asymmetry and a predominance of negative corrected asymmetry are consistent with substantial correction sensitivity. The 144 synthetic pilot cases showed noise-dependent raw A and negative corrected-A bias for symmetric scenes; they did not calibrate DES-specific uncertainties. <a href="#ref7">Thorp et al. (2021)</a> provide related noise/resolution robustness context; their corrections are not applied to these DES values. The chosen minimum-scatter corner and the integer sky rotation versus subpixel galaxy interpolation remain potential bias sources. Blank-region alternatives, matched reference-code settings, injected galaxies with realistic sky and PSF, and mask perturbations are required before physical morphology classification.</p>
<p>Clumpiness profiles retain their different definitions. Strong rank agreement does not show that either profile is preferable, and normalization-adjusted overlays do not make their excluded regions or smoothing scales equivalent. Aperture Gini/M20 and positive-flux shapes are further affected by aperture choice and signed background noise. No merger labels, Hubble types or published segmentation-based thresholds are assigned here.</p>
<p>All summary statistics are descriptive. The multiple associations share the same galaxies and derived quantities; this report does not use p-values to claim new independent effects. The interior-fit comparison and the radius-order-1 / constant-offset comparison quantify selection/model sensitivity on their stated cohorts rather than resolving all systematic uncertainty.</p></section>
<section id="appendix"><h2>Appendices: definitions, coverage and reproducibility</h2>{appendix}</section>
<section id="references"><h2>References</h2>{refs}</section><p><a href="#summary">Back to executive summary ↑</a></p></main><script>document.querySelectorAll('[data-search]').forEach(input=>input.addEventListener('input',()=>{{const q=input.value.toLowerCase();document.querySelectorAll('#'+input.dataset.search+' tbody tr').forEach(row=>row.hidden=!row.textContent.toLowerCase().includes(q));}}));</script></body></html>'''
    htmlpath=root/'report_v2.html';htmlpath.write_text(html)
    pd.DataFrame(plot_index).to_csv(data/'figure_index.csv',index=False)
    provenance_meta['html_sha256']=sha256_file(htmlpath);provenance_meta['compact_table_unchanged']=sha256_file(morph_path)==provenance['compact_table_sha256']
    provenance_meta['files_sha256']={str(p.relative_to(root)):sha256_file(p) for p in sorted(data.iterdir()) if p.is_file()}
    (root/'report_v2_manifest.json').write_text(json.dumps(provenance_meta,indent=2,default=str))
    if not (root/'report_v1.html').exists():shutil.copy2(root/'index.html',root/'report_v1.html')
    temporary=root/'.index_v2.tmp';temporary.write_text(html);os.replace(temporary,root/'index.html')
    LOG.info('Published scientific report revision 2: %s; %d main figures; %d range-audit rows',htmlpath,len(plot_index),len(book.records))


if __name__=='__main__':main()
