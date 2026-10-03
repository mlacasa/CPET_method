"""E20 secondary summaries of immutable E08 metrics; never reruns simulations.

Run with the project Anaconda interpreter. Creates a fresh timestamped run.
No manuscript files are edited by this script.
"""
from pathlib import Path
from datetime import datetime
import json
import hashlib
import shutil
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
J = HERE.parent
SOURCE = J / 'punto5/runs/smoothing_20260930T104226+0200'
METHODS = ['raw', 'mean5_p1', 'median5_p1', 'mean5_p2', 'median5_p2']
NAMES = dict(zip(METHODS, ['No pre-smoothing', 'Mean, one pass', 'Median, one pass', 'Mean, two passes', 'Median, two passes']))
FAMILIES = ['clean', 'gaussian', 'impulse', 'burst', 'saw', 'true_peak', 'true_oscillation', 'border', 'slope', 'missing']
FLABELS = ['Clean', 'Gaussian', 'Impulses', 'Bursts', 'Saw noise', 'True peaks', 'True oscillations', 'Edge peak', 'Slope change', 'Missingness']
SEED = 20261001
B = 2000

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def dump(p, obj):
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf-8')

def main():
    run = HERE / 'runs' / datetime.now().astimezone().strftime('validation_%Y%m%dT%H%M%S%z')
    run.mkdir(parents=True, exist_ok=False)
    for sub in ['data', 'tables', 'figures', 'source_snapshot']:
        (run/sub).mkdir()
    config = dict(stage='E20', source=str(SOURCE), created_at=datetime.now().astimezone().isoformat(),
        purpose='Secondary engineering validation summaries using E08 known-truth metrics',
        primary_outputs=['final-curve RMSE', 'paired-difference RMSE', 'signed bias', 'coverage'],
        methods=METHODS, adopted='mean5_p1', main_channel='VCO2', supplementary_channels='all 17',
        family_order=FAMILIES, aggregation='replicates already averaged in source; equal sessions and scenarios within family; equal family means within participant; equal participants',
        uncertainty=dict(seed=SEED, resamples=B, generator='NumPy PCG64', unit='paired participant template',
            same_indices_for_all_endpoints=True, interval='nominal 95% percentile, numpy linear quantiles',
            interpretation='conditional on fixed synthetic design and stored replicate averages; not clinical or Monte Carlo uncertainty; not multiplicity adjusted'),
        missing='E08 finite-value aggregation retained; report per-endpoint n and unavailable source counts. Bootstrap all 60 clusters and average available endpoints; retain zero coverage. Missing family means must occur for all families of a participant endpoint, not silently change family weights.',
        normalized_errors='in E08 channel amplitude scale, not relative to true paired change or percentages',
        contrasts='paired participant mean5_p1 minus raw; no p-values, significance tally or clinical equivalence claim',
        status='revision-stage analysis; E08 aggregate results already known; configuration saved before new summaries/intervals',
        new_simulations=False, change_filter_selection=False)
    dump(run/'analysis_config.json', config)
    shutil.copy2(__file__, run/'source_snapshot/summarize_validation.py')
    inputs = sorted((SOURCE/'synthetic').glob('metrics_template_*.parquet')) + sorted((SOURCE/'synthetic').glob('paired_template_*.parquet'))
    for name in ['analysis_config.json', 'channel_scales.csv', 'truth_model_dictionary.json', 'data_dictionary.json',
                 'synthetic/scenario_manifest.csv', 'tables/T-SR4_known_paired_change.csv', 'tables/T-SR4_by_phase_variable.csv', 'tables/T-SR4_global.csv']:
        inputs.append(SOURCE/name)
    frozen = [dict(path=str(p), sha256=sha(p), bytes=p.stat().st_size) for p in inputs]
    dump(run/'input_manifest.json', frozen)
    print('RUN', run, flush=True)
    curve = pd.concat([pd.read_parquet(p) for p in inputs if p.name.startswith('metrics_template_')], ignore_index=True)
    pair = pd.concat([pd.read_parquet(p) for p in inputs if p.name.startswith('paired_template_')], ignore_index=True)
    assert len(curve)==60*2*2*17*31*5 and len(pair)==60*2*17*31*5
    keys=['template_id','phase','channel','scenario','method']
    assert not curve.duplicated(keys+['session']).any() and not pair.duplicated(keys).any()
    assert set(curve.family)==set(FAMILIES) and set(pair.method)==set(METHODS)
    assert (curve.replicates==5).all() and (pair.replicates==5).all()
    scales=pd.read_csv(SOURCE/'channel_scales.csv').set_index('channel').amplitude
    channels=list(scales.index)
    for f in [curve,pair]:
        assert f.template_id.nunique()==60
    ccols=['grid_rmse','grid_rmse_units','grid_mae','grid_bias','grid_coverage']
    pcols=['paired_grid_rmse','paired_grid_mae','paired_grid_bias','paired_grid_coverage']
    pair['paired_grid_rmse_units']=pair.paired_grid_rmse*pair.channel.map(scales)
    pair['paired_grid_bias_units']=pair.paired_grid_bias*pair.channel.map(scales)
    curve['grid_bias_units']=curve.grid_bias*curve.channel.map(scales)
    ccols+=['grid_bias_units']; pcols+=['paired_grid_rmse_units','paired_grid_bias_units']
    audits=[]
    for name,f,cols in [('curve',curve,ccols),('paired',pair,pcols)]:
        for k,g in f.groupby(['phase','channel','family','method'],sort=True):
            for col in cols:
                audits.append(dict(kind=name,phase=k[0],channel=k[1],family=k[2],method=k[3],metric=col,
                    rows=len(g),finite=int(np.isfinite(g[col]).sum()),missing=int(g[col].isna().sum()),zero_coverage=int((g['grid_coverage' if name=='curve' else 'paired_grid_coverage']==0).sum())))
    pd.DataFrame(audits).to_csv(run/'data/source_availability.csv',index=False)
    # All methods restore the same original masks and share binning/interpolation support.
    coverage_spreads=[]
    for f,col,extra in [(curve,'grid_coverage',['session']),(pair,'paired_grid_coverage',[])]:
        grouped=f.groupby(keys[:-1]+extra)[col]
        spread=grouped.max()-grouped.min(); assert np.allclose(spread,0,atol=1e-14)
        coverage_spreads.append(float(spread.max()))
    fkeys=['template_id','phase','channel','family','method']
    cf=curve.groupby(fkeys)[ccols].mean().reset_index()
    pf=pair.groupby(fkeys)[pcols].mean().reset_index()
    fam=cf.merge(pf,on=fkeys,validate='one_to_one')
    cols=ccols+pcols
    assert fam.groupby(['template_id','phase','channel','method']).size().eq(10).all()
    family_counts=fam.groupby(['template_id','phase','channel','method'])[cols].count()
    assert family_counts.isin([0,10]).all().all(), 'Partial family availability changes weights: inspect before proceeding'
    participant=fam.groupby(['template_id','phase','channel','method'])[cols].mean().reset_index()
    fam.to_parquet(run/'data/participant_family_metrics.parquet',index=False)
    participant.to_csv(run/'data/participant_metrics.csv',index=False)
    # Reproduce prior aggregates before interpreting new confidence intervals.
    previous=pd.read_csv(SOURCE/'tables/T-SR4_known_paired_change.csv')
    reproduced=participant.groupby(['channel','phase','method'])[['paired_grid_rmse','paired_grid_coverage']].mean().reset_index()
    cmp=previous.merge(reproduced,on=['channel','phase','method'],suffixes=('_old','_new'),validate='one_to_one')
    regression={col:float((cmp[col+'_old']-cmp[col+'_new']).abs().max()) for col in ['paired_grid_rmse','paired_grid_coverage']}
    assert max(regression.values())<1e-12
    previous_curve=pd.read_csv(SOURCE/'tables/T-SR4_by_phase_variable.csv')
    reproduced_curve=participant.groupby(['channel','phase','method'])[['grid_rmse','grid_rmse_units']].mean().reset_index()
    cmpc=previous_curve.merge(reproduced_curve,on=['channel','phase','method'],suffixes=('_old','_new'),validate='one_to_one')
    for col in ['grid_rmse','grid_rmse_units']:
        regression[col]=float((cmpc[col+'_old']-cmpc[col+'_new']).abs().max()); assert regression[col]<1e-10
    ids=sorted(participant.template_id.unique())
    ix=np.random.default_rng(SEED).integers(0,60,size=(B,60))
    np.save(run/'data/bootstrap_indices.npy',ix,allow_pickle=False)
    rows=[]
    def summarize(frame, grouping, level):
        for k,g in frame.groupby(grouping,sort=True):
            k=k if isinstance(k,tuple) else (k,)
            g=g.set_index('template_id').reindex(ids)
            x=g[cols].to_numpy(); assert x.shape==(60,len(cols))
            assert np.isfinite(x[ix]).sum(axis=1).min()>0
            boots=np.nanmean(x[ix],axis=1); lo,hi=np.quantile(boots,[.025,.975],axis=0)
            row=dict(zip(grouping,k));row.update(level=level,n_templates=60)
            for z,col in enumerate(cols): row.update({col:float(np.nanmean(x[:,z])),col+'_low':float(lo[z]),col+'_high':float(hi[z]),col+'_n':int(np.isfinite(x[:,z]).sum())})
            rows.append(row)
    summarize(participant,['phase','channel','method'],'balanced')
    summarize(fam,['phase','channel','family','method'],'family')
    summary=pd.DataFrame(rows)
    balanced=summary[summary.level=='balanced'].drop(columns='family')
    family=summary[summary.level=='family']
    balanced.to_csv(run/'data/validation_by_channel.csv',index=False,float_format='%.15g')
    family.to_csv(run/'data/validation_by_family.csv',index=False,float_format='%.15g')
    contrasts=[]
    for level,f,grouping in [('balanced',participant,['phase','channel']),('family',fam,['phase','channel','family'])]:
        for k,g in f.groupby(grouping,sort=True):
            a=g[g.method=='mean5_p1'].set_index('template_id').reindex(ids)
            b=g[g.method=='raw'].set_index('template_id').reindex(ids)
            assert np.array_equal(np.isfinite(a[cols]),np.isfinite(b[cols]))
            delta=a[cols].to_numpy()-b[cols].to_numpy(); boots=np.nanmean(delta[ix],axis=1)
            lo,hi=np.quantile(boots,[.025,.975],axis=0)
            for z,col in enumerate(cols):
                contrasts.append(dict(zip(grouping,k))|dict(level=level,metric=col,contrast='mean5_p1 minus raw',n_templates=60,n_finite=int(np.isfinite(delta[:,z]).sum()),mean=float(np.nanmean(delta[:,z])),low=float(lo[z]),high=float(hi[z])))
    pd.DataFrame(contrasts).to_csv(run/'data/validation_contrasts.csv',index=False,float_format='%.15g')
    # Main table: show all candidates, not just the adopted method and raw control.
    def ci(row,key): return f"{row[key]:.4f} [{row[key+'_low']:.4f}, {row[key+'_high']:.4f}]"
    tex=[r'\begin{table}[tbp]',r'\centering\scriptsize\color{revisionblue}',
        r'\caption{Reconstruction of synthetic VCO$_2$ curves and known paired differences. Values are means of participant-level, equally family-weighted RMSEs with nominal 95\% bootstrap intervals. Errors are normalized by the channel amplitude scale, not by the true paired change. The control omits pre-smoothing but retains binning and interpolation.}',
        r'\label{tab:known_truth_validation}',r'\begin{tabular}{llll}',r'\toprule',
        r'Method & Curve RMSE [95\% CI] & Difference RMSE [95\% CI] & Joint support (\%) \\',r'\midrule']
    for phase in ['ramp','recovery']:
        tex.append(r'\multicolumn{4}{l}{\textbf{'+phase.title()+r'; 60 paired templates}} \\')
        for method in METHODS:
            row=balanced[(balanced.channel=='VCO2')&(balanced.phase==phase)&(balanced.method==method)].iloc[0]
            tex.append(f"{NAMES[method]} & {ci(row,'grid_rmse')} & {ci(row,'paired_grid_rmse')} & {100*row.paired_grid_coverage:.2f}"+r' \\')
        tex.append(r'\midrule')
    tex+= [r'\bottomrule',r'\end{tabular}',r'\par\smallskip\parbox{\textwidth}{\scriptsize Mean, one pass is the adopted method. Intervals resample the 60 paired participant templates; they are descriptive and conditional on the fixed simulation design. Equal method coverage does not mean complete coverage.}',r'\end{table}']
    (run/'tables/known_truth_validation.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')
    # Main figure: all family outcomes of adopted versus raw, plus shared coverage.
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,2,figsize=(8.2,7.2),sharex='col',layout='constrained')
    for col,phase in enumerate(['ramp','recovery']):
        ax=axs[0,col]
        for method,offset,color in [('raw',-.11,'#777777'),('mean5_p1',.11,'#0046a0')]:
            d=family[(family.channel=='VCO2')&(family.phase==phase)&(family.method==method)].set_index('family').reindex(FAMILIES)
            v=d.paired_grid_rmse.to_numpy();lo=d.paired_grid_rmse_low.to_numpy();hi=d.paired_grid_rmse_high.to_numpy()
            ax.errorbar(np.arange(10)+offset,v,yerr=np.maximum(np.vstack([v-lo,hi-v]),0),fmt='o',ms=4,capsize=2,color=color,label=NAMES[method])
        ax.set_title(('A' if col==0 else 'B')+'  '+phase.title()+'\nPaired-difference error',loc='left',fontweight='bold',fontsize=11)
        ax.set_ylabel('Mean normalized RMSE (95% CI)',fontsize=10);ax.set_ylim(bottom=0);ax.grid(axis='y',alpha=.2);ax.legend(fontsize=9)
        ax=axs[1,col]
        d=family[(family.channel=='VCO2')&(family.phase==phase)&(family.method=='mean5_p1')].set_index('family').reindex(FAMILIES)
        for key,label,style in [('grid_coverage','Single curve','s--'),('paired_grid_coverage','Joint pair','o-')]:
            ax.plot(np.arange(10),100*d[key],style,label=label,ms=4)
        ax.set_ylim(0,105);ax.set_ylabel('Available grid locations (%)');ax.grid(axis='y',alpha=.2)
        ax.set_title(('C' if col==0 else 'D')+'  '+phase.title()+': support',loc='left',fontweight='bold',fontsize=11)
        ax.set_xticks(np.arange(10),FLABELS,rotation=55,ha='right',fontsize=9);ax.legend(fontsize=9,loc='lower left')
        ax.text(.02,.30,'Same support across\nall five methods',transform=ax.transAxes,fontsize=9)
    fig.savefig(run/'figures/known_truth_validation.png',dpi=300)
    plt.close(fig)
    main_rows=balanced[balanced.channel=='VCO2'].copy()
    main_rows.to_csv(run/'data/main_table_values.csv',index=False,float_format='%.15g')
    (family[family.channel=='VCO2']).to_csv(run/'data/main_figure_values.csv',index=False,float_format='%.15g')
    summary_note=dict(regression_max_absolute_differences=regression,coverage_max_method_spreads=coverage_spreads,
        source_rows=dict(curve=len(curve),paired=len(pair)),n_family_rows=len(fam),n_participant_rows=len(participant),
        source_missing_counts={col:int(curve[col].isna().sum()) for col in ccols}|{col:int(pair[col].isna().sum()) for col in pcols},
        point_estimates_only_not_recomputed_simulations=True,bootstrap_resamples=B,seed=SEED,
        conditional_uncertainty=True,inputs_unchanged=all(sha(Path(x['path']))==x['sha256'] for x in frozen))
    assert summary_note['inputs_unchanged']
    dump(run/'quality_assurance.json',summary_note)
    dump(run/'environment.json',dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,matplotlib=matplotlib.__version__))
    (HERE/'LATEST_VALIDATION_RUN.txt').write_text(run.name+'\n',encoding='utf-8')
    print(main_rows[['phase','method','grid_rmse','paired_grid_rmse','paired_grid_rmse_low','paired_grid_rmse_high','paired_grid_coverage']].to_string(index=False),flush=True)
    print('QA',summary_note,flush=True)

if __name__=='__main__':
    main()
