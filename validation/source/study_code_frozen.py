"""F5-03: crossed mean/median x one/two passes; synthetic fidelity and real impact.

No source manuscript, notebook, or historical matrix is modified. All calculations
use a single implementation with original NaNs restored after EACH pass.
"""
from __future__ import annotations
from datetime import datetime
from pathlib import Path
import argparse
import json
import sys
import warnings
import numpy as np
import pandas as pd
from audit_reproduction import HERE,ROOT,M,D,S,UNITS,sha,dump

SEED=20260930
METHODS=['raw','mean5_p1','median5_p1','mean5_p2','median5_p2']
PRIMARY=METHODS[1:]
CHANNELS=list(UNITS)
REFERENCE='20260930T101337+0200'

def apply_filter(values,method):
    """Columns independent; windows overlap, stay within phase, restore original mask."""
    values=np.asarray(values,float); original=np.isfinite(values)
    if method=='raw': return values.copy()
    kind,passes=method.split('_p'); result=pd.DataFrame(values)
    for _ in range(int(passes)):
        rolling=result.rolling(5,center=True,min_periods=1)
        result=(rolling.mean() if kind=='mean5' else rolling.median()).where(original)
    return result.to_numpy()

def grid_matrix(x,y,grid,h):
    """Same median-by-bin operator for EVERY candidate, with at least two finite bins."""
    valid=np.isfinite(x)&(x>=grid[0])&(x<=grid[-1])
    f=pd.DataFrame(y[valid]);f['__bin']=np.round(x[valid]/h)*h
    med=f.groupby('__bin').median()
    out=np.full((len(grid),y.shape[1]),np.nan)
    for j in range(y.shape[1]):
        b=med[j].dropna()
        if len(b)>=2:
            inside=(grid>=b.index.min())&(grid<=b.index.max())
            out[inside,j]=np.interp(grid[inside],b.index,b.to_numpy())
    return out

def scenarios():
    rows=[]
    def add(name,family,**kw): rows.append(dict(scenario=name,family=family,**kw))
    add('clean','clean')
    for sd in [.05,.15]: add(f'gaussian_{sd}','gaussian',sd=sd)
    for amp in [.5,1.,3.]:
        for sign in [-1.,1.]: add(f'impulse_{sign*amp:g}','impulse',artifact=amp*sign,length=1)
    for n in [2,3]: add(f'burst_{n}','burst',artifact=1.,length=n)
    for amp in [.05,.15]: add(f'saw_{amp}','saw',saw=amp)
    for amp in [.1,.3]:
        for duration in [5.,15.,30.]:
            add(f'peak_{amp}_{duration:g}','true_peak',feature='peak',amp=amp,duration=duration)
            add(f'oscillation_{amp}_{duration:g}','true_oscillation',feature='oscillation',amp=amp,duration=duration)
        add(f'border_peak_{amp}','border',feature='peak',amp=amp,duration=15.,border=True)
    add('slope_change','slope',feature='slope',amp=.3)
    for mode in ['random10','gap15','truncated10']: add(mode,'missing',missing=mode,sd=.05)
    return rows

def baseline(channel,u,phase):
    """Explicit mathematical templates, not fitted latent clinical trajectories."""
    if phase=='recovery':
        if channel in ['HRR_L','SpO2']: return 1-np.exp(-3*u)
        if channel=='Pdia': return .55+.15*np.exp(-3*u)
        if channel=='RPM': return .2+.1*np.exp(-6*u)
        return np.exp(-3*u)
    if channel in ['VO2kg','VO2','Psist']: return u
    if channel in ['VCO2','VE','BF','RER']: return .6*u+np.maximum(u-.55,0)*(.4/.45)
    if channel in ['HR','O2pulse']: return (1-np.exp(-2*u))/(1-np.exp(-2))
    if channel in ['EqO2','EqCO2','PETO2']: return ((u-.35)**2)/(.65**2)
    if channel=='PETCO2': return (1-np.exp(-5*u))-.8*np.maximum(u-.65,0)
    if channel=='HRR_L': return 1-u
    if channel=='SpO2': return .8-.3*u
    if channel=='Pdia': return .4+.2*u
    if channel=='RPM': return .6+.1*np.sin(2*np.pi*u)
    raise KeyError(channel)

def make_truth(channel,u,phase,scenario,center,equiv_duration,session):
    b=baseline(channel,u,phase)
    # Known session change; synthetic, independent of clinical CPET1/2 differences.
    if session==1: b=b+.05*(u if phase=='ramp' else np.exp(-3*u))
    f=np.zeros_like(u);sign=-1 if channel in ['SpO2','HRR_L','PETCO2'] else 1
    if scenario.get('feature')=='peak':
        width=scenario['duration']/max(equiv_duration,1.)
        f=sign*scenario['amp']*np.exp(-4*np.log(2)*((u-center)/width)**2)
    elif scenario.get('feature')=='oscillation':
        period=scenario['duration']/max(equiv_duration,1.)
        f=scenario['amp']*np.sin(2*np.pi*(u-center)/period)
    elif scenario.get('feature')=='slope': f=scenario['amp']*np.maximum(u-center,0)
    return b+f,b,f

def metric_columns(pred,truth):
    ok=np.isfinite(pred)&np.isfinite(truth);n=ok.sum(axis=0)
    d=np.where(ok,pred-truth,0.)
    denom=np.where(n>0,n,np.nan)
    return dict(rmse=np.sqrt((d*d).sum(axis=0)/denom),mae=abs(d).sum(axis=0)/denom,bias=d.sum(axis=0)/denom,coverage=n/len(pred))

def collapse_reps(meta,metrics,reps=5):
    f=meta.copy()
    for name,arr in metrics.items():
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',RuntimeWarning)
            f[name]=np.nanmean(np.asarray(arr).reshape(len(f),reps),axis=1)
    f['replicates']=reps
    return f

def prepare():
    run=HERE/'runs'/('smoothing_'+datetime.now().astimezone().strftime('%Y%m%dT%H%M%S%z'))
    run.mkdir(parents=True,exist_ok=False)
    for d in ['synthetic','real','figures','tables','text','selected']: (run/d).mkdir()
    source=json.loads((HERE/'runs'/REFERENCE/'input_hash_audit.json').read_text())
    assert all(sha(ROOT/r['path'])==r['actual'] for r in source),'Historical inputs changed'
    cfg=dict(run_id=run.name,created_at=datetime.now().astimezone().isoformat(),root_seed=SEED,replicates=5,methods=METHODS,primary_candidates=PRIMARY,channels=CHANNELS,window=5,center=True,stride=1,min_periods=1,restore_original_mask_after_each_pass=True,phase_local=True,bin_operator='median; nearest ties-to-even; >=2 occupied bins; no boundary extrapolation',reference_run=REFERENCE,truth='analytic channel-specific mathematical functions of phase coordinate; not fitted or validated physiological ground truth',scale='per-channel median of record-phase empirical q95-q05, floor only if entire channel scale zero; used for units, not filter tuning',sampling='all 60 paired empirical timestamp/workload/missingness patterns; 120 records; two phases',real_reference='raw QC export for impact only; not noise-free truth',decision=dict(primary='mean normalized RMSE: 25% native and 75% gridded; mean within each scenario family, equal family weights, then equal phases/channels, participant means',candidates=PRIMARY,tie='exact numerical tie within 1e-12: fewer passes, then mean before median (literature precedent); not equivalence',uncertainty='2000 participant-cluster bootstrap resamples; sessions/replicates stay together',sensitivity='native-only, grid-only, clean-feature families, leave-one-family-out, core7, exclude BP and RPM',feature_checks='report peak amplitude/area/width/location; do not discard adverse cases',interpretation='single global operational choice among four specified candidates; no universal clinical-optimality claim'),scope='User explicitly requests mean/median then 2-pass comparison, per-variable evaluation, one selected method and supplementary justification',scenarios=scenarios())
    dump(run/'analysis_config.json',cfg)
    dump(run/'input_manifest.json',source)
    # Freeze code and rule BEFORE assessing candidate outcomes.
    (run/'study_code_frozen.py').write_bytes(Path(__file__).read_bytes())
    (HERE/'LATEST_SMOOTHING_RUN.txt').write_text(run.name,encoding='utf-8')
    return run

def tests(run):
    checks=[]
    x=np.arange(12,dtype=float)[:,None];x[4]=np.nan
    for m in METHODS:
        y=apply_filter(x,m)
        assert y.shape==x.shape and np.array_equal(np.isnan(x),np.isnan(y))
        checks.append(dict(test='shape_and_original_mask',method=m,passed=True))
    y=np.zeros((11,1));y[5]=100
    assert apply_filter(y,'mean5_p1')[5,0]==20
    assert apply_filter(y,'median5_p1')[5,0]==0
    # Two means are a triangular effective window, not a second independent estimate.
    assert np.isclose(apply_filter(y,'mean5_p2')[5,0],20)
    assert np.count_nonzero(apply_filter(y,'mean5_p2'))==9
    for n in [0,1,2]:
        a=grid_matrix(np.arange(n,dtype=float),np.ones((n,1)),np.arange(5,dtype=float),1)
        assert np.isfinite(a).sum()==(2 if n==2 else 0)
        checks.append(dict(test=f'{n}_bin_contract',method='all',passed=True))
    for m in PRIMARY:
        constant=apply_filter(np.ones((12,3))*7,m)
        assert np.allclose(constant,7);checks.append(dict(test='constant_preservation',method=m,passed=True))
    pd.DataFrame(checks).to_csv(run/'contract_tests.csv',index=False)

def load_data():
    data={}
    for qc,suffix in [('main','mainHardNA'),('sens','sensHardPlusSusNA')]:
        f=pd.read_csv(D/f'cpet_long_all_patients_v3_{suffix}.csv')
        data[qc]={k:g.sort_values('t') for k,g in f.groupby(['patient_id','test'])}
    master=pd.read_csv(S/'step3_master_by_test_v3.csv').set_index(['patient_id','test'])
    return data,master

def segment(g,m,phase):
    if phase=='ramp':
        s=g.loc[(g.t>=m.t_ramp_start)&(g.t<=m.t_peak)].copy();x=100*s.W.cummax().to_numpy()/m.W_peak;hi=100.;h=1.
    else:
        s=g.loc[(g.t>=m.rec_anchor_t)&(g.t<=m.rec_anchor_t+180)].copy();x=s.t.to_numpy()-m.rec_anchor_t;hi=180.;h=5.
    return s,x,np.arange(0,hi+.1,h),h

def run_real(run,data,master):
    rows=[];audit=[];samples=[];scales=[]
    ids=sorted({k[0] for k in data['main']})
    for qc in data:
        for phase in ['ramp','recovery']:
            arrays={method:np.full((60,2,17,101 if phase=='ramp' else 37),np.nan) for method in METHODS}
            for pi,pid in enumerate(ids):
                for si,test in enumerate(['CPET1','CPET2']):
                    s,x,grid,h=segment(data[qc][pid,test],master.loc[pid,test],phase);y=s[CHANNELS].to_numpy(float);t=s.t.to_numpy(float)
                    for method in METHODS:
                        filtered=apply_filter(y,method);mapped=grid_matrix(x,filtered,grid,h);arrays[method][pi,si]=mapped.T
                        assert np.array_equal(np.isfinite(y),np.isfinite(filtered))
                        for ci,c in enumerate(CHANNELS):
                            a=filtered[:,ci];b=y[:,ci];valid=np.isfinite(b);q=np.nanquantile(b,[.05,.95]) if valid.any() else [np.nan,np.nan]
                            if qc=='main' and method=='raw': scales.append(dict(channel=c,phase=phase,patient_id=pid,test=test,offset=q[0],amplitude=q[1]-q[0]))
                            rows.append(dict(qc=qc,phase=phase,patient_id=pid,test=test,channel=c,unit=UNITS[c],method=method,n_native=len(s),n_valid=int(valid.sum()),peak=float(np.nanmax(a)) if valid.any() else np.nan,trough=float(np.nanmin(a)) if valid.any() else np.nan,peak_time_s=float(t[np.nanargmax(a)]) if valid.any() else np.nan,rmse_vs_raw=float(np.sqrt(np.mean((a[valid]-b[valid])**2))) if valid.any() else np.nan,area_on_valid_times=float(np.trapz(a[valid],t[valid])) if valid.sum()>1 else np.nan,defined_grid=int(np.isfinite(mapped[:,ci]).sum())))
                        audit.append(dict(qc=qc,phase=phase,patient_id=pid,test=test,method=method,rows_in=len(y),rows_out=len(filtered),missing_in=int(np.isnan(y).sum()),missing_out=int(np.isnan(filtered).sum()),original_masks_equal=True))
                        if pi==0 and si==0 and qc=='main':
                            ci=CHANNELS.index('VCO2')
                            samples.extend(dict(phase=phase,method=method,t=float(tt),x=float(xx),value_original=float(bb),value_filtered=float(aa)) for tt,xx,bb,aa in zip(t,x,y[:,ci],filtered[:,ci]))
            for method,a in arrays.items():
                delta=a[:,1]-a[:,0]
                np.savez_compressed(run/'real'/f'{qc}_{phase}_{method}.npz',X_cpet1=a[:,0],X_cpet2=a[:,1],X_delta=delta,M_delta=np.isfinite(delta).astype('uint8'),grid=grid,patient_ids=np.array(ids,dtype=str),var_names=np.array(CHANNELS,dtype=str))
            print('Real matrices',qc,phase,flush=True)
    pd.DataFrame(rows).to_parquet(run/'real/variable_metrics.parquet',index=False)
    pd.DataFrame(audit).to_csv(run/'real/row_mask_audit.csv',index=False)
    pd.DataFrame(samples).to_csv(run/'real/first_template_VCO2_internal.csv',index=False)
    sf=pd.DataFrame(scales);sf.to_csv(run/'real/channel_scale_inputs.csv',index=False)
    scale=sf.groupby('channel')[['offset','amplitude']].median().reindex(CHANNELS)
    assert (scale.amplitude>0).all(),'Cannot normalize constant/absent channel; inspect before proceeding'
    scale.to_csv(run/'channel_scales.csv')
    return scale

def run_synthetic(run,data,master,scale):
    scenario_list=scenarios();ns=len(scenario_list);reps=5;ncol=ns*reps
    metadata=pd.DataFrame([dict(scenario=s['scenario'],family=s['family']) for s in scenario_list])
    config=json.loads((run/'analysis_config.json').read_text());assert config['scenarios']==scenario_list
    ids=sorted({k[0] for k in data['main']});manifest=[];examples=[];allrows=[];pairrows=[];featureframes=[]
    for pi,pid in enumerate(ids):
        for phase in ['ramp','recovery']:
            paired={}
            for si,test in enumerate(['CPET1','CPET2']):
                s,x,grid,h=segment(data['main'][pid,test],master.loc[pid,test],phase)
                u=x/grid[-1];ug=grid/grid[-1];t=s.t.to_numpy(float)
                duration=float(t[-1]-t[0]) if len(t)>1 else 0.;span=float(u.max()-u.min()) if len(u) else 0.
                equiv=duration/span if span>0 else 180.
                center0=float(u.min()+.55*span) if len(u) else .55
                for ci,c in enumerate(CHANNELS):
                    valid=np.isfinite(s[c].to_numpy(float));A=float(scale.loc[c,'amplitude']);offset=float(scale.loc[c,'offset'])
                    truth=np.empty((len(u),ncol));tg=np.empty((len(grid),ncol));base=np.empty_like(truth);bg=np.empty_like(tg);obs=np.empty_like(truth)
                    col=0;centers=[]
                    for sj,scenario in enumerate(scenario_list):
                        center=float(u.max()-7.5/max(equiv,1)) if scenario.get('border') else center0;centers.append(center)
                        tr,ba,fe=make_truth(c,u,phase,scenario,center,equiv,si);gt,gb,gf=make_truth(c,ug,phase,scenario,center,equiv,si)
                        for rep in range(reps):
                            seed=[SEED,pi,si,0 if phase=='ramp' else 1,ci,sj,rep]
                            rng=np.random.default_rng(np.random.SeedSequence(seed));noisy=tr.copy()
                            if scenario.get('sd'): noisy+=rng.normal(0,scenario['sd'],len(u))
                            if scenario.get('artifact'):
                                loc=int(np.argmin(abs(u-center))) if len(u) else 0;noisy[loc:loc+scenario['length']]+=scenario['artifact']
                            if scenario.get('saw'): noisy+=scenario['saw']*np.where(np.arange(len(u))%2,1,-1)
                            mask=valid.copy()
                            if scenario.get('missing')=='random10': mask&=rng.random(len(u))>=.1
                            if scenario.get('missing')=='gap15': mask&=abs(u-center)*equiv>7.5
                            if scenario.get('missing')=='truncated10': mask&=u<=u.max()-.1*span
                            noisy[~mask]=np.nan
                            truth[:,col]=tr;base[:,col]=ba;tg[:,col]=gt;bg[:,col]=gb;obs[:,col]=noisy;col+=1
                    context=dict(template_id=f'T{pi+1:03d}',session=test,phase=phase,channel=c,unit=UNITS[c],amplitude_scale=A)
                    manifest.append(dict(**context,n_native=len(u),n_observed=int(valid.sum()),coordinate_min=float(u.min()) if len(u) else np.nan,coordinate_max=float(u.max()) if len(u) else np.nan,equivalent_duration_s=equiv,seed_components_prefix=f'{SEED},{pi},{si},{0 if phase=="ramp" else 1},{ci}',scenarios=ns,replicates=reps))
                    for method in METHODS:
                        filtered=apply_filter(obs,method);mapped=grid_matrix(x,filtered,grid,h)
                        assert np.array_equal(np.isfinite(filtered),np.isfinite(obs))
                        nm=metric_columns(filtered,truth);gm=metric_columns(mapped,tg)
                        met={f'native_{k}':v for k,v in nm.items()};met.update({f'grid_{k}':v for k,v in gm.items()})
                        met['native_rmse_units']=nm['rmse']*A;met['grid_rmse_units']=gm['rmse']*A
                        fm=collapse_reps(metadata,met)
                        for k,v in context.items(): fm[k]=v
                        fm['method']=method;allrows.append(fm)
                        # Feature estimates in the true-peak scenarios; do not drop poorly sampled peaks.
                        feature_rows=[]
                        for sj,scenario in enumerate(scenario_list):
                            if scenario.get('feature')!='peak': continue
                            center=centers[sj];width=scenario['duration']/max(equiv,1);local=abs(u-center)<=2*width
                            sign=-1 if c in ['SpO2','HRR_L','PETCO2'] else 1
                            for rep in range(reps):
                                j=sj*reps+rep;ok=local&np.isfinite(filtered[:,j]);n=int(ok.sum())
                                residual=sign*(filtered[ok,j]-base[ok,j]);true_residual=sign*(truth[ok,j]-base[ok,j]);xx=u[ok]*equiv
                                height=float(np.max(residual)) if n else np.nan
                                time=float(xx[np.argmax(residual)]-center*equiv) if n else np.nan
                                above=residual>=scenario['amp']/2
                                # Width uses the connected above-half-height segment around the estimated maximum.
                                width_est=np.nan
                                if n>=2 and np.any(above):
                                    peak=int(np.argmax(residual));left=right=peak
                                    while left>0 and above[left-1]:left-=1
                                    while right+1<n and above[right+1]:right+=1
                                    width_est=float(xx[right]-xx[left])
                                feature_rows.append(dict(**context,method=method,scenario=scenario['scenario'],replicate=rep,n_local=n,peak_height_ratio=height/scenario['amp'],peak_time_error_s=abs(time),feature_area_error_norm=float(abs(np.trapz(residual-true_residual,xx))/(scenario['amp']*scenario['duration'])) if n>=2 else np.nan,width_error_s=abs(width_est-scenario['duration']) if np.isfinite(width_est) else np.nan,adequately_sampled=bool(n>=3 and np.count_nonzero(ok&(abs(u-center)<=width/2))>=2)))
                        if feature_rows:
                            featureframes.append(pd.DataFrame(feature_rows))
                        paired[si,c,method]=(mapped,tg)
                        if pi==0 and si==0 and phase=='ramp' and c=='VCO2':
                            for name in ['impulse_1','peak_0.3_15','gaussian_0.15','burst_3']:
                                sj=next(k for k,v in enumerate(scenario_list) if v['scenario']==name);j=sj*reps
                                examples.extend(dict(scenario=name,method=method,coordinate=float(xx),truth=float(yy),observed=float(oo),filtered=float(ff)) for xx,yy,oo,ff in zip(x,truth[:,j],obs[:,j],filtered[:,j]))
            for c in CHANNELS:
                for method in METHODS:
                    p1,t1=paired[0,c,method];p2,t2=paired[1,c,method];pm=metric_columns(p2-p1,t2-t1)
                    f=collapse_reps(metadata,{f'paired_grid_{k}':v for k,v in pm.items()})
                    f['template_id']=f'T{pi+1:03d}';f['phase']=phase;f['channel']=c;f['method']=method;pairrows.append(f)
        if (pi+1)%5==0:print('Synthetic templates',pi+1,'/60',flush=True)
        # Bounded memory and consolidated per-template metrics; no million-row Python dict list.
        pd.concat(allrows,ignore_index=True).to_parquet(run/'synthetic'/f'metrics_template_{pi:02d}.parquet',index=False);allrows.clear()
        pd.concat(pairrows,ignore_index=True).to_parquet(run/'synthetic'/f'paired_template_{pi:02d}.parquet',index=False);pairrows.clear()
        pd.concat(featureframes,ignore_index=True).to_parquet(run/'synthetic'/f'features_template_{pi:02d}.parquet',index=False);featureframes.clear()
    pd.DataFrame(manifest).to_csv(run/'synthetic/template_manifest.csv',index=False)
    pd.DataFrame(scenario_list).to_csv(run/'synthetic/scenario_manifest.csv',index=False)
    pd.DataFrame(examples).to_csv(run/'synthetic/example_curves.csv',index=False)
    dump(run/'simulation_complete.json',dict(completed_at=datetime.now().astimezone().isoformat(),templates=60,records=120,phases=2,channels=17,methods=METHODS,scenarios=ns,replicates=5,replicate_aggregation='mean per record/channel/scenario; not independent participants',seed=SEED))

def main():
    p=argparse.ArgumentParser();p.add_argument('--run');p.add_argument('--stage',choices=['all','real','synthetic'],default='all');args=p.parse_args()
    run=Path(args.run) if args.run else prepare();tests(run);data,master=load_data()
    if args.stage in ['all','real']:scale=run_real(run,data,master)
    else:scale=pd.read_csv(run/'channel_scales.csv',index_col='channel')
    if args.stage in ['all','synthetic']:run_synthetic(run,data,master,scale)
    print(run,flush=True)

if __name__=='__main__':main()
