from pathlib import Path
import os,sys,csv,json,hashlib
root=Path(__file__).resolve().parents[1];base=root.parent;ev=root/'evidence';run=root/'run'
sys.path.insert(0,str(base/'Schutz_2018_phase1'/'deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
protocol=json.loads((ev/'protocol.json').read_text(encoding='utf8'));ph=hashlib.sha256((ev/'protocol.json').read_bytes()).hexdigest()
assert ph==(ev/'protocol_sha256.txt').read_text(encoding='utf8')
for name,h in protocol['inputs_sha256'].items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==h
assert hashlib.sha256((ev/'frozen_variants.json').read_bytes()).hexdigest()==protocol['variant_file_sha256']
variants=json.loads((ev/'frozen_variants.json').read_text(encoding='utf8'));models=list(protocol['models']);baselines=protocol['baselines'];estimators=models+baselines
digest=lambda ids:hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
def metric(yy,pp,ww):
 ww=ww/ww.sum();ee=yy-pp;mu=np.sum(ww*yy);sst=np.sum(ww*(yy-mu)**2);mse=np.sum(ww*ee**2)
 return {'MAE':float(np.sum(ww*np.abs(ee))),'RMSE':float(np.sqrt(mse)),'R2':float(1-mse/sst) if sst>1e-15 else 'NA'}
metric_rows=[];result={};fitaudits=[]
for variant,v in variants.items():
 df=pd.DataFrame(v['rows']).set_index('observation_id',drop=False);df.logRR=df.logRR.astype(float)
 ids=df.index.tolist();y=df.logRR.to_numpy(float);N=len(ids);studies=sorted(df.study_id.unique(),key=int);S=len(studies)
 inv=np.array([studies.index(s) for s in df.study_id]);ns=np.bincount(inv)
 tasks=sorted(set(f['task'] for f in v['splits']));pred={t:{e:np.full(N,np.nan) for e in estimators} for t in tasks}
 for i,f in enumerate(v['splits']):
  r=json.loads((run/variant/f'fold_{i:04d}.json').read_text(encoding='utf8'));tr,te,pu=[f[role+'_ids'] for role in ('train','test','purge')]
  assert r['protocol_sha256']==ph and r['predictors']==v['predictors'] and r['fold_id']==f['fold_id'] and r['test_ids']==te and r['model_seed']==20261006+f['seed_index'] and r['model_fits']==3
  assert r['preprocessor_fit_ids_sha256']==digest(tr) and not set(df.loc[tr,'study_id'])&set(df.loc[te+pu,'study_id'])
  assert len(set(tr+te+pu))==len(tr+te+pu)==N
  for role,lst in zip(('train','test','purge'),(tr,te,pu)):assert digest(lst)==f[role+'_ids_sha256']==r[role+'_ids_sha256']
  train=df.loc[tr];a=train[v['numeric_predictors']].replace('NA',np.nan).astype(float).to_numpy();med=np.nanmedian(a,axis=0);np.testing.assert_allclose(med,r['numeric_train_medians'],rtol=0,atol=1e-14)
  ai=np.where(np.isnan(a),med,a);_,ii,nn=np.unique(train.study_id,return_inverse=True,return_counts=True);w=len(tr)/(len(nn)*nn[ii])
  mu=np.average(ai,axis=0,weights=w);var=np.average((ai-mu)**2,axis=0,weights=w);sc=np.sqrt(var);sc[sc==0]=1
  np.testing.assert_allclose(mu,r['numeric_train_scaler_mean'],rtol=1e-12,atol=1e-12);np.testing.assert_allclose(sc,r['numeric_train_scaler_scale'],rtol=1e-12,atol=1e-12)
  assert set(r['categorical_training_vocabulary'])==set(v['categorical_predictors'])
  for c in v['categorical_predictors']:
   lv=sorted(set(train[c].astype(str).replace('NA','__MISSING__'))|{'__MISSING__'});assert lv==r['categorical_training_vocabulary'][c]
   assert int((~df.loc[te,c].astype(str).replace('NA','__MISSING__').isin(lv)).sum())==r['test_unknown_counts'][c]
  bs={'zero_response':0.,'training_observation_mean':float(train.logRR.mean()),'training_study_equal_mean':float(train.groupby('study_id').logRR.mean().mean())}
  ix=df.index.get_indexer(te)
  for e in estimators:
   pp=np.array(r['predictions'][e]);assert np.isnan(pred[f['task']][e][ix]).all() and len(pp)==len(ix) and np.isfinite(pp).all();pred[f['task']][e][ix]=pp
   if e in bs:
    np.testing.assert_allclose(r['baseline_constants'][e],bs[e],rtol=1e-14,atol=1e-14);np.testing.assert_allclose(pp,bs[e],rtol=0,atol=1e-14)
  fitaudits.append({'variant':variant,'fold_id':f['fold_id'],'preprocessor_fit_ids_sha256':r['preprocessor_fit_ids_sha256'],'study_overlap':0,'frozen_predictors_only':True,'training_only_statistics_reconstructed':True})
 assert all(np.isfinite(pp).all() for t in tasks for pp in pred[t].values())
 np.savez_compressed(run/variant/'oof_predictions.npz',observation_ids=np.array(ids),study_ids=np.array(df.study_id),logRR=y,**{t+'__'+e:pp for t in tasks for e,pp in pred[t].items()})
 rng=np.random.default_rng(20261006);sampled=rng.integers(0,S,size=(2000,S));C=np.stack([np.bincount(row,minlength=S) for row in sampled]);assert np.all(C.sum(axis=1)==S)
 if S==169:
  original=np.load(base/'Schutz_2018_phase1'/'run'/'bootstrap_draws.npz');assert np.array_equal(C,original['study_multiplicity']) and np.array_equal(np.array(studies),original['study_ids'])
 sy=np.bincount(inv,weights=y);sy2=np.bincount(inv,weights=y*y);points={};draws={};archive={'study_multiplicity':C,'sampled_study_index':sampled,'study_ids':np.array(studies)};support={}
 for t in tasks:
  points[t]={};draws[t]={}
  for e,pp in pred[t].items():
   points[t][e]={};draws[t][e]={};ae=np.bincount(inv,weights=np.abs(y-pp));se=np.bincount(inv,weights=(y-pp)**2)
   for wt,w in (('observation',np.ones(N)),('study_equal',1/ns[inv])):
    mt=metric(y,pp,w);points[t][e][wt]=mt
    expected=[mean_absolute_error(y,pp,sample_weight=w),np.sqrt(mean_squared_error(y,pp,sample_weight=w)),r2_score(y,pp,sample_weight=w)]
    np.testing.assert_allclose([mt[m] for m in ('MAE','RMSE','R2')],expected,rtol=1e-12,atol=1e-12)
    if wt=='observation':den=C@ns;mae=C@ae/den;mse=C@se/den;mu=C@sy/den;vv=C@sy2/den-mu*mu
    else:den=S;mae=C@(ae/ns)/den;mse=C@(se/ns)/den;mu=C@(sy/ns)/den;vv=C@(sy2/ns)/den-mu*mu
    val={'MAE':mae,'RMSE':np.sqrt(mse),'R2':1-mse/vv};draws[t][e][wt]=val
    for m,arr in val.items():assert np.isfinite(arr).all();archive[f'{t}__{e}__{wt}__{m}']=arr
    for k in (0,17,1999):
     ix=np.concatenate([np.flatnonzero(inv==s) for s in sampled[k]]);ww=np.ones(len(ix)) if wt=='observation' else 1/ns[inv[ix]];mm=metric(y[ix],pp[ix],ww)
     np.testing.assert_allclose([val[m][k] for m in ('MAE','RMSE','R2')],[mm[m] for m in ('MAE','RMSE','R2')],rtol=1e-12,atol=1e-12)
  if t in ('A','B'):
   support[t]={}
   for e in models:
    pt=points[t][e];checks={}
    for b in baselines[1:]:
     upper=float(np.quantile(draws[t][e]['study_equal']['RMSE']-draws[t][b]['study_equal']['RMSE'],1-.05/12));checks['adjusted_delta_RMSE_upper_vs_'+b]=upper<0
     checks['study_equal_MAE_no_worse_vs_'+b]=pt['study_equal']['MAE']<=points[t][b]['study_equal']['MAE']
     for m in ('MAE','RMSE'):checks['observation_'+m+'_no_worse_vs_'+b]=pt['observation'][m]<=points[t][b]['observation'][m]
    checks['study_equal_R2_positive']=pt['study_equal']['R2']>0
    for wt in ('observation','study_equal'):checks[wt+'_RMSE_no_worse_vs_zero']=pt[wt]['RMSE']<=points[t]['zero_response'][wt]['RMSE']
    support[t][e]={'supported':all(checks.values()),'checks':checks}
 for t in tasks:
  nf=sum(f['task']==t for f in v['splits'])
  for e in estimators:
   for wt in ('observation','study_equal'):
    for m in ('MAE','RMSE','R2'):
     arr=draws[t][e][wt][m];lo,hi=np.quantile(arr,[.025,.975])
     metric_rows.append({'scenario':variant,'task':t,'kind':'absolute','estimator':e,'baseline':'NA','weighting':wt,'metric':m,'point':points[t][e][wt][m],'ci95_lower':float(lo),'ci95_upper':float(hi),'adjusted_one_sided_upper':'NA','n_observations':N,'n_studies':S,'n_predictors':len(v['predictors']),'outer_folds':nf,'bootstrap_n':2000,'seed':20261006,'phase1_gate_supported':str(support[t][e]['supported']).lower() if t in support and e in models else 'NA'})
  for e in models:
   for b in baselines:
    for wt in ('observation','study_equal'):
     for m in ('MAE','RMSE','R2'):
      arr=draws[t][e][wt][m]-draws[t][b][wt][m];lo,hi=np.quantile(arr,[.025,.975])
      upper=float(np.quantile(arr,1-.05/12)) if t in ('A','B') and b!='zero_response' and wt=='study_equal' and m=='RMSE' else 'NA'
      metric_rows.append({'scenario':variant,'task':t,'kind':'paired_delta_model_minus_baseline','estimator':e,'baseline':b,'weighting':wt,'metric':m,'point':float(points[t][e][wt][m]-points[t][b][wt][m]),'ci95_lower':float(lo),'ci95_upper':float(hi),'adjusted_one_sided_upper':upper,'n_observations':N,'n_studies':S,'n_predictors':len(v['predictors']),'outer_folds':nf,'bootstrap_n':2000,'seed':20261006,'phase1_gate_supported':str(support[t][e]['supported']).lower() if t in support else 'NA'})
 np.savez_compressed(run/variant/'bootstrap_draws.npz',**archive)
 result[variant]={'points':points,'support':support,'cohort_n':N,'studies_n':S,'bootstrap_archive_sha256':hashlib.sha256((run/variant/'bootstrap_draws.npz').read_bytes()).hexdigest(),'predictions_archive_sha256':hashlib.sha256((run/variant/'oof_predictions.npz').read_bytes()).hexdigest()}

# Carry original summary for side-by-side comparison without any Phase1 refits.
p1=base/'Schutz_2018_phase1';original_boot=list(csv.DictReader((p1/'phase1_bootstrap.tsv').open(encoding='utf8',newline=''),delimiter='\t'))
for r in original_boot:
 metric_rows.append({'scenario':'phase1_main_1723_12_reference','task':r['task'],'kind':r['kind'],'estimator':r['estimator'],'baseline':r['baseline'],'weighting':r['weighting'],'metric':r['metric'],'point':float(r['point']),'ci95_lower':float(r['ci95_lower']),'ci95_upper':float(r['ci95_upper']),'adjusted_one_sided_upper':float(r['adjusted_one_sided_upper']) if r['adjusted_one_sided_upper']!='NA' else 'NA','n_observations':1723,'n_studies':169,'n_predictors':12,'outer_folds':{'A':169,'B':40,'C':6}[r['task']],'bootstrap_n':2000,'seed':20261006,'phase1_gate_supported':'false' if r['task'] in ('A','B') and r['estimator'] in models else 'NA'})
assert len(metric_rows)==720

# Descriptive audit uses unchanged Phase1 OOF only; strata were frozen from metadata.
main=variants['main_1723_crop_blind_11']['rows'];df=pd.DataFrame(main).set_index('observation_id',drop=False);y=df.logRR.astype(float).to_numpy();N=len(df);studies=sorted(df.study_id.unique(),key=int);inv=np.array([studies.index(s) for s in df.study_id]);ns=np.bincount(inv)
oldpred={t:{e:np.full(N,np.nan) for e in estimators} for t in ('A','B','C')};seen=set()
for r in csv.DictReader((p1/'phase1_predictions.tsv').open(encoding='utf8',newline=''),delimiter='\t'):
 i=df.index.get_loc(r['observation_id']);key=(r['task'],r['estimator'],i);assert key not in seen;seen.add(key)
 assert float(r['logRR'])==y[i];oldpred[r['task']][r['estimator']][i]=float(r['prediction'])
assert len(seen)==31014 and all(np.isfinite(pp).all() for t in oldpred for pp in oldpred[t].values())
crop_bands={}
for crop,rs in df.groupby('crop'):
 s=rs.study_id.nunique();n=len(rs);crop_bands[crop]='singleton' if s==1 else 'high_support' if s>=3 and n>=20 else 'low_support'
failure_rows=[]
def auditrow(scope,t,e,grouping,g,wt,m,stat,value,idx):
 failure_rows.append({'scope':scope,'task':t,'estimator':e,'grouping':grouping,'group_id':g,'weighting':wt,'metric':m,'statistic':stat,'value':value,'n_observations':len(idx),'n_studies':int(len(set(inv[idx]))),'n_crops':int(df.iloc[idx].crop.nunique()),'source':'unchanged_Phase1_main12_OOF','descriptive_only':'true'})
for t in oldpred:
 for e,pp in oldpred[t].items():
  groups={'study_size':{'large_ge20':np.flatnonzero(ns[inv]>=20),'small_lt20':np.flatnonzero(ns[inv]<20)},'crop_support':{g:np.flatnonzero(df.crop.map(crop_bands).to_numpy()==g) for g in ('high_support','low_support','singleton')}}
  for grouping,gg in groups.items():
   for g,ix in gg.items():
    nsub=np.bincount(inv[ix],minlength=len(studies))
    for wt,w in (('observation',np.ones(len(ix))),('study_equal',1/nsub[inv[ix]])):
     for m,val in metric(y[ix],pp[ix],w).items():auditrow('metadata_stratum',t,e,grouping,g,wt,m,'point',val,ix)
  perstudy={'MAE':[],'RMSE':[]}
  for s,study in enumerate(studies):
   ix=np.flatnonzero(inv==s);mt=metric(y[ix],pp[ix],np.ones(len(ix)))
   for m in perstudy:perstudy[m].append(mt[m]);auditrow('individual_study',t,e,'study_id',study,'within_study',m,'point',mt[m],ix)
  for m,values in perstudy.items():
   allidx=np.arange(N)
   auditrow('study_error_distribution',t,e,'all_studies','ALL','study_equal',m,'mean_of_study_metric',float(np.mean(values)),allidx)
   auditrow('study_error_distribution',t,e,'all_studies','ALL','study_equal',m,'std_of_study_metric',float(np.std(values,ddof=1)),allidx)
   for q in protocol['failure_audit']['study_error_quantiles']:auditrow('study_error_distribution',t,e,'all_studies','ALL','study_equal',m,'quantile_'+str(q),float(np.quantile(values,q)),allidx)
assert len(failure_rows)==6984
for name,h in protocol['inputs_sha256'].items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==h
data={'phase2_sensitivity_metrics.tsv':metric_rows,'phase2_failure_audit.tsv':failure_rows}
(ev/'delivery_tables.json').write_text(json.dumps(data,ensure_ascii=True,allow_nan=False),encoding='utf8')
summary={'protocol_sha256':ph,'variants':result,'sensitivity_metric_rows':len(metric_rows),'failure_audit_rows':len(failure_rows),'fit_QC':fitaudits,'all_inputs_unchanged':True,'independent_sklearn_point_metrics':True,'bootstrap_expanded_cluster_draws_checked':[0,17,1999],'source_Phase1_OOF_count':31014,'failure_no_refitting':True}
(ev/'evaluation_summary.json').write_text(json.dumps(summary,ensure_ascii=True,indent=2,allow_nan=False),encoding='utf8')
console_summary={'rows':{k:len(v) for k,v in data.items()},'sensitivity_task_supported':{v:{t:any(a['supported'] for a in ss.values()) for t,ss in r['support'].items()} for v,r in result.items()},'sensitivity_study_equal_points':{v:{t:{e:rr['study_equal'] for e,rr in r.items()} for t,r in z['points'].items()} for v,z in result.items()}}
print(json.dumps(console_summary),flush=True)
