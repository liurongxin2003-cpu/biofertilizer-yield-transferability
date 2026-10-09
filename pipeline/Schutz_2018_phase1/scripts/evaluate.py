from pathlib import Path
import os,sys,csv,json,hashlib
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
import pandas as pd

ev=root/'evidence'; run=root/'run'; old=root.parent/'Schutz_2018_phase0'
protocol=json.loads((ev/'protocol.json').read_text(encoding='utf8'))
ph=hashlib.sha256((ev/'protocol.json').read_bytes()).hexdigest()
assert ph==(ev/'protocol_sha256.txt').read_text(encoding='utf8')
for f,h in protocol['input_sha256'].items():assert hashlib.sha256((old/f).read_bytes()).hexdigest()==h
df=pd.read_csv(old/'analysis_matrix.tsv',sep='\t',keep_default_na=False,dtype=str).set_index('observation_id',drop=False)
df['logRR']=df['logRR'].astype(float)
df.study_id=df.study_id.astype(str)
splits=list(csv.DictReader((old/'split_manifest.tsv').open(encoding='utf8',newline=''),delimiter='\t'))
ids=df.index.tolist(); y=df.logRR.to_numpy(float)
studies=sorted(df.study_id.unique(),key=int); inv=np.array([studies.index(s) for s in df.study_id])
ns=np.bincount(inv); S=len(ns); N=len(ids)
models=list(protocol['models']); baselines=protocol['baselines']; estimators=models+baselines
tasks=protocol['designs']; design_task={v:k for k,v in tasks.items()}
pred={d:{e:np.full(N,np.nan) for e in estimators} for d in tasks.values()}
folds={d:np.empty(N,dtype=object) for d in tasks.values()}
prediction_rows=[]; fold_provenance=[]
unpack=lambda s:[] if s=='NA' else s.split(',')
def digest(x):return hashlib.sha256(('\n'.join(x)+'\n').encode()).hexdigest()
for i,f in enumerate(splits):
 r=json.loads((run/f'fold_{i:04d}.json').read_text(encoding='utf8'))
 tr,te,pu=[unpack(f[k+'_ids']) for k in ('train','test','purge')]
 assert r['protocol_sha256']==ph and r['fold_id']==f['fold_id'] and r['test_ids']==te
 assert r['model_seed']==20261006+i and r['model_fits']==3
 assert r['preprocessor_fit_ids_sha256']==digest(tr)
 for role,lst in zip(('train','test','purge'),(tr,te,pu)):assert digest(lst)==f[role+'_ids_sha256']==r[role+'_ids_sha256']
 assert not set(df.loc[tr,'study_id'])&set(df.loc[te+pu,'study_id'])
 assert len(set(tr+te+pu))==1723 and len(tr+te+pu)==1723
 # Independently reconstruct training-only preprocessing and all baseline constants.
 train=df.loc[tr]; a=train[protocol['numeric_predictors']].replace('NA',np.nan).astype(float).to_numpy()
 med=np.nanmedian(a,axis=0); np.testing.assert_allclose(med,r['numeric_train_medians'],rtol=0,atol=1e-14)
 ai=np.where(np.isnan(a),med,a)
 _,ii,nn=np.unique(train.study_id,return_inverse=True,return_counts=True); w=len(tr)/(len(nn)*nn[ii])
 mean=np.average(ai,axis=0,weights=w); var=np.average((ai-mean)**2,axis=0,weights=w)
 scale=np.sqrt(var);scale[scale==0]=1
 np.testing.assert_allclose(mean,r['numeric_train_scaler_mean'],rtol=1e-12,atol=1e-12)
 np.testing.assert_allclose(scale,r['numeric_train_scaler_scale'],rtol=1e-12,atol=1e-12)
 for c in protocol['categorical_predictors']:
  lv=sorted(set(train[c].astype(str).replace('NA','__MISSING__'))|{'__MISSING__'})
  assert lv==r['categorical_training_vocabulary'][c]
  un=int((~df.loc[te,c].astype(str).replace('NA','__MISSING__').isin(lv)).sum())
  assert un==r['test_unknown_counts'][c]
 b={'zero_response':0.,'training_observation_mean':float(train.logRR.mean()),'training_study_equal_mean':float(train.groupby('study_id').logRR.mean().mean())}
 for name,v in b.items():np.testing.assert_allclose(r['baseline_constants'][name],v,rtol=1e-14,atol=1e-14)
 idx=df.index.get_indexer(te);d=f['design'];folds[d][idx]=f['fold_id']
 for e in estimators:
  yp=np.asarray(r['predictions'][e]);assert len(yp)==len(te) and np.isfinite(yp).all()
  assert np.isnan(pred[d][e][idx]).all();pred[d][e][idx]=yp
  if e in b:np.testing.assert_allclose(yp,b[e],rtol=0,atol=1e-14)
  for j,oid in enumerate(te):
   row=df.loc[oid];prediction_rows.append({'task':design_task[d],'estimator':e,'observation_id':oid,'study_id':int(row.study_id),'logRR':float(row.logRR),'prediction':float(yp[j]),'residual_observed_minus_predicted':float(row.logRR-yp[j]),'crop':row.crop,'inoculant_function':str(row.inoculant_function),'fold_id':f['fold_id'],'primary_or_secondary':f['primary_or_secondary'],'train_n':int(f['train_n']),'purge_n':int(f['purge_n']),'train_ids_sha256':r['train_ids_sha256'],'preprocessor_fit_ids_sha256':r['preprocessor_fit_ids_sha256'],'protocol_sha256':ph})
 fold_provenance.append({k:r[k] for k in ('fold_id','test_unknown_counts','transformed_columns','train_n','test_n','purge_n')})
for d in pred:
 for e in pred[d]:assert np.isfinite(pred[d][e]).all()
assert len(prediction_rows)==31014

def metric(yy,pp,ww):
 ww=ww/ww.sum();err=yy-pp; mu=np.sum(ww*yy);sst=np.sum(ww*(yy-mu)**2)
 mse=np.sum(ww*err**2)
 return {'MAE':float(np.sum(ww*np.abs(err))),'RMSE':float(np.sqrt(mse)),'R2':float(1-mse/sst) if sst>1e-15 else 'NA'}
metric_rows=[];point={}
for d in pred:
 point[d]={}
 for e,pp in pred[d].items():
  point[d][e]={}
  for weighting,w in (('observation',np.ones(N)),('study_equal',1/ns[inv])):
   mt=metric(y,pp,w);point[d][e][weighting]=mt
   metric_rows.append({'scope':'overall_oof','task':design_task[d],'estimator':e,'weighting':weighting,'group_id':'ALL','n_observations':N,'n_studies':S,'MAE':mt['MAE'],'RMSE':mt['RMSE'],'R2':mt['R2'],'precision_stratum':'all_frozen_rows'})
 # No averaging per-fold R2: pooled held-out predictions define primary metrics.
for f in splits:
 d=f['design'];idx=df.index.get_indexer(unpack(f['test_ids']));gi=inv[idx];nc=np.bincount(gi,minlength=S);st=int((nc>0).sum())
 for e,pp in pred[d].items():
  for weighting,w in (('observation',np.ones(len(idx))),('study_equal',1/nc[gi])):
   mt=metric(y[idx],pp[idx],w)
   metric_rows.append({'scope':'fold_descriptive','task':design_task[d],'estimator':e,'weighting':weighting,'group_id':f['fold_id'],'n_observations':len(idx),'n_studies':st,**mt,'precision_stratum':'n_studies_ge3_and_obs_ge20' if st>=3 and len(idx)>=20 else 'small_group'})
for d in pred:
 for s,study in enumerate(studies):
  idx=np.flatnonzero(inv==s)
  for e,pp in pred[d].items():
   mt=metric(y[idx],pp[idx],np.ones(len(idx)));mt['R2']='NA'
   metric_rows.append({'scope':'study_descriptive','task':design_task[d],'estimator':e,'weighting':'within_study','group_id':study,'n_observations':len(idx),'n_studies':1,**mt,'precision_stratum':'descriptive_only_R2_not_defined_as_study_score'})

rng=np.random.default_rng(protocol['bootstrap']['seed'])
sampled=rng.integers(0,S,size=(2000,S)); C=np.stack([np.bincount(row,minlength=S) for row in sampled])
assert C.shape==(2000,169) and np.all(C.sum(axis=1)==169)
sy=np.bincount(inv,weights=y);sy2=np.bincount(inv,weights=y*y)
draws={};archive={'study_multiplicity':C,'sampled_study_index':sampled,'study_ids':np.array(studies)}
bootstrap_rows=[]
for d in pred:
 draws[d]={}
 for e,pp in pred[d].items():
  ae=np.bincount(inv,weights=np.abs(y-pp));se=np.bincount(inv,weights=(y-pp)**2);draws[d][e]={}
  for wt in ('observation','study_equal'):
   if wt=='observation':den=C@ns;mae=C@ae/den;mse=C@se/den;mu=C@sy/den;vv=C@sy2/den-mu*mu
   else:den=S;mae=C@(ae/ns)/den;mse=C@(se/ns)/den;mu=C@(sy/ns)/den;vv=C@(sy2/ns)/den-mu*mu
   vals={'MAE':mae,'RMSE':np.sqrt(mse),'R2':1-mse/vv};draws[d][e][wt]=vals
   for m,v in vals.items():
    assert np.isfinite(v).all();archive[f'{design_task[d]}__{e}__{wt}__{m}']=v
    lo,hi=np.quantile(v,[.025,.975])
    bootstrap_rows.append({'kind':'absolute','task':design_task[d],'estimator':e,'baseline':'NA','weighting':wt,'metric':m,'point':point[d][e][wt][m],'ci95_lower':float(lo),'ci95_upper':float(hi),'adjusted_one_sided_upper':'NA','bootstrap_n':2000,'seed':20261006,'direction':'lower_better' if m!='R2' else 'higher_better'})
 for e in models:
  for b in baselines:
   for wt in ('observation','study_equal'):
    for m in ('MAE','RMSE','R2'):
     v=draws[d][e][wt][m]-draws[d][b][wt][m];lo,hi=np.quantile(v,[.025,.975])
     upper=float(np.quantile(v,1-.05/12)) if design_task[d] in ('A','B') and wt=='study_equal' and m=='RMSE' and b!='zero_response' else 'NA'
     bootstrap_rows.append({'kind':'paired_delta_model_minus_baseline','task':design_task[d],'estimator':e,'baseline':b,'weighting':wt,'metric':m,'point':float(point[d][e][wt][m]-point[d][b][wt][m]),'ci95_lower':float(lo),'ci95_upper':float(hi),'adjusted_one_sided_upper':upper,'bootstrap_n':2000,'seed':20261006,'direction':'negative_improvement' if m!='R2' else 'positive_improvement'})
assert len(bootstrap_rows)==270
np.savez_compressed(run/'bootstrap_draws.npz',**archive)

# Independent expanded-cluster verification of multiplicity and weighted point metrics.
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
for d in pred:
 for e,pp in pred[d].items():
  for wt,w in (('observation',np.ones(N)),('study_equal',1/ns[inv])):
   expected=[mean_absolute_error(y,pp,sample_weight=w),np.sqrt(mean_squared_error(y,pp,sample_weight=w)),r2_score(y,pp,sample_weight=w)]
   np.testing.assert_allclose([point[d][e][wt][m] for m in ('MAE','RMSE','R2')],expected,rtol=1e-12,atol=1e-12)
  for draw in (0,17,1999):
   ix=np.concatenate([np.flatnonzero(inv==s) for s in sampled[draw]])
   for wt,w in (('observation',np.ones(len(ix))),('study_equal',1/ns[inv[ix]])):
    mt=metric(y[ix],pp[ix],w)
    np.testing.assert_allclose([draws[d][e][wt][m][draw] for m in ('MAE','RMSE','R2')],[mt[m] for m in ('MAE','RMSE','R2')],rtol=1e-12,atol=1e-12)
support={}
for task in ('A','B'):
 d=tasks[task];support[task]={}
 for e in models:
  p=point[d][e];checks={}
  for b in baselines[1:]:
   upper=float(np.quantile(draws[d][e]['study_equal']['RMSE']-draws[d][b]['study_equal']['RMSE'],1-.05/12))
   checks[f'adjusted_delta_RMSE_upper_vs_{b}']=upper<0
   checks[f'study_equal_MAE_no_worse_vs_{b}']=p['study_equal']['MAE']<=point[d][b]['study_equal']['MAE']
   for m in ('MAE','RMSE'):checks[f'observation_{m}_no_worse_vs_{b}']=p['observation'][m]<=point[d][b]['observation'][m]
  checks['study_equal_R2_positive']=p['study_equal']['R2']>0
  for wt in ('observation','study_equal'):checks[f'{wt}_RMSE_no_worse_vs_zero']=p[wt]['RMSE']<=point[d]['zero_response'][wt]['RMSE']
  support[task][e]={'supported':all(checks.values()),'checks':checks}

out={'phase1_predictions.tsv':prediction_rows,'phase1_metrics.tsv':metric_rows,'phase1_bootstrap.tsv':bootstrap_rows}
(ev/'delivery_tables.json').write_text(json.dumps(out,ensure_ascii=True,allow_nan=False),encoding='utf8')
(ev/'classification_evidence.json').write_text(json.dumps(support,indent=2),encoding='utf8')
(ev/'evaluation_summary.json').write_text(json.dumps({'protocol_sha256':ph,'metrics':point,'support':support,'fold_provenance':fold_provenance,'prediction_rows':len(prediction_rows),'metric_rows':len(metric_rows),'bootstrap_summary_rows':len(bootstrap_rows),'bootstrap_archive_sha256':hashlib.sha256((run/'bootstrap_draws.npz').read_bytes()).hexdigest(),'verification':{'all_folds_training_only_provenance':True,'all_baselines_recomputed_training_only':True,'independent_sklearn_point_metrics':True,'expanded_cluster_draw_checks':[0,17,1999],'no_input_changes':True}},ensure_ascii=True,indent=2,allow_nan=False),encoding='utf8')
print(json.dumps({'prediction_rows':len(prediction_rows),'metric_rows':len(metric_rows),'bootstrap_rows':len(bootstrap_rows),'task_supported':{t:any(v['supported'] for v in support[t].values()) for t in support},'overall_metrics':point},ensure_ascii=True),flush=True)
