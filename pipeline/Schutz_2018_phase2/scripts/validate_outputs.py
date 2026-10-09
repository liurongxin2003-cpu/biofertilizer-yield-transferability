from pathlib import Path
import os,sys,csv,json,hashlib,math
root=Path(__file__).resolve().parents[1];base=root.parent;ev=root/'evidence';sys.path.insert(0,str(base/'Schutz_2018_phase1'/'deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import numpy as np
import pandas as pd
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
read=lambda p:list(csv.DictReader(p.open(encoding='utf8',newline=''),delimiter='\t'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p=json.loads((ev/'protocol.json').read_text(encoding='utf8'));es=json.loads((ev/'evaluation_summary.json').read_text(encoding='utf8'));v=json.loads((ev/'frozen_variants.json').read_text(encoding='utf8'))
expected=json.loads((ev/'delivery_tables.json').read_text(encoding='utf8'));out={};checks=[];hashes={}
for name,rows in expected.items():
 delivered=read(root/name);assert len(rows)==len(delivered)
 for a,b in zip(rows,delivered):
  assert list(a)==list(b)
  for k,val in a.items():
   if isinstance(val,(int,float)):assert float(b[k])==float(val),(name,k,val,b[k])
   else:assert b[k]==val
 out[name]=delivered;hashes[name]=sha(root/name)
 checks.append('complete_value_roundtrip:'+name)
for name,h in p['inputs_sha256'].items():assert sha(base/name)==h
assert sha(ev/'frozen_variants.json')==p['variant_file_sha256']
checks.append('all_Phase0_Phase1_and_original_source_hashes_unchanged')
npred=0
for name,x in v.items():
 assert len(x['rows'])==x['cohort_n'] and len(set(r['study_id'] for r in x['rows']))==x['studies_n']
 for r in x['rows']:
  yy=float(r['logRR']);a=float(r['treatment_dry_yield_kg_ha']);b=float(r['control_dry_yield_kg_ha']);assert math.isfinite(yy) and a>0 and b>0 and yy==math.log(a/b)
  assert (yy>0)==(a>b) and (yy<0)==(a<b)
 z=np.load(root/'run'/name/'bootstrap_draws.npz');S=x['studies_n'];ss=np.random.default_rng(20261006).integers(0,S,size=(2000,S))
 assert np.array_equal(z['sampled_study_index'],ss) and np.array_equal(z['study_multiplicity'],np.stack([np.bincount(r,minlength=S) for r in ss]))
 assert sha(root/'run'/name/'bootstrap_draws.npz')==es['variants'][name]['bootstrap_archive_sha256']
 assert sha(root/'run'/name/'oof_predictions.npz')==es['variants'][name]['predictions_archive_sha256']
 preds=np.load(root/'run'/name/'oof_predictions.npz');assert np.array_equal(preds['observation_ids'],np.array([r['observation_id'] for r in x['rows']]))
 for key in preds.files:
  if '__' in key:assert np.isfinite(preds[key]).all();npred+=len(preds[key])
 for r in out['phase2_sensitivity_metrics.tsv']:
  if r['scenario']!=name:continue
  arr=z[f"{r['task']}__{r['estimator']}__{r['weighting']}__{r['metric']}"]
  if r['kind'].startswith('paired'):arr=arr-z[f"{r['task']}__{r['baseline']}__{r['weighting']}__{r['metric']}"]
  np.testing.assert_allclose([float(r['ci95_lower']),float(r['ci95_upper'])],np.quantile(arr,[.025,.975]),rtol=0,atol=1e-15)
  if r['adjusted_one_sided_upper']!='NA':assert float(r['adjusted_one_sided_upper'])==float(np.quantile(arr,1-.05/12))
assert npred==52320
checks.extend(['all52320_sensitivity_predictions_finite_complete','same_logRR_definition_all1758_rows_valid_no_changes','2000_seeded_cluster_draws_and_all_sensitivity_CIs_reproduced'])
reference=read(base/'Schutz_2018_phase1'/'phase1_bootstrap.tsv');refs=[r for r in out['phase2_sensitivity_metrics.tsv'] if r['scenario']=='phase1_main_1723_12_reference'];assert len(refs)==len(reference)==270
for a,b in zip(reference,refs):
 for k in ('kind','task','estimator','baseline','weighting','metric','point','ci95_lower','ci95_upper','adjusted_one_sided_upper','bootstrap_n','seed'):
  if k in ('point','ci95_lower','ci95_upper','adjusted_one_sided_upper') and a[k]!='NA':assert float(a[k])==float(b[k])
  else:assert a[k]==b[k]
checks.append('all270_Phase1_reference_rows_preserved')
main=pd.DataFrame(v['main_1723_crop_blind_11']['rows']).set_index('observation_id',drop=False);main.logRR=main.logRR.astype(float)
preds=read(base/'Schutz_2018_phase1'/'phase1_predictions.tsv');estimators=list(p['models'])+p['baselines'];yp={t:{e:pd.Series(index=main.index,dtype=float) for e in estimators} for t in ('A','B','C')}
for r in preds:
 assert float(r['logRR'])==main.loc[r['observation_id'],'logRR'];yp[r['task']][r['estimator']].loc[r['observation_id']]=float(r['prediction'])
cb=np.load(root/'run'/'main_1723_crop_blind_11'/'oof_predictions.npz')
for t in ('A','B'):
 for b in p['baselines']:assert np.array_equal(cb[t+'__'+b],yp[t][b].to_numpy())
assert all('crop' not in r['predictors'] and 'crop' not in r['categorical_predictors'] for n,r in v.items() if n.endswith('_11'))
checks.append('cropblind_exact_same_baselines_and_original_primary_members_no_crop_input')
size=main.groupby('study_id').size();cstat=main.groupby('crop').agg(n=('observation_id','size'),studies=('study_id','nunique'))
cbands={c:'singleton' if r.studies==1 else 'high_support' if r.studies>=3 and r.n>=20 else 'low_support' for c,r in cstat.iterrows()}
cache={}
for t in yp:
 for e,pp in yp[t].items():
  er=main.logRR-pp;tmp=pd.DataFrame({'ae':er.abs(),'se':er**2,'study':main.study_id});st=tmp.groupby('study').mean(numeric_only=True);st['RMSE']=np.sqrt(st.se)
  for m,col in (('MAE','ae'),('RMSE','RMSE')):
   values=st[col].to_numpy();cache[(t,e,m)]={'mean_of_study_metric':float(np.mean(values)),'std_of_study_metric':float(np.std(values,ddof=1)),**{'quantile_'+str(q):float(np.quantile(values,q)) for q in p['failure_audit']['study_error_quantiles']}}
  for sid,r in st.iterrows():cache[(t,e,'study',sid)]={'MAE':float(r.ae),'RMSE':float(r.RMSE)}
  for grouping,gg in (('study_size',{'large_ge20':main.study_id.map(size)>=20,'small_lt20':main.study_id.map(size)<20}),('crop_support',{g:main.crop.map(cbands)==g for g in ('high_support','low_support','singleton')})):
   for g,mask in gg.items():
    mm=main.loc[mask];pc=pp.loc[mask];ss=mm.groupby('study_id').size()
    for wt,w in (('observation',np.ones(len(mm))),('study_equal',1/mm.study_id.map(ss).to_numpy())):
     cache[(t,e,grouping,g,wt)]={'MAE':float(mean_absolute_error(mm.logRR,pc,sample_weight=w)),'RMSE':float(np.sqrt(mean_squared_error(mm.logRR,pc,sample_weight=w))),'R2':float(r2_score(mm.logRR,pc,sample_weight=w)),'n_observations':len(mm),'n_studies':mm.study_id.nunique(),'n_crops':mm.crop.nunique()}
for r in out['phase2_failure_audit.tsv']:
 assert r['source']=='unchanged_Phase1_main12_OOF' and r['descriptive_only']=='true'
 t,e,m=r['task'],r['estimator'],r['metric']
 if r['scope']=='study_error_distribution':expected_value=cache[(t,e,m)][r['statistic']]
 elif r['scope']=='individual_study':expected_value=cache[(t,e,'study',r['group_id'])][m]
 else:
  q=cache[(t,e,r['grouping'],r['group_id'],r['weighting'])];expected_value=q[m]
  for k in ('n_observations','n_studies','n_crops'):assert int(r[k])==q[k]
 np.testing.assert_allclose(float(r['value']),expected_value,rtol=1e-12,atol=1e-12)
checks.append('all6984_failure_descriptions_independently_recomputed_from_Phase1_only')
fit=json.loads((ev/'fit_summary.json').read_text());assert fit['model_fits']==1278 and fit['outer_folds']==426 and fit['tuning_fits']==fit['SHAP']==0
checks.append('426_fold_training_provenance_and1278_fixed_fits_no_tuning_no_SHAP')
qc={'passed':True,'failed':[],'checks':checks,'output_sha256':hashes,'row_counts':{k:len(r) for k,r in out.items()},'formal_decision_written':False}
(ev/'output_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf8');print(json.dumps(qc))
print(json.dumps([r for r in out['phase2_failure_audit.tsv'] if r['scope']=='metadata_stratum' and r['task'] in ('A','B') and r['weighting']=='study_equal' and r['metric']=='RMSE'],ensure_ascii=True))
