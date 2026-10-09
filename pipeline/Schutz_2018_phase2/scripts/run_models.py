from pathlib import Path
import os,sys
root=Path(__file__).resolve().parents[1];base=root.parent
sys.path.insert(0,str(base/'Schutz_2018_phase1'/'deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,hashlib,time,concurrent.futures,platform
import numpy as np
import pandas as pd
import scipy,sklearn,xgboost
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from threadpoolctl import threadpool_limits
ev=root/'evidence';run=root/'run'
protocol=json.loads((ev/'protocol.json').read_text(encoding='utf8'));ph=hashlib.sha256((ev/'protocol.json').read_bytes()).hexdigest()
assert ph==(ev/'protocol_sha256.txt').read_text(encoding='utf8')
for name,h in protocol['inputs_sha256'].items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==h
assert hashlib.sha256((ev/'frozen_variants.json').read_bytes()).hexdigest()==protocol['variant_file_sha256']
variants=json.loads((ev/'frozen_variants.json').read_text(encoding='utf8'))
old=json.loads((base/'Schutz_2018_phase1'/'evidence'/'protocol.json').read_text(encoding='utf8'))
assert protocol['models']==old['models'] and protocol['libraries']==old['libraries']
versions={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scikit-learn':sklearn.__version__,'xgboost':xgboost.__version__}
for k,v in versions.items():
 if k!='python':assert v==protocol['libraries'][k]
(ev/'runtime_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf8')
digest=lambda ids:hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
frames={}
for name,v in variants.items():
 df=pd.DataFrame(v['rows']).set_index('observation_id',drop=False)
 assert len(df)==v['cohort_n'] and df.study_id.nunique()==v['studies_n'] and df.index.is_unique
 for n in v['numeric_predictors']:df[n]=df[n].replace('NA',np.nan).astype(float)
 for c in v['categorical_predictors']:df[c]=df[c].replace('NA','__MISSING__').astype(str)
 df.logRR=df.logRR.astype(float);assert np.isfinite(df.logRR).all()
 assert set(v['predictors']).isdisjoint({'logRR','study_id','pair_id','treatment_dry_yield_kg_ha','control_dry_yield_kg_ha','treatment_excel_row','control_excel_row'})
 for f in v['splits']:
  tr,te,pu=[f[r+'_ids'] for r in ('train','test','purge')]
  assert not set(df.loc[tr,'study_id'])&set(df.loc[te+pu,'study_id'])
  assert len(set(tr+te+pu))==len(tr+te+pu)==len(df)
  for r,ids in zip(('train','test','purge'),(tr,te,pu)):assert digest(ids)==f[r+'_ids_sha256']
 frames[name]=df;(run/name).mkdir(exist_ok=True)
def fit_fold(item):
 name,i,f=item;v=variants[name];df=frames[name];nums=v['numeric_predictors'];cats=v['categorical_predictors'];features=v['predictors']
 dest=run/name/f'fold_{i:04d}.json'
 if dest.exists():
  r=json.loads(dest.read_text(encoding='utf8'));assert r['protocol_sha256']==ph and r['fold_id']==f['fold_id'];return r
 start=time.monotonic();tr=f['train_ids'];te=f['test_ids'];pu=f['purge_ids']
 train=df.loc[tr,features];test=df.loc[te,features];y=df.loc[tr,'logRR'].to_numpy();groups=df.loc[tr,'study_id'].to_numpy()
 _,inv,counts=np.unique(groups,return_inverse=True,return_counts=True)
 weights=len(tr)/(len(counts)*counts[inv]);assert np.isclose(weights.mean(),1)
 assert np.allclose(np.bincount(inv,weights=weights),len(tr)/len(counts))
 a=train[nums].to_numpy(float);b=test[nums].to_numpy(float);med=np.nanmedian(a,axis=0);assert np.isfinite(med).all()
 am=np.isnan(a).astype(float);bm=np.isnan(b).astype(float);ai=np.where(np.isnan(a),med,a);bi=np.where(np.isnan(b),med,b)
 scaler=StandardScaler().fit(ai,sample_weight=weights)
 ac=train[cats].to_numpy(object);bc=test[cats].to_numpy(object)
 levels=[sorted(set(ac[:,j])|{'__MISSING__'}) for j in range(len(cats))]
 enc=OneHotEncoder(categories=levels,handle_unknown='ignore',sparse_output=False,dtype=np.float64).fit(ac)
 ao=enc.transform(ac);bo=enc.transform(bc)
 au=np.column_stack([~np.isin(ac[:,j],lv) for j,lv in enumerate(levels)]).astype(float);bu=np.column_stack([~np.isin(bc[:,j],lv) for j,lv in enumerate(levels)]).astype(float)
 assert au.sum()==0
 ar=np.column_stack([ai,am,ao,au]);br=np.column_stack([bi,bm,bo,bu]);az=np.column_stack([scaler.transform(ai),am,ao,au]);bz=np.column_stack([scaler.transform(bi),bm,bo,bu])
 assert all(np.isfinite(x).all() for x in (ar,br,az,bz))
 seed=20261006+f['seed_index']
 estimators={'Ridge':Ridge(**protocol['models']['Ridge']),'RandomForest':RandomForestRegressor(**protocol['models']['RandomForest'],random_state=seed),'XGBoost':XGBRegressor(**protocol['models']['XGBoost'],random_state=seed)}
 predictions={}
 for model,est in estimators.items():
  xx,zz=(az,bz) if model=='Ridge' else (ar,br);est.fit(xx,y,sample_weight=weights)
  yp=est.predict(zz).astype(float);assert np.isfinite(yp).all();predictions[model]=yp.tolist()
 means=np.bincount(inv,weights=y)/counts
 baselines={'zero_response':0.,'training_observation_mean':float(np.mean(y)),'training_study_equal_mean':float(np.mean(means))}
 for model,vv in baselines.items():predictions[model]=[vv]*len(te)
 out={'variant':name,'fold_id':f['fold_id'],'fold_index':i,'seed_index':f['seed_index'],'task':f['task'],'test_ids':te,'protocol_sha256':ph,'predictors':features,'train_ids_sha256':digest(tr),'test_ids_sha256':digest(te),'purge_ids_sha256':digest(pu),'preprocessor_fit_ids_sha256':digest(tr),'train_n':len(tr),'test_n':len(te),'purge_n':len(pu),'train_studies_n':len(counts),'study_overlap':0,'model_seed':seed,'model_fits':3,'numeric_train_medians':med.tolist(),'numeric_train_scaler_mean':scaler.mean_.tolist(),'numeric_train_scaler_scale':scaler.scale_.tolist(),'categorical_training_vocabulary':dict(zip(cats,levels)),'test_unknown_counts':dict(zip(cats,bu.sum(axis=0).astype(int).tolist())),'transformed_columns':ar.shape[1],'baseline_constants':baselines,'predictions':predictions,'elapsed_seconds':time.monotonic()-start}
 dest.write_text(json.dumps(out,ensure_ascii=True,allow_nan=False),encoding='utf8');return out
if __name__=='__main__':
 jobs=[(name,i,f) for name,v in variants.items() for i,f in enumerate(v['splits'])];assert len(jobs)==426
 start=time.monotonic();results=[]
 with threadpool_limits(limits=1),concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
  futures=[executor.submit(fit_fold,j) for j in jobs]
  for fut in concurrent.futures.as_completed(futures):
   r=fut.result();results.append(r)
   if len(results)%25==0 or len(results)==len(jobs):print(json.dumps({'completed_folds':len(results),'total_folds':len(jobs),'latest_variant':r['variant'],'latest_fold':r['fold_id'],'elapsed_seconds':round(time.monotonic()-start,1)},ensure_ascii=True),flush=True)
 assert len(results)==426 and sum(r['model_fits'] for r in results)==1278
 summary={'protocol_sha256':ph,'model_fits':1278,'outer_folds':426,'folds_by_variant':{v:len(x['splits']) for v,x in variants.items()},'tuning_fits':0,'SHAP':0,'all_inputs_unchanged':True,'runtime':versions,'elapsed_seconds':time.monotonic()-start}
 (ev/'fit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8');print(json.dumps(summary),flush=True)
