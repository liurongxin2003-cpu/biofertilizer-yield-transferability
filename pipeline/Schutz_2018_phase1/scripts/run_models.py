from pathlib import Path
import os,sys
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'deps'))
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import csv,json,hashlib,time,concurrent.futures,warnings,platform
import numpy as np
import pandas as pd
import scipy,sklearn,xgboost
from sklearn.preprocessing import OneHotEncoder,StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from threadpoolctl import threadpool_limits

ev=root/'evidence';run=root/'run';old=root.parent/'Schutz_2018_phase0'
protocol=json.loads((ev/'protocol.json').read_text(encoding='utf8'))
protocol_sha=hashlib.sha256((ev/'protocol.json').read_bytes()).hexdigest()
assert protocol_sha==(ev/'protocol_sha256.txt').read_text(encoding='utf8')
for name,sha in protocol['input_sha256'].items():assert hashlib.sha256((old/name).read_bytes()).hexdigest()==sha
records=list(csv.DictReader((old/'analysis_matrix.tsv').open(encoding='utf8',newline=''),delimiter='\t'))
splits=list(csv.DictReader((old/'split_manifest.tsv').open(encoding='utf8',newline=''),delimiter='\t'))
features=protocol['predictors'];nums=protocol['numeric_predictors'];cats=protocol['categorical_predictors']
df=pd.DataFrame(records).set_index('observation_id',drop=False)
assert len(df)==1723 and df.index.is_unique and df.study_id.nunique()==169
assert set(features).isdisjoint({'logRR','study_id','pair_id','treatment_dry_yield_kg_ha','control_dry_yield_kg_ha','treatment_excel_row','control_excel_row'})
# Parsing, not fitting. No global imputation, vocabularies, scaling or feature screening.
for n in nums:df[n]=df[n].replace('NA',np.nan).astype(float)
for c in cats:df[c]=df[c].replace('NA','__MISSING__').astype(str)
df['logRR']=df['logRR'].astype(float)
assert np.isfinite(df.logRR).all()
def unpack(s):return [] if s=='NA' else s.split(',')
def digest(ids):return hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
for f in splits:
 tr,te,pu=[unpack(f[r+'_ids']) for r in ('train','test','purge')]
 assert not set(df.loc[tr,'study_id'])&set(df.loc[te,'study_id'])
 assert not set(df.loc[tr,'study_id'])&set(df.loc[pu,'study_id'])
 assert set(tr)|set(te)|set(pu)==set(df.index) and len(tr)+len(te)+len(pu)==1723
 for role,ids in zip(('train','test','purge'),(tr,te,pu)):assert digest(ids)==f[role+'_ids_sha256']
versions={'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'scikit-learn':sklearn.__version__,'xgboost':xgboost.__version__}
for k in versions:
 if k!='python':assert versions[k]==protocol['libraries'][k],(k,versions[k])
(ev/'runtime_versions.json').write_text(json.dumps(versions,indent=2),encoding='utf8')
def fit_fold(item):
 i,f=item;dest=run/f'fold_{i:04d}.json'
 if dest.exists():
  saved=json.loads(dest.read_text(encoding='utf8'));assert saved['protocol_sha256']==protocol_sha and saved['fold_id']==f['fold_id'];return saved
 start=time.monotonic();tr=unpack(f['train_ids']);te=unpack(f['test_ids']);pu=unpack(f['purge_ids'])
 train=df.loc[tr,features];test=df.loc[te,features]
 y=df.loc[tr,'logRR'].to_numpy();groups=df.loc[tr,'study_id'].to_numpy()
 _,inv,counts=np.unique(groups,return_inverse=True,return_counts=True)
 weights=len(tr)/(len(counts)*counts[inv])
 assert np.isclose(weights.mean(),1) and np.allclose(np.bincount(inv,weights=weights),len(tr)/len(counts))
 a=train[nums].to_numpy(dtype=float);b=test[nums].to_numpy(dtype=float)
 med=np.nanmedian(a,axis=0);assert np.isfinite(med).all()
 am=np.isnan(a).astype(float);bm=np.isnan(b).astype(float)
 ai=np.where(np.isnan(a),med,a);bi=np.where(np.isnan(b),med,b)
 scaler=StandardScaler().fit(ai,sample_weight=weights)
 ac=train[cats].to_numpy(dtype=object);bc=test[cats].to_numpy(dtype=object)
 levels=[sorted(set(ac[:,j])|{'__MISSING__'}) for j in range(len(cats))]
 enc=OneHotEncoder(categories=levels,handle_unknown='ignore',sparse_output=False,dtype=np.float64).fit(ac)
 ao=enc.transform(ac);bo=enc.transform(bc)
 au=np.column_stack([~np.isin(ac[:,j],lv) for j,lv in enumerate(levels)]).astype(float)
 bu=np.column_stack([~np.isin(bc[:,j],lv) for j,lv in enumerate(levels)]).astype(float)
 assert au.sum()==0
 ar=np.column_stack([ai,am,ao,au]);br=np.column_stack([bi,bm,bo,bu])
 az=np.column_stack([scaler.transform(ai),am,ao,au]);bz=np.column_stack([scaler.transform(bi),bm,bo,bu])
 assert np.isfinite(ar).all() and np.isfinite(br).all() and np.isfinite(az).all() and np.isfinite(bz).all()
 seed=protocol['base_seed']+i
 estimators={'Ridge':Ridge(**protocol['models']['Ridge']),'RandomForest':RandomForestRegressor(**protocol['models']['RandomForest'],random_state=seed),'XGBoost':XGBRegressor(**protocol['models']['XGBoost'],random_state=seed)}
 predictions={}
 for name,model in estimators.items():
  xx,zz=(az,bz) if name=='Ridge' else (ar,br)
  model.fit(xx,y,sample_weight=weights)
  pred=model.predict(zz).astype(float);assert np.isfinite(pred).all()
  predictions[name]=pred.tolist()
 means=np.bincount(inv,weights=y)/counts
 baseline={'zero_response':0.0,'training_observation_mean':float(np.mean(y)),'training_study_equal_mean':float(np.mean(means))}
 for name,mean in baseline.items():predictions[name]=[mean]*len(te)
 out={'fold_id':f['fold_id'],'fold_index':i,'design':f['design'],'held_out_level':f['held_out_level'],'primary_or_secondary':f['primary_or_secondary'],'test_ids':te,'protocol_sha256':protocol_sha,'train_ids_sha256':digest(tr),'test_ids_sha256':digest(te),'purge_ids_sha256':digest(pu),'train_n':len(tr),'test_n':len(te),'purge_n':len(pu),'train_studies_n':len(counts),'study_overlap':0,'model_seed':seed,'model_fits':3,'transformed_columns':ar.shape[1],'preprocessor_fit_ids_sha256':digest(tr),'numeric_train_medians':med.tolist(),'numeric_train_scaler_mean':scaler.mean_.tolist(),'numeric_train_scaler_scale':scaler.scale_.tolist(),'categorical_training_vocabulary':dict(zip(cats,levels)),'test_unknown_counts':dict(zip(cats,bu.sum(axis=0).astype(int).tolist())),'baseline_constants':baseline,'predictions':predictions,'elapsed_seconds':time.monotonic()-start}
 dest.write_text(json.dumps(out,ensure_ascii=True,allow_nan=False),encoding='utf8')
 return out

if __name__=='__main__':
 start=time.monotonic();results=[]
 with threadpool_limits(limits=1),concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
  futures=[executor.submit(fit_fold,item) for item in enumerate(splits)]
  for future in concurrent.futures.as_completed(futures):
   result=future.result();results.append(result)
   if len(results)%10==0 or len(results)==215:print(json.dumps({'completed_folds':len(results),'total_folds':215,'latest_fold':result['fold_id'],'elapsed_seconds':round(time.monotonic()-start,1)},ensure_ascii=True),flush=True)
 assert len(results)==215 and sum(r['model_fits'] for r in results)==645
 summary={'protocol_sha256':protocol_sha,'model_fits':645,'outer_folds':215,'preprocessing_fit_scope':'each outer train only','tuning_fits':0,'SHAP':0,'elapsed_seconds':time.monotonic()-start,'runtime':versions,'all_input_hashes_verified':True}
 (ev/'fit_summary.json').write_text(json.dumps(summary,ensure_ascii=True,indent=2),encoding='utf8')
 print(json.dumps(summary,ensure_ascii=True),flush=True)
