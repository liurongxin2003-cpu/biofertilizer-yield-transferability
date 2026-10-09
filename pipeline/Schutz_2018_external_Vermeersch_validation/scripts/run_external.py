from pathlib import Path
import os, sys, csv, json, hashlib, platform, time
from datetime import datetime, timezone
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT.parent
sys.path.insert(0, str(BASE / 'Schutz_2018_phase1/deps'))
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
import scipy, sklearn, xgboost, joblib, threadpoolctl
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
from threadpoolctl import threadpool_limits

EV = ROOT / 'evidence'
RUN = ROOT / 'run'
EV.mkdir(parents=True, exist_ok=True)
RUN.mkdir(exist_ok=True)
def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding='utf8')
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def digest(ids):
    return hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()

assert not (RUN/'fit_started.json').exists(), 'Refuse automatic refitting: a fit has already started.'
old = BASE/'Schutz_2018_phase0'
prpath = BASE/'Schutz_2018_phase1/evidence/protocol.json'
pr = json.loads(prpath.read_text(encoding='utf8'))
assert sha(prpath) == '1420e0b80930256655aa60e145bc8b8ae93823a080ffa4c4a40fe0e64dfcbb0d'
sources = {str(prpath.relative_to(BASE)): sha(prpath)}
for name, expected in pr['input_sha256'].items():
    assert sha(old/name) == expected
    sources[str((old/name).relative_to(BASE))] = expected
audroot = BASE/'Schutz_2018_external_Sato_Vermeersch'
for name in ('structural_audit.json','vermeersch_pairing_audit.json','source_hashes.json','EXTERNAL_SCREEN_REPORT.md'):
    sources[str((audroot/name).relative_to(BASE))] = sha(audroot/name)
versions = {'python':platform.python_version(), 'numpy':np.__version__, 'pandas':pd.__version__, 'scipy':scipy.__version__, 'scikit-learn':sklearn.__version__, 'xgboost':xgboost.__version__, 'joblib':joblib.__version__, 'threadpoolctl':threadpoolctl.__version__}
assert versions['python']=='3.12.14'
for key, version in pr['libraries'].items():
    assert versions[key] == version, (key, versions[key], version)

# Freeze the full-cohort extension and descriptive interpretation before fitting/evaluation.
ep = {'frozen_at_utc':datetime.now(timezone.utc).isoformat(), 'phase1_protocol_sha256':sha(prpath), 'source_sha256':sources,
      'predictors':pr['predictors'], 'models':pr['models'], 'training_n':1723, 'training_studies':169,
      'training_weights':pr['training_weights'], 'preprocessing':pr['preprocessing'],
      'seed':pr['base_seed'], 'seed_extension':'Full-cohort fit uses base_seed=20261006; no held-out fold index and no seed search.',
      'fits_per_model':1, 'external_outcome_in_fit_or_tuning':False, 'SHAP':False,
      'outcome':'ln(treatment mean yield at 15% moisture / paired control mean at 15% moisture); common dry factor 0.85 cancels',
      'missing_external_predictors':['soil_organic_matter_pct','soil_ph_water','available_p_kg_ha','n_applied_kg_ha','p_applied_kg_ha'],
      'baselines':pr['baselines'], 'metrics':['MAE','RMSE','prediction_mean','observed_mean','mean_bias'],
      'weighting':['comparison equal: 1/57','field-year equal: 1/(5*n_in_field_year)'],
      'dependence':'13 controls nested in 5 field-years; no iid comparison CI or significance claim. Merelbeke recurs across years.',
      'decision_rule':{
        'scope':'Descriptive external comparison, no statistical significance or universal transferability claim; report all models without selecting a deployable winner.',
        'model_concordant_gain':'Strictly lower MAE and RMSE than BOTH frozen training-mean baselines under BOTH weighting schemes; RMSE no worse than zero under both schemes.',
        'EXTERNAL_SUGGESTS_TRANSFERABILITY':'All three prespecified models show concordant gain.',
        'EXTERNAL_SUPPORTS_LIMITED_GENERALIZABILITY':'No model shows concordant gain AND each model has field-year-equal RMSE no better than at least one frozen training-mean baseline.',
        'EXTERNAL_INCONCLUSIVE':'Other, discordant patterns or failed structural checks.',
        'numeric_tolerance':1e-12}}
save(EV/'external_protocol.json',ep)
save(EV/'runtime_versions.json',versions)

records = list(csv.DictReader((old/'analysis_matrix.tsv').open(encoding='utf8', newline=''), delimiter='\t'))
df = pd.DataFrame(records).set_index('observation_id',drop=False)
assert len(df)==1723 and df.index.is_unique and df.study_id.nunique()==169
features, nums, cats = pr['predictors'], pr['numeric_predictors'], pr['categorical_predictors']
for n in nums: df[n]=df[n].replace('NA',np.nan).astype(float)
for c in cats: df[c]=df[c].replace('NA','__MISSING__').astype(str)
y = df.logRR.astype(float).to_numpy()
assert np.isfinite(y).all()

# Only the frozen predictors and provenance enter the external input construction.
pairs = json.loads((audroot/'vermeersch_pairing_audit.json').read_text(encoding='utf8'))
structure = json.loads((audroot/'structural_audit.json').read_text(encoding='utf8'))['vermeersch']
assert len(pairs)==57==structure['comparisons']
assert len({p['experiment_id'] for p in pairs})==5
assert len({p['control_cluster_id'] for p in pairs})==13
assert len({(p['experiment_id'],p['cultivar'],p['strain']) for p in pairs})==57
external = pd.DataFrame([p['frozen_predictors'] for p in pairs],columns=features)
for p in pairs:
    assert set(p['frozen_predictors']) == set(features)
    assert sum(v is not None for v in p['frozen_predictors'].values()) == 7
    assert {k for k,v in p['frozen_predictors'].items() if v is None} == set(ep['missing_external_predictors'])
    for k, v in structure['strict_predictors'].items(): assert p['frozen_predictors'][k] == v
for n in nums: external[n]=external[n].astype(float)
# Use original audited categorical scalar codes, not pandas float string coercion.
for c in cats: external[c]=[str(p['frozen_predictors'][c]) if p['frozen_predictors'][c] is not None else '__MISSING__' for p in pairs]
assert external.fillna('NA').drop_duplicates().shape[0]==1
train=df[features]
_,inv,counts=np.unique(df.study_id.to_numpy(),return_inverse=True,return_counts=True)
weights=len(df)/(len(counts)*counts[inv])
assert np.isclose(weights.mean(),1) and np.allclose(np.bincount(inv,weights=weights),len(df)/len(counts))
a=train[nums].to_numpy(dtype=float); b=external[nums].to_numpy(dtype=float)
med=np.nanmedian(a,axis=0); assert np.isfinite(med).all()
am=np.isnan(a).astype(float); bm=np.isnan(b).astype(float)
ai=np.where(np.isnan(a),med,a); bi=np.where(np.isnan(b),med,b)
scaler=StandardScaler().fit(ai,sample_weight=weights)
ac=train[cats].to_numpy(dtype=object); bc=external[cats].to_numpy(dtype=object)
levels=[sorted(set(ac[:,j])|{'__MISSING__'}) for j in range(len(cats))]
enc=OneHotEncoder(categories=levels,handle_unknown='ignore',sparse_output=False,dtype=np.float64).fit(ac)
ao=enc.transform(ac); bo=enc.transform(bc)
au=np.column_stack([~np.isin(ac[:,j],lv) for j,lv in enumerate(levels)]).astype(float)
bu=np.column_stack([~np.isin(bc[:,j],lv) for j,lv in enumerate(levels)]).astype(float)
assert au.sum()==0
ar=np.column_stack([ai,am,ao,au]); br=np.column_stack([bi,bm,bo,bu])
az=np.column_stack([scaler.transform(ai),am,ao,au]); bz=np.column_stack([scaler.transform(bi),bm,bo,bu])
assert all(np.isfinite(m).all() for m in (ar,br,az,bz))
basevalues={'zero_response':0.0,'training_observation_mean':float(np.mean(y)), 'training_study_equal_mean':float(np.mean(np.bincount(inv,weights=y)/counts))}
prep={'numeric_train_medians':med.tolist(), 'numeric_train_scaler_mean':scaler.mean_.tolist(), 'numeric_train_scaler_scale':scaler.scale_.tolist(),
      'numeric_order':nums, 'categorical_order':cats, 'categorical_training_vocabulary':dict(zip(cats,levels)),
      'external_unknown_counts':dict(zip(cats,bu.sum(axis=0).astype(int).tolist())), 'transformed_columns':ar.shape[1],
      'preprocessor_fit_ids_sha256':digest(df.index.tolist()), 'training_study_ids':sorted(df.study_id.unique().tolist()),
      'baseline_constants':basevalues}
save(EV/'preprocessing_audit.json',prep)
joblib.dump({'medians':med,'scaler':scaler,'encoder':enc,'levels':levels,'numeric_order':nums,'categorical_order':cats},RUN/'preprocessor.joblib')
save(RUN/'fit_started.json',{'utc':datetime.now(timezone.utc).isoformat(), 'protocol_sha256':sha(EV/'external_protocol.json')})
models={'Ridge':Ridge(**pr['models']['Ridge']), 'RandomForest':RandomForestRegressor(**pr['models']['RandomForest'],random_state=pr['base_seed']), 'XGBoost':XGBRegressor(**pr['models']['XGBoost'],random_state=pr['base_seed'])}
preds={}; fitlog=[]
with threadpool_limits(limits=1):
    for name,model in models.items():
        start=time.monotonic()
        xx,zz=(az,bz) if name=='Ridge' else (ar,br)
        model.fit(xx,y,sample_weight=weights)
        pred=model.predict(zz).astype(float)
        assert np.isfinite(pred).all() and np.ptp(pred)<1e-12
        preds[name]=pred.tolist()
        joblib.dump(model,RUN/(name+'.joblib'))
        fitlog.append({'model':name,'fits':1,'train_n':1723,'train_studies':169,'seed':pr['base_seed'],'seconds':time.monotonic()-start})
        save(EV/'fit_log.json',fitlog)
        print(json.dumps({'fitted':name,'fits':1,'seconds':fitlog[-1]['seconds']}),flush=True)
for name,constant in basevalues.items(): preds[name]=[constant]*57
# Predictions are written before computing/inspecting any external outcome.
save(RUN/'predictions_before_evaluation.json',{'external_ids':[f'VERM2026_{i+1:03d}' for i in range(57)],'predictions':preds,'protocol_sha256':sha(EV/'external_protocol.json')})

# Evaluation stage: source means are now used solely to calculate fixed logRR.
obs=[]
for i,p in enumerate(pairs):
    t=float(p['treatment_mean_kg_ha_at_15pct_moisture']); c=float(p['control_mean_kg_ha_at_15pct_moisture'])
    assert np.isfinite([t,c]).all() and t>0 and c>0
    target=float(np.log(t/c))
    assert np.isclose(target,np.log((t*.85)/(c*.85)),rtol=0,atol=1e-15)
    row={'comparison_id':f'VERM2026_{i+1:03d}', 'field_year_id':p['experiment_id'],'control_cluster_id':p['control_cluster_id'],
         'year':p['year'],'site':p['site'],'cultivar_provenance_only':p['cultivar'],'strain_provenance_only':p['strain'],
         'treatment_yield_kg_ha_15pct_moisture':t,'control_yield_kg_ha_15pct_moisture':c,
         'treatment_source_docx_row':p['treatment_source_docx_row'],'control_source_docx_row':p['control_source_docx_row'],'external_logRR':target}
    row.update({k:('NA' if p['frozen_predictors'][k] is None else p['frozen_predictors'][k]) for k in features})
    row.update({'prediction_'+k:values[i] for k,values in preds.items()})
    obs.append(row)
targets=np.array([r['external_logRR'] for r in obs])
g=np.array([r['field_year_id'] for r in obs]); cg=np.array([r['control_cluster_id'] for r in obs])
_,gi,gcount=np.unique(g,return_inverse=True,return_counts=True)
metricweights={'comparison':np.ones(57)/57,'field_year_equal':1/(5*gcount[gi])}
for cc in set(cg):
    members=[r for r in obs if r['control_cluster_id']==cc]
    assert len({r['control_yield_kg_ha_15pct_moisture'] for r in members})==1
    assert len({r['field_year_id'] for r in members})==1
assert np.allclose(np.bincount(gi,weights=metricweights['field_year_equal']),.2)
metrics=[]
def metricsrow(scope,weighting,field_id,idx,w,name):
    yy=targets[idx];pp=np.array(preds[name])[idx]; w=w/w.sum();error=pp-yy
    return {'scope':scope,'weighting':weighting,'field_year_id':field_id,'method':name,
            'n_comparisons':int(len(yy)),'n_control_clusters':len(set(cg[idx])),'n_field_years':len(set(g[idx])),
            'MAE':float(np.sum(w*np.abs(error))),'RMSE':float(np.sqrt(np.sum(w*error**2))),
            'prediction_mean':float(np.sum(w*pp)),'observed_mean':float(np.sum(w*yy)),'mean_bias':float(np.sum(w*error))}
for wt,w in metricweights.items():
    for name in preds: metrics.append(metricsrow('overall',wt,'ALL',np.arange(57),w,name))
for field in sorted(set(g)):
    idx=np.where(g==field)[0]
    for name in preds: metrics.append(metricsrow('field_year','within_field_comparison_equal',field,idx,np.ones(len(idx)),name))
lookup={(r['weighting'],r['method']):r for r in metrics if r['scope']=='overall'}
gain={}
for name in models:
    gain[name]=all(lookup[(wt,name)][m] < lookup[(wt,baseline)][m]-1e-12 for wt in metricweights for m in ('MAE','RMSE') for baseline in ('training_observation_mean','training_study_equal_mean')) and all(lookup[(wt,name)]['RMSE']<=lookup[(wt,'zero_response')]['RMSE']+1e-12 for wt in metricweights)
limited=not any(gain.values()) and all(any(lookup[('field_year_equal',name)]['RMSE'] >= lookup[('field_year_equal',b)]['RMSE']-1e-12 for b in ('training_observation_mean','training_study_equal_mean')) for name in models)
decision='EXTERNAL_SUGGESTS_TRANSFERABILITY' if all(gain.values()) else ('EXTERNAL_SUPPORTS_LIMITED_GENERALIZABILITY' if limited else 'EXTERNAL_INCONCLUSIVE')
save(EV/'delivery_tables.json',{'external_predictions.tsv':obs,'external_metrics.tsv':metrics})
save(EV/'evaluation_summary.json',{'candidate_conclusion':decision,'model_concordant_gain':gain,'field_year_counts':dict(zip(*[v.tolist() for v in np.unique(g,return_counts=True)])),
      'outcome_n':57,'outcome_missing':0,'outcome_nonfinite':0,'outcome_min':float(targets.min()),'outcome_max':float(targets.max()),
      'unique_external_predictor_profiles':1,'external_known_predictors':7,'external_missing_predictors':5,'model_fits':3,'fit_log':fitlog,
      'source_hashes_unchanged':all(sha(BASE/name)==expected for name,expected in sources.items()),'overall_metrics':[r for r in metrics if r['scope']=='overall']})
assert all(sha(BASE/name)==expected for name,expected in sources.items())
print(json.dumps({'metrics':[r for r in metrics if r['scope']=='overall'],'candidate_conclusion':decision},indent=2),flush=True)
