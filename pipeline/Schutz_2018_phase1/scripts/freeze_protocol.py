from pathlib import Path
import csv,json,hashlib,datetime
root=Path(__file__).resolve().parents[1];ev=root/'evidence';old=root.parent/'Schutz_2018_phase0'
expected={'analysis_matrix.tsv':'66fd6e7ddd65426225167bfa3016f518cad26ff836e0748decdd3afa39209094','split_manifest.tsv':'0c9c0881db2fc520af24aec352f2f9066801cd7cdd0e57ea61c3270b4e0c6583','data_dictionary.tsv':'f2e49b119be38ac84214f33cb6bc7abac0808402fb7fa9aa815f9b420fb72cfa'}
for name,digest in expected.items():assert hashlib.sha256((old/name).read_bytes()).hexdigest()==digest,name
assert not list((root/'run').glob('fold_*.json')),'Cannot freeze after predictions exist'
dictionary=list(csv.DictReader((old/'data_dictionary.tsv').open(encoding='utf8'),delimiter='\t'))
features=[r['column'] for r in dictionary if r['role']=='main_predictor']
categorical=[r['column'] for r in dictionary if r['role']=='main_predictor' and r['dtype'].startswith('category')]
numeric=[r['column'] for r in dictionary if r['role']=='main_predictor' and r['dtype']=='float']
assert len(features)==12 and len(categorical)==len(numeric)==6
protocol={
 'frozen_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'input_sha256':expected,'predictors':features,'categorical_predictors':categorical,'numeric_predictors':numeric,
 'target':'logRR','cohort_n':1723,'studies_n':169,
 'designs':{'A':'LOSO','B':'LOCO_whole_study_purge','C':'function_secondary_whole_study_purge'},
 'models':{
  'Ridge':{'alpha':10.0,'solver':'cholesky','fit_intercept':True},
  'RandomForest':{'n_estimators':300,'max_depth':6,'min_samples_leaf':8,'max_features':1.0,'bootstrap':True,'oob_score':False,'n_jobs':1},
  'XGBoost':{'n_estimators':400,'max_depth':3,'learning_rate':0.03,'min_child_weight':5.0,'reg_lambda':10.0,'reg_alpha':0.0,'subsample':1.0,'colsample_bytree':1.0,'tree_method':'hist','objective':'reg:squarederror','n_jobs':1,'verbosity':0}
 },
 'base_seed':20261006,'model_seed_policy':'base_seed+zero_based_frozen_split_manifest_index',
 'hyperparameter_search':False,'feature_selection':False,'SHAP':False,
 'training_weights':'N_train/(S_train*n_rows_in_training_study); mean row weight1; each study equal total weight',
 'preprocessing':'training-only numeric median; six fixed numeric missing indicators; Ridge-only training sample-weighted numeric StandardScaler; train-only categorical vocabulary with fixed MISSING token; OHE handle_unknown=ignore plus six fixed unknown indicators; all derive only from12 source predictors',
 'baselines':['zero_response','training_observation_mean','training_study_equal_mean'],
 'metrics':['MAE','RMSE','R2'],'metric_weights':['observation','study_equal'],
 'bootstrap':{'resamples':2000,'seed':20261006,'unit':'original study_id','resampled_clusters_per_draw':169,'paired_draws_shared_across_models_baselines_designs':True,'CI':'95% percentile; no refitting; conditional on fixed out-of-fold predictions','draw_storage':'run/bootstrap_draws.npz; study multiplicity preserves all observations/folds within each sampled study'},
 'decision_rule':{
  'primary_loss':'study_equal RMSE; negative model-minus-baseline delta means improvement',
  'multiplicity':'one-sided upper percentile bound at 1-0.05/12 for3models*2primarydesigns*2mean_baselines; secondary function not in conclusion gate',
  'model_design_supported':'upper adjusted paired deltaRMSE<0 vs both training mean baselines; point study_equal R2>0; point study_equal MAE and observation MAE/RMSE no worse than both mean baselines; point RMSE under both weights no worse than zero-response baseline',
  'GENERALIZABLE_SIGNAL':'A and B each have at least one supported prespecified model; can be different models, no deployable winner selected',
  'STUDY_GENERALIZES_CROP_TRANSFER_WEAK':'A has a supported model, B does not',
  'LIMITED_GENERALIZABLE_SIGNAL':'A does not have a supported model; any B-only support reported descriptively',
  'interpretation':'conservative prespecified conditional uncertainty classification, not formal causal inference or proof that any possible algorithm cannot generalize'
 },
 'out_of_scope':['grouped5fold repeat sensitivity','unweighted-training sensitivity','crop-blind11 sensitivity','external_validation','final full-cohort model','SHAP'],
 'libraries':{'numpy':'2.3.5','pandas':'3.0.1','scipy':'1.16.3','scikit-learn':'1.7.2','xgboost':'3.0.5','joblib':'1.5.2','threadpoolctl':'3.6.0'}
}
data=json.dumps(protocol,ensure_ascii=True,indent=2).encode()
(ev/'protocol.json').write_bytes(data)
(ev/'protocol_sha256.txt').write_text(hashlib.sha256(data).hexdigest(),encoding='utf8')
print(json.dumps({'protocol_sha256':hashlib.sha256(data).hexdigest(),'models':protocol['models'],'inputs_verified':True,'frozen_before_any_fit':True,'predictors':features},ensure_ascii=True))
