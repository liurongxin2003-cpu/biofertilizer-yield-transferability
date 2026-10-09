from pathlib import Path
import csv,json,hashlib,sys
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'deps'))
import numpy as np
ev=root/'evidence';old=root.parent/'Schutz_2018_phase0'
p=json.loads((ev/'protocol.json').read_text(encoding='utf8'))
t=json.loads((ev/'delivery_tables.json').read_text(encoding='utf8'))
checks=[];hashes={};loaded={}
for name,rows in t.items():
 out=list(csv.DictReader((root/name).open(encoding='utf8',newline=''),delimiter='\t'));assert len(out)==len(rows)
 for r,q in zip(rows,out):
  assert list(r)==list(q)
  for k,v in r.items():
   if isinstance(v,(int,float)):assert float(q[k])==float(v),(name,k,q[k],v)
   else:assert q[k]==v
 loaded[name]=out;hashes[name]=hashlib.sha256((root/name).read_bytes()).hexdigest()
 checks.append('complete_value_roundtrip:'+name)
for name,h in p['input_sha256'].items():assert hashlib.sha256((old/name).read_bytes()).hexdigest()==h
checks.append('all_frozen_input_hashes_unchanged')
source={r['observation_id']:r for r in csv.DictReader((old/'analysis_matrix.tsv').open(encoding='utf8',newline=''),delimiter='\t')}
splits={r['fold_id']:r for r in csv.DictReader((old/'split_manifest.tsv').open(encoding='utf8',newline=''),delimiter='\t')}
seen=set()
for r in loaded['phase1_predictions.tsv']:
 key=(r['task'],r['estimator'],r['observation_id']);assert key not in seen;seen.add(key)
 assert float(r['logRR'])==float(source[r['observation_id']]['logRR'])
 assert r['observation_id'] in splits[r['fold_id']]['test_ids'].split(',')
 assert r['train_ids_sha256']==r['preprocessor_fit_ids_sha256']==splits[r['fold_id']]['train_ids_sha256']
 assert np.isfinite(float(r['prediction']))
assert len(seen)==3*6*1723
checks.extend(['all_31014_predictions_unique_finite_exact_target','all_predictions_are_frozen_test_members'])
z=np.load(root/'run'/'bootstrap_draws.npz');ss=np.random.default_rng(20261006).integers(0,169,size=(2000,169))
assert np.array_equal(ss,z['sampled_study_index'])
assert np.array_equal(np.stack([np.bincount(r,minlength=169) for r in ss]),z['study_multiplicity'])
for r in loaded['phase1_bootstrap.tsv']:
 arr=z[f"{r['task']}__{r['estimator']}__{r['weighting']}__{r['metric']}"]
 if r['kind'].startswith('paired'):arr=arr-z[f"{r['task']}__{r['baseline']}__{r['weighting']}__{r['metric']}"]
 np.testing.assert_allclose([float(r['ci95_lower']),float(r['ci95_upper'])],np.quantile(arr,[.025,.975]),rtol=0,atol=1e-15)
 if r['adjusted_one_sided_upper']!='NA':assert float(r['adjusted_one_sided_upper'])==float(np.quantile(arr,1-.05/12))
assert sum(r['adjusted_one_sided_upper']!='NA' for r in loaded['phase1_bootstrap.tsv'])==12
checks.extend(['2000_seeded_paired_study_draws_reproduced','all270_bootstrap_CIs_and12_adjusted_bounds_verified'])
fit=json.loads((ev/'fit_summary.json').read_text());assert fit['model_fits']==645 and fit['outer_folds']==215 and fit['tuning_fits']==0 and fit['SHAP']==0
checks.append('645_fits_215_outer_folds_zero_tuning_zero_SHAP')
summary={'passed':True,'checks':checks,'failed':[],'output_sha256':hashes,'row_counts':{k:len(v) for k,v in loaded.items()},'formal_decision_written':False}
(ev/'output_QC.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
print(json.dumps(summary))
print(json.dumps([r for r in loaded['phase1_bootstrap.tsv'] if r['kind'].startswith('paired') and r['task'] in ('A','B') and r['weighting']=='study_equal' and r['metric']=='RMSE' and r['baseline']!='zero_response']))
