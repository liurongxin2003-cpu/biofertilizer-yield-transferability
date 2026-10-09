from pathlib import Path
from collections import Counter
from decimal import Decimal, localcontext
import json,csv,hashlib,math,ast
root=Path(__file__).resolve().parents[1];ev=root/'evidence'
expected=json.loads((ev/'delivery_tables.json').read_text(encoding='utf8'))
checks={};errors=[]
def check(name,condition):
 checks[name]=bool(condition)
 if not condition: errors.append(name)
def matches(a,v):
 if isinstance(v,bool): return a==str(v).lower()
 if isinstance(v,float):
  try:return math.isclose(float(a),v,rel_tol=1e-14,abs_tol=1e-15)
  except ValueError:return False
 return a==str(v)
tables={}
for name,rs in expected.items():
 actual=list(csv.DictReader((root/name).open(encoding='utf8',newline=''),delimiter='\t'));tables[name]=actual
 check(name+'_row_count',len(actual)==len(rs))
 check(name+'_preserved_values',all(set(a)==set(r) and all(matches(a[k],v) for k,v in r.items()) for a,r in zip(actual,rs)))
matrix=tables['analysis_matrix.tsv'];ids={r['observation_id'] for r in matrix};lookup={r['observation_id']:r for r in matrix}
check('unique1723',len(ids)==len(matrix)==1723)
check('169studies',len({r['study_id'] for r in matrix})==169)
check('dictionary_all20_columns',set(r['column'] for r in tables['data_dictionary.tsv'])==set(matrix[0]) and len(matrix[0])==20)
features=[r['column'] for r in tables['data_dictionary.tsv'] if r['role']=='main_predictor']
check('exact_frozen12',features==json.loads((ev/'predictor_columns.json').read_text(encoding='utf8')) and len(features)==12)
max_error=0
with localcontext() as ctx:
 ctx.prec=50
 for r in matrix:
  o=r['observation_id'];yt=float(r['treatment_dry_yield_kg_ha']);yc=float(r['control_dry_yield_kg_ha']);v=float(r['logRR'])
  check(o+'_positive_finite_sources',math.isfinite(yt) and math.isfinite(yc) and yt>0 and yc>0)
  oracle=float((Decimal.from_float(yt)/Decimal.from_float(yc)).ln())
  error=abs(v-oracle);max_error=max(max_error,error)
  check(o+'_decimal50_logRR_verified',math.isfinite(v) and math.isclose(v,oracle,rel_tol=2e-14,abs_tol=2e-15))
  check(o+'_direction',(yt>yc and v>0) or (yt<yc and v<0) or (yt==yc and v==0))
audit=tables['outcome_audit.tsv']
check('outcome_audit_one_per_observation',len(audit)==1723 and {r['observation_id'] for r in audit}==ids)
check('all_retained_unmodified',all(r['retained_in_matrix']=='true' and r['value_modified']=='false' for r in audit))
check('audit_matrix_outcome_consistent',all(r['logRR']==lookup[r['observation_id']]['logRR'] for r in audit))
def unpack(s):return [] if s=='NA' else s.split(',')
def ihash(a):return hashlib.sha256(('\n'.join(a)+'\n').encode()).hexdigest()
folds=tables['split_manifest.tsv'];old=json.loads((root.parent/'Schutz_2018_phase_minus1'/'evidence'/'frozen_fold_keys.json').read_text(encoding='utf8'))
check('215_folds',len(folds)==215 and {r['fold_id'] for r in folds}==set(old))
for r in folds:
 f=r['fold_id'];roles={role:unpack(r[role+'_ids']) for role in ('train','test','purge')}; sets={k:set(v) for k,v in roles.items()}
 check(f+'_counted_unique_ids',all(len(a)==len(sets[role])==int(r[role+'_n']) for role,a in roles.items()))
 check(f+'_complete_disjoint',set.union(*sets.values())==ids and not sets['train']&sets['test'] and not sets['train']&sets['purge'] and not sets['test']&sets['purge'])
 studies={k:{lookup[o]['study_id'] for o in a} for k,a in roles.items()}
 check(f+'_study_no_leakage',not studies['train']&studies['test'] and not studies['train']&studies['purge'] and studies['purge']<=studies['test'] and r['study_overlap_train_test']=='0' and r['study_overlap_train_purge']=='0')
 for role,a in roles.items():
  check(f+'_'+role+'_id_hash',ihash(a)==r[role+'_ids_sha256'])
  check(f+'_'+role+'_study_count',len(studies[role])==int(r[role+'_studies_n']))
  key='purged' if role=='purge' else role
  converted={(int(lookup[o]['study_id']),int(lookup[o]['pair_id'])) for o in a}
  check(f+'_'+role+'_unchanged_old_fold',converted=={tuple(p) for p in old[f][key]})
 for o in sets['test']:
  row=lookup[o]
  check(f+'_'+o+'_test_definition',row['study_id']==r['held_out_level'] if r['design']=='LOSO' else row['crop']==r['held_out_level'] if r['design']=='LOCO_whole_study_purge' else row['inoculant_function']==r['held_out_function_code'])
for design,n in [('LOSO',169),('LOCO_whole_study_purge',40),('function_secondary_whole_study_purge',6)]:
 ds=[r for r in folds if r['design']==design];tested=[o for r in ds for o in unpack(r['test_ids'])]
 check(design+'_fold_count',len(ds)==n)
 check(design+'_each_test_once',Counter(tested)==Counter(ids))
source_summary=json.loads((ev/'summary.json').read_text(encoding='utf8'))
for relative,digest in source_summary['inputs_sha256'].items():check('source_unchanged_'+relative,hashlib.sha256((root.parent/relative).read_bytes()).hexdigest()==digest)
for script in ('prepare.py','deliver.mjs'):
 text=(root/'scripts'/script).read_text(encoding='utf8')
 check(script+'_no_ml_imports',not any(k in text for k in ('import sklearn','from sklearn','import shap','import xgboost','import lightgbm','import torch','import tensorflow')))
 if script.endswith('.py'):
  tree=ast.parse(text);calls={n.func.attr for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute)}
  check(script+'_no_estimator_calls',not calls&{'fit','fit_transform','predict','predict_proba'})
check('no_terminal_state_before_finalization',not (root/'PHASE0_REPORT.md').exists() and not (root/'PHASE0_STATE.md').exists())
check('all_four_previews_reviewed',json.loads((ev/'visual_review.json').read_text(encoding='utf8')).get('all_tables_legible') is True)
out={'all_output_checks_passed':all(checks.values()),'all_structural_checks_passed':source_summary['structural_checks_passed'] and all(checks.values()),'checks_count':len(checks),'failed_checks':errors,'decimal50_max_absolute_error':max_error,'file_sha256':{name:hashlib.sha256((root/name).read_bytes()).hexdigest() for name in expected}}
(ev/'predecision_validation.json').write_text(json.dumps(out,ensure_ascii=True,indent=2),encoding='utf8')
print(json.dumps(out,ensure_ascii=True,indent=2))
