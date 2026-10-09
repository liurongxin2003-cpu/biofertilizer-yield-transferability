from pathlib import Path
import json,csv,hashlib,collections,re,ast,math
root=Path(__file__).resolve().parents[1]; ev=root/'evidence'
expected=json.loads((ev/'delivery_tables.json').read_text(encoding='utf8'))
validation=json.loads((ev/'delivery_validation.json').read_text(encoding='utf8'))
checks={}
def check(name,condition):
 checks[name]=bool(condition)
 assert condition,name
def same_value(a,v):
 if isinstance(v,bool): return a==str(v).lower()
 if isinstance(v,float):
  try: return math.isclose(float(a),v,rel_tol=1e-14,abs_tol=1e-15)
  except ValueError: return False
 return a==str(v)
for name,rs in expected.items():
 actual=list(csv.DictReader((root/name).open(encoding='utf8',newline=''),delimiter='\t'))
 check(name+'_rows',len(actual)==len(rs))
 check(name+'_values',all(set(a)==set(r) and all(same_value(a[k],v) for k,v in r.items()) for a,r in zip(actual,rs)))
manifest=expected['phase_minus1_cohort_manifest.tsv']; keys={(r['study_id'],r['pair_id']) for r in manifest}
check('cohort_unique1723',len(manifest)==len(keys)==1723)
check('studies169',len({s for s,p in keys})==169)
check('study_counts_reconcile',sum(r['observations'] for r in expected['phase_minus1_study_audit.tsv'])==1723)
for dim in ('study','crop','function','country','climate','environment'):
 check('distribution_'+dim,sum(r['observations'] for r in expected['phase_minus1_distribution.tsv'] if r['dimension']==dim)==1723)
for r in expected['phase_minus1_predictor_availability.tsv']:
 check('missing_count_'+r['predictor'],r['observed']+r['missing']==1723)
folds=json.loads((ev/'frozen_fold_keys.json').read_text(encoding='utf8'))
fold_rows=expected['phase_minus1_fold_audit.tsv']; check('folds215',len(folds)==len(fold_rows)==215)
for design,n in [('LOSO',169),('LOCO_whole_study_purge',40),('function_secondary_whole_study_purge',6)]:
 visited=[]; rows=[r for r in fold_rows if r['design']==design]; check(design+'_fold_count',len(rows)==n)
 for r in rows:
  f=folds[f"{design}:{r['held_out']}"]; train={tuple(k) for k in f['train']}; test={tuple(k) for k in f['test']}; purge={tuple(k) for k in f['purged']}
  check(f"{design}:{r['held_out']}:partition",train|test|purge==keys and not train&test and not train&purge and not test&purge)
  check(f"{design}:{r['held_out']}:isolation",not {k[0] for k in train}&{k[0] for k in test})
  check(f"{design}:{r['held_out']}:counts",len(train)==r['training_observations'] and len(test)==r['test_observations'] and len(purge)==r['purged_observations'])
  visited.extend(test)
 check(design+'_every_row_once',len(visited)==1723 and set(visited)==keys)
check('all_main_predictors_train_covered',all(r['all_missing_training_folds']==0 for r in expected['phase_minus1_coverage_summary.tsv']))
check('no_outcome_values_in_manifest',all(r['outcome_computed'] is False for r in manifest) and all(not re.search('yield|logrr|lnrr|response_ratio',k,re.I) for k in manifest[0]))
check('external_raw_read',all((ev/f'lutz_biomass{y}.csv').exists() for y in (2018,2019,2020)))
check('primary_source_unchanged',hashlib.sha256((root.parent/'Schutz_2018_phase_minus2'/'evidence'/'Supplementary_Data_Sheet_S2.xlsx').read_bytes()).hexdigest()=='a534b4b18a1d4fd235d7ae4a62d415bbb3c821221a94b826590cc8aa95fec2b2')
check('terminal_files_not_yet_written',not (root/'PHASE_MINUS1_REPORT.md').exists() and not (root/'PHASE_MINUS1_STATE.md').exists())
check('15_tsv',len(list(root.glob('*.tsv')))==15)
check('all_previews_reviewed',len(list(ev.glob('preview_*.png')))==15 and (ev/'visual_review.json').exists())
view=json.loads((ev/'visual_review.json').read_text(encoding='utf8'));check('visual_pass',view.get('all_tables_legible') is True)
report=(ev/'report.md.in').read_text(encoding='utf8'); state=(ev/'state.md.in').read_text(encoding='utf8')
check('draft_terminal_placeholder',report.count('__DECISION__')==1 and state.count('__DECISION__')==1)
check('no_computations',json.loads((ev/'audit_summary.json').read_text(encoding='utf8'))['outcome_computations']==0)
# Check our executable metadata audit contains no estimator or target transformations.
for name in ('audit.py','external_audit.py','complete_evidence.py','prepare_report.py'):
 tree=ast.parse((root/'scripts'/name).read_text(encoding='utf8'))
 forbidden={'log','log1p','log2','log10','fit','fit_transform','predict','cross_val_score','train_test_split','shap_values'}
 calls=[n.func.attr if isinstance(n.func,ast.Attribute) else n.func.id if isinstance(n.func,ast.Name) else '' for n in ast.walk(tree) if isinstance(n,ast.Call)]
 check('no_ml_or_log_calls_'+name,not set(calls)&forbidden)
out={'all_checks_verified':all(checks.values()),'checks':checks,'formal_decision_written':False}
(ev/'predecision_validation.json').write_text(json.dumps(out,ensure_ascii=True,indent=2),encoding='utf8')
print(json.dumps({'checks_passed':len(checks),'all_checks_verified':out['all_checks_verified'],'formal_decision_written':False},ensure_ascii=True))
