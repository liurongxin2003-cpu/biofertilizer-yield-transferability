from pathlib import Path
from collections import defaultdict, Counter
import csv, json, hashlib, math
import openpyxl

root=Path(__file__).resolve().parents[1]; ev=root/'evidence'
previous=root.parent/'Schutz_2018_phase_minus1'
source=root.parent/'Schutz_2018_phase_minus2'/'evidence'/'Supplementary_Data_Sheet_S2.xlsx'
cohort_file=previous/'phase_minus1_cohort_manifest.tsv'
old_splits_file=previous/'evidence'/'frozen_fold_keys.json'
inputs={str(p.relative_to(root.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,cohort_file,old_splits_file,previous/'phase_minus1_predictor_availability.tsv',previous/'PHASE_MINUS1_STATE.md')}
checks={}; issues=[]
def check(name,condition):
 checks[name]=bool(condition)
 if not condition: issues.append(name)
check('source_hash_unchanged',inputs[str(source.relative_to(root.parent))]=='a534b4b18a1d4fd235d7ae4a62d415bbb3c821221a94b826590cc8aa95fec2b2')
cohort=list(csv.DictReader(cohort_file.open(encoding='utf8',newline=''),delimiter='\t'))
keys=sorted((int(r['study_id']),int(r['pair_id'])) for r in cohort)
check('cohort_size1723',len(keys)==1723)
check('cohort_unique1723',len(set(keys))==1723)
check('study_count169',len({s for s,p in keys})==169)
check('cohort_key_hash_frozen',hashlib.sha256(json.dumps(keys,separators=(',',':')).encode()).hexdigest()=='09c2c2da671d53504413548198c869fe1d94cbca37bf964eedd84e2679c2592f')
w=openpyxl.load_workbook(source,read_only=True,data_only=True)
sh=w['Sheet3']; raw=list(sh.iter_rows(values_only=True)); fields=list(raw[0])
check('source_shape',len(raw)==3525 and len(fields)==66)
rows=[dict(zip(fields,vs),excel_row=i+2) for i,vs in enumerate(raw[1:])]
paired=defaultdict(list)
for r in rows: paired[(r['Studynr'],r['Pair'])].append(r)
function_field='N fixing=0, P solubilizing=1, both=2, unclear=3, S or K solubilizer=4, AMF=5'
function_labels={0:'N_fixing',1:'P_solubilizing',2:'N_and_P_traits',3:'unclear',4:'S_or_K_solubilizer',5:'AMF'}
aliases={'Paddy':'rice','rice':'rice','springwheat':'wheat','winter wheat':'wheat','durum wheat':'wheat','Triticum turgidum':'wheat','maize':'maize','silage maize':'maize','kidney bean':'Phaseolus_beans','Snap bean':'Phaseolus_beans','runner bean(Phaseolus coccineus)':'Phaseolus_beans'}
spec=[
 ('crop','Crop','category','crop group',False),
 ('climate_zone','Climatic zone (after Koeppen)','category','Koeppen code',True),
 ('soil_texture','Soil texture','category','source texture label',True),
 ('soil_organic_matter_pct','OM % ','float','percent',True),
 ('soil_ph_water','pH (water)','float','pH',True),
 ('available_p_kg_ha','Available P merged kg/ha','float','kg P/ha; author merged conversion',True),
 ('inoculant_function',function_field,'category_code','code0-5; nominal, not ordinal',False),
 ('n_applied_kg_ha','N applied kg/ha','float','kg N/ha',True),
 ('k_applied_kg_ha','K applied kg/ha','float','kg K/ha',True),
 ('p_applied_kg_ha','P applied kg/ha','float','kg P/ha',True),
 ('irrigation_status','Rainfed/irrigated (0/1)','category_code','0 rainfed;1 irrigated',True),
 ('fertilization_type','Type of fertilization','category','mineral;mix;no;organic',False)
]
def clean(v): return v.strip() if isinstance(v,str) else v
def missing(v): return v is None or isinstance(v,str) and (v.strip().casefold() in ('','na','n/a','nan','none','not specified') or v.strip().startswith('#'))
def numeric(v): return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
def oid(s,p): return f'S{s:03d}_P{p:03d}'
def join_ids(rs): return ','.join(r['observation_id'] for r in sorted(rs,key=lambda r:(r['study_id'],r['pair_id'])))
def keyhash(rs): return hashlib.sha256(('\n'.join(r['observation_id'] for r in sorted(rs,key=lambda r:(r['study_id'],r['pair_id'])))+'\n').encode()).hexdigest()
matrix=[]; audit=[]
for m in sorted(cohort,key=lambda r:(int(r['study_id']),int(r['pair_id']))):
 s,p=int(m['study_id']),int(m['pair_id']); arms=paired[(s,p)]
 ts=[r for r in arms if r['With biofertilizers yes=1/no=0']==1]
 cs=[r for r in arms if r['With biofertilizers yes=1/no=0']==0]
 check(f'pair_{s}_{p}_structure',len(arms)==2 and len(ts)==len(cs)==1)
 if len(ts)!=1 or len(cs)!=1: raise RuntimeError('Unrecoverable source arm structure; see issue list')
 t,c=ts[0],cs[0]
 check(f'pair_{s}_{p}_frozen_source_rows',t['excel_row']==int(m['treatment_excel_row']) and c['excel_row']==int(m['control_excel_row']))
 yt,yc=t['Yield dry (kg/ha)'],c['Yield dry (kg/ha)']
 ym=missing(yt) or missing(yc); valid=numeric(yt) and numeric(yc) and yt>0 and yc>0
 logrr=math.log(yt/yc) if valid else None
 finite=logrr is not None and math.isfinite(logrr)
 direction='increase' if valid and yt>yc else 'decrease' if valid and yt<yc else 'equal' if valid else 'invalid'
 direction_ok=finite and ((yt>yc and logrr>0) or (yt<yc and logrr<0) or (yt==yc and logrr==0))
 check(f'outcome_{s}_{p}_valid',not ym and valid and finite and direction_ok)
 row={'observation_id':oid(s,p),'study_id':s,'pair_id':p,'logRR':logrr if finite else 'NA'}
 for name,field,dtype,unit,allows_na in spec:
  v=clean(t[field])
  if name=='crop':
   v=aliases.get(str(v),str(v)); check(f'crop_{s}_{p}_unchanged',v==m['crop_group'])
  if missing(v):
   v='NA'; check(f'{name}_{s}_{p}_missing_allowed',allows_na)
  elif dtype=='float': check(f'{name}_{s}_{p}_numeric',numeric(v))
  elif dtype=='category_code': check(f'{name}_{s}_{p}_known_code',v in function_labels if name=='inoculant_function' else v in (0,1))
  row[name]=v
 check(f'function_{s}_{p}_unchanged',function_labels.get(row['inoculant_function'])==m['function_group'])
 row.update(treatment_dry_yield_kg_ha=yt if numeric(yt) else 'NA',control_dry_yield_kg_ha=yc if numeric(yc) else 'NA',treatment_excel_row=t['excel_row'],control_excel_row=c['excel_row'])
 matrix.append(row)
 audit.append({'observation_id':row['observation_id'],'study_id':s,'pair_id':p,'treatment_dry_yield_kg_ha':row['treatment_dry_yield_kg_ha'],'control_dry_yield_kg_ha':row['control_dry_yield_kg_ha'],'logRR':row['logRR'],'source_yield_missing':ym,'source_yield_valid_positive_finite':valid,'logRR_missing':logrr is None,'logRR_nonfinite':logrr is not None and not finite,'direction':direction,'direction_consistent':direction_ok})
# Flags are audit-only. Their thresholds never define inclusion, preprocessing or labels.
values=[r['logRR'] for r in matrix if isinstance(r['logRR'],(int,float)) and math.isfinite(r['logRR'])]
def quantile(v,q):
 a=sorted(v); z=(len(a)-1)*q; lo=math.floor(z); hi=math.ceil(z)
 return a[lo]+(a[hi]-a[lo])*(z-lo)
q1=quantile(values,.25); q3=quantile(values,.75); iqr=q3-q1; low=q1-3*iqr; high=q3+3*iqr
threshold=math.log(2)
for r in audit:
 v=r['logRR']; isfinite=isinstance(v,(int,float)) and math.isfinite(v)
 r.update(abs_logRR_gt_ln2=isfinite and abs(v)>threshold,below_Q1_minus3IQR=isfinite and v<low,above_Q3_plus3IQR=isfinite and v>high,retained_in_matrix=True,value_modified=False)
old=json.loads(old_splits_file.read_text(encoding='utf8'))
split_rows=[]
for design,key in [('LOSO','study_id'),('LOCO_whole_study_purge','crop'),('function_secondary_whole_study_purge','inoculant_function')]:
 levels=sorted({r[key] for r in matrix},key=lambda v:str(function_labels[v]) if key=='inoculant_function' else str(v))
 for i,group in enumerate(levels,1):
  test=[r for r in matrix if r[key]==group]; test_studies={r['study_id'] for r in test}
  train=[r for r in matrix if r[key]!=group and r['study_id'] not in test_studies]
  purge=[r for r in matrix if r[key]!=group and r['study_id'] in test_studies]
  train_studies={r['study_id'] for r in train}; purge_studies={r['study_id'] for r in purge}
  held=function_labels[group] if key=='inoculant_function' else group
  fold=f'{design}:{held}'
  for role,rs in [('train',train),('test',test),('purged',purge)]:
   actual={(r['study_id'],r['pair_id']) for r in rs}
   check(f'{fold}_{role}_matches_previous_freeze',fold in old and actual=={tuple(k) for k in old[fold][role]})
  check(f'{fold}_study_isolation',not train_studies&test_studies and not train_studies&purge_studies and purge_studies<=test_studies)
  check(f'{fold}_nonempty_train_test',bool(train) and bool(test))
  trainids={r['observation_id'] for r in train}; testids={r['observation_id'] for r in test}; purgeids={r['observation_id'] for r in purge}
  check(f'{fold}_complete_disjoint_partition',trainids|testids|purgeids=={r['observation_id'] for r in matrix} and not trainids&testids and not trainids&purgeids and not testids&purgeids)
  split_rows.append({'fold_id':fold,'design':design,'fold_number':i,'primary_or_secondary':'secondary' if key=='inoculant_function' else 'primary','held_out_level':held,'held_out_function_code':group if key=='inoculant_function' else 'NA','train_n':len(train),'test_n':len(test),'purge_n':len(purge),'train_studies_n':len(train_studies),'test_studies_n':len(test_studies),'purge_studies_n':len(purge_studies),'study_overlap_train_test':len(train_studies&test_studies),'study_overlap_train_purge':len(train_studies&purge_studies),'train_study_ids':','.join(map(str,sorted(train_studies))),'test_study_ids':','.join(map(str,sorted(test_studies))),'purge_study_ids':','.join(map(str,sorted(purge_studies))) or 'NA','train_ids':join_ids(train),'test_ids':join_ids(test),'purge_ids':join_ids(purge) or 'NA','train_ids_sha256':keyhash(train),'test_ids_sha256':keyhash(test),'purge_ids_sha256':keyhash(purge)})
for design,expected in [('LOSO',169),('LOCO_whole_study_purge',40),('function_secondary_whole_study_purge',6)]:
 ds=[r for r in split_rows if r['design']==design]
 check(design+'_fold_count',len(ds)==expected)
 tested=[i for r in ds for i in r['test_ids'].split(',')]
 check(design+'_every_observation_tested_once',len(tested)==len(set(tested))==1723 and set(tested)=={r['observation_id'] for r in matrix})
check('all215_frozen_folds_preserved',set(r['fold_id'] for r in split_rows)==set(old) and len(split_rows)==215)
dictionary=[]
sources={name:(field,dtype,unit,allows_na) for name,field,dtype,unit,allows_na in spec}
for index,name in enumerate(matrix[0],1):
 if name in sources:
  field,dtype,unit,allows_na=sources[name]; role='main_predictor'
  observed={str(r[name]) for r in matrix if r[name]!='NA'}
  definition='Frozen Phase-1 crop groups; source aliases mapped conservatively' if name=='crop' else '0=N_fixing;1=P_solubilizing;2=N_and_P_traits;3=unclear;4=S_or_K_solubilizer;5=AMF' if name=='inoculant_function' else '0=rainfed;1=irrigated' if name=='irrigation_status' else 'Source treatment-arm metadata, whitespace trimmed only'
  future='No imputation/encoding now; fit any imputation, scaling and vocabulary only in future training folds'
 elif name=='logRR':
  field='Yield dry (kg/ha) from flag1 treatment and flag0 coded comparator';dtype='float';unit='natural logarithm;dimensionless';allows_na=False;role='outcome';definition='log(treatment_dry_yield_kg_ha/control_dry_yield_kg_ha); positive=higher dry yield vs coded comparator; no relabeling';future='Only outcome; cannot be a predictor'
 elif name in ('study_id','pair_id','observation_id'):
  field='Studynr' if name=='study_id' else 'Pair' if name=='pair_id' else 'Studynr + Pair';dtype='integer' if name!='observation_id' else 'string';unit='identifier';allows_na=False;role='grouping_id' if name=='study_id' else 'row_identifier';definition='Entire numeric Studynr is minimum isolation unit' if name=='study_id' else 'Pair ID is unique only within Studynr' if name=='pair_id' else 'S{study_id:03d}_P{pair_id:03d}; sorted numeric (study_id,pair_id)';future='Exclude from predictor input'
 else:
  field='Yield dry (kg/ha)' if 'yield' in name else 'Sheet3 Excel row (header=1)';dtype='float' if 'yield' in name else 'integer';unit='kg/ha dry harvested product' if 'yield' in name else '1-based Excel row';allows_na=False;role='outcome_source_exclude' if 'yield' in name else 'provenance_exclude';definition='flag1=treatment;flag0=coded control; source value unchanged; some comparators contain active inoculants';future='Exclude from predictor input'
 dictionary.append({'column_index':index,'column':name,'role':role,'dtype':dtype,'unit_or_code':unit,'source_field':field,'missing_allowed':allows_na,'missing_token':'NA','missing_n':sum(r[name]=='NA' for r in matrix),'definition':definition,'future_handling':future})
check('exact12_main_predictors',sum(r['role']=='main_predictor' for r in dictionary)==12 and len(matrix[0])==20)
prior_features=list(csv.DictReader((previous/'phase_minus1_predictor_availability.tsv').open(encoding='utf8',newline=''),delimiter='\t'))
for d in dictionary:
 if d['role']=='main_predictor':
  oldfeature=next(r for r in prior_features if r['predictor']==d['source_field'])
  check(d['column']+'_missingness_matches_frozen',d['missing_n']==int(oldfeature['missing']))
tables={'analysis_matrix.tsv':matrix,'split_manifest.tsv':split_rows,'outcome_audit.tsv':audit,'data_dictionary.tsv':dictionary}
summary={'observations':len(matrix),'studies':len({r['study_id'] for r in matrix}),'predictors':12,'matrix_columns':20,'outcome':'logRR','outcome_missing':sum(r['logRR_missing'] for r in audit),'outcome_nonfinite':sum(r['logRR_nonfinite'] for r in audit),'invalid_source_yield':sum(not r['source_yield_valid_positive_finite'] for r in audit),'direction_mismatches':sum(not r['direction_consistent'] for r in audit),'direction_counts':dict(Counter(r['direction'] for r in audit)),'outcome_quantiles':{str(q):quantile(values,q) for q in (0,.01,.05,.25,.5,.75,.95,.99,1)},'extreme_rules':{'absolute_logRR_threshold':threshold,'lower_Q1_minus3IQR':low,'upper_Q3_plus3IQR':high,'quantile_method':'linear interpolation at (n-1)*q','action':'flag only; retain all values and observations'},'extreme_counts':{'abs_logRR_gt_ln2':sum(r['abs_logRR_gt_ln2'] for r in audit),'below_Q1_minus3IQR':sum(r['below_Q1_minus3IQR'] for r in audit),'above_Q3_plus3IQR':sum(r['above_Q3_plus3IQR'] for r in audit)},'no_deletion':True,'no_winsorization':True,'no_preprocessing_fitted':True,'models_trained':0,'prediction_performance_read':0,'inputs_sha256':inputs,'structural_checks_passed':all(checks.values()),'structural_issues':issues,'fold_counts':dict(Counter(r['design'] for r in split_rows))}
(ev/'delivery_tables.json').write_text(json.dumps(tables,ensure_ascii=True,allow_nan=False),encoding='utf8')
(ev/'summary.json').write_text(json.dumps(summary,ensure_ascii=True,indent=2),encoding='utf8')
(ev/'structural_checks.json').write_text(json.dumps(checks,ensure_ascii=True),encoding='utf8')
(ev/'predictor_columns.json').write_text(json.dumps([s[0] for s in spec]),encoding='utf8')
print(json.dumps({k:v for k,v in summary.items() if k!='inputs_sha256'},ensure_ascii=True,indent=2))
