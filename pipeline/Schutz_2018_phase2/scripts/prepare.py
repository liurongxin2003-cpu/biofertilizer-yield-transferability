from pathlib import Path
from collections import defaultdict,Counter
import csv,json,math,hashlib,datetime
import openpyxl
root=Path(__file__).resolve().parents[1];ev=root/'evidence'; ev.mkdir(parents=True,exist_ok=True);(root/'run').mkdir(exist_ok=True)
base=root.parent;p0=base/'Schutz_2018_phase0';p1=base/'Schutz_2018_phase1';pm1=base/'Schutz_2018_phase_minus1';pm2=base/'Schutz_2018_phase_minus2'
assert not list((root/'run').glob('*/fold_*.json')),'Cannot freeze after Phase2 fits'
read=lambda p:list(csv.DictReader(p.open(encoding='utf8',newline=''),delimiter='\t'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
old_protocol=json.loads((p1/'evidence'/'protocol.json').read_text(encoding='utf8'))
assert sha(p1/'evidence'/'protocol.json')=='1420e0b80930256655aa60e145bc8b8ae93823a080ffa4c4a40fe0e64dfcbb0d'
for n,h in old_protocol['input_sha256'].items():assert sha(p0/n)==h
old_qc=json.loads((p1/'evidence'/'output_QC.json').read_text(encoding='utf8'))
for n,h in old_qc['output_sha256'].items():assert sha(p1/n)==h
assert 'decision: `LIMITED_GENERALIZABLE_SIGNAL`' in (p1/'PHASE1_STATE.md').read_text(encoding='utf8')
source=pm2/'evidence'/'Supplementary_Data_Sheet_S2.xlsx';assert sha(source)=='a534b4b18a1d4fd235d7ae4a62d415bbb3c821221a94b826590cc8aa95fec2b2'
files=[p0/n for n in old_protocol['input_sha256']]+[p1/n for n in old_qc['output_sha256']]+[source,pm2/'phase_minus2_pairing_audit.tsv',pm1/'phase_minus1_validation_freeze.tsv',p1/'evidence'/'protocol.json',p1/'PHASE1_STATE.md',p1/'evidence'/'evaluation_summary.json']
inputs={str(p.relative_to(base)):sha(p) for p in files}
main=read(p0/'analysis_matrix.tsv');original_splits=read(p0/'split_manifest.tsv')
dictionary=read(p0/'data_dictionary.tsv')
spec={r['column']:r for r in dictionary if r['role']=='main_predictor'}
assert list(spec)==old_protocol['predictors'] and len(spec)==12
raw=list(openpyxl.load_workbook(source,read_only=True,data_only=True)['Sheet3'].iter_rows(values_only=True))
fields=list(raw[0]);assert len(raw)==3525 and len(fields)==66
arms=defaultdict(list)
for i,values in enumerate(raw[1:],2):
 r=dict(zip(fields,values));r['excel_row']=i;arms[(int(r['Studynr']),int(r['Pair']))].append(r)
pair_audit=read(pm2/'phase_minus2_pairing_audit.tsv')
approved=[r for r in pair_audit if r['included_as_yield_comparison']=='true']
assert len(approved)==1758 and len(set(r['study_id'] for r in approved))==171
aliases={'Paddy':'rice','rice':'rice','springwheat':'wheat','winter wheat':'wheat','durum wheat':'wheat','Triticum turgidum':'wheat','maize':'maize','silage maize':'maize','kidney bean':'Phaseolus_beans','Snap bean':'Phaseolus_beans','runner bean(Phaseolus coccineus)':'Phaseolus_beans'}
function_labels={0:'N_fixing',1:'P_solubilizing',2:'N_and_P_traits',3:'unclear',4:'S_or_K_solubilizer',5:'AMF'}
def missing(v):return v is None or isinstance(v,str) and (v.strip().casefold() in ('','na','n/a','nan','none','not specified') or v.strip().startswith('#'))
full=[];extra_audit=[];byid={r['observation_id']:r for r in main}
for audit in sorted(approved,key=lambda r:(int(r['study_id']),int(r['pair_id']))):
 s,p=int(audit['study_id']),int(audit['pair_id']);rs=arms[(s,p)]
 assert len(rs)==2
 t=next(r for r in rs if r['With biofertilizers yes=1/no=0']==1);c=next(r for r in rs if r['With biofertilizers yes=1/no=0']==0)
 yt,yc=t['Yield dry (kg/ha)'],c['Yield dry (kg/ha)'];assert math.isfinite(yt) and math.isfinite(yc) and yt>0 and yc>0
 assert {t['excel_row'],c['excel_row']}==set(map(int,audit['excel_rows'].split(',')))
 assert yt==float(audit['treatment_yield_dry_kg_ha']) and yc==float(audit['control_yield_dry_kg_ha'])
 oid=f'S{s:03d}_P{p:03d}'
 r={'observation_id':oid,'study_id':str(s),'pair_id':str(p),'logRR':repr(math.log(yt/yc))}
 for name,d in spec.items():
  v=t[d['source_field']];v=v.strip() if isinstance(v,str) else v
  if name=='crop':v=aliases.get(str(v),str(v))
  if missing(v):v='NA'
  elif d['dtype']=='float':assert isinstance(v,(float,int)) and math.isfinite(v)
  elif d['dtype']=='category_code':v=int(v);assert v in (function_labels if name=='inoculant_function' else (0,1))
  r[name]=str(v)
 r.update(treatment_dry_yield_kg_ha=repr(float(yt)),control_dry_yield_kg_ha=repr(float(yc)),treatment_excel_row=str(t['excel_row']),control_excel_row=str(c['excel_row']))
 if oid in byid:
  old=byid[oid];assert set(r)==set(old)
  for k in r:
   if k=='logRR' or k in old_protocol['numeric_predictors'] or k in ('treatment_dry_yield_kg_ha','control_dry_yield_kg_ha'):
    assert (r[k]==old[k]=='NA') or float(r[k])==float(old[k]),(oid,k,r[k],old[k])
   else:assert r[k]==old[k],(oid,k,r[k],old[k])
  r=old.copy() # Preserve shared rows and frozen logRR literally.
 else:
  assert audit['metadata_consistent_pair']=='false'
  extra_audit.append({'observation_id':oid,'study_id':s,'pair_id':p,'mismatch_fields':audit['background_mismatch_fields'],'logRR':float(r['logRR']),'treatment_excel_row':t['excel_row'],'control_excel_row':c['excel_row'],'metadata_source':'treatment arm; unchanged source values'})
 full.append(r)
assert len(extra_audit)==35 and set(byid)<=set(r['observation_id'] for r in full)
def digest(ids):return hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()
def splits_for(rows):
 out=[]
 for task,design in old_protocol['designs'].items():
  key={'A':'study_id','B':'crop','C':'inoculant_function'}[task]
  levels=sorted(set(r[key] for r in rows),key=lambda v:function_labels[int(v)] if task=='C' else str(v))
  for i,g in enumerate(levels,1):
   test=[r for r in rows if r[key]==g];studies={r['study_id'] for r in test}
   train=[r for r in rows if r['study_id'] not in studies];purge=[r for r in rows if r['study_id'] in studies and r[key]!=g]
   assert train and test and len(train+test+purge)==len(rows)
   held=function_labels[int(g)] if task=='C' else g
   f={'fold_id':design+':'+held,'design':design,'task':task,'fold_number':i,'held_out_level':held,'primary_or_secondary':'secondary' if task=='C' else 'primary'}
   for role,rr in (('train',train),('test',test),('purge',purge)):
    ids=[r['observation_id'] for r in rr];f[role+'_ids']=ids;f[role+'_ids_sha256']=digest(ids);f[role+'_n']=len(ids);f[role+'_studies_n']=len(set(r['study_id'] for r in rr))
   assert not {r['study_id'] for r in train}&{r['study_id'] for r in test+purge}
   out.append(f)
 return out
fullsplits=splits_for(full)
main_splits=[]
for i,f in enumerate(original_splits):
 if f['design']==old_protocol['designs']['C']:continue
 out={k:f[k] for k in ('fold_id','design','held_out_level','primary_or_secondary')};out['task']='A' if f['design']=='LOSO' else 'B';out['fold_number']=int(f['fold_number']);out['seed_index']=i
 for role in ('train','test','purge'):
  out[role+'_ids']=[] if f[role+'_ids']=='NA' else f[role+'_ids'].split(',')
  for tail in ('n','studies_n'):out[role+'_'+tail]=int(f[role+'_'+tail])
  out[role+'_ids_sha256']=f[role+'_ids_sha256'];assert digest(out[role+'_ids'])==f[role+'_ids_sha256']
 main_splits.append(out)
for i,f in enumerate(fullsplits):f['seed_index']=i
for fs,rows in ((fullsplits,full),(main_splits,main)):
 for task in set(f['task'] for f in fs):
  tested=[oid for f in fs if f['task']==task for oid in f['test_ids']]
  assert len(tested)==len(set(tested))==len(rows) and set(tested)==set(r['observation_id'] for r in rows)
assert Counter(f['task'] for f in fullsplits)=={'A':171,'B':40,'C':6}
assert len(main_splits)==209
variants={'expanded_1758_12':{'rows':full,'splits':fullsplits,'predictors':old_protocol['predictors'],'categorical_predictors':old_protocol['categorical_predictors'],'numeric_predictors':old_protocol['numeric_predictors'],'cohort_n':1758,'studies_n':171},'main_1723_crop_blind_11':{'rows':main,'splits':main_splits,'predictors':[v for v in old_protocol['predictors'] if v!='crop'],'categorical_predictors':[v for v in old_protocol['categorical_predictors'] if v!='crop'],'numeric_predictors':old_protocol['numeric_predictors'],'cohort_n':1723,'studies_n':169}}
data=json.dumps(variants,ensure_ascii=True,indent=2,allow_nan=False).encode();(ev/'frozen_variants.json').write_bytes(data)
nc=Counter(r['study_id'] for r in main);cropgroups=defaultdict(list)
for r in main:cropgroups[r['crop']].append(r)
cropbands={c:'singleton' if len(set(r['study_id'] for r in rs))==1 else 'high_support' if len(set(r['study_id'] for r in rs))>=3 and len(rs)>=20 else 'low_support' for c,rs in cropgroups.items()}
metadata={'study_size_large_ge20':{'studies':sum(n>=20 for n in nc.values()),'observations':sum(n for n in nc.values() if n>=20)},'study_size_small_lt20':{'studies':sum(n<20 for n in nc.values()),'observations':sum(n for n in nc.values() if n<20)},'crop_support':{g:{'crops':sum(v==g for v in cropbands.values()),'observations':sum(len(cropgroups[c]) for c,v in cropbands.items() if v==g),'studies':len(set(r['study_id'] for c,v in cropbands.items() if v==g for r in cropgroups[c]))} for g in ('high_support','low_support','singleton')}}
newprotocol={'frozen_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'phase1_protocol_sha256':sha(p1/'evidence'/'protocol.json'),'inputs_sha256':inputs,'variant_file_sha256':hashlib.sha256(data).hexdigest(),'variants':{v:{k:x[k] for k in ('cohort_n','studies_n','predictors','categorical_predictors','numeric_predictors')} for v,x in variants.items()},'models':old_protocol['models'],'training_weights':old_protocol['training_weights'],'preprocessing':old_protocol['preprocessing'],'base_seed':20261006,'model_seed_policy':'Phase1:base_seed+frozen variant split manifest index; cropblind preserves original main manifest indices','baselines':old_protocol['baselines'],'metrics':old_protocol['metrics'],'metric_weights':old_protocol['metric_weights'],'bootstrap':dict(old_protocol['bootstrap'],resampled_clusters_per_draw='variant original study count169 or171'),'decision_rule':old_protocol['decision_rule'],'libraries':old_protocol['libraries'],'phase2_decision':'ROBUST_LIMITED_GENERALIZABILITY iff neither A nor B has a Phase1-supported model in either sensitivity; otherwise SENSITIVITY_DEPENDENT. Small point changes or descriptive subgroup gains alone cannot reverse the conclusion.','failure_audit':{'source':'unchanged Phase1 main12 OOF only; all3models+3baselines; A/B/C','large_study':'fullmaincohort study size>=20 comparisons; small<20; fixed before Phase2 fitting/audit','crop_support':'high>=3 studies AND>=20 comparisons; singleton=1study; low=allother; thresholds inheritedPhaseMinus1','study_error_quantiles':[0,.1,.25,.5,.75,.9,.95,1],'subset_metrics':'descriptive only, no subgroup bootstrap/gates; study-equal weights recomputed within each descriptive subset; no retraining or selection'},'failure_audit_metadata_counts':metadata,'out_of_scope':['expanded_cohort_x_cropblind_factorial','unweighted training','new splits','new models','hyperparameter search','SHAP','feature hunting','deletion/winsorization','external_validation']}
payload=json.dumps(newprotocol,ensure_ascii=True,indent=2).encode();(ev/'protocol.json').write_bytes(payload);(ev/'protocol_sha256.txt').write_text(hashlib.sha256(payload).hexdigest(),encoding='utf8')
(ev/'cohort_reconstruction_QC.json').write_text(json.dumps({'all1723_shared_rows_exactly_preserved':True,'expanded_cohort_n':1758,'expanded_studies_n':171,'added_comparisons':extra_audit,'added_studies':sorted(set(r['study_id'] for r in full)-set(nc),key=int),'crop_count':len(cropgroups),'fold_counts':{v:dict(Counter(f['task'] for f in x['splits'])) for v,x in variants.items()},'analysis_X_role_validation':True,'source_xlsx_actual_read':True,'structural_errors':[]},ensure_ascii=True,indent=2),encoding='utf8')
print(json.dumps({'phase2_protocol_sha256':hashlib.sha256(payload).hexdigest(),'variants':newprotocol['variants'],'fold_counts':{v:dict(Counter(f['task'] for f in x['splits'])) for v,x in variants.items()},'failure_audit_metadata_counts':metadata,'frozen_before_phase2_fit':True},ensure_ascii=True))
