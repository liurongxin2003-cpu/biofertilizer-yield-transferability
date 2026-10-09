from pathlib import Path
from collections import Counter,defaultdict
import json,math,re,hashlib
root=Path(__file__).resolve().parents[1]; ev=root/'evidence'
raw=json.loads((ev/'explore.json').read_text(encoding='utf8')); rows=raw['rows_data']; fields=[s['field'] for s in raw['schema']]
flag='With biofertilizers yes=1/no=0'; func=fields[38]; y='Yield dry (kg/ha)'; sd='SD dry'; n='Replicates'
def clean(v): return v.strip() if isinstance(v,str) else v
def missing(v): return v is None or isinstance(v,str) and (v.strip().casefold() in ('','na','n/a','nan','none','not specified') or v.startswith('#'))
def numeric(v): return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
crop_aliases={'Paddy':'rice','rice':'rice','springwheat':'wheat','winter wheat':'wheat','durum wheat':'wheat','Triticum turgidum':'wheat','maize':'maize','silage maize':'maize','kidney bean':'Phaseolus_beans','Snap bean':'Phaseolus_beans','runner bean(Phaseolus coccineus)':'Phaseolus_beans'}
pairs=defaultdict(list)
for r in rows: pairs[(r['Studynr'],r['Pair'])].append(r)
valid=[]; pair_audit=[]; mismatch_counts=Counter()
background=['Year','Crop','Variety','Location','Country','OM % ','pH (water)','N applied kg/ha','P applied kg/ha','K applied kg/ha','Type of fertilization','Rainfed/irrigated (0/1)']
for (study,pair),rs in sorted(pairs.items()):
 ts=[r for r in rs if r[flag]==1]; cs=[r for r in rs if r[flag]==0]; reason=[]; t=ts[0] if len(ts)==1 else None; c=cs[0] if len(cs)==1 else None
 if len(rs)!=2: reason.append('not_two_rows')
 if len(ts)!=1 or len(cs)!=1: reason.append('not_exactly_one_treatment_one_control')
 if t and c:
  if clean(t['Study'])!=clean(c['Study']) or clean(t['Crop'])!=clean(c['Crop']): reason.append('study_or_crop_mismatch')
  if not all(numeric(r[y]) and r[y]>0 for r in (t,c)): reason.append('invalid_yield')
  if not all(numeric(r[n]) and r[n]>0 for r in (t,c)): reason.append('invalid_n')
  if not all(numeric(r[sd]) and r[sd]>=0 for r in (t,c)): reason.append('invalid_sd')
  mismatches=[k for k in background if clean(t[k])!=clean(c[k])]
  mismatch_counts.update(mismatches)
 else: mismatches=[]
 accepted=not reason
 out={'study_id':study,'study_label':rs[0]['Study'],'pair_id':pair,'raw_rows':len(rs),'excel_rows':','.join(str(r['excel_row']) for r in rs),'flags':','.join(str(r[flag]) for r in rs),'included_as_yield_comparison':accepted,'exclusion_reason':';'.join(reason) or 'NA','treatment_yield_dry_kg_ha':t[y] if t else 'NA','control_yield_dry_kg_ha':c[y] if c else 'NA','treatment_SD_dry':t[sd] if t else 'NA','control_SD_dry':c[sd] if c else 'NA','treatment_n':t[n] if t else 'NA','control_n':c[n] if c else 'NA','background_mismatch_fields':' | '.join(mismatches) or 'NA','effect_size_computed':False}
 pair_audit.append(out)
 if accepted:
  env='greenhouse' if t['Study'] in {'Gravel 2007','Luna 2012','Bernabeu 2015'} else 'field'
  original_crop=str(clean(t['Crop']))
  valid.append({'study':study,'label':t['Study'],'pair':pair,'t':t,'c':c,'environment':env,'crop':crop_aliases.get(original_crop,original_crop),'raw_crop':original_crop,'group':t[func],'mismatches':mismatches})
studies=defaultdict(list)
for p in valid: studies[p['study']].append(p)
assert len(set(r['Studynr'] for r in rows))==171
assert all(len({clean(r['Study']) for r in rows if r['Studynr']==i})==1 for i in set(r['Studynr'] for r in rows))
refs=json.loads((ev/'s1_probe.json').read_text(encoding='utf8'))['paragraphs']; numbered={int(m.group(1)):m.group(2) for text in refs if (m:=re.match(r'^(\d+)\.\s*(.*)',text))}
assert set(numbered)==set(range(1,172))
# No unseen-control ID is asserted. These signatures identify repeated *reported* controls only.
control_fields=['Studynr','Year','Crop','Variety','Location','Country','Data summarized from 2 years (no=0/yes=1)','First, second, third or forth year data (0/1/2/3)','N applied kg/ha','P applied kg/ha','K applied kg/ha','Type of fertilization','Type of organic fertilizer','Nutrient content organic amendment N kg/ha','Nutrient content organic amendment P kg/ha','Nutrient content organic amendment K (Kg/ha)','Total P input - mineral + organic','Total N input - mineral + organic','Rainfed/irrigated (0/1)','Treatments','Biofertilizer species','Strain',y,sd,n]
signatures=defaultdict(list)
for p in valid: signatures[tuple(clean(p['c'][k]) for k in control_fields)].append(p)
dependencies=[]
for i,(sig,ps) in enumerate(signatures.items(),1):
 if len(ps)>1:
  dependencies.append({'candidate_control_cluster':f'S{ps[0]["study"]}_C{i}','study_id':ps[0]['study'],'study_label':ps[0]['label'],'n_comparisons':len(ps),'pair_ids':','.join(str(p['pair']) for p in ps),'control_excel_rows':','.join(str(p['c']['excel_row']) for p in ps),'control_yield_dry_kg_ha':ps[0]['c'][y],'control_SD_dry':ps[0]['c'][sd],'control_n':ps[0]['c'][n],'status':'suspected_shared_control_not_proven','safe_split_unit':'whole_study','signature_fields':' | '.join(control_fields)})
dropdup_fields=[k for k in fields if k not in ('Pair','Hilfszeilen')]
duplicate_pairs=defaultdict(list)
for p in valid: duplicate_pairs[(tuple(clean(p['t'][k]) for k in dropdup_fields),tuple(clean(p['c'][k]) for k in dropdup_fields))].append(p)
duplicate_sets=[ps for ps in duplicate_pairs.values() if len(ps)>1]
duplicate_surplus=sum(len(ps)-1 for ps in duplicate_sets)
# Retain source observations in audit; conservative lower bound removes observational duplicates.
conservative=[ps[0] for ps in duplicate_pairs.values()]
core={'Crop','Climatic zone (after Koeppen)','Soil texture','OM % ','pH (water)','Available P merged kg/ha',func,'N applied kg/ha','K applied kg/ha','P applied kg/ha','Rainfed/irrigated (0/1)','Type of fertilization'}
grouping={'Study','Studynr','Pair'}
postoutcome={'NUE Ntot','SD NUE Ntot','Hilfszeilen','Water content % yield','Yield dry (kg/ha)','SD dry','Yield fresh (Kg/ha)','SD','SD in %','SD averages and real in % merged','SD averages and real merged','SD type original type/modelled (0/1)','PUE','PUESD'}
exclude={'With biofertilizers yes=1/no=0'}
sensitivity={'Variety','Year','Replicates','Location','Country','Latitude','Longitude','Soil type','EC (dS/m)','CEC (cmolc/kg)','Available N kg/ha merged','Total N kg/ha merged','Available K merged kg/ha','Method','0-30cm topsoil sampling (yes/no)','Bulk density (kg/dm3)','Bulk density calculated (Post and Kwon 2000)','Bulk density merged (kg/dm3)','Data summarized from 2 years (no=0/yes=1)','First, second, third or forth year data (0/1/2/3)','Sowing density (plants or seeds/ha)','Plot size (m2)','Origin','Application','Propagule count (cfu/g)','Treatments','Biofertilizer species','Strain','P source','N source','Type of organic fertilizer','Nutrient content organic amendment N kg/ha','Nutrient content organic amendment P kg/ha','Nutrient content organic amendment K (Kg/ha)','Total P input - mineral + organic','Total N input - mineral + organic'}
assert core|grouping|postoutcome|exclude|sensitivity==set(fields)
availability=[]; schema=[]
for i,k in enumerate(fields,1):
 vals=[p['t'][k] for p in valid]; observed=[v for v in vals if not missing(v)]; cats=set(str(clean(v)) for v in observed); covered={p['study'] for p in valid if not missing(p['t'][k])}
 role='primary_predictor' if k in core else 'grouping_only' if k in grouping else 'outcome_or_outcome_derived_exclude' if k in postoutcome else 'pairing_only_exclude' if k in exclude else 'sensitivity_only'
 note='Use missing category/indicator; learn transformations only within training studies.' if k in core else 'NA'
 if k=='Available N kg/ha merged': note='102 source cells are #REF!; classify as missing, never zero.'
 if k in {'Available P merged kg/ha','Bulk density merged (kg/dm3)','Bulk density calculated (Post and Kwon 2000)'}: note='Source author-derived values; density mostly estimated; no reconstruction in this audit.'
 if k=='Sowing density (plants or seeds/ha)': note='Mixed plants/ha and kg seed/ha units; sensitivity-only pending unit QC.'
 if k=='Variety': note='31 pairs have cultivar mismatches, some apparent sequential autofill; sensitivity-only pending source QC.'
 if k in {'Biofertilizer species','Strain','Application','Propagule count (cfu/g)','Origin','Treatments'}: note='Free text/mixed syntax; sensitivity-only, do not require unpublished joins.'
 if k in postoutcome: note='Not a predictor. Yield-related fields must not enter predictor matrix.'
 availability.append({'column_index':i,'field':k,'role':role,'comparison_denominator':len(valid),'treatment_nonmissing':len(observed),'treatment_missing':len(valid)-len(observed),'treatment_missing_fraction':(len(valid)-len(observed))/len(valid),'covered_studies':len(covered),'distinct_observed_values':len(cats),'treatment_non_numeric_nonmissing':sum(not numeric(v) for v in observed),'raw_excel_errors':sum(isinstance(r[k],str) and r[k].startswith('#') for r in rows),'note':note})
 schema.append({'column_index':i,'original_field':k,'worksheet':'Sheet3','stored_rows_excluding_header':len(rows),'missing_tokens':'blank;NA;N/A;Excel error','role':role,'native_effect_size_column':False})
study_audit=[]
for i,ps in sorted(studies.items()):
 study_audit.append({'study_id':i,'study_label':ps[0]['label'],'public_s1_reference':numbered[i],'valid_yield_comparisons':len(ps),'conservative_nonduplicate_comparisons':sum(p['study']==i for p in conservative),'field_comparisons':sum(p['environment']=='field' for p in ps),'greenhouse_comparisons':sum(p['environment']=='greenhouse' for p in ps),'unique_crop_labels':len({p['crop'] for p in ps}),'unique_inoculant_function_groups':len({p['group'] for p in ps}),'countries':' | '.join(sorted({str(p['t']['Country']) for p in ps})),'publication_year':ps[0]['t']['Year'],'multi_observation_study':len(ps)>1,'suspected_shared_control_clusters':sum(d['study_id']==i for d in dependencies),'raw_pairing_issues':sum(a['study_id']==i and not a['included_as_yield_comparison'] for a in pair_audit)})
labels={0:'N_fixing',1:'P_solubilizing',2:'N_and_P_traits',3:'unclear',4:'S_or_K_solubilizer',5:'AMF'}
categories=[]; splits=[]
for kind in ('crop','inoculant_function'):
 key='crop' if kind=='crop' else 'group'; groups=sorted({p[key] for p in valid},key=str)
 for group in groups:
  test=[p for p in valid if p[key]==group]; ts={p['study'] for p in test}; train=[p for p in valid if p[key]!=group and p['study'] not in ts]; purged=[p for p in valid if p[key]!=group and p['study'] in ts]
  name=labels[group] if kind=='inoculant_function' else group
  categories.append({'category_type':kind,'category':name,'raw_code':group,'yield_comparisons':len(test),'conservative_nonduplicate_comparisons':sum(p[key]==group for p in conservative),'studies':len(ts),'field_comparisons':sum(p['environment']=='field' for p in test),'field_studies':len({p['study'] for p in test if p['environment']=='field'})})
  assert ts.isdisjoint({p['study'] for p in train})
  splits.append({'design':'leave_one_'+kind+'_out_with_study_purge','held_out_group':name,'test_comparisons':len(test),'test_studies':len(ts),'training_comparisons_after_purge':len(train),'training_studies_after_purge':len({p['study'] for p in train}),'purged_same_study_other_group_comparisons':len(purged),'study_overlap_after_purge':0,'technically_executable':bool(train and test),'meaning':'unseen crop context' if kind=='crop' else 'unseen functional category, not unseen microbial strain/species'})
for i,ps in sorted(studies.items()):
 splits.append({'design':'leave_one_study_out','held_out_group':str(i),'test_comparisons':len(ps),'test_studies':1,'training_comparisons_after_purge':len(valid)-len(ps),'training_studies_after_purge':len(studies)-1,'purged_same_study_other_group_comparisons':0,'study_overlap_after_purge':0,'technically_executable':True,'meaning':'new published study; institution independence not guaranteed'})
field=[p for p in valid if p['environment']=='field']; gh=[p for p in valid if p['environment']=='greenhouse']
sdtype='SD type original type/modelled (0/1)'
metrics={'xlsx_stored_rows_with_header':3525,'xlsx_columns':66,'raw_arm_records':len(rows),'raw_pair_keys':len(pairs),'structurally_valid_yield_comparisons':len(valid),'conservative_nonduplicate_yield_comparisons':len(conservative),'pairing_issues':sum(not a['included_as_yield_comparison'] for a in pair_audit),'yield_studies':len(studies),'study_id_missing_raw_records':sum(missing(r['Studynr']) or missing(r['Study']) for r in rows),'field_yield_comparisons':len(field),'field_studies':len({p['study'] for p in field}),'greenhouse_yield_comparisons':len(gh),'greenhouse_studies':len({p['study'] for p in gh}),'raw_crop_labels':len({r['Crop'] for r in rows}),'trimmed_crop_labels':len({p['crop'] for p in valid}),'inoculant_function_groups':len({p['group'] for p in valid}),'primary_predictors':len(core),'suspected_shared_control_clusters':len(dependencies),'comparisons_in_suspected_shared_control_clusters':sum(d['n_comparisons'] for d in dependencies),'studies_with_suspected_shared_controls':len({d['study_id'] for d in dependencies}),'duplicate_pair_signature_sets':len(duplicate_sets),'duplicate_pair_signature_surplus':duplicate_surplus,'multi_observation_studies':sum(len(ps)>1 for ps in studies.values()),'native_lnRR_or_effect_size_columns':0,'pairs_with_positive_treatment_control_yield':sum(p['t'][y]>0 and p['c'][y]>0 for p in valid),'pairs_with_complete_yield_SD_n':sum(all(numeric(r[k]) for r in (p['t'],p['c']) for k in (y,sd,n)) for p in valid),'original_SD_flag_both_arms':sum(str(p['t'][sdtype])=='0' and str(p['c'][sdtype])=='0' for p in valid),'modelled_SD_flag_any_arm':sum(any(str(p[r][sdtype])=='1' for r in ('t','c')) for p in valid),'SD_flag_missing_any_arm':sum(any(missing(p[r][sdtype]) for r in ('t','c')) for p in valid),'pairs_with_background_field_mismatch':sum(bool(p['mismatches']) for p in valid),'control_species_nonmissing_pairs':sum(not missing(p['c']['Biofertilizer species']) for p in valid),'outcome_recomputed':0,'ML_models_trained':0}
counts=[{'metric':k,'value':v,'scope':'yield comparisons only; no NUE/PUE observations counted'} for k,v in metrics.items()]
metrics['trimmed_crop_labels']=len({p['raw_crop'] for p in valid})
metrics['conservative_crop_groups_for_holdout']=len({p['crop'] for p in valid})
for a in counts:
 if a['metric']=='trimmed_crop_labels': a['value']=metrics['trimmed_crop_labels']
counts.append({'metric':'conservative_crop_groups_for_holdout','value':metrics['conservative_crop_groups_for_holdout'],'scope':'alias-merged crop groups; no outcome changes'})
for original_crop in sorted({p['raw_crop'] for p in valid}):
 ps=[p for p in valid if p['raw_crop']==original_crop]
 categories.append({'category_type':'original_crop_label_to_holdout_mapping','category':original_crop,'raw_code':ps[0]['crop'],'yield_comparisons':len(ps),'conservative_nonduplicate_comparisons':len(ps),'studies':len({p['study'] for p in ps}),'field_comparisons':sum(p['environment']=='field' for p in ps),'field_studies':len({p['study'] for p in ps if p['environment']=='field'})})
strict=[p for p in valid if not p['mismatches']]
strict_metrics={'metadata_consistent_comparisons':len(strict),'metadata_consistent_studies':len({p['study'] for p in strict}),'metadata_consistent_field_comparisons':sum(p['environment']=='field' for p in strict),'metadata_consistent_field_studies':len({p['study'] for p in strict if p['environment']=='field'})}
metrics.update(strict_metrics); counts.extend({'metric':k,'value':v,'scope':'conservative bound after withholding 35 metadata-mismatch pairs'} for k,v in strict_metrics.items())
for a in pair_audit: a['metadata_consistent_pair']=a['included_as_yield_comparison'] and a['background_mismatch_fields']=='NA'
for kind in ('crop','inoculant_function'):
 key='crop' if kind=='crop' else 'group'
 for group in sorted({p[key] for p in strict},key=str):
  test=[p for p in strict if p[key]==group]; ts={p['study'] for p in test}; train=[p for p in strict if p[key]!=group and p['study'] not in ts]; purged=[p for p in strict if p[key]!=group and p['study'] in ts]
  assert ts.isdisjoint({p['study'] for p in train})
  splits.append({'design':'metadata_consistent_leave_one_'+kind+'_out_with_study_purge','held_out_group':labels[group] if kind=='inoculant_function' else group,'test_comparisons':len(test),'test_studies':len(ts),'training_comparisons_after_purge':len(train),'training_studies_after_purge':len({p['study'] for p in train}),'purged_same_study_other_group_comparisons':len(purged),'study_overlap_after_purge':0,'technically_executable':bool(train and test),'meaning':'technical split audit only; no ML'})
tables={'phase_minus2_schema.tsv':schema,'phase_minus2_counts.tsv':counts,'phase_minus2_missingness.tsv':availability,'phase_minus2_study_audit.tsv':study_audit,'phase_minus2_category_audit.tsv':categories,'phase_minus2_pairing_audit.tsv':pair_audit,'phase_minus2_dependency_audit.tsv':dependencies,'phase_minus2_split_audit.tsv':splits}
summary={'metrics':metrics,'background_mismatch_counts':dict(mismatch_counts),'invalid_pairs':[a for a in pair_audit if not a['included_as_yield_comparison']],'duplicate_pairs':[{'study_id':ps[0]['study'],'pair_ids':[p['pair'] for p in ps]} for ps in duplicate_sets],'yield_min':min(r[y] for r in rows),'yield_max':max(r[y] for r in rows),'min_training_studies_crop':min(a['training_studies_after_purge'] for a in splits if a['design'].startswith('leave_one_crop')),'min_training_studies_function':min(a['training_studies_after_purge'] for a in splits if a['design'].startswith('leave_one_inoculant')),'numeric_study_id_unambiguous_label':True,'author_year_label_not_unique':True,'source_sha256':hashlib.sha256((ev/'Supplementary_Data_Sheet_S2.xlsx').read_bytes()).hexdigest()}
(ev/'audit_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf8')
(ev/'delivery_tables.json').write_text(json.dumps(tables,indent=2),encoding='utf8')
print(json.dumps(summary,ensure_ascii=True))
