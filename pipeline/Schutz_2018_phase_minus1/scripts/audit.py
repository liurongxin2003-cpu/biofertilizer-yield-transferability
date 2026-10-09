from pathlib import Path
from collections import Counter, defaultdict
import csv, json, hashlib, math, re
import openpyxl

root = Path(__file__).resolve().parents[1]
old = root.parent / 'Schutz_2018_phase_minus2'
ev = root / 'evidence'
source = old / 'evidence' / 'Supplementary_Data_Sheet_S2.xlsx'
sha = hashlib.sha256(source.read_bytes()).hexdigest()
assert sha == 'a534b4b18a1d4fd235d7ae4a62d415bbb3c821221a94b826590cc8aa95fec2b2'
wb = openpyxl.load_workbook(source, read_only=True, data_only=True)
sheet = wb['Sheet3']; data = list(sheet.iter_rows(values_only=True))
fields = list(data[0]); rows = [dict(zip(fields, values), excel_row=i+2) for i, values in enumerate(data[1:])]
assert len(rows) == 3524 and len(fields) == 66
pairing = list(csv.DictReader((old/'phase_minus2_pairing_audit.tsv').open(encoding='utf8'), delimiter='\t'))
keys = {(int(r['study_id']), int(r['pair_id'])) for r in pairing if r['metadata_consistent_pair'] == 'true'}
sens = {(int(r['study_id']), int(r['pair_id'])) for r in pairing if r['included_as_yield_comparison'] == 'true'}
assert len(keys)==1723 and len({k[0] for k in keys})==169
assert len(sens)==1758 and len({k[0] for k in sens})==171
flag='With biofertilizers yes=1/no=0'; func=fields[38]
def clean(v): return v.strip() if isinstance(v,str) else v
def missing(v): return v is None or isinstance(v,str) and (v.strip().casefold() in ('','na','n/a','nan','none','not specified') or v.strip().startswith('#'))
def num(v): return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
def level(v): return str(clean(v))
aliases={'Paddy':'rice','rice':'rice','springwheat':'wheat','winter wheat':'wheat','durum wheat':'wheat','Triticum turgidum':'wheat','maize':'maize','silage maize':'maize','kidney bean':'Phaseolus_beans','Snap bean':'Phaseolus_beans','runner bean(Phaseolus coccineus)':'Phaseolus_beans'}
labels={0:'N_fixing',1:'P_solubilizing',2:'N_and_P_traits',3:'unclear',4:'S_or_K_solubilizer',5:'AMF'}
paired=defaultdict(list)
for r in rows: paired[(r['Studynr'],r['Pair'])].append(r)
primary=[]; manifest=[]
for k in sorted(keys):
    rs=paired[k]; t=next(r for r in rs if r[flag]==1); c=next(r for r in rs if r[flag]==0)
    assert len(rs)==2
    # Validate source means only. No ratios, effect sizes, logarithms or outcome distribution are computed.
    assert all(num(r['Yield dry (kg/ha)']) and r['Yield dry (kg/ha)']>0 and num(r['SD dry']) and r['SD dry']>=0 and num(r['Replicates']) and r['Replicates']>0 for r in (t,c))
    raw=level(t['Crop']); crop=aliases.get(raw,raw)
    env='greenhouse' if clean(t['Study']) in ('Gravel 2007','Luna 2012','Bernabeu 2015') else 'field'
    p={'study':k[0],'pair':k[1],'t':t,'c':c,'crop':crop,'function':labels[t[func]],'country':level(t['Country']),'climate':level(t['Climatic zone (after Koeppen)']) if not missing(t['Climatic zone (after Koeppen)']) else 'MISSING','environment':env}
    primary.append(p)
    manifest.append({'study_id':k[0],'study_label':clean(t['Study']),'pair_id':k[1],'treatment_excel_row':t['excel_row'],'control_excel_row':c['excel_row'],'crop_raw':raw,'crop_group':crop,'function_group':p['function'],'country':p['country'],'climate':p['climate'],'environment':env,'outcome_computed':False})
studies=defaultdict(list)
for p in primary: studies[p['study']].append(p)
reference_rows=list(csv.DictReader((old/'phase_minus2_study_audit.tsv').open(encoding='utf8'),delimiter='\t'))
refs={int(r['study_id']):r['public_s1_reference'] for r in reference_rows}
distribution=[]
for kind in ['study','crop','function','country','climate','environment']:
    for group in sorted({p[kind] for p in primary},key=str):
        ps=[p for p in primary if p[kind]==group]
        distribution.append({'dimension':kind,'level':group,'observations':len(ps),'studies':len({p['study'] for p in ps}),'observation_share':len(ps)/1723,'field_observations':sum(p['environment']=='field' for p in ps),'study_label':clean(ps[0]['t']['Study']) if kind=='study' else 'NA'})
study_rows=[{'study_id':s,'study_label':clean(ps[0]['t']['Study']),'observations':len(ps),'observation_share':len(ps)/1723,'crops':' | '.join(sorted({p['crop'] for p in ps})),'functions':' | '.join(sorted({p['function'] for p in ps})),'countries':' | '.join(sorted({p['country'] for p in ps})),'climates':' | '.join(sorted({p['climate'] for p in ps})),'publication_year':ps[0]['t']['Year'],'public_s1_reference':refs[s]} for s,ps in sorted(studies.items())]
features=['Crop','Climatic zone (after Koeppen)','Soil texture','OM % ','pH (water)','Available P merged kg/ha',func,'N applied kg/ha','K applied kg/ha','P applied kg/ha','Rainfed/irrigated (0/1)','Type of fertilization']
cats=set([features[i] for i in (0,1,2,6,10,11)])
def value(p,f): return p['crop'] if f=='Crop' else clean(p['t'][f])
availability=[]; level_rows=[]
for f in features:
    observed=[value(p,f) for p in primary if not missing(value(p,f))]
    varying=sum(len({level(value(p,f)) for p in ps if not missing(value(p,f))})>1 for ps in studies.values())
    ls={level(v) for v in observed}
    sole=0
    if f in cats:
        for l in sorted(ls):
            ps=[p for p in primary if not missing(value(p,f)) and level(value(p,f))==l]; ss={p['study'] for p in ps}
            sole+=len(ss)==1
            level_rows.append({'predictor':f,'level':l,'observations':len(ps),'studies':len(ss),'study_ids':','.join(map(str,sorted(ss))),'single_study_level':len(ss)==1})
    else: assert all(num(v) for v in observed)
    availability.append({'predictor':f,'type':'categorical' if f in cats else 'numeric','observations':1723,'observed':len(observed),'missing':1723-len(observed),'missing_fraction':1-len(observed)/1723,'covered_studies':len({p['study'] for p in primary if not missing(value(p,f))}),'within_study_varying_studies':varying,'observed_levels':len(ls) if f in cats else 'NA','single_study_levels':sole if f in cats else 'NA','observed_min':min(observed) if f not in cats else 'NA','observed_max':max(observed) if f not in cats else 'NA','policy':'training-only encoding; separate MISSING and UNKNOWN' if f in cats else 'training-only imputation + missing indicator; retain out-of-range test values'})
folds=[]; fold_features=[]; unseen=[]; fold_keys={}
for design,key in [('LOSO','study'),('LOCO_whole_study_purge','crop'),('function_secondary_whole_study_purge','function')]:
    for group in sorted({p[key] for p in primary},key=str):
        test=[p for p in primary if p[key]==group]; ts={p['study'] for p in test}
        train=[p for p in primary if p[key]!=group and p['study'] not in ts]
        purge=[p for p in primary if p[key]!=group and p['study'] in ts]
        tr={p['study'] for p in train}; assert not ts&tr
        fid=f'{design}:{group}'
        folds.append({'design':design,'held_out':group,'test_observations':len(test),'test_studies':len(ts),'training_observations':len(train),'training_studies':len(tr),'purged_observations':len(purge),'study_overlap':0,'test_study_ids':','.join(map(str,sorted(ts))),'training_study_ids':','.join(map(str,sorted(tr)))})
        fold_keys[fid]={'test':[(p['study'],p['pair']) for p in test],'train':[(p['study'],p['pair']) for p in train],'purged':[(p['study'],p['pair']) for p in purge]}
        for f in features:
            a=[value(p,f) for p in train if not missing(value(p,f))]; b=[value(p,f) for p in test if not missing(value(p,f))]
            new=set(); novel=below=above=0
            if f in cats:
                la={level(v) for v in a}; new={level(v) for v in b}-la
                novel=sum(level(v) in new for v in b)
                for l in sorted(new):
                    matches=[p for p in test if not missing(value(p,f)) and level(value(p,f))==l]
                    unseen.append({'design':design,'held_out':group,'predictor':f,'unseen_level':l,'test_observations':len(matches),'test_studies':len({p['study'] for p in matches})})
            elif a:
                below=sum(v<min(a) for v in b); above=sum(v>max(a) for v in b)
            fold_features.append({'design':design,'held_out':group,'predictor':f,'training_observed':len(a),'test_observed':len(b),'test_missing':len(test)-len(b),'test_observations':len(test),'unseen_level_count':len(new),'unseen_observations':novel,'unseen_fraction_all_test':novel/len(test),'test_below_training_min':below,'test_above_training_max':above,'training_min':min(a) if a and f not in cats else 'NA','training_max':max(a) if a and f not in cats else 'NA','test_all_missing':len(b)==0,'training_all_missing':len(a)==0})
coverage=[]
for design in ['LOSO','LOCO_whole_study_purge','function_secondary_whole_study_purge']:
    for f in features:
        fs=[r for r in fold_features if r['design']==design and r['predictor']==f]
        denom=sum(r['test_observations'] for r in fs); assert denom==1723
        coverage.append({'design':design,'predictor':f,'folds':len(fs),'all_missing_training_folds':sum(r['training_all_missing'] for r in fs),'all_missing_test_folds':sum(r['test_all_missing'] for r in fs),'folds_with_unseen_levels':sum(r['unseen_observations']>0 for r in fs),'unseen_test_observations':sum(r['unseen_observations'] for r in fs),'unseen_fraction':sum(r['unseen_observations'] for r in fs)/denom,'out_of_training_range_observations':sum(r['test_below_training_min']+r['test_above_training_max'] for r in fs),'out_of_training_range_fraction':sum(r['test_below_training_min']+r['test_above_training_max'] for r in fs)/denom})
control_fields=['Studynr','Year','Crop','Variety','Location','Country','Data summarized from 2 years (no=0/yes=1)','First, second, third or forth year data (0/1/2/3)','N applied kg/ha','P applied kg/ha','K applied kg/ha','Type of fertilization','Type of organic fertilizer','Nutrient content organic amendment N kg/ha','Nutrient content organic amendment P kg/ha','Nutrient content organic amendment K (Kg/ha)','Total P input - mineral + organic','Total N input - mineral + organic','Rainfed/irrigated (0/1)','Treatments','Biofertilizer species','Strain','Yield dry (kg/ha)','SD dry','Replicates']
sig=defaultdict(list)
for p in primary: sig[tuple(clean(p['c'][f]) for f in control_fields)].append(p)
dep=[]
for i,ps in enumerate(sorted(sig.values(),key=lambda ps:(ps[0]['study'],ps[0]['pair'])),1):
    if len(ps)>1: dep.append({'candidate_cluster':f'S{ps[0]["study"]}_C{i}','study_id':ps[0]['study'],'observations':len(ps),'pair_ids':','.join(str(p['pair']) for p in ps),'signature_sha256':hashlib.sha256(json.dumps([clean(ps[0]['c'][f]) for f in control_fields],ensure_ascii=True).encode()).hexdigest(),'status':'suspected repeated control; not proven unique experimental control','minimum_isolation_unit':'study'})
sdflag='SD type original type/modelled (0/1)'
summary={'primary_observations':1723,'primary_studies':169,'sensitivity_observations':1758,'sensitivity_studies':171,'source_sha256':sha,'cohort_keys_sha256':hashlib.sha256(json.dumps(sorted(keys),separators=(',',':')).encode()).hexdigest(),'main_predictors':12,'crop_groups':len({p['crop'] for p in primary}),'countries':len({p['country'] for p in primary}),'climate_groups_observed':len({p['climate'] for p in primary}-{'MISSING'}),'function_groups':6,'largest_study_n':max(map(len,studies.values())),'largest_study_share':max(map(len,studies.values()))/1723,'top5_study_share':sum(sorted(map(len,studies.values()),reverse=True)[:5])/1723,'multi_observation_studies':sum(len(ps)>1 for ps in studies.values()),'single_study_crop_groups':sum(r['dimension']=='crop' and r['studies']==1 for r in distribution),'suspected_control_clusters':len(dep),'observations_in_suspected_control_clusters':sum(d['observations'] for d in dep),'studies_with_suspected_controls':len({d['study_id'] for d in dep}),'control_species_nonmissing':sum(not missing(p['c']['Biofertilizer species']) for p in primary),'original_sd_both':sum(all(level(p[a][sdflag])=='0' for a in ('t','c')) for p in primary),'modelled_sd_any':sum(any(level(p[a][sdflag])=='1' for a in ('t','c')) for p in primary),'sd_provenance_missing_any':sum(any(missing(p[a][sdflag]) for a in ('t','c')) for p in primary),'outcome_computations':0,'models_trained':0}
summary['folds']={d:{'folds':sum(f['design']==d for f in folds),'min_training_studies':min(f['training_studies'] for f in folds if f['design']==d),'max_training_studies':max(f['training_studies'] for f in folds if f['design']==d),'min_test_studies':min(f['test_studies'] for f in folds if f['design']==d),'max_purged_observations':max(f['purged_observations'] for f in folds if f['design']==d)} for d in ('LOSO','LOCO_whole_study_purge','function_secondary_whole_study_purge')}
tables={'phase_minus1_cohort_manifest.tsv':manifest,'phase_minus1_study_audit.tsv':study_rows,'phase_minus1_distribution.tsv':distribution,'phase_minus1_predictor_availability.tsv':availability,'phase_minus1_predictor_levels.tsv':level_rows,'phase_minus1_fold_audit.tsv':folds,'phase_minus1_fold_predictor_coverage.tsv':fold_features,'phase_minus1_unseen_levels.tsv':unseen,'phase_minus1_coverage_summary.tsv':coverage,'phase_minus1_dependency_audit.tsv':dep}
(ev/'audit_summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=True),encoding='utf8')
(ev/'delivery_tables.json').write_text(json.dumps(tables,ensure_ascii=True,allow_nan=False),encoding='utf8')
(ev/'frozen_fold_keys.json').write_text(json.dumps(fold_keys,ensure_ascii=True),encoding='utf8')
print(json.dumps(summary,indent=2,ensure_ascii=True))
