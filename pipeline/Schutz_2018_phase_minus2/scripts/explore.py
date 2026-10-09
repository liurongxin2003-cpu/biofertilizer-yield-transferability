from pathlib import Path
from collections import Counter,defaultdict
import json,openpyxl,math
root=Path(__file__).resolve().parents[1]
wb=openpyxl.load_workbook(root/'evidence/Supplementary_Data_Sheet_S2.xlsx',data_only=True)
s=wb.active; header=[c.value for c in s[1]]
rows=[dict(zip(header,values),excel_row=i) for i,values in enumerate(s.iter_rows(min_row=2,values_only=True),2) if any(v is not None for v in values)]
def missing(v): return v is None or isinstance(v,str) and v.strip().lower() in ('','na','n/a','nan','none','not specified')
schema=[]
for i,k in enumerate(header,1):
 vals=[r[k] for r in rows if not missing(r[k])]; counts=Counter(str(v) for v in vals)
 schema.append({'column_index':i,'field':k,'nonmissing':len(vals),'missing':len(rows)-len(vals),'unique':len(counts),'top10':counts.most_common(10)})
pairs=defaultdict(list)
for r in rows: pairs[(r['Studynr'],r['Pair'])].append(r)
issues=[{'key':k,'n':len(v),'flags':Counter(str(r['With biofertilizers yes=1/no=0']) for r in v),'rows':[r['excel_row'] for r in v]} for k,v in pairs.items() if len(v)!=2 or set(r['With biofertilizers yes=1/no=0'] for r in v)!={0,1}]
out={'rows':len(rows),'dimensions':[s.max_row,s.max_column],'schema':schema,'unique_study_labels':len(set(r['Study'] for r in rows)),'unique_studynr':len(set(r['Studynr'] for r in rows)),'pair_count':len(pairs),'pair_size_distribution':Counter(len(v) for v in pairs.values()),'pair_issues':issues,'last_rows':rows[-4:],'rows_data':rows}
(root/'evidence/explore.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps({k:v for k,v in out.items() if k not in ('rows_data','schema','last_rows')},ensure_ascii=True))
for x in schema: print(json.dumps(x,ensure_ascii=True))
print('LAST_ROWS',json.dumps(out['last_rows'],ensure_ascii=True))
