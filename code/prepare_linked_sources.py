#!/usr/bin/env python3
"""Public source-acquisition adapter; SPDX-License-Identifier: MIT.
Reads author-acquired authentic attachments; never downloads or redistributes.
--check-only changes nothing. --write-workspace creates a NEW private workspace,
adds the accepted attachments and the minimum S1 paragraph probe for old audit.py.
No research script, model, or bootstrap is executed.
"""
from pathlib import Path
import argparse,csv,json,hashlib,zipfile,xml.etree.ElementTree as ET,shutil
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def extract_s1(p):
 with zipfile.ZipFile(p) as z:
  if z.testzip() is not None:raise ValueError('DOCX ZIP CRC failure')
  body=ET.fromstring(z.read('word/document.xml'))
 ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
 return [''.join(n.text or '' for n in q.findall('.//w:t',ns)).strip() for q in body.findall('.//w:p',ns)]
def validate(source_root):
 manifest=list(csv.DictReader((ROOT/'provenance/linked_sources.tsv').open(encoding='utf8'),delimiter='\t'));out=[]
 for r in manifest:
  rel=Path(r['expected_relative_path']);p=(source_root/rel).resolve()
  assert p.is_relative_to(source_root) and p.is_file(),r['expected_relative_path']
  assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['expected_relative_path']
  out.append(r)
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument('--sources-root',type=Path,required=True);g=p.add_mutually_exclusive_group(required=True);g.add_argument('--check-only',action='store_true');g.add_argument('--write-workspace',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
 report=a.report.resolve() if a.report else None
 if report and report.is_relative_to(ROOT):raise ValueError('Report must be outside sealed release')
 source=a.sources_root.resolve();records=validate(source);s1=next(r for r in records if r['expected_relative_path'].endswith('Sheet_S1.docx'));paragraphs=extract_s1(source/s1['expected_relative_path'])
 result={'status':'PASS','source_files_checked':len(records),'source_sha256_checked':True,'s1_paragraphs':len(paragraphs),'scope':'attachment identity and deterministic reference-text extraction only','research_scripts_executed':False,'workspace_created':False}
 if a.write_workspace:
  dest=a.write_workspace.resolve()
  if dest.is_relative_to(ROOT) or dest.exists():raise ValueError('Workspace must be NEW and outside sealed release')
  if report and report.is_relative_to(dest):raise ValueError('Report must be outside workspace')
  shutil.copytree(ROOT,dest)
  for r in records:
   q=dest/r['expected_relative_path'];q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source/r['expected_relative_path'],q)
  target=dest/'pipeline/Schutz_2018_phase_minus2/evidence/s1_probe.json'
  target.write_text(json.dumps({'_public_derivation':{'source_sha256':s1['sha256'],'scope':'minimum paragraphs required by historical pairing audit, no forensic metadata claims'},'paragraphs':paragraphs},ensure_ascii=False,indent=2),encoding='utf8')
  result['workspace_created']=True
 if report:report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(result,indent=2),encoding='utf8')
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
