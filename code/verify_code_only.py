#!/usr/bin/env python3
"""MIT. Verify final code-only inventory, source identity, syntax and publication boundary."""
from pathlib import Path
import argparse,ast,csv,hashlib,json,re
ROOT=Path(__file__).resolve().parents[1]
FORBIDDEN={".pdf",".png",".svg",".jpg",".jpeg",".tif",".tiff",".xlsx",".xls",".docx",".npz",".npy",".rds",".rda",".rdata",".parquet",".pkl",".pickle",".joblib",".zip",".bundle",".csv"}
def verify():
    rows=list(csv.DictReader((ROOT/"manifest_sha256.tsv").open(encoding="utf8"),delimiter="\t"))
    expected={r["file"] for r in rows}
    actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file()}
    if len(rows)!=len(expected) or actual!=expected|{"manifest_sha256.tsv"}:raise ValueError("Manifest file set mismatch")
    py=0
    for r in rows:
        p=ROOT/r["file"]
        if p.is_symlink() or not p.resolve().is_relative_to(ROOT):raise ValueError("Unsafe file path")
        b=p.read_bytes()
        if len(b)!=int(r["bytes"]) or hashlib.sha256(b).hexdigest()!=r["sha256"]:raise ValueError("Hash mismatch: "+r["file"])
        if p.suffix.lower() in FORBIDDEN:raise ValueError("Withheld payload in package")
        s=b.decode("utf-8-sig")
        if p.suffix==".py":ast.parse(s,filename=r["file"]);py+=1
        if re.search(r"(?<![A-Za-z])[A-Za-z]:[\\/](?!/)",s) or any(x in s for x in ["data"+":image/","iVBOR"+"w0KGgo","JVBE"+"Ri0"]):raise ValueError("Private path or embedded image: "+r["file"])
        if re.search(r"gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9]{20,}|-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----",s):raise ValueError("Credential-like content")
        if r["file"].startswith("pipeline/") and ("_private_literal(" in s or ".private_literal(" in s or "EJA_PRIVATE_LITERALS" in s):raise ValueError("Obsolete literal-sidecar dependency")
    sources=list(csv.DictReader((ROOT/"provenance/source_code_map.tsv").open(encoding="utf8"),delimiter="\t"))
    for r in sources:
        if hashlib.sha256((ROOT/r["file"]).read_bytes()).hexdigest()!=r["final_sha256"]:raise ValueError("Source identity mismatch")
    return {"status":"PASS","scope":"inventory/hash/source-identity/Python syntax only",
            "files":len(rows)+1,"python_parsed":py,"scientific_sources":len(sources),
            "byte_identical_originals":sum(r["status"]=="BYTE_IDENTICAL_ORIGINAL" for r in sources),
            "scientific_reproduction":"NOT_VERIFIED","training":"NOT_RUN","figures":"NOT_RUN"}
def main():
    p=argparse.ArgumentParser(description=__doc__);p.parse_args();print(json.dumps(verify(),indent=2))
if __name__=="__main__":main()
