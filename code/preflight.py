#!/usr/bin/env python3
"""MIT. Read-only input-presence and dependency-discovery check, without loading scientific data."""
from pathlib import Path
import argparse,json,importlib.util as iu,importlib.metadata as md
ROOT=Path(__file__).resolve().parents[1]
DEPS={"numpy":"numpy","pandas":"pandas","scipy":"scipy","sklearn":"scikit-learn","xgboost":"xgboost","joblib":"joblib","threadpoolctl":"threadpoolctl","openpyxl":"openpyxl","docx":"python-docx"}
def inspect(inputs_root,stage=None):
    c=json.loads((ROOT/"config/input_contracts.json").read_text(encoding="utf8"))
    scripts=[r for r in c["scripts"] if stage is None or r["stage"]==stage]
    if not scripts:raise ValueError("Unknown stage")
    paths=sorted({p for r in scripts for p in r["required_paths_relative_to_pipeline"]})
    absent=[p for p in paths if not (inputs_root/p).exists()]
    deps={}
    for module,dist in DEPS.items():
        present=iu.find_spec(module) is not None
        try:version=md.version(dist) if present else None
        except md.PackageNotFoundError:version=None
        deps[dist]={"discoverable":present,"installed_version":version}
    return {"status":"INPUTS_WITHHELD_OR_MISSING" if absent else "KEY_INPUTS_PRESENT_ONLY",
            "scripts_checked":len(scripts),"representative_input_paths":len(paths),"missing":absent,
            "python_dependencies":deps,"input_list_exhaustive":False,
            "data_values_read":False,"scientific_scripts_imported":False,"scientific_execution":"NOT_RUN"}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs-root",type=Path,default=ROOT/"pipeline")
    p.add_argument("--stage");p.add_argument("--require-inputs",action="store_true")
    a=p.parse_args()
    try:r=inspect(a.inputs_root.resolve(),a.stage)
    except ValueError as e:p.error(str(e))
    print(json.dumps(r,indent=2))
    if a.require_inputs and r["missing"]:raise SystemExit(2)
if __name__=="__main__":main()
