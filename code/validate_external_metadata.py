#!/usr/bin/env python3
"""MIT. Validate the format of an AUTHENTIC private metadata file; never fit or analyse."""
from pathlib import Path
import argparse,json,math
ROOT=Path(__file__).resolve().parents[1]
PREDICTORS={"crop","climate_zone","soil_texture","inoculant_function","irrigation_status","fertilization_type"}
def number(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def validate(d):
    required={"field_sites","mapped_predictors","potassium_product_rate_kg_ha","potassium_k2o_mass_fraction"}
    if not isinstance(d,dict) or set(d)!=required:raise ValueError("Wrong metadata fields")
    if not isinstance(d["field_sites"],list) or not d["field_sites"]:raise ValueError("Authentic sites required")
    seen=set()
    for r in d["field_sites"]:
        if not isinstance(r,dict) or set(r)!={"year","site","latitude","longitude"}:raise ValueError("Wrong site fields")
        if not isinstance(r["year"],int) or isinstance(r["year"],bool):raise ValueError("year must be integer")
        if not isinstance(r["site"],str) or not r["site"].strip():raise ValueError("site must be nonempty")
        if not number(r["latitude"]) or not -90<=r["latitude"]<=90:raise ValueError("Invalid latitude")
        if not number(r["longitude"]) or not -180<=r["longitude"]<=180:raise ValueError("Invalid longitude")
        key=(r["year"],r["site"])
        if key in seen:raise ValueError("Duplicate site/year")
        seen.add(key)
    m=d["mapped_predictors"]
    if not isinstance(m,dict) or set(m)!=PREDICTORS:raise ValueError("Wrong mapped predictor fields")
    for k,v in m.items():
        if k in {"inoculant_function","irrigation_status"}:
            if not isinstance(v,int) or isinstance(v,bool):raise ValueError("Nominal code must be integer")
        elif not isinstance(v,str) or not v.strip():raise ValueError("Mapped label must be nonempty")
    for k in ["potassium_product_rate_kg_ha","potassium_k2o_mass_fraction"]:
        if not number(d[k]) or d[k]<0:raise ValueError("Invalid product input")
    if d["potassium_k2o_mass_fraction"]>1:raise ValueError("Mass fraction exceeds 1")
    return {"status":"FORMAT_VALID","scope":"structure/types only; source authenticity and frozen identity NOT VERIFIED"}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group(required=True);g.add_argument("--schema",action="store_true");g.add_argument("--check",type=Path)
    a=p.parse_args()
    if a.schema:
        print((ROOT/"config/external_source_metadata.schema.json").read_text(encoding="utf8"));return
    if not a.check.is_file():p.error("Private metadata file missing")
    try:print(json.dumps(validate(json.loads(a.check.read_text(encoding="utf8"))),indent=2))
    except (ValueError,TypeError) as e:p.error(str(e))
if __name__=="__main__":main()
