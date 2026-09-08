#!/usr/bin/env python3
"""Exact dyadic replay of the new leading full-pressure initial-jet certificate."""
from pathlib import Path
from fractions import Fraction as F
import json
here=Path(__file__).resolve().parent
data=json.loads((here/"general-envelope-initial-jet-certificate.json").read_text())
def ep(key,side):
    x=data[key]["exact"][side]
    return (-1 if x["sign"] else 1)*F(int(x["mantissa"]))*F(2)**x["exponent"]
assert data["status"]=="PASS"
assert all(data["exact_rational_predicates"].values())
assert data["parameters"]["k"]==-20
assert data["parameters"]["radial_center"]=="7/50"
assert data["parameters"]["radial_width"]=="1/10"
assert ep("Gt_actual","lower")>11 and ep("Gt_actual","upper")<106
assert ep("Gtt_actual_full_pressure","lower")>-101000
assert ep("Gtt_actual_full_pressure","upper")<-13000
assert ep("moving_radial_value_Gtt_actual","lower")>-101000
assert ep("moving_radial_value_Gtt_actual","upper")<-12000
assert ep("NEW_full_pressure_radial_error_bound","upper")<F(13,10)
assert ep("NEW_weighted_residual_norm","upper")<F(855612,100000)
assert ep("S_actual","upper")<0
assert ep("inner_pressure_tail","lower")>=0
assert ep("outer_pressure_tail","lower")>=0
assert "upper-bound expression" in data["pressure_error_interval_semantics"]
assert "lower endpoint is not a lower bound" in data["pressure_error_interval_semantics"]
for relative in data["dependencies"].values():
    assert (here/relative).is_file(),relative
print("PASS: new leading data, complete residual/tails and dependencies, "
      "Gt>11, actual Gtt<-13000, tracked Gtt<-12000, and new pressure error<1.3.")
