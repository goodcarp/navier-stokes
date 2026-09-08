#!/usr/bin/env python3
"""Falsifiable prediction for the N=7 quadrupling time, from the N=4,5,6 fits only.
Run again against dep_n7_partial.json / dep_n7.json when the in-flight N=7 run crosses 4 M0."""
import os, json, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))

runs = json.load(open(os.path.join(HERE, "dep_main.json")))["runs"]
P = {r["N"]: r for r in runs}
use = [4, 5, 6]
x = np.array([P[N]["logReE"] for N in use])
y4 = np.array([1.0 / P[N]["cross"]["T4"]["t"] for N in use])
A = np.vstack([x, np.ones(len(x))]).T
cf, *_ = np.linalg.lstsq(A, y4, rcond=None)
n7 = json.load(open(os.path.join(HERE, "dep_n7_partial.json")))["runs"][0]
lr7 = n7["logReE"]
pred_fit = 1.0 / (cf[0] * lr7 + cf[1])

# independent route: the ratio T4/T32, linear in 1/N over N=4,5,6, evaluated at N=7
rat = np.array([P[N]["cross"]["T4"]["t"] / P[N]["cross"]["T32"]["t"] for N in use])
inv = 1.0 / np.array(use, dtype=float)
B = np.vstack([inv, np.ones(3)]).T
cr, *_ = np.linalg.lstsq(B, rat, rcond=None)
rat7 = cr[0] / 7.0 + cr[1]
pred_rat = rat7 * n7["cross"]["T32"]["t"]

print(f"N=7: log Re_E = {lr7:.4f}, measured T32 = {n7['cross']['T32']['t']:.4f}, "
      f"T2 = {n7['cross'].get('T2',{}).get('t')}")
print(f"  route A  1/T4 = A logReE + B fitted on N=4,5,6 :  T4(N=7) = {pred_fit:.4f}")
print(f"  route B  T4/T32 linear in 1/N on N=4,5,6       :  ratio(7) = {rat7:.4f}, "
      f"T4(N=7) = {pred_rat:.4f}")
print(f"  PREDICTION: T4(N=7)*M0 = {0.5*(pred_fit+pred_rat):.4f} "
      f"(two routes differ by {abs(pred_fit-pred_rat):.4f})")
print(f"  DEPLETION would instead require T4(N=7) ~ T4(N=6) = "
      f"{P[6]['cross']['T4']['t']:.4f} or larger.")
t4 = n7["cross"].get("T4", {}).get("t")
print("  in-flight N=7 measured T4 =", t4 if t4 else
      f"NOT YET REACHED (run at t={n7['t_final']:.4f}, ||omega||={n7['ommax_final']:.4f})")
