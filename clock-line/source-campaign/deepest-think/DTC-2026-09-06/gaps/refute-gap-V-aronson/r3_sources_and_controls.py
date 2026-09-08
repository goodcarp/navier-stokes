#!/usr/bin/env python3
"""
REFUTER r3 -- (a) what the files in gap-V-aronson/sources/ ACTUALLY are, and
                  whether the three verbatim quotes in section 3 are at source;
              (b) whether the pre-declared refuter rules R1 and R3 have any failure mode.
"""
import json, os, hashlib, re
G = "~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/gaps/gap-V-aronson"
OUT = {}

# ---------- (a1) classify every file the note lists as a source ----------
print("(a1) sources/ file audit")
rows=[]
for fn in sorted(os.listdir(os.path.join(G,'sources'))):
    p=os.path.join(G,'sources',fn)
    if os.path.isdir(p): continue
    b=open(p,'rb').read()
    if b[:5]==b'%PDF-':      kind='real PDF'
    elif b.lstrip()[:9].lower()==b'<!doctype' or b.lstrip()[:5].lower()==b'<html': kind='HTML PAGE (not a paper)'
    else:                    kind='text'
    marker=''
    if kind.startswith('HTML'):
        t=re.sub(r'<[^>]+>',' ',b.decode('utf-8','replace'))
        t=re.sub(r'\s+',' ',t)
        m=re.search(r"No document for '([^']+)'",t)
        if m: marker="arXiv: no document for %s"%m.group(1)
        elif '404' in t[:400]: marker="404 page (%d bytes)"%len(b)
        else: marker=t[:90].strip()
    rows.append(dict(file=fn,bytes=len(b),kind=kind,marker=marker,
                     sha256=hashlib.sha256(b).hexdigest()[:16]))
    print("   %-26s %9d  %-24s %s" % (fn,len(b),kind,marker))
OUT['sources_audit']=rows

# ---------- (a2) the three section-3 quotes, checked against the REAL survey ----------
real = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "aronson_survey_1707.04620v1.txt")).read()
flat = re.sub(r'\s+',' ',real)
quotes = {
 "FS  (used to dismiss Fabes-Stroock)":
   "Neither Nash nor Fabes & Stroock consider the full equation (1)",
 "NS  (used to dismiss Norris-Stroock)":
   "they are forced to assume the uniform continuity of A and E",
 "PE  (used to dismiss Porper-Eidel'man)":
   "a slight generalization of equations (1) and",
}
print("\n(a2) verbatim check of the three section-3 quotes against arXiv:1707.04620v1")
qres={}
for k,q in quotes.items():
    hit = q in flat
    near = None
    if not hit:
        # nearest variant
        for cand in ["Neither Nash or Fabes & Stroock consider the full equation (1)"]:
            if cand in flat: near=cand
    qres[k]=dict(quoted=q,verbatim_at_source=hit,source_reads=near)
    print("   %-40s verbatim=%-5s %s" % (k,hit, "" if hit else "SOURCE READS: '%s'"%near))
OUT['quote_check']=qres
OUT['survey_retrieved_by_this_refuter']=True
OUT['survey_present_in_gapV_sources']=False

# ---------- (b) do R1 and R3 have a failure mode? ----------
print("\n(b) failure-mode audit of the pre-declared refuter rules")
v4=json.load(open(os.path.join(G,'v4_results.json')))
bnds=[r['bound'] for r in v4['A_rows']]
print("   R1: rule = 'refuted if P_emp - 3se > min(B1,B2)'.  P_emp is a probability, so P_emp <= 1.")
print("       min(B1,B2) over the five TEST A rows = %s" % ["%.3f"%b for b in bnds])
print("       min over all rows = %.4f  -> %s"
      % (min(bnds), "ALL BOUNDS EXCEED 1: R1 has NO failure mode" if min(bnds)>1 else "R1 could have fired"))
OUT['R1_min_bound']=float(min(bnds)); OUT['R1_vacuous']=bool(min(bnds)>1)

drift=abs(v4['B_Xtau']-v4['B_x0'])
ds=[r['frac']*(r['nu']*0.5)**0.5 for r in v4['B_rows']]
print("   R3: rule = 'non-diagnostic unless the Eulerian control violates the bound'.")
print("       drift displacement of x0 over the window = %.4f ; the datum hole has radius d <= %.4f"
      % (drift, max(ds)))
print("       the hole is advected clear of x0 by %.1fx its own radius -> the control is FORCED"
      % (drift/max(ds)))
print("       -> R3 also has no failure mode as constructed (it is a demonstration, not a test).")
OUT['R3_drift']=float(drift); OUT['R3_dmax']=float(max(ds)); OUT['R3_ratio']=float(drift/max(ds))

bb=[r['bound'] for r in v4['B_rows']]; dv=[r['dev_material'] for r in v4['B_rows']]
print("   R2: bounds %s ; measurements %s -> R2 CAN fire (margin %.1fx at the tightest row)"
      % (["%.3f"%x for x in bb], ["%.2e"%x for x in dv], min(b/d for b,d in zip(bb,dv))))
OUT['R2_can_fire']=True; OUT['R2_tightest_margin']=float(min(b/d for b,d in zip(bb,dv)))

json.dump(OUT,open("r3_results.json","w"),indent=1)
print("\nwrote r3_results.json")
