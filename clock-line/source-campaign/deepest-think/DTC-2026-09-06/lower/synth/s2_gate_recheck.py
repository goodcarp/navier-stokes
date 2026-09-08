# Independent re-check of the two FL-043 gate failures reported by refuters.
import re, json, subprocess, os
out={}
p="~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/lower/prove-lagrangian/check_constants.py"
src=open(p).read()
out['star_zero_mask_lines']=[l.strip() for l in src.splitlines() if '*0+' in l or '*0 +' in l]
out['n_check_calls']=len(re.findall(r'^\s*ck\(', src, re.M))
m=re.findall(r'ALL %d CHECKS PASS.*', src)
out['summary_line']=[l.strip() for l in src.splitlines() if 'CHECKS PASS' in l]
log="~/Desktop/Solve Navier Stokes/campaign/deepest-think/DTC-2026-09-06/lower/prove-lagrangian/check_log.txt"
if os.path.exists(log):
    L=open(log).read()
    out['log_pass_lines']=L.count('PASS')
    out['log_fail_lines']=L.count('FAIL')
    out['log_tail']=L.strip().splitlines()[-3:]
print(json.dumps(out,indent=1))
json.dump(out,open('s2_results.json','w'),indent=1)
