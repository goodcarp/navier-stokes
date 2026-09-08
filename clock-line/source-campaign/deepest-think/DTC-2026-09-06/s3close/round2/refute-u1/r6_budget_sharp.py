import sys, math
sys.path.insert(0, "/private/tmp/claude-501/-Users-spaceman-Desktop/e34d31b2-1bd6-4688-bd1c-674bd591026e/scratchpad/refute-u1")
import u5_budget as B
CST = B.CSTAR
FS = [0.05,0.1,0.25,0.5,1.0,2.0,4.0,8.0,16.0,32.0]
def Lstar(col, target):
    def eps(L): return min(B.assemble(L, CST, col, f=f)["eps"] for f in FS)
    lo, hi = 2.0, 1e9
    if eps(hi) > target: return None
    for _ in range(45):
        mid = math.sqrt(lo*hi)
        if eps(mid) <= target: hi = mid
        else: lo = mid
        if hi/lo < 1+1e-7: break
    return hi
for nm, cr in [("u1 headline  phi0=30    C_R=29.0242", 29.02416532565612),
               ("sharp        phi0=30    C_R=17.4007", 17.400705),
               ("sharp cone   sup>=delta C_R=40.2293", 40.229279),
               ("sharp+taper  sup ALL    C_R=100.6122", 100.612238)]:
    B.C_R_PROVED = cr
    a = Lstar("proved", 0.5)
    print("%-40s L_*(proved,eps<=1/2) = %10.1f   log Lambda_* = %11.1f" % (nm, a, 2*a+B.LOGREE_SHIFT))
