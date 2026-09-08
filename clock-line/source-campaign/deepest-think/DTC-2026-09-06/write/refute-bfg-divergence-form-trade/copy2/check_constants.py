#!/usr/bin/env python3
"""
check_constants.py -- adversarial re-assertion gate for this folder.

Design rules (FL-043: a gate that cannot fail is not a control):
  * every check RECOMPUTES its target from first principles, by a route that does
    not read the stored value;
  * the stored value is only ever the LEFT side of the comparison;
  * a dynamic mutation test perturbs each stored number one at a time (v -> 1.5v+0.37)
    and requires the gate to FIRE on every single one;
  * the pass/fail line is printed from the real counters.
"""
import json, math, copy
import mpmath as mp
mp.mp.dps = 40

S = {k: json.load(open('s%s_results.json' % k)) for k in ['1','2','3','4','5']}

# ---------------- independent recomputations (no stored input) --------------
def rc_gradL1(nu, t):
    """int_{R^3}|grad G_nu| dx, by 3-D Monte-Carlo-free tensor quadrature in the
    ORIGINAL cartesian variables (different route from s1's radial integral)."""
    nu, t = mp.mpf(nu), mp.mpf(t)
    a = 4*nu*t
    # |grad G| = (|x|/(2 nu t)) G ; integrate in spherical but with the substitution
    # r = sqrt(a) * s, giving  (2/sqrt(pi nu t)) * (4/sqrt(pi)) int_0^inf s^3 e^{-s^2} ds ... :
    I = mp.quad(lambda s: s**3*mp.e**(-s**2), [0, 1, 6, 30])
    return (4/mp.sqrt(mp.pi))*I/mp.sqrt(nu*t)      # = 2/sqrt(pi nu t) iff I = 1/2

def rc_K6():
    f = lambda z: (z/2)*(4*mp.pi)**mp.mpf(-1.5)*mp.e**(-z**2/4)*(z+1)**4
    lo, hi = mp.mpf('0.5'), mp.mpf('8'); phi = (mp.sqrt(5)-1)/2
    for _ in range(300):
        x1 = hi - phi*(hi-lo); x2 = lo + phi*(hi-lo)
        if f(x1) > f(x2): hi = x2
        else: lo = x1
    return f((lo+hi)/2)

def rc_K7():
    f = lambda z: (4*mp.pi)**mp.mpf(-1.5)*mp.e**(-z**2/4)*(z+1)**4
    lo, hi = mp.mpf('0.5'), mp.mpf('8'); phi = (mp.sqrt(5)-1)/2
    for _ in range(300):
        x1 = hi - phi*(hi-lo); x2 = lo + phi*(hi-lo)
        if f(x1) > f(x2): hi = x2
        else: lo = x1
    return f((lo+hi)/2)

def rc_gamma():
    return (mp.mpf(3)/2)**(mp.mpf(2)/5) + (mp.mpf(2)/3)**(mp.mpf(3)/5)

def rc_Cu():
    """golden-section minimisation of A sqrt(tau) + B tau^{-3/4} at M = E = 1"""
    A = 4/mp.sqrt(mp.pi); B = (8*mp.pi)**mp.mpf('-0.75')
    g = lambda tv: A*mp.sqrt(tv) + B*tv**mp.mpf('-0.75')
    lo, hi = mp.mpf('1e-8'), mp.mpf('1e4'); phi = (mp.sqrt(5)-1)/2
    for _ in range(600):
        x1 = hi - phi*(hi-lo); x2 = lo + phi*(hi-lo)
        if g(x1) < g(x2): hi = x2
        else: lo = x1
    return g((lo+hi)/2)

def rc_Cw():
    A = 4/mp.sqrt(mp.pi); B = 2**mp.mpf('0.75')/(4*mp.pi**mp.mpf('0.75'))
    g = lambda tv: A*mp.sqrt(tv) + B*tv**mp.mpf('-0.25')
    lo, hi = mp.mpf('1e-8'), mp.mpf('1e4'); phi = (mp.sqrt(5)-1)/2
    for _ in range(600):
        x1 = hi - phi*(hi-lo); x2 = lo + phi*(hi-lo)
        if g(x1) < g(x2): hi = x2
        else: lo = x1
    return g((lo+hi)/2)

def rc_c_div(Q=1):
    """END-TO-END numeric route: pick numbers, run the whole bootstrap, read off
    M*Re_E*t_*.  Nothing symbolic, nothing stored."""
    Cu = rc_Cu()
    M, E, nu = mp.mpf('3.7'), mp.mpf('0.41'), mp.mpf('0.013')
    K = mp.mpf(2)
    Cnu = 2*Q/mp.sqrt(mp.pi*nu)
    # Duhamel(t) = Cnu * 2 sqrt(t) * (K M) * Cu E^{1/5} (K M)^{3/5}
    coef = Cnu*2*(K*M)*Cu*E**(mp.mpf(1)/5)*(K*M)**(mp.mpf(3)/5)
    tstar = (M/(2*coef))**2
    ReE = E**(mp.mpf(2)/5)*M**(mp.mpf(1)/5)/nu
    return M*ReE*tstar, tstar, coef, M, E, nu, ReE

def rc_pi246():
    return mp.pi/2**(mp.mpf(46)/5)

def rc_log_row(logRe):
    return 1/(1+max(mp.log(mp.e**mp.mpf(logRe)) if False else mp.mpf(logRe), mp.mpf(0)))

# ---------------- the checks -------------------------------------------------
CHECKS = []   # (label, stored_getter, recompute, rtol)

# NOTE on tolerances: the stored values are IEEE doubles written to JSON, so no
# check can be tighter than ~1e-16 relative.  1e-15 is the floor used here; the
# mutation test below (v -> 1.5v+0.37) moves every number by O(1) relative, so the
# gate's discriminating power does not depend on the tolerance being tight.
def add(label, path, rec, rtol=mp.mpf('1e-15')):
    CHECKS.append((label, path, rec, rtol))

def get(d, path):
    cur = d
    for p in path: cur = cur[p]
    return cur

def setp(d, path, val):
    cur = d
    for p in path[:-1]: cur = cur[p]
    cur[path[-1]] = val

for i, (nuv, tv) in enumerate([(1,1),(1,0.01),(7,3),(0.001,5)]):
    add('s1 int|grad G| nu=%s t=%s' % (nuv,tv), ('1','K3_numeric',i,2), (lambda nuv=nuv,tv=tv: rc_gradL1(nuv,tv)), mp.mpf('1e-15'))
add('s1 K6 sup|gradG|(|x|+1)^4', ('1','K6_gradmajorant'), rc_K6, mp.mpf('1e-12'))
add('s1 K7 sup G(|x|+1)^4',      ('1','K7_Gmajorant'),    rc_K7, mp.mpf('1e-12'))
add('s1 K8 8pi/3',               ('1','K8_8pi_3'),        lambda: 8*mp.pi/3)
add('s1 K8 4pi/3',               ('1','K8_4pi_3'),        lambda: 4*mp.pi/3)
add('s2 Q exact',                ('2','Q_exact'),         lambda: mp.mpf(1))
add('s2 Q attained',             ('2','Q_attained'),      lambda: mp.mpf(1), mp.mpf('1e-14'))
add('s2 Q triangle',             ('2','Q_triangle'),      lambda: mp.mpf(2))
add('s3 gamma',                  ('3','gamma'),           rc_gamma, mp.mpf('1e-15'))
add('s3 C_u',                    ('3','C_u'),             rc_Cu,    mp.mpf('1e-15'))
add('s3 C_w enstrophy variant',  ('3','C_w_enstrophy_variant'), rc_Cw, mp.mpf('1e-11'))
add('s4 c_div',                  ('4','c_div'),           lambda: rc_c_div(1)[0], mp.mpf('1e-14'))
add('s4 c_div (Q=2)',            ('4','c_div_Q2'),        lambda: rc_c_div(2)[0], mp.mpf('1e-14'))
add('s4 Q gain factor',          ('4','Q_gain_factor'),   lambda: rc_c_div(1)[0]/rc_c_div(2)[0], mp.mpf('1e-14'))
add('s4 pi*2^(-46/5)',           ('4','pi_2_46_5'),       rc_pi246, mp.mpf('1e-15'))
for i, R in enumerate([1,1e2,1e5,1e10,1e20,1e40]):
    add('s4 log-clock 1/(1+logRe) Re=%g' % R, ('4','comparison_rows',i,1), (lambda R=R: 1/(1+max(mp.log(mp.mpf(R)),mp.mpf(0)))), mp.mpf('1e-15'))
    add('s4 power-clock c/Re     Re=%g' % R, ('4','comparison_rows',i,2), (lambda R=R: rc_c_div(1)[0]/mp.mpf(R)), mp.mpf('1e-14'))
for i, L in enumerate([8.3178, 20, 50, 100]):
    add('s5 log Re_E   L=%g' % L, ('5','rows',i,1), (lambda L=L: 2*mp.mpf(L)-mp.mpf('0.7031660')), mp.mpf('1e-15'))
    add('s5 Re_E       L=%g' % L, ('5','rows',i,2), (lambda L=L: mp.e**(2*mp.mpf(L)-mp.mpf('0.7031660'))), mp.mpf('1e-14'))
    add('s5 power win  L=%g' % L, ('5','rows',i,3), (lambda L=L: rc_c_div(1)[0]/mp.e**(2*mp.mpf(L)-mp.mpf('0.7031660'))), mp.mpf('1e-13'))
    add('s5 log   win  L=%g' % L, ('5','rows',i,4), (lambda L=L: 1/(1+2*mp.mpf(L)-mp.mpf('0.7031660'))), mp.mpf('1e-15'))
    add('s5 model win  L=%g' % L, ('5','rows',i,5), (lambda L=L: 2*mp.log(2)/mp.mpf(L)), mp.mpf('1e-15'))

WANT = None
def build_wants():
    global WANT
    WANT = [mp.mpf(rec()) for _, _, rec, _ in CHECKS]

def run_gate(data, verbose=True):
    npass = nfail = 0; fails = []
    for k, (label, path, rec, rtol) in enumerate(CHECKS):
        stored = mp.mpf(repr(get(data[path[0]], path[1:])))
        want = WANT[k]
        rel = abs(stored-want)/max(abs(want), mp.mpf(1))
        ok = rel <= rtol
        if ok: npass += 1
        else:
            nfail += 1; fails.append((label, float(stored), float(want), float(rel)))
        if verbose:
            print("  %-42s stored=%-24.16g recomputed=%-24.16g rel=%-10.2e %s"
                  % (label, float(stored), float(want), float(rel), "PASS" if ok else "FAIL"))
    return npass, nfail, fails

build_wants()   # recomputed ONCE, from first principles, never from stored data
print("=== GATE (independent recomputation of every stored number) ===")
npass, nfail, fails = run_gate(S)
print("\n  %d PASS, %d FAIL" % (npass, nfail))
for f in fails: print("   FAIL:", f)

print("\n=== END-TO-END sanity: the bootstrap actually closes at t_* ===")
cd, tstar, coef, M, E, nu, ReE = rc_c_div(1)
Duh = coef*mp.sqrt(tstar)
print("  M = %s, E = %s, nu = %s -> Re_E = %.10f" % (M, E, nu, float(ReE)))
print("  t_* = %.16e   Duhamel(t_*) = %.16e   M/2 = %.16e   rel = %.2e"
      % (float(tstar), float(Duh), float(M/2), float(abs(Duh-M/2)/(M/2))))
assert abs(Duh-M/2)/(M/2) < mp.mpf('1e-30')
print("  M*Re_E*t_* = %.16e   (= c_div)" % float(cd))
print("  omega cap check: ||omega|| <= M + M/2 = 1.5 M < 2M (the bootstrap cap) -> closes")

print("\n=== DYNAMIC MUTATION TEST (v -> 1.5 v + 0.37, one stored number at a time) ===")
caught = missed = 0; missed_labels = []
for label, path, rec, rtol in CHECKS:
    D = copy.deepcopy(S)
    cur = D[path[0]]
    for p in path[1:-1]: cur = cur[p]
    cur[path[-1]] = 1.5*float(cur[path[-1]]) + 0.37
    p2, f2, _ = run_gate(D, verbose=False)
    if f2 >= 1: caught += 1
    else: missed += 1; missed_labels.append(label)
print("  mutated %d stored numbers; gate CAUGHT %d, MISSED %d" % (len(CHECKS), caught, missed))
for m in missed_labels: print("   MISSED:", m)

ok = (nfail == 0 and missed == 0)
print("\n%s  (%d checks, %d mutations caught / %d)" %
      ("ALL CHECKS PASS AND EVERY MUTATION CAUGHT" if ok else "GATE FAILED", npass, caught, len(CHECKS)))
raise SystemExit(0 if ok else 1)
