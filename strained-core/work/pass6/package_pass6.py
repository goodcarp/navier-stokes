from pathlib import Path
import hashlib,json,shutil,sys

src=Path(__file__).resolve().parent
base=src.parents[1]/'outputs'/'euler-ns-transfer'
dst=base/'pass6'
checks=dst/'checks'
checks.mkdir(parents=True,exist_ok=True)
for p in src.glob('*.md'):
    text=p.read_text().replace('python3 work/pass6/','python3 checks/')
    (dst/p.name).write_text(text)
for p in src.glob('*.py'):
    if p.name!='package_pass6.py':
        shutil.copy2(p,checks/p.name)
for p in src.glob('*.json'):
    shutil.copy2(p,checks/p.name)
(checks/'run_checks.py').write_text('''#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
here=Path(__file__).resolve().parent
programs=sorted(here.glob('verify_*.py'))
for p in programs:
    print('CHECK '+p.name,flush=True)
    r=subprocess.run([sys.executable,str(p)],cwd=here,capture_output=True,text=True)
    print(r.stdout,end='',flush=True)
    if r.stderr: print(r.stderr,end='',flush=True)
    if r.returncode: raise SystemExit(r.returncode)
print(f'PASS: {len(programs)} finite checks. No return or NS blowup certified.',flush=True)
''')
(dst/'REPORT.md').write_text(r'''# Sixth pass: reject the axisymmetric return and quantify three-dimensional transfer

8 September 2026. **The proposed indefinitely repeating axisymmetric rotating-core class is ruled out. A concrete three-dimensional initial-data family retains favorable local feedback and supplies positive mean torque, but a usable repeated stage remains unproved.**

## The obstruction changes the next task

The actual axisymmetric circulation maximum principle gives

\[
 G(V_{n+1})\le q_n^\alpha G(V_n),\qquad
 G(V)=\|rV_\theta\|_\infty .
\]

Under indefinitely shrinking scales and fixed \(\alpha>0\), this tends to zero. The pass5 return class's positive normalized rotation and uniform local profile bounds instead force \(G\ge g_*>0\). The [circulation proof](circulation-return-obstruction.md) derives the contradiction, a finite return-count bound, and a global ceiling for the rotational receiving Reynolds observable. It uses [Lei–Zhang's established scalar maximum principle](https://msp.org/pjm/2017/289-1/pjm-v289-n1-p06-s.pdf). Finite smooth axisymmetric forcing with a finite weighted torque budget cannot remove this obstruction.

Adding a bounded cubic core coordinate to the old class is therefore insufficient. This is a stronger restriction than the first-order shape defects found in pass5. The earlier finite local certificate and receiving gain remain valid.

The [escape audit](escape-route-audit.md) shows that allowing rotation to vanish does not rescue a class with uniform global C² bounds over complete normalized stages: its meridional profiles converge locally to zero. This conclusion needs those stronger global assumptions. At the critical exponent \(\alpha=0\), a uniform full-stage velocity bound implies the axisymmetric Type I bound excluded by [Chen–Strain–Tsai–Yau, Theorem 1.1](https://arxiv.org/pdf/0709.4230). Neither statement proves regularity for arbitrary axisymmetric flows. A genuine second scale or loss of those uniform bounds remains a separate research possibility.

## A concrete three-dimensional initial transfer

The [full angular-momentum equation](nonaxisymmetric-torque-budget.md) identifies the terms absent under axisymmetry. After azimuthal averaging, direct pressure torque cancels; transfer into mean rotation comes from actual quadratic velocity correlations.

An explicit compact solenoidal fourfold perturbation supplies strictly positive initial torque into a fixed receiver. The [receiver calculation](rotational-receiver-torque.md) proves the exact sign and frequency-independent bound

\[
 |Q_L(v)|\le C_\chi\|v\|_2^2.
\]

Thus a required receiver increment \(D\) from this stress term must pay
\(\int_0^\tau\|v\|_2^2dt\ge2I_LD/C_\chi\).
Increasing carrier frequency cannot provide unlimited mean transfer at fixed fluctuation energy.

The [three-dimensional embedding](three-dimensional-initial-gate.md) retunes the remote swirl exactly and gives an explicit norm-form interval of nonzero perturbation amplitudes for which

\[
 p_{zz}'(0)+32\le-\frac{290649}{17500},\qquad
 (b/\Omega)''(0)\ge\frac{290649}{35000}>0.
\]

These estimates use the complete three-dimensional pressure. The [functional audit](three-dimensional-gate-audit.md) verifies the Sobolev constants and retuning requirement, and the [full-family audit](nonaxisymmetric-family-audit.md) checks the actual seed, central symmetry, pressure tuning and receiver transfer together. No numerical useful-size perturbation amplitude or duration has been certified.

## The remaining decisive test

Initial compatibility and positive local torque do not prove enough sustained torque for a return. The compact receiver seed can initially draw on the existing outer circulation budget; it does not immediately increase that global mean maximum.

A separate [outer-supported family](outer-mean-maximum-target.md) addresses that distinction. At seed amplitudes \(3/4\le\lambda\le1\), exact central-pressure retuning remains possible, and a cutoff interval certificate proves an actual initial mean-maximum growth rate exceeding \(626/35\), including viscosity. **The full central pressure-derivative gate has not been checked at those finite amplitudes.** This is distinct from the small-amplitude family that retains the gate.

The [remaining pressure test](outer-pressure-compatibility-test.md) reduces the full calculation to two new interaction coefficients. At the already torque-admissible amplitude \(\lambda=3/4\), a sufficient condition is

\[
 D_\nu+A_{3/4}E < \frac{5093339}{210000}\approx24.254 .
\]

Neither interaction coefficient has been evaluated or enclosed. The [independent compatibility audit](outer-compatibility-independent-audit.md) checks the expansion and inequality directions. This is the immediate finite-amplitude compatibility test, not an asserted inequality.

The next target is one actual three-dimensional stage with a common amplitude/time window that simultaneously preserves the full feedback margin, supplies the required evolved stress, pays viscosity, and enters a controlled inherited successor state. The endpoint includes the old field, new angular modes, pressure tails, core deformation and outer geometry. No fresh packet or affine-core reset is allowed.

The previous pass's twelve checks and interval certificate remain unchanged. Run python3 checks/run_checks.py for this pass's finite identities. [Verification](VERIFICATION.md) and [sources](SOURCES.json) state the scope. No same-solution return, infinite compatible sequence, singular solution, or formal proof-kernel verification has been obtained.
''')
(dst/'VERIFICATION.md').write_text('''# Verification and limits

Run python3 checks/run_checks.py from this pass directory. Python with SymPy and mpmath is required. CHECK_RESULTS.txt records the packaged execution. The outer mean-maximum check imports the unchanged interval engine from the sibling pass5/certificate directory and recomputes its 64-panel cutoff-jet bounds.

The finite algebra covers the circulation equation and scaling, receiver lower bounds, forcing scaling, the conditional meridional recurrence and normalization, the fully three-dimensional angular equation, exact Fourier-mode covariance, selected-sign receiver torque, its frequency-independent energy bound, the pressure-functional Lipschitz constants, exact neutral retuning, and discrete-symmetry central identities.

The proofs additionally use scalar comparison, the earlier explicitly bounded Sobolev estimates, compactness and Liouville under the stated global bounds, smooth local existence/uniqueness, and the inherited pass5 interval certificate. These are analytical dependencies, not statements mechanically verified by the algebra programs. The maximum-principle equations and the Type I theorem were checked against the primary papers linked in the notes.

This pass supplies rigorous obstructions to specified normalized axisymmetric classes and a norm-form family of genuinely three-dimensional favorable initial data with positive local transfer. It does not enclose a useful perturbation amplitude, sustained torque, full terminal-state error, invariant return class, or blowup. The weighted/local profile assumptions and global C2 assumptions are kept distinct. A published conditional regularity theorem is not a general regularity claim.

The outer finite-amplitude family has a separate initial mean-maximum certificate, with viscosity and exact central pressure tuning. Its full pressure derivative is untested. The amplitude-polynomial checker reduces that pending test without evaluating the two new interaction coefficients. The two certified initial properties must not be combined across the distinct small-core and finite-outer seed families.

There is no new NS time-stepping computation or formal proof-kernel replay. The earlier pass2 through pass5 artifacts and their hashes are retained unchanged. Hashes identify files; finite checker exits do not prove the missing dynamical construction.
''')
source={
 'date':'2026-09-08',
 'inherited_manifests':[{'path':'../'+str(p.relative_to(base)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [base/'SOURCES.json',base/'pass5'/'SOURCES.json']],
 'primary_sources':[
 {'title':'Lei and Zhang, Criticality of the axially symmetric Navier–Stokes equations (PJM 2017)','url':'https://msp.org/pjm/2017/289-1/pjm-v289-n1-p06-s.pdf','checked':'Page 171, equations (1-3)–(1-4); scalar circulation maximum principle; also xi equation.'},
 {'title':'Chen, Strain, Tsai and Yau, Lower bounds on the blow-up rate of the axisymmetric Navier–Stokes equations II','url':'https://arxiv.org/pdf/0709.4230','checked':'Theorem 1.1, PDF page 2; full Type I velocity bound, initial H1/2 and bounded r utheta, spacetime pressure L5/3.'}],
 'new_derivations':'Exact normalized circulation obstruction, conditional meridional recurrence, compact Fourier seed and receiver torque budget, quantitative C4-symmetric embedding of the pass5 initial gate.',
 'verification':{'python':sys.version,'formal_kernel_replay':False,'NS_time_evolution':False,'same_solution_return_proved':False,'NS_blowup_proved':False}}
(dst/'SOURCES.json').write_text(json.dumps(source,indent=2)+'\n')
print(dst)
