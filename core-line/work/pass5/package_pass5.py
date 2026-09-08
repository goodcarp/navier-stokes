from pathlib import Path
import hashlib,json,shutil,sys

src=Path(__file__).resolve().parent
base=src.parents[1]/'outputs'/'euler-ns-transfer'
dst=base/'pass5'
dst.mkdir(parents=True,exist_ok=True)
cert=dst/'certificate'
cert.mkdir(exist_ok=True)
for path in src.glob('*.md'):
    text=path.read_text().replace('`work/pass3/quantitative-core-stage.md`','[the pass3 derivative estimates](../pass3/quantitative-core-stage.md)')
    (dst/path.name).write_text(text)
for path in src.glob('*.py'):
    if path.name not in {'package_pass5.py','evaluate_radial_pump.py'}:
        shutil.copy2(path,cert/path.name)
for name in ['radial-strain-interval-1024.json','core-twb-interval.json','local-interval-results.json','FULL_INITIAL_GATE_CERTIFICATE.json']:
    if name != 'FULL_INITIAL_GATE_CERTIFICATE.json' or not (cert/name).exists():
        shutil.copy2(src/name,cert/name)

(cert/'run_checks.py').write_text('''#!/usr/bin/env python3
"""Run the finite symbolic, interval-consistency and exact arithmetic checks."""
from pathlib import Path
import subprocess,sys
here=Path(__file__).resolve().parent
programs=sorted(here.glob('verify_*.py'))
for program in programs:
    print('CHECK '+program.name,flush=True)
    result=subprocess.run([sys.executable,str(program)],cwd=here,capture_output=True,text=True)
    print(result.stdout,end='',flush=True)
    if result.stderr:
        print(result.stderr,end='',flush=True)
    if result.returncode:
        raise SystemExit(result.returncode)
print(f'PASS: {len(programs)} finite checks. No return map or blowup is certified.',flush=True)
''')

(dst/'REPORT.md').write_text('''# Fifth pass: a certified initial feedback inequality

8 September 2026. **The full initial pressure inequality is now certified for an explicit compact smooth, unforced ordinary Navier–Stokes datum. Same-solution return and blowup remain open.**

Replacing the thin cylindrical pump with a radial annular strain makes the entire initial pressure-feedback test reduce exactly to angular degrees zero, two and four and bounded radial integrals. This reduction includes the nonlocal pressure response and all viscosity terms. It applies to this initial scalar functional; it is not a finite-dimensional evolution model.

For the exact datum with central strain and rotation both one, annular strength 7/5, and viscosity 1/1000, outward-rounded interval integration and exact analytical bounds give

\\[
p_{zz}'(0)+32\le-\\frac{290649}{8750}<-33,
\qquad (b/\\Omega)''(0)\ge\\frac{290649}{17500}>16.
\\]

The [main proof](CERTIFIED_INITIAL_FEEDBACK.md) specifies the datum, exact amplitude tuning and complete component identity. The [independent reconstruction audit](full-gate-certificate-audit.md) checks that no interaction or viscous contribution was omitted. Multiplying velocity by 1000 and accelerating time by 1000 converts the same construction to viscosity one.

The [finite-time corollary](finite-feedback-corollary.md) proves that one actual local solution has a positive interval with increasing strain/rotation ratio, strain growth faster than 2b², increasing normalized exterior swirl stress, and actual receiving-velocity/Reynolds gain. It gives a sufficient interval in terms of the new datum's H10 norm. That norm has not been numerically enclosed, so this package claims no decimal duration or useful-size stage gain.

A [robustness corollary](robust-initial-data.md) also gives strict favorable initial signs with the direct rational swirl amplitude A=1100, without exact neutral tuning. This is a separate local datum; it does not provide an invariant class.

This supersedes the immediate initial-sign task for a **new radial datum**. It does not certify the older cylindrical candidate's approximate value of −241.35. The old data and reports remain unchanged.

The next mathematical obligation is to carry the actual terminal velocity and inherited exterior geometry into a reusable class with retained gain. The [return calculation](actual-return-next-lemma.md) proves an unavoidable first-order cubic core deformation: p_zzzz(0) >= 3693184/105861 > 34. Receiver rescaling cancels the central matrix's first-order change but cannot erase this higher derivative. The [inherited-geometry calculation](inherited-geometry-first-variation.md) supplies exact packet-moment rates and a necessary aspect-ratio return budget, and retains the packet-induced meridional deformation. The [independent finite-time audit](finite-feedback-independent-audit.md) reviews the local corollary and higher-jet calculation.

Favorable initial derivatives do not establish persistence to a singular time, and choosing smaller receivers does not supply an infinite sequence of compatible stages. A proposed larger return class is stated explicitly; its invariance is unproved.

To reproduce, run `python3 certificate/run_checks.py`, then `python3 certificate/certify_full_initial_gate.py --recompute` from this directory. The latter recomputes all four interval integrals. The saved [certificate](certificate/FULL_INITIAL_GATE_CERTIFICATE.json) contains the exact dyadic endpoints. [Verification and limitations](VERIFICATION.md) describe the trust base: analytical arguments and mpmath interval arithmetic, without a formal proof-kernel replay or an NS blowup proof.
''')

(dst/'VERIFICATION.md').write_text('''# Verification and limitations

The certificate combines exact analytic reductions with inclusion-preserving range quadrature using mpmath.iv at 35 decimal digits. Rational panel partitions, endpoint cutoff estimates, cumulative-integral interval propagation, and exact dyadic endpoint comparison are documented in the mathematical notes. Approximate priorities choose which panels to subdivide; they do not determine the enclosure error.

`python3 certificate/run_checks.py` runs the finite symbolic identities, support/kernel bounds, cutoff consistency tests, central-jet and receiver arithmetic, and exact certificate aggregation. CHECK_RESULTS.txt records the packaged run. These checks supplement the analytic proofs and audits; they do not mechanically prove those entire arguments.

`python3 certificate/certify_full_initial_gate.py --recompute` independently regenerates the four integral enclosures used by the assembler: the combined radial meridional cubic, core strain/swirl coefficient, positive outer swirl stress and its viscous coefficient. FULL_INITIAL_GATE_CERTIFICATE.json records the recomputation and all exact endpoints. Omitting --recompute only aggregates stored component bounds.

The analytical dependencies include the earlier radial pressure/cubic identity, extended to the present sign-changing radial profile; the degree-four cutoff for this scalar functional; the core and annular adjoint kernels; support and viscosity cancellations; and local smooth NS theory. Separate audits check the full decomposition, cumulative quadrature, and finite-time deductions. mpmath's interval elementary functions and Python's exact integer/rational operations remain trusted software. No Lean or other formal proof-kernel replay has been performed.

The stated inequality concerns one actual initial-time pressure derivative. Local existence and strict margins also give a positive interval of the specified improving observables. There was no NS time-stepping simulation. No numerical lower bound for a useful stage duration, invariant profile class, inherited-state return, indefinite amplification, assembled all-order smooth forcing, or singular solution is certified.

The exploratory cylindrical and coarse radial pressure evaluations are not inputs to this certificate. Earlier pass2, pass3 and pass4 packages are unchanged. File hashes identify the delivered artifacts and do not enlarge the mathematical claim.
''')

sources={
 'date':'2026-09-08',
 'inherited_manifests':[{'path':str(p.relative_to(base.parent)).replace('euler-ns-transfer/','../',1),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in [base/'SOURCES.json',base/'pass3'/'SOURCES.json',base/'pass4'/'SOURCES.json']],
 'new_external_sources':[],
 'derivations':'Explicit radial datum, full initial NS pressure identity, adjoint kernels, interval range integration and local smooth-solution corollary.',
 'verification':{'python':sys.version,'formal_kernel_replay':False,'NS_time_evolution':False,'full_initial_pressure_inequality_interval_certified':True,'local_positive_interval_proved':True,'numerical_useful_duration_certified':False,'terminal_profile_return_proved':False,'NS_blowup_proved':False}
}
(dst/'SOURCES.json').write_text(json.dumps(sources,indent=2)+'\n')
print(dst)
