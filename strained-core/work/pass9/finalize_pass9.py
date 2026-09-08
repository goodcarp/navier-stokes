#!/usr/bin/env python3
from pathlib import Path
import ast,hashlib,json,re,shutil
root=Path(__file__).resolve().parents[2];top=root/'outputs/euler-ns-transfer';p9=top/'pass9'

def verify_manifest(folder):
    lines=(folder/'SHA256SUMS.txt').read_text().splitlines()
    for line in lines:
        wanted,relative=line.split('  ',1)
        assert hashlib.sha256((folder/relative).read_bytes()).hexdigest()==wanted,(folder,relative)
    return len(lines)

history={f'pass{i}':verify_manifest(top/f'pass{i}') for i in range(2,9)}
assert list(history.values())==[17,16,33,39,23,21,35]
log=(p9/'CHECK_RESULTS.txt').read_text()
runs=re.findall(r'^RUN (\S+)$',log,re.M);passes=re.findall(r'^PASS (\S+)$',log,re.M)
assert len(runs)==8 and runs==passes
assert 'PASS: 8 finite programs; two interval certificates recomputed' in log
for name in ['general-envelope-initial-jet-certificate.json','affine-envelope-budget-certificate.json','fixed-amplitude-initialization-certificate.json']:
    assert json.loads((p9/'checks'/name).read_text())['status']=='PASS'

(top/'STATUS.md').write_text('''# Current research state

Updated 8 September 2026, after the ninth transfer pass.
**The full Navier–Stokes goal remains open. No singular solution, complete
blowup proof or formal proof-kernel verification has been obtained.**

The preceding eighth pass made concrete progress by constructing and auditing
three compatible initial gains. This ninth pass also makes concrete progress:
it encloses the new seed's actual initial receiver deceleration, evaluates
a resolved nonlinear planar gain pulse, proves a finite gain in an exact
affine linear envelope model, and derives a full-field validation route.
An explicit rational amplitude now removes implicit pressure retuning from
the selected next initial datum.

## Selected next full-field test

Use `u0=M+1020v+(1/4)w_L`, with the unchanged compact leading-envelope
profiles and viscosity `1/1000`. This is an explicit initial datum, not a
regenerated endpoint. The [fixed-amplitude calculation](pass9/fixed-amplitude-initialization.md)
and its independent audit give actual initial bounds:

- Core ratio derivative `beta'(0)>0.02484` and acceleration `beta''(0)>6.09`.
- Receiving mean derivative `G_t(0,r_*,4)>102.49`.
- Fractional nonaxisymmetric energy derivative `K'(0)/K(0)>5.829`.

The amplitude intentionally gives a strict initial core cone. Its correct
off-neutral identity is used; the old neutral shortcut is not propagated.
The pressure-neutral pass8 family remains a valid historical result.

## What the new evolution evidence establishes

The [full-pressure new-envelope jet](pass9/general-envelope-jet-evaluation.md)
retains the steep envelope, viscosity and axial covariance flux. Its fresh
whole-space pressure-gradient error is below 1.3. Uniformly over the stated
independent amplitude ranges, the initial fixed-circle acceleration is
below -13000; tracking its radial maximum at fixed packet-center height
still gives acceleration below -12000. These are initial statements, not
turning-time bounds or acceleration bounds for the global spatial maximum.

In the nonlinear planar model, the selected amplitude shows about 0.175%
gain in the tracked mean branch and 14.6% fluctuation-energy gain at time
0.001. Radial/time refinement supports that observation. A longer A=950
comparison run shows a later torque reversal and loss of fixed-receiver
gain. These are finite-radius slice calculations; they omit the compact
axial pressure, evolving exterior strain, central core and older field.
They are not time evolution of the full compact 3D datum.

The separate affine linear model has a certified finite-time energy gain
above 15% for the entire envelope and a finite positive shear-work budget.
Its energy normalization, local torque and missing nonlinear terms are
kept separate from both the slice and the actual NS solution.

## Current next action

Build and validate one complete compact three-dimensional approximation
from the explicit datum, carrying all mean, axial, core and older-field
effects. The [H4 endpoint bridge](pass9/evolution-validation-bridge.md)
can certify retained mean, energy, core-cone and full-velocity gains without
requiring every intermediate derivative to stay positive. It must also
enclose the actual inherited endpoint in a successor class stable in its
named norm. H7 curvature validation remains an optional stronger route.

No full-field residual, useful actual duration or terminal-class margin
has yet been enclosed. Numerical reconstruction, exterior tails, domain
errors and every time-slab error still need to be carried. A fresh packet
or a reduced-model pulse cannot replace this step.
[NEXT_TARGET.md](NEXT_TARGET.md) retains the full-stage and infinite-assembly
obligations.

Eight packaged finite programs passed, including two recomputed interval
certificates and the fixed-amplitude exact algebra. Pass2–pass8 manifests
remain unchanged. See the [ninth-pass report](pass9/REPORT.md) and
[verification](pass9/VERIFICATION.md) for scope.

The pass6 bounded axisymmetric rotating-core return obstruction remains
valid under its stated assumptions; no general axisymmetric regularity
theorem is asserted. There is no missing authorization or external blocker.
The remaining gap is mathematical, and the full goal remains active.
''')

f=top/'NEXT_TARGET.md';s=f.read_text()
a=re.search(r'## Immediate target after the (?:eighth|ninth) pass',s).start()
b=s.index('The [pass6 circulation theorem]',a)
s=s[:a]+'''## Immediate target after the ninth pass

**The next test has an explicit initial field.** Choose
`u0=M+1020v+(1/4)w_L`, viscosity `1/1000`, with the unchanged compact
leading envelope centered at 7/50, width 1/10, angular mode four and radial
carrier -20. The [fixed-amplitude theorem](pass9/fixed-amplitude-initialization.md)
preserves initial mean and fluctuation-energy gains and gives a strict
initial core cone. No implicitly retuned amplitude has to be represented.
This field is not a regenerated endpoint of an earlier solution.

**The next required result is a validated full 3D endpoint with retained
gain and inherited geometry.** The [H4 bridge](pass9/evolution-validation-bridge.md)
supplies sufficient residual, error-majorant and endpoint tests. Construct
a smooth, divergence-free whole-space approximation including the actual
axial envelopes, generated vertical motion, evolved mean, central core and
older field. Bound its complete vector residual and initial reconstruction
error. Every subsequent time slab inherits the previous error radius.

Validate prescribed endpoint mean, fluctuation-energy, core-cone and
full-field velocity/Reynolds gains as needed by the successor mechanism.
Prove actual endpoint membership in a named stable successor class after
charging rescaling and coordinate errors. H4 is sufficient only if that
successor lemma is stable in H4. Positive derivative or curvature signs at
every intermediate time are optional stronger conditions, not automatic
requirements of retained gain. Use the H7 alternative if those conditions
or a stronger successor class are needed.

The finite planar simulations select an early time around 0.001 for testing,
not a certified compact-3D duration. At the selected A=1020 they show about
0.175% tracked mean-branch gain and 14.6% finite-radius fluctuation-energy
gain at that time, with refinement checks. The longer A=950 comparison
shows a later torque reversal. Its timing cannot be assigned to the actual
selected solution. The exact affine-envelope finite gain is also a separate
linear model, with a finite shear-work budget.

The actual full-pressure initial-jet certificate now shows deceleration
for the leading seed too, including radial recentering at fixed axial height.
It does not supply an actual turning time or rule out a useful endpoint gain.
Do not substitute those initial jets for a time-dependent error enclosure.

Raw slice data are not already a full-field H4 approximation. They omit
three-dimensional effects and exterior energy tails; residual numerical
mean circulation gives a nonintegrable whole-plane mean-energy tail unless
it is removed in a controlled reconstruction. A small finite-grid energy
balance diagnostic does not bound the complete 3D residual.

The pass8 three-sign result, pass9 explicit-amplitude inequalities and
reduced-model finite pulse are closed diagnostic work. Do not repeat them
as the main target. No affine-core reset, fresh outer seed, prescribed mean
stress or pressure closure is a substitute for the evolved endpoint.
Uniform repeated compatibility and admissible all-order force assembly
remain separate after one actual stage.

'''+s[b:];f.write_text(s)

f=top/'REPORT.md';s=f.read_text();a=s.index('**Latest:**');b=s.index('\n\n',a)
s=s[:a]+'''**Latest:** [The ninth-pass report](pass9/REPORT.md) adds the leading seed's full-pressure initial deceleration, a resolved finite pulse in a nonlinear planar model, an exact affine-envelope finite-gain theorem and an H4 full-field endpoint-validation route. The selected next initial field now has explicit rational amplitude and an audited strict core cone. No compact 3D persistence interval, inherited return or blowup has been proved; [STATUS.md](STATUS.md) records the current state.'''+s[b:];f.write_text(s)
f=top/'VERIFICATION.md';s=f.read_text();a=s.index('This file records the first pass.');b=s.index('\n\n',a)
s=s[:a]+'''This file records the first pass. The latest [ninth-pass verification](pass9/VERIFICATION.md) records eight successful finite programs, including two recomputed interval certificates, explicit-amplitude exact algebra and scoped nonlinear slice checks. Pass2–pass8 remain unchanged historical packages. Actual initial identities, exact reduced-model statements, numerical slice diagnostics and conditional full-field validation are distinguished. No full compact 3D time integration, useful actual duration, inherited return or blowup is established.'''+s[b:];f.write_text(s)

f=p9/'VERIFICATION.md';s=f.read_text()
if 'All eight programs passed.' not in s:
    s=s.replace('The packaged runner executes eight finite programs:', 'All eight programs passed. The packaged output is saved in\n[CHECK_RESULTS.txt](CHECK_RESULTS.txt).\n\nThe runner executes eight finite programs:')
f.write_text(s)
f=p9/'SOURCES.json';d=json.loads(f.read_text());d['verification'].update(finite_programs_passed=8,new_interval_certificates_recomputed=True,earlier_packages_unchanged=True,earlier_manifest_entries=history)
f.write_text(json.dumps(d,indent=2)+'\n')

for cache in p9.rglob('__pycache__'):shutil.rmtree(cache)
for f in p9.rglob('*.py'):ast.parse(f.read_text(),filename=str(f))
for f in list(p9.rglob('*.md'))+[top/x for x in ['STATUS.md','NEXT_TARGET.md','REPORT.md','VERIFICATION.md']]:
    s=f.read_text();assert not any(ord(c)<32 and c not in '\n\r\t' for c in s),f
    prose=re.sub(r'`[^`]*`','',s)
    for link in re.findall(r'\]\(([^)]+)\)',prose):
        if '://' in link or link.startswith('#'):continue
        assert (f.parent/link.split('#')[0]).exists(),(f,link)
def write_manifest(folder,files):
    (folder/'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(folder)}\n' for f in sorted(files)))
write_manifest(p9,[f for f in p9.rglob('*') if f.is_file() and f.name!='SHA256SUMS.txt' and '__pycache__' not in f.parts])
write_manifest(top,[f for f in top.iterdir() if f.is_file() and f.name!='SHA256SUMS.txt']+list((top/'checks').glob('*.py')))
assert verify_manifest(top)==10
n=verify_manifest(p9)
assert {f'pass{i}':verify_manifest(top/f'pass{i}') for i in range(2,9)}==history
print(json.dumps(dict(status='PASS',finite_programs=8,interval_certificates_recomputed=2,pass9_manifest_entries=n,top_manifest_entries=10,historical_manifest_entries=history,links='resolved',Python_syntax='valid',figure='visually reviewed'),indent=2))
