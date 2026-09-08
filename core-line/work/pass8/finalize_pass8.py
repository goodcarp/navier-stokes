from pathlib import Path
import ast, hashlib, json, re, shutil

root = Path(__file__).resolve().parents[2]
top = root / 'outputs/euler-ns-transfer'
p8 = top / 'pass8'

def verify_manifest(folder):
    lines = (folder / 'SHA256SUMS.txt').read_text().splitlines()
    for line in lines:
        wanted, rel = line.split('  ', 1)
        assert hashlib.sha256((folder / rel).read_bytes()).hexdigest() == wanted, (folder, rel)
    return len(lines)

historical = {f'pass{i}': verify_manifest(top / f'pass{i}') for i in range(2, 8)}
assert list(historical.values()) == [17, 16, 33, 39, 23, 21]
log = (p8 / 'CHECK_RESULTS.txt').read_text()
runs = re.findall(r'^RUN (\S+)$', log, re.M)
passes = re.findall(r'^PASS (\S+)$', log, re.M)
assert len(runs) == 11 and runs == passes
assert 'PASS: 11 finite programs; all three new interval certificates recomputed.' in log
certs = {f.name: json.loads(f.read_text()) for f in (p8 / 'checks').glob('*certificate.json')}
assert len(certs) == 3 and all(x['status'] == 'PASS' for x in certs.values())
leading = certs['leading-envelope-certificate.json']
assert leading['mean_derivative_lower'] == '8871/800'
assert leading['fluctuation_energy_fractional_derivative_lower'] == '12853/2205'
pressure = certs['leading-pressure-certificate.json']
assert pressure['common_T_negative_margin'] == '12493177/1120000'
assert pressure['beta_second_lower'] == '12493177/2240000'

(top / 'STATUS.md').write_text('''# Current research state

Updated 8 September 2026, after the eighth transfer pass.
**The full Navier–Stokes goal remains open. No singular solution, complete
blowup proof or formal proof-kernel verification has been obtained.**

The preceding seventh pass made concrete progress by closing the common
initial pressure and mean-torque test. This eighth pass also makes concrete
progress: it diagnoses energy depletion and initial receiver deceleration
for that family, then constructs a separate family passing three compatible
initial tests, including net perturbation-energy growth after viscosity.

## Latest result and selected candidate

The [eighth-pass report](pass8/REPORT.md) selects the leading-envelope seed
with radial center 7/50, width 1/10, angular mode four, radial carrier -20,
amplitude 1/8 <= lambda <= 1/4, and exact pressure retuning. Viscosity is
1/1000. The [same-datum pressure theorem](pass8/leading-full-pressure-feedback.md)
and [envelope calculation](pass8/leading-envelope-alternative.md) prove:

- Initial receiving mean angular-momentum derivative G_t >8871/800.
- Initial fractional perturbation-energy derivative K'/K >12853/2205,
  including the mean strain and viscous losses.
- Complete initial pressure-feedback test
  p_zz'(0)+32 <-12493177/1120000, hence
  (b/Omega)''(0)>12493177/2240000.

All three inequalities apply to the same new family. The reversed chirality
draws energy from the mean; the declining envelope supplies favorable
torque at the receiving maximum. Total kinetic energy still dissipates.
Strict uniform margins and smooth local dependence give a common positive
interval of simultaneous gain. No useful numerical duration or prescribed
finite gain has been enclosed.

Eleven packaged finite programs passed, including recomputation of three
new interval certificates. Independent analytical audits accepted the
scoped calculations. These are not formal proof-kernel replays or numerical
Navier–Stokes time evolution; see [verification](pass8/VERIFICATION.md).

## What the previous family teaches us

The pass7 trailing family still passes its initial pressure and mean-growth
tests. The [actual energy identity](pass8/actual-fluctuation-energy.md) now
shows K'(0)/K(0)<-43.66 for those data. The [full-pressure initial jets](pass8/initial-receiver-deceleration.md)
also show strong receiver deceleration: the value along its radial maximum
at fixed packet-center axial height has initial second derivative <-240000,
and that radius initially moves outward. This is not a bound for the
unrestricted global maximum or a proof of a turning time.

These old-family conclusions are not transferred to the new leading seed.
The latter is a new initial datum, not an evolved endpoint or a regeneration
of the old solution. The separate Kelvin calculation is a model diagnostic;
the new steep envelope must be carried in the actual dynamics.

## Preserved obstruction and current next action

The [pass6 circulation obstruction](pass6/circulation-return-obstruction.md)
still rules out the former indefinitely repeating bounded axisymmetric
rotating-core class. Earlier local certificates remain valid. No general
axisymmetric regularity theorem is asserted.

The next task is a useful quantitative interval for the actual fully coupled
leading-envelope solution, with retained gain and a controlled inherited
endpoint. Carry the spatial covariance gradient, pressure-generated vertical
velocity, changing mean geometry, viscosity, core curvature and older field.
Then establish a compatible continuation within the same solution.
[NEXT_TARGET.md](NEXT_TARGET.md) preserves the obligations; initial signs,
a fresh replacement packet and qualitative continuity do not supply them.

There is no missing authorization or external tool blocker. The remaining
gap is mathematical, and the full goal remains active.
''')

f = top / 'NEXT_TARGET.md'
s = f.read_text()
a = re.search(r'## Immediate target after the (?:seventh|eighth) pass', s).start()
b = s.index('The [pass6 circulation theorem]', a)
s = s[:a] + '''## Immediate target after the eighth pass

**The selected leading-envelope family now passes three initial tests.**
The [eighth-pass report](pass8/REPORT.md) combines initial receiving
mean-angular-momentum growth, net nonaxisymmetric kinetic-energy growth,
and favorable complete central pressure feedback on the same exactly
retuned data, with viscosity 1/1000. Its radial envelope is centered at
7/50 with width 1/10, angular mode four and radial carrier -20; the seed
amplitude is 1/8 <= lambda <= 1/4.

The certified bounds are G_t>8871/800, K'/K>12853/2205, and
p_zz'(0)+32<-12493177/1120000. The full pressure calculation includes
viscosity and the nonlocal strong-swirl interaction. These give a common
qualitative positive interval, without an enclosed useful duration.
Do not repeat these already closed initial tests as the main target.

The earlier trailing seed still passes the seventh-pass initial test, but
its perturbation energy initially declines at fractional rate below -43.66.
Its packet-center maximum strongly decelerates even with radial recentering.
Those initial jets neither prove its eventual failure nor describe the new
leading seed. The new seed is a separate datum, not same-solution regeneration.

**The next action is quantitative persistence and inherited gain for the
leading-envelope solution.** Select one amplitude in [1/8,1/4] and bound
the actual fully coupled evolution over a useful interval. The steep
envelope is essential to the torque; it cannot be discarded as a small
correction to the old constant-amplitude Kelvin model. Carry the spatial
covariance gradient, pressure-generated vertical velocity, radial and axial
deformation, viscosity, core curvature and older field. The exact covariance
and whole-space pressure identities in pass8 are available for this step.
The former order-10^-3 diagnostic window is not a certified duration or
an automatic scale choice for these new data.

Prove a finite retained full-field gain and a named norm bound on the actual
endpoint that enables a successor stage. No affine-core reset, fresh outer
packet, prescribed mean stress or differentiated radial-pressure closure is
permitted. Qualitative continuity alone does not give a prescribed useful gain.
Uniform repeated compatibility remains a separate requirement after one stage.

''' + s[b:]
f.write_text(s)

f = top / 'REPORT.md'
s = f.read_text(); a = s.index('**Latest:**'); b = s.index('\n\n', a)
s = s[:a] + '''**Latest:** [The eighth-pass report](pass8/REPORT.md) diagnoses initial energy depletion and receiver deceleration in the previous candidate, then constructs a separate leading-envelope family with three compatible initial gains: receiving mean growth, net perturbation-energy growth and favorable full central pressure feedback, including viscosity. The [sixth-pass axisymmetric return obstruction](pass6/circulation-return-obstruction.md) and earlier local certificates remain valid. A useful quantitative stage, inherited return and blowup remain open; [STATUS.md](STATUS.md) records the current state.''' + s[b:]
f.write_text(s)
f = top / 'VERIFICATION.md'
s = f.read_text(); a = s.index('This file records the first pass.'); b = s.index('\n\n', a)
s = s[:a] + '''This file records the first pass. The latest [eighth-pass verification](pass8/VERIFICATION.md) records eleven successful finite programs, including three recomputed interval certificates for the previous family's initial receiver jets and the new leading family's compatible initial gains. The pass2–pass7 packages remain unchanged historical records. No actual NS time integration, useful numerical stage duration, inherited return or blowup has been proved.''' + s[b:]
f.write_text(s)

f = p8 / 'REPORT.md'; s = f.read_text()
s = s.replace('Run python3 checks/run_checks.py to recompute the three new interval\ncertificates and the eleven finite programs.', 'Run `python3 checks/run_checks.py` to execute all eleven programs, including\nrecomputation of three new interval certificates.')
f.write_text(s)
f = p8 / 'VERIFICATION.md'; s = f.read_text()
s = s.replace('the separate\nKelvin and endcap models', 'the separate\nKelvin model and exact endcap estimate')
s = s.replace('Run python3 checks/run_checks.py. It regenerates:', 'All eleven programs passed; the saved run is in [CHECK_RESULTS.txt](CHECK_RESULTS.txt).\nRun `python3 checks/run_checks.py`. It regenerates:')
f.write_text(s)

f = p8 / 'SOURCES.json'; d = json.loads(f.read_text())
d['verification'].update(finite_check_programs_passed=11, new_interval_certificates_recomputed=True, earlier_packages_unchanged=True, earlier_package_manifest_entries=historical)
f.write_text(json.dumps(d, indent=2) + '\n')

for cache in p8.rglob('__pycache__'):
    shutil.rmtree(cache)
for f in p8.rglob('*.py'):
    ast.parse(f.read_text(), filename=str(f))
for f in list(p8.rglob('*.md')) + [top / x for x in ['STATUS.md', 'NEXT_TARGET.md', 'REPORT.md', 'VERIFICATION.md']]:
    s = f.read_text()
    assert not any(ord(c) < 32 and c not in '\n\r\t' for c in s), f
    prose = re.sub(r'`[^`]*`', '', s)
    for link in re.findall(r'\]\(([^)]+)\)', prose):
        if '://' in link or link.startswith('#'):
            continue
        target = link.split('#')[0]
        assert (f.parent / target).exists(), (f, target)

def write_manifest(folder, files):
    text = ''.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(folder)}\n' for f in sorted(files))
    (folder / 'SHA256SUMS.txt').write_text(text)

write_manifest(p8, [f for f in p8.rglob('*') if f.is_file() and f.name != 'SHA256SUMS.txt' and '__pycache__' not in f.parts])
write_manifest(top, [f for f in top.iterdir() if f.is_file() and f.name != 'SHA256SUMS.txt'] + list((top / 'checks').glob('*.py')))
assert verify_manifest(top) == 10
new_count = verify_manifest(p8)
assert {f'pass{i}': verify_manifest(top / f'pass{i}') for i in range(2, 8)} == historical
print(json.dumps(dict(status='PASS', programs=11, new_interval_certificates=3, pass8_manifest_entries=new_count, top_manifest_entries=10, unchanged_historical_entries=historical, local_links='resolved', python_syntax='valid'), indent=2))
