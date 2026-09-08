#!/usr/bin/env python3
from pathlib import Path
root=Path(__file__).resolve().parents[2]
top=root/'outputs/euler-ns-transfer'
(top/'STATUS.md').write_text("""# Current research state

Updated 8 September 2026, after the seventh transfer pass.
**The full Navier–Stokes goal remains open. No singular solution, complete
blowup proof or formal proof-kernel verification has been obtained.**

The immediately preceding status-only turn added no new progress. This turn
makes concrete progress: the two previously unevaluated full-pressure
interaction coefficients are now bounded, and the finite outer family
passes both initial requirements on the same actual datum.

## Latest result

The [seventh-pass common-family theorem](pass7/common-family-feedback.md)
uses the unchanged outer chiral seed, exact pressure retuning, viscosity
1/1000, and every amplitude 3/4 <= lambda <= 1. It proves:

- D_nu < 4-(96/5)Cw, including all core/pump pressure and viscosity.
- E < 1/60, retaining the strong-swirl cross pressure and its full tails.
- p_zz'(0)+32 < -435297/70000 and (b/Omega)''(0)>435297/140000.
- The initial global positive azimuthal-mean angular-momentum maximum
  grows at a fixed maximizing circle at rate >626/35, by the earlier
  unchanged calculation for this very same family.

Smooth local existence and uniform strict margins give a common short
interval with simultaneous actual core and mean-maximum gains.
No useful numerical duration or prescribed finite gain has been enclosed.
The inner infinitesimal perturbation from pass6 is not substituted here.

The [full D proof](pass7/D-pressure-kernel.md), [full E reduction](pass7/E-interaction-reduction.md),
interval certificates and independent audits are packaged in the
[seventh-pass report](pass7/REPORT.md). Both decisive interval calculations
are recomputed by its verification runner.

## Preserved obstruction and current next action

The [pass6 circulation obstruction](pass6/circulation-return-obstruction.md)
still rules out the former indefinitely repeating bounded axisymmetric
rotating-core class. The earlier pass5 local certificate remains valid.
No general axisymmetric regularity theorem is asserted.

The next task is a useful quantitative interval for the actual fully coupled
three-dimensional solution that retains gain and supplies a controlled
inherited endpoint. The evolving angular correlation, pressure, viscosity,
core curvature, older field and outer shape must all be carried.
An initial sign or passive phase model does not supply that estimate.
After one useful stage, a true return or other infinitely compatible
continuation remains separate. [NEXT_TARGET.md](NEXT_TARGET.md) preserves
the full obligations and forbids packet or affine-profile resets.

There is no missing authorization or external tool blocker. The remaining
gap is mathematical, and the full goal remains active.
""")
next_file=top/'NEXT_TARGET.md'
s=next_file.read_text()
start=s.index('## Immediate target after the sixth pass')
end=s.index('## Progress and remaining gap after the third pass')
new="""## Immediate target after the seventh pass

**The common finite-amplitude initial test is now passed.**
The [seventh-pass theorem](pass7/common-family-feedback.md) joins the full
central pressure-feedback inequality and initial growth of the positive
global mean-angular-momentum maximum on the same outer family, for every
3/4 <= lambda <= 1. The exact retuning is unchanged.

The formerly unknown coefficients satisfy D_nu<4-(96/5)Cw and E<1/60.
The latter retains the complete nonlocal strong-swirl cross pressure by an
exact horizontal Green trial and an explicit whole-space residual bound.
Consequently p_zz'(0)+32<-435297/70000. The same data already have initial
mean-maximum derivative >626/35, including viscosity. The strict uniform
margins give simultaneous short-time gains, but no useful numerical duration.
Do not repeat the already closed two-coefficient test as the main target.

**The next action is quantitative persistence and inherited gain.**
At one selected amplitude (lambda=1 is admissible), bound the actual fully
coupled evolution over a useful interval. Track the evolving nonaxisymmetric
correlation, radial and axial wave covector, pressure polarization, viscosity,
core curvature, outer deformation, and the old field. A frozen-background
phase estimate may diagnose scales but cannot replace these dynamics.

Prove a finite retained full-field gain and a named norm bound on the actual
endpoint that enables a successor stage. No affine-core reset, fresh outer
packet, prescribed mean stress or differentiated radial-pressure closure is
permitted. Qualitative continuity alone does not give a prescribed useful gain.
Uniform repeated compatibility remains a separate requirement after one stage.

The [pass6 circulation theorem](pass6/circulation-return-obstruction.md)
continues to rule out the former bounded axisymmetric rotating-core return.
The conditional escape audit also remains valid under its stated stronger
global assumptions. These do not exclude all second-scale mechanisms or
prove general axisymmetric regularity.

The fixed-receiver stress budget remains |Q_L(v)|<=C_chi ||v||2². A required
increment D supplied by this stress costs integrated perturbation energy
at least 2 I_L D/C_chi. Frequency alone cannot evade that cost.

"""
next_file.write_text(s[:start]+new+s[end:])
p=top/'REPORT.md';s=p.read_text();a=s.index('**Latest:**');b=s.index('\n\n',a)
s=s[:a]+"""**Latest:** [The seventh-pass report](pass7/REPORT.md) closes the finite-amplitude three-dimensional initial compatibility test: the same exactly retuned outer family has favorable full pressure feedback and growing mean angular momentum, including viscosity. The [sixth-pass axisymmetric return obstruction](pass6/circulation-return-obstruction.md) and [fifth-pass local certificate](pass5/REPORT.md) remain valid. A useful quantitative stage, inherited return and blowup remain open; [STATUS.md](STATUS.md) records the current state."""+s[b:];p.write_text(s)
p=top/'VERIFICATION.md';s=p.read_text();a=s.index('This file records');b=s.index('\n\n',a)
s=s[:a]+"""This file records the first pass. The latest [seventh-pass verification](pass7/VERIFICATION.md) recomputes the two full-pressure interval bounds and checks their common finite-amplitude consequence. The [sixth-pass verification](pass6/VERIFICATION.md) and [fifth-pass verification](pass5/VERIFICATION.md) remain unchanged historical records. No useful numerical stage duration, inherited return or blowup has been proved."""+s[b:];p.write_text(s)
print('Updated STATUS, NEXT_TARGET, top REPORT and top VERIFICATION.')
