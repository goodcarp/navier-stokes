<!-- DRAFT for the author's pass. -->

# What happened on 8 September, and what we did about it

## The announcement

Levent Alpöge and Tristan Buckmaster posted three papers and a statement. The papers prove that the equations of an ideal fluid, of a stratified fluid, and of flow through porous rock can all be driven to a singularity in finite time by a smooth force, starting from smooth data, and that a computer has checked the proofs in Lean. The statement says two more things. They believe they can do it with a little viscosity too. And an outside lab had told them, days earlier, that it had done it for the real thing: the Navier–Stokes equations with full viscosity, with a smooth force, which is one of the two ways to win the Clay prize as Fefferman wrote it.

Nobody outside that lab has seen that proof. So the first useful thing to do was not to guess whether it exists. It was to work out what any such proof would have to look like.

## Two fences

A blowup is a race between a fluid's tendency to pile up gradients and viscosity's tendency to smooth them. Write viscosity as a fractional operator of order α, so that real Navier–Stokes is α = 2. Now ask: what does a bounded velocity buy you?

The answer is scale. Rescale a solution so that lengths shrink by λ; its velocity grows like λ to the power α minus one. At α = 1 that power is zero, so bounded velocity says nothing about the small scales. Above α = 1 it says everything: an energy estimate closes, and the solution continues. We proved this in full for every order between one and two. It is Lemma A, and it is not deep. What it does is decisive: the released construction, and the whole family it belongs to, keeps the velocity bounded by design. Only the gradients blow up. So that family cannot reach any order above one, force or no force, and Navier–Stokes is at two.

The second fence is geometric. In an axisymmetric flow, rotate a singular point about the axis and you get a whole circle of singular points. Caffarelli, Kohn and Nirenberg proved in 1982 that the singular set is too thin to contain a circle. So an axisymmetric singularity sits on the axis. The released construction lives in a doughnut away from the axis. That is fine for Euler and for weak viscosity, where the thinness theorem degrades, and not fine at full viscosity.

Put the two together with two older theorems about the rate, and you get a corridor. Any axisymmetric solution of the Clay class that blows up does so on the axis, faster than the self-similar rate, with swirl. We wrote that up as a theorem, in the Clay class exactly, with the one hypothesis it needs stated in the open.

## What their certificates say

Their Lean repository is real and it is careful. Three projects, each stating its theorem over Mathlib alone and proving it with the standard axioms and nothing else. We put thirty-three independent readers on the statement files, two refuters on every finding, and found no loophole. What we did find is that the certificates are quieter than the papers. The Euler certificate drops the axisymmetry, the swirl, the doughnut and the circulation blowup that the paper's theorem states; the library proves them, the bridge to the statement discards them, and the fidelity note does not say so. The Boussinesq certificate certifies a solution the paper itself says is not the printed one. And no hash anywhere binds a paper to the repository commit. None of this makes anything false. It means the phrase "Lean-verified" is doing less than a reader assumes, and we wrote down exactly how much less.

We also started compiling their Euler certificate on a 2019 laptop, Mathlib from source, because nobody had rebuilt it independently. It is running on a 2019 laptop as this is written, Mathlib from source; the axiom printout goes into the repository the moment it lands.

## What their mechanism can and cannot be pushed to

The released papers add layers at ever higher frequency and ever smaller amplitude, so the envelope never grows and the slope always does. We took the amplitude system straight from their paper, added the damping that viscosity would contribute, and asked how much damping it survives. The answer, reproduced to four digits by a toy model: at most a dissipation order of one half, and only after re-optimizing everything. As written, the schedules survive none at all, because they widen the frequency gap at every stage, and the correction budget that keeps the force smooth is the same budget that forbids a fixed damping. So one prediction is on record: their unreleased viscous paper cannot be the released construction plus damping.

## The clock

Everything above is a fence around someone else's idea. Our own line asks a different question. Never mind whether a blowup exists; how fast can the quantity a blowup needs, the vorticity, actually double?

The lower bound was on the record from an earlier sitting: doubling takes at least a constant over the amplitude times the logarithm of the Reynolds number. On 8 September we assembled the matching upper bound, for an explicit family of vortex rings, with an explicit constant. Two things had to be true for the assembly to close. One was a second-derivative bound that the first proof of the day established with no logarithm at all. The other was a bootstrap that had carried a factor of e to the thirty-four in every previous attempt. That factor turned out to be self-inflicted: the argument compared the true flow to a reference strain fixed in advance, and priced the difference with a crude constant that fed back into itself. Define the reference from the true solution's own strain, and the error term is zero, not small. The factor becomes 3.4.

The assembly is written end to end. Its own note lists what remains, and a blind review from the other line added five bookkeeping items nobody on this side had caught; they are being cleared. When it closes, a published theorem that claims a logarithm-free doubling bound is false as stated, by a family that satisfies every hypothesis its authors wrote.

## The strained core

The other line spent the day building. Ordinary Navier–Stokes, no force, viscosity one over a thousand. A rotating core with a little strain, a pump outside it, and a receiver inside it. Eleven passes, each with an independent audit and check scripts. What it has: certified initial inequalities, one finite gain theorem on one real solution, a numerical lead reproduced to the digit, an explicit next datum. What is still open: a stage that persists, a return to a state it can restart from, and the cascade they would make. Its most important result so far is a negative one, and it is the kind that saves a year. Under the scaling any cascade of this kind must use, the normalized circulation shrinks at every stage, so a rotating core cannot be what drives the cascade forever. The line has already moved in response, to a non-axisymmetric envelope, which is the one direction the fences leave open.

## What Lean does here, and does not

Our own Lean is ninety-one small theorems: the exponent algebra behind Lemma A, the measure of a circle behind Lemma B, the exponent window, and, more than we expected, the two coefficients of the core pressure, which turned out to be theorems for every cutoff profile in the class. No partial differential equation is formalized on our side. Their certificate is a different object: a real formalization of a real theorem, which we are rebuilding rather than trusting. Keeping those two things straight is most of what "Lean-verified" should mean.

## What would change everything

A resolved three-dimensional run of the strained-core datum past its first torque reversal. A written return lemma. Or that manuscript, which will have to pass through both fences on its first page.
