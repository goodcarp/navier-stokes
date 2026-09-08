# Bounded primary-source check

Date: 2026-09-08. This lookup concerns publicly readable mathematical statements, not unpublished claims or corporate attribution. It is not an exhaustive literature ranking or a proof audit.

## Newly announced ordinary/hypodissipative Navier–Stokes

Targeted searches of arXiv, Buckmaster's author site, the public `tristanbuckmaster/fluid_lean` repository and OpenAI's official domain did **not locate a newly public ordinary forced NS proof or an Alpöge–Buckmaster hypodissipative paper**. The inspected release remains the Euler/IPM/Boussinesq material already recorded by the parent task. This is a bounded negative lookup, not proof that no new file exists.

Limitations: the author homepage returned limited indexed content; GitHub's web extraction omitted its file listing. A direct GitHub API check failed through both the shell network path and web fetch. The browser tool reported that no browser was available. Therefore this lookup did not establish the latest repository HEAD or certify absence of newly added files. The parent's previously pinned commit remains the available inspected provenance. Relevant landing pages: [Buckmaster](https://cims.nyu.edu/~tristanb/), [public repository](https://github.com/tristanbuckmaster/fluid_lean).

An existing, distinct result is [Córdoba–Martínez-Zoroa–Zheng's earlier hypodissipative blowup theorem](https://arxiv.org/abs/2407.06776). Their dissipation notation is `|nabla|^alpha`, with `0<=alpha<(22-8sqrt(7))/9`; ordinary NS is alpha=2 in that convention. The force is in `L1_t C1,epsilon_x intersect Linfinity_t L2_x`, weaker than smoothness in all space-time derivatives. This is not the newly announced smooth-forcing extension.

## Most useful quantitative actual-NS stage reference found

Jeong–Yoneda, *Vortex stretching and enhanced dissipation for the incompressible 3D Navier–Stokes equations*, Math. Ann. 380 (2021), 2041–2072. The arXiv title retains “anomalous” in its metadata; the inspected v2 PDF states “Enhanced dissipation.” [Primary PDF](https://arxiv.org/pdf/2001.02333), Theorem 1.1, equations (1.5)–(1.6), and §2.3, especially (2.33), (2.37).

The theorem constructs smooth finite-energy periodic data and viscosities `nu_n -> 0` in a 2.5-dimensional stretching setup. For every `0<a0<1` it obtains a viscosity-weighted lower bound on averaged gradient energy relative to initial energy. The horizontal torus can stay fixed when `a0<=1/2`. These are actual NS solutions. Equation (2.37) quantitatively selects viscosity so the stretching survives the comparison error. The small-scale component does not feed back into the planar base.

**Research inference:** §2.3 provides a quantitative benchmark for our proposed frequency window. It does not supply fixed-viscosity regeneration. A three-dimensional coupled packet still needs its own feedback estimate.

## Most directly relevant whole-space finite-energy amplification reference found

Choi–Jeong, *On vortex stretching for anti-parallel axisymmetric flows*, American Journal of Mathematics 147 (2025), 1251–1284. [Publisher primary PDF](https://preprint.press.jhu.edu/ajm/sites/default/files/AJM-choi-jeong-FINAL.pdf), §6, Proposition 6.1, printed pp. 23–24; [arXiv record](https://arxiv.org/abs/2110.09079).

Proposition 6.1 assumes compactly supported bounded axisymmetric no-swirl Euler vorticity with finite kinetic energy, global locally Lipschitz velocity, and unbounded long-time vorticity L2 norm. For finite-energy viscous initial vorticities converging strongly in L2 to that datum, it concludes

\[
\liminf_{\nu\downarrow0}\sup_{t\ge0}\|\omega^\nu(t)\|_{L^2(\mathbb R^3)}=\infty.
\]

The paper supplies Euler examples meeting the premise and recalls global regularity of these no-swirl NS solutions. The argument is an inviscid-limit contradiction and supplies no explicit viscosity-to-gain or time-to-gain bound. It establishes large finite amplification in a family, not a finite-time NS singularity or iterative one-parent regeneration.

## A tempting source that must remain a model reference

Moffatt–Kimura, *Towards a finite-time singularity of the Navier–Stokes equations. Part 3. Maximal vorticity amplification* (JFM, 2023), gives an exact solution to its reduced vortex-ring dynamical system and estimates the vortex Reynolds number needed for a prescribed gain. It explicitly reports eventual violation of the model assumptions and makes no NS-singularity claim. It is useful for circulation/reconnection/saturation diagnostics, but its exact ODE solution is not a rigorously validated actual-NS stage. [Primary paper](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/towards-a-finitetime-singularity-of-the-navierstokes-equations-part-3-maximal-vorticity-amplification/EB8BDB2188734A7FF788E44C1BD7ACD3).

## Decision

Use Jeong–Yoneda §2.3 as the first quantitative comparison reference. Use Choi–Jeong Proposition 6.1 as an actual whole-space amplification baseline and a warning about missing quantitative parameters. No newly public theorem found in this lookup closes the current viscous stage or forcing-summability obligations.
