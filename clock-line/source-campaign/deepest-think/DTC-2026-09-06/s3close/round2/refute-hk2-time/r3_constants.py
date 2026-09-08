#!/usr/bin/env python3
"""r3 -- closed-form checks: (1) hk2 sec 7's Gaussian tail; (2) the growth factor e^{int Gamma}
per V-b (inherited) vs per the record's proved Theorem Gamma; (3) diffusion length vs the
strained equatorial-layer half-thickness; (4) the global sup of ||Hess eta_0||_F for (D-B), which
multiplies the tail; (5) the TV-growth exponent of a transported scalar's gradient measure;
(6) the assembly's consumed constant vs hk2's window sup."""
import json, math, numpy as np
from scipy.stats import chi2
from rlib import *
out = {}
delta = math.radians(7.5); sd = math.sin(delta)
# (1) tail.  hk2: nu tau = (rho0 delta)^2 theta_max / L  (V-b campaign normalisation, delta in radians? or sin delta?)
for lbl, dd in (('delta_rad', delta), ('sin_delta', sd)):
    for L in (10.0, 40.0, 160.0):
        nutau = dd**2*THETA_MAX/L
        row = {}
        for d in (0.38, 0.2029, 0.1672, 0.1662):
            x = d*d/(8*nutau)                        # hk2's 'exponent' = (d/2)^2 / (2 nu tau) = chi^2_5 threshold
            # 5-D Brownian with generator nu Lap: per-coordinate variance 2 nu t.  |Z|^2/(2 nu tau) ~ chi2_5.
            p = float(chi2.sf(x, 5))
            row[str(d)] = dict(hk2_exponent_x=x, exp_minus_x=math.exp(-x), P_chi2_5=p, exp_minus_x_over_2=math.exp(-x/2))
        out[f'tail_{lbl}_L{int(L)}'] = row
nutau10 = sd**2*THETA_MAX/10.0
out['hk2_formula_at_d0.38_L10'] = 0.38**2*10/(8*sd**2*THETA_MAX)
out['hk2_formula_at_d0.38_L10_delta_rad'] = 0.38**2*10/(8*delta**2*THETA_MAX)
print("hk2's own formula d^2 L/(8 (rho0 delta)^2 theta_max) at d=0.38, L=10:",
      f"{out['hk2_formula_at_d0.38_L10']:.4f} (sin delta)   {out['hk2_formula_at_d0.38_L10_delta_rad']:.4f} (delta rad)   vs printed 1.444")
for k in ('tail_sin_delta_L10', 'tail_delta_rad_L10'):
    print(k)
    for d, v in out[k].items():
        print(f"   d={d}: x={v['hk2_exponent_x']:.3f}  e^-x={v['exp_minus_x']:.3e}  e^-x/2={v['exp_minus_x_over_2']:.3e}  P(chi2_5>x)={v['P_chi2_5']:.3e}")
# (2) growth factors
g = {}
for Cp in (0.0, 76.11, 119.33, 151.15):
    g[str(Cp)] = {str(L): 2.25*math.exp(Cp*THETA_MAX/L) for L in (10, 40, 160, 640, 2560, 10240)}
out['e_int_Gamma'] = g
print("e^{int_0^tau Gamma} = (9/4) e^{C' theta_max/L}:")
for Cp, row in g.items():
    print(f"  C'={Cp:>7}: " + "  ".join(f"L={L}:{v:.4g}" for L, v in row.items()))
# (3) diffusion length vs strained layer thickness at lam = 3/2
lam = 1.5; r, z = strain_point(30.0, lam); rho = math.hypot(r, z)
layer_half = 0.20*rho/lam**2          # tanh(cos phi / w): |z| <~ w rho in the label; z shrinks by lam^-2
out['layer'] = dict(sqrt_nu_tau_L10=math.sqrt(nutau10), sqrt_nu_tau_L40=math.sqrt(nutau10/4),
                    strained_layer_half_thickness_at_3over2=layer_half, rho_at_3over2=rho,
                    ratio_L10=math.sqrt(nutau10)/layer_half)
print("diffusion length sqrt(nu tau) at L=10:", f"{math.sqrt(nutau10):.4f}", " strained layer half-thickness at lam=3/2:",
      f"{layer_half:.4f}", " ratio", f"{out['layer']['ratio_L10']:.3f}")
# (4) global sup of ||Hess eta_0||_F for (D-B), dense grid over the support
F = Field(7.5, 0.20, 10.0, 1.0)
rr = np.linspace(0.005, 3.0, 1200); zz = np.linspace(-3.0, 3.0, 2400)
RR, ZZ = np.meshgrid(rr, zz, indexing='ij')
H = F.hessF(RR, ZZ); G = F.grad(RR, ZZ)
i = np.unravel_index(np.argmax(H), H.shape); j = np.unravel_index(np.argmax(G), G.shape)
out['global_sup_hessF_DB'] = dict(value=float(H[i]), r=float(RR[i]), z=float(ZZ[i]),
                                  phi_deg=math.degrees(math.atan2(RR[i], ZZ[i])))
out['global_sup_grad_DB'] = dict(value=float(G[j]), r=float(RR[j]), z=float(ZZ[j]),
                                 phi_deg=math.degrees(math.atan2(RR[j], ZZ[j])))
print("global sup ||Hess eta_0||_F (D-B) ~", f"{H[i]:.2f} at r={RR[i]:.3f} z={ZZ[i]:.3f} (phi={out['global_sup_hessF_DB']['phi_deg']:.1f} deg)",
      "; sup|grad eta_0| ~", f"{G[j]:.2f} at phi={out['global_sup_grad_DB']['phi_deg']:.1f} deg")
# tail term in Lemma 7.2 at the d actually used: ||Hess||_inf * P(|Z|>d/2) vs the plateau-ish ball value
for d in (0.1662, 0.1672, 0.2029):
    x = d*d/(8*nutau10); p = float(chi2.sf(x, 5))
    out[f'tail_term_L10_d{d}'] = dict(P=p, hess_tail=p*float(H[i]), grad_tail=p*float(G[j]))
    print(f"  L=10, d={d}: P(|Z|>d/2)={p:.3e}; ||Hess||_inf*P = {p*H[i]:.3f}; |grad|_inf*P = {p*G[j]:.4f}")
# (5) TV growth of the gradient measure: d/dt int|grad eta| <= int (||grad b|| + div b)|grad eta|, div_5 b = 2a
out['TV_growth_exponent_pure_strain'] = dict(formula='exp(int (||grad b|| + div_5 b)) = exp(4 int a) = lam^4', at_3over2=1.5**4)
# (6) assembly mismatch
out['assembly_consumed_CK'] = 66.6622; out['hk2_window_sup'] = 161.7735
out['assembly_mismatch_ratio'] = 161.7735/66.6622
print("assembly a8 used C_K = 66.6622; hk2 sec 0 window sup = 161.7735; ratio", f"{out['assembly_mismatch_ratio']:.4f}")
json.dump(out, open('r3_results.json','w'), indent=1); print("WROTE r3_results.json")
