#!/usr/bin/env python3
"""x8 -- two minor checks.
 (i) PROOF.md writes the datum's finest scale as 'rho0 delta = sqrt(nu/M)' with delta = 7.5
     deg, but b5_sizing.py uses sin(delta).  How much does the 7-digit headline constant move?
 (ii) reconcile 'conservative by 1.4401' with gap-V sec 5(iii)'s 'optimistic by 6.2':
     gap-V's 6.2 = 5(e^{2c}-1)/(c(1+f)^2) compares prove-lagrangian's nu tau/r^2 with gap-V's
     UPPER BOUND ||Hess||_op E|Z|^2, not with the true value; the true (affine) value is the
     trace contraction.  Compute all three."""
import json, math
R={}
TH=4*(1-math.sqrt(2/3)); I2=4/5-16*math.sqrt(6)/135; I4=-4/7+27*math.sqrt(6)/28
c=math.sqrt(1.5)*TH; n=5
d_deg=7.5; d_rad=math.radians(d_deg); d_sin=math.sin(d_rad)
phi0=math.radians(30.0); r0_over_rho0=math.sin(phi0)
R['eps_bulk_L_using_sin_delta']=(d_sin/r0_over_rho0)**2*I2
R['eps_bulk_L_using_delta_radians']=(d_rad/r0_over_rho0)**2*I2
R['ratio']=R['eps_bulk_L_using_delta_radians']/R['eps_bulk_L_using_sin_delta']
print(f"(i) eps_bulk*L with sin(delta): {R['eps_bulk_L_using_sin_delta']:.9f}   "
      f"with delta in radians: {R['eps_bulk_L_using_delta_radians']:.9f}   "
      f"ratio {R['ratio']:.6f}")
print("    PROOF.md displays 0.03473454, i.e. the sin(delta) convention, while its prose says "
      "'rho0 delta'.")
# (ii)
# true affine leading term / (nu tau / r0^2)  =  I2/theta_max
R['true_over_lagrangian']=I2/TH
R['lagrangian_over_true']=TH/I2
# gap-V's upper bound: (1/2)||Hess eta_P||_op E|Z|^2 / |eta_0| with E|Z|^2 <= n nu tau (e^{2c}-1)/c
# = (1/2)(2M/r0^3) * n nu tau (e^{2c}-1)/c / (M/r0) = n nu tau (e^{2c}-1)/(c r0^2)
R['gapV_upper_over_lagrangian']=n*(math.exp(2*c)-1)/c
R['gapV_upper_over_true']=R['gapV_upper_over_lagrangian']*TH/I2
c734=TH   # gap-V uses c = 0.734 (Gamma(0) tau)
R['gapV_quoted_factor_at_c_0p734']=n*(math.exp(2*c734)-1)/c734
print(f"(ii) true affine s_C / (nu tau/r0^2)          = I2/theta_max = {R['true_over_lagrangian']:.6f}")
print(f"     prove-lagrangian's pricing / true         = {R['lagrangian_over_true']:.6f}  (conservative)")
print(f"     gap-V's operator-norm upper bound / prove-lagrangian's pricing, at c={c:.6f}: "
      f"{R['gapV_upper_over_lagrangian']:.4f}")
print(f"     ... at gap-V's own c = Gamma(0)tau = {c734:.6f}: {R['gapV_quoted_factor_at_c_0p734']:.4f}"
      f"   (gap-V sec 5(iii) quotes 6.2 with a (1+f)^2)")
print(f"     gap-V's upper bound / the true affine value = {R['gapV_upper_over_true']:.4f}")
json.dump(R,open('x8_results.json','w'),indent=1)
print("WROTE x8_results.json")
