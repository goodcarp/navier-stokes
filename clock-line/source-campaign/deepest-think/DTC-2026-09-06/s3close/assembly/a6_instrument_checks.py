"""
a6 -- instrument checks of my Theorem-V.4 sizing pipeline against the V-b refuter's
      published numbers at THEIR design point, and of the feedback exponent algebra.

Refuter's design point (write/refute-V-b-bulk-viscous-loss/PROOF.md sec.1, sec.3):
  f = 0, phi0 = 30 deg, delta = 7.5 deg, d = rho0 sin(delta), R_- = r0 - d,
  window theta_max = 4(1 - sqrt(2/3)), Grownwall constant c_G = sqrt(3/2)*theta_max,
  rho0 sin(delta) = sqrt(nu/M), n = 5, hatK2 = 1, lambda(theta) = (1 - kappa theta/2)^-2,
  kappa = 1/2, I2 = 4/5 - 16 sqrt6/135, I4 = -4/7 + 27 sqrt6/28.
Their published values:  E_4 = 1.9523e-2, Eh1 = 1.1572e-2, Eh2 = 3.7089e-1,
  Eh2/Eh1 = 32.050, (Eh1+Eh2)/(seat's) = 87.166, eps_bulk(L=10) = 3.4735e-3.
"""
import json, math
OUT = {}

theta_max = 4.0*(1.0 - math.sqrt(2.0/3.0))
c_G = math.sqrt(1.5)*theta_max
I2 = 4.0/5.0 - 16.0*math.sqrt(6.0)/135.0
I4 = -4.0/7.0 + 27.0*math.sqrt(6.0)/28.0
delta = 7.5*math.pi/180.0
s = 1.0/math.sin(delta)
phi0 = 30.0*math.pi/180.0
n = 5.0
rho0 = 1.0; M = 1.0
r0 = math.sin(phi0)              # f = 0
d = math.sin(delta)              # = rho0 sin delta = sqrt(nu/M)
R_m = r0 - d
N = r0/math.sin(delta)
hatK2 = 1.0

OUT["theta_max"] = theta_max; OUT["c_G"] = c_G; OUT["I2"] = I2; OUT["I4"] = I4
OUT["s"] = s; OUT["r0"] = r0; OUT["d"] = d; OUT["R_minus"] = R_m; OUT["N"] = N

rows = {}
for L in [10.0, 40.0, 160.0, 640.0]:
    nu = M*rho0**2/s**2
    tau = theta_max/(M*L)
    sy = nu*I2/(M*L)               # nu * int lam^-2 dt
    sz = nu*I4/(M*L)
    eps_bulk = sy/r0**2
    V1 = n*nu*tau*(math.exp(2*c_G)-1.0)/c_G
    K2 = hatK2*M/rho0
    EW = 0.5*K2*math.exp(c_G)*tau*V1
    Eh1 = EW*(r0/R_m**2)
    Eh2 = 4.0*N*EW/d
    trC = 8*sy + 2*sz
    trC2 = 4*(2*sy)**2 + (2*sz)**2
    E4 = (r0/R_m**5)*(trC**2 + 2*trC2)
    seat = 0.5*math.exp(c_G)*n*(math.exp(2*c_G)-1.0)/c_G*K2*tau*r0   # PROOF.md sec.5 formula
    rows["L=%g" % L] = {"eps_bulk": eps_bulk, "Eh1": Eh1, "Eh2": Eh2,
                        "Eh2_over_Eh1": Eh2/Eh1, "E_4": E4,
                        "seat_sized_E_hess": seat*eps_bulk,
                        "correction_factor": (Eh1+Eh2)/(seat*eps_bulk)}
OUT["rows"] = rows
OUT["published"] = {"E_4_L10": 1.9523e-2, "Eh1_L10": 1.1572e-2, "Eh2_L10": 3.7089e-1,
                    "Eh2_over_Eh1": 32.050, "correction_factor": 87.166,
                    "eps_bulk_L10": 3.4735e-3, "theta_max_over_I2": 1.440118}
p = OUT["published"]; r = rows["L=10"]
OUT["relerr"] = {"E_4": abs(r["E_4"]-p["E_4_L10"])/p["E_4_L10"],
                 "Eh1": abs(r["Eh1"]-p["Eh1_L10"])/p["Eh1_L10"],
                 "Eh2": abs(r["Eh2"]-p["Eh2_L10"])/p["Eh2_L10"],
                 "Eh2_over_Eh1": abs(r["Eh2_over_Eh1"]-p["Eh2_over_Eh1"])/p["Eh2_over_Eh1"],
                 "correction_factor": abs(r["correction_factor"]-p["correction_factor"])/p["correction_factor"],
                 "eps_bulk": abs(r["eps_bulk"]-p["eps_bulk_L10"])/p["eps_bulk_L10"],
                 "theta_max_over_I2": abs(theta_max/I2-p["theta_max_over_I2"])/p["theta_max_over_I2"]}
OUT["theta_max_over_I2"] = theta_max/I2

# ---- feedback exponent algebra, symbolic-then-numeric --------------------------
import json as _j
A1 = _j.load(open("a1_results.json"))
kappa = A1["kappa_delta"]["7.5"]; r_h = A1["r_h"]["7.5"]
geom_hi, geom_lo = 1.5, 4.0/9.0
p_inf = 1.5 + geom_hi*2.0*kappa*(15.0*math.pi/4.0)/(r_h*geom_lo)
cstar = math.log(1.5)/kappa
sens = 0.34913942333376546/2.1181705106873974
OUT["feedback"] = {"p_over_ML_proved_Linf": p_inf, "pc_proved": p_inf*cstar,
                   "exp_pc_proved": math.exp(p_inf*cstar),
                   "p_over_ML_sharp_Linf": 1.5 + (p_inf-1.5)*sens,
                   "pc_sharp": (1.5 + (p_inf-1.5)*sens)*cstar,
                   "exp_pc_sharp": math.exp((1.5 + (p_inf-1.5)*sens)*cstar),
                   "sens_factor": sens, "c_star": cstar}
# ---- small-delta deficit law  1 - 2 kappa_delta ~ delta^3/4 --------------------
lawchk = {}
for dd in [3.0, 5.0, 7.5, 10.0]:
    dr = dd*math.pi/180.0
    lawchk["%g" % dd] = {"one_minus_2kappa": 1.0-2.0*A1["kappa_delta"]["%g" % dd],
                         "delta_cubed_over_4": dr**3/4.0}
OUT["small_delta_law"] = lawchk

with open("a6_results.json", "w") as f:
    json.dump(OUT, f, indent=1, sort_keys=True)
print(json.dumps({"relerr": OUT["relerr"], "feedback": OUT["feedback"],
                  "small_delta_law": OUT["small_delta_law"],
                  "rows_L10": rows["L=10"]}, indent=1, sort_keys=True))
