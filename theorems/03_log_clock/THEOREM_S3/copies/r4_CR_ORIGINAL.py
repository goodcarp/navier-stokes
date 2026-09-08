"""Independent re-derivation of the velocity remainder constant C_R of u1/PROOF.md sec.3,
keeping the geometric factors u1 drops:

   R_y = (a(X) - frak_a(|X|)) Y      => |R_y| <= C_a(phi) M_s * (rho sin phi)     [NOT rho]
   R_z = -INT_0^z 2[a - frak_a(rho_z)] - INT_0^z 2[frak_a(rho_z)-frak_a(rho)] - INT_0^z r d_r a

Exact auxiliary integrals (sympy-checked below):
   INT_0^|z| rho_zeta dzeta = (|z| rho + r^2 log((|z|+rho)/r))/2
   INT_0^|z| log(rho/rho_zeta) dzeta = |z| - r arctan(|z|/r)      <=  |z|      (NO log(1/sin))
Collar constant: min( 3.999218/sin phi ,  6.281958/delta )   -- the second is the
taper-repaired collar of write/L3v-and-gamma-bound sec.B5 (PROVED there, 47.99 at delta=7.5deg).
"""
import math
import numpy as np
import sympy as sp

# --- verify the two auxiliary integrals exactly
r, zz, ze = sp.symbols('r z zeta', positive=True)
rho = sp.sqrt(r**2 + zz**2); rz = sp.sqrt(r**2 + ze**2)
I1 = sp.integrate(rz, (ze, 0, zz)); I1c = (zz*rho + r**2*sp.log((zz+rho)/r))/2
print("INT rho_zeta  residual:", sp.simplify(I1 - I1c))
I2 = sp.integrate(sp.log(rho/rz), (ze, 0, zz)); I2c = zz - r*sp.atan(zz/r)
print("INT log(rho/rho_zeta) residual:", sp.simplify(sp.expand_log(sp.simplify(I2 - I2c), force=True)))

CF, CI, CC, CAX = 0.291999, 0.014754, 3.999218, math.pi/8
G1 = 0.5 + 1.6049285 + 0.1324254          # G_collar = 0
DELTA = 7.5*math.pi/180.0
COLLAR_TAPER = 6.281958/DELTA             # 47.99, L3v sec.B5

def C_a(sphi, taper=False):
    coll = CC/sphi
    if taper: coll = min(coll, COLLAR_TAPER)
    return CF + CI + coll + CAX

def CR_u1(phi, taper=False):
    """u1's own assembly: 3 C_a + log(1/sin) + G1 (all geometric factors dropped)."""
    s = math.sin(phi)
    return 3*C_a(s, taper) + math.log(1.0/s) + G1

def CR_sharp(phi, taper=False):
    """same three pieces, geometric factors KEPT, exact zeta-integrals."""
    s, c = math.sin(phi), abs(math.cos(phi))       # rho = 1
    rr, z = s, c
    Ry = C_a(s, taper)*s
    # term 1 : 2 INT_0^z C_a(phi_zeta) dzeta,  sin phi_zeta = r/rho_zeta
    n = 20000
    zg = np.linspace(0.0, z, n+1)
    rz = np.sqrt(rr**2 + zg**2)
    coll = CC*rz/rr
    if taper: coll = np.minimum(coll, COLLAR_TAPER)
    t1 = 2*np.trapz(CF + CI + coll + CAX, zg)
    # term 2 : 2*(M/2)*INT log(rho/rho_zeta) = z - r arctan(z/r)
    t2 = z - rr*math.atan(z/rr) if rr > 0 else z
    # term 3 : G1 * z
    t3 = G1*z
    return Ry + t1 + t2 + t3

print("\n phi   C_a(phi)   C_R u1      C_R sharp   C_R sharp+taper")
for d in [90,60,45,30,20,15,10,7.5,5,2,1,0.5,0.2]:
    p = d*math.pi/180
    print("%5.1f  %9.4f  %10.4f  %10.4f  %10.4f" %
          (d, C_a(math.sin(p)), CR_u1(p), CR_sharp(p), CR_sharp(p, True)))

gr = np.linspace(1e-5, math.pi/2, 4000)
print("\nsup over ALL phi of C_R sharp+taper =", max(CR_sharp(p, True) for p in gr))
print("sup over the cone [delta,pi/2] of C_R sharp (no taper repair) =",
      max(CR_sharp(p) for p in gr if p >= DELTA))
print("sup_phi C_a(phi) sin phi =", max(C_a(math.sin(p))*math.sin(p) for p in gr),
      " (= C_a(90deg) =", C_a(1.0), ")")

# what these do to mu(c)L and hence to L_*
LM = math.exp(0.75*0.8113825936219724)
def muL(CR): return LM**3*(CR + 2*math.log(LM))*(LM**2 - 1.0)
for nm, CR in [("u1 headline  phi0=30", CR_u1(math.pi/6)),
               ("u1 cone      phi=delta", CR_u1(DELTA)),
               ("sharp        phi0=30", CR_sharp(math.pi/6)),
               ("sharp cone   sup>=delta", max(CR_sharp(p) for p in gr if p >= DELTA)),
               ("sharp+taper  sup all phi", max(CR_sharp(p, True) for p in gr))]:
    print("%-26s C_R=%9.4f  mu(c)L=%9.3f  L_* scale vs headline = %.4f"
          % (nm, CR, muL(CR), muL(CR)/muL(CR_u1(math.pi/6))))
