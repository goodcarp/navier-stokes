"""
f2_datum_s3 -- THEOREM_S3's OWN datum as a field object, in the interface hk2's instrument
               expects (hk2lib.DatumB / StrainedB), so that hk2's k3 assembly can be run on it
               without editing anything in s3close/hk2.

THE TWO DATA ARE NOT THE SAME OBJECT.

  hk2 (D-B), the campaign datum, radial ramp in RHO:
      Theta_hk2(rho) = (1/2)[ tanh((rho-rho_0)/w_0) - tanh((rho-R)/w_1) ] ,  w_0 = 0.25 rho_0
      angular:  tanh(sin phi / sin delta) * tanh(cos phi / w) ,  w = 0.20      (real-analytic)

  THEOREM_S3 sec.1.1, radial ramp in LOG RHO:
      Theta_S3(rho)  = (1/2)[ tanh((u-eps_r)/eps_r) - tanh((u-L+eps_r)/eps_r) ] , u = log(rho/rho_0)
      angular:  sgn(cos phi) min(1, phi_ax/delta) min(1, |phi-pi/2|/delta_m)     (Lipschitz, kinked)

Units: rho_0 = 1, M = 1.  eta = omega^theta/r = G(rho) W(phi),
      G(rho) = -M Theta(rho)/rho ,   W(phi) = A(phi)/sin phi ,
      A(phi) = sgn(cos phi) min(1, phi_ax/delta) min(1, |phi - pi/2|/delta_m) .

A is piecewise LINEAR in phi (A'' = 0 on every piece) and is C^1 across the equator
phi = pi/2 (both one-sided slopes are -1/delta_m), so the only kinks are at
phi = delta, pi/2 - delta_m, pi/2 + delta_m, pi - delta.

Derivatives are taken in (rho, phi) and pushed to (r,z) by the exact chain rule -- NOT in the
t = cos phi variable hk2lib uses, because W' = (A' sin phi - A cos phi)/sin^2 phi suffers
cancellation on the axis in that variable.
"""
import math
import numpy as np

DEG = math.pi/180.0


class DatumS3:
    def __init__(self, delta_deg=7.5, dm_deg=5.0, eps_r=0.25, L=10.0, M=1.0):
        self.delta = delta_deg*DEG
        self.dm = dm_deg*DEG
        self.eps_r = eps_r
        self.L = L
        self.R = math.exp(L)
        self.M = M
        self.delta_deg = delta_deg
        self.dm_deg = dm_deg

    # ---------------- radial profile, in u = log rho
    def Th_u(self, u, k=0):
        e = self.eps_r
        a = (u - e)/e
        b = (u - self.L + e)/e
        if k == 0:
            return 0.5*(np.tanh(a) - np.tanh(b))
        ca = 1.0/np.cosh(np.clip(a, -350, 350))**2
        cb = 1.0/np.cosh(np.clip(b, -350, 350))**2
        if k == 1:                                   # dTheta/du
            return 0.5*(ca - cb)/e
        ta, tb = np.tanh(a), np.tanh(b)
        if k == 2:                                   # d^2Theta/du^2
            return (tb*cb - ta*ca)/e**2
        raise ValueError

    def Th(self, rho, k=0):
        """Theta and its rho-derivatives (k = 0,1,2)."""
        u = np.log(rho)
        T0 = self.Th_u(u, 0)
        if k == 0:
            return T0
        T1 = self.Th_u(u, 1)
        if k == 1:
            return T1/rho
        T2 = self.Th_u(u, 2)
        if k == 2:
            return (T2 - T1)/rho**2
        raise ValueError

    # ---------------- angular profile A(phi) and W = A/sin phi
    def A_of_phi(self, phi, k=0):
        s_sign = np.where(np.cos(phi) >= 0.0, 1.0, -1.0)
        pax = np.minimum(phi, math.pi - phi)
        dpax = np.where(phi < math.pi/2, 1.0, -1.0)
        eqd = np.abs(phi - math.pi/2)
        deqd = np.where(phi > math.pi/2, 1.0, -1.0)
        g1 = np.minimum(1.0, pax/self.delta)
        dg1 = np.where(pax < self.delta, dpax/self.delta, 0.0)
        g2 = np.minimum(1.0, eqd/self.dm)
        dg2 = np.where(eqd < self.dm, deqd/self.dm, 0.0)
        if k == 0:
            return s_sign*g1*g2
        if k == 1:
            return s_sign*(dg1*g2 + g1*dg2)
        if k == 2:
            return np.zeros_like(np.asarray(phi, dtype=float))   # A is piecewise linear
        raise ValueError

    def W(self, phi, k=0):
        s = np.sin(phi)
        c = np.cos(phi)
        A0 = self.A_of_phi(phi, 0)
        if k == 0:
            return A0/s
        A1 = self.A_of_phi(phi, 1)
        if k == 1:
            return (A1*s - A0*c)/s**2
        if k == 2:                       # A'' = 0
            return -2.0*A1*c/s**2 + A0/s + 2.0*A0*c*c/s**3
        raise ValueError

    # ---------------- eta and its (r,z) derivatives
    def eta_rz(self, r, z):
        r = np.asarray(r, dtype=float)
        z = np.asarray(z, dtype=float)
        rho = np.hypot(r, z)
        phi = np.arctan2(np.maximum(r, 1e-300), z)
        s, c = np.sin(phi), np.cos(phi)
        M = self.M
        Th0 = self.Th(rho, 0)
        Tu1 = self.Th_u(np.log(rho), 1)
        Tu2 = self.Th_u(np.log(rho), 2)
        G0 = -M*Th0/rho
        G1 = -M*(Tu1 - Th0)/rho**2
        G2 = -M*(Tu2 - 3.0*Tu1 + 2.0*Th0)/rho**3
        W0 = self.W(phi, 0)
        W1 = self.W(phi, 1)
        W2 = self.W(phi, 2)
        rr_, rz_ = s, c
        pr_, pz_ = c/rho, -s/rho
        rrr, rrz, rzz = c*c/rho, -s*c/rho, s*s/rho
        prr = -2.0*s*c/rho**2
        prz = (s*s - c*c)/rho**2
        pzz = 2.0*s*c/rho**2
        eta = G0*W0
        eta_r = G1*rr_*W0 + G0*W1*pr_
        eta_z = G1*rz_*W0 + G0*W1*pz_
        eta_rr = (G2*rr_**2*W0 + G1*rrr*W0 + 2*G1*rr_*W1*pr_ + G0*W2*pr_**2 + G0*W1*prr)
        eta_rz = (G2*rr_*rz_*W0 + G1*rrz*W0 + G1*(rr_*W1*pz_ + rz_*W1*pr_)
                  + G0*W2*pr_*pz_ + G0*W1*prz)
        eta_zz = (G2*rz_**2*W0 + G1*rzz*W0 + 2*G1*rz_*W1*pz_ + G0*W2*pz_**2 + G0*W1*pzz)
        return eta, eta_r, eta_z, eta_rr, eta_rz, eta_zz

    def grad_norm(self, r, z):
        _, er, ez, *_ = self.eta_rz(r, z)
        return np.sqrt(er*er + ez*ez)

    def grad_norm_fast(self, r, z):
        return self.grad_norm(r, z)

    def hess_F(self, r, z):
        _, er, ez, err, erz, ezz = self.eta_rz(r, z)
        return np.sqrt(err**2 + 2*erz**2 + ezz**2 + 3*(er/np.maximum(r, 1e-300))**2)

    def hess_op(self, r, z):
        _, er, ez, err, erz, ezz = self.eta_rz(r, z)
        tr = err + ezz
        disc = np.sqrt(np.maximum((err-ezz)**2 + 4*erz**2, 0.0))
        lam = np.maximum(np.abs(0.5*(tr+disc)), np.abs(0.5*(tr-disc)))
        return np.maximum(lam, np.abs(er/np.maximum(r, 1e-300)))


class StrainedS3:
    """eta_lam = eta_0 o T_lam^{-1},  T_lam(y,z) = (lam y, lam^{-2} z)."""
    def __init__(self, base, lam):
        self.b = base
        self.lam = float(lam)
        self.M = base.M
        self.R = base.R

    def _pull(self, r, z):
        return np.asarray(r, dtype=float)/self.lam, (self.lam**2)*np.asarray(z, dtype=float)

    def eta_rz(self, r, z):
        lam = self.lam
        e, er, ez, err, erz, ezz = self.b.eta_rz(*self._pull(r, z))
        return (e, er/lam, ez*lam**2, err/lam**2, erz*lam, ezz*lam**4)

    def grad_norm(self, r, z):
        _, er, ez, *_ = self.eta_rz(r, z)
        return np.sqrt(er*er + ez*ez)

    def grad_norm_fast(self, r, z):
        return self.grad_norm(r, z)

    def hess_F(self, r, z):
        _, er, ez, err, erz, ezz = self.eta_rz(r, z)
        return np.sqrt(err**2 + 2*erz**2 + ezz**2 + 3*(er/np.maximum(r, 1e-300))**2)


def dist_to_cone(r, z, ang):
    """distance in the (r,z) half plane to the cone phi = ang about the +z axis."""
    return abs(r*math.cos(ang) - z*math.sin(ang))


def geometry_s3(phi0_deg=30.0, f=1.0, delta_deg=7.5, dm_deg=5.0, lam=1.0):
    """the tracked point ((1+f)rho_0, phi_0) pushed by T_lam, and its distances to the FOUR
    kink cones of the theorem's datum (delta, pi/2 - delta_m, pi/2 + delta_m, pi - delta),
    each transported by the same law tan -> lam^3 tan.  The radial ramp is smooth, so unlike
    hk2's (D-A) it imposes NO exclusion."""
    phi0 = phi0_deg*DEG
    rho0 = 1.0 + f
    r0, z0 = rho0*math.sin(phi0), rho0*math.cos(phi0)
    r, z = lam*r0, z0/lam**2
    rho = math.hypot(r, z)
    cones = {}
    for name, ang in [("taper", delta_deg*DEG), ("eq_lo", math.pi/2 - dm_deg*DEG),
                      ("eq_hi", math.pi/2 + dm_deg*DEG), ("taper_mirror", math.pi - delta_deg*DEG)]:
        ta = math.tan(ang)
        angl = math.atan2(lam**3*ta, 1.0) if ta >= 0 else math.pi + math.atan(lam**3*ta)
        cones[name] = {"cone_deg": math.degrees(angl), "dist": dist_to_cone(r, z, angl)}
    dmin = min(v["dist"] for v in cones.values())
    return {"r": r, "z": z, "rho": rho, "phi_deg": math.degrees(math.atan2(r, z)),
            "lam": lam, "f": f, "cones": cones, "d_layers": dmin}


if __name__ == "__main__":
    # self-test: the analytic (rho,phi) chain rule against 4th-order central differences,
    # scaled by |eta|/rho^k as hk2's k2_datum.py does (a relative test at the plateau values
    # of the second derivatives is meaningless: they are O(1e-5) there).
    import json
    D = DatumS3(delta_deg=7.5, dm_deg=5.0, eps_r=0.25, L=10.0)
    def f(rr, zz):
        return float(np.asarray(D.eta_rz(np.array([rr]), np.array([zz]))[0]).reshape(-1)[0])
    def d1(g, h):
        return (-g(2*h) + 8*g(h) - 8*g(-h) + g(-2*h))/(12*h)
    def d2(g, h):
        return (-g(2*h) + 16*g(h) - 30*g(0.0) + 16*g(-h) - g(-2*h))/(12*h*h)
    rows, worst = [], 0.0
    for (r, z) in [(0.5, 0.866), (2.5, 4.330), (1.2, 0.3), (0.9, -1.4), (1.5, 2.0), (0.2, 0.05)]:
        e, er, ez, err, erz, ezz = [float(np.asarray(v).reshape(-1)[0])
                                    for v in D.eta_rz(np.array([r]), np.array([z]))]
        rho = math.hypot(r, z)
        h = 1e-3
        fr = d1(lambda t: f(r+t, z), h)
        fz = d1(lambda t: f(r, z+t), h)
        frr = d2(lambda t: f(r+t, z), h)
        fzz = d2(lambda t: f(r, z+t), h)
        frz = d1(lambda t: d1(lambda v_: f(r+t, z+v_), h), h)
        scale = abs(e)
        ana = [er, ez, err, erz, ezz]
        num = [fr, fz, frr, frz, fzz]
        sc = [scale/rho, scale/rho, scale/rho**2, scale/rho**2, scale/rho**2]
        rel = max(abs(a-n)/max(abs(n), 1e-3*s_) for a, n, s_ in zip(ana, num, sc))
        worst = max(worst, rel)
        rows.append({"r": r, "z": z, "scaled_worst_rel": rel})
        print("(r,z)=(%.2f,%.3f)  worst scaled err = %.3e" % (r, z, rel))
    out = {"points": rows, "worst_scaled_rel": worst,
           "note": "scaling convention copied from hk2/k2_datum.py section 1"}
    json.dump(out, open("f2_datum_check.json", "w"), indent=1)
    print("worst =", worst)
