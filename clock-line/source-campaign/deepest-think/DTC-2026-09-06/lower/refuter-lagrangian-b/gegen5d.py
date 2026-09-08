#!/usr/bin/env python3
"""x2 -- MY OWN 5D Gegenbauer instrument for a(x) at an ARBITRARY point, applied to the
STRAINED configuration.  Tests the step the attempt did NOT test:
LEMMA 1 is proved only for a(0) (the origin, where there is no vorticity).  The integro-ODE
needs the strain at the MATERIAL POINT in the DEFORMED configuration.  Here I compute it.

eta = sum_l e_l(rho) C_l^{3/2}(t),  e_l = q_l(rho)/rho ;  Delta_5 psi = -eta ;
 psi_l(rho)  = rho/(2l+3) [ int_0^inf e^{(1-l)w} q_l(rho e^w) dw + int_{-inf}^0 e^{(l+4)w} q_l dw ]
 psi_l'(rho) =   1/(2l+3) [ l int_0^inf e^{(1-l)w} q_l dw - (l+3) int_{-inf}^0 e^{(l+4)w} q_l dw ]
 a = -d_z psi = -sum_l [ psi_l' t C_l + (psi_l/rho)(1-t^2) C_l' ].
Strained field: eta_lam(rho,t) = Xi(t)/rho on {rho_0 < rho G(t) < R},
 G = sqrt((1-t^2)/lam^2 + lam^4 t^2), t'=lam^2 t/G, Xi = -M lam sgn(t) h(phi') / sqrt(1-t^2).
"""
import numpy as np, json, sys
M = 1.0

def gegen_and_deriv(lmax, t):
    C = np.zeros((lmax+2, t.size)); D = np.zeros((lmax+2, t.size))
    # Gegenbauer alpha=3/2:  n C_n = 2 t (n+alpha-1) C_{n-1} - (n+2alpha-2) C_{n-2}
    C[0]=1.0; C[1]=3*t; D[0]=0.0; D[1]=3.0
    for n in range(2, lmax+1):
        C[n] = (2*t*(n+0.5)*C[n-1] - (n+1)*C[n-2])/n
        D[n] = (2*(n+0.5)*(t*D[n-1]+C[n-1]) - (n+1)*D[n-2])/n
    return C[:lmax+1], D[:lmax+1]
def Nl(l): return (l+1)*(l+2)/(l+1.5)

class Field:
    def __init__(self, lam, hfun, L, lmax=801, nt=20001, nu=8001, pad=4.0):
        self.lam, self.L, self.lmax = lam, L, lmax
        rho0, R = 1.0, np.exp(L)
        t = np.linspace(-1+1e-13, 1-1e-13, nt)
        G = np.sqrt((1-t**2)/lam**2 + lam**4*t**2)
        g = np.log(G)
        tp = np.clip(lam**2*t/G, -1.0, 1.0); php = np.arccos(tp)
        Xi = -M*lam*np.sign(t)*hfun(php)/np.sqrt(1-t**2)
        C, _ = gegen_and_deriv(lmax, t)
        wt = np.gradient(t)
        ls = np.arange(lmax+1)
        W = ((1-t**2)*Xi*wt)[None,:]*C / np.array([Nl(l) for l in ls])[:,None]   # (lmax+1, nt)
        # support in u=log rho:  -g(t) < u < L - g(t)   <=>   g(t) > -u  and  g(t) < L-u
        order = np.argsort(g); gs = g[order]; Ws = W[:, order]
        cum = np.concatenate([np.zeros((lmax+1,1)), np.cumsum(Ws, axis=1)], axis=1)  # (lmax+1, nt+1)
        umin = -gs.max() - pad; umax = L - gs.min() + pad
        u = np.linspace(umin, umax, nu)
        lo = np.searchsorted(gs, -u, side='right')      # count of g <= -u
        hi = np.searchsorted(gs, L-u, side='left')      # count of g <  L-u
        hi = np.maximum(hi, lo)
        q = cum[:, hi] - cum[:, lo]                     # (lmax+1, nu)
        self.t, self.u, self.du, self.q = t, u, u[1]-u[0], q
        self.Hbulk = cum[:, -1]                          # all-t coefficients (bulk)
    def a(self, rho_ev, t_ev, lmax=None):
        lmax = self.lmax if lmax is None else lmax
        u, du, q = self.u, self.du, self.q
        i0 = int(np.argmin(np.abs(u - np.log(rho_ev))))
        du_p = u[i0:] - u[i0]; du_m = u[:i0+1] - u[i0]
        Ce, De = gegen_and_deriv(lmax, np.array([t_ev])); Ce=Ce[:,0]; De=De[:,0]
        tot = 0.0; parts=[]
        for l in range(1, lmax+1, 2):
            Ip = np.sum(np.exp((1-l)*du_p)*q[l, i0:])*du
            Im = np.sum(np.exp((l+4)*du_m)*q[l, :i0+1])*du
            psil  = rho_ev/(2*l+3.0)*(Ip+Im)
            psilp = (l*Ip - (l+3)*Im)/(2*l+3.0)
            c = -(psilp*t_ev*Ce[l] + (psil/rho_ev)*(1-t_ev**2)*De[l])
            tot += c; parts.append(c)
        return tot, parts

h1 = lambda p: np.ones_like(p)
def htap(d): return lambda p: np.minimum(1.0, np.minimum(p, np.pi-p)/d)


