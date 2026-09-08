#!/usr/bin/env python3
"""
Audit of the Kim-Jeong (arXiv:2111.14078, JFA 283 (2022) 109673) dyadic ring
stack, purely as BOOKKEEPING of the quantities the paper itself defines.

Their data (eq. (2.1) / eq. (13) of Bang-Cheskidov):
    omega_0 = sum_{k=n0}^{m} k^{-alpha} phi( r/(1/8)^{k-1}, z/(1/8)^k ),  1/q < alpha < 1
so ring k sits at radius ~ 8^{-(k-1)} with amplitude k^{-alpha}.

Their Lemma 3.3 (p.8): scale-dependent window  T_k = min{T, c1 (1-alpha) k^{-1+alpha}},
on which the outer rings stay put and  I_n(t) ~ I_n(0) = C n^{-alpha}.
Their stretching rate at ring k is  a_k = (1/((d-1)|B_d|)) sum_{n<k} I_n  ~  sum_{n<k} n^{-alpha}.

Questions checked numerically (exact rational/float arithmetic on the stated formulas):
 (Q1) a_k / (M0 * log(R_outer/R_inner))  ->  what constant?   [M0 = ||omega_0||_inf = n0^{-alpha}]
 (Q2) a_k * T_k  ->  is it O(1)?  (i.e. is T_k the doubling time 1/a_k?)
 (Q3) For Prop 4.1 (p.10): M0*T(m) and A(m)*T(m) where
      T(m) = c1(1-alpha) m^{(-1+alpha)/2},  A(m) = (1/4) m^{c2-alpha}.
NO claim is made here about any PDE; this only re-derives the paper's own exponents.
"""
import math

def a_k(k, n0, alpha):
    """sum_{n=n0}^{k-1} n^{-alpha}  (their I_n(0) ~ C n^{-alpha}, C set to 1)."""
    return sum(n**(-alpha) for n in range(n0, k))

def report():
    LOG8 = math.log(8.0)
    print("Q1/Q2: ring-k stretching rate vs M0*log(R/rho0), and a_k*T_k")
    print(f"{'alpha':>7} {'n0':>4} {'k':>7} {'M0':>10} {'a_k':>12} "
          f"{'efolds':>8} {'a_k/(M0*efolds)':>16} {'a_k*T_k/c1':>11}")
    for alpha in (0.02, 0.05, 0.10, 0.25, 0.50):
        for n0 in (2,):
            for k in (10, 100, 1000, 10000):
                M0 = n0**(-alpha)                      # max initial amplitude
                A  = a_k(k, n0, alpha)                 # stretching rate at ring k
                # radii: ring j at 8^{-(j-1)};  outermost j=n0, innermost j=k
                efolds = (k - n0) * LOG8               # log(R_outer/R_inner)
                Tk = (1.0 - alpha) * k**(-1.0 + alpha) # = T_k / c1
                print(f"{alpha:>7.2f} {n0:>4d} {k:>7d} {M0:>10.4f} {A:>12.4f} "
                      f"{efolds:>8.1f} {A/(M0*efolds):>16.4f} {A*Tk:>11.4f}")
    print()
    print("Q3: Prop 4.1 (p.10) clock, c2 < 1/4, alpha < c2. T(m)=c1(1-a)m^{(a-1)/2}, A(m)=(1/4)m^{c2-a}")
    print(f"{'alpha':>7} {'c2':>6} {'m':>9} {'M0':>8} {'M0*T/c1':>12} {'A(m)':>12} {'A*T/c1':>12} {'growth A/M0':>12}")
    for alpha, c2 in ((0.05, 0.20), (0.10, 0.20), (0.20, 0.24)):
        for m in (10**3, 10**6, 10**9, 10**12):
            M0 = 2**(-alpha)
            T  = (1-alpha) * m**((alpha-1)/2.0)     # = T(m)/c1
            A  = 0.25 * m**(c2-alpha)
            print(f"{alpha:>7.2f} {c2:>6.2f} {m:>9d} {M0:>8.4f} {M0*T:>12.3e} "
                  f"{A:>12.4e} {A*T:>12.3e} {A/M0:>12.4e}")

if __name__ == "__main__":
    report()

def q4():
    """
    Q4: read Kim-Jeong Prop 4.1 (p.10) as a CLOCK.
    N = number of dyadic rings ~ m ; Lambda = log(R_outer/R_inner) = (m-n0)*log 8.
    Growth  G = A(m)/M0 ~ m^{c2} ;  time  T = c1(1-alpha) m^{(alpha-1)/2}.
    Effective exponent beta with  T = c/(M0 * Lambda^beta)  is beta = (1-alpha)/2
    (up to the log-8 conversion, checked numerically below by a log-log fit).
    Constraint of Prop 4.1: 0 < alpha < c2 and (from (3.10)) c2 < 1/4.
    """
    import math
    print()
    print("Q4: effective clock exponent beta in  T ~ c/(M0 * Lambda^beta)   [Lambda = log(R/rho0)]")
    print(f"{'alpha':>7} {'beta=(1-a)/2':>13} {'fit beta (m=1e6..1e12)':>24} {'growth exponent c2':>20}")
    LOG8 = math.log(8.0)
    for alpha, c2 in ((0.02, 0.24), (0.05, 0.20), (0.10, 0.20), (0.24, 0.245)):
        m1, m2 = 1e6, 1e12
        L1, L2 = (m1-2)*LOG8, (m2-2)*LOG8
        T1, T2 = (1-alpha)*m1**((alpha-1)/2), (1-alpha)*m2**((alpha-1)/2)
        beta_fit = -(math.log(T2)-math.log(T1))/(math.log(L2)-math.log(L1))
        print(f"{alpha:>7.2f} {(1-alpha)/2:>13.4f} {beta_fit:>24.4f} {c2:>20.3f}")
    print("  -> beta stays strictly below 1/2; the upper-bound clock c/(M0(1+log Re)) has beta = 1.")
    print("  -> at FIXED growth factor (say 2) Prop 4.1 needs only m = 2^(1/c2) = O(1) rings,")
    print("     so it gives NO Lambda-dependence at all for a fixed-factor doubling clock:")
    for c2 in (0.20, 0.24, 0.245):
        print(f"     c2={c2:.3f} -> m needed for growth factor 2 is 2^(1/c2) = {2**(1/c2):.1f} rings")
