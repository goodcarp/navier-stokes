#!/usr/bin/env python3
"""Check arithmetic in quantitative-core-stage.md; no PDE simulation."""

from fractions import Fraction as F
from math import comb
import sympy as sp


CE = 65536
PRODUCT = {4: 512, 6: 2048, 8: 16384}


def bounds(viscosity, velocity_bound):
    V = velocity_bound
    M1 = 3 * viscosity * V + PRODUCT[8] * V * V
    M2 = (3 * viscosity + 2 * PRODUCT[6] * V) * M1
    M3 = ((3 * viscosity + 2 * PRODUCT[4] * V) * M2
          + 2 * PRODUCT[4] * M1 * M1)
    return M1, M2, M3


def run():
    assert 9 * (2**10 - 1)**2 * comb(13, 3) < CE**2
    for r, C in PRODUCT.items():
        assert 9 * 2**(2 * r) * comb(r + 3, 3) < C**2
    print("PASS: explicit Sobolev product and energy counting constants.")

    S, lam, V = sp.symbols("S lam V", positive=True)
    normal = bounds(lam, V)
    scaled = bounds(S * lam, S * V)
    for degree, original, result in zip((2, 3, 4), normal, scaled):
        assert sp.expand(result - S**degree * original) == 0
    print("PASS: time-derivative bounds scale as S^2, S^3, S^4.")

    # Test the interval restrictions and integral remainders using exact rationals.
    # Inputs below are positive model norm/margin values, not computed bump norms.
    for nu in [F(0), F(1, 10), F(3)]:
        for X0 in [F(1), F(5), F(20)]:
            for omega in [F(1, 4), F(2)]:
                for k0 in [F(1, 10), F(1), F(7)]:
                    M1, M2, M3 = bounds(nu, 2 * X0)
                    TE = 1 / (4 * CE * X0)
                    tau = min(TE, k0 / M2, omega * k0 / M3)
                    assert tau > 0
                    assert X0 / (1 - CE * X0 * TE) == F(4, 3) * X0
                    for fraction in [F(0), F(1, 3), F(1, 2), F(1)]:
                        t = fraction * tau
                        q_lower = k0 * t - M2 * t*t / 2
                        h_prime_lower = 2 * omega * k0 * t - M3 * t*t
                        h_taylor_lower = (2 * omega + omega * k0 * t*t
                                          - M3 * t**3 / 3)
                        desired = 2 * omega * (1 + k0 * t*t / 4)
                        modal_lower = 1 + k0 * t*t / 2 - M3 * t**3 / (4 * omega)
                        assert q_lower >= k0 * t / 2
                        assert h_prime_lower >= omega * k0 * t
                        assert h_taylor_lower >= desired
                        assert modal_lower >= 1 + k0 * t*t / 4
                    gain = k0 * tau*tau / 4
                    q_squared = (1 + gain / 2) / (1 + gain)
                    assert 0 < q_squared < 1
                    assert q_squared * (1 + gain) == 1 + gain / 2
    print("PASS: interval restrictions imply strain, vorticity, modal and defined Reynolds gains.")

    # Verify the derivative-norm H2 comparison polynomial coefficientwise.
    x, y, z = sp.symbols("x y z", real=True)
    radius_squared = x*x + y*y + z*z
    norm_symbol = (1 + radius_squared + x**4 + y**4 + z**4
                   + x*x*y*y + x*x*z*z + y*y*z*z)
    assert sp.expand(2 * norm_symbol - (1 + radius_squared)**2) == (
        1 + x**4 + y**4 + z**4)
    print("PASS: order-two derivative norm dominates the stated Bessel norm.")
    print("All arithmetic checks passed. No numerical PDE evolution was performed.")


if __name__ == "__main__":
    run()
