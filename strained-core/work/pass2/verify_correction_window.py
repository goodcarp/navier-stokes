#!/usr/bin/env python3
"""Exact and asymptotic checks for correction-window.md; no PDE existence test.

Uses only the standard library. The scalar recurrence checked here is the
displayed finite source majorant, not imported source code or a Lean proof.
"""

from fractions import Fraction as F
from math import comb, factorial, log


CORE_RE_EXP = F(1, 4)
CARRIER_EXP = F(1, 16)
TIME_EXP = 2 + CORE_RE_EXP
FORCE_EXP = 3 + 2 * CORE_RE_EXP


def depth(order):
    return 38 * order + 87


def physical_exponent(levels, spatial, temporal):
    return (CARRIER_EXP * (levels + 1 - spatial - 2 * temporal)
            - FORCE_EXP - spatial - TIME_EXP * temporal)


def catalan(n):
    return comb(2 * n, n) // (n + 1)


def finite_source_recursion(E, H, D, J):
    """Actual finite-depth scalar recursion, with coefficients zero above J."""
    a = [F(0)] * (2 * J + 5)
    b = [F(0)] * (2 * J + 5)
    s = [F(0)] * (2 * J + 5)
    a[0] = F(E)

    def at(array, n):
        return array[n] if 0 <= n < len(array) else F(0)

    b[0] = F(E)
    for n in range(1, 2 * J + 5):
        conv1 = sum(a[i] * (a[n - 1 - i] + b[n - 1 - i])
                    for i in range(n))
        conv2 = sum(a[i] * (a[n - 2 - i] + b[n - 2 - i])
                    for i in range(n - 1))
        s[n] = H * (6 * at(a, n - 1) + 3 * at(a, n - 2)
                    + 3 * conv1 + 3 * conv2)
        if n <= J:
            a[n] = 4 * D * s[n]
        b[n] = a[n] + H * (at(a, n - 1) + at(a, n - 2))
    return a, b, s


def run():
    assert TIME_EXP == F(9, 4)
    assert FORCE_EXP == F(7, 2)
    assert CORE_RE_EXP - 2 * CARRIER_EXP == F(1, 8)

    for K in range(101):
        J = depth(K)
        assert K + 4 * J + 13 == 153 * K + 361
        assert physical_exponent(J, 0, K) == 2
        for k in range(K + 1):
            for j in range(K + 1 - k):
                exact = F(J + 1 - 56 - 17 * k - 38 * j, 16)
                assert physical_exponent(J, k, j) == exact
                assert exact >= 2
                # A slow mean with this normalized bound also gives rho^2.
                mean_power = F(11, 2) + F(9, 4) * K
                assert mean_power - FORCE_EXP - k - TIME_EXP * j >= 2
    print("PASS: exact dimensional exponents and depth J=38K+87.")

    for E, H, D in [(1, 1, 1), (2, 3, 5), (5, 2, 3),
                     (F(3, 2), F(5, 4), F(7, 3))]:
        B = 256 * D * E * H * H
        assert 4 * D * H * (9 + 6 * E * (1 + 3 * H)) <= B
        for J in [0, 1, 2, 5, 9]:
            a, b, source = finite_source_recursion(E, H, D, J)
            for n in range(0, len(a)):
                assert a[n] <= E * catalan(n) * B ** n
                assert b[n] <= 3 * H * E * catalan(n) * B ** n
                if n > 0:
                    assert source[n] <= E * (4 * B) ** n
            delta = F(1, 16) / B
            assert 4 * B * delta <= F(1, 2)
            remainder = sum(source[g] * delta ** g
                            for g in range(J + 1, 2 * J + 5))
            assert remainder <= 2 * E * (4 * B * delta) ** (J + 1)
    print("PASS: finite scalar recursion, Catalan coefficient bound, terminal sum.")

    # Exact Reynolds-loss thresholds. Equality needs constants, not this test.
    for p, s, expected in [(F(0), F(1), True),
                            (F(1, 8), F(1), True),
                            (F(1, 4), F(1), False),
                            (F(1, 3), F(1), False)]:
        theta = s * CARRIER_EXP - p * CORE_RE_EXP
        assert (theta > 0) == expected
        assert (theta > 0) == (p < s / 4)
    assert F(1, 3) < F(1, 2)  # A different M could fit p/s=1/3.
    assert not F(3, 4) < F(1, 2)
    print("PASS: proposed p/s<1/4 and general p/s<1/2 power-gap tests.")

    # A model order constant, NOT measured source/PDE constants.
    # Use m=n^3, so floor(m^(1/3))=n exactly. Work in log space.
    print("Illustrative constant-budget ratios for log A=(J+1)(rF+1)log(rF+2):")
    ratios = []
    for n in [100, 1000, 10000, 100000, 1000000, 10000000]:
        m = n ** 3
        J = depth(n)
        rF = 153 * n + 361
        log_A = (J + 1) * (rF + 1) * log(rF + 2)
        ratio = log_A / (m * log(2))
        ratios.append(ratio)
        # Polynomial-in-J harmonic overhead, e.g. n_harm <= J^2.
        log_max_damping_ratio = 4 * log(J) - (m / 8) * log(2)
        assert log_max_damping_ratio < 0
        print(f"  K={n:g}, m=K^3: log(A)/log(1/rho)={ratio:.6g}")
    assert all(a > b for a, b in zip(ratios, ratios[1:]))
    assert ratios[-1] < 1
    print("PASS: illustrative constants eventually fit; polynomial harmonics stay below ceiling.")

    # An explicit nonzero Gaussian profile has nonzero even derivatives at 0.
    # If a mean of size rho^(2c) survives with this profile, a sufficiently high
    # derivative diverges. This is a force-regularity model, not a PDE example.
    for c in [F(0), F(1), F(10), F(100)]:
        n = 1
        while 2 * n <= 2 * c - FORCE_EXP:
            n += 1
        k = 2 * n
        assert 2 * c - FORCE_EXP - k < 0
        gaussian_derivative_nonzero = (-1) ** n * factorial(2 * n) // factorial(n)
        assert gaussian_derivative_nonzero != 0
    print("PASS: fixed-power nonzero mean floor fails some high derivative in the model.")
    print("All checks passed. No Navier–Stokes solution or full-force estimate was constructed.")


if __name__ == "__main__":
    run()
