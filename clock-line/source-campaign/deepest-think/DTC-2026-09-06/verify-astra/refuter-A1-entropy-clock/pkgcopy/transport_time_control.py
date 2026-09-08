#!/usr/bin/env python3
"""Exact nonautonomous scalar-kernel orientation control.

For v_t + t*v_x1 - nu*Delta(v)=0,
K(t,x;s,y)=G_{nu*(t-s)}(y-x+((t*t-s*s)/2)*e1).
At fixed terminal t, tau=t-s gives the Gaussian density with mean
x-(t*tau-tau*tau/2)*e1 and covariance 2*nu*tau*I.
Its unit mass follows from Gaussian normalization. Its elapsed-time
drift is -(t-tau)*e1, so its density equation has +(t-tau)*p_y1.

This spatially constant drift has infinite physical kinetic energy.
It is an auxiliary scalar kernel control, not a finite-energy NS example.
No root implementation or other project code is read or imported.
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def run():
    t, tau, nu = sp.symbols("t tau nu", positive=True)
    x = sp.symbols("x1:4", real=True)
    y = sp.symbols("y1:4", real=True)
    displacement = t*tau - tau**2/2
    q = (y[0]-x[0]+displacement, y[1]-x[1], y[2]-x[2])
    p = (4*sp.pi*nu*tau)**(-sp.Rational(3, 2)) * sp.exp(
        -sum(component**2 for component in q)/(4*nu*tau))
    time_derivative = sp.diff(p, tau)
    laplacian = sum(sp.diff(p, coordinate, 2) for coordinate in y)
    axial_derivative = sp.diff(p, y[0])

    def normalized_residual(drift_coefficient):
        # Equation being tested: p_tau=nu*Delta(p)+coefficient*p_y1.
        return sp.factor(sp.simplify((time_derivative-nu*laplacian
                                     -drift_coefficient*axial_derivative)/p))

    correct = normalized_residual(t-tau)
    wrong_sign = normalized_residual(-(t-tau))
    wrong_elapsed_time = normalized_residual(tau)
    wrong_frozen_time = normalized_residual(t)
    assert correct == 0, correct
    assert sp.simplify(wrong_sign + (t-tau)*q[0]/(nu*tau)) == 0
    assert sp.simplify(wrong_elapsed_time
                       + (t-2*tau)*q[0]/(2*nu*tau)) == 0
    assert sp.simplify(wrong_frozen_time-q[0]/(2*nu)) == 0

    # One strictly interior point: terminal t=3, elapsed tau=1, source s=2.
    point = {t: 3, tau: 1, nu: 1,
             x[0]: 0, x[1]: 0, x[2]: 0,
             y[0]: 1, y[1]: 0, y[2]: 0}
    witness = {
        "correct": correct.subs(point),
        "wrong_sign": wrong_sign.subs(point),
        "wrong_elapsed_time_argument": wrong_elapsed_time.subs(point),
        "wrong_frozen_terminal_time": wrong_frozen_time.subs(point),
    }
    assert witness == {
        "correct": 0, "wrong_sign": -7,
        "wrong_elapsed_time_argument": -sp.Rational(7, 4),
        "wrong_frozen_terminal_time": sp.Rational(7, 4),
    }, witness

    # Check that the source-time displacement becomes the expression used.
    s = sp.symbols("s", real=True)
    assert sp.expand(((t*t-s*s)/2).subs(s, t-tau)-displacement) == 0

    source = Path(__file__).resolve()
    return {
        "status": "pass",
        "source": str(source),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "sympy_version": sp.__version__,
        "domain": "nu>0; 0<tau<t; all x,y in physical R^3",
        "kernel": "G_{nu*tau}(y-x+(t*tau-tau^2/2)e1)",
        "correct_equation": "p_tau=nu*Delta_y p+(t-tau)*partial_y1 p",
        "normalized_residuals": {
            "correct": str(correct),
            "wrong_sign": str(wrong_sign),
            "wrong_elapsed_time_argument": str(wrong_elapsed_time),
            "wrong_frozen_terminal_time": str(wrong_frozen_time),
        },
        "exact_interior_witness": {name: str(value)
                                   for name, value in witness.items()},
        "scope": "Auxiliary scalar transport kernel with infinite-energy drift; no finite-energy NS or new lifespan claim.",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path)
    options = parser.parse_args()
    rendered = json.dumps(run(), indent=2) + "\n"
    if options.receipt:
        options.receipt.parent.mkdir(parents=True, exist_ok=True)
        options.receipt.write_text(rendered)
    print(rendered, end="")
