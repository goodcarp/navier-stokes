#!/usr/bin/env python3
"""Exact algebra/constant checks for the moving-cylinder a priori lemma.

This does not verify a PDE solution, existence, a NS bootstrap, or S3.
It reads the named live proof only to record its current revision hash.
"""
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sympy as s


def main():
    checks = []
    def exact(label, expression):
        assert s.simplify(expression) == 0, (label, expression)
        checks.append(label)

    d, n, r, nu, L, h = s.symbols('d n r nu L h', positive=True)
    zeta = 1-d
    exact('cutoff drift combination', 2*zeta**2+4*zeta*d-2*(1-d**2))
    # For radial f(d), Delta f = (4d f_dd + 2n f_d)/r^2.
    lap_zeta2 = (4*d*s.diff(zeta**2,d,2)+2*n*s.diff(zeta**2,d))/r**2
    grad_zeta_sq = 4*d/r**2
    exact('diffusion cutoff constant', -lap_zeta2+8*grad_zeta_sq
          -(4*n*(1-d)+24*d)/r**2)
    Z, H, G, Q = s.symbols('Z H G Q', real=True)
    exact('Young inequality perfect square',
          2*Z**2*H**2+8*G**2*Q**2-8*Z*G*H*Q-2*(Z*H-2*G*Q)**2)
    for dimension in range(1,13):
        cn = max(4*dimension,24)
        # C_n - linear interpolant is a nonnegative combination of 1-d,d.
        exact(f'cutoff maximum n={dimension}',
              cn-(4*dimension*(1-d)+24*d)
              -((cn-4*dimension)*(1-d)+(cn-24)*d))
        assert cn-4*dimension >= 0 and cn-24 >= 0
    assert max(4*5,24) == 24 and 1+2+24 == 27
    checks.append('five-dimensional constant 27')
    # The spatial/time coefficient bound follows termwise on 0<=d<=1,
    # 0<=sigma<=h. These exact positive decompositions certify each step.
    exact('one minus cutoff square', 1-zeta**2-d*(2-d))
    A = (1+2*L*h+s.Max(4*n,24)*nu*h/r**2)/(2*nu)
    exact('choice of Bernstein coefficient',
          2*nu*A-(1+2*L*h+s.Max(4*n,24)*nu*h/r**2))
    # Minkowski's first weight inequality follows after squaring:
    # (sqrt(1+a^2)+c)^2 - [1+(a+c)^2]
    # = 2c/(sqrt(1+a^2)+a), which is nonnegative.
    a,c = s.symbols('a c', nonnegative=True)
    exact('tail weight comparison',
          (s.sqrt(1+a*a)+c)**2-(1+(a+c)**2)
          -2*c/(s.sqrt(1+a*a)+a))
    # Integral of Xi_9 on R^5 via the beta integral in radial coordinates.
    beta_integral = s.gamma(s.Rational(5,2))*s.gamma(2)/(2*s.gamma(s.Rational(9,2)))
    sphere_area = 2*s.pi**s.Rational(5,2)/s.gamma(s.Rational(5,2))
    exact('five-dimensional tail integral',sphere_area*beta_integral-16*s.pi**2/105)
    # Exact bounded solution for a time-dependent linear drift; alpha need
    # only be continuous, so its time-Hölder modulus is not manufactured.
    x,t = s.symbols('x t', real=True)
    At,St = s.Function('A')(t),s.Function('Sigma')(t)
    al = s.Function('alpha')(t)
    u = s.exp(-nu*St)*s.cos(s.exp(-At)*x)
    residual = s.diff(u,t)+al*x*s.diff(u,x)-nu*s.diff(u,x,2)
    exact('linear-drift example PDE',residual.subs({s.diff(At,t):al,
                                               s.diff(St,t):s.exp(-2*At)}))
    # Initial-trace extension: exact finite algebra only. The note separately
    # proves the pathwise comparison and states the representation premise.
    G0,B0,D0 = s.symbols('G0 B0 D0', nonnegative=True)
    exact('quadratic weighted-gradient closure',
          G0+2*s.sqrt(G0*B0)*D0+B0*D0**2-(s.sqrt(G0)+s.sqrt(B0)*D0)**2)
    exact('Brownian quadratic moment',
          2*s.gamma(s.Rational(5+2,2))/s.gamma(s.Rational(5,2))-5)
    exact('Brownian tenth moment',
          2**5*s.gamma(s.Rational(5+10,2))/s.gamma(s.Rational(5,2))-45045)
    assert 45045 < 3**10
    checks.append('c_5_10 less than 10/3')
    # q=2 Doob factor and sigma=sqrt(2nu): (2sqrt(nt)sigma)^2=8n nu t.
    exact('quadratic noise constant', (2*s.sqrt(n*t)*s.sqrt(2*nu))**2-8*n*nu*t)
    source = Path('~/Desktop/Solve Navier Stokes/campaign/'
        'deepest-think/DTC-2026-09-06/s3close/round2/pmax-h3v/PROOF.md')
    pin = None
    if source.is_file():
        data=source.read_bytes()
        pin=dict(path=str(source),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),
                 mtime_utc=datetime.fromtimestamp(source.stat().st_mtime,timezone.utc).isoformat())
    result=dict(status='PASS',checks=len(checks),labels=checks,
                scope='Exact cutoff, Young, weight, tail-integral, Gaussian-moment and example-PDE algebra only; '
                      'not existence, solution regularity, NS bootstrap or S3 certification.',
                source_revision_observed=pin)
    out=Path(__file__).with_name('comoving-gradient-decay-check.json')
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
