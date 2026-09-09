"""
r0_iv -- outward-rounded interval arithmetic on numpy arrays, plus truncated-Taylor series whose
coefficients are polynomials in the two direction parameters (a, c).

Written from scratch for the round-3 referee sitting (L-14: no code shared with fix5/x1, whose
arithmetic is plain float and whose search is a grid).  Outward rounding is one ulp for +,-,*,/
(a single IEEE operation is correct to half an ulp) and four ulps for sqrt/log/exp/arctan/tanh
(libm accuracy on this platform is well inside four ulps for these functions).
"""
import math
import numpy as np


def _dn(x, n=1):
    for _ in range(n):
        x = np.nextafter(x, -np.inf)
    return x


def _up(x, n=1):
    for _ in range(n):
        x = np.nextafter(x, np.inf)
    return x


class IV(object):
    """interval [lo, hi], componentwise over a numpy array"""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, dtype=float)
        hi = lo if hi is None else np.asarray(hi, dtype=float)
        self.lo = lo
        self.hi = hi

    def mag(self):
        return np.maximum(np.abs(self.lo), np.abs(self.hi))

    def __repr__(self):
        return "IV[%r, %r]" % (self.lo, self.hi)


def iv(x, y=None):
    return IV(x, y)


def const(x, shape):
    a = np.full(shape, float(x))
    return IV(a, a)


def add(A, B):
    return IV(_dn(A.lo + B.lo), _up(A.hi + B.hi))


def sub(A, B):
    return IV(_dn(A.lo - B.hi), _up(A.hi - B.lo))


def neg(A):
    return IV(-A.hi, -A.lo)


def mul(A, B):
    p1 = A.lo*B.lo
    p2 = A.lo*B.hi
    p3 = A.hi*B.lo
    p4 = A.hi*B.hi
    lo = np.minimum(np.minimum(p1, p2), np.minimum(p3, p4))
    hi = np.maximum(np.maximum(p1, p2), np.maximum(p3, p4))
    return IV(_dn(lo), _up(hi))


def smul(s, A):
    """scalar (python float) times interval"""
    if s >= 0:
        return IV(_dn(s*A.lo), _up(s*A.hi))
    return IV(_dn(s*A.hi), _up(s*A.lo))


def inv(A):
    if np.any((A.lo <= 0.0) & (A.hi >= 0.0)):
        raise ValueError("interval contains zero in inv")
    lo = np.minimum(1.0/A.lo, 1.0/A.hi)
    hi = np.maximum(1.0/A.lo, 1.0/A.hi)
    return IV(_dn(lo), _up(hi))


def div(A, B):
    return mul(A, inv(B))


def ipow(A, n):
    """integer power of an interval; A must be positive for odd/even alike here"""
    if n == 0:
        return IV(np.ones_like(A.lo), np.ones_like(A.hi))
    out = A
    for _ in range(n-1):
        out = mul(out, A)
    return out


def sqrt(A):
    if np.any(A.lo < 0.0):
        raise ValueError("sqrt of a negative interval")
    return IV(_dn(np.sqrt(A.lo), 4), _up(np.sqrt(A.hi), 4))


def log(A):
    if np.any(A.lo <= 0.0):
        raise ValueError("log of a non-positive interval")
    return IV(_dn(np.log(A.lo), 4), _up(np.log(A.hi), 4))


def exp(A):
    return IV(_dn(np.exp(A.lo), 4), _up(np.exp(A.hi), 4))


def arctan(A):
    return IV(_dn(np.arctan(A.lo), 4), _up(np.arctan(A.hi), 4))


def tanh(A):
    return IV(_dn(np.tanh(A.lo), 4), _up(np.tanh(A.hi), 4))


def sin_mono(A):
    """sin on an interval contained in [0, pi/2]: monotone increasing"""
    return IV(_dn(np.sin(A.lo), 4), _up(np.sin(A.hi), 4))


def cos_mono(A):
    """cos on an interval contained in [0, pi/2]: monotone decreasing"""
    return IV(_dn(np.cos(A.hi), 4), _up(np.cos(A.lo), 4))


# ------------------------------------------------------------------ polynomials in (a, c)
NORD = 5                       # Taylor orders s^0 .. s^4
MON = []
for d in range(NORD):
    for i in range(d, -1, -1):
        MON.append((i, d - i))
IDX = {m: k for k, m in enumerate(MON)}
DEG = [i + j for (i, j) in MON]
NMON = len(MON)                # 15

# max over {a^2 + c^2 <= 1} of |a|^i |c|^j
MAXMON = []
for (i, j) in MON:
    if i + j == 0:
        MAXMON.append(1.0)
    else:
        k = i + j
        v = (i**i if i else 1)*(j**j if j else 1)/float(k**k)
        MAXMON.append(math.sqrt(v))


def pzero():
    return [None]*NMON


def pconst(A):
    p = pzero()
    p[0] = A
    return p


def padd(P, Q):
    out = pzero()
    for k in range(NMON):
        if P[k] is None:
            out[k] = Q[k]
        elif Q[k] is None:
            out[k] = P[k]
        else:
            out[k] = add(P[k], Q[k])
    return out


def psub(P, Q):
    out = pzero()
    for k in range(NMON):
        if Q[k] is None:
            out[k] = P[k]
        elif P[k] is None:
            out[k] = neg(Q[k])
        else:
            out[k] = sub(P[k], Q[k])
    return out


def pmul(P, Q):
    out = pzero()
    for ka in range(NMON):
        A = P[ka]
        if A is None:
            continue
        for kb in range(NMON):
            B = Q[kb]
            if B is None:
                continue
            if DEG[ka] + DEG[kb] > NORD - 1:
                raise ValueError("polynomial degree overflow: the grading is wrong")
            kk = IDX[(MON[ka][0] + MON[kb][0], MON[ka][1] + MON[kb][1])]
            t = mul(A, B)
            out[kk] = t if out[kk] is None else add(out[kk], t)
    return out


def pscale(s, P):
    out = pzero()
    for k in range(NMON):
        if P[k] is not None:
            out[k] = smul(s, P[k])
    return out


def pmulIV(A, P):
    out = pzero()
    for k in range(NMON):
        if P[k] is not None:
            out[k] = mul(A, P[k])
    return out


def pmag_disc(P):
    """upper bound on |P(a,c)| over the closed disc a^2 + c^2 <= 1"""
    tot = None
    for k in range(NMON):
        if P[k] is None:
            continue
        t = P[k].mag()*MAXMON[k]
        tot = t if tot is None else _up(tot + t)
    return tot


def peval_box(P, A, C):
    """interval evaluation of P at interval a in A, c in C"""
    tot = None
    for k in range(NMON):
        if P[k] is None:
            continue
        i, j = MON[k]
        t = P[k]
        if i:
            t = mul(t, ipow_signed(A, i))
        if j:
            t = mul(t, ipow_signed(C, j))
        tot = t if tot is None else add(tot, t)
    if tot is None:
        return IV(np.zeros_like(A.lo), np.zeros_like(A.lo))
    return tot


def ipow_signed(A, n):
    """integer power of an interval that may straddle zero"""
    if n == 0:
        return IV(np.ones_like(A.lo), np.ones_like(A.lo))
    if n % 2 == 1:
        return IV(_dn(A.lo**n), _up(A.hi**n))
    m = np.maximum(np.abs(A.lo), np.abs(A.hi))
    lo = np.where((A.lo <= 0) & (A.hi >= 0), 0.0, np.minimum(np.abs(A.lo), np.abs(A.hi))**n)
    return IV(_dn(lo), _up(m**n))


# ------------------------------------------------------------------ series (list of NORD polys)
def szero():
    return [pzero() for _ in range(NORD)]


def sadd(S, T):
    return [padd(S[k], T[k]) for k in range(NORD)]


def ssub(S, T):
    return [psub(S[k], T[k]) for k in range(NORD)]


def smul_series(S, T):
    out = szero()
    for i in range(NORD):
        if all(x is None for x in S[i]):
            continue
        for j in range(NORD - i):
            if all(x is None for x in T[j]):
                continue
            out[i+j] = padd(out[i+j], pmul(S[i], T[j]))
    return out


def scompose(fd, S):
    """f(S(s)) with fd[k] an interval enclosure of f^{(k)} at the base point S[0]."""
    du = [pzero()] + [S[k] for k in range(1, NORD)]
    out = szero()
    out[0] = pconst(fd[0])
    pw = szero()
    pw[0] = pconst(IV(np.ones_like(fd[0].lo), np.ones_like(fd[0].lo)))
    fac = 1.0
    for k in range(1, NORD):
        pw = smul_series(pw, du)
        fac *= k
        fk = smul(1.0/fac, fd[k])
        for m in range(NORD):
            if all(x is None for x in pw[m]):
                continue
            out[m] = padd(out[m], pmulIV(fk, pw[m]))
    return out


def s_sqrt(S):
    a = S[0][0]
    h = sqrt(a)
    d = [h,
         smul(0.5, inv(h)),
         smul(-0.25, inv(mul(a, h))),
         smul(0.375, inv(mul(ipow(a, 2), h))),
         smul(-0.9375, inv(mul(ipow(a, 3), h)))]
    return scompose(d, S)


def s_log(S):
    a = S[0][0]
    d = [log(a), inv(a), smul(-1.0, inv(ipow(a, 2))), smul(2.0, inv(ipow(a, 3))),
         smul(-6.0, inv(ipow(a, 4)))]
    return scompose(d, S)


def s_inv(S):
    a = S[0][0]
    d = [inv(a), smul(-1.0, inv(ipow(a, 2))), smul(2.0, inv(ipow(a, 3))),
         smul(-6.0, inv(ipow(a, 4))), smul(24.0, inv(ipow(a, 5)))]
    return scompose(d, S)


def s_arctan(S):
    a = S[0][0]
    one = IV(np.ones_like(a.lo), np.ones_like(a.lo))
    d1 = inv(add(one, ipow(a, 2)))
    d = [arctan(a),
         d1,
         smul(-2.0, mul(a, ipow(d1, 2))),
         mul(add(smul(6.0, ipow(a, 2)), smul(-2.0, one)), ipow(d1, 3)),
         mul(smul(24.0, sub(a, ipow(a, 3))), ipow(d1, 4))]
    return scompose(d, S)
