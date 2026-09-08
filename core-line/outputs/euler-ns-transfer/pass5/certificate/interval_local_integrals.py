#!/usr/bin/env python3
"""Outward-rounded range quadrature for the actual smooth cutoff and local integrals.

Uses mpmath.iv elementary arithmetic/exponential/square-root enclosures. Every
leaf contributes area times a range enclosure; no sampled-derivative or
floating-point convergence estimate is used as an error bound. Exact dyadic
endpoint tuples are emitted alongside human-readable intervals.
"""
import argparse
import heapq
import itertools
import json
import time
from fractions import Fraction as F
from mpmath import iv


def interval(lo, hi=None):
    if hi is None and hasattr(lo,'_mpi_'):
        return lo
    if hi is None:
        hi = lo
    def endpoint(x):
        x = F(x)
        return iv.mpf(x.numerator)/x.denominator
    return iv.mpf([endpoint(lo).a, endpoint(hi).b])


def hull(xs):
    return iv.mpf([min(x.a for x in xs), max(x.b for x in xs)])


def logistic(t):
    return 1/(1+iv.exp(-1/t+1/(1-t)))


def cutoff_jets(X, order=2):
    """Enclose f,f',f'' with derivatives in the original cutoff argument.

    f(x)=1 for x<=1/4, =0 for x>=1; transition t=(4x-1)/3.
    Endpoint tails use exp(-1/t) polynomial bounds, never infinite derivatives.
    """
    zero, one = iv.mpf(0), iv.mpf(1)
    quarter, delta = interval(F(1,4)), interval(F(1,32))
    pieces = []
    if X.a <= quarter:
        pieces.append([one]+[zero]*order)
    if X.b >= one:
        pieces.append([zero]*(order+1))
    if X.b <= quarter or X.a >= one:
        return tuple(hull([p[j] for p in pieces]) for j in range(order+1))
    ta = max(zero, (4*X.a-1)/3)
    tb = min(one, (4*X.b-1)/3)
    bounds = [zero, delta, one-delta, one]
    for k in range(3):
        lo, hi = max(ta,bounds[k]), min(tb,bounds[k+1])
        if lo > hi:
            continue
        T = iv.mpf([lo.a, hi.b])
        if k in (0,2):
            # Reflect the upper endpoint tail into 0<=tau<=delta.
            reflected = (k == 2)
            tau = 1-T if reflected else T
            top = tau.b
            if top == 0:
                vals = [one]+[zero]*order
            else:
                low_value = logistic(top)
                high_value = one if tau.a == 0 else logistic(tau.a)
                vals = [iv.mpf([low_value.a,high_value.b])]
                if order:
                    e = iv.exp(-1/top+1/(1-top))
                    m1 = (iv.mpf(4)/3)*e*(top**-2+(1-top)**-2)
                    vals.append(iv.mpf([-m1.b,zero]))
                if order >= 2:
                    m2 = (iv.mpf(16)/9)*e*((top**-2+(1-top)**-2)**2
                                          +2*top**-3+2*(1-top)**-3)
                    vals.append(iv.mpf([-m2.b,m2.b]))
            if reflected:
                vals[0] = 1-vals[0]
                if order >= 2:
                    vals[2] = -vals[2]
            pieces.append(vals)
        else:
            lower, upper = logistic(T.b), logistic(T.a)
            p = iv.mpf([lower.a,upper.b])
            vals = [p]
            if order:
                h = p*(1-p)
                # p(1-p) never exceeds 1/4; retaining this bound improves range.
                h = iv.mpf([max(zero,h.a),min(quarter,h.b)])
                z1 = -T**-2-(1-T)**-2
                vals.append((iv.mpf(4)/3)*h*z1)
            if order >= 2:
                z2 = 2*T**-3-2*(1-T)**-3
                vals.append((iv.mpf(16)/9)*h*((1-2*p)*z1*z1+z2))
            pieces.append(vals)
    return tuple(hull([p[j] for p in pieces]) for j in range(order+1))


def exact_endpoints(value):
    def encode(t):
        sign,mantissa,exponent,bits = t
        return dict(sign=int(sign),mantissa=str(mantissa),exponent=int(exponent))
    return dict(lower=encode(value._mpi_[0]),upper=encode(value._mpi_[1]))


def range_quadrature(fn, partitions, max_panels=4096, relative_width=.1):
    """Validated 1D or tensor-box range integration on rational partitions.

    Adaptive priorities and the stopping rule are heuristics; they do not enter
    enclosure validity. Memory is linear in max_panels.
    """
    serial = itertools.count()
    heap = []
    total_lo, total_hi = iv.mpf(0), iv.mpf(0)
    base_lengths = [F(p[-1])-F(p[0]) for p in partitions]
    def leaf(box):
        value = fn(*[interval(a,b) for a,b in box])
        volume = F(1)
        for a,b in box:
            volume *= b-a
        value *= interval(volume)
        return (-float(value.delta),next(serial),box,value)
    intervals = [[(F(a),F(b)) for a,b in zip(p[:-1],p[1:])] for p in partitions]
    for box in itertools.product(*intervals):
        entry = leaf(box)
        heapq.heappush(heap,entry)
        total_lo += entry[3].a
        total_hi += entry[3].b
    while len(heap) < max_panels:
        enclosure = iv.mpf([total_lo.a,total_hi.b])
        scale = max(abs(float(enclosure.a)),abs(float(enclosure.b)),1.e-300)
        if float(enclosure.delta) <= relative_width*scale:
            break
        _,_,box,old = heapq.heappop(heap)
        axis = max(range(len(box)),key=lambda k: (box[k][1]-box[k][0])/base_lengths[k])
        midpoint = (box[axis][0]+box[axis][1])/2
        children = []
        for segment in [(box[axis][0],midpoint),(midpoint,box[axis][1])]:
            child = list(box); child[axis] = segment
            entry = leaf(tuple(child)); children.append(entry)
            heapq.heappush(heap,entry)
        # Subtract endpoints, not the old interval (which would inflate widths).
        total_lo = total_lo-old.a+sum((e[3].a for e in children),iv.mpf(0))
        total_hi = total_hi-old.b+sum((e[3].b for e in children),iv.mpf(0))
    result = iv.mpf([total_lo.a,total_hi.b])
    return dict(interval=str(result),exact_dyadic_bounds=exact_endpoints(result),
                panels=len(heap),strictly_positive=bool(result.a>0),
                strictly_negative=bool(result.b<0),
                width_display=float(result.delta))


def swirl_integrand(kind, radial_scale=F(63,200), axial_scale=F(9,20), distance=F(4)):
    """C_v or d_v for w=r f(r²/sr²) f((z-d)²/sz²), mirrored at -d.

    Coordinates xi=r²/sr² in [0,1], eta=(z-d)/sz in [-1,1].
    The full 3D Jacobian, including both packets, is 2 pi sr² sz.
    """
    sr,sz,d = interval(radial_scale),interval(axial_scale),interval(distance)
    jacobian = 2*iv.pi*sr**2*sz
    def integrand(xi,eta):
        order = 0 if kind == 'Cv' else 1
        radial = cutoff_jets(xi,order)
        axial = cutoff_jets(eta**2,order)
        q = sr**2*xi
        z = d+sz*eta
        R2 = q+z*z
        Qtt = 3*(4*z*z-q)/(4*iv.pi*R2**(iv.mpf(7)/2))
        W = radial[0]*axial[0]
        if kind == 'Cv':
            return jacobian*Qtt*q*W**2
        Qrr = (12*R2**2-105*q*z*z)/(4*iv.pi*R2**(iv.mpf(9)/2))
        wr = (radial[0]+2*xi*radial[1])*axial[0]
        wz_over_r = radial[0]*axial[1]*(2*eta/sz)
        return jacobian*2*(Qrr*W**2+Qtt*(wr**2+q*wz_over_r**2))
    return integrand


def pump_integrand(kind, sigma=F(7,10),height=F(1),distance=F(4),ratio=F(9,20),
                   pump_amplitude=F(1),swirl_amplitude=F(1)):
    """Range functions for Cp, dP and the local outer convective contribution.

    Domain xi=r²/sigma² in [0,1], eta=(z-d)/height in [-1,1].
    'local_pzz_convection' is +2 int Q_ij u_i ((u.grad)u)_j, i.e. the
    contribution from acceleration -(u.grad)u to pzz'. It omits the separate
    nonlocal pressure-acceleration contribution and is not the full derivative.
    """
    sig,h,d,tau = map(interval,(sigma,height,distance,ratio))
    c,A = interval(pump_amplitude),interval(swirl_amplitude)
    jacobian = 2*iv.pi*sig**2*h
    def integrand(xi,eta):
        alpha,a1,a2 = cutoff_jets(xi,2)
        beta,b1,b2 = cutoff_jets(eta**2,2)
        q,z = sig**2*xi,d+h*eta
        bz = b1*2*eta/h
        bzz = (2*b1+4*eta**2*b2)/h**2
        K = beta+z*bz
        H,Hq,Hz = alpha*K,a1*K/sig**2,alpha*(2*bz+z*bzz)
        Z = -2*z*(alpha+xi*a1)*beta
        Zr_over_r = -4*z*(2*a1+xi*a2)*beta/sig**2
        Zz = -2*(alpha+xi*a1)*K
        R2 = q+z*z
        den = 4*iv.pi*R2**(iv.mpf(9)/2)
        Qrr = (12*R2**2-105*q*z*z)/den
        Qtt = (15*z*z-3*R2)/(4*iv.pi*R2**(iv.mpf(7)/2))
        Qzz = -(105*z**4-90*R2*z*z+9*R2**2)/den
        Qrz_over_r = -(105*z**3-45*R2*z)/den
        if kind == 'Cp':
            return jacobian*(Qrr*q*H**2+Qzz*Z**2+2*Qrz_over_r*q*H*Z)
        if kind == 'dP':
            vr = H+2*q*Hq
            return jacobian*2*(Qrr*(vr**2+q*Hz**2)+Qtt*H**2
                +Qzz*(q*Zr_over_r**2+Zz**2)
                +2*Qrz_over_r*q*(vr*Zr_over_r+Hz*Zz))
        radial = cutoff_jets(xi/tau**2,1)
        axial = cutoff_jets(eta**2/tau**2,1)
        W = A*radial[0]*axial[0]
        Wq = A*radial[1]*axial[0]/(sig*tau)**2
        Wz = A*radial[0]*axial[1]*2*eta/(h*tau**2)
        H,Hq,Hz,Z,Zq,Zz = c*H,c*Hq,c*Hz,c*Z,c*Zr_over_r/2,c*Zz
        nr = H*(H+2*q*Hq)+Z*Hz-W**2
        nt = 2*H*(W+q*Wq)+Z*Wz
        nz = 2*q*H*Zq+Z*Zz
        return jacobian*2*(Qrr*q*H*nr+Qtt*q*W*nt+Qzz*Z*nz
                           +Qrz_over_r*q*(H*nz+Z*nr))
    if kind not in ('Cp','dP','local_pzz_convection'):
        raise ValueError(kind)
    return integrand


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--kind',choices=['Cv','dv','Cp','dP','local_pzz_convection'],default='Cv')
    parser.add_argument('--panels',type=int,default=4096)
    parser.add_argument('--relative-width',type=float,default=.1)
    parser.add_argument('--dps',type=int,default=35)
    args = parser.parse_args()
    iv.dps = args.dps
    partitions = [[F(0),F(1,4),F(1,2),F(3,4),F(1)],
                  [F(-1),F(-3,4),F(-1,2),F(0),F(1,2),F(3,4),F(1)]]
    started = time.monotonic()
    fn = swirl_integrand(args.kind) if args.kind in ('Cv','dv') else pump_integrand(args.kind)
    result = range_quadrature(fn,partitions,args.panels,args.relative_width)
    result.update(kind=args.kind,precision_decimal_digits=args.dps,
                  seconds=time.monotonic()-started,
                  profile='actual smooth logistic cutoff; sigma=7/10, h=1, ratio=9/20, d=4',
                  method='outward-rounded range inclusion on every rectangle; exact dyadic bounds emitted')
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()
