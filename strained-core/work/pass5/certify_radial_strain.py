#!/usr/bin/env python3
"""Validated range enclosure of the radial meridional cubic coefficient.

The formula is the pass4 exact radial identity, applied homogeneously to
Phi=psi+c eta with Phi(0)=1. Every panel uses outward-rounded interval ranges;
the accumulated Q,K integrals are also intervals, including on each panel.
No numerical pressure solve or floating-point convergence bound is used.
"""
import argparse,json,time
from fractions import Fraction as F
from mpmath import iv
from interval_local_integrals import interval,cutoff_jets,exact_endpoints


def certify(c=F(7,5),n=1024,dps=35):
    iv.dps=dps
    C=interval(c); zero=interval(0)
    # In [0,1/4], Phi=1 and Phi'=0 exactly.
    Q=interval(F(1,4)); K=zero
    coefficient=interval(F(340,49))
    started=time.monotonic()
    jets=[]
    for j in range(n):
        jets.append(cutoff_jets(interval(F(1,4)+F(3*j,4*n),F(1,4)+F(3*(j+1),4*n)),1))

    def transition(left,right,kind):
        nonlocal Q,K,coefficient
        h=(right-left)/n
        H=interval(h); partial=interval(0,h)
        for j,(f,fp) in enumerate(jets):
            s=interval(left+j*h,left+(j+1)*h)
            if kind=='core':
                phi,der=f,fp
            elif kind=='inner':
                phi,der=-C*(1-f),C*interval(F(4,25))*fp
            elif kind=='outer':
                phi,der=-C*f,-C*interval(F(3,44))*fp
            else:
                raise ValueError(kind)
            ks=s**(iv.mpf(7)/2)*(phi*der+interval(F(7,9))*s*der**2)
            qrange=Q+partial*phi
            krange=K+partial*ks
            integrand=(-interval(F(18368,735))*s*phi*der**2
                       -interval(F(448,105))*s**2*der**3
                       -interval(F(256,105))*qrange*der**2
                       +interval(F(1536,49))*s**(-iv.mpf(7)/2)*der*krange)
            coefficient+=H*integrand
            Q+=H*phi
            K+=H*ks
        return dict(kind=kind,coefficient=str(coefficient),Q=str(Q),K=str(K))

    stages=[transition(F(1,4),F(1),'core')]
    # Phi=0 between s=1 and25/16: no increment in Q or K.
    stages.append(transition(F(25,16),F(25,4),'inner'))
    # Phi=-c on [25/4,25]. K and the cubic integrand vanish there.
    Q-=C*interval(F(25)-F(25,4))
    stages.append(transition(F(25),F(36),'outer'))
    return dict(c=str(c),panels_per_transition=n,precision_decimal_digits=dps,
                coefficient_interval=str(coefficient),exact_dyadic_bounds=exact_endpoints(coefficient),
                stages=stages,seconds=time.monotonic()-started,
                method='outward-rounded range integration with cumulative interval inclusions',
                scope='meridional radial cubic only; full NS feedback and iteration separate')


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--c',default='7/5')
    ap.add_argument('--panels',type=int,default=1024)
    ap.add_argument('--dps',type=int,default=35)
    args=ap.parse_args()
    print(json.dumps(certify(F(args.c),args.panels,args.dps)),flush=True)
