#!/usr/bin/env python3
"""Small independent manufactured audit; no production endpoint is loaded."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import hashlib,json,math,time
import numpy as np
import sympy as sp
from scipy.linalg import lstsq
from compact_potential_reconstruction import RadialBasis,fit_mode,evaluate_increment,lower_banded,FastMAC

HERE=Path(__file__).resolve().parent

def manual_window(z,L=16.):
    z=np.asarray(z);t=(z*z-36)/(L*L/4-36)
    q=np.zeros_like(z,float);qp=q.copy();q[t<=0]=1
    mask=(t>0)&(t<1);a=t[mask]
    H=1/a-1/(1-a);v=1/(1+np.exp(-H))
    q[mask]=v
    qp[mask]=v*(1-v)*(-1/a**2-1/(1-a)**2)*2*z[mask]/(L*L/4-36)
    return q,qp

def analytic_radial(r,m,R):
    """Explicit global polynomial r^m*(1-r²/R²)^9, not spline evaluation."""
    r=np.asarray(r,float);s=r*r/R**2;inside=r<=R
    t=np.maximum(1-s,0);fac=1.
    if m:
        peak=m/(m+18);fac=1/(peak**(m/2)*(1-peak)**9)
    power=(r/R)**m*fac;lower=np.zeros_like(r) if m==0 else m/R*(r/R)**(m-1)*fac
    B=power*t**9;Dr=lower*t**9-18*r/R**2*power*t**8
    Lap=4*power/R**2*(-9*(m+1)*t**8+72*s*t**7)
    mb=lower*t**9
    return tuple(np.where(inside,x,0.) for x in (B,Dr,Lap,mb))

def analytic_velocity(r,z,m,R,k,origin,p,t,window):
    B,Dr,Lap,mb=analytic_radial(r,m,R)
    e=np.exp(1j*k[:,None]*(np.asarray(z)[None,:]-origin))
    P=p@e;Pz=(1j*k*p)@e;T=t@e
    u=np.stack((Dr[:,None]*Pz+1j*mb[:,None]*T,
                1j*mb[:,None]*Pz-Dr[:,None]*T,-Lap[:,None]*P))
    if window:
        q,qp=manual_window(z)
        u=q*u+np.stack((Dr[:,None]*qp*P,1j*mb[:,None]*qp*P,np.zeros_like(u[2])))
    return u

def weighted_relative(g,x,y):
    return float(np.sqrt(np.sum(g.mass[:,None]*abs(x-y)**2)/max(np.sum(g.mass[:,None]*abs(y)**2),1e-300)))

def manufactured_fit(m,nz):
    g=FastMAC(nr=80,nz=nz,R=4.,L=16.,mapping=4.,J=4,nu=.001);g.mapping=4.
    j=m//4;rng=np.random.default_rng(8130+m+nz)
    p=np.zeros(nz,complex);t=p.copy()
    selected=(0,1,-1,2,-3,nz//2,-nz//2+1)
    for kval in selected:
        where=np.flatnonzero(g.ki==kval)
        if len(where):
            p[where[0]]=.02*(rng.normal()+1j*rng.normal())
            t[where[0]]=.03*(rng.normal()+1j*rng.normal())
    p[~g.keep]=0;t[~g.keep]=0
    if m==0:
        p[0]=p[0].real;t[0]=t[0].real
        for kval in g.ki[g.ki>0]:
            pos=np.flatnonzero(g.ki==kval)[0];neg=np.flatnonzero(g.ki==-kval)[0]
            p[neg]=p[pos].conj();t[neg]=t[pos].conj()
    _,Dr,Lap,mb=analytic_radial(g.r,m,g.R)
    _,Df,_,mbf=analytic_radial(g.rf[1:-1],m,g.R)
    data=np.concatenate((1j*Df[:,None]*g.kz*p+1j*mbf[:,None]*t,
                         -mb[:,None]*g.kz*p-Dr[:,None]*t,-Lap[:,None]*p))
    basis,P,T,fit,backward=fit_mode(g,data,j,intervals=12)
    fit_error=weighted_relative(g,fit,data)
    # Arbitrary off-grid radii include the axis, joins and the exterior.
    rr=np.r_[0.,.013,.071,.27,.63,1.13,1.91,2.73,3.48,3.99,4.,4.2]
    zz=np.array([-9.,-7.4,-6.7,-2.1,0.,1.4,6.7,7.4,9.])
    eval_errors={}
    for win in (False,True):
        got=evaluate_increment(basis,P,T,g.kz,g.z[0],g.L,rr,zz,window=win)
        want=analytic_velocity(rr,zz,m,g.R,g.kz,g.z[0],p,t,win)
        eval_errors[str(win)]=float(np.linalg.norm(got-want)/max(np.linalg.norm(want),1e-300))
    # Independently stack the actual P/T velocity columns and use dense QR.
    _,Br,Lc,mc=basis.matrices(g.r);_,Bfr,_,mf=basis.matrices(g.rf[1:-1]);n=basis.n
    dense=[];outside_band=[];condition=[]
    for col in (0,1,nz-1):
        k=g.kz[col]
        A=np.zeros((g.nv,2*n),complex)
        A[g.slr,::2]=1j*k*Bfr;A[g.slr,1::2]=1j*mf
        A[g.slt,::2]=-k*mc;A[g.slt,1::2]=-Br;A[g.slz,::2]=-Lc
        Aw=np.sqrt(g.mass)[:,None]*A;bw=np.sqrt(g.mass)*data[:,col]
        x,_,rank,_=lstsq(Aw,bw,lapack_driver='gelsy')
        assert rank==2*n
        production=np.empty(2*n,complex);production[::2]=P[:,col];production[1::2]=T[:,col]
        dense.append(float(np.linalg.norm(Aw@(x-production))/max(np.linalg.norm(bw),1e-300)))
        G=Aw.conj().T@Aw
        assert np.max(abs(G.imag))<1e-12*max(1.,np.max(abs(G)))
        band=lower_banded(G.real);restored=np.zeros_like(G.real)
        for d in range(len(band)):
            idx=np.arange(2*n-d);restored[idx+d,idx]=band[d,:len(idx)]
            if d:restored[idx,idx+d]=band[d,:len(idx)]
        outside_band.append(float(np.linalg.norm(G-restored)/np.linalg.norm(G)))
        condition.append(float(np.linalg.cond(Aw)))
    # Endpoint indices: the final retained function has order exactly nine.
    endjets=[float(np.max(abs(basis.spline(g.R**2,nu=d)))) for d in range(9)]
    order9=float(abs(basis.spline(g.R**2,nu=9)[-1]))
    assert max(endjets)==0. and order9>0.
    out=dict(m=m,nz=nz,fit_weighted_relative=fit_error,offgrid_relative=eval_errors,
        dense_QR_velocity_relative=max(dense),band_truncation_relative=max(outside_band),
        design_condition_max=max(condition),normal_equation_backward=backward,
        outer_endpoint_derivatives_0_to_8=endjets,last_retained_derivative9=order9)
    assert fit_error<2e-8,out
    assert max(eval_errors.values())<2e-7,out
    assert max(dense)<2e-8,out
    assert max(outside_band)<1e-14,out
    return out

def symbolic_checks():
    x,y,z,r=sp.symbols('x y z r',real=True)
    I=sp.I;s=x*x+y*y
    def curlpot(P,T):return sp.Matrix([sp.diff(P,x,z)+sp.diff(T,y),sp.diff(P,y,z)-sp.diff(T,x),-sp.diff(P,x,2)-sp.diff(P,y,2)])
    results=[]
    for m in (-4,0,4,16):
        n=abs(m);H=(x+(I if m>=0 else -I)*y)**n
        P=H*(1+s)*(z+z*z);T=H*(2+s*s)*(1+I*z)
        u=curlpot(P,T);pr=r**n*(1+r*r)*(z+z*z);tr=r**n*(2+r**4)*(1+I*z)
        ur=sp.diff(pr,r,z)+I*m*tr/r;ut=I*m*sp.diff(pr,z)/r-sp.diff(tr,r)
        uz=-(sp.diff(pr,r,2)+sp.diff(pr,r)/r-m*m*pr/r**2)
        axis=[sp.simplify(v) for v in (ur+I*ut,ur-I*ut,uz)]
        for actual,expected in zip(u.subs({x:r,y:0}),[ur,ut,uz]):assert sp.simplify(actual-expected)==0
        for value,power in zip(axis,(abs(m+1),abs(m-1),abs(m))):
            reduced=sp.Poly(sp.cancel(value/r**power),r,z)
            assert all(mon[0]%2==0 for mon,_ in reduced.terms())
        assert sp.expand(sp.diff(u[0],x)+sp.diff(u[1],y)+sp.diff(u[2],z))==0
        results.append(dict(mode=m,axis_orders=[abs(m+1),abs(m-1),abs(m)]))
    # Axial window commutator is checked in Cartesian derivatives with a
    # nonconstant polynomial q, so it does not reuse cylindrical code.
    P=(x+I*y)**4*(1+s)*(z+z*z);T=(x+I*y)**4*(2+s)*(1+I*z)
    q=1+2*z+3*z*z+z**3;u=curlpot(P,T);uq=curlpot(q*P,q*T)
    assert all(sp.expand(v)==0 for v in uq-q*u-sp.diff(q,z)*sp.Matrix([sp.diff(P,x),sp.diff(P,y),0]))
    # Ordered Cartesian tensor norm vs complex raising/lowering derivatives.
    f=(x+I*y)**4*(1+s)*(1+z)+I*z*(x-I*y)
    normchecks=[]
    def sq(a):return sp.expand(a*sp.conjugate(a))
    def plus(a):return sp.diff(a,x)+I*sp.diff(a,y)
    def minus(a):return sp.diff(a,x)-I*sp.diff(a,y)
    for j in range(5):
        left=0;right=0
        for a in range(j+1):
            for b in range(j-a+1):
                c=j-a-b
                left+=math.factorial(j)//(math.factorial(a)*math.factorial(b)*math.factorial(c))*sq(sp.diff(f,x,a,y,b,z,c))
        for h in range(j+1):
            for a in range(h+1):
                g=sp.diff(f,z,j-h)
                for _ in range(a):g=plus(g)
                for _ in range(h-a):g=minus(g)
                right+=sp.Rational(math.comb(j,h)*math.comb(h,a),2**h)*sq(g)
        assert sp.expand(left-right)==0
        normchecks.append(j)
    # A wrong radial parity can look small at the axis but fails H5.
    badT=(x**4-6*x*x*y*y+y**4)*sp.sqrt(s)
    bad_derivative=sp.simplify(-sp.diff(badT,x,6).subs({x:1,y:1}))
    assert bad_derivative!=0
    return dict(axis_checks=results,axial_window_commutator=True,ordered_norm_orders=normchecks,
        bad_m4_potential='r^5 cos(4 theta): velocity degree4; a nonzero fifth derivative is homogeneous degree -1, hence not locally L2 across axis.',
        bad_velocity_d5_at_1_1=str(bad_derivative))

def main():
    start=time.time();source=HERE/'compact_potential_reconstruction.py'
    before=hashlib.sha256(source.read_bytes()).hexdigest()
    fit=[manufactured_fit(m,nz) for m,nz in ((0,17),(4,18),(16,17))]
    symbolic=symbolic_checks()
    after=hashlib.sha256(source.read_bytes()).hexdigest()
    out=dict(status='PASS',scope='Manufactured finite algebra and floating-point reconstruction tests; no production fit or residual certificate.',
        reconstruction_sha256_before=before,reconstruction_sha256_after=after,source_changed_during_check=before!=after,
        checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),fits=fit,symbolic=symbolic,seconds=time.time()-start)
    target=HERE/'compact-potential-independent-checks.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',seconds=out['seconds'],output=str(target),source_changed=before!=after)))

if __name__=='__main__':main()
