"""Independent matrix specification/check for a harmonic-exterior MAC model.

This is a small algebra reference, not a time integrator or a modification
of evolve_full_mac.py. Its finite-difference Gram is a discrete model.
"""
from pathlib import Path
import json
import numpy as np
from check_exterior_harmonic_mass import dtn


def assemble(m,k,n=8,R=3.):
    rf=R*np.linspace(0,1,n+1)**1.25
    r=(rf[:-1]+rf[1:])/2
    h=np.diff(rf);dc=np.diff(r);delta=R-r[-1]
    vol=np.diff(rf*rf)/2
    zero=(m==0 and k==0)
    nq=n-1 if zero else n
    nv=nq+2*n
    Q=np.zeros((n,nv),complex);Q[:nq,:nq]=np.eye(nq)
    T=np.zeros((n,nv),complex);T[:,nq:nq+n]=np.eye(n)
    Z=np.zeros((n,nv),complex);Z[:,nq+n:]=np.eye(n)
    qR=Q[-1]
    kap=0. if zero else dtn(m,k,R)
    lt=0. if zero else -1j*m/(R*kap)
    lz=0. if zero else -1j*k/kap
    dext=0. if zero else 1+2*R*(m*m/R**2+k*k)/kap+m*m/(R*kap)**2
    w=np.r_[rf[1:-1]*dc,R*delta][:nq]
    interior_mass=np.r_[w,vol,vol]
    M=np.diag(interior_mass).astype(complex)
    if not zero:M+=R/kap*np.outer(qR.conj(),qR)
    prev=np.vstack((np.zeros((1,nv)),Q[:-1]))
    Iq=(Q+prev)/2
    Br=(Q-prev)/h[:,None]
    Dt=(T[1:]-T[:-1])/dc[:,None]
    Dz=(Z[1:]-Z[:-1])/dc[:,None]
    Bt=np.vstack((Dt,(lt*qR-T[-1])/delta))
    Bz=np.vstack((Dz,(lz*qR-Z[-1])/delta))
    wg=np.r_[rf[1:-1]*dc,R*delta]
    angular=(1j*m*Iq-T,Iq+1j*m*T,1j*m*Z)
    def gram(A,weights):return A.conj().T@(weights[:,None]*A)
    K=gram(Br,vol)+gram(Bt,wg)+gram(Bz,wg)
    K+=sum((gram(A,vol/r**2) for A in angular),np.zeros_like(K))
    K+=k*k*np.diag(interior_mass)+dext*np.outer(qR.conj(),qR)
    D=(rf[1:,None]*Q-rf[:-1,None]*prev)/vol[:,None]+1j*m*T/r[:,None]+1j*k*Z
    G=-np.linalg.solve(M,D.conj().T*vol[None,:])
    gather=np.vstack((Iq,T,Z))
    return dict(rf=rf,r=r,vol=vol,M=M,K=K,D=D,G=G,Q=Q,T=T,Z=Z,
                qR=qR,kap=kap,lt=lt,lz=lz,dext=dext,delta=delta,
                Br=Br,Bt=Bt,Bz=Bz,angular=angular,wg=wg,
                interior_mass=interior_mass,gather=gather,nq=nq,nv=nv)


def projector(M,D,vol):
    sq=np.sqrt(vol)
    coupling=D.conj().T*sq[None,:]
    Minv=np.linalg.solve(M,coupling)
    A=sq[:,None]*(D@Minv)
    eig,U=np.linalg.eigh((A+A.conj().T)/2)
    keep=eig>max(1.,eig[-1])*1.e-12
    inv=(U[:,keep]/eig[keep])@U[:,keep].conj().T
    return np.eye(len(M))-Minv@inv@(sq[:,None]*D)


def check_case(m,k):
    b=assemble(m,k);M,K,D,G=b['M'],b['K'],b['D'],b['G']
    vol=b['vol'];P=projector(M,D,vol)
    rng=np.random.default_rng(12039)
    u=rng.normal(size=len(M))+1j*rng.normal(size=len(M))
    v=P@u
    require=[np.max(abs(D@P)),np.max(abs(P@P-P)),
             np.max(abs(P.conj().T@M-M@P)),np.max(abs(K-K.conj().T))]
    assert max(require)<2.e-11,(m,k,require)
    assert np.linalg.eigvalsh(K)[0]>-1.e-11
    direct=np.sum(vol*abs(b['Br']@u)**2)
    direct+=np.sum(b['wg']*(abs(b['Bt']@u)**2+abs(b['Bz']@u)**2))
    direct+=sum(np.sum(vol/b['r']**2*abs(A@u)**2) for A in b['angular'])
    direct+=k*k*np.sum(b['interior_mass']*abs(u)**2)+b['dext']*abs(b['qR']@u)**2
    gram_error=abs(direct-np.vdot(u,K@u).real)/max(1.,direct)
    assert gram_error<1.e-14
    # The weak scatter is adjoint to gather in the FULL modal mass.
    F=rng.normal(size=3*len(vol))+1j*rng.normal(size=3*len(vol))
    scatter=np.linalg.solve(M,b['gather'].conj().T@(np.tile(vol,3)*F))
    weak_defect=abs(np.vdot(v,M@(P@scatter))-np.vdot(b['gather']@v,np.tile(vol,3)*F))
    assert weak_defect<1.e-10
    if b['kap']:
        p=rng.normal(size=len(vol))+1j*rng.normal(size=len(vol))
        q=b['qR']@(G@p)
        pR=p[-1]/(1+b['kap']*b['delta'])
        assert abs(q+b['kap']*pR)<1.e-12
        assert abs(b['lt']*q-1j*m*pR/b['rf'][-1])<1.e-12
        assert abs(b['lz']*q-1j*k*pR)<1.e-12
        # The new outer derivative rows realize those same traces.
        assert abs((b['T']@u)[-1]+b['delta']*(b['Bt']@u)[-1]-b['lt']*(b['qR']@u))<1.e-12
        assert abs((b['Z']@u)[-1]+b['delta']*(b['Bz']@u)[-1]-b['lz']*(b['qR']@u))<1.e-12
    else:
        assert np.max(abs(v[:b['nq']]))<1.e-12
        assert np.linalg.norm(G@np.ones(len(vol)))<1.e-12
    # CN on the exactly constrained finite-dimensional space dissipates.
    alpha=.01
    sq=np.sqrt(vol);C=D.conj().T*sq[None,:]
    saddle=np.block([[M+alpha*K,C],[C.conj().T,np.zeros((len(vol),len(vol)))]])
    rhs=np.r_[(M-alpha*K)@v,np.zeros(len(vol))]
    sol=np.linalg.lstsq(saddle,rhs,rcond=1.e-13)[0][:len(M)]
    E0=np.vdot(v,M@v).real;E1=np.vdot(sol,M@sol).real
    cn=E1-E0+alpha*np.vdot(sol+v,K@(sol+v)).real
    assert abs(cn)<1.e-9*max(1.,E0)
    assert E1<=E0+1.e-10 and np.max(abs(D@sol))<1.e-10
    return dict(m=m,k=k,operator_defects=list(map(float,require)),
                positive_gram_relative_error=float(gram_error),
                nonlinear_weak_adjoint_defect=float(weak_defect),
                CN_energy_identity_relative_defect=float(abs(cn)/max(1.,E0)))


def axis_injection_check():
    """Example strict mode-four axis patch with two even Taylor powers.

    First three q-face, theta-center and z-center values come from six
    regular helical coefficients. Other values remain independent.
    This checks congruence reduction, not convergence of this tiny patch.
    """
    m=4;b=assemble(m,.7);nv=b['nv'];nq=b['nq'];n=len(b['r'])
    chosen=list(range(3))+list(range(nq,nq+3))+list(range(nq+n,nq+n+3))
    free=[i for i in range(nv) if i not in chosen]
    C=np.zeros((nv,6+len(free)),complex)
    scale=b['rf'][3]
    for j,r in enumerate(b['rf'][1:4]/scale):
        C[j,:4]=[r**5/2,r**7/2,r**3/2,r**5/2]
    for j,r in enumerate(b['r'][:3]/scale):
        C[nq+j,:4]=[r**5/(2j),r**7/(2j),-r**3/(2j),-r**5/(2j)]
        C[nq+n+j,4:6]=[r**4,r**6]
    C[free,6+np.arange(len(free))]=1
    M=C.conj().T@b['M']@C;K=C.conj().T@b['K']@C;D=b['D']@C
    P=projector(M,D,b['vol'])
    div=np.max(abs(D@P));orth=np.max(abs(P.conj().T@M-M@P))
    assert div<2.e-10 and orth<2.e-10
    assert np.linalg.eigvalsh(M)[0]>0 and np.linalg.eigvalsh(K)[0]>0
    return {'mode':m,'reduced_dimension':len(M),'full_dimension':nv,
            'divergence':float(div),'kinetic_self_adjoint_defect':float(orth)}


if __name__=='__main__':
    cases=[check_case(*p) for p in [(0,0),(0,.7),(4,0),(4,.7),(4,-.7),(8,2.3)]]
    plus=assemble(4,.7);minus=assemble(4,-.7)
    S=np.ones(plus['nv']);S[-len(plus['r']):]=-1
    reflection=np.max(abs(minus['K']-S[:,None]*plus['K']*S[None,:]))
    assert reflection<1.e-13
    out={'status':'PASS','scope':'Independent finite-dimensional harmonic-exterior MAC model; no NS evolution or residual enclosure',
         'cases':cases,'axial_reflection_stiffness_defect':float(reflection),
         'axis_constraint_congruence':axis_injection_check()}
    Path(__file__).with_name('exterior-mac-block-check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
