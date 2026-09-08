#!/usr/bin/env python3
"""Independent small-grid structural tests of evolve_full_mac.MAC.

Floating-point identities and manufactured consistency, not NS validation.
No production solver file is modified.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft
from scipy.signal import resample
from evolve_full_mac import MAC


def norm(g,U):
    return np.sqrt(max(0.,g.inner(U,U)))


def relative(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a),np.linalg.norm(b),1e-300))


def random_state(g,rng):
    U=rng.normal(size=(g.J+1,g.nv,g.nz))+1j*rng.normal(size=(g.J+1,g.nv,g.nz))
    U[0]=U[0].real
    return U


def explicit_gram(g,a,j):
    ar,at,az=a[g.slr],a[g.slt],a[g.slz]
    uc=g.I@ar;m=g.modes[j]
    radial=np.sum(g.V*abs(g.Df@ar)**2)+np.sum(g.Wg*(abs(g.E@at)**2+abs(g.E@az)**2))
    angular=np.sum(g.V/g.r**2*(abs(1j*m*uc-at)**2+abs(1j*m*at+uc)**2+m*m*abs(az)**2))
    return float(radial+angular)


def direct_stokes(g,ahat,j,col,h):
    # Independent dense saddle solve, including the actual SIGNED axial k.
    m=g.modes[j];k=g.kz[col];n=g.nr;alpha=h*g.nu/2
    D=np.concatenate((g.Divr.toarray(),np.diag(1j*m/g.r),1j*k*np.eye(n)),axis=1)
    BM=-D.conj().T*g.V[None,:]
    A=np.diag(g.mass*(1+alpha*k*k))+alpha*g.K[j].toarray()
    saddle=np.block([[A,BM],[BM.conj().T,np.zeros((n,n))]])
    rhs=np.zeros(g.nv+n,complex)
    rhs[:g.nv]=g.mass*ahat-alpha*(g.K[j]@ahat+g.mass*k*k*ahat)
    take=np.arange(g.nv+n)
    if j==0 and g.ki[col]==0:take=take[:-1]
    out=np.zeros_like(rhs);out[take]=np.linalg.solve(saddle[np.ix_(take,take)],rhs[take])
    return out[:g.nv]


def one_grid(nr,nz,dealias='padding'):
    rng=np.random.default_rng(903+nr)
    g=MAC(nr=nr,nz=nz,R=4.,L=8.,mapping=2.,J=2,nu=.001,dealias=dealias)
    raw=random_state(g,rng);U=g.project(raw);V=g.project(random_state(g,rng))
    P2=g.project(U);complement=raw-U
    answer=dict(nr=nr,nz=nz,J=g.J,dealias=g.dealias,nzpad=g.nzpad)
    adjoint=[];gather=[];grams=[];lowest=[]
    for j in range(g.J+1):
        a=fft(raw[j],axis=-1)
        p=rng.normal(size=(nr,nz))+1j*rng.normal(size=(nr,nz))
        da=g.div_hat(a,j);gp=g.grad_hat(p,j)
        l=np.sum(g.V[:,None]*p.conj()*da)
        rr=np.sum(g.mass[:,None]*gp.conj()*a)
        adjoint.append(float(abs(l+rr)/max(abs(l),abs(rr),1e-300)))
        f=rng.normal(size=(3,nr,nz))+1j*rng.normal(size=(3,nr,nz))
        l=np.sum(g.mass[:,None]*raw[j].conj()*g.scatter(f))
        rr=np.sum(g.V[None,:,None]*g.gather(raw[j]).conj()*f)
        gather.append(float(abs(l-rr)/max(abs(l),abs(rr),1e-300)))
        vec=rng.normal(size=g.nv)+1j*rng.normal(size=g.nv)
        form=np.vdot(vec,g.K[j]@vec)
        expected=explicit_gram(g,vec,j)
        grams.append(float(abs(form-expected)/expected))
        scaled=g.K[j].toarray()/np.sqrt(g.mass[:,None]*g.mass[None,:])
        lowest.append(float(np.linalg.eigvalsh(scaled)[0]))
    answer.update(gradient_negative_adjoint_relative=max(adjoint),
        gather_scatter_adjoint_relative=max(gather),
        independent_gradient_gram_relative=max(grams),
        minimum_generalized_viscous_eigenvalue=min(lowest),
        projection_idempotence_relative=norm(g,P2-U)/norm(g,U),
        projected_divergence_L2=g.divergence_norm(U),
        projected_divergence_over_input_divergence=g.divergence_norm(U)/g.divergence_norm(raw),
        projection_orthogonality_relative=abs(g.inner(U,complement))/(norm(g,U)*norm(g,complement)),
        projected_norm_ratio=norm(g,U)/norm(g,raw))
    # Weighted self-adjointness, including filtering and real mode-zero space.
    W=random_state(g,rng)
    lhs=g.inner(g.project(raw),W);rhs=g.inner(raw,g.project(W))
    answer['projection_self_adjoint_relative']=abs(lhs-rhs)/max(norm(g,raw)*norm(g,W),1e-300)

    F=g.nonlinear(U,project=False);PF=g.project(F)
    answer['raw_lamb_energy_work_relative']=abs(g.inner(U,F))/(norm(g,U)*norm(g,F))
    answer['projected_lamb_energy_work_relative']=abs(g.inner(U,PF))/(norm(g,U)*norm(g,PF))
    # Arbitrary angular rotations, not merely a collocation-node shift.
    phase=np.exp(1j*g.modes*.173)[:,None,None]
    answer['nonlinear_rotation_equivariance_relative']=norm(g,g.nonlinear(phase*U)-phase*PF)/norm(g,PF)
    shift=3
    answer['nonlinear_axial_translation_relative']=norm(g,g.nonlinear(np.roll(U,shift,axis=-1))-np.roll(PF,shift,axis=-1))/norm(g,PF)

    h=.007
    after=g.diffuse(U,h);mid=(after+U)/2
    E0=.5*g.inner(U,U);E1=.5*g.inner(after,after)
    loss=h*g.nu*g.dissipation(mid)
    answer['CN_energy_ratio']=E1/E0
    answer['CN_midpoint_energy_identity_relative']=abs(E1-E0+loss)/(E0+E1+loss)
    answer['CN_divergence_L2']=g.divergence_norm(after)
    direct_errors=[]
    if nr==16:
        ah=fft(U,axis=-1);bh=fft(after,axis=-1)
        indices=[(0,0),(0,3),(0,-3),(1,3),(1,-3),(2,1),(2,-1)]
        if dealias=='padding':indices.extend([(1,nz//2-1),(1,-nz//2+1)])
        for j,kindex in indices:
            col=int(np.flatnonzero(g.ki==kindex)[0])
            direct_errors.append(relative(bh[j,:,col],direct_stokes(g,ah[j,:,col],j,col,h)))
    answer['direct_signed_k_stokes_relative']=max(direct_errors,default=None)

    # Cartesian oddness in this C4 class: horizontal components even in z,
    # vertical component odd in z. The pressure and Stokes maps must preserve it.
    rev=(-np.arange(nz))%nz
    odd=raw.copy();ref=odd[...,rev]
    odd[:,:g.slz.start]=(odd[:,:g.slz.start]+ref[:,:g.slz.start])/2
    odd[:,g.slz]=(odd[:,g.slz]-ref[:,g.slz])/2
    odd=g.project(odd)
    def odd_defect(a):
        transformed=a[...,rev].copy();transformed[:,g.slz]*=-1
        return norm(g,a-transformed)/max(norm(g,a),1e-300)
    answer['projection_odd_symmetry_relative']=odd_defect(odd)
    answer['nonlinear_odd_symmetry_relative']=odd_defect(g.nonlinear(odd))
    answer['CN_odd_symmetry_relative']=odd_defect(g.diffuse(odd,h))
    # A single small split step tests the implementation path, not stability
    # over a proposed research interval.
    small=.01*U
    stepped=g.step(small,1e-4)
    answer['one_split_step_divergence_L2']=g.divergence_norm(stepped)
    answer['one_split_step_energy_ratio']=g.inner(stepped,stepped)/g.inner(small,small)
    answer['mode_zero_imaginary_max']=float(np.max(abs(stepped[0].imag)))
    return answer


def manufactured_curl(nr):
    g=MAC(nr=nr,nz=64,R=4.,L=8.,mapping=2.,J=1,dealias='padding')
    r=g.r[:,None];z=g.z[None,:];R=g.R;k=2*np.pi/g.L
    q=1+.2*np.cos(k*z);qz=-.2*k*np.sin(k*z)
    ss=r*r;f=(1-ss/R**2)**4;fs=-4/R**2*(1-ss/R**2)**3
    fss=12/R**4*(1-ss/R**2)**2
    result={}
    for j in (0,1):
        a=np.zeros((g.nv,g.nz),complex)
        if j==0:
            a[g.slt]=r*f*q
            exact=np.stack((-r*f*qz,np.zeros_like(r*q),(2*f+2*ss*fs)*q))
        else:
            rf=g.rf[1:-1,None];sf=rf*rf;ff=(1-sf/R**2)**4
            a[g.slr]=4j*rf**3/R**4*ff*q
            a[g.slt]=-r**3/R**4*(4*f+2*ss*fs)*q
            exact=np.stack((r**3/R**4*(4*f+2*ss*fs)*qz,
                4j*r**3/R**4*f*qz,-4*r**4/R**4*(5*fs+ss*fss)*q))
        computed=g.curl_cells(a,j)
        weight=g.V[None,:,None]*g.L/g.nz
        error=float(np.sqrt(np.sum(weight*abs(computed-exact)**2)/np.sum(weight*abs(exact)**2)))
        result[str(g.modes[j])]=dict(curl_relative_L2=error,
            first_four_rows_absolute_error=float(np.max(abs(computed[:,:4]-exact[:,:4]))),
            max_absolute_error=float(np.max(abs(computed-exact))))
    rigid=np.zeros((g.nv,g.nz),complex);rigid[g.slt]=g.r[:,None]
    result['rigid_core_first_four_cells_vorticity_error']=float(np.max(abs(g.curl_cells(rigid,0)[2,:4]-2)))
    return dict(nr=nr,dealias=g.dealias,results=result,
        scope='Smooth wall-compatible manufactured m0 and m4 toroidal fields; sampled before projection. Rigid-core test excludes the incompatible outer wall.')


def lift(U,g,G):
    src=fft(U,axis=-1,norm='forward')
    dst=np.zeros((G.J+1,G.nv,G.nz),complex)
    for col,k in enumerate(g.ki):
        if g.keep[col]:dst[:g.J+1,:,int(k)%G.nz]=src[:,:,col]
    return ifft(dst,axis=-1,norm='forward')


def dealias_check():
    g=MAC(nr=12,nz=32,R=4.,L=8.,mapping=2.,J=2,dealias='padding')
    G=MAC(nr=12,nz=96,R=4.,L=8.,mapping=2.,J=4,dealias='padding')
    rng=np.random.default_rng(815);U=g.project(random_state(g,rng))
    F=g.nonlinear(U);FF=G.nonlinear(lift(U,g,G))
    big=fft(FF,axis=-1,norm='forward');small=np.zeros_like(U)
    for col,k in enumerate(g.ki):
        if g.keep[col]:small[:,:,col]=big[:g.J+1,:,int(k)%G.nz]
    restricted=ifft(small,axis=-1,norm='forward')
    # Explicitly verify that this test retains substantial input beyond the
    # old two-thirds state filter, and that up/down sampling are adjoints on
    # the retained state space in the physical L/nz norm.
    ah=fft(U,axis=-1,norm='forward');high=ah.copy()
    high[:,:,abs(g.ki)<g.nz/3]=0
    high_state=ifft(high,axis=-1,norm='forward')
    up=resample(U,g.nzpad,axis=-1)
    force=rng.normal(size=up.shape)+1j*rng.normal(size=up.shape);force[0]=force[0].real
    factors=np.r_[1.,np.full(g.J,2.)]
    def arbitrary_grid_inner(a,b):
        return float(2*np.pi*g.L/a.shape[-1]*np.real(np.sum(factors[:,None,None]*g.mass[None,:,None]*a.conj()*b)))
    lhs=arbitrary_grid_inner(up,force);rhs=arbitrary_grid_inner(U,resample(force,g.nz,axis=-1))
    adjoint_error=abs(lhs-rhs)/max(np.sqrt(arbitrary_grid_inner(up,up)*arbitrary_grid_inner(force,force)),1e-300)
    return dict(dealias='padding',angular_J_low=2,angular_J_reference=4,nz_low=32,nz_reference=96,
        nz_product_grid=g.nzpad,input_energy_fraction_beyond_old_filter=g.inner(high_state,high_state)/g.inner(U,U),
        largest_retained_axial_integer=int(max(abs(g.ki[g.keep]))),
        resample_adjoint_relative=adjoint_error,
        retained_nonlinearity_relative=norm(g,F-restricted)/norm(g,F),
        scope='Same radial operators and input coefficients, independent larger angular/axial product grid.')


def run():
    start=time.monotonic()
    structural=[one_grid(n,n*2,'padding') for n in (16,32,64)]
    structural.append(one_grid(16,32,'filter'))
    curls=[manufactured_curl(n) for n in (16,32,64)]
    alias=dealias_check()
    for a in structural:
        for key in ['gradient_negative_adjoint_relative','gather_scatter_adjoint_relative',
            'independent_gradient_gram_relative','projection_idempotence_relative',
            'projection_orthogonality_relative','projection_self_adjoint_relative',
            'raw_lamb_energy_work_relative','projected_lamb_energy_work_relative',
            'nonlinear_rotation_equivariance_relative','nonlinear_axial_translation_relative',
            'CN_midpoint_energy_identity_relative','projection_odd_symmetry_relative',
            'nonlinear_odd_symmetry_relative','CN_odd_symmetry_relative']:
            assert a[key]<2e-10,(a['nr'],key,a[key])
        assert a['minimum_generalized_viscous_eigenvalue']>=-1e-9
        assert a['projected_norm_ratio']<=1+1e-12
        assert a['CN_energy_ratio']<=1+1e-12
        assert a['one_split_step_energy_ratio']<=1+1e-10
        assert a['projected_divergence_over_input_divergence']<1e-10
        if a['direct_signed_k_stokes_relative'] is not None:
            assert a['direct_signed_k_stokes_relative']<2e-10
    for m in ('0','4'):
        assert curls[-1]['results'][m]['curl_relative_L2']<curls[-2]['results'][m]['curl_relative_L2']
    assert max(v['results']['rigid_core_first_four_cells_vorticity_error'] for v in curls)<1e-12
    assert alias['retained_nonlinearity_relative']<2e-10
    assert alias['resample_adjoint_relative']<2e-10
    assert alias['input_energy_fraction_beyond_old_filter']>.1
    return dict(status='PASS bounded floating-point structural/manufactured MAC checks',
        solver_sha256=hashlib.sha256(Path(__file__).with_name('evolve_full_mac.py').read_bytes()).hexdigest(),
        structural=structural,manufactured_curl=curls,dealiasing=alias,seconds=time.monotonic()-start,
        scope='Discrete adjoints, projection, Gram viscosity, constrained CN, Lamb energy work, symmetries, and radial curl consistency only. No whole-space reconstruction, interval residual, axis high-derivative enclosure, or research endpoint is certified.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(Path(__file__).with_name('full-mac-structural-checks.json')))
    args=ap.parse_args();result=run();Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
