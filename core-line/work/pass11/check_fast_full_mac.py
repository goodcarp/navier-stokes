#!/usr/bin/env python3
"""Independent equivalence/structure tests of the FastMAC implementation.

These are floating-point discrete-operator tests, not a continuum NS certificate.
The reference solver, fast solver and this checker are hashed before and after.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import hashlib,json,time
from pathlib import Path
import numpy as np
from scipy.fft import fft,ifft
from evolve_full_mac import MAC
from fast_full_mac import FastMAC
from check_full_mac import random_state,norm,relative,direct_stokes

HERE=Path(__file__).resolve().parent
FILES=('evolve_full_mac.py','fast_full_mac.py','check_fast_full_mac.py','check_full_mac.py')
def hashes():return {s:hashlib.sha256((HERE/s).read_bytes()).hexdigest() for s in FILES}
def rel(g,x,y):return norm(g,x-y)/max(norm(g,x),norm(g,y),1e-300)
def energy_work(g,x,f):return abs(g.inner(x,f))/max(norm(g,x)*norm(g,f),1e-300)
def restrict(g,U):
    a=fft(U,axis=-1);a[:,:,~g.keep]=0
    out=ifft(a,axis=-1);out[0]=out[0].real
    return out

def exact_template(g,j,kindex,h):
    # Independently assemble the full signed-frequency saddle matrix in the
    # ORIGINAL unknown order. Compare after the optimized permutation.
    n=g.nr;alpha=h*g.nu/2;k=2*np.pi*kindex/g.L
    D=np.concatenate((g.Divr.toarray(),np.diag(1j*g.modes[j]/g.r),1j*k*np.eye(n)),axis=1)
    BM=-D.conj().T*g.V[None,:]
    A=np.diag(g.mass*(1+alpha*k*k))+alpha*g.K[j].toarray()
    return np.block([[A,BM],[BM.conj().T,np.zeros((n,n))]])

def one(nr,nz,J,dealias):
    kw=dict(nr=nr,nz=nz,J=J,dealias=dealias,R=4.,L=8.,mapping=2.,nu=.001)
    r=MAC(**kw);g=FastMAC(**kw);rng=np.random.default_rng(20260908+nr*13+nz+J)
    raw=random_state(g,rng);U=.01*r.project(raw)
    projection=rel(g,r.project(raw),g.project(raw))
    f=r.nonlinear(U);ff=g.nonlinear(U)
    fr=r.nonlinear(U,project=False);frf=g.nonlinear(U,project=False)
    h=.003
    vr=r.diffuse(U,h);vf=g.diffuse(U,h)
    mid=(U+vf)/2;e0=.5*g.inner(U,U);e1=.5*g.inner(vf,vf);loss=h*g.nu*g.dissipation(mid)
    step_r=r.step(U,1e-4);step_f=g.step(U,1e-4)
    answer=dict(parameters=kw,nz_product_grid=g.nzpad if dealias=='padding' else nz,
        nphi_product_grid=g.nphi,projection_equivalence=projection,
        nonlinear_projected_equivalence=rel(g,f,ff),
        nonlinear_raw_retained_equivalence=rel(g,restrict(g,fr),frf),
        nonlinear_raw_unfiltered_difference=rel(g,fr,frf),
        raw_fast_energy_work=energy_work(g,U,frf),projected_fast_energy_work=energy_work(g,U,ff),
        nonlinear_divergence_over_field_norm=g.divergence_norm(ff)/max(norm(g,ff),1e-300),
        CN_equivalence=rel(g,vr,vf),CN_energy_ratio=e1/e0,
        CN_midpoint_energy_identity=abs(e1-e0+loss)/(e0+e1+loss),
        CN_divergence_over_field_norm=g.divergence_norm(vf)/norm(g,vf),
        full_step_equivalence=rel(g,step_r,step_f),
        full_step_energy_ratio=g.inner(step_f,step_f)/g.inner(U,U),
        full_step_divergence_over_field_norm=g.divergence_norm(step_f)/norm(g,step_f),
        mode_zero_imaginary_max=float(abs(step_f[0].imag).max()))
    # Both signs of the HIGHEST retained axial frequency; nonzero m modes
    # have independent complex +k and -k coefficients, not false Hermitian pairs.
    high=int(max(abs(g.ki[g.keep])))
    ah=fft(U,axis=-1);bh=fft(vf,axis=-1)
    checks=[]
    for j in sorted(set((0,J))):
        for kval in sorted(set((0,1,-1,high,-high))):
            col=int(np.flatnonzero(g.ki==kval)[0])
            dense=direct_stokes(g,ah[j,:,col],j,col,h)
            checks.append(dict(m=int(g.modes[j]),k=kval,dense_saddle_relative=relative(bh[j,:,col],dense)))
    answer['signed_frequency_dense_CN']=checks
    # Template entries and index positions, including alpha=0 and zero mode
    # gauge removal. For nonzero k use magnitude for the optimized template.
    template=[]
    for hh in (0.,h):
        alpha=hh*g.nu/2
        for j in sorted(set((0,J))):
            for kval in (0,1,high):
                mat,perm,diag,linear=g._stokes_template(j,alpha,j==0 and kval==0)
                a=mat.copy();k=2*np.pi*kval/g.L
                a.data[diag]+=alpha*k*k*g.mass;a.data[linear]*=k
                full=exact_template(g,j,kval,hh)
                err=relative(a.toarray(),full[np.ix_(perm,perm)])
                assert len(set(diag))==g.nv
                assert len(set(linear))==len(linear)
                assert not set(diag)&set(linear)
                template.append(err)
    answer['maximum_independent_template_matrix_relative']=max(template)
    answer['CN_h_zero_is_projection']=rel(g,g.diffuse(U,0.),U)
    # Arbitrary angular phase and exact fractional axial shift catch failure
    # to reverse BOTH m and k when constructing the physical reality partner.
    phase=np.exp(1j*g.modes*.137)[:,None,None]
    answer['angular_equivariance']=rel(g,g.nonlinear(phase*U),phase*ff)
    axialphase=np.exp(1j*g.kz*.217)[None,None,:]
    shifted=ifft(fft(U,axis=-1)*axialphase,axis=-1);shifted[0]=shifted[0].real
    expected=ifft(fft(ff,axis=-1)*axialphase,axis=-1);expected[0]=expected[0].real
    answer['axial_equivariance']=rel(g,g.nonlinear(shifted),expected)
    # Reflection around the actual z=0 grid point, valid on odd and even grids.
    rev=(2*(nz//2)-np.arange(nz))%nz
    O=raw.copy();ref=raw[...,rev]
    O[:,:g.slz.start]=(O[:,:g.slz.start]+ref[:,:g.slz.start])/2
    O[:,g.slz]=(O[:,g.slz]-ref[:,g.slz])/2;O=g.project(.01*O)
    def odd_error(X):
        Y=X[...,rev].copy();Y[:,g.slz]*=-1
        return rel(g,X,Y)
    answer['odd_C4_full_step_defect']=odd_error(g.step(O,1e-4))
    for name in ('projection_equivalence','nonlinear_projected_equivalence',
        'nonlinear_raw_retained_equivalence','raw_fast_energy_work','projected_fast_energy_work',
        'CN_equivalence','CN_midpoint_energy_identity','full_step_equivalence',
        'maximum_independent_template_matrix_relative','CN_h_zero_is_projection',
        'angular_equivariance','axial_equivariance','odd_C4_full_step_defect'):
        assert answer[name]<2e-10,(kw,name,answer[name])
    for name in ('nonlinear_divergence_over_field_norm','CN_divergence_over_field_norm','full_step_divergence_over_field_norm'):
        assert answer[name]<2e-10,(kw,name,answer[name])
    assert max(v['dense_saddle_relative'] for v in checks)<2e-10
    assert answer['CN_energy_ratio']<=1+1e-12
    assert answer['full_step_energy_ratio']<=1+1e-10
    return answer

def run():
    start=time.time();before=hashes()
    # nr=28 crosses the fused FFT's 24-row block boundary.
    specs=[(12,31,0,'padding'),(12,32,1,'padding'),(16,33,2,'padding'),
        (28,64,3,'padding'),(12,31,2,'filter'),(16,32,2,'filter')]
    results=[one(*spec) for spec in specs]
    after=hashes();assert before==after,'Solver/checker mutated during test'
    return dict(status='PASS',scope='Floating-point equivalence and discrete structure only; no full-space NS validation.',
        tested_domain='Real m=0, C4 Fourier representation, projected retained axial state; even Nyquist removed.',
        raw_api_difference='Fast raw nonlinearity already filters excluded axial output; reference raw branch does not.',
        source_sha256=before,results=results,seconds=time.time()-start)

if __name__=='__main__':
    result=run();path=HERE/'fast-full-mac-checks.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],seconds=result['seconds'],output=str(path),
        maximum_projected_nonlinearity_error=max(a['nonlinear_projected_equivalence'] for a in result['results']),
        maximum_full_step_error=max(a['full_step_equivalence'] for a in result['results']))))
