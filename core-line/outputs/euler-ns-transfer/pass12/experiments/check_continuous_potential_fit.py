#!/usr/bin/env python3
"""Independent small dense-QR / polynomial-quadrature checks."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import hashlib,json
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import lstsq
from numpy.polynomial.legendre import leggauss
import continuous_potential_fit as fit
from compact_potential_reconstruction import FastMAC

HERE=Path(__file__).resolve().parent

def gauss(breaks,n):
    x,w=leggauss(n);s=np.unique(breaks);a=s[:-1,None];b=s[1:,None]
    return ((a+b)/2+(b-a)*x/2).ravel(),((b-a)*w/2).ravel()

def main():
    answers=[];rng=np.random.default_rng(122093)
    for m in (0,4,16):
        g=FastMAC(nr=32,nz=15,R=4,L=16,mapping=3,J=4);g.mapping=3.
        data=.01*(rng.normal(size=(g.nv,g.nz))+1j*rng.normal(size=(g.nv,g.nz)))
        basis,P,T,fitted,backward=fit.fit_mode(g,data,m//4,intervals=10)
        breaks=np.unique(np.r_[np.sqrt(np.unique(basis.knots)),g.rf[1:-1],g.r])
        # Five extra Gauss nodes per span, independent dense target evaluation.
        r,w=gauss(breaks,m+23);wr=w*r
        B,D,L,M=basis.matrices(r);n=basis.n
        def target(col):
            parts=[]
            for nodes,values,sign in ((g.rf[1:-1],data[g.slr,col],-1),
                    (g.r,data[g.slt,col],-1),(g.r,data[g.slz,col],1)):
                xx=np.r_[-g.R,-nodes[::-1],nodes,g.R]
                vv=np.r_[0.,sign*values[::-1],values,0.]
                parts.append(CubicSpline(xx,vv)(r))
            return np.concatenate(parts)
        errors=[];orthogonality=[];ratios=[];crosses=[]
        for col in (0,1,g.nz-1):
            k=g.kz[col]
            A=np.block([[1j*k*D,1j*M],[-k*M,-D],[-L,np.zeros_like(L)]])
            W=np.tile(wr,3);Aw=np.sqrt(W)[:,None]*A;bw=np.sqrt(W)*target(col)
            qr,_,rank,_=lstsq(Aw,bw,lapack_driver='gelsy');assert rank==2*n
            coeff=np.r_[P[:,col],T[:,col]];v=Aw@coeff;res=v-bw
            errors.append(float(np.linalg.norm(Aw@(coeff-qr))/max(np.linalg.norm(bw),1e-300)))
            orthogonality.append(float(np.linalg.norm(Aw.conj().T@res)/max(np.linalg.norm(Aw)*np.linalg.norm(bw),1e-300)))
            ratios.append(float(np.linalg.norm(v)/np.linalg.norm(bw)))
            G=Aw.conj().T@Aw
            crosses.append(float(np.linalg.norm(G[:n,n:])/np.linalg.norm(G)))
        out=dict(m=m,mapping=g.mapping,dense_QR_relative=max(errors),normal_orthogonality=max(orthogonality),
            maximum_continuous_projected_norm_ratio=max(ratios),cross_Gram_relative=max(crosses),
            reported_backward=backward)
        assert max(errors)<1e-8,out
        assert max(orthogonality)<1e-10,out
        assert max(ratios)<=1+1e-10,out
        assert max(crosses)<1e-13,out
        answers.append(out)
    # Exact counterexample to sample norm controlling either core rotation
    # or continuous energy: a toroidal polynomial bump below the first cell.
    # T=(eps²/20)*(1-r²/eps²)^10_+ q(z), so utheta=r(1-r²/eps²)^9_+ q(z).
    eps=g.r[0]/2
    sampled=g.r*np.maximum(1-g.r*g.r/eps**2,0)**9
    assert np.array_equal(sampled,np.zeros_like(sampled))
    # At fixed z with q=1, integral_R2 |u|² = pi eps^4 Beta(2,19)
    # = pi eps^4/(19*20), while central rotation is exactly one.
    hidden=dict(epsilon=float(eps),sampled_velocity_norm=0.,central_rotation=1.,
        exact_horizontal_L2_squared_factor='pi*epsilon^4/380',
        interpretation='Multiplying by any amplitude keeps all MAC velocity samples zero and scales actual energy by amplitude²; H4 and higher norms are also uncontrolled.',
        scope='An admissible smooth-enough compact potential family, not a claim this exact bump lies in the fixed twelve-interval degree-nine trial space.')
    result=dict(status='PASS',scope='Small dense-QR equivalence, exact polynomial Gram degree checks and a sample-norm counterexample; no production fit or residual enclosure.',
        source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in
            ('continuous_potential_fit.py','compact_potential_reconstruction.py','observe_compact_reconstruction.py','check_continuous_potential_fit.py')},
        results=answers,invisible_axis_bump=hidden)
    path=HERE/'continuous-potential-independent-checks.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',output=str(path),largest_dense_QR_error=max(a['dense_QR_relative'] for a in answers))))

if __name__=='__main__':main()
