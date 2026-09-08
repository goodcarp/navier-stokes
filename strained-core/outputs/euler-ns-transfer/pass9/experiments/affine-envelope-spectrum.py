"""Exploratory quadrature of an exact linear affine envelope model.

No finite-amplitude NS or circular-vortex time integration is performed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,argparse
from pathlib import Path
import numpy as np
from scipy.special import gammainc,gamma
here=Path(__file__).resolve().parent
for dependency_dir in (here.parent/'pass4', here.parents[1]/'pass4/experiments'):
    if (dependency_dir/'evaluate_pressure_gate.py').is_file():
        sys.path.insert(0,str(dependency_dir))
        break
else:
    raise FileNotFoundError('Expected sibling pass4/evaluate_pressure_gate.py or pass4/experiments/evaluate_pressure_gate.py')
from evaluate_pressure_gate import cutoff

def run(n=8192,length=2.,A=960.,end=.02,steps=2000):
    rs=.21547557533348247
    B=float(cutoff(np.array([(rs/.315)**2]))[0])
    h=4/rs; Om=A*B; shear=2*Om; S=shear*h
    c=1.4; nu=.001; k0=-20.
    dx=length/n
    x=(np.arange(n)-n//2)*dx
    e=cutoff(((rs+x-.14)/.1)**2)
    psi=e*np.exp(1j*k0*x)
    xi=2*np.pi*np.fft.fftfreq(n,d=dx)
    dk=2*np.pi/length
    psihat=dx/np.sqrt(2*np.pi)*np.fft.fft(np.fft.ifftshift(psi))
    zeta=(xi*xi+h*h)*psihat
    z2=abs(zeta)**2
    Z0=float(np.sum(z2)*dk/4)
    times=np.linspace(0,end,steps+1)
    data=[]
    for t in times:
        if t==0: J0=J1=J2=0.
        else:
            J0=-np.expm1(-2*c*t)/(2*c)
            J1=gamma(2)*gammainc(2,2*c*t)/(2*c)**2
            J2=gamma(3)*gammainc(3,2*c*t)/(2*c)**3
        K=xi+S*t; Q=K*K+h*h
        I=(xi*xi+h*h)*J0+2*xi*S*J1+S*S*J2
        damp=np.exp(-nu*I)
        fhat=zeta*damp/Q
        p=np.sum(fhat)*dk/np.sqrt(2*np.pi)
        px=np.sum(1j*K*fhat)*dk/np.sqrt(2*np.pi)
        pxx=np.sum(-K*K*fhat)*dk/np.sqrt(2*np.pi)
        dil=np.exp(-2*c*t)
        local_R=-h*dil*np.imag(np.conj(p)*px)/2
        local_force=h*np.exp(-3*c*t)*np.imag(np.conj(p)*pxx)/2
        energy=float(dil*np.sum(z2*damp*damp/Q)*dk/4)
        stress=float(-h*dil*np.sum(K*z2*damp*damp/(Q*Q))*dk/2)
        edot=-2*c*energy+shear*stress-nu*dil*dil*np.sum(z2*damp*damp)*dk/2
        data.append((float(t),energy,float(edot),float(local_R),float(rs*local_force),stress))
    a=np.array(data)
    first_break=np.flatnonzero((a[:,2]<=0)|(a[:,4]<=0))
    peak=int(np.argmax(a[:,1]))
    checks=[]
    for tt in [0,.0001,.00025,.0005,.001,.002,.005,.01,.02]:
        if tt<=end:
            row=a[np.argmin(abs(a[:,0]-tt))]
            checks.append(dict(t=row[0],energy_ratio=row[1]/a[0,1],energy_rate=row[2],
                               local_covariance=row[3],initial_radius_times_flat_force=row[4],integrated_covariance=row[5]))
    return dict(status='EXPLORATORY Fourier quadrature of exact linear rotating-affine envelope model; not a rigorous enclosure',
                n=n,length=length,A=A,B=B,h=h,Omega=Om,c=c,nu=nu,
                initial_energy=float(a[0,1]),initial_enstrophy=Z0,
                enstrophy_energy_cap=Z0/(h*h*a[0,1]),
                maximum_sampled_energy_ratio=float(a[peak,1]/a[0,1]),sampled_peak_time=float(a[peak,0]),
                first_sample_break_of_joint_positive_energy_rate_and_flat_force=(None if not len(first_break) else float(a[first_break[0],0])),
                samples=checks,
                scope='flat momentum-flux divergence at the central model ray; not actual cylindrical angular-momentum torque')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=8192)
    ap.add_argument('--length',type=float,default=2.)
    ap.add_argument('--A',type=float,default=960.)
    ap.add_argument('--end',type=float,default=.02)
    ap.add_argument('--steps',type=int,default=2000)
    ap.add_argument('--output',default=str(here/'affine-envelope-spectrum-8192.json'))
    args=ap.parse_args();output=args.output;del args.output
    result=run(**vars(args))
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
