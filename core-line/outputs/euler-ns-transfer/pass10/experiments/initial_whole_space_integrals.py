"""Independent whole-space quadrature references for the compact initial datum.

Computes initial core p_zz by source/contact and stress kernels, and actual
3D fluctuation K,K' by full-field tensors and separated profile integrals.
Finite compact supports, both axial endcaps, no finite-cylinder pressure solve.
Gauss comparisons estimate numerical consistency; they are NOT error bounds.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse,json,time
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from initial_full_field import fields,cutoff,_quadratic_composition,_axial_cutoff


def quadrature(order,points):
    nodes,weights=leggauss(order)
    grids=[];wgts=[]
    for lo,hi in zip(points[:-1],points[1:]):
        grids.append((lo+hi)/2+(hi-lo)*nodes/2)
        wgts.append((hi-lo)*weights/2)
    return np.concatenate(grids),np.concatenate(wgts)


def evaluate(order=64,A=1020.,lam=.25,c=1.4,nu=.001):
    started=time.monotonic()
    # Full radial M, integrated in spherical coordinates. Twenty-four angular
    # nodes exactly integrate its polynomial angular dependence, up to roundoff.
    R,wR=quadrature(order,[0.,.5,1.,1.25,2.5,5.,6.])
    mu,wmu=leggauss(24)
    rr=R[:,None]*np.sqrt(1-mu[None,:]**2);zz=R[:,None]*mu[None,:]
    M=fields(rr,zz,A=0.,lam=0.,c=c)
    M_velocity_square=np.sum(abs(M['u'][0])**2,axis=0)
    M_gradient_square=np.sum(abs(M['gradient'][0])**2,axis=(0,1))
    E_M=float(np.pi*np.sum(wR*R**2*(M_velocity_square@wmu)))
    D_M=float(2*np.pi*np.sum(wR*R**2*(M_gradient_square@wmu)))
    # Constant g=4 in the affine ball has zero angular PV contribution.
    # Subtract it there exactly rather than integrate numerical 0/R near R=0.
    gM=M['source'][0]-4*(R[:,None]<.5)
    pzz_M_source=float(-4/3+.5*np.sum(wR/R*(gM@((3*mu**2-1)*wmu))))

    rpoints=[0.,.04,.09,.1575,.19,.24,.315]
    zpoints=[3.4,3.55,3.7,3.775,4.,4.225,4.3,4.45,4.6]
    r,wr=quadrature(order,rpoints);z,wz=quadrature(order,zpoints)
    totals=dict(pzz_packet_source=0.,Cv_stress=0.,CL_stress=0.,
                E_v=0.,D_v=0.,K_tensor=0.,production_tensor=0.,D_seed_tensor=0.)
    for start in range(0,len(z),32):
        Z=z[None,start:start+32];Rho=r[:,None]
        weight=wr[:,None]*wz[None,start:start+32]
        full=fields(Rho,Z,A=A,lam=lam,c=c)
        mean_M=fields(Rho,Z,A=0.,lam=0.,c=c)
        gpacket=full['source'][0]-mean_M['source'][0]
        norm2=Rho**2+Z**2
        # Includes the full 2pi angular integral and two symmetric z packets.
        source_kernel=Rho*(2*Z**2-Rho**2)/norm2**2.5
        totals['pzz_packet_source']+=float(np.sum(weight*source_kernel*gpacket))
        Qr=(12*norm2**2-105*Rho**2*Z**2)/(4*np.pi*norm2**4.5)
        Qt=3*(4*Z**2-Rho**2)/(4*np.pi*norm2**3.5)

        f,fr,_=_quadratic_composition(Rho,0.,63/200,2)
        qv,qvz,_=_axial_cutoff(Z,9/20)
        V=Rho*f*qv;Vr=(f+Rho*fr)*qv;Vz=Rho*f*qvz
        totals['Cv_stress']+=float(np.sum(weight*4*np.pi*Rho*Qt*V**2))
        totals['E_v']+=float(np.sum(weight*2*np.pi*Rho*V**2))
        totals['D_v']+=float(np.sum(weight*4*np.pi*Rho*((f*qv)**2+Vr**2+Vz**2)))
        e,er,_=_quadratic_composition(Rho,7/50,1/10,2)
        qs,_,_=_axial_cutoff(Z,3/5)
        totals['CL_stress']+=float(np.sum(weight*2*np.pi*Rho*qs**2
            *(Qr*16*e**2/Rho**2+Qt*(er**2+400*e**2))))
        u4=full['u'][4];G4=full['gradient'][4];G0=full['gradient'][0]
        totals['K_tensor']+=float(np.sum(weight*4*np.pi*Rho*np.sum(abs(u4)**2,axis=0)))
        work=np.real(np.einsum('i...,ij...,j...->...',u4.conj(),G0,u4))
        totals['production_tensor']-=float(np.sum(weight*8*np.pi*Rho*work))
        totals['D_seed_tensor']+=float(np.sum(weight*8*np.pi*Rho*np.sum(abs(G4)**2,axis=(0,1))))

    # Independent separated identities for the same actual fluctuation energy.
    e,er,err=_quadratic_composition(r,7/50,1/10,2)
    Cprime=2*r/(63/200)**2*cutoff((r/(63/200))**2,1)
    R0=float(np.sum(wr*r*(er**2+(400+16/r**2)*e**2)))
    R1=float(np.sum(wr*r*((err+er/r-(400+16/r**2)*e)**2+400*(2*er+e/r)**2)))
    Ir=float(-np.sum(wr*r*Cprime*e**2))
    qs,qsz,_=_axial_cutoff(z,3/5);qv,_,_=_axial_cutoff(z,9/20)
    Z0=float(2*np.sum(wz*qs**2));Z1=float(2*np.sum(wz*qsz**2));Iz=float(2*np.sum(wz*qs**2*qv))
    K_separated=np.pi*lam**2*R0*Z0/2
    D_seed_separated=np.pi*lam**2*(R1*Z0+R0*Z1)
    production_separated=np.pi*lam**2*(80*A*Ir*Iz-c*R0*Z0)
    Kprime_separated=production_separated-nu*D_seed_separated
    pzz_source=pzz_M_source+totals['pzz_packet_source']
    pzz_stress=-76/35-A*A*totals['Cv_stress']-lam*lam*totals['CL_stress']
    total_energy=E_M+A*A*totals['E_v']+totals['K_tensor']
    total_D=D_M+A*A*totals['D_v']+totals['D_seed_tensor']
    Kprime_tensor=totals['production_tensor']-nu*totals['D_seed_tensor']
    return dict(status='EXPLORATORY WHOLE-SPACE QUADRATURE REFERENCE; NOT AN INTERVAL CERTIFICATE',
        parameters=dict(order_per_piece=order,A=A,lam=lam,c=c,nu=nu),
        pzz_source=pzz_source,pzz_stress=pzz_stress,pzz_M_source=pzz_M_source,
        pzz_M_exact=-76/35,pzz_source_stress_difference=pzz_source-pzz_stress,
        beta_prime_from_stress=-4-pzz_stress/2,
        E_M=E_M,D_M=D_M,total_energy=total_energy,total_velocity_gradient_L2_squared=total_D,
        total_energy_prime=-nu*total_D,
        Kprime_tensor=Kprime_tensor,Kprime_separated=Kprime_separated,
        K_separated=K_separated,D_seed_separated=D_seed_separated,
        production_separated=production_separated,
        Kprime_over_K=Kprime_tensor/totals['K_tensor'],
        separated_moments=dict(R0=R0,R1=R1,Ir=Ir,Z0=Z0,Z1=Z1,Iz=Iz),
        seconds=time.monotonic()-started,**totals)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--orders',type=int,nargs='+',default=[32,48,72])
    ap.add_argument('--output',default=str(Path(__file__).with_name('initial-whole-space-integrals.json')))
    args=ap.parse_args()
    results=[]
    for order in args.orders:
        result=evaluate(order)
        results.append(result)
        print(json.dumps(result),flush=True)
    Path(args.output).write_text(json.dumps(dict(status='QUADRATURE COMPARISON, NOT VALIDATED ERROR BOUNDS',results=results),indent=2)+'\n')
