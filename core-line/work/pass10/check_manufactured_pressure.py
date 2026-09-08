#!/usr/bin/env python3
"""Independent compact manufactured-pressure tests for both cylinder solvers.

Analytic source and gradients do not use the initial-field cutoff evaluator.
These are numerical consistency/convergence tests, not continuum enclosures.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse,json,time
from pathlib import Path
import numpy as np
from full_initial_evolution import Cylinder
from high_order_pressure import HighOrderCylinder


def bump_jets(s):
    # f(s)=exp(-s/(1-s)) below 1, zero above; independent C-infinity bump.
    f=np.zeros_like(s);fp=f.copy();fpp=f.copy()
    mask=s<1
    den=1-s[mask]
    v=np.exp(-s[mask]/den)
    f[mask]=v
    fp[mask]=-v/den**2
    fpp[mask]=v*(den**-4-2*den**-3)
    return f,fp,fpp


def manufactured(grid,m,a=1.4):
    r=grid.r[:,None];z=grid.z[None,:];k=2*np.pi/grid.L
    q=1+.3*np.cos(k*z)+.2*np.sin(2*k*z)
    qz=-.3*k*np.sin(k*z)+.4*k*np.cos(2*k*z)
    qzz=-.3*k*k*np.cos(k*z)-.8*k*k*np.sin(2*k*z)
    f,fp,fpp=bump_jets((r/a)**2)
    phase=1 if m==0 else (1+.5j)
    h=(r/a)**m*f
    hr=(m*r**(m-1)*f+2*r**(m+1)*fp/a**2)/a**m
    lapr=(4*(m+1)*r**m*fp/a**2+4*r**(m+2)*fpp/a**4)/a**m
    p=phase*h*q
    source=-phase*(lapr*q+h*qzz)
    grad=phase*np.stack((hr*q,1j*m*h*q/r,h*qz))
    return p,source,grad


def l2(grid,a):
    return float(np.sqrt(np.sum(abs(a)**2*grid.vol[:,None])*grid.L/grid.nz))


def fv_operator_interior(grid,p,m):
    flux=grid.rf[1:-1,None]*(p[1:]-p[:-1])/(grid.r[1:]-grid.r[:-1])[:,None]
    left=np.concatenate((np.zeros_like(flux[:1]),flux[:-1]),axis=0)
    radial=-(flux-left)/grid.vol[:-1,None]
    return radial-grid.dz(p,2)[:-1]+m*m*p[:-1]/grid.r[:-1,None]**2


def one_test(cls,nr,nz=64):
    grid=cls(nr=nr,nz=nz,R=4.,L=6.,mapping=4.)
    out={}
    for m in (0,4,8):
        exact,source,exact_grad=manufactured(grid,m)
        p=grid.poisson(source,m,tag=f'manufactured_m{m}')
        grad=grid.pressure_gradient(p,m)
        divgrad=np.einsum('ii...->...',grid.gradient(grad,m))
        # The source is analytic, not produced by applying the tested operator.
        source_adjusted=source.astype(complex,copy=True)
        if cls is Cylinder and m==0:
            repair=grid.repairs[-1]['uniform_z_compatibility_repair_amplitude']
            # Only this bump is the solver's specified repair profile.
            from initial_full_field import cutoff
            source_adjusted-=(repair[0]+1j*repair[1])*cutoff((grid.r/(grid.R/2))**2)[:,None]
        if cls is Cylinder:
            algebra=fv_operator_interior(grid,p,m)-source_adjusted[:-1]
        else:
            algebra=(-grid.eDrr@p-(grid.eDr@p)/grid.r[:,None]
                     -grid.dz(p,2)+m*m*p/grid.r[:,None]**2-source_adjusted)[:-1]
        # Report both gradient-divergence and the separately stored scalar Laplacian.
        mismatch=source+divgrad
        scalar_mismatch=source+grid.scalar_lap(p,m)
        derivative_only=grid.pressure_gradient(exact,m)-exact_grad
        z0=grid.nz//2
        axis_truth=1.3 if m==0 else 0.
        axis_value=grid.axis(p,0.)
        axis_sampled=grid.axis(exact,0.)
        gradient_scale=float(np.max(abs(exact_grad)))
        radial_scale=(grid.r[0]/1.4)**m*1.3*(1 if m==0 else (1+.5j))
        exact_regular=bump_jets(np.array([(grid.r[0]/1.4)**2]))[0][0]
        out[str(m)]=dict(
            pressure_relative_L2=l2(grid,p-exact)/l2(grid,exact),
            pressure_relative_Linf=float(np.max(abs(p-exact))/np.max(abs(exact))),
            gradient_relative_L2=l2(grid,grad-exact_grad)/l2(grid,exact_grad),
            gradient_relative_Linf=float(np.max(abs(grad-exact_grad))/gradient_scale),
            analytic_sample_gradient_relative_L2=l2(grid,derivative_only)/l2(grid,exact_grad),
            solver_operator_relative_Linf=float(np.max(abs(algebra))/np.max(abs(source))),
            div_pressure_gradient_plus_source_relative_L2=l2(grid,mismatch)/l2(grid,source),
            div_pressure_gradient_plus_source_relative_Linf=float(np.max(abs(mismatch))/np.max(abs(source))),
            scalar_laplacian_plus_source_relative_L2=l2(grid,scalar_mismatch)/l2(grid,source),
            axis_extrapolated_pressure_error=float(abs(axis_value-axis_truth)),
            analytic_axis_extrapolation_error=float(abs(axis_sampled-axis_truth)),
            first_four_rows_gradient_max_error_over_global_scale=float(np.max(abs(grad[:,:4]-exact_grad[:,:4]))/gradient_scale),
            first_cell_regularized_pressure=[float((p[0,z0]/radial_scale).real),float((p[0,z0]/radial_scale).imag)],
            exact_first_cell_regularized_pressure=float(exact_regular),
        )
    return dict(nr=nr,nz=nz,R=grid.R,L=grid.L,mapping=4.,modes=out,repairs_or_flux=grid.repairs)


def run(sizes=(96,192,384,768)):
    started=time.time();results={}
    for cls in (Cylinder,HighOrderCylinder):
        runs=[one_test(cls,n) for n in sizes]
        rates={}
        for mode in ('0','4','8'):
            rates[mode]={}
            for key in ('pressure_relative_L2','gradient_relative_L2',
                        'div_pressure_gradient_plus_source_relative_L2'):
                errors=[run['modes'][mode][key] for run in runs]
                rates[mode][key]=[float(np.log(errors[i]/errors[i+1])/np.log(sizes[i+1]/sizes[i]))
                                  for i in range(len(sizes)-1)]
        results[cls.__name__]=dict(runs=runs,observed_orders=rates)
    checks={}
    for name,res in results.items():
        last=res['runs'][-1];prior=res['runs'][-2]
        checks[name]=dict(
            all_modes_pressure_and_gradient_decrease=all(last['modes'][m][key]<prior['modes'][m][key]
                  for m in ('0','4','8') for key in ('pressure_relative_L2','gradient_relative_L2')),
            all_modes_operator_residual_below_1e_minus_8=all(last['modes'][m]['solver_operator_relative_Linf']<1e-8
                                                               for m in ('0','4','8')),
        )
    return dict(status='COMPLETED numerical manufactured consistency tests',
        exact_profile='p_m=(r/1.4)^m exp(-s/(1-s)) q(z), s=(r/1.4)^2<1, zero outside; q=1+.3 cos(2pi z/6)+.2 sin(4pi z/6)',
        nonzero_mode_complex_multiplier='1+0.5i',
        source='analytic -Delta_m p; no discrete-operator manufacture',
        results=results,checks=checks,seconds=time.time()-started,
        scope='Periodic-z compact-radial scalar pressure tests, m=0,4,8. Numerical convergence only; not a whole-space NS validation, no rigorous truncation-error enclosure. A small solve residual does not imply a commuting divergence/gradient projection.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(Path(__file__).with_name('manufactured-pressure-checks.json')))
    ap.add_argument('--sizes',nargs='+',type=int,default=[96,192,384,768]);args=ap.parse_args()
    result=run(tuple(args.sizes))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(checks=result['checks'],orders={k:v['observed_orders'] for k,v in result['results'].items()},seconds=result['seconds']),indent=2))
