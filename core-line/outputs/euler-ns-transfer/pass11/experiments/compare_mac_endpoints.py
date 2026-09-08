#!/usr/bin/env python3
"""Compare completed, same-grid run_resolved_mac checkpoints.

Retains complex angular phases. All norms/observables are finite-cylinder
quantities; reference initial values below are numerical, not interval bounds.
"""
import argparse,json,math
from pathlib import Path
import numpy as np
from scipy.fft import fft
from scipy.interpolate import CubicSpline

ALLOWED_DIFFERENCES={'dt','steps','fft_workers','workers'}
GRID_KEYS={'nr','nz','R','L','mapping','J','nu','T','dealias','backend'}
REFERENCE={'K':6.4437717665751,'mean_G':41.87025495048,'core_b':1.,'core_Omega':1.}
OBS_KEYS=('energy','K','core_b','core_Omega','core_beta','mean_G_fixed','mean_G_branch','mean_G_branch_radius')


def fail(message):raise ValueError(message)
def equal_parameters(a,b,label,allowed=frozenset()):
    if set(a)!=set(b):fail(f'{label}: parameter key sets differ')
    bad={k:(a[k],b[k]) for k in a if k not in allowed and a[k]!=b[k]}
    if bad:fail(f'{label}: incompatible parameters {bad}')

def load_endpoint(json_path,npz_path):
    report=json.loads(Path(json_path).read_text())
    if report.get('complete') is not True:fail(f'{json_path}: run is not complete')
    p=report.get('parameters',{})
    if not GRID_KEYS<=set(p):fail(f'{json_path}: required parameters missing')
    if not {'dt','steps'}<=set(p):fail('Missing time-step metadata')
    if not isinstance(p['steps'],int) or p['steps']<1 or not p['dt']>0:
        fail('Invalid time-step metadata')
    if not math.isclose(p['dt']*p['steps'],p['T'],rel_tol=2e-14,abs_tol=2e-14):
        fail('Time step times step count differs from T')
    if p['dealias']!='padding' or p['backend']!='fast':fail('Expected run_resolved_mac fast/padding metadata')
    with np.load(npz_path,allow_pickle=False) as saved:
        for name in ('U','r','rf','z','parameters'):
            if name not in saved:fail(f'{npz_path}: missing {name}')
        pp=json.loads(str(saved['parameters'].item()))
        equal_parameters(p,pp,'JSON/checkpoint',allowed=frozenset())
        U=np.array(saved['U'],dtype=np.complex128);r=np.array(saved['r']);rf=np.array(saved['rf']);z=np.array(saved['z'])
    nr,nz,J=int(p['nr']),int(p['nz']),int(p['J'])
    if nr<4 or nz<7 or U.shape!=(J+1,3*nr-1,nz):fail('Invalid field shape/grid size')
    if r.shape!=(nr,) or rf.shape!=(nr+1,) or z.shape!=(nz,):fail('Invalid coordinate shapes')
    if not all(np.isfinite(x).all() for x in (U,r,rf,z)):fail('Nonfinite checkpoint')
    if np.any(np.diff(rf)<=0) or rf[0]!=0 or rf[-1]!=p['R']:fail('Invalid radial faces')
    if not np.array_equal(r,(rf[:-1]+rf[1:])/2):fail('Cell radii disagree with saved faces')
    expected_z=(np.arange(nz)-nz//2)*p['L']/nz
    if not np.array_equal(z,expected_z):fail('Physical axial origin/grid does not match the driver')
    expected_rf=p['R']*np.sinh(p['mapping']*np.arange(nr+1)/nr)/np.sinh(p['mapping'])
    if not np.array_equal(rf,expected_rf):fail('Radial coordinates disagree with grid parameters')
    if abs(U[0].imag).max()>1e-13*max(1.,abs(U[0]).max()):fail('Mode zero is not a real physical axial field')
    history=report.get('history',[])
    if len(history)<2 or history[0].get('time')!=0:fail('Missing initial/final history')
    if not math.isclose(history[-1].get('time',math.nan),p['T'],rel_tol=2e-14,abs_tol=2e-14):fail('Final history time differs from T')
    V=(rf[1:]**2-rf[:-1]**2)/2;W=rf[1:-1]*np.diff(r)
    mass=np.r_[W,V,V]
    obs=observables(U,p,r,z,mass)
    for key in OBS_KEYS:
        stored=history[-1].get(key)
        if not isinstance(stored,(int,float)) or not math.isfinite(stored):fail(f'Invalid saved endpoint {key}')
        if obs[key] is None or not math.isclose(stored,obs[key],rel_tol=2e-10,abs_tol=2e-10):
            fail(f'Checkpoint/JSON endpoint mismatch in {key}: {obs[key]} versus {stored}')
    return dict(parameters=p,U=U,r=r,rf=rf,z=z,mass=mass,report=report,
        observables=obs,json_path=str(Path(json_path)),npz_path=str(Path(npz_path)))


def per_mode_squared_norm(U,p,mass):
    factor=np.r_[1.,np.full(int(p['J']),2.)]
    return 2*np.pi*p['L']/p['nz']*factor*np.sum(mass[None,:,None]*abs(U)**2,axis=(1,2))


def observables(U,p,r,z,mass):
    nr,nz=int(p['nr']),int(p['nz']);a=U[0].real
    st=slice(nr-1,2*nr-1);sz=slice(2*nr-1,3*nr-1)
    xx=r[:4]**2;scale=xx[-1]
    axis=np.linalg.solve(np.array([(xx/scale)**j for j in range(4)]),[1.,0,0,0])
    zc=nz//2;omega=float(axis@(a[st][:4,zc]/r[:4]))
    h=p['L']/nz;uz=a[sz]
    dz=(uz[:,zc+3]-9*uz[:,zc+2]+45*uz[:,zc+1]-45*uz[:,zc-1]+9*uz[:,zc-2]-uz[:,zc-3])/(60*h)
    b=.5*float(axis@dz[:4])
    ki=np.rint(np.fft.fftfreq(nz)*nz).astype(int);kz=2*np.pi*ki/p['L']
    phase=np.exp(1j*kz*(4-z[0]))
    receiver=np.real(fft(a[st],axis=-1,norm='forward')@phase)
    spline=CubicSpline(r,r*receiver);rstar=.21547557533348247
    candidates=[float(x) for x in spline.derivative().roots(extrapolate=False) if .18<x<.3]
    rm=max(candidates,key=lambda x:float(spline(x))) if candidates else rstar
    norms=per_mode_squared_norm(U,p,mass)
    return dict(energy=float(norms.sum()/2),K=float(norms[1:].sum()/2),core_b=b,core_Omega=omega,
        core_beta=b/omega if omega!=0 else None,mean_G_fixed=float(spline(rstar)),
        mean_G_branch=float(spline(rm)),mean_G_branch_radius=rm)


def safe_ratio(a,b):return float(a/b) if b>0 else (0. if a==0 else None)

def gain_margins(e):
    o=e['observables'];initial=e['report']['history'][0]
    discrete={key:o[key]-initial[key] for key in ('K','mean_G_fixed','mean_G_branch','core_b','core_Omega','core_beta')}
    return dict(same_run_discrete_initial_differences=discrete,
        numerical_whole_space_initial_comparisons=dict(K_minus_reference=o['K']-REFERENCE['K'],
            fixed_mean_minus_reference=o['mean_G_fixed']-REFERENCE['mean_G'],
            branch_mean_minus_reference=o['mean_G_branch']-REFERENCE['mean_G'],
            b_minus_one=o['core_b']-1.,Omega_minus_one=o['core_Omega']-1.,
            b_minus_Omega=o['core_b']-o['core_Omega'],beta_minus_one=o['core_beta']-1.),
        reference_values=REFERENCE,scope='Numerical diagnostic differences only; no interval lower bounds or whole-space endpoint validation.')


def compare(json_a,npz_a,json_b,npz_b):
    a=load_endpoint(json_a,npz_a);b=load_endpoint(json_b,npz_b)
    equal_parameters(a['parameters'],b['parameters'],'Run A/run B',ALLOWED_DIFFERENCES)
    for key in ('r','rf','z'):
        if not np.array_equal(a[key],b[key]):fail(f'Incompatible saved {key} coordinates')
    if a['report'].get('source_sha256')!=b['report'].get('source_sha256'):
        fail('Source hashes differ; this is not an isolated time-step comparison')
    # The driver deterministically samples the same initial datum on this grid.
    # Different initial observables indicate an incompatible comparison.
    for key in OBS_KEYS:
        ia=a['report']['history'][0].get(key);ib=b['report']['history'][0].get(key)
        if ia!=ib:fail(f'Initial history differs in {key}')
    p=a['parameters'];mass=a['mass']
    na=per_mode_squared_norm(a['U'],p,mass);nb=per_mode_squared_norm(b['U'],p,mass)
    nd=per_mode_squared_norm(b['U']-a['U'],p,mass)
    differences=[dict(m=4*j,L2_absolute=float(np.sqrt(d)),L2_relative_to_B=safe_ratio(np.sqrt(d),np.sqrt(bb)),
        energy_A=float(aa/2),energy_B=float(bb/2),energy_B_minus_A=float((bb-aa)/2))
        for j,(aa,bb,d) in enumerate(zip(na,nb,nd))]
    result=dict(status='SAME-GRID DISCRETE ENDPOINT COMPARISON; NOT A CONTINUUM ERROR BOUND',
        orientation='Signed observable differences are B minus A; relative norms use B as denominator.',
        definitions={'full_L2':'sqrt(2*pi*L/Nz * sum_j angular_factor_j * sum mass * |B_j-A_j|^2)',
            'seed_L2':'The same norm restricted to all nonzero angular modes, including generated harmonics.',
            'angular_factor':'1 at m=0, 2 for positive m; complex coefficient differences retain phase.',
            'receiver':'The driver cubic radial interpolation at physical z=4; branch candidates have .18<r<.3, not a global maximum.'},
        A=dict(json=a['json_path'],checkpoint=a['npz_path'],parameters=a['parameters'],observables=a['observables'],gain_margins=gain_margins(a)),
        B=dict(json=b['json_path'],checkpoint=b['npz_path'],parameters=b['parameters'],observables=b['observables'],gain_margins=gain_margins(b)),
        full_L2_absolute=float(np.sqrt(nd.sum())),seed_L2_absolute=float(np.sqrt(nd[1:].sum())),
        full_L2_A=float(np.sqrt(na.sum())),full_L2_B=float(np.sqrt(nb.sum())),
        seed_L2_A=float(np.sqrt(na[1:].sum())),seed_L2_B=float(np.sqrt(nb[1:].sum())),
        full_L2_relative_to_B=safe_ratio(np.sqrt(nd.sum()),np.sqrt(nb.sum())),
        seed_L2_relative_to_B=safe_ratio(np.sqrt(nd[1:].sum()),np.sqrt(nb[1:].sum())),
        per_mode=differences,observables_B_minus_A={key:b['observables'][key]-a['observables'][key] for key in OBS_KEYS},
        source_sha256=a['report'].get('source_sha256'),
        limitation='Same-grid time-step agreement does not bound spatial/domain errors or the whole-space residual; the numerical initial references are not intervals.')
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--a-json',required=True);ap.add_argument('--a-npz',required=True)
    ap.add_argument('--b-json',required=True);ap.add_argument('--b-npz',required=True)
    ap.add_argument('--output')
    args=ap.parse_args()
    try:out=compare(args.a_json,args.a_npz,args.b_json,args.b_npz)
    except (ValueError,KeyError,OSError) as exc:ap.exit(2,f'Endpoint comparison rejected: {exc}\n')
    serialized=json.dumps(out,indent=2,allow_nan=False)+'\n'
    if args.output:Path(args.output).write_text(serialized)
    print(serialized,end='')
