#!/usr/bin/env python3
"""Small manufactured checks; does not load or rerun a research trajectory."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import hashlib,json,tempfile
from pathlib import Path
import numpy as np
from evolve_full_mac import MAC
from compare_mac_endpoints import compare,per_mode_squared_norm

HERE=Path(__file__).resolve().parent

def main():
    results=[]
    for nz in (32,33):
        g=MAC(nr=16,nz=nz,R=1.,L=10.,mapping=1.2,J=2,nu=.001,dealias='padding')
        r=g.r[:,None];z=g.z[None,:];k=2*np.pi/g.L
        U=np.zeros((3,g.nv,nz),complex)
        U[0,g.slt]=r;U[0,g.slz]=2*np.sin(k*z)/k
        radial=g.rf[1:-1,None]
        U[1,g.slr]=.07*(1+1j)*radial**3*np.exp(2j*k*z)
        U[1,g.slt]=.05*(1-2j)*r**3*np.exp(-3j*k*z)
        U[1,g.slz]=.04*r**4*np.exp(1j*k*z)
        V=U.copy();theta=.7;eps=.031
        V[1]*=np.exp(1j*theta)
        V[2,g.slt]=.02*(1+3j)*r**7*np.exp(-2j*k*z)
        V[0,g.slt]+=eps*r*np.cos(k*z)
        pa=dict(nr=g.nr,nz=g.nz,R=g.R,L=g.L,mapping=1.2,J=g.J,nu=g.nu,T=.01,
            dt=.002,steps=5,dealias='padding',backend='fast',fft_workers=1)
        pb=dict(pa,dt=.001,steps=10,fft_workers=2)
        initial=g.observe(U,0.)
        source={'manufactured':'No evolved trajectory; same source identifier on both fixtures.'}
        def report(X,p):return dict(complete=True,parameters=p,source_sha256=source,
            history=[initial,g.observe(X,p['T'])])
        with tempfile.TemporaryDirectory(prefix='mac-endpoint-check-') as tmp:
            tmp=Path(tmp);aj=tmp/'a.json';az=tmp/'a.npz';bj=tmp/'b.json';bz=tmp/'b.npz'
            a=report(U,pa);b=report(V,pb)
            def save_b(rep,arr=V,coords=None,params=None):
                bj.write_text(json.dumps(rep))
                c=dict(r=g.r,rf=g.rf,z=g.z) if coords is None else coords
                np.savez(bz,U=arr,parameters=json.dumps(rep['parameters'] if params is None else params),**c)
            aj.write_text(json.dumps(a));np.savez(az,U=U,r=g.r,rf=g.rf,z=g.z,parameters=json.dumps(pa))
            save_b(b);out=compare(aj,az,bj,bz)
            # Direct independent mode-weighted sums and phase identity.
            pref=2*np.pi*g.L/g.nz
            e_seed_A=pref*2*np.sum(g.mass[:,None]*abs(U[1])**2)
            e_mode8=pref*2*np.sum(g.mass[:,None]*abs(V[2])**2)
            axis_error=eps*r*np.cos(k*z)
            e_axis=pref*np.sum(g.V[:,None]*abs(axis_error)**2)
            seed_expected=(2-2*np.cos(theta))*e_seed_A+e_mode8
            total_expected=seed_expected+e_axis
            errors=dict(full_L2=abs(out['full_L2_absolute']**2-total_expected)/total_expected,
                seed_L2=abs(out['seed_L2_absolute']**2-seed_expected)/seed_expected,
                seed_relative=abs(out['seed_L2_relative_to_B']-np.sqrt(seed_expected/(e_seed_A+e_mode8))),
                core_Omega_difference=abs(out['observables_B_minus_A']['core_Omega']-eps),
                physical_z_receiver_difference=abs(out['observables_B_minus_A']['mean_G_fixed']-
                    eps*.21547557533348247**2*np.cos(k*4)),
                K_difference=abs(out['observables_B_minus_A']['K']-e_mode8/2))
            assert max(errors.values())<2e-12,errors
            # Same coefficient magnitudes do NOT imply zero state difference.
            assert out['per_mode'][1]['L2_absolute']>.01
            assert abs(out['per_mode'][1]['energy_B_minus_A'])<1e-14
            assert out['per_mode'][2]['L2_relative_to_B']==1.
            rejections=[]
            def reject(name,rep=b,arr=V,coords=None,params=None):
                save_b(rep,arr,coords,params)
                try:compare(aj,az,bj,bz)
                except ValueError as exc:rejections.append(dict(case=name,reason=str(exc)))
                else:raise AssertionError(f'{name} was accepted')
            def changed(name,value):
                rep=dict(b);rep['parameters']=dict(pb);rep['parameters'][name]=value;return rep
            reject('viscosity',changed('nu',.002))
            reject('final_time',changed('T',.02))
            reject('angular_truncation',changed('J',3))
            reject('dealias',changed('dealias','filter'))
            reject('unknown_extra_parameter',changed('datum_amplitude',1021))
            reject('JSON_npz_dt_disagreement',params=dict(pb,dt=.003))
            reject('physical_z_origin',coords=dict(r=g.r,rf=g.rf,z=g.z+.001))
            rr=g.r.copy();rr[3]+=.001
            reject('radial_coordinates',coords=dict(r=rr,rf=g.rf,z=g.z))
            reject('incomplete',rep=dict(b,complete=False))
            reject('different_source',rep=dict(b,source_sha256={'manufactured':'different'}))
            bad=V.copy();bad[0,0,0]=np.nan
            reject('nonfinite',arr=bad)
            bad=V.copy();bad[0,0,0]=1j
            reject('nonreal_mode_zero',arr=bad)
            bad=V.copy();bad[0,g.slt]+=r
            reject('checkpoint_history_mismatch',arr=bad)
            initial_bad=dict(initial);initial_bad['K']+=1
            reject('different_initial_observable',rep=dict(b,history=[initial_bad,b['history'][-1]]))
            results.append(dict(nz=nz,analytical_errors=errors,rejection_tests=rejections,
                full_L2=out['full_L2_absolute'],seed_L2=out['seed_L2_absolute'],
                seed_relative=out['seed_L2_relative_to_B'],phase_preserved=True))
    hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
        ('compare_mac_endpoints.py','check_compare_mac_endpoints.py','evolve_full_mac.py')}
    result=dict(status='PASS',scope='Manufactured saved fields and metadata only; no NS run used.',
        source_sha256=hashes,results=results)
    path=HERE/'compare-mac-endpoints-checks.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status='PASS',manufactured_grid_count=len(results),rejection_count=sum(len(r['rejection_tests']) for r in results),output=str(path))))

if __name__=='__main__':main()
