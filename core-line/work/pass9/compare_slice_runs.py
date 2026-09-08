#!/usr/bin/env python3
"""Summarize numerical sensitivity; this does not certify a PDE solution."""
import argparse,json
from pathlib import Path

def run(folder):
    names=['slice-768.json','slice-1536.json','slice-angular-check.json','slice-boundary-check.json',
           'slice-fixed-A1020.json','slice-fixed-A1020-coarse.json']
    data={name:json.loads((folder/name).read_text()) for name in names}
    out={}
    for name,d in data.items():
        records=d['records']; first=records[0]
        target=min(records,key=lambda x:abs(x['t']-.001)); assert abs(target['t']-.001)<1e-12
        peak=max(records,key=lambda x:x['mean_maximum'])
        out[name]=dict(parameters=d['parameters'],at_t_001=dict(
            branch_value_increment=target['mean_maximum']-first['mean_maximum'],
            branch_value_fractional_gain=target['mean_maximum']/first['mean_maximum']-1,
            branch_radius_displacement=target['mean_maximum_radius']-first['mean_maximum_radius'],
            fixed_receiver_increment=target['mean_fixed_initial_receiver']-first['mean_fixed_initial_receiver'],
            fluctuation_energy_fractional_gain=target['fluctuation_energy_per_axial_length']/first['fluctuation_energy_per_axial_length']-1,
            material_slice_weighted_energy_fractional_gain=target['material_slice_weighted_fluctuation_energy']/first['fluctuation_energy_per_axial_length']-1,
            torque=target['torque_at_radial_maximum']),
            sampled_peak_time=peak['t'],sampled_peak_increment=peak['mean_maximum']-first['mean_maximum'],
            last_fixed_receiver_increment=records[-1]['mean_fixed_initial_receiver']-first['mean_fixed_initial_receiver'],
            final_torque=records[-1]['torque_at_radial_maximum'],
            relative_energy_balance_diagnostic=d['total_energy_balance_relative_residual'])
    coarse=out['slice-768.json']['at_t_001'];fine=out['slice-1536.json']['at_t_001']
    diffs={key:abs(fine[key]-coarse[key]) for key in fine}
    assert 0.0019<fine['branch_value_fractional_gain']<.0021
    assert .12<fine['fluctuation_energy_fractional_gain']<.15
    assert out['slice-1536.json']['final_torque']<0
    assert out['slice-1536.json']['last_fixed_receiver_increment']<0
    assert diffs['branch_value_increment']<.0005
    assert diffs['fluctuation_energy_fractional_gain']<.001
    for name in ['slice-angular-check.json','slice-boundary-check.json']:
        test=out[name]['at_t_001']
        assert abs(test['branch_value_increment']-coarse['branch_value_increment'])<1e-6
        assert abs(test['fluctuation_energy_fractional_gain']-coarse['fluctuation_energy_fractional_gain'])<1e-5
    chosen=out['slice-fixed-A1020.json']['at_t_001']
    chosen_coarse=out['slice-fixed-A1020-coarse.json']['at_t_001']
    chosen_diffs={key:abs(chosen[key]-chosen_coarse[key]) for key in chosen}
    assert .0017<chosen['branch_value_fractional_gain']<.0018
    assert .14<chosen['fluctuation_energy_fractional_gain']<.15
    assert chosen_diffs['branch_value_increment']<.0005
    assert chosen_diffs['fluctuation_energy_fractional_gain']<.001
    return dict(status='PASS numerical sensitivity checks',scope='Tracked radial spline critical branch and finite-radius planar energies; no compact 3D evolution or rigorous continuum error enclosure',runs=out,radial_and_time_refinement_differences=diffs,selected_amplitude_refinement_differences=chosen_diffs)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--folder',default=str(Path(__file__).resolve().parent));ap.add_argument('--output')
    args=ap.parse_args();result=run(Path(args.folder));text=json.dumps(result,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    print(text)
