"""Compare stored production diagnostics against independent physical integrals.

These checks do not replay a full production grid or certify its continuum
error. Deliberately retain and quantify the failed projection/core checks.
"""
import json
from pathlib import Path

here=Path(__file__).resolve().parent
def read(name):return json.loads((here/name).read_text())

whole=read('initial-whole-space-integrals.json')['results'][-1]
assert abs(whole['pzz_source']-whole['pzz_stress'])<1e-9
assert abs(whole['Kprime_tensor']-whole['Kprime_separated'])<1e-8
images=read('initial-periodic-images.json')
image_by_period={x['period']:x for x in images['results'] if x['positive_image_pairs']==32}
comparisons=[]
for period,name in [(16,'full-coefficients-fd768.json'),(32,'full-coefficients-fd768-L32.json')]:
    data=read(name);reference=image_by_period[period]
    init=data['initial']
    error=init['pressure_zz_at_core']-reference['periodic_pzz_reference']
    assert abs(error)<1e-4
    assert abs(init['fluctuation_energy_derivative']-whole['Kprime_tensor'])<1e-3
    assert abs(init['energy_derivative']-whole['total_energy_prime'])<.2
    # A tiny scalar solve residual would not justify discarding this failure.
    assert init['numerical_div_N0_L2']>100
    comparisons.append(dict(period=period,pressure_reference_discrepancy=error,
        periodic_image_contribution=reference['correction_stress'],
        rigorous_image_tail_bound=reference['rigorous_image_tail_bound_float'],
        Kprime_reference_discrepancy=init['fluctuation_energy_derivative']-whole['Kprime_tensor'],
        unresolved_divergence_L2=init['numerical_div_N0_L2']))
coarse=read('full-initial-jets-audit384.json')
fine=read('full-coefficients-fd768-L32.json')
assert abs(fine['direct_core_jet_checks']['b2_pressure_discrepancy'])>1
manufactured=read('manufactured-pressure-checks.json')
assert all(all(checks.values()) for checks in manufactured['checks'].values())
print(json.dumps(dict(status='PASS: independent numerical reference comparisons and explicit failure retention',
    comparisons=comparisons,
    unresolved_fine_b2_discrepancy=fine['direct_core_jet_checks']['b2_pressure_discrepancy'],
    scope='Stored runs compared; no evolved endpoint, commuting projection, residual bound or continuum accuracy certificate.'),indent=2))
