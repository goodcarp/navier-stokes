#!/usr/bin/env python3
from pathlib import Path
import subprocess,sys
here=Path(__file__).resolve().parent
checks=['derive_affine_core_pressure.py', 'verify_affine_pressure_independent.py', 'verify_outer_reservoir.py', 'verify_iteration_audit.py', 'verify_full_feedback_gate.py', 'verify_coefficient_expansion.py', 'derive_core_pressure_derivative.py', 'verify_pressure_viscosity.py', 'verify_pressure_quadrature.py', 'verify_angular_residual_certificate.py']
for name in checks:
    print('CHECK:',name,flush=True)
    args=[sys.executable,str(here/name)]
    if name=='verify_pressure_viscosity.py': args.append('--symbolic-only')
    subprocess.run(args,check=True)
print('All ten packaged checks passed. See VERIFICATION.md for exact scope.')
