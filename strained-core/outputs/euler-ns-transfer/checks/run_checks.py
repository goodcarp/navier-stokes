"""Run the finite algebra checks with one Python interpreter.

Requires Python 3.9+ and sympy==1.14.0 for the recorded environment.
No source code from the external proof repository is executed.
"""
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
for filename in ['phase_heat_check.py','viscosity_check.py','direction_check.py']:
    print(f'\nRunning {filename}',flush=True)
    subprocess.run([sys.executable,str(HERE/filename)],check=True,cwd=HERE)
print('\nALL FINITE ALGEBRA CHECKS PASSED. No Navier–Stokes existence or blowup is certified.')
