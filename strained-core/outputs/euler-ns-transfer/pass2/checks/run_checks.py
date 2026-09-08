#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
here=Path(__file__).resolve().parent
for name in ['verify_active_core_pressure.py', 'verify_local_pressure_audit.py', 'verify_affine_stretch_gate.py', 'verify_correction_window.py']:
    print("CHECK:",name,flush=True)
    subprocess.run([sys.executable,str(here/name)],check=True)
print("All packaged checks passed; see VERIFICATION.md for scope.")
