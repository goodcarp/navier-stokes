#!/usr/bin/env python3
from pathlib import Path
import subprocess, sys
here=Path(__file__).resolve().parent
for name in ['verify_quantitative_core_stage.py', 'verify_rotational_receiver.py', 'verify_material_core_transfer.py', 'verify_receiver_flux.py']:
    print("CHECK:",name,flush=True)
    subprocess.run([sys.executable,str(here/name)],check=True)
print("All four packaged checks passed. See VERIFICATION.md for their scope.")
