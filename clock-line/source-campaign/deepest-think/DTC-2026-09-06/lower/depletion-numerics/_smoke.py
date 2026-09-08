import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_dep import run
from scipy.fft import set_workers
with set_workers(8):
    r = run('shell', 3, label='smoke')
print('T32 frozen 0.6666902575799509 vs', r['cross'].get('T32',{}).get('t'))
print('T2  frozen 1.1305544037651174 vs', r['cross'].get('T2',{}).get('t'))
h = r['hist']
for tt in [0,0.25,0.5,1.0,1.5,2.0,3.0,4.0]:
    print('C(%.2f)=%.4f  om=%.4f' % (tt, np.interp(tt,h['t'],h['C']), np.interp(tt,h['t'],h['ommax'])))
