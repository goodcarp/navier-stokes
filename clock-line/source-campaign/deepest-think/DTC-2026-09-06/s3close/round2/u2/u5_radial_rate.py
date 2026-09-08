#!/usr/bin/env python3
"""
u5 -- the radial-rate lemma and the sharpness of the two maximum principles.

(a) Gamma_rad := sup_x Lambda_max(sym grad_5 b) and the closed form
        Lambda_max = max( a , ( -a + sqrt((3a+2Q)^2 + (2P-om)^2) )/2 ) ,
    Q = r d_r a, P = r d_z a, om = omega^theta   -- verified numerically.
(b) the envelope
        Gamma_rad <= max( A , A + 2 Ghat + lam M/2 , 2 Chat_a M + 2 Ghat + lam M/2 )
    under |a| <= A, a >= -Chat_a M, |Q|,|P| <= Ghat, |om| <= lam M -- verified by
    random search over the admissible box.
(c) SHARPNESS: for the pure axisymmetric strain b = (a r, -2 a z) the two maximum
    principles P1 (factor e^{INT Gamma_rad}) and P2 (factor e^{INT (Gamma + 2 Gamma_rad)})
    reproduce the EXACT transported values lambda and lambda^4.  So neither weight
    is lossy on the reference field; all the loss is in Gamma_rad <= Gamma.
"""
import json, math
import numpy as np

OUT = {}
rng = np.random.default_rng(7761298)

def grad5b(a, Q, P, om):
    Mx = a*np.diag([1., 1., 1., 1., -2.])
    Mx[0, 0] += Q; Mx[0, 4] += P; Mx[4, 0] += P - om; Mx[4, 4] += -Q
    return Mx

# (a) closed form for the largest eigenvalue of the symmetric part
worst = 0.0
for _ in range(200000):
    a, Q, P, om = rng.normal(size=4)*rng.choice([0.1, 1.0, 10.0])
    S = 0.5*(grad5b(a, Q, P, om) + grad5b(a, Q, P, om).T)
    lam_max = np.linalg.eigvalsh(S)[-1]
    closed = max(a, 0.5*(-a + math.sqrt((3*a+2*Q)**2 + (2*P-om)**2)))
    worst = max(worst, abs(lam_max - closed)/max(1.0, abs(lam_max)))
OUT['a_closed_form_max_abs_rel_err'] = worst

# (b) the envelope over the admissible box
def envelope(A, Ca, Gh, lamM):
    return max(A, A + 2*Gh + lamM/2, 2*Ca + 2*Gh + lamM/2)
viol = 0; worst_ratio = 0.0
for _ in range(200000):
    A = 10.0**rng.uniform(-1, 3); Ca = 10.0**rng.uniform(-1, 2)
    Gh = 10.0**rng.uniform(-1, 3); lamM = 10.0**rng.uniform(-1, 1)
    a = rng.uniform(-Ca, A); Q = rng.uniform(-Gh, Gh)
    P = rng.uniform(-Gh, Gh); om = rng.uniform(-lamM, lamM)
    lam_max = max(a, 0.5*(-a + math.sqrt((3*a+2*Q)**2 + (2*P-om)**2)))
    env = envelope(A, Ca, Gh, lamM)
    worst_ratio = max(worst_ratio, lam_max/env)
    if lam_max > env*(1+1e-12):
        viol += 1
OUT['b_envelope_violations'] = viol
OUT['b_envelope_worst_ratio'] = worst_ratio

# (c) sharpness on the pure strain.  b = (a r, -2 a z), a > 0 constant on [0,s].
# Gamma = ||grad_5 b||_op = 2a ; sym grad_5 b = diag(a,a,a,a,-2a) so Gamma_rad = a.
for lam in (1.25, 1.5):
    s_int_a = math.log(lam)                    # INT a = log lambda
    Gamma_int = 2*s_int_a                      # c_G
    P1_factor = math.exp(s_int_a)              # e^{INT Gamma_rad}
    P2_factor = math.exp(Gamma_int + 2*s_int_a)
    OUT[f'c_strain_lam{lam}'] = dict(
        c_G=Gamma_int,
        P1_maxprinciple=P1_factor, P1_exact=lam, P1_rel=abs(P1_factor-lam)/lam,
        P2_maxprinciple=P2_factor, P2_exact=lam**4, P2_rel=abs(P2_factor-lam**4)/lam**4,
        crude_P1_eGamma=math.exp(Gamma_int), crude_P2_e3Gamma=math.exp(3*Gamma_int))
    # direct check that Gamma_rad = a and Gamma = 2a for the strain
    S = 0.5*(grad5b(1.0, 0.0, 0.0, 0.0) + grad5b(1.0, 0.0, 0.0, 0.0).T)
    OUT['c_strain_Gamma_rad_over_Gamma'] = float(np.linalg.eigvalsh(S)[-1]
                                                 / np.linalg.norm(grad5b(1., 0., 0., 0.), 2))

# (d) the exact transported factors for eta_lambda = eta_0 o T_lam^{-1}:
#     |T_lam alpha| <= lam |alpha| and |grad eta_lam| <= lam^2 |grad eta_0|
#     => sup rho|eta_lam| <= lam sup |alpha||eta_0| ,  sup rho^2|grad| <= lam^4 sup
for lam in (1.25, 1.5):
    A = np.diag([1/lam, 1/lam, 1/lam, 1/lam, lam*lam])   # (T_lam)^{-T} = T_lam^{-1}
    OUT[f'd_Tlam_inv_op_lam{lam}'] = float(np.linalg.norm(A, 2))
    OUT[f'd_Tlam_op_lam{lam}'] = float(np.linalg.norm(np.diag([lam]*4+[1/lam**2]), 2))

with open('u5_results.json', 'w') as f:
    json.dump(OUT, f, indent=1)
for k, v in OUT.items():
    print(f'{k:36s} {v}')
