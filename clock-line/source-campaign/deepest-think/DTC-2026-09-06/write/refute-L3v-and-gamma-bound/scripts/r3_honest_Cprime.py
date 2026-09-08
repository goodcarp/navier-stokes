"""r3 - what C' actually is when the ONLY input for ||f||_inf is the PROVED envelope
(Corollary A1), instead of the COMPUTED value of the L3v sum.

The target seat's column labelled "C' (all inputs proved)" (JSON key C_prime_fully_proved)
substitutes  m['L3v_sum']  -- the g5 COMPUTED value 0.4280...0.7911, whose tail contains a
fitted residual-decay extrapolation -- into (B.5).  Corollary A1, the thing Part A actually
PROVES, gives ||f||_inf <= 0.899159 sqrt(2/pi) (2 lam M + V_lam), i.e. 4.30...6.46 lam M.
This script recomputes C' with that envelope, changing nothing else.
"""
import json, numpy as np

g8 = json.load(open('../copy/g8_results.json'))
g9 = json.load(open('../copy/g9_results.json'))
g6 = json.load(open('../copy/g6_results.json'))

Kf = 3 * np.pi / 4
Kw = 2.0 / 3.0
tailsum = g9['sum_l_ge_3_l^-1.5']
ratio = g9['C_ratio_sup']
A1coef = tailsum * ratio          # 0.899159

out = {}
rows = {}
for k, v in g8['rows'].items():
    lam = float(k.split('_')[0][3:]); dd = k.split('del')[1]
    m, pr = v['measured'], v['proved_pieces']
    # the PROVED envelope for ||f||_inf, from Corollary A1 with the seat's own computed V_lam
    key9 = k
    V = g9['rows'][key9]['V'] if key9 in g9['rows'] else None
    Cb_meas_V = np.sqrt(2/np.pi) * (2*lam + V)                  # uses the computed V
    Cb_env_V  = np.sqrt(2/np.pi) * (2*lam + 6*lam)              # uses only the PROVED V <= 6 lam M
    f_A1_measV = A1coef * Cb_meas_V
    f_A1_envV  = A1coef * Cb_env_V
    for tag, fnorm in [('A1_with_computed_V', f_A1_measV), ('A1_with_proved_V', f_A1_envV)]:
        adev = Kf * fnorm + Kw * pr['sup_omega_eff_env']
        rgrad = 3 * fnorm + 3 * adev + pr['sup_s2_Wdev_env']
        Cp = (2 * (pr['a1_offset'] + adev + m['edge_a'] + pr['collar_ellip_a'])
              + (pr['r_grad_a1'] + rgrad + m['edge_grad'] + pr['collar_ellip_grad']) + lam)
        rows.setdefault(k, {})[tag] = {'f_norm': float(fnorm), 'adev': float(adev),
                                       'r_grad_adev': float(rgrad), 'C_prime': float(Cp),
                                       'crossover_L_10pct': float(Cp / (0.2 * m['kappa']))}
    rows[k]['seat_C_prime_fully_proved'] = v['C_prime_fully_proved']
    rows[k]['seat_L3v_sum_used_as_f_norm'] = m['L3v_sum']
    rows[k]['kappa'] = m['kappa']
    print('%-18s  seat"fully proved" C\'=%6.2f (uses ||f||=%.4f)   HONEST A1 C\'=%7.2f (||f||=%.4f)  '
          'A1+provedV C\'=%7.2f (||f||=%.4f)'
          % (k, v['C_prime_fully_proved'], m['L3v_sum'],
             rows[k]['A1_with_computed_V']['C_prime'], rows[k]['A1_with_computed_V']['f_norm'],
             rows[k]['A1_with_proved_V']['C_prime'], rows[k]['A1_with_proved_V']['f_norm']))

out['rows'] = rows
out['seat_C_full_max'] = max(v['C_prime_fully_proved'] for v in g8['rows'].values())
out['honest_A1_computedV_C_max'] = max(r['A1_with_computed_V']['C_prime'] for r in rows.values())
out['honest_A1_provedV_C_max'] = max(r['A1_with_proved_V']['C_prime'] for r in rows.values())
out['honest_A1_computedV_crossover_max'] = max(r['A1_with_computed_V']['crossover_L_10pct'] for r in rows.values())
out['honest_A1_provedV_crossover_max'] = max(r['A1_with_proved_V']['crossover_L_10pct'] for r in rows.values())
out['honest_A1_computedV_crossover_min'] = min(r['A1_with_computed_V']['crossover_L_10pct'] for r in rows.values())
out['inflation_factor_max'] = out['honest_A1_computedV_C_max'] / out['seat_C_full_max']
print()
print(json.dumps({k: out[k] for k in out if k != 'rows'}, indent=1))
json.dump(out, open('r3_results.json', 'w'), indent=1)
