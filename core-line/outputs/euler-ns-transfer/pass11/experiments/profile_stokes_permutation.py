#!/usr/bin/env python3
"""Read-only comparison of orderings for the exact MAC Stokes saddle matrix.

No evolution, no source solver edits, and no discretization change. Timings are
diagnostics under concurrent load, not performance guarantees or PDE validation.
"""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from scipy.sparse import diags, bmat
from scipy.sparse.linalg import splu
from evolve_full_mac import MAC


def main():
    source = Path(__file__).with_name('evolve_full_mac.py')
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    g = MAC(nr=96, nz=256, R=8., L=16., mapping=4., J=4,
            nu=.001, dealias='padding')
    alpha = (2e-5 / 2) * g.nu / 2
    rng = np.random.default_rng(730129)
    records = []
    for j in range(g.J + 1):
        for ki in [0, 1, 2, 4, 8, 32, 64, 127]:
            n = g.nr
            k = 2 * np.pi * ki / g.L
            D = bmat([[g.Divr, diags(1j*g.modes[j]/g.r),
                       1j*k*diags(np.ones(n))]], format='csc')
            BM = -D.conj().T @ diags(g.V)
            A = diags(g.mass*(1+alpha*k*k)) + alpha*g.K[j]
            full = bmat([[A, BM], [BM.conj().T, None]], format='csc')
            old, right = [], []
            for i in range(n):
                if i > 0:
                    old.append(i-1)
                old.extend((n-1+i, 2*n-1+i))
                right.extend((n-1+i, 2*n-1+i))
                if i < n-1:
                    right.append(i)
                if not (j == 0 and ki == 0 and i == n-1):
                    old.append(g.nv+i)
                    right.append(g.nv+i)
            permutations = {'left_face_before_pressure': np.array(old),
                            'right_face_before_pressure': np.array(right)}
            assert sorted(old) == sorted(right)
            active = np.array(sorted(old))
            rhs = np.zeros((g.nv+g.nr, 2), complex)
            rhs[:g.nv] = rng.normal(size=(g.nv,2))+1j*rng.normal(size=(g.nv,2))
            solutions = {}
            for name, perm in permutations.items():
                mat = full[perm][:,perm]
                t0 = time.perf_counter()
                lu = splu(mat, permc_spec='NATURAL')
                factor_seconds = time.perf_counter()-t0
                # Warm the solve and then time 40 identical two-column solves.
                bp = rhs[perm].copy(order='F')
                xp = lu.solve(bp)
                t0 = time.perf_counter()
                for _ in range(40):
                    xp = lu.solve(bp)
                solve_seconds = (time.perf_counter()-t0)/40
                x = np.zeros_like(rhs)
                x[perm] = xp
                solutions[name] = x
                residual = full[active][:,active] @ x[active]-rhs[active]
                relres = np.linalg.norm(residual)/np.linalg.norm(rhs[active])
                assert relres < 1e-8, (j, ki, name, relres)
                records.append(dict(j=j,m=int(g.modes[j]),k_index=ki,
                    ordering=name, dimension=mat.shape[0],matrix_nnz=mat.nnz,
                    L_nnz=lu.L.nnz,U_nnz=lu.U.nnz,
                    fill_ratio=(lu.L.nnz+lu.U.nnz)/mat.nnz,
                    nnz_per_dimension=(lu.L.nnz+lu.U.nnz)/mat.shape[0],
                    factor_seconds=factor_seconds,
                    two_rhs_solve_seconds=solve_seconds,
                    relative_residual=relres))
            error = np.linalg.norm(solutions['left_face_before_pressure']-
                                   solutions['right_face_before_pressure'])
            scale = np.linalg.norm(solutions['left_face_before_pressure'])
            assert error/scale < 1e-7, (j,ki,error/scale)
            for row in records[-2:]:
                row['relative_ordering_solution_difference'] = error/scale
    summaries = {}
    for name in permutations:
        subset = [r for r in records if r['ordering']==name]
        summaries[name] = dict(
            max_total_LU_nnz=max(r['L_nnz']+r['U_nnz'] for r in subset),
            max_fill_ratio=max(r['fill_ratio'] for r in subset),
            sum_factor_seconds=sum(r['factor_seconds'] for r in subset),
            mean_two_rhs_solve_seconds=float(np.mean([r['two_rhs_solve_seconds'] for r in subset])),
            max_relative_residual=max(r['relative_residual'] for r in subset))
    result = dict(status='PASS: exact matrix comparison; performance diagnostic only',
        parameters=dict(nr=96,nz=256,J=4,dt=2e-5,alpha=alpha),
        source_sha256_before=source_hash,
        source_sha256_after=hashlib.sha256(source.read_bytes()).hexdigest(),
        summaries=summaries,records=records)
    Path(__file__).with_name('stokes-permutation-timing96.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))
    for row in records:
        if row['k_index'] <= 2:
            print(json.dumps(row))


if __name__ == '__main__':
    main()
