# NOTE — SYNTHESIZER seat — DTC-2026-09-06 / lower

Deliverable: `../SYNTHESIS.md`. This note records only what this seat computed itself.

## Independence declaration
Read: the four attempt records and five refuter records as supplied in the brief;
`../refuter-corollary-blowup-rate/NOTE.md`; `../prove-lagrangian/check_constants.py` and
`check_log.txt`; `../refuter-lagrangian/r1_results.json`;
`../refuter-prove-duhamel/r3_results.json`, `r4_results.json`;
`~/Desktop/TORMENT NEXUS/LAWS.md` (L-53…L-62);
`~/Desktop/Tricritical Exploration/ESTATE_GRADING_2026-08-04.md` §0 (rubric, quoted verbatim).
Did NOT read, import or execute: any solver (`nsring.py`, `gegen5d.py`, `ns5d.py`), any viscous
run JSON, any PDF in `../corollary-blowup-rate/pdf/` or `../refuter2-*/pdf/`, theorem (i)'s
`LOGARITHMIC_RECORD_CLOCK.md`. Nothing outside `lower/` was written; no other seat's files were
modified.

## s1_constants.py — every clock constant recomputed from scratch
Reproduces, to all printed digits, numbers the seats report separately:
- 4log(3/2) = 1.6218604324326575; 8(1−√(2/3)) = 1.4680273525781917; ratio 0.9051502356317253.
- θ: frozen log(3/2) = 0.4054651081081644; accelerated 2(1−√(2/3)) = 0.36700683814454793; Riccati 0.5.
- c₂ = 2θ/κ at κ=1/2: **2.0000000** (Riccati), 1.6218604 (frozen), 1.4680274 (accelerated).
- c₂ at the datum's own κ: 2.001121 (κ=0.49972), **2.077404** (κ=0.48137, the datum the 16 runs
  used), 2.206677 (κ=0.45317).
- Q′ against viscous-numerics 0.988 ± 0.018: frozen z = −3.342/−0.610/**−0.081** at
  κ = 0.437/0.415/0.411; accelerated z = −8.232/−5.758/**−5.280**.
- advertised gain 9.484976 % vs accelerated shortfall at κ=0.411 9.619365 % — the two agree to
  0.13 percentage points, confirming refuter-lagrangian-b's find (G).
- c_E: sharp shell (2/5)log(0.172403978) = −0.7031659387; mollified (2/5)log(0.1418961) =
  −0.7810640717; **mismatch 0.0778981**.
- Corollary form (c) effective coefficient M·s·log(1/s)/c₁ at T−t = 1e-6/1e-12/1e-30/1e-100:
  **3.6713090 / 4.2338579 / 4.6625158 / 4.8937339** — matches refuter2's 3.671/4.234/4.663/4.894
  exactly, by an instrument written here with no sight of theirs.

## s2_gate_recheck.py — the Route-A FL-043 finding verified at source
Static + log read of `../prove-lagrangian/check_constants.py` and `check_log.txt`:
- exactly one line contains a `*0+` mask, and it is the one refuter-lagrangian-b names:
  `ck('sum |H| ||C||/l^2 (lam=1)', s5['sum_absH_over_l2_tail_from_3']*0+0.3958, 0.3958, 1e-9)`;
- the summary is literally `print(("ALL %d CHECKS PASS" % 0) if not fails else ...)`;
- `check_log.txt`: 88 PASS lines, 0 FAIL lines, final line **`ALL 0 CHECKS PASS`**.
Both halves of the finding confirmed independently. (The NOTE claims "88 assertions"; the
source contains 55 `ck(` call sites, several in loops.)

## Deviation
DEV-S1: this seat re-derived no PDE identity and re-ran no solver. Its independent contribution
is arithmetic and static analysis only; every mathematical status in `../SYNTHESIS.md` is
attributed to the seat that established it, and where two seats disagree the disagreement is
printed rather than averaged (L-11).
