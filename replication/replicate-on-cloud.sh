#!/bin/bash
# Independent replication + comparator replay of the two 2026-09-08 Lean releases on a rented CPU box.
# Needs: Ubuntu 22.04/24.04, >= 128 GB RAM, ~32 cores, ~200 GB disk, no GPU. Run as a normal user with sudo.
#   scp this file to the box, then:   bash replicate-on-cloud.sh 2>&1 | tee replicate.log
# Everything lands under ~/replication/<project>/ with build.log, AXIOMS.txt and comparator.log; copy that folder back.
set -u
JOBS=$(( $(nproc) - 2 ))
mkdir -p ~/replication && cd ~/replication
sudo apt-get update -qq && sudo apt-get install -y -qq git curl build-essential python3 zstd >/dev/null
# elan (Lean toolchain manager); each project pins its own toolchain via lean-toolchain
curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none >/dev/null
export PATH="$HOME/.elan/bin:$PATH"
export LEAN_NUM_THREADS=$JOBS

record () {  # $1 = project dir, $2 = axioms script (relative), $3 = label
  local d="$1" ax="$2" label="$3"
  cd "$d" || return 1
  echo "=== $label: $(cat lean-toolchain) start $(date -u +%FT%TZ)" | tee -a build.log
  lake exe cache get >> build.log 2>&1 || echo "(no cache for this toolchain; building Mathlib from source)" | tee -a build.log
  lake build >> build.log 2>&1; echo "=== lake build exit=$? $(date -u +%FT%TZ)" | tee -a build.log
  [ -n "$ax" ] && lake env lean "$ax" > AXIOMS.txt 2>&1 && echo "=== axioms: $(sort AXIOMS.txt | uniq -c | sort -rn | head -3 | tr '\n' ' ')" | tee -a build.log
  cd ~/replication
}

# ---------- 1. Alpöge–Buckmaster: fluid_lean (three projects; Lean 4.32.2; Mathlib from source once, then reused per project) ----------
[ -d fluid_lean ] || git clone -q https://github.com/tristanbuckmaster/fluid_lean.git
( cd fluid_lean && git log -1 --format='fluid_lean commit %H %ad' --date=iso )
record fluid_lean/euler-blowup      scripts/PrintAxioms.lean "AB euler-blowup"
record fluid_lean/boussinesq-blowup scripts/PrintAxioms.lean "AB boussinesq-blowup"
( cd fluid_lean/affinecore && lake build Challenge Solution >> build.log 2>&1 )   # its defaultTargets omit these (audit finding P4)
record fluid_lean/affinecore        scripts/PrintAxioms.lean "AB affinecore"

# ---------- 2. OpenAI: NavierStokesAndEuler (one project; Lean 4.34.0-rc2; README says the Mathlib cache exists) ----------
[ -d NavierStokesAndEuler ] || git clone -q https://github.com/openai/NavierStokesAndEuler.git
( cd NavierStokesAndEuler && git log -1 --format='NavierStokesAndEuler commit %H %ad' --date=iso )
record NavierStokesAndEuler "" "OA NavierStokesAndEuler"
# their axiom printout: see ComparatorChallenges/README.md for the exact command; the challenge files are
# ComparatorChallenges/NavierStokes.lean and ComparatorChallenges/Euler.lean

# ---------- 3. Comparator (statement match + axiom restriction + kernel replay under nanoda) ----------
# Follow each repository's own instructions: fluid_lean/<project>/README.md ("Building and checking") and
# NavierStokesAndEuler/ComparatorChallenges/README.md. The comparator lives at https://github.com/leanprover/comparator;
# its README states the sandboxing it assumes (audit finding P-02/P-03: run the check BEFORE building solution libraries
# unsandboxed, and use the systemd-run RestrictAddressFamilies wrapper it documents). Record each verdict in comparator.log.
[ -d comparator ] || git clone -q https://github.com/leanprover/comparator.git
echo "comparator cloned; run it per each project's README and append verdicts to ~/replication/<project>/comparator.log"
echo "=== done $(date -u +%FT%TZ); copy ~/replication/*/{build.log,AXIOMS.txt,comparator.log} back to the repo's replication/ folder"
