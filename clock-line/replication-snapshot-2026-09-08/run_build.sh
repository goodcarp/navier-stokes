#!/bin/bash
# Independent replication of the Alpöge–Buckmaster Euler Lean certificate (euler-blowup), laptop-throttled.
export PATH="$HOME/.elan/bin:$PATH"
export LEAN_NUM_THREADS=4
cd "$(dirname "$0")"
LOG=build.log
echo "=== replication start $(date -u +%FT%TZ) host=$(hostname) threads=$LEAN_NUM_THREADS" >> $LOG
# memory watchdog: kill the build if free+inactive memory stays under 1.2 GB for two checks in a row
( low=0; while true; do
    sleep 30
    pgrep -f "lake build" >/dev/null || break
    pg=$(vm_stat | awk '/page size of/{gsub("[^0-9]","",$8); print $8}')
    free=$(vm_stat | awk -v pg="$pg" '/Pages free/{f=$3} /Pages inactive/{i=$3} /Pages speculative/{s=$3} END{gsub("\\.","",f);gsub("\\.","",i);gsub("\\.","",s); printf "%d", (f+i+s)*pg/1e9}')
    if [ "$free" -lt 2 ]; then low=$((low+1)); else low=0; fi
    if [ $low -ge 2 ]; then echo "WATCHDOG: free memory ${free} GB twice, killing lake $(date -u +%FT%TZ)" >> $LOG; pkill -f "lake build"; pkill -x lean; break; fi
  done ) &
WD=$!
nice -n 19 lake build >> $LOG 2>&1
RC=$?
echo "=== lake build exit=$RC $(date -u +%FT%TZ)" >> $LOG
kill $WD 2>/dev/null
if [ $RC -eq 0 ]; then
  echo "=== axioms:" >> $LOG
  nice -n 19 lake env lean scripts/PrintAxioms.lean >> $LOG 2>&1
  echo "=== axioms exit=$? $(date -u +%FT%TZ)" >> $LOG
fi
