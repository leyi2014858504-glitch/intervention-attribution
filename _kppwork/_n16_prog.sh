#!/bin/bash
# Progress reporter for the n=16 sweep.  Counts recorded rows across workers.
TAG=${1:-r_sat}
W=/mnt/d/atmos/_kppwork
tot=0
for K in $(seq 1 12); do
  f=$W/n16_v_${TAG}_$K.csv
  if [ -f "$f" ]; then
    c=$(awk 'END{print NR-1}' "$f")
    tot=$((tot + c))
  fi
done
alive=$(pgrep -c -f 'sweep_n16.sh' || echo 0)
pct=$(awk -v a=$tot -v b=65536 'BEGIN{printf "%.2f", 100*a/b}')
echo "[$TAG] $(date +%H:%M:%S)  recorded=$tot / 65536  (${pct}%)  workers_alive=$alive"
