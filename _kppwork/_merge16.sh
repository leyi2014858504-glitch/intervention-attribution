#!/bin/bash
# Merge the 24 per-worker n=16 shards into one truth table per scenario, the way
# scale_v.csv is kept.  The shards themselves stay ignored (they are per-worker
# partials); the merged table is what gets versioned.
W=/mnt/d/atmos/_kppwork
for TAG in r_sat r_lim1; do
  OUT=$W/n16_v_${TAG}.csv
  echo "name,v_ppb,nfail,rc,nrow" > $OUT
  for K in $(seq 1 12); do
    f=$W/n16_v_${TAG}_$K.csv
    [ -f "$f" ] || continue
    awk 'NR>1' "$f" >> $OUT
  done
  n=$(awk 'END{print NR-1}' $OUT)
  u=$(awk -F, 'NR>1{print $1}' $OUT | sort -u | wc -l)
  echo "$TAG: $n rows, $u distinct masks (expect 65536 / 65536)"
done
ls -l $W/n16_v_r_sat.csv $W/n16_v_r_lim1.csv
