#!/bin/bash
echo "WSL_OK"
W=/mnt/d/atmos/_kppwork
for TAG in r_sat r_lim1; do
  tot=0
  for K in $(seq 1 12); do
    f=$W/n16_v_${TAG}_$K.csv
    if [ -f "$f" ]; then
      c=$(awk 'END{print NR-1}' "$f")
      tot=$((tot + c))
    fi
  done
  echo "$TAG: $tot / 65536 rows on disk"
done
echo "--- r_sat anchors (must match the frozen V-gate) ---"
grep -h '^m00000,' $W/n16_v_r_sat_*.csv
grep -h '^m00012,' $W/n16_v_r_sat_*.csv
echo "--- n=8 / n=12 oracles still present ---"
awk 'END{print "scen_v.csv rows:", NR-1}' $W/scen_v.csv
awk 'END{print "co_v.csv rows:", NR-1}' $W/co_v.csv
