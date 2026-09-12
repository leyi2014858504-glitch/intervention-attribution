#!/bin/bash
# Merge the 12 sham worker shards (dedup by coalition name) and run the
# exact-equality control checks.
set -u
W=/mnt/d/atmos/_kppwork
first=$(ls $W/sham_v_1.csv 2>/dev/null)
head -1 $first > $W/sham_v.csv
for K in $(seq 1 12); do tail -n +2 $W/sham_v_$K.csv 2>/dev/null; done \
  | awk -F, '!seen[$1]++' >> $W/sham_v.csv
{
  echo "merged rows: $(wc -l < $W/sham_v.csv)  (expect 8193)"
  echo "rows per scenario:"
  awk -F, 'NR>1 {split($1,a,"_m"); print a[1]}' $W/sham_v.csv | sort | uniq -c
  awk -F, 'NR>1 && ($3!="0" || $4!="0" || $5!="864") {bad++} END {print "BAD_ROWS="bad+0}' $W/sham_v.csv
}
echo "=== sham control checks ==="
cd $W && python3 exp_sham.py
