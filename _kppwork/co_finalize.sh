#!/bin/bash
# Final merge (dedup by coalition name) + 12-player analysis.
set -u
W=/mnt/d/atmos/_kppwork
head -1 $W/co_v_1.csv > $W/co_v.csv
for K in 1 2 3 4 5 6; do tail -n +2 $W/co_v_$K.csv; done \
  | awk -F, '!seen[$1]++' >> $W/co_v.csv
{
  echo "merged rows: $(wc -l < $W/co_v.csv)"
  echo "rows per scenario:"
  awk -F, 'NR>1 {split($1,a,"_m"); print a[1]}' $W/co_v.csv | sort | uniq -c
  awk -F, 'NR>1 && ($3!="0" || $4!="0" || $5!="864") {bad++} END {print "BAD_ROWS="bad+0}' $W/co_v.csv
}
echo "=== 12-player analysis ==="
cd $W && python3 exp_co12.py
