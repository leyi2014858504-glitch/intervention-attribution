#!/bin/bash
# Rebuild co_v.csv from the clean co_v_m8 shards (v48 column), preserving the
# original header.  The pre-2026-09-16 file had 18 corrupt rows (interleaved
# appends after a WSL crash); the m8 re-run reproduced the other 4094 rows
# bit-for-bit and re-judge confirmed the new values are correct.
set -u
W=/mnt/d/atmos/_kppwork
cp $W/co_v.csv $W/co_v_corrupt_backup.csv
echo "name,v_ppb,nfail,rc,nrow" > $W/co_v.csv
for K in $(seq 1 6); do
  awk -F, 'NR>1{printf "%s,%s,%s,%s,%s\n", $1,$2,$4,$5,$6}' $W/co_v_m8_$K.csv >> $W/co_v.csv
done
n=$(awk 'END{print NR-1}' $W/co_v.csv)
u=$(awk -F, 'NR>1{print $1}' $W/co_v.csv | sort -u | wc -l)
echo "co_v.csv rebuilt: $n rows, $u distinct (expect 8192 / 8192)"
