#!/bin/bash
# Progress probe for the parallel 12-player sweep.
date
tot=0
for K in 1 2 3 4 5 6; do
  n=$(wc -l < /mnt/d/atmos/_kppwork/co_v_$K.csv 2>/dev/null || echo 0)
  echo "worker $K: $n rows"
  tot=$((tot + n))
done
echo "TOTAL rows (incl. 6 headers) = $tot"
echo "live cb05.exe processes: $(pgrep -c cb05.exe || echo 0)"
cat /mnt/d/atmos/_kppwork/co_prog_*.txt 2>/dev/null | tail -8
cat /mnt/d/atmos/_kppwork/co_done.txt 2>/dev/null
