#!/bin/bash
W=/mnt/d/atmos/_kppwork
date
tot=0
for K in 1 2 3 4 5 6; do
  n=$(wc -l < $W/scale_v_$K.csv 2>/dev/null || echo 0)
  echo "worker $K: $n rows"
  tot=$((tot + n))
done
echo "TOTAL rows (incl. 6 headers) = $tot  / 6970 expected"
echo "workers alive: $(pgrep -c -f 'sweep_scale.sh' || echo 0)"
echo "cb05 running: $(pgrep -c cb05.exe || echo 0)"
tail -4 $W/scale_prog_*.txt 2>/dev/null
[ -f $W/scale_done.txt ] && echo "DONE MARKER PRESENT"
