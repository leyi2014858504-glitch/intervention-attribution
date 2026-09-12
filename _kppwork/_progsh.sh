#!/bin/bash
W=/mnt/d/atmos/_kppwork
date
tot=0
for K in $(seq 1 12); do
  n=$(wc -l < $W/sham_v_$K.csv 2>/dev/null || echo 0)
  printf "w%-3s %5s   " "$K" "$n"
  [ $((K % 4)) -eq 0 ] && echo ""
  tot=$((tot + n))
done
echo ""
echo "TOTAL rows = $tot  / 8192 expected"
echo "workers alive: $(pgrep -c -f 'sweep_sham.sh' || echo 0)"
echo "cb05 running:  $(pgrep -c cb05.exe || echo 0)"
cat $W/sham_worker_*.out 2>/dev/null | grep -v '^$' | tail -5
