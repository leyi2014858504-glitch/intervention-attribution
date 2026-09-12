#!/bin/bash
W=/mnt/d/atmos/_kppwork
for K in 1 2 3 4 5 6; do
  echo "===== sham_worker_$K.out ====="
  cat $W/sham_worker_$K.out 2>/dev/null | tail -20
  echo "----- sham_prog_$K.txt -----"
  cat $W/sham_prog_$K.txt 2>/dev/null | tail -5
done
echo "===== last rows written ====="
for K in 1 2 3 4 5 6; do echo "w$K: $(tail -1 $W/sham_v_$K.csv 2>/dev/null)"; done
echo "===== any error in the per-run logs? ====="
for K in 1 2 3 4 5 6; do echo "w$K run.log: $(tail -2 /root/kpp6s$K/run.log 2>/dev/null | tr '\n' ' ')"; done
echo "===== worker dirs exist? ====="
for K in 1 2 3 4 5 6; do echo "kpp6s$K: $(ls /root/kpp6s$K/cb05.exe 2>/dev/null | wc -l)"; done
