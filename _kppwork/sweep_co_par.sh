#!/bin/bash
# Parallel driver for the 12-player (NOx/VOC/CO x 4 windows) coalition sweep.
# Splits the 8192 namelists over NW workers, each with its own build dir, then
# merges the per-worker CSVs.  Run this FOREGROUND (a detached parent lets the
# WSL VM idle out and the sweep dies).
set -u
W=/mnt/d/atmos/_kppwork
N=$W/nmldir_co
NW=6

rm -f $W/co_v.csv $W/co_v_*.csv $W/co_prog_*.txt $W/co_done.txt

ls $N > $W/co_files.txt
for K in $(seq 1 $NW); do
  awk -v k=$K -v n=$NW 'NR % n == (k % n)' $W/co_files.txt > $W/chunk_$K.txt
done

for K in $(seq 1 $NW); do
  rm -rf /root/kpp6w$K
  cp -r /root/kpp6 /root/kpp6w$K
done

for K in $(seq 1 $NW); do
  bash $W/sweep_co_w.sh $K &
done
wait

cat $W/co_v_1.csv > $W/co_v.csv
for K in $(seq 2 $NW); do tail -n +2 $W/co_v_$K.csv >> $W/co_v.csv; done
{
  echo "PAR_DONE $(date)"
  wc -l $W/co_v.csv
  awk -F, 'NR>1 && ($3!="0" || $4!="0" || $5!="864") {bad++} END {print "BAD_ROWS="bad+0}' $W/co_v.csv
} > $W/co_done.txt
