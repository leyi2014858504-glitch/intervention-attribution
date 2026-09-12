#!/bin/bash
# Sweep the 12-player (NOx/VOC/CO x 4 windows) coalition space through
# /root/kpp6.  Pure bash so it survives Windows-side process reaping;
# resumable (skips rows already present in the output csv).
set -u
D=/root/kpp6
N=/mnt/d/atmos/_kppwork/nmldir_co
O=/mnt/d/atmos/_kppwork/co_v.csv
P=/mnt/d/atmos/_kppwork/co_progress.txt
[ -f $O ] || echo "name,v_ppb,nfail,rc,nrow" > $O
done_list=$(mktemp)
grep -o '^[a-z0-9_]*' $O | sort -u > $done_list 2>/dev/null || true
n=0
for f in $N/*.nml; do
  b=$(basename $f .nml)
  grep -qx "$b" $done_list && continue
  cp $f $D/cb05.nml
  cd $D
  rm -f traj.dat run.log
  timeout 120 ./cb05.exe > run.log 2>&1
  rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
  nrow=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
  v=$(awk 'NR>1 && $1>=48 {s+=$2; n++} END {if(n>0) printf "%.6f", s/n; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "$b,${v:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  n=$((n+1))
  [ $((n % 512)) -eq 0 ] && echo "completed $n at $b ($(date +%H:%M:%S))" >> $P
done
echo "SWEEP_CO_DONE $(date)" >> $P
rm -f $done_list
