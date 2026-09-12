#!/bin/bash
# One parallel worker of the 12-player sweep.  $1 = worker index.
# RESUME-SAFE: appends to its own csv and skips coalitions already recorded,
# so a reaped worker can simply be relaunched.
set -u
K=$1
W=/mnt/d/atmos/_kppwork
N=$W/nmldir_co
D=/root/kpp6w$K
O=$W/co_v_$K.csv
P=$W/co_prog_$K.txt
[ -f $O ] || echo "name,v_ppb,nfail,rc,nrow" > $O
DONE=$W/done_$K.txt
grep -o '^[a-z0-9_]*' $O | sort -u > $DONE 2>/dev/null || true
n=0
while read -r b; do
  nm=${b%.nml}
  grep -qx "$nm" $DONE && continue
  cp $N/$b $D/cb05.nml
  cd $D
  rm -f traj.dat run.log
  timeout 120 ./cb05.exe > run.log 2>&1
  rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
  nrow=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
  v=$(awk 'NR>1 && $1>=48 {s+=$2; n++} END {if(n>0) printf "%.6f", s/n; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "$nm,${v:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  echo "$nm" >> $DONE
  n=$((n+1))
  [ $((n % 256)) -eq 0 ] && echo "$K: +$n at $b $(date +%H:%M:%S)" >> $P
done < $W/chunk_$K.txt
echo "$K: DONE $(date)" >> $P
