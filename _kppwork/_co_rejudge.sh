#!/bin/bash
# Adjudicate the 18 disagreeing co12 masks: re-run each fresh in kpp6w1 and
# print old / new / rejudge side by side.
set -u
W=/mnt/d/atmos/_kppwork
D=/root/kpp6w1
cd $D
for M in "r_sat 191" "r_sat 294" "r_sat 705" "r_sat 709" "r_sat 711" "r_lim1 3956" "r_lim1 4044"; do
  set -- $M
  scen=$1; m=$2
  nm=$(printf "%s_m%04d" $scen $m)
  cp $W/nmldir_co/$nm.nml cb05.nml
  rm -f traj.dat run.log
  ./cb05.exe > run.log 2>&1
  v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END{printf "%.6f", s/c}' traj.dat)
  old=$(grep "^$nm," $W/co_v.csv | cut -d, -f2)
  new=$(grep -h "^$nm," $W/co_v_m8_*.csv | cut -d, -f2)
  echo "$nm  old=$old  new=$new  rejudge=$v"
done
