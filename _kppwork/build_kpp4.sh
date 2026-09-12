#!/bin/bash
# Build the scenario backend (/root/kpp4) from the validated pure-chem dir.
set -e
rm -rf /root/kpp4
cp -r /root/kpp2 /root/kpp4
cd /root/kpp4
cp /mnt/d/atmos/_kppwork/cb05s_Main.f90 cb05_Main.f90
dos2unix cb05_Main.f90
rm -f cb05.nml
make -f Makefile_cb05 > mk1.log 2>&1
echo "MAKE_RC=$?"
{
  echo "=== smoke: default weekday 72h, dt=300s ==="
  t0=$(date +%s)
  timeout 300 ./cb05.exe > run_smoke.log 2>&1
  echo RUN_RC=$?
  echo wall_s=$(( $(date +%s) - t0 ))
  grep -c WARN run_smoke.log || true
  tail -3 run_smoke.log
  echo "--- traj: day stats (o3 mean/max per 24h) ---"
  awk 'NR>1 { d=int(($1)/24); n[d]++; s[d]+=$2; if($2>mx[d]) mx[d]=$2 }
       END { for (k=0;k<=2;k++) printf "day%d meanO3=%.2f maxO3=%.2f\n", k+1, s[k]/n[k], mx[k] }' traj.dat
  wc -l traj.dat
} > /mnt/d/atmos/_smoke.txt 2>&1
echo SMOKE_WRITTEN
