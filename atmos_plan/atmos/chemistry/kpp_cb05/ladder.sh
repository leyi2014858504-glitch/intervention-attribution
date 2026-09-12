#!/bin/bash
# G1/G2 convergence ladder for the KPP-3.5 CB05 standalone box.
# Same scenario each run; only RTOL in cb05_Main.f90 changes. ATOL = 1
# (WRF-Chem module_wkppc_constants: rtols=1.E-3, atols=1.).
cd /root/kpp2 || exit 1
[ -f cb05_Main.f90.tpl ] || cp cb05_Main.f90 cb05_Main.f90.tpl
OUT=/mnt/d/atmos/_ladder.txt
: > "$OUT"

for cfg in "L1_rtol1e-3 1.0d-3" "L2_rtol1e-4 1.0d-4" "L3_rtol1e-6 1.0d-6"; do
  set -- $cfg
  name=$1; rtol=$2
  cp cb05_Main.f90.tpl cb05_Main.f90
  sed -i "s/RTOL(i) = 1.0d-4/RTOL(i) = $rtol/" cb05_Main.f90
  sed -i "s/ATOL(i) = 1.0d-3/ATOL(i) = 1.0d+0/" cb05_Main.f90
  make -f Makefile_cb05 > build_$name.log 2>&1 || { echo "$name BUILD_FAIL" >> "$OUT"; continue; }
  t0=$(date +%s)
  timeout 900 ./cb05.exe > run_$name.log 2>&1
  rc=$?
  t1=$(date +%s)
  nf=$(grep -ci -e warn -e error run_$name.log); [ -z "$nf" ] && nf=0
  awk -v name="$name" -v nf="$nf" -v rc="$rc" -v wall="$((t1-t0))" '
    match($0, /o3= *[0-9.Ee+-]+/) {
      s = substr($0, RSTART+3, RLENGTH-3); gsub(/[ ;]/, "", s); v[++n] = s
    }
    END {
      if (n == 0) { printf "%s rc=%s NO_DATA\n", name, rc; exit }
      sum = 0; for (i = 1; i <= n; i++) sum += v[i]
      printf "%s steps=%d fails=%s wall_s=%s meanO3_ppb=%.4f lastO3_ppb=%.4f\n", \
             name, n, nf, wall, sum/n/2.45e10, v[n]/2.45e10
    }
  ' run_$name.log >> "$OUT"
done
echo LADDER_DONE >> "$OUT"
