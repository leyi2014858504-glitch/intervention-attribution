#!/bin/bash
# Refresh the 12-player backend from a clean copy, smoke-test it, then run the
# isoprene probe.  Guards against the transient "empty run.log / no traj.dat"
# state that /root/kpp6 fell into after the parallel sweep touched it.
set -u
W=/mnt/d/atmos/_kppwork
rm -rf /root/kpp6
cp -r /root/kpp5 /root/kpp6
cp $W/cb05d_Main.f90 /root/kpp6/cb05_Main.f90
touch /root/kpp6/cb05_Main.f90
cd /root/kpp6
make -f Makefile_cb05 > mk7.log 2>&1 || true
grep -qi error mk7.log && { echo "REBUILD FAILED"; grep -i -B2 -A4 error mk7.log | head -30; exit 1; }

echo "=== smoke on the refreshed /root/kpp6 ==="
cat > cb05.nml <<'NML'
&RUNCFG
  TEND_S = 259200.0,
  DT_S = 300.0,
  EMIS_NOX = 16.0, EMIS_VOC = 12.0, EMIS_ISO = 1.0, EMIS_CO = 40.0,
  NOX_FACTOR = 1.0, VOC_FACTOR = 1.0,
  NOXW = 1.0,1.0,1.0,1.0, VOCW = 1.0,1.0,1.0,1.0, COW = 1.0,1.0,1.0,1.0,
  WKEND = .FALSE., DUMMYF = .TRUE., CONSTJ = .FALSE., C_H2O = 5.0e18,
  KDIL_DAY = 1500.0, KDIL_NIGHT = 86400.0,
  BG_O3 = 30.0, BG_NO2 = 5.0, BG_NO = 1.0, BG_CO = 120.0,
  BG_VOC = 20.0, BG_ISO = 0.2,
  TMIN_K = 288.0, TMAX_K = 301.0,
/
NML
rm -f traj.dat run.log
timeout 120 ./cb05.exe > run.log 2>&1
rc=$?
n=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
echo "smoke rc=$rc traj=$n"
[ "$rc" = "0" ] && [ "$n" = "864" ] || { echo "SMOKE FAILED"; cat run.log; exit 1; }
v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END {printf "%.6f", s/c}' traj.dat)
echo "smoke v48=$v  (expect 21.182293)"

echo ""
echo "=== probe ==="
cd $W && python3 probe_isop.py
