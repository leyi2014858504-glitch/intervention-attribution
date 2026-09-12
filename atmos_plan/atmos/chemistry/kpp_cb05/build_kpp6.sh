#!/bin/bash
# Build /root/kpp6 = the dummy-reaction backend PLUS the third intervention
# sector (COW), i.e. the 12-player game.  The mechanism is untouched, so this
# copies the frozen /root/kpp5 tree and only recompiles the driver.
#
# BACKWARD-COMPAT GATE: the namelist below does NOT set COW, so it must
# reproduce the frozen 8-player values bit for bit:
#     r_sat  mask 0 = 21.182293 ppb      (all windows on)
#     r_sat  mask 3 = 30.240258 ppb      (NOx W1+W2 off)
set -e
W=/mnt/d/atmos/_kppwork
R=/root/kpp6

rm -rf $R
cp -r /root/kpp5 $R
cd $R
rm -f cb05.nml traj.dat state_out.dat run.log traj.dat.keep

echo "=== driver swap + recompile (mechanism unchanged) ==="
cp $W/cb05d_Main.f90 cb05_Main.f90
dos2unix -q cb05_Main.f90
touch cb05_Main.f90
make -f Makefile_cb05 > mk6.log 2>&1 || true
if grep -qi "error" mk6.log; then
  echo MAKE_FAILED; grep -i -B2 -A4 "error" mk6.log | head -40; exit 1
fi
echo MAKE_OK

run_one () {   # $1 = label, $2 = NOXW, $3 = VOCW  (COW line appended by caller)
  cat > cb05.nml <<NML
&RUNCFG
  TEND_S = 259200.0,
  DT_S = 300.0,
  EMIS_NOX = 16.0, EMIS_VOC = 12.0, EMIS_ISO = 1.0, EMIS_CO = 40.0,
  NOX_FACTOR = 1.0, VOC_FACTOR = 1.0,
  NOXW = $2, VOCW = $3,
$4
  WKEND = .FALSE., DUMMYF = .TRUE., CONSTJ = .FALSE., C_H2O = 5.0e18,
  KDIL_DAY = 1500.0, KDIL_NIGHT = 86400.0,
  BG_O3 = 30.0, BG_NO2 = 5.0, BG_NO = 1.0, BG_CO = 120.0,
  BG_VOC = 20.0, BG_ISO = 0.2,
  TMIN_K = 288.0, TMAX_K = 301.0,
/
NML
  timeout 300 ./cb05.exe > run_v6.log 2>&1; rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run_v6.log | tr -dc '0-9' || true)
  v=$(awk 'NR>1 && $1>=48 {s+=$2; n++} END {if(n>0) printf "%.6f", s/n; else printf "NaN"}' traj.dat 2>/dev/null || printf "NaN")
  echo "  $1  rc=$rc nfail=${nf:-NA}  v48=$v"
  if [ "$rc" != "0" ] || [ "$v" = "NaN" ]; then
    echo "  ---- run_v6.log ----"; tail -6 run_v6.log
    echo "  ---- cb05.nml ----"; cat cb05.nml
  fi
}

{
  echo "=== backward-compat gate (COW omitted -> defaults all ON) ==="
  run_one "r_sat mask0 (expect 21.182293)" "1.0,1.0,1.0,1.0" "1.0,1.0,1.0,1.0" ""
  run_one "r_sat mask3 (expect 30.240258)" "0.0,0.0,1.0,1.0" "1.0,1.0,1.0,1.0" ""

  echo ""
  echo "=== COW switch actually does something ==="
  run_one "CO all off"  "1.0,1.0,1.0,1.0" "1.0,1.0,1.0,1.0" "  COW = 0.0,0.0,0.0,0.0,"
  run_one "CO W3 off"   "1.0,1.0,1.0,1.0" "1.0,1.0,1.0,1.0" "  COW = 1.0,1.0,0.0,1.0,"
} > /mnt/d/atmos/_build6.txt 2>&1

echo BUILD6_DONE
