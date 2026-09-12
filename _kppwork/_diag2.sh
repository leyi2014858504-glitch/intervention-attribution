#!/bin/bash
# Deeper diagnosis: fresh dir, explicit rc, file listing, and a control run
# with the known-good backward-compat namelist.
set -u
rm -rf /root/isotest
cp -r /root/kpp6 /root/isotest
cd /root/isotest

echo "=== control: EMIS_ISO = 1.0 (known good) ==="
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
./cb05.exe > run.log 2>&1
echo "rc=$?  traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)  logbytes=$(wc -c < run.log)"

echo ""
echo "=== suspect: EMIS_ISO = 0.0 ==="
sed -i 's/EMIS_ISO = 1.0/EMIS_ISO = 0.0/' cb05.nml
rm -f traj.dat run.log
./cb05.exe > run.log 2>&1
echo "rc=$?  traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)  logbytes=$(wc -c < run.log)"
cat run.log

echo ""
echo "=== control2: EMIS_ISO = 3.0 ==="
sed -i 's/EMIS_ISO = 0.0/EMIS_ISO = 3.0/' cb05.nml
rm -f traj.dat run.log
./cb05.exe > run.log 2>&1
echo "rc=$?  traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)  logbytes=$(wc -c < run.log)"
