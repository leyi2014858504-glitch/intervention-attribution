#!/bin/bash
# Generate + build the DUMMY-REACTION backend in /root/kpp5.
#
#   /root/kpp4 : scenario driver with operator splitting (cb05s_Main.f90)
#   /root/kpp5 : scenario driver with dummy-reaction forcing (cb05d_Main.f90)
# The two exist side by side for the A/B comparison.
set -e
export KPP_HOME=/opt/KPP        # without this KPP cannot find int/rosenbrock.def
W=/mnt/d/atmos/_kppwork
R=/root/kpp5

rm -rf $R
mkdir -p $R
cd $R

echo "=== stage mechanism inputs ==="
for f in cb05.kpp cb05_sorg_aq.def cb05_sorg_aq.eqn cb05_sorg_aq.spc atoms_red; do
  cp $W/$f .
  dos2unix -q $f
done
# the archived .def still defines EP2/EP3, which clash with KPP's own generic
# interfaces in the generated Rates module -> rename to EP2X/EP3X (see ren2.sed)
cp $W/ren2.sed .
dos2unix -q ren2.sed
sed -i -f ren2.sed cb05_sorg_aq.def
grep -c "FUNCTION EP2X\|FUNCTION EP3X" cb05_sorg_aq.def

echo "=== run KPP ==="
/opt/KPP/bin/kpp cb05.kpp > kpp.log 2>&1 || true
# KPP exits 0 even on fatal errors -> judge by the success banner
if ! grep -q "successfully created" kpp.log; then
  echo "KPP_FAILED"
  tail -30 kpp.log
  exit 1
fi
echo "KPP_OK"
grep -E "NSPEC|NVAR|NFIX|NREACT" cb05_Parameters.f90 | head -4

echo "=== compile ==="
cp /root/kpp2/Makefile_cb05 .
cp $W/cb05d_Main.f90 cb05_Main.f90
dos2unix -q cb05_Main.f90
make -f Makefile_cb05 > mk.log 2>&1 || true
if grep -qi "error" mk.log; then
  echo "MAKE_FAILED"
  grep -i -B2 -A4 "error" mk.log | head -40
  exit 1
fi
echo "MAKE_OK"

{
  echo "=== smoke A: dummy forcing ON, weekday, 72h, dt=300s ==="
  rm -f cb05.nml
  t0=$(date +%s)
  timeout 300 ./cb05.exe > run_smoke.log 2>&1; echo RUN_RC=$?
  echo wall_s=$(( $(date +%s) - t0 ))
  grep -c WARN run_smoke.log || true
  tail -2 run_smoke.log
  awk 'NR>1 { d=int(($1)/24); n[d]++; s[d]+=$2; if($2>mx[d]) mx[d]=$2 }
       END { for (k=0;k<=2;k++) printf "day%d meanO3=%.2f maxO3=%.2f\n", k+1, s[k]/n[k], mx[k] }' traj.dat

  echo ""
  echo "=== smoke B: closed box (DUMMYF=F, CONSTJ=T), 24h, dt=1800s ==="
  cat > cb05.nml <<'NML'
&RUNCFG
  TEND_S = 86400.0,
  DT_S = 1800.0,
  DUMMYF = .FALSE.,
  CONSTJ = .TRUE.,
  C_H2O = 0.0,
/
NML
  timeout 300 ./cb05.exe > run_closed.log 2>&1; echo RUN_RC=$?
  grep -c WARN run_closed.log || true
  tail -2 run_closed.log
  cp state_out.dat closed_out.dat
  grep -E "^(o3|no2|oh|H2O|M|FSRC|FBCK|FSNK) " closed_out.dat
} > /mnt/d/atmos/_smoke5.txt 2>&1

echo SMOKE_WRITTEN
