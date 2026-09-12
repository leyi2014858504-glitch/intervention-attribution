#!/bin/bash
# V2: re-run the G3 conservation gate on the DUMMY-FORCING mechanism.
#
# Uses the unchanged pure-chemistry driver (the G1/G2/G3 driver of the old
# build).  It never touches the forcing globals, so every dummy rate is zero
# and the box is closed by construction.
set -e
W=/mnt/d/atmos/_kppwork
R=/root/kpp5g3

rm -rf $R
cp -r /root/kpp5 $R
cd $R
cp $W/cb05_Main.f90 cb05_Main.f90
dos2unix -q cb05_Main.f90
make -f Makefile_cb05 > mk_g3.log 2>&1 || true
if grep -qi "error" mk_g3.log; then echo "MAKE_FAILED"; grep -i -B2 -A4 error mk_g3.log | head -30; exit 1; fi

timeout 600 ./cb05.exe > run_g3.log 2>&1 || true
{
  echo "=== V2: G3 on the dummy-forcing mechanism ==="
  echo -n "solver failures (IERR<0): "
  grep -c WARN run_g3.log || true
  tail -1 run_g3.log
  echo "--- conservation (chemistry reactions only) ---"
  python3 $W/g3_check.py $R 191
  echo "--- O3 at t_end vs the old-build reference ---"
  grep "^o3 " cb05_tend.dat
  echo "old-build reference:     0.576169439164043E+13"
  cp cb05_t0.dat cb05_tend.dat cb05_Function.f90 /mnt/d/atmos/_kppwork/src/ 2>/dev/null || true
  mv /mnt/d/atmos/_kppwork/src/cb05_t0.dat /mnt/d/atmos/_kppwork/src/new_t0.dat
  mv /mnt/d/atmos/_kppwork/src/cb05_tend.dat /mnt/d/atmos/_kppwork/src/new_tend.dat
} > /mnt/d/atmos/_v2.txt 2>&1
echo V2_WRITTEN
