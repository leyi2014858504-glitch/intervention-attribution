#!/bin/bash
cd /root/kpp2 || exit 1
cp /mnt/d/atmos/_kppwork/cb05_Main.f90 .
dos2unix cb05_Main.f90
cp cb05_Main.f90 cb05_Main.f90.tpl
make -f Makefile_cb05 > mk12.log 2>&1
{
  echo "=== BUILD ==="
  grep -c Error mk12.log
  timeout 600 ./cb05.exe > run12.log 2>&1
  echo "=== RUN ==="
  echo RUN_RC=$?
  echo -n 'true_failures(IERR<0)='
  grep -c WARN run12.log
  echo -n 'success_notes(IERR=1 not counted)='
  tail -1 run12.log
  echo "=== G3 ==="
  python3 g3_check.py
} > /mnt/d/atmos/_final.txt 2>&1
