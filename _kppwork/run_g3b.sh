#!/bin/bash
cd /root/kpp2 || exit 1
cp /mnt/d/atmos/_kppwork/cb05_Main.f90 .
dos2unix cb05_Main.f90
cp cb05_Main.f90 cb05_Main.f90.tpl
apt-get install -y python3-numpy > /dev/null 2>&1
make -f Makefile_cb05 > mk11.log 2>&1
{
  echo MAKE_RC=$?
  grep -c Error mk11.log
  timeout 600 ./cb05.exe > run11.log 2>&1
  echo RUN_RC=$?
  echo -n 'fails='; grep -c WARN run11.log
  echo -n 't0_lines='; wc -l < cb05_t0.dat
  echo -n 'tend_lines='; wc -l < cb05_tend.dat
  head -3 cb05_tend.dat
  cp cb05_t0.dat cb05_tend.dat cb05_Function.f90 /mnt/d/atmos/_kppwork/src/
} > /mnt/d/atmos/_g3d.txt 2>&1
