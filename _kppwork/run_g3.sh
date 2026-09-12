#!/bin/bash
cd /root/kpp2 || exit 1
cp cb05_Main.f90.tpl cb05_Main.f90
sed -i 's/ATOL(i) = 1.0d-3/ATOL(i) = 1.0d+0/' cb05_Main.f90
make -f Makefile_cb05 > /dev/null 2>&1
timeout 600 ./cb05.exe > run_g3.log 2>&1
{
  echo RUN_RC=$?
  echo -n 'fails='; grep -ci -e warn -e error run_g3.log
  echo -n 'dat_lines='; wc -l < cb05.dat
  echo '--- dat head ---'; head -2 cb05.dat | cut -c1-160
  echo '--- dat tail ---'; tail -1 cb05.dat | cut -c1-160
  echo '--- stoich files ---'; ls cb05_Stoichiom*
} > /mnt/d/atmos/_g3b.txt 2>&1
