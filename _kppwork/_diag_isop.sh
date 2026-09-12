#!/bin/bash
# Diagnose why the ISOP probe run produced an empty run.log.
set -u
cd /root/kpp6
ls -la cb05.exe || echo "NO EXE"
cp /mnt/d/atmos/_kppwork/cb05.nml ./cb05.nml
echo "--- namelist written by kpp_box ---"
cat cb05.nml
rm -f traj.dat run.log
timeout 60 ./cb05.exe > run.log 2>&1
echo "rc=$?"
echo "--- run.log ---"
cat run.log
echo "--- traj rows: $(awk 'END{print NR-1}' traj.dat 2>/dev/null) ---"
