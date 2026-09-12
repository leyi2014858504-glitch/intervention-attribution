#!/bin/bash
# Isolate: is it the directory (/root/kpp6) or the kpp_box-written namelist?
echo "=== A) working dir /root/twice + kpp_box namelist ==="
cd /root/twice
cp /mnt/d/atmos/_kppwork/cb05.nml ./cb05.nml
rm -f traj.dat run.log state_out.dat
timeout 60 ./cb05.exe > run.log 2>&1
echo "rc=$? traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null) logbytes=$(wc -c < run.log 2>/dev/null)"
cat run.log

echo ""
echo "=== B) /root/kpp6 + isotest namelist ==="
cd /root/kpp6
cp /root/isotest/cb05.nml ./cb05.nml
rm -f traj.dat run.log state_out.dat
timeout 60 ./cb05.exe > run.log 2>&1
echo "rc=$? traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null) logbytes=$(wc -c < run.log 2>/dev/null)"
cat run.log

echo ""
echo "=== C) file inventory of /root/kpp6 ==="
ls -la /root/kpp6 | head -12
echo "count=$(ls /root/kpp6 | wc -l)"

echo ""
echo "=== D) kpp_box namelist content ==="
cat /mnt/d/atmos/_kppwork/cb05.nml
