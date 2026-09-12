#!/bin/bash
cd /root/kpp6
cp /root/isotest/cb05.nml ./cb05.nml
rm -f traj.dat run.log
timeout 60 ./cb05.exe > run.log 2>&1
echo "kpp6_rc=$?"
echo "traj=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)"
echo "logbytes=$(wc -c < run.log)"
cat run.log
# refresh kpp6 from the known-good copy either way
cd /
rm -rf /root/kpp6
cp -r /root/isotest /root/kpp6
echo "kpp6 refreshed from isotest"
