#!/bin/bash
# 1) generate the 2 scenarios x 65536 namelists
# 2) SMOKE TEST: the all-players-ON config of the n=16 ladder must reproduce the
#    frozen V-gate value 21.182293.  If it does not, the ladder is not the same
#    model family and the whole comparison would be void.
set -u
W=/mnt/d/atmos/_kppwork
echo "=== generating ==="
python3 $W/gen_n16.py
echo "=== sizes ==="
du -sh /root/nml16/r_sat /root/nml16/r_lim1 2>/dev/null
echo "=== V-gate smoke test on /root/kpp7 ==="
cp /root/nml16/vgate.nml /root/kpp7/cb05.nml
cd /root/kpp7
rm -f traj.dat run.log
timeout 200 ./cb05.exe > run.log 2>&1
echo "rc=$?"
v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END{if(c>0) printf "%.6f", s/c; else printf "NaN"}' traj.dat)
echo "v48 = $v    (frozen reference 21.182293)"
awk 'NR==1{print "header:", $0}' traj.dat
