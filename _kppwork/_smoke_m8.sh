#!/bin/bash
# Smoke test for the dual-metric extraction added to sweep_n16.sh.
# Runs one coalition, prints new (d3mean, mda8) and the original table row;
# d3mean must match the original v_ppb bit-for-bit.
set -u
W=/mnt/d/atmos/_kppwork
cd /root/kpp7n1
b=$(head -1 /root/n16_chunk_r_sat_1.txt)
nm=${b%.nml}
cp /root/nml16/r_sat/$b cb05.nml
rm -f traj.dat run.log
./cb05.exe > run.log 2>&1
echo "mask=$nm"
awk 'NR>1{n++;t[n]=$1;o[n]=$2; if($1>=48){s+=$2;c++}}
END{
  if(c>0) printf "d3mean=%.6f ", s/c; else {print "NaN"; exit}
  if(n<2){printf "mda8=NaN\n"; exit}
  dt=t[2]-t[1]; w=int(8.0/dt+0.5); if(w<1)w=1
  rs=0; best=-1
  for(i=1;i<=n;i++){rs+=o[i]; if(i>w)rs-=o[i-w]
    if(i>=w && t[i]>=48){a=rs/w; if(a>best)best=a}}
  if(best>=0)printf "mda8=%.6f\n",best; else printf "mda8=NaN\n"
}' traj.dat
echo -n "original row: "; grep "^$nm," $W/n16_v_r_sat.csv
# second mask from a different worker chunk for a bit more coverage
b2=$(sed -n '2500p' /root/n16_chunk_r_sat_1.txt)
nm2=${b2%.nml}
cp /root/nml16/r_sat/$b2 cb05.nml
rm -f traj.dat run.log
./cb05.exe > run.log 2>&1
echo "mask=$nm2"
awk 'NR>1{n++;t[n]=$1;o[n]=$2; if($1>=48){s+=$2;c++}}
END{
  if(c>0) printf "d3mean=%.6f ", s/c; else {print "NaN"; exit}
  dt=t[2]-t[1]; w=int(8.0/dt+0.5); if(w<1)w=1
  rs=0; best=-1
  for(i=1;i<=n;i++){rs+=o[i]; if(i>w)rs-=o[i-w]
    if(i>=w && t[i]>=48){a=rs/w; if(a>best)best=a}}
  if(best>=0)printf "mda8=%.6f\n",best; else printf "mda8=NaN\n"
}' traj.dat
echo -n "original row: "; grep "^$nm2," $W/n16_v_r_sat.csv
