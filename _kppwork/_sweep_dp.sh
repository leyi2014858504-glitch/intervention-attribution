#!/bin/bash
# Deer Park Taylor-2 campaign: 512 derivative + 150 truth runs, 12 workers.
#   $1 = runall | worker
#   $2 = worker index (worker mode only)
set -u
W=/mnt/d/atmos/_kppwork
N=/root/nmlhddm/dp
NW=12

if [ "${1:-}" = "runall" ]; then
  ls $N > /root/dp_files.txt
  echo "[dp] tasks: $(wc -l < /root/dp_files.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/dp_files.txt > /root/dp_chunk_$K.txt
    [ -d /root/kpp7n$K ] || cp -r /root/kpp7 /root/kpp7n$K
  done
  rm -f $W/dp_v_*.csv /root/dp_done_*.txt
  for K in $(seq 1 $NW); do
    bash "$0" worker $K > $W/dp_worker_$K.out 2>&1 &
  done
  wait
  n=$(cat $W/dp_v_*.csv 2>/dev/null | grep -c '^[a-z]')
  bad=$(cat $W/dp_v_*.csv 2>/dev/null | awk -F, '$2=="NaN"||$3=="NaN"' | wc -l)
  echo "[dp] ALL WORKERS EXITED $(date) -- rows $n, NaN $bad"
  exit 0
fi

if [ "${1:-}" = "worker" ]; then
  K=$2
  D=/root/kpp7n$K
  O=$W/dp_v_$K.csv
  [ -d $D ] || cp -r /root/kpp7 $D
  [ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc" > $O
  DONE=/root/dp_done_$K.txt
  grep -o '^[a-z0-9_]*' $O | sort -u > $DONE 2>/dev/null || true
  while read -r b; do
    nm=${b%.nml}
    grep -qx "$nm" $DONE && continue
    cp $N/$b $D/cb05.nml
    cd $D
    rm -f traj.dat run.log
    timeout 120 ./cb05.exe > run.log 2>&1
    rc=$?
    nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
    read -r v v8 <<EOF2
$(awk 'NR>1{n++;t[n]=$1;o[n]=$2; if($1>=48){s+=$2;c++}}
END{
  if(c>0) printf "%.6f ", s/c; else printf "NaN ";
  if(n<2){printf "NaN"; exit}
  dt=t[2]-t[1]; w=int(8.0/dt+0.5); if(w<1)w=1;
  rs=0; best=-1;
  for(i=1;i<=n;i++){rs+=o[i]; if(i>w)rs-=o[i-w];
    if(i>=w && t[i]>=48){a=rs/w; if(a>best)best=a}}
  if(best>=0)printf "%.6f",best; else printf "NaN"
}' traj.dat 2>/dev/null)
EOF2
    echo "$nm,${v:-NaN},${v8:-NaN},${nf:-NA},$rc" >> $O
    echo "$nm" >> $DONE
  done < /root/dp_chunk_$K.txt
  echo "$K: DONE $(date)" >> $W/dp_prog_$K.txt
  exit 0
fi

echo "usage: _sweep_dp.sh runall | worker K"
