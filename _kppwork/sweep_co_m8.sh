#!/bin/bash
# co12 (3 sectors x 4 windows) MDA8 campaign: dual-metric full enumeration,
# 2 scenarios x 4096 coalitions.  Run AFTER the n16 m8 campaign (CPU).
# Resumable: workers append to co_v_m8_$K.csv and skip recorded names.
set -u
W=/mnt/d/atmos/_kppwork
N=$W/nmldir_co
NW=6

_prep () {
  ls $N > $W/co_files_m8.txt
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' $W/co_files_m8.txt > $W/co_chunk_m8_$K.txt
    [ -d /root/kpp6w$K ] || cp -r /root/kpp6 /root/kpp6w$K
  done
}

if [ "${1:-}" = "runall" ]; then
  _prep
  for K in $(seq 1 $NW); do
    bash "$0" $K > $W/co_worker_m8_$K.out 2>&1 &
  done
  wait
  echo "[co-m8] ALL WORKERS EXITED $(date)"
  exit 0
fi

K=$1
D=/root/kpp6w$K
O=$W/co_v_m8_$K.csv
P=$W/co_prog_m8_$K.txt
[ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc,nrow" > $O
DONE=/root/co_done_m8_$K.txt
grep -o '^[a-z0-9_]*' $O | sort -u > $DONE 2>/dev/null || true
n=0
while read -r b; do
  nm=${b%.nml}
  grep -qx "$nm" $DONE && continue
  cp $N/$b $D/cb05.nml
  cd $D
  rm -f traj.dat run.log
  timeout 120 ./cb05.exe > run.log 2>&1
  rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
  nrow=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
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
  echo "$nm,${v:-NaN},${v8:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  echo "$nm" >> $DONE
  n=$((n+1))
  [ $((n % 256)) -eq 0 ] && echo "$K: +$n $(date +%H:%M:%S)" >> $P
done < $W/co_chunk_m8_$K.txt
echo "$K: DONE $(date)" >> $P
