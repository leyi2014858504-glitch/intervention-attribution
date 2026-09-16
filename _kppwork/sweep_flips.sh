#!/bin/bash
# Scan A worker loop: 1825 dual-metric box runs from /root/flipscan, 12 workers
# on the existing /root/kpp7n$K builds (same pattern as sweep_n16.sh).
#   $1 = "runall"  -> prep + foreground + wait
#   $1 = N         -> worker N
# Resumable: workers append to flip_v_$K.csv and skip recorded names.
set -u
W=/mnt/d/atmos/_kppwork
N=/root/flipscan
NW=12

_prep () {
  ls $N > /root/flip_files.txt
  echo "[flip] total tasks: $(wc -l < /root/flip_files.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/flip_files.txt > /root/flip_chunk_$K.txt
    [ -d /root/kpp7n$K ] || cp -r /root/kpp7 /root/kpp7n$K
  done
}

if [ "${1:-}" = "runall" ]; then
  _prep
  for K in $(seq 1 $NW); do
    bash "$0" $K > $W/flip_worker_$K.out 2>&1 &
  done
  wait
  echo "[flip] ALL WORKERS EXITED $(date)"
  exit 0
fi

K=$1
D=/root/kpp7n$K
O=$W/flip_v_$K.csv
[ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc,nrow" > $O
DONE=/root/flip_done_$K.txt
grep -o '^[a-z0-9]*' $O | sort -u > $DONE 2>/dev/null || true
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
done < /root/flip_chunk_$K.txt
echo "$K: DONE $(date)" >> $W/flip_prog_$K.txt
