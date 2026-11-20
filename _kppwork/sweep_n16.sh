#!/bin/bash
# n=16 sweep: 2 sectors x 8 windows = 2**16 = 65536 coalitions per scenario.
#   $1 = "runall"  -> prep + run all workers in the FOREGROUND and wait
#   $1 = "launch"  -> same but detached (workers get reaped ~45 s after the
#                     launching session ends, so prefer runall)
#   $1 = N         -> run worker N
#   $2 = scenario tag (r_sat | r_lim1), default r_sat
#
# Full enumeration means there is nothing to sample here: every one of the 65536
# coalitions is a real box run, and the resulting table IS the oracle.
# Resumable: workers append and skip what they already recorded.
set -u
W=/mnt/d/atmos/_kppwork
ROOT=/root/nml16
NW=12
TAG=${2:-r_sat}
SUF=${SUF:-}
N=$ROOT/$TAG

_prep () {
  rm -f $W/n16_done_$TAG$SUF.txt
  ls $N > /root/n16_files_$TAG.txt
  echo "[$TAG] total tasks: $(wc -l < /root/n16_files_$TAG.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/n16_files_$TAG.txt > /root/n16_chunk_${TAG}_$K.txt
  done
  for K in $(seq 1 $NW); do
    if [ "${REFRESH:-0}" = "1" ]; then rm -rf /root/kpp7n$K; fi
    [ -d /root/kpp7n$K ] || cp -r /root/kpp7 /root/kpp7n$K
  done
}

if [ "${1:-}" = "runall" ]; then
  _prep
  for K in $(seq 1 $NW); do
    bash "$0" $K $TAG > $W/n16_worker_${TAG}_m8_$K.out 2>&1 &
  done
  wait
  echo "[$TAG] ALL WORKERS EXITED $(date)"
  exit 0
fi

if [ "${1:-}" = "launch" ]; then
  _prep
  for K in $(seq 1 $NW); do
    setsid bash "$0" $K $TAG < /dev/null > $W/n16_worker_${TAG}_m8_$K.out 2>&1 &
  done
  sleep 6
  echo "[$TAG] LAUNCHED $(date)  workers alive: $(pgrep -c -f 'sweep_n16.sh' || echo 0)"
  exit 0
fi

K=$1
D=/root/kpp7n$K
# MDA8 campaign: outputs get an _m8 suffix so the original d3_mean-only
# tables stay untouched -- v_ppb must reproduce them bit-for-bit (free
# regression check on the driver), v_mda8 is the new primary metric.
O=$W/n16_v_${TAG}_m8${SUF}_$K.csv
P=$W/n16_prog_${TAG}_m8${SUF}_$K.txt
[ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc,nrow" > $O
DONE=/root/n16_done_${TAG}_m8${SUF}_$K.txt
grep -o '^m[0-9]*' $O | sort -u > $DONE 2>/dev/null || true
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
  # one pass: d3_mean (t>=48) and MDA8 (max rolling 8h mean over windows
  # ending at t>=48) -- same rule as exp_mda8.mean_and_mda8(t0=48)
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
  [ $((n % 512)) -eq 0 ] && echo "$K: +$n $(date +%H:%M:%S)" >> $P
done < /root/n16_chunk_${TAG}_$K.txt
echo "$K: DONE $(date)" >> $P
