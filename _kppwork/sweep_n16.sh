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
N=$ROOT/$TAG

_prep () {
  rm -f $W/n16_done_$TAG.txt
  ls $N > /root/n16_files_$TAG.txt
  echo "[$TAG] total tasks: $(wc -l < /root/n16_files_$TAG.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/n16_files_$TAG.txt > /root/n16_chunk_${TAG}_$K.txt
  done
  for K in $(seq 1 $NW); do
    [ -d /root/kpp7n$K ] || cp -r /root/kpp7 /root/kpp7n$K
  done
}

if [ "${1:-}" = "runall" ]; then
  _prep
  for K in $(seq 1 $NW); do
    bash "$0" $K $TAG > $W/n16_worker_${TAG}_$K.out 2>&1 &
  done
  wait
  echo "[$TAG] ALL WORKERS EXITED $(date)"
  exit 0
fi

if [ "${1:-}" = "launch" ]; then
  _prep
  for K in $(seq 1 $NW); do
    setsid bash "$0" $K $TAG < /dev/null > $W/n16_worker_${TAG}_$K.out 2>&1 &
  done
  sleep 6
  echo "[$TAG] LAUNCHED $(date)  workers alive: $(pgrep -c -f 'sweep_n16.sh' || echo 0)"
  exit 0
fi

K=$1
D=/root/kpp7n$K
O=$W/n16_v_${TAG}_$K.csv
P=$W/n16_prog_${TAG}_$K.txt
[ -f $O ] || echo "name,v_ppb,nfail,rc,nrow" > $O
DONE=/root/n16_done_${TAG}_$K.txt
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
  v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END{if(c>0) printf "%.6f", s/c; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "$nm,${v:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  echo "$nm" >> $DONE
  n=$((n+1))
  [ $((n % 512)) -eq 0 ] && echo "$K: +$n $(date +%H:%M:%S)" >> $P
done < /root/n16_chunk_${TAG}_$K.txt
echo "$K: DONE $(date)" >> $P
