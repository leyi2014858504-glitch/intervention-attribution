#!/bin/bash
# SHAM-SWITCH sweep: 2 scenarios x 4096 coalitions, EMIS_CO = 0, /root/kpp6.
#   $1 = "launch"  -> partition (into /root, fast) + copy worker dirs + detach
#   $1 = "runall"  -> same, but run the workers in the FOREGROUND and wait.
#                    Detached workers get reaped ~45 s after the launching
#                    session ends; a foreground `wait` keeps the session busy
#                    so the whole batch lives as long as the sandbox allows.
#   $1 = N         -> run worker N
# Resumable: workers append and skip what they already recorded.
set -u
W=/mnt/d/atmos/_kppwork
N=/root/sham_nml
NW=12

_prep () {
  rm -f $W/sham_done.txt
  ls $N > /root/sham_files.txt
  echo "total tasks: $(wc -l < /root/sham_files.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/sham_files.txt > /root/sham_chunk_$K.txt
  done
  for K in $(seq 1 $NW); do
    [ -d /root/kpp6s$K ] || cp -r /root/kpp6 /root/kpp6s$K
  done
}

if [ "${1:-}" = "runall" ]; then
  _prep
  for K in $(seq 1 $NW); do
    bash "$0" $K > $W/sham_worker_$K.out 2>&1 &
  done
  wait
  echo "ALL WORKERS EXITED $(date)"
  exit 0
fi

if [ "${1:-}" = "launch" ]; then
  _prep
  for K in $(seq 1 $NW); do
    setsid bash "$0" $K < /dev/null > $W/sham_worker_$K.out 2>&1 &
  done
  sleep 6
  echo "LAUNCHED $(date)  workers alive: $(pgrep -c -f 'sweep_sham.sh' || echo 0)"
  exit 0
fi

K=$1
D=/root/kpp6s$K
O=$W/sham_v_$K.csv
P=$W/sham_prog_$K.txt
[ -f $O ] || echo "name,v_ppb,nfail,rc,nrow" > $O
DONE=$W/sham_done_$K.txt
grep -o '^[a-z0-9_]*' $O | sort -u > $DONE 2>/dev/null || true
n=0
while read -r b; do
  nm=${b%.nml}
  grep -qx "$nm" $DONE && continue
  cp $N/${nm}.nml $D/cb05.nml
  cd $D
  rm -f traj.dat run.log
  timeout 120 ./cb05.exe > run.log 2>&1
  rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
  nrow=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
  v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END {if(c>0) printf "%.6f", s/c; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "$nm,${v:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  echo "$nm" >> $DONE
  n=$((n+1))
  [ $((n % 512)) -eq 0 ] && echo "$K: +$n $(date +%H:%M:%S)" >> $P
done < /root/sham_chunk_$K.txt
echo "$K: DONE $(date)" >> $P
