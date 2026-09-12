#!/bin/bash
# Scale test sweep over nmldir_scale (n=24, 3 sectors x 8 three-hour windows).
#   $1 = worker index 1..NW   -> run that worker
#   $1 = "launch"             -> partition, copy build dirs, detach workers
# Pure bash + setsid, same pattern as the 8192-run sweep.  Resumable.
set -u
W=/mnt/d/atmos/_kppwork
N=$W/nmldir_scale
NW=6

if [ "${1:-}" = "launch" ]; then
  # NOTE: per-worker csvs are deliberately NOT deleted -- the worker resumes
  # from them, so relaunching after a reaped session loses no work.
  rm -f $W/scale_done.txt
  ls $N > $W/scale_files.txt
  echo "total tasks: $(wc -l < $W/scale_files.txt)"
  for K in $(seq 1 $NW); do
    awk -v k=$K -v n=$NW 'NR % n == (k % n)' $W/scale_files.txt > $W/scale_chunk_$K.txt
    echo "chunk $K: $(wc -l < $W/scale_chunk_$K.txt)"
  done
  for K in $(seq 1 $NW); do
    rm -rf /root/kpp7w$K
    cp -r /root/kpp7 /root/kpp7w$K
  done
  for K in $(seq 1 $NW); do
    setsid bash "$0" $K < /dev/null > $W/scale_worker_$K.out 2>&1 &
  done
  sleep 6
  echo "LAUNCHED $(date)  workers alive: $(pgrep -c -f 'sweep_scale.sh' || echo 0)"
  exit 0
fi

K=$1
D=/root/kpp7w$K
O=$W/scale_v_$K.csv
P=$W/scale_prog_$K.txt
[ -f $O ] || echo "name,v_ppb,nfail,rc,nrow" > $O
DONE=$W/scale_done_$K.txt
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
  v=$(awk 'NR>1 && $1>=48 {s+=$2; c++} END {if(c>0) printf "%.6f", s/c; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "$nm,${v:-NaN},${nf:-NA},$rc,${nrow:-0}" >> $O
  echo "$nm" >> $DONE
  n=$((n+1))
  [ $((n % 256)) -eq 0 ] && echo "$K: +$n at $b $(date +%H:%M:%S)" >> $P
done < $W/scale_chunk_$K.txt
echo "$K: DONE $(date)" >> $P
