#!/bin/bash
# Relaunch the 12-player sweep with FULLY DETACHED workers.
#
# Why setsid: a plain `&` + wait put the workers in the sandbox's process
# group, and the whole group was reaped mid-run.  setsid gives each worker its
# own session so it survives the launcher exiting.
#
# Chunks were previously partitioned with `NR % n == k`, which is empty for
# k == n -> worker 6 idled while 1..5 did everything.  Fixed to `k % n`, so
# worker 6 now owns positions 6, 12, 18, ...  Workers 1..5 keep their old
# chunks and resume from their own csv.
set -u
W=/mnt/d/atmos/_kppwork
N=$W/nmldir_co
NW=6

ls $N > $W/co_files.txt
for K in $(seq 1 $NW); do
  awk -v k=$K -v n=$NW 'NR % n == (k % n)' $W/co_files.txt > $W/chunk_$K.txt
done

for K in $(seq 1 $NW); do
  rm -rf /root/kpp6w$K
  cp -r /root/kpp6 /root/kpp6w$K
done
rm -f $W/co_done.txt

for K in $(seq 1 $NW); do
  setsid bash $W/sweep_co_w.sh $K < /dev/null > $W/worker_$K.out 2>&1 &
done
sleep 6
echo "LAUNCHED $(date)"
pgrep -c -f sweep_co_w.sh || echo "no workers alive"
