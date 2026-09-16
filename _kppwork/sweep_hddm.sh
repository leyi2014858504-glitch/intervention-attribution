#!/bin/bash
# Taylor-2 coefficient runs for the HDDM check: 2 scenarios x 512 runs.
# Same worker pattern as sweep_n16.sh (kpp7 builds, resumable, dual metric).
#   $1 = "runall" | worker index
set -u
W=/mnt/d/atmos/_kppwork
N=/root/nmlhddm
NW=12
TAG=${2:-r_sat}

if [ "${1:-}" = "runall" ]; then
  # Scenarios run SEQUENTIALLY: both use the same /root/kpp7n$K build dirs,
  # so overlapping them makes workers clobber each other's cb05.nml/traj.dat.
  # The first attempt did exactly that and wrote plausible-looking but WRONG
  # values (a clobbered namelist still runs), so a partial purge is not enough
  # -- start from clean tables.
  #   $2 = "resume" keeps existing tables (use only after a clean run)
  if [ "${2:-}" != "resume" ]; then
    rm -f $W/hddm_v_*.csv /root/hddm_done_*.txt
    echo "[hddm] tables cleared"
  fi
  for T in r_sat r_lim1; do
    ls $N/$T > /root/hddm_files_$T.txt
    echo "[hddm $T] tasks: $(wc -l < /root/hddm_files_$T.txt)"
    for K in $(seq 1 $NW); do
      awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/hddm_files_$T.txt > /root/hddm_chunk_${T}_$K.txt
    done
    for K in $(seq 1 $NW); do
      bash "$0" $K $T > $W/hddm_worker_${T}_$K.out 2>&1 &
    done
    wait
    n=$(cat $W/hddm_v_${T}_*.csv 2>/dev/null | grep -c '^[a-z]')
    bad=$(cat $W/hddm_v_${T}_*.csv 2>/dev/null | awk -F, '$2=="NaN"||$3=="NaN"' | wc -l)
    echo "[hddm $T] workers done $(date) -- rows $n, NaN $bad"
  done
  echo "[hddm] ALL WORKERS EXITED $(date)"
  exit 0
fi

K=$1
D=/root/kpp7n$K
O=$W/hddm_v_${TAG}_$K.csv
[ -d $D ] || cp -r /root/kpp7 $D
[ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc" > $O
DONE=/root/hddm_done_${TAG}_$K.txt
grep -o '^[a-z0-9_]*' $O | sort -u > $DONE 2>/dev/null || true
while read -r b; do
  nm=${b%.nml}
  grep -qx "$nm" $DONE && continue
  cp $N/$TAG/$b $D/cb05.nml
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
done < /root/hddm_chunk_${TAG}_$K.txt
echo "$K: DONE $(date)" >> $W/hddm_prog_${TAG}_$K.txt
