#!/bin/bash
# Hourly-MDA8 Taylor-2 coefficient runs for the CITY-ANCHORED HDDM check.
# Same 512 derivative runs as sweep_hddm_city.sh, with the regulatory hourly
# MDA8 metric (8-h means of hourly values ending at h = 48..71).
set -u
W=/mnt/d/atmos/_kppwork
N=/root/nmlhddm
NW=12
LABELS="dallas_2016w deerpark_2017s"

if [ "${1:-}" = "runall" ]; then
  for T in $LABELS; do
    ls $N/$T > /root/hddmh_files_$T.txt
    echo "[hddmh $T] tasks: $(wc -l < /root/hddmh_files_$T.txt)"
    for K in $(seq 1 $NW); do
      awk -v k=$K -v n=$NW 'NR % n == (k % n)' /root/hddmh_files_$T.txt > /root/hddmh_chunk_${T}_$K.txt
    done
    for K in $(seq 1 $NW); do
      bash "$0" $K $T > $W/hddmh_worker_${T}_$K.out 2>&1 &
    done
    wait
    n=$(cat $W/hddm_h_v_${T}_*.csv 2>/dev/null | grep -c '^[a-z]')
    echo "[hddmh $T] workers done $(date) -- rows $n"
  done
  echo "[hddmh] ALL WORKERS EXITED $(date)"
  exit 0
fi

K=$1
TAG=$2
D=/root/kpp7n$K
O=$W/hddm_h_v_${TAG}_$K.csv
[ -d $D ] || cp -r /root/kpp7 $D
[ -f $O ] || echo "name,v_ppb,v_mda8,nfail,rc" > $O
DONE=/root/hddmh_done_${TAG}_$K.txt
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
  for(h=40;h<=71;h++){
    s=0;c=0;
    for(i=1;i<=n;i++) if(t[i]>=h && t[i]<h+1){s+=o[i];c++}
    hm[h]=(c>0?s/c:0)
  }
  best=-1;
  for(h=48;h<=71;h++){
    a=0; for(j=h-7;j<=h;j++) a+=hm[j]; a/=8;
    if(a>best)best=a
  }
  if(best>=0)printf "%.6f",best; else printf "NaN"
}' traj.dat 2>/dev/null)
EOF2
  echo "$nm,${v:-NaN},${v8:-NaN},${nf:-NA},$rc" >> $O
  echo "$nm" >> $DONE
done < /root/hddmh_chunk_${TAG}_$K.txt
echo "$K: DONE $(date)" >> $W/hddmh_prog_${TAG}_$K.txt
