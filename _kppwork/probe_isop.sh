#!/bin/bash
# Isoprene-leverage probe, done the SAME way as the 8192-run sweep: pure bash
# driving ./cb05.exe directly (no Python<->WSL per-run round trip, which proved
# flaky once the parallel sweep had churned /root).
#
# Configs (2 scenarios x 6):
#   mask 0    (all windows on)             x EMIS_ISO in {0, 1, 3, 6}
#   mask 240  (anthropogenic VOC all off)  x EMIS_ISO in {1, 0}
# mask bit w = NOx window w+1 off, bit 4+w = VOC window w+1 off.
set -u
D=/root/kpp6
O=/mnt/d/atmos/_kppwork/isop_probe.csv

write_nml () {  # $1 = NOX_FACTOR, $2 = VOC_FACTOR, $3 = VOCW, $4 = ISO
  cat > $D/cb05.nml <<NML
&RUNCFG
  TEND_S = 259200.0,
  DT_S = 300.0,
  STATE_IN = 'none',
  EMIS_NOX = 16.0, EMIS_VOC = 12.0, EMIS_ISO = $4, EMIS_CO = 40.0,
  NOX_FACTOR = $1, VOC_FACTOR = $2,
  NOXW = 1.0,1.0,1.0,1.0, VOCW = $3, COW = 1.0,1.0,1.0,1.0,
  WKEND = .FALSE., DUMMYF = .TRUE., CONSTJ = .FALSE., C_H2O = 5.0e18,
  KDIL_DAY = 1500.0, KDIL_NIGHT = 86400.0,
  BG_O3 = 30.0, BG_NO2 = 5.0, BG_NO = 1.0, BG_CO = 120.0,
  BG_VOC = 20.0, BG_ISO = 0.2,
  TMIN_K = 288.0, TMAX_K = 301.0,
/
NML
}

metrics () {  # $1 = traj path -> "mean24,mda8"
  awk '
    NR>1 { t[NR]=$1; o[NR]=$2 }
    END {
      n=0; s=0
      for (i=2;i<=NR;i++) if (t[i]>=48) { s+=o[i]; n++ }
      mean = (n>0) ? s/n : 0
      W=96; best=0
      for (i=2;i<=NR;i++) {
        if (t[i]<48 || i-W+1<2) continue
        ss=0; for (j=i-W+1;j<=i;j++) ss+=o[j]
        if (ss/W>best) best=ss/W
      }
      printf "%.6f,%.6f", mean, best
    }' "$1"
}

echo "scenario,mask,iso,mean24,mda8,nfail,nrow" > $O
run_one () {  # $1 label, $2 noxf, $3 vocf, $4 vocw, $5 iso, $6 mask
  write_nml "$2" "$3" "$4" "$5"
  cd $D
  rm -f traj.dat run.log
  timeout 120 ./cb05.exe > run.log 2>&1
  rc=$?
  nf=$(grep -o 'NFAIL=[ ]*[0-9]*' run.log | tr -dc '0-9')
  nrow=$(awk 'END{print NR-1}' traj.dat 2>/dev/null)
  m=$(metrics $D/traj.dat)
  echo "$1,$6,$5,$m,${nf:-NA},${nrow:-0}" >> $O
  echo "$1 mask=$6 iso=$5 -> $m rc=$rc nfail=${nf:-NA} nrow=${nrow:-0}" \
    >> /mnt/d/atmos/_kppwork/isop_probe_progress.txt
}

: > /mnt/d/atmos/_kppwork/isop_probe_progress.txt
echo "=== A: ISO sweep, all windows on ==="
for iso in 0.0 1.0 3.0 6.0; do
  run_one r_sat  1.0 1.0 "1.0,1.0,1.0,1.0" $iso 0
  run_one r_lim1 0.5 1.5 "1.0,1.0,1.0,1.0" $iso 0
done
echo "=== B: 2x2 anthro-VOC(off) x biogenic-ISO ==="
for iso in 1.0 0.0; do
  run_one r_sat  1.0 1.0 "0.0,0.0,0.0,0.0" $iso 240
  run_one r_lim1 0.5 1.5 "0.0,0.0,0.0,0.0" $iso 240
done
echo "PROBE_DONE $(date)" >> /mnt/d/atmos/_kppwork/isop_probe_progress.txt
