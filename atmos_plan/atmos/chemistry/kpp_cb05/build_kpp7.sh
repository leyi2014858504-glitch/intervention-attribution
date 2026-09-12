#!/bin/bash
# Build /root/kpp7 = the 12-player backend with EIGHT 3-hour windows per
# sector instead of four 6-hour ones (3 sectors x 8 windows = 24 players).
#
# The 4-window game must be an exact SUB-LATTICE of the 8-window game: turning
# off both halves of a 6-hour block has to reproduce the 4-window result.  That
# is the V-gate below, and it is what makes the finer grid trustworthy.
#
# Only the driver changes; the mechanism is untouched (so KPP is not re-run and
# species indices cannot move).
set -e
W=/mnt/d/atmos/_kppwork
R=/root/kpp7

rm -rf $R
cp -r /root/kpp6 $R
cd $R
rm -f cb05.nml traj.dat state_out.dat run.log
rm -f cb05_Main.f90 cb05_Main.o

cp $W/cb05d_Main.f90 cb05_Main.f90
dos2unix -q cb05_Main.f90

python3 - <<'PY'
import re, io
p = 'cb05_Main.f90'
s = open(p).read()
subs = [
 ("REAL(kind=dp)    :: NOXW(4), VOCW(4), COW(4)",
  "REAL(kind=dp)    :: NOXW(8), VOCW(8), COW(8)"),
 ("!    NOXW/VOCW/COW = per-window (4 x 6h) multipliers on the NOx / VOC / CO",
  "!    NOXW/VOCW/COW = per-window (8 x 3h) multipliers on the NOx / VOC / CO"),
 ("        IWIN = MIN(4, INT(HH/6.0d0) + 1)",
  "        IWIN = MIN(8, INT(HH/3.0d0) + 1)"),
 ("!       6-hour intervention windows: [0,6) [6,12) [12,18) [18,24)",
  "!       3-hour intervention windows: [0,3) [3,6) ... [21,24)"),
]
for a, b in subs:
    assert a in s, "pattern not found: %r" % a[:60]
    s = s.replace(a, b, 1)
open(p, 'w').write(s)
print("driver patched: 8 x 3h windows")
PY

make -f Makefile_cb05 > mk8.log 2>&1 || true
grep -qi "error" mk8.log && { echo MAKE_FAILED; grep -i -B2 -A4 error mk8.log | head -40; exit 1; }
echo MAKE_OK

run_gate () {  # $1 label, $2 expected, $3 NOXW(8), $4 VOCW(8), $5 COW(8)
  cat > cb05.nml <<NML
&RUNCFG
  TEND_S = 259200.0,
  DT_S = 300.0,
  STATE_IN = 'none',
  EMIS_NOX = 16.0, EMIS_VOC = 12.0, EMIS_ISO = 1.0, EMIS_CO = 40.0,
  NOX_FACTOR = 1.0, VOC_FACTOR = 1.0,
  NOXW = $3,
  VOCW = $4,
  COW  = $5,
  WKEND = .FALSE., DUMMYF = .TRUE., CONSTJ = .FALSE., C_H2O = 5.0e18,
  KDIL_DAY = 1500.0, KDIL_NIGHT = 86400.0,
  BG_O3 = 30.0, BG_NO2 = 5.0, BG_NO = 1.0, BG_CO = 120.0,
  BG_VOC = 20.0, BG_ISO = 0.2,
  TMIN_K = 288.0, TMAX_K = 301.0,
/
NML
  rm -f traj.dat run.log
  timeout 300 ./cb05.exe > run.log 2>&1
  v=$(awk 'NR>1 && $1>=48 {s+=$2; n++} END {if(n>0) printf "%.6f", s/n; else printf "NaN"}' traj.dat 2>/dev/null)
  echo "  $1  v48=$v   expect=$2"
}

ALL="1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0"
# NOx 6h block W2 = 3h windows #3,#4 (0-based 2,3); W4 = #7,#8 (0-based 6,7)
NOXW_W2OFF="1.0,1.0,0.0,0.0,1.0,1.0,1.0,1.0"
NOXW_W4OFF="1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0"
NOX_ALL_OFF="0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0"
# VOC 6h block W4 = 3h #7,#8
VOCW_W4OFF="1.0,1.0,1.0,1.0,1.0,1.0,0.0,0.0"

{
  echo "=== V-gate: the 4-window game is a sub-lattice of the 8-window game ==="
  run_gate "all 8 windows ON          (m000)" "21.182293" "$ALL" "$ALL" "$ALL"
  run_gate "NOx 6h-block W2 off       (m002)" "27.371391" "$NOXW_W2OFF" "$ALL" "$ALL"
  run_gate "NOx 6h-block W4 off       (m008)" "34.342422" "$NOXW_W4OFF" "$ALL" "$ALL"
  run_gate "NOx all 8 off             (m015)" "46.762071" "$NOX_ALL_OFF" "$ALL" "$ALL"
  run_gate "VOC all 8 off             (m240)" "10.247975" "$ALL" "$NOX_ALL_OFF" "$ALL"
  run_gate "NOx+VOC all off           (m255)" "43.814426" "$NOX_ALL_OFF" "$NOX_ALL_OFF" "$ALL"
  echo ""
  echo "=== half-block asymmetry: split W2 into its two 3h halves ==="
  run_gate "NOx 3h#3 off only         (no 4-window counterpart)" "-" \
           "1.0,1.0,0.0,1.0,1.0,1.0,1.0,1.0" "$ALL" "$ALL"
  run_gate "NOx 3h#4 off only         (no 4-window counterpart)" "-" \
           "1.0,1.0,1.0,0.0,1.0,1.0,1.0,1.0" "$ALL" "$ALL"
} > /mnt/d/atmos/_gate7.txt 2>&1
echo GATE7_DONE
