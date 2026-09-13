#!/bin/bash
# Build /root/kppa = the dummy backend with an EXTENDED trajectory output.
#
# Why: to anchor the box to AQS we must fit the observed diurnal cycles of O3,
# NOx, CO and the lumped VOC. The stock traj.dat only carries O3 and NOx, so the
# CO and PAR comparison has nothing to fit against.
#
# The new columns are APPENDED, so the existing 7-column parsers (kpp_box.metrics)
# keep working untouched, and the V-gate value must stay bit-for-bit 21.182293 --
# a change in that number would mean we perturbed the physics.
set -e
W=/mnt/d/atmos/_kppwork
R=/root/kppa

rm -rf $R
cp -r /root/kpp5 $R
cd $R
rm -f cb05.nml traj.dat state_out.dat run.log

python3 - <<'PY'
p = 'cb05_Main.f90'
s = open(p).read()

old_hdr = "      WRITE(20,'(A)') '# t_h  o3_ppb  nox_ppb  temp_K  Sday  kdil_per_h  emis_nox_ppbh'"
new_hdr = ("      WRITE(20,'(A)') '# t_h  o3_ppb  nox_ppb  temp_K  Sday  kdil_per_h  emis_nox_ppbh"
           "  co_ppb  no_ppb  no2_ppb  par  tol  xyl  isop_ppb'")
assert old_hdr in s, "header pattern not found"
s = s.replace(old_hdr, new_hdr, 1)

old_w = """        WRITE(20,'(7ES13.5)') T/3600.0d0, &
              C(ind_o3)/PPB2N, &
              (C(ind_no)+C(ind_no2))/PPB2N, &
              TEMP, SDAY, K_DIL*3600.0d0, ENOX"""
new_w = """        WRITE(20,'(14ES13.5)') T/3600.0d0, &
              C(ind_o3)/PPB2N, &
              (C(ind_no)+C(ind_no2))/PPB2N, &
              TEMP, SDAY, K_DIL*3600.0d0, ENOX, &
              C(ind_co)/PPB2N, &
              C(ind_no)/PPB2N, &
              C(ind_no2)/PPB2N, &
              C(ind_par)/PPB2N, &
              C(ind_tol)/PPB2N, &
              C(ind_xyl)/PPB2N, &
              C(ind_isop)/PPB2N"""
assert old_w in s, "write pattern not found"
s = s.replace(old_w, new_w, 1)

open(p, 'w').write(s)
print("driver patched: traj now 14 columns")
PY

make -f Makefile_cb05 > mka.log 2>&1 || true
grep -qi "error" mka.log && { echo MAKE_FAILED; grep -i -B2 -A4 error mka.log | head -40; exit 1; }
echo MAKE_OK

# --- regression: the frozen reference must not move ------------------------
cat > cb05.nml <<'NML'
&RUNCFG
  TEND_S = 259200.0,
  DT_S = 300.0,
  STATE_IN = 'none',
  EMIS_NOX = 16.0, EMIS_VOC = 12.0, EMIS_ISO = 1.0, EMIS_CO = 40.0,
  NOX_FACTOR = 1.0, VOC_FACTOR = 1.0,
  NOXW = 1.0, VOCW = 1.0,
  WKEND = .FALSE., DUMMYF = .TRUE., CONSTJ = .FALSE., C_H2O = 5.0e18,
  KDIL_DAY = 1500.0, KDIL_NIGHT = 86400.0,
  BG_O3 = 30.0, BG_NO2 = 5.0, BG_NO = 1.0, BG_CO = 120.0,
  BG_VOC = 20.0, BG_ISO = 0.2,
  TMIN_K = 288.0, TMAX_K = 301.0,
/
NML
rm -f traj.dat run.log
timeout 300 ./cb05.exe > run.log 2>&1
v=$(awk 'NR>1 && $1>=48 {s+=$2; n++} END {if(n>0) printf "%.6f", s/n; else printf "NaN"}' traj.dat)
echo "V-gate check: v48=$v   expect=21.182293"
awk 'NR<=3 || (NR>192 && NR<197) {print "  cols:", NF, $0}' traj.dat | head -8
echo KPPA_DONE
