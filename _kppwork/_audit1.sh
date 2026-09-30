#!/bin/bash
# Audit verification, part 1: do the claimed artifacts exist in the repo?
set -u
cd /mnt/d/atmos

echo "=== 1. n16 oracle tables tracked?"
git ls-files _kppwork/n16_v_* || true

echo
echo "=== 2. row counts, header, NaN counts"
for f in _kppwork/n16_v_r_sat.csv _kppwork/n16_v_r_lim1.csv \
         _kppwork/n16_v_r_sat_m8.csv _kppwork/n16_v_r_lim1_m8.csv; do
  [ -f "$f" ] || { echo "-- $f: MISSING ON DISK"; continue; }
  rows=$(wc -l < "$f")
  size=$(du -k "$f" | cut -f1)
  echo "-- $f: rows=$rows size=${size}K"
  head -2 "$f" | sed 's/^/     /'
  awk -F, 'NR>1{for(i=1;i<=NF;i++) if($i=="NaN"||$i=="") nan[i]++}
           END{printf "     NaN/empty by col:"; for(i in nan) printf " col%d=%d", i, nan[i]; print ""}' "$f"
done

echo
echo "=== 3. audit scripts tracked?"
git ls-files | grep -E '_mda8_wincheck|_grid_check' || echo "  NONE (wincheck/grid scripts not tracked)"

echo
echo "=== 4. anchored-fit result files tracked?"
git ls-files | grep -E 'results_fit_' || echo "  NONE"

echo
echo "=== 5. what .gitignore says about the dallas fit:"
grep -n 'results_fit_dallas' .gitignore || echo "  (no rule)"

echo
echo "=== 6. which metric do the benchmark sources use?"
for s in _kppwork/attrib_bench.txt _kppwork/bench_n16.txt _kppwork/scaling_n16.txt _kppwork/cost_frontier.txt; do
  echo "-- $s:"
  grep -m3 -nE 'v_ppb|v_mda8|d3|MDA8|metric' "$s" 2>/dev/null | sed 's/^/     /' || echo "     (no metric mention)"
done