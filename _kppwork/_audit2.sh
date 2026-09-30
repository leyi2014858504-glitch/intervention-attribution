#!/bin/bash
# Audit verification, part 2: metrics used by each analysis, source of the
# 1.6-1.8x cost law, and whether O3 enters the anchoring fit cost.
set -u
cd /mnt/d/atmos

echo "=== A. which metric column does each analysis read?"
for s in _kppwork/exp_attrib.py _kppwork/exp_n16_bench.py _kppwork/exp_scale.py \
         _kppwork/exp_calib.py _kppwork/exp_selfdiag.py _kppwork/_n16_spectrum.py \
         _kppwork/_m8_spectrum.py _kppwork/exp_sham.py; do
  [ -f "$s" ] || continue
  hits=$(grep -nE 'v_mda8|v_ppb|d3_mean|\["v_|header|fieldnames|row\[' "$s" | head -4)
  echo "-- $s"
  if [ -n "$hits" ]; then echo "$hits" | sed 's/^/     /'; else echo "     (no explicit metric string)"; fi
done

echo
echo "=== B. cost_frontier.txt (first 30 lines)"
head -30 _kppwork/cost_frontier.txt

echo
echo "=== C. bench_n16.txt (first 30 lines)"
head -30 _kppwork/bench_n16.txt

echo
echo "=== D. does the anchor fit put O3 in the cost? (make_cost in _fit_anchor.py)"
grep -nA6 -E 'def make_cost|def cost' _fit_anchor.py | head -40

echo
echo "=== E. d3_mean definition in the campaign awk (spin-up convention)"
grep -n 'NR>1 && \$1' _kppwork/sweep_n16.sh | head -3
grep -n 'NR>1' _kppwork/_sweep_dp.sh | head -3