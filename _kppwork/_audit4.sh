#!/bin/bash
# Audit verification, part 4: does the anchor fit include O3? where do the
# 620/1370/1850 run-count numbers come from?
set -u
cd /mnt/d/atmos

echo "=== A. TARGETS in _fit_anchor.py (what enters the cost)"
grep -nB2 -A12 '^TARGETS' _fit_anchor.py

echo
echo "=== B. where do 620 / 1370 / 1850 appear?"
grep -rn '\b620\b' _kppwork/*.txt _kppwork/*.py *.txt 2>/dev/null | grep -v Binary | head -8
echo "--"
grep -rn '\b1370\b\|\b1850\b' _kppwork/*.txt _kppwork/*.py 2>/dev/null | head -8

echo
echo "=== C. exp_n16.py: how does it define runs-for-eps?"
grep -n 'def probe\|eps\|budgets\|runs\b' _kppwork/exp_n16.py | head -25

echo
echo "=== D. any output file of exp_n16.py?"
ls -la _kppwork/*.txt | grep -iE 'n16|spectr|scale'