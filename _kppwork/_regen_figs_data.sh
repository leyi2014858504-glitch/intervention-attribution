#!/bin/bash
# Regenerate results_quota.txt (Fig 7) and results_co12.txt (Table 3 / Fig 3)
# -- pure post-processing on the released oracle tables, zero chemistry runs.
cd /mnt/d/atmos/_kppwork
python3 exp_quota.py > /mnt/d/atmos/results_quota.txt 2>&1
echo "QUOTA EXIT $?"
python3 exp_co12.py > _co12_run.log 2>&1
echo "CO12 EXIT $?"
cp co12_analysis.txt /mnt/d/atmos/results_co12.txt 2>/dev/null || echo "no co12_analysis.txt"
tail -3 _co12_run.log
grep -n "0.10\|0.31\|max" /mnt/d/atmos/results_co12.txt | head -8