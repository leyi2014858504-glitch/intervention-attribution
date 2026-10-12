#!/bin/bash
# Re-derive scaling_dual_metrics.txt with q1/q3 columns (zero chemistry runs).
# The script writes the output file itself; stdout goes to a log.
cd /mnt/d/atmos
python3 _kppwork/_dual_bench.py > _kppwork/_dual_bench_run.log 2>&1
echo "DUAL BENCH EXIT $?"
grep -m8 "q1\|  8  \|  16 " _kppwork/scaling_dual_metrics.txt