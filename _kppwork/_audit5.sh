#!/bin/bash
# Audit verification, part 5: provenance of the Dallas anchor signature
# (Table 9 claims "7/8 NOx windows raise O3") and of the lvg baselines.
set -u
cd /mnt/d/atmos

echo "=== A. what scenario/base does _local_vs_global.py use?"
grep -nE 'SCEN|scen|ANCHOR|EMIS_NOX|dallas|deerpark|BASE' _local_vs_global.py | head -20

echo
echo "=== B. any script that runs a Dallas-anchored lvg/regime test?"
grep -rlnE 'dallas' --include='*.py' . 2>/dev/null | head -12

echo
echo "=== C. search results files for a Dallas per-window NOx signature"
grep -rlnE '7/8|NOx.*windows' results_*.txt _kppwork/*.txt 2>/dev/null | head -8

echo
echo "=== D. is the dallas fit file on disk, and what does it contain?"
ls -la results_fit_dallas_2016w.txt 2>/dev/null && cat results_fit_dallas_2016w.txt

echo
echo "=== E. scan_anchor / anchor-related outputs tracked?"
git ls-files | grep -iE 'scan_anchor|anchor|kdil'