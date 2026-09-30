#!/bin/bash
# Release the spin-up steady-state diagnostic as a result file
# (3-day vs 10-day last-day composite, at the Dallas anchored slice).
cd /mnt/d/atmos
python3 _diag_fit.py > results_spinup_dallas.txt 2>&1
echo "SPINUP EXIT $?"
cat results_spinup_dallas.txt