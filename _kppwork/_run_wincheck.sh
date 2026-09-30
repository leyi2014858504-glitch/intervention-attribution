#!/bin/bash
# Release the MDA8-window check output (was run ad hoc; now saved as a result file).
cd /mnt/d/atmos
python3 _kppwork/_mda8_wincheck.py > _kppwork/results_mda8_wincheck.txt 2>&1
echo "WINCHECK EXIT $?"
cat _kppwork/results_mda8_wincheck.txt