#!/bin/bash
# Top dividends at n=16 (MDA8) with member names -> results_top_dividends.txt
cd /mnt/d/atmos
bash _kppwork/_top_dividends.sh > results_top_dividends.txt 2>&1
echo "EXIT $?"
cat results_top_dividends.txt