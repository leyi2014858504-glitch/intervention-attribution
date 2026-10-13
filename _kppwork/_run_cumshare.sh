#!/bin/bash
# R_q cumulative spectral-mass share, both metrics, all n (zero chemistry runs).
cd /mnt/d/atmos
python3 _kppwork/_cum_share.py > _kppwork/_cum_share.log 2>&1
echo "EXIT $?"
cat _kppwork/_cum_share.log