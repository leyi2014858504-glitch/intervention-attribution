#!/bin/bash
# Dallas 2016 winter: 8-knob anchoring under the current protocol
# (EMIS x4 + KDIL x2 + BG_O3 + BG_CO, met-bounded tau_BLH and TMIN/TMAX).
# Overwrites the stale 6-knob results_fit_dallas_2016w.txt on success.
cd /mnt/d/atmos
python3 _fit_anchor.py dallas_2016w 60 > _fit_dallas.log 2>&1
echo "FIT EXIT $?"
tail -5 _fit_dallas.log