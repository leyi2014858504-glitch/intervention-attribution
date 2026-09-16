#!/bin/bash
# detached runner for the MDA8 self-diagnostic validation
cd /mnt/d/atmos/_kppwork
python3 exp_selfdiag.py "" v_mda8 > /tmp/sd_m8.log 2>&1
echo "EXIT $?" >> /tmp/sd_m8.log
echo DONE >> /tmp/sd_m8.log
