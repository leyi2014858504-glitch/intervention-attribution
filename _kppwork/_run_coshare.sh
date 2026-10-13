#!/bin/bash
# f_CO(k) = 1 - W_8(k)/W_12(k), both metrics (zero chemistry runs).
cd /mnt/d/atmos
python3 _kppwork/_co_share.py > _kppwork/_co_share.log 2>&1
echo "EXIT $?"
cat _kppwork/_co_share.log