#!/bin/bash
# Calibration runner.  Log goes to the Windows-side mount (persistent across
# WSL restarts; /tmp is tmpfs and vanishes).  Threads capped so the VM survives
# the n=16 dense lstsq.
set -u
cd /mnt/d/atmos/_kppwork
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
python3 -u exp_calib.py > /mnt/d/atmos/_kppwork/calib_run.log 2>&1
echo "EXIT $?" >> /mnt/d/atmos/_kppwork/calib_run.log
