#!/bin/bash
# Rebuild+gate the 12-player backend, then generate the sham namelists.
set -u
W=/mnt/d/atmos/_kppwork
rm -rf $W/nmldir_sham          # drop the half-written Windows-side copy
cd $W
dos2unix -q gen_sham.py sweep_sham.sh
echo "=== rebuild + backward-compat gate (/root/kpp6) ==="
bash build_kpp6.sh
echo ""
echo "=== gate output ==="
cat $W/../_build6.txt 2>/dev/null || true
echo ""
echo "=== generate sham namelists into /root/sham_nml ==="
python3 gen_sham.py
echo "namelists in /root/sham_nml: $(ls /root/sham_nml 2>/dev/null | wc -l)"
