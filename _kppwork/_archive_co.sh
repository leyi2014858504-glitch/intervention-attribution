#!/bin/bash
# Archive the 12-player (CO-sector) experiment.
set -eu
W=/mnt/d/atmos/_kppwork
A=/mnt/d/atmos/atmos_plan/atmos/chemistry/kpp_cb05
cd $W
for f in exp_co12.py gen_co_sweep.py sweep_co_w.sh sweep_co_par2.sh \
         co_watch.sh co_finalize.sh build_kpp6.sh cb05d_Main.f90 kpp_box.py; do
  cp $f $A/$f
done
cp co12_analysis.txt $A/results_co12.txt
cp co_v.csv $A/co_v.csv
cp /mnt/d/atmos/_co12.txt $A/results_co12_raw.txt
echo "archived. co-related files:"
ls $A | grep -i "co" || true
echo "total files in archive: $(ls $A | wc -l)"
