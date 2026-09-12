#!/bin/bash
# Archive the time-resolution scale stress test.
set -eu
W=/mnt/d/atmos/_kppwork
A=/mnt/d/atmos/atmos_plan/atmos/chemistry/kpp_cb05
cd $W
for f in build_kpp7.sh gen_scale_masks.py sweep_scale.sh scale_finalize.sh \
         exp_scale.py; do
  cp $f $A/$f
done
cp scale_v.csv      $A/scale_v.csv
cp scale_result.txt $A/results_scale.txt
cp /mnt/d/atmos/_gate7.txt $A/results_gate7.txt
echo archived
ls $A | grep -Ei "scale|gate7"
