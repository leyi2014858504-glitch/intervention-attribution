#!/bin/bash
# Archive the sham-switch control.
set -eu
W=/mnt/d/atmos/_kppwork
A=/mnt/d/atmos/atmos_plan/atmos/chemistry/kpp_cb05
cd $W
for f in gen_sham.py sweep_sham.sh exp_sham.py sham_finalize.sh; do
  cp $f $A/$f
done
cp sham_v.csv      $A/sham_v.csv
cp sham_result.txt $A/results_sham.txt
echo archived
ls $A | grep -i sham
