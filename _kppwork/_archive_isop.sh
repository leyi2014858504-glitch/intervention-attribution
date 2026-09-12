#!/bin/bash
# Archive the CO-coefficient evidence + the isoprene probe.
set -eu
W=/mnt/d/atmos/_kppwork
A=/mnt/d/atmos/atmos_plan/atmos/chemistry/kpp_cb05
cd $W
cp exp_co12_coeff.py probe_isop.py probe_isop.sh run_isop_probe.sh $A/
cp co12_coeff.txt $A/results_co12_coeff.txt
cp isop_probe.csv $A/isop_probe.csv
echo "archived:"
ls $A | grep -Ei "isop|coeff"
