#!/bin/bash
# Archive the cost-accuracy frontier (zero-new-run analysis).
set -eu
W=/mnt/d/atmos/_kppwork
A=/mnt/d/atmos/atmos_plan/atmos/chemistry/kpp_cb05
cp $W/exp_cost_frontier.py $A/exp_cost_frontier.py
cp $W/cost_frontier.txt $A/results_cost_frontier.txt
echo archived
