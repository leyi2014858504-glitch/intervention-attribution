#!/bin/bash
# Deeper probe: are the worker shells and the parent sweep still alive?
date
echo "--- worker bash processes ---"
pgrep -af "sweep_co_w.sh" || echo "NO sweep_co_w.sh processes"
echo "--- parent ---"
pgrep -af "sweep_co_par.sh" || echo "NO sweep_co_par.sh"
echo "--- cb05 ---"
pgrep -af "cb05.exe" || echo "NO cb05.exe"
echo "--- row counts (repeat twice, 8s apart) ---"
for K in 1 2 3 4 5 6; do echo "w$K $(wc -l < /mnt/d/atmos/_kppwork/co_v_$K.csv 2>/dev/null)"; done
sleep 8
echo "--- after 8s ---"
for K in 1 2 3 4 5 6; do echo "w$K $(wc -l < /mnt/d/atmos/_kppwork/co_v_$K.csv 2>/dev/null)"; done
cat /mnt/d/atmos/_kppwork/co_done.txt 2>/dev/null && echo "(sweep finished)"
ls -la /mnt/d/atmos/_kppwork/co_v.csv 2>/dev/null
