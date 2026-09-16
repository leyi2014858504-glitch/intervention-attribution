#!/bin/bash
# MDA8 campaign: full n=16 enumeration, both scenarios, dual-metric output.
# Resumable (workers skip names already in their _m8 csv).
S=/mnt/d/atmos/_kppwork/sweep_n16.sh
echo "[m8] start r_sat  $(date)"
bash $S runall r_sat
echo "[m8] start r_lim1 $(date)"
bash $S runall r_lim1
echo "[m8] CAMPAIGN COMPLETE $(date)"
