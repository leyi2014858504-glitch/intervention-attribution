#!/bin/bash
cd /mnt/d/atmos/_kppwork
echo "namelists: $(ls nmldir_sham 2>/dev/null | wc -l)  (expect 8192)"
echo "--- sham_sat_m0000.nml ---"
cat nmldir_sham/sham_sat_m0000.nml
