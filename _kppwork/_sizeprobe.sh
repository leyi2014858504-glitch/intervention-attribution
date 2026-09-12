#!/bin/bash
cd /mnt/d/atmos
echo "=== total files (tracked candidates) ==="
find . -type f | wc -l
echo "=== files per top-level dir ==="
for d in _kppwork atmos_plan; do
  printf "%-14s %s\n" "$d" "$(find $d -type f | wc -l)"
done
echo "=== root-level scratch files (_*.txt / _*.sh) ==="
ls -1 _*.txt _*.sh _*.html 2>/dev/null | wc -l
echo "=== generated namelists ==="
for d in _kppwork/nmldir _kppwork/nmldir_scen _kppwork/nmldir_co _kppwork/nmldir_scale; do
  printf "%-26s %s\n" "$d" "$(ls $d 2>/dev/null | wc -l)"
done
echo "=== g4runs artifacts ==="
find _kppwork/g4runs -type f 2>/dev/null | wc -l
echo "=== sizes ==="
du -sh _kppwork atmos_plan 2>/dev/null
du -sh . 2>/dev/null
