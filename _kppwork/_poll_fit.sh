#!/bin/bash
# poll for the Deer Park anchor fit to finish (max ~25 min)
F=/mnt/d/atmos/results_fit_deerpark_2017s.txt
for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30; do
  if [ -f "$F" ]; then
    echo "FIT_DONE"
    cat "$F"
    exit 0
  fi
  if ! pgrep -f _fit_anchor >/dev/null; then
    echo "PROC_EXITED_NO_FILE"
    tail -5 /mnt/d/atmos/_kppwork/fit_deerpark_2017s.log
    exit 1
  fi
  sleep 30
done
echo "STILL_RUNNING_AFTER_15MIN"
tail -4 /mnt/d/atmos/_kppwork/fit_deerpark_2017s.log
