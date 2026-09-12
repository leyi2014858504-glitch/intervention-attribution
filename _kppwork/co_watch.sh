#!/bin/bash
# Wait for all 6 sweep workers to report DONE, merge their CSVs (dedup by
# coalition name), run the 12-player Moebius analysis, then stop.
# Bounded so a reaped worker cannot hang this forever.
set -u
W=/mnt/d/atmos/_kppwork
for i in $(seq 1 400); do
  c=$(grep -l "DONE" $W/co_prog_*.txt 2>/dev/null | wc -l)
  [ "$c" -ge 6 ] && break
  sleep 10
done
echo "watcher: $c/6 workers reported DONE at $(date)" > $W/co_watch.txt
{
  echo "working dir set: $W"
  head -1 $W/co_v_1.csv > $W/co_v.csv
  for K in 1 2 3 4 5 6; do tail -n +2 $W/co_v_$K.csv; done \
    | awk -F, '!seen[$1]++' >> $W/co_v.csv
  echo "merged rows: $(wc -l < $W/co_v.csv)"
  echo "distinct scenarios:"
  awk -F, 'NR>1 {split($1,a,"_m"); print a[1]}' $W/co_v.csv | sort | uniq -c
  echo "bad rows (nfail!=0 or rc!=0 or nrow!=864):"
  awk -F, 'NR>1 && ($3!="0" || $4!="0" || $5!="864") {print; bad++} END {print "BAD="bad+0}' $W/co_v.csv | head -20
} >> $W/co_watch.txt 2>&1
cd $W && python3 exp_co12.py >> $W/co_watch.txt 2>&1
echo "WATCHER_FINISHED $(date)" >> $W/co_watch.txt
