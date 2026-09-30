#!/bin/bash
# Release slimming: delete WeChat screenshots (untracked), remove tracked
# junk logs/bookkeeping, report repo stats.
set -u
cd /mnt/d/atmos

rm -f _kppwork/mmexport*.png _kppwork/mmexport*.jpg
echo "mmexport deleted from disk (they were untracked)"

git rm -q out35.log _kppwork/calib_run.log _kppwork/fit_deerpark_2017s.log \
  _kppwork/co_list.txt _kppwork/scen_list.txt
echo "tracked junk removed (5 files)"

echo "--- stray build outputs check:"
git ls-files | grep -E 'cb05\.nml$|\.dat$|atoms_red$|\.mod$|\.exe$' || echo "  none"

echo "--- repo stats:"
echo "  tracked files: $(git ls-files | wc -l)"
git count-objects -vH | grep size-pack
echo "--- status:"
git status --short
