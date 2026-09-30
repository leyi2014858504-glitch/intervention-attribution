#!/bin/bash
# Final wrap: commit the README licence note, verify repo cleanliness.
set -u
cd /mnt/d/atmos

git add README.md
git -c user.name=unknown -c user.email=leyi2014858504@gmail.com commit -q -m 'README: note RELEASE_NOTES.md and the pending licence decision'
echo "latest commit: $(git log --oneline -1)"

echo "--- screenshots on disk:"
ls _kppwork/mmexport* 2>/dev/null || echo "  none (deleted)"

echo "--- tracked files: $(git ls-files | wc -l)"
echo "--- working tree:"
git status --short | head -10
echo "--- branch vs origin:"
git status -sb | head -1