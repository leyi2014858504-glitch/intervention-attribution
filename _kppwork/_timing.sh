#!/bin/bash
# One-off: time a single 3-day box run to size the Dallas anchor job.
cd /mnt/d/atmos
python3 - <<'EOF'
import time, sys
sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box
t = time.time()
kpp_box.run_scenario({}, "timing", "anchor")
print("one 3-day run: %.1f s" % (time.time() - t))
EOF