"""Diagnose the co12 v48 regression failure: raw offset vs game change."""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))

old, new = {}, {}
for r in csv.DictReader(open(os.path.join(HERE, "co_v.csv"))):
    if r["name"].startswith("r_sat_m") and r["v_ppb"] not in ("", "NaN"):
        old[int(r["name"].rsplit("_m", 1)[1])] = float(r["v_ppb"])
for f in glob.glob(os.path.join(HERE, "co_v_m8_*.csv")):
    for r in csv.DictReader(open(f)):
        if r["name"].startswith("r_sat_m") and r["v_ppb"] not in ("", "NaN"):
            new[int(r["name"].rsplit("_m", 1)[1])] = float(r["v_ppb"])

print("rows old=%d new=%d" % (len(old), len(new)))
o0, n0 = old[0], new[0]
print("raw baseline (all ON): old=%.6f new=%.6f  offset=%.6f" % (o0, n0, n0 - o0))
d_raw = max(abs(new[m] - old[m]) for m in old)
d_game = max(abs((new[m] - n0) - (old[m] - o0)) for m in old)
print("max|raw diff|   = %.6f" % d_raw)
print("max|game diff|  = %.6f" % d_game)
# a few probes
for m in (1, 16, 256, 4095):
    print("  m=%-5d old=%9.4f new=%9.4f  game_old=%8.4f game_new=%8.4f"
          % (m, old[m], new[m], old[m] - o0, new[m] - n0))
