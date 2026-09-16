"""Cross-validate the OLD n16 v48 shards against the NEW m8 campaign's
v_ppb column, mask by mask, both scenarios.  If they agree bit-for-bit the
old n=16 tables were NOT corrupted and the v48 divergence story stands;
if not, the v48 n=16 numbers must be re-derived from the new tables."""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))

for scen in ("r_sat", "r_lim1"):
    old, new = {}, {}
    for K in range(1, 13):
        f = os.path.join(HERE, "n16_v_%s_%d.csv" % (scen, K))
        if os.path.exists(f):
            for r in csv.DictReader(open(f)):
                if r["v_ppb"] not in ("", "NaN"):
                    old[int(r["name"][1:])] = float(r["v_ppb"])
        f = os.path.join(HERE, "n16_v_%s_m8_%d.csv" % (scen, K))
        if os.path.exists(f):
            for r in csv.DictReader(open(f)):
                if r["v_ppb"] not in ("", "NaN"):
                    new[int(r["name"][1:])] = float(r["v_ppb"])
    common = set(old) & set(new)
    bad = [(m, old[m], new[m]) for m in common if abs(old[m] - new[m]) > 1e-6]
    print("%s: old=%d new=%d common=%d disagree=%d" %
          (scen, len(old), len(new), len(common), len(bad)))
    for m, o, n in sorted(bad)[:10]:
        print("   m=%-6d old=%10.4f new=%10.4f" % (m, o, n))
