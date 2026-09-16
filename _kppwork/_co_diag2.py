"""Which masks disagree between co_v.csv and co_v_m8?"""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))

for scen in ("r_sat", "r_lim1"):
    old, new = {}, {}
    for r in csv.DictReader(open(os.path.join(HERE, "co_v.csv"))):
        if r["name"].startswith(scen + "_m") and r["v_ppb"] not in ("", "NaN"):
            old[int(r["name"].rsplit("_m", 1)[1])] = float(r["v_ppb"])
    for f in glob.glob(os.path.join(HERE, "co_v_m8_*.csv")):
        for r in csv.DictReader(open(f)):
            if r["name"].startswith(scen + "_m") and r["v_ppb"] not in ("", "NaN"):
                new[int(r["name"].rsplit("_m", 1)[1])] = float(r["v_ppb"])
    bad = [(m, old[m], new[m]) for m in old if abs(new[m] - old[m]) > 1e-6]
    print("%s: %d/%d masks disagree" % (scen, len(bad), len(old)))
    for m, o, n in sorted(bad)[:12]:
        print("   m=%-5d old=%9.4f new=%9.4d pop=%d" % (m, o, n, bin(m).count("1")))
