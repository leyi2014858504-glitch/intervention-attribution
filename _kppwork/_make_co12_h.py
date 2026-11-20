import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 12
FULL = 1 << N

hm = {}
for k in range(1, 7):
    with open(os.path.join(W, "co_hm_v_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            scen, ms = row["name"].rsplit("_m", 1)
            h = np.array([float(row["h%02d" % hh]) for hh in range(40, 72)])
            A = np.array([h[i - 7:i + 1].mean() for i in range(7, 32)])
            hm.setdefault(scen, {})[int(ms)] = A.max()

for scen in hm:
    assert len(hm[scen]) == FULL, (scen, len(hm[scen]))

shards = [[] for _ in range(6)]
for mask in range(FULL):
    shards[mask % 6].append(mask)
for k, masks in enumerate(shards, 1):
    with open(os.path.join(W, "co_vh_m8_%d.csv" % k), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "v_ppb", "v_mda8", "nfail", "rc"])
        for scen in ("r_sat", "r_lim1"):
            for mask in masks:
                w.writerow(["%s_m%04d" % (scen, mask), "0",
                            "%.6f" % hm[scen][mask], 0, 0])
for scen in ("r_sat", "r_lim1"):
    print(scen, "all-on %.4f" % hm[scen][0])
print("wrote co_vh_m8_1..6.csv")
