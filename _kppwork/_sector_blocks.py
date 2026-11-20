"""Sector-level sufficiency check, kernel-free, from the canonical tables (MDA8).

The allocation answers "how much of the joint effect is charged to each sector".
This script answers the complementary policy question: "can the benefit be
obtained by acting on one sector alone, and do the two sectors add up?"

  v(NOx block off), v(VOC block off), v(all off), v(none)
  sum of the eight single NOx (VOC) effects
  block interaction = joint - (NOx block + VOC block)

Writes results_sector_blocks.txt
"""
import csv
import os

import numpy as np

D = os.environ.get("ATMOS_DIR", "d:/atmos")
W = os.path.join(D, "_kppwork")
N = 16
GAMES = {"dallas_2016w": "city", "deerpark_2017s": "city", "r_sat": "synth", "r_lim1": "synth"}
NOX = (1 << 8) - 1
VOC = NOX << 8


def load(tag, kind):
    v = np.zeros(1 << N)
    if kind == "city":
        for k in range(1, 13):
            p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (tag, k))
            if not os.path.exists(p):
                continue
            for row in csv.DictReader(open(p)):
                val = row.get("v_mda8", "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    else:
        for row in csv.DictReader(open(os.path.join(W, "n16_v_%s_m8.csv" % tag))):
            val = row.get("v_mda8", "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


L = []
add = L.append
add("Sector-level sufficiency check (MDA8, ppb), canonical n=16 tables")
add("block = all eight windows of one sector switched off together")
add("")
add("%-14s %9s %9s %9s %9s %10s %9s" %
    ("game", "NOx block", "VOC block", "sum", "joint", "interaction", "singles"))
for g, kind in GAMES.items():
    v = load(g, kind)
    v0 = v[0]
    jn = v[NOX] - v0
    jv = v[VOC] - v0
    jb = v[NOX | VOC] - v0
    sing = sum(v[1 << i] - v0 for i in range(N))
    sing_n = sum(v[1 << i] - v0 for i in range(8))
    sing_v = sum(v[1 << (8 + i)] - v0 for i in range(8))
    add("%-14s %+9.3f %+9.3f %+9.3f %+9.3f %+10.3f %+9.3f" %
        (g, jn, jv, jn + jv, jb, jb - jn - jv, sing))
    add("   single sums: NOx %+.3f   VOC %+.3f   all %+.3f   (joint minus all-singles %+.3f)"
        % (sing_n, sing_v, sing, jb - sing))
add("")
txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_sector_blocks.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_sector_blocks.txt")
