#!/usr/bin/env python3
"""Single-window effects versus the joint effect.

If the 16 single-window effects were a decomposition, their sum would equal
v(all off) - v(all on).  It does not, and the gap measures the interaction that
the allocation has to redistribute.

Writes results_single_vs_joint.txt
"""
import csv
import os

import numpy as np

D = os.environ.get("ATMOS_DIR", "/mnt/d/atmos")
N = 16
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")


def load(tag):
    v = np.full(1 << N, np.nan)
    for k in range(1, 13):
        p = os.path.join(D, "_kppwork", "n16_v_%s_m8_v2_%d.csv" % (tag, k))
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get("v_mda8", "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    return v


L = []
add = L.append
add("Single-window effects vs the joint effect (MDA8, ppb)")
add("%-16s %10s %10s %10s %10s" %
    ("game", "sum single", "joint", "gap", "gap/joint"))
for g in GAMES:
    v = load(g)
    base = v[0]
    single = np.array([v[1 << i] - base for i in range(N)])
    joint = v[(1 << N) - 1] - base
    gap = joint - single.sum()
    add("%-16s %+10.3f %+10.3f %+10.3f %9.0f%%" %
        (g, single.sum(), joint, gap, 100 * gap / joint))
add("")
add("Per-window single effects (ppb), NOx then VOC")
add("%-16s %s" % ("game", " ".join("%7s" % w for w in WIN)))
for g in GAMES:
    v = load(g)
    base = v[0]
    nx = [v[1 << i] - base for i in range(8)]
    vc = [v[1 << (8 + i)] - base for i in range(8)]
    add("%-16s %s" % (g + " NOx", " ".join("%+7.2f" % x for x in nx)))
    add("%-16s %s" % (g + " VOC", " ".join("%+7.2f" % x for x in vc)))

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_single_vs_joint.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_single_vs_joint.txt")
