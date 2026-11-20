#!/usr/bin/env python3
"""MLMC two-level test for BOTH metrics on the same coarse/fine probe tables.

fine = DT_S 300 s, coarse = DT_S 1800 s, same 1,019 coalitions per city.
Reports the level-difference variance ratio and the speed-up under
  (a) step-proportional cost  (c1 = 300/1800, a coupled driver or a CTM)
  (b) measured wall-clock cost (c1 = 0.25/0.35, the current sweep harness)

Writes results_mlmc.txt
"""
import csv
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results_mlmc.txt")
LABELS = ("dallas_2016w", "deerpark_2017s")
METRICS = (("v_mda8", "MDA8"), ("v_ppb", "d3_mean"))
C_STEP = 300.0 / 1800.0
C_WALL = 0.25 / 0.35


def load(pattern, col):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        for row in csv.DictReader(open(f)):
            try:
                d[int(row["name"][1:])] = float(row[col])
            except (ValueError, KeyError):
                continue
    return d


def steps(label, fine, coarse):
    ch = json.load(open(os.path.join(HERE, "probe_chains_%s.json" % label)))["chains"]
    df, dc = [], []
    for c in ch:
        ms = c["masks"]
        for k in range(1, len(ms)):
            prev = ms[k - 1]
            if prev in fine and ms[k] in fine and prev in coarse and ms[k] in coarse:
                df.append(fine[ms[k]] - fine[prev])
                dc.append(coarse[ms[k]] - coarse[prev])
    return np.array(df), np.array(dc)


def main():
    L, add = [], None
    add = L.append
    add("MLMC two-level test, both metrics, city n=16 games")
    add("fine DT_S 300 s / coarse DT_S 1800 s, 1,019 coalitions per city")
    add("")
    for label in LABELS:
        for col, mname in METRICS:
            fine = load("probe_probe_%s_*.csv" % label, col)
            coarse = load("probe_probedt1800_%s_*.csv" % label, col)
            df, dc = steps(label, fine, coarse)
            Vf = float(np.var(df))
            V1 = float(np.var(dc))
            V2 = float(np.var(df - dc))
            if not (Vf > 0):
                continue
            sp_step = Vf / (np.sqrt(V1 * C_STEP) + np.sqrt(V2)) ** 2
            sp_wall = Vf / (np.sqrt(V1 * C_WALL) + np.sqrt(V2)) ** 2
            add("=" * 76)
            add("%s   metric = %s   steps %d" % (label, mname, len(df)))
            add("=" * 76)
            add("  Var(fine) %.5f   Var(coarse) %.5f   Var(fine-coarse) %.5f  (%.2f%% of fine)"
                % (Vf, V1, V2, 100 * V2 / Vf))
            add("  MLMC speed-up, step-scaled cost : %.2fx" % sp_step)
            add("  MLMC speed-up, wall-clock cost  : %.2fx" % sp_wall)
            add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
