#!/usr/bin/env python3
"""Score DDM-1 and Taylor-2 on the Deer Park truth sample.

Inputs (dp_v_*.csv shards): base, p/m single-axis, c cross-corner runs
(coefficients) and t<mask> truth runs (150 stratified coalitions).

Models (same as _hddm_coalitions.py):
    DDM-1 : Delta_hat(S) = -sum_{i in S} D1_i
    HD-2  : Delta_hat(S) = -sum D1_i + 1/2 sum_{i,j in S} D2_ij
Truth : v(t<mask>) - v(base), both metrics.

Reports per k-group (small <=4, mid 5-9, large >=10) and overall:
median/max |error| relative to max|truth|, sign-flip rate (floor 1% of
max|truth|), for both metrics.

Writes results_deerpark_taylor.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_deerpark_taylor.txt")
DELTA = 0.05
NPLAY = 16
FLIP_FLOOR = 0.01


def load():
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "dp_v_*.csv"))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def derivatives(d, j):
    cen = d["base"][j]
    d1 = np.zeros(NPLAY)
    d2 = np.zeros((NPLAY, NPLAY))
    for i in range(NPLAY):
        vp = d["p%02d" % i][j]
        vm = d["m%02d" % i][j]
        d1[i] = (vp - vm) / (2 * DELTA)
        d2[i, i] = (vp - 2 * cen + vm) / (DELTA ** 2)
    for i in range(NPLAY):
        for k in range(i + 1, NPLAY):
            v = (d["c%02d_%02d_pp" % (i, k)][j] - d["c%02d_%02d_pm" % (i, k)][j]
                 - d["c%02d_%02d_mp" % (i, k)][j] + d["c%02d_%02d_mm" % (i, k)][j]) \
                / (4 * DELTA ** 2)
            d2[i, k] = d2[k, i] = v
    return d1, d2


def main():
    d = load()
    L = []
    add = L.append
    add("Taylor-2 vs DDM-1 at the DEER PARK anchored slice (2017s, LIMITED side)")
    add("150 stratified truth coalitions (small k<=4: 40, mid 5-9: 40, large >=10: 70)")
    add("errors normalised by max|truth|; flips need both sides > 1% of max|truth|")
    add("")
    for col, mname in (("v_mda8", "MDA8 (primary)"), ("v_ppb", "d3_mean (secondary)")):
        j = 1 if col == "v_mda8" else 0
        d1, d2 = derivatives(d, j)
        truth = {}
        for name, v in d.items():
            if name.startswith("t") and name[1:].isdigit():
                truth[int(name[1:])] = v[j] - d["base"][j]
        masks = np.array(sorted(truth))
        tp = np.array([truth[m] for m in masks])
        M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)
        ks = M.sum(axis=1).astype(int)
        scale = np.abs(tp).max() or 1.0
        p1 = -(M @ d1)
        p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
        e1 = np.abs(p1 - tp) / scale
        e2 = np.abs(p2 - tp) / scale
        fl1 = (np.sign(p1) != np.sign(tp)) & (np.abs(p1) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)
        fl2 = (np.sign(p2) != np.sign(tp)) & (np.abs(p2) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)

        add("=" * 84)
        add("METRIC = %s   max|truth| = %.3f ppb   (base %.3f)" %
            (mname, scale, d["base"][j]))
        add("=" * 84)
        add("%-10s %5s | %9s %9s | %9s %9s | %8s %8s" %
            ("k-group", "n", "DDM1 med", "DDM1 max", "HD2 med", "HD2 max",
             "flip1%", "flip2%"))
        for gname, sel in (("small k<=4", ks <= 4), ("mid 5<=k<=9", (ks >= 5) & (ks <= 9)),
                           ("large k>=10", ks >= 10), ("all", ks >= 1)):
            add("%-10s %5d | %8.1f%% %8.1f%% | %8.1f%% %8.1f%% | %7.1f%% %7.1f%%" %
                (gname, int(sel.sum()), 100 * np.median(e1[sel]), 100 * e1[sel].max(),
                 100 * np.median(e2[sel]), 100 * e2[sel].max(),
                 100 * fl1[sel].mean(), 100 * fl2[sel].mean()))
        worst1 = int(np.argmax(e1))
        worst2 = int(np.argmax(e2))
        add("  worst DDM-1 : mask %016x  k=%d  err %.1f%%  (pred %+.2f truth %+.2f)" %
            (masks[worst1], ks[worst1], 100 * e1[worst1], p1[worst1], tp[worst1]))
        add("  worst HD-2  : mask %016x  k=%d  err %.1f%%  (pred %+.2f truth %+.2f)" %
            (masks[worst2], ks[worst2], 100 * e2[worst2], p2[worst2], tp[worst2]))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("wrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
