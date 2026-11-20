#!/usr/bin/env python3
"""Is the AMPLITUDE error structured, in the same way the sign error is?

Sign errors live in one 2x2 cell (no dominant window, with the morning NOx
window).  This script asks the companion question for |error|:

  absolute error  |p - truth|            -- follows effect size, mostly trivial
  relative error  |p - truth| / max|truth|
  ratio-to-truth  |p| / |truth|          -- under/over-prediction factor

Tabulated on a 2-D grid of composition (number of NOx players, number of VOC
players), which keeps k = n_NOx + n_VOC explicit.  Same protocol and
coefficients as _hddm_cities.py.

Writes results_hddm_amplitude.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_hddm_amplitude.txt")
DELTA = 0.05
NPLAY = 16
LABELS = ("dallas_2016w", "deerpark_2017s")
FLIP_FLOOR = 0.01


def load_frac(label):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_v_%s_*.csv" % label))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def load_truth(label):
    raw = {}
    for K in range(1, 13):
        p = os.path.join(HERE, "n16_v_%s_m8_%d.csv" % (label, K))
        if not os.path.exists(p):
            continue
        for row in csv.DictReader(open(p)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            raw[int(row["name"][1:])] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return raw, raw[0]


def derivatives(frac, cen, j):
    d1 = np.zeros(NPLAY)
    d2 = np.zeros((NPLAY, NPLAY))
    for i in range(NPLAY):
        vp = frac["p%02d" % i][j]
        vm = frac["m%02d" % i][j]
        d1[i] = (vp - vm) / (2 * DELTA)
        d2[i, i] = (vp - 2 * cen + vm) / (DELTA ** 2)
    for i in range(NPLAY):
        for k in range(i + 1, NPLAY):
            d2[i, k] = d2[k, i] = (
                frac["c%02d_%02d_pp" % (i, k)][j]
                - frac["c%02d_%02d_pm" % (i, k)][j]
                - frac["c%02d_%02d_mp" % (i, k)][j]
                + frac["c%02d_%02d_mm" % (i, k)][j]) / (4 * DELTA ** 2)
    return d1, d2


def main():
    L, add = [], None
    add = L.append
    add("Amplitude error by coalition composition (number of NOx vs VOC players)")
    add("CITY-ANCHORED n=16 games, all 65535 coalitions")
    add("")
    add("r2 = |pred - truth| / max|truth|      (MD = median)")
    add("abs = |pred - truth| in ppb             (MD = median)")
    add("U   = |pred| / |truth|, >1 over-prediction, <1 under-prediction")
    add("      computed only where |truth| > 1% of max|truth|")
    add("")

    for label in LABELS:
        frac = load_frac(label)
        raw, base_abs = load_truth(label)
        if len(frac) != 512 or len(raw) != 1 << NPLAY:
            add("!! %s incomplete -- skipped" % label)
            continue
        masks = np.arange(1, 1 << NPLAY)
        M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)
        n_nox = M[:, :8].sum(axis=1).astype(int)
        n_voc = M[:, 8:].sum(axis=1).astype(int)

        for j, mname in ((1, "MDA8"), (0, "d3_mean")):
            b = base_abs[j]
            tp = np.array([raw[m][j] - b for m in masks])
            d1, d2 = derivatives(frac, b, j)
            p1 = -(M @ d1)
            p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
            scale = np.abs(tp).max() or 1.0
            big = np.abs(tp) > FLIP_FLOOR * scale

            add("=" * 100)
            add("%s   METRIC = %s   max|truth| = %.3f ppb" % (label, mname, scale))
            add("=" * 100)

            for tag, p in (("DDM-1", p1), ("HDDM-2", p2)):
                e = np.abs(p - tp)
                r = e / scale
                add("-- %s --  median r2 = %.1f%%   max r2 = %.1f%%   "
                    "median abs = %.3f ppb   max abs = %.3f ppb"
                    % (tag, 100 * np.median(r), 100 * r.max(),
                       np.median(e), e.max()))
                add("   correlation: |error| vs |truth|  r = %.3f   "
                    "|error| vs n_NOx  r = %.3f   |error| vs k  r = %.3f"
                    % (np.corrcoef(e, np.abs(tp))[0, 1],
                       np.corrcoef(e, n_nox)[0, 1],
                       np.corrcoef(e, n_nox + n_voc)[0, 1]))
                add("   median |error| by number of NOx players (all VOC counts pooled):")
                row_n = []
                for a in range(0, 9):
                    s = n_nox == a
                    row_n.append("%d:%.2f" % (a, np.median(e[s])) if s.any() else "%d:--" % a)
                add("     " + "  ".join(row_n))
                add("   median |error| by number of VOC players (all NOx counts pooled):")
                row_v = []
                for a in range(0, 9):
                    s = n_voc == a
                    row_v.append("%d:%.2f" % (a, np.median(e[s])) if s.any() else "%d:--" % a)
                add("     " + "  ".join(row_v))
                # under/over prediction on the big-effect coalitions
                u = np.abs(p[big]) / np.abs(tp[big])
                add("   |pred|/|truth| on 1%%-floor coalitions: median %.2f   "
                    "p10 %.2f   p90 %.2f   share under-predicting %.1f%%"
                    % (np.median(u), np.percentile(u, 10), np.percentile(u, 90),
                       100 * (u < 1).mean()))
                add("")

            # composition grid for HDDM-2 relative error
            add("HDDM-2 median r2 (%%) on the composition grid  [rows = n_NOx, cols = n_VOC]")
            hdr = "     " + "".join("%7d" % b_ for b_ in range(9))
            add(hdr)
            for a in range(9):
                cells = []
                for bb in range(9):
                    if a + bb == 0:
                        cells.append("      -")
                        continue
                    s = (n_nox == a) & (n_voc == bb)
                    cells.append("%7.1f" % (100 * np.median(np.abs(p2[s] - tp[s]) / scale))
                                 if s.any() else "      .")
                add("%4d " % a + "".join(cells))
            add("")

            add("DDM-1 median r2 (%%) on the composition grid")
            add(hdr)
            for a in range(9):
                cells = []
                for bb in range(9):
                    if a + bb == 0:
                        cells.append("      -")
                        continue
                    s = (n_nox == a) & (n_voc == bb)
                    cells.append("%7.1f" % (100 * np.median(np.abs(p1[s] - tp[s]) / scale))
                                 if s.any() else "      .")
                add("%4d " % a + "".join(cells))
            add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
