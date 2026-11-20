#!/usr/bin/env python3
"""Rebuild manuscript Table `error-structure` (Dallas MDA8) on the refit anchor.

(a) sign-error rate by coalition membership of the dominant window (NOx 12-15)
    and the morning NOx window (6-9);
(b) median absolute DDM-1 / HDDM-2 error (ppb) against the number of NOx
    windows off, VOC pooled, plus the span produced by varying the VOC count
    from 0 to 8 at fixed NOx count.

Same 512 fractional runs / delta / flip floor as _hddm_cities.py.
Writes results_error_structure.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_error_structure.txt")
DELTA = 0.05
NPLAY = 16
LABEL = "dallas_2016w"
COL = "v_mda8"
FLIP_FLOOR = 0.01
DOM, MOR = 4, 2          # NOx 12-15 and NOx 6-9


def load_frac():
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_h_v_%s_*.csv" % LABEL))):
        for row in csv.DictReader(open(f)):
            if row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = float(row["v_mda8"])
    return d


def load_truth():
    raw = {}
    for K in range(1, 13):
        p = os.path.join(HERE, "n16_vh_%s_m8_%d.csv" % (LABEL, K))
        with open(p) as fh:
            for row in csv.DictReader(fh):
                if row["v_mda8"] in ("", "NaN"):
                    continue
                raw[int(row["name"][1:])] = float(row["v_mda8"])
    return raw


def derivatives(frac, cen):
    d1 = np.zeros(NPLAY)
    d2 = np.zeros((NPLAY, NPLAY))
    for i in range(NPLAY):
        vp, vm = frac["p%02d" % i], frac["m%02d" % i]
        d1[i] = (vp - vm) / (2 * DELTA)
        d2[i, i] = (vp - 2 * cen + vm) / (DELTA ** 2)
    for i in range(NPLAY):
        for k in range(i + 1, NPLAY):
            d2[i, k] = d2[k, i] = (
                frac["c%02d_%02d_pp" % (i, k)] - frac["c%02d_%02d_pm" % (i, k)]
                - frac["c%02d_%02d_mp" % (i, k)] + frac["c%02d_%02d_mm" % (i, k)]
            ) / (4 * DELTA ** 2)
    return d1, d2


def main():
    frac, raw = load_frac(), load_truth()
    b = raw[0]
    masks = np.arange(1, 1 << NPLAY)
    M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)
    tp = np.array([raw[m] - b for m in masks])
    d1, d2 = derivatives(frac, b)
    p1 = -(M @ d1)
    p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
    scale = np.abs(tp).max()
    big = np.abs(tp) > FLIP_FLOOR * scale
    fl1 = (np.sign(p1) != np.sign(tp)) & big & (np.abs(p1) > FLIP_FLOOR * scale)
    fl2 = (np.sign(p2) != np.sign(tp)) & big & (np.abs(p2) > FLIP_FLOOR * scale)
    e1 = np.abs(p1 - tp)
    e2 = np.abs(p2 - tp)

    has_dom = M[:, DOM] > 0
    has_mor = M[:, MOR] > 0
    cls = [
        ("contains dominant window", has_dom),
        ("omits dominant, contains morning NOx", (~has_dom) & has_mor),
        ("omits dominant, omits morning NOx", (~has_dom) & (~has_mor)),
    ]
    L = []
    add = L.append
    add("Error structure, Dallas %s, MDA8, refit anchor (max|truth| = %.3f ppb)" % (LABEL, scale))
    add("")
    add("(a) sign error by coalition class")
    add("  %-42s %8s %9s %9s" % ("class", "n", "DDM-1", "HDDM-2"))
    for name, flag in cls:
        add("  %-42s %8d %8.2f%% %8.2f%%" %
            (name, int(flag.sum()), 100 * fl1[flag].mean(), 100 * fl2[flag].mean()))

    add("")
    add("(b) median absolute error (ppb) by NOx windows off, VOC pooled")
    nnox = M[:, :8].sum(axis=1).astype(int)
    nvoc = M[:, 8:].sum(axis=1).astype(int)
    add("  %-8s" % "metric" + "".join("%8d" % k for k in range(1, 9)) + "%10s" % "span")
    for name, e in (("DDM-1", e1), ("HDDM-2", e2)):
        med = [np.median(e[nnox == k]) for k in range(1, 9)]
        # same definition as the manuscript column: median error pooled over all
        # NOx counts at VOC=8 divided by the same at VOC=0
        span = np.median(e[nvoc == 8]) / np.median(e[nvoc == 0])
        add("  %-8s" % name + "".join("%8.2f" % v for v in med)
            + "%10.2fx" % span)

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
