#!/usr/bin/env python3
"""Are DDM-1 / HDDM-2 sign errors concentrated near particular KINDS of
coalitions?  Structural stratification of the city-anchored n=16 audit.

Same protocol and coefficients as _hddm_cities.py (512 fractional-multiplier
runs, delta = 0.05, 152 Taylor coefficients, flip floor = 1% of max|truth|).
Truth is the exact enumeration, 65535 coalitions.

Because broad features correlate with coalition size k (a 12-window coalition
is far more likely to contain the 6-9 h window), a raw rate comparison between
"contains X" and "does not contain X" is confounded by k.  Every stratified rate
below is therefore DIRECTLY STANDARDISED over k:

    rate_F = sum_k  P(k) * rate_F(k),     rate_F(k) = flips{F = 1, size k} / n{F = 1, size k}

with P(k) taken from the full coalition population.  That answers "at the same
intervention size, does containing X change the flip rate?".

Player indices (bit i of the mask):  i = 0..7   NOx windows 0-3,3-6,...,21-24
                                    i = 8..15  VOC windows 0-3,3-6,...,21-24

Writes results_hddm_strata.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_hddm_strata.txt")
DELTA = 0.05
NPLAY = 16
LABELS = ("dallas_2016w", "deerpark_2017s")
METRICS = (("v_mda8", "MDA8"), ("v_ppb", "d3_mean"))
FLIP_FLOOR = 0.01
WINS = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24"]


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


def standardise(flag, flips, ks):
    """k-standardised rate for a binary stratum flag, plus the raw rate."""
    tot = len(ks)
    pk = {k: (ks == k).sum() / tot for k in range(1, NPLAY + 1)}
    num = den = 0.0
    for k in range(1, NPLAY + 1):
        sel = (ks == k) & flag
        n = int(sel.sum())
        if n == 0:
            continue
        num += pk[k] * flips[sel].mean()
        den += pk[k]
    return (num / den if den else float("nan")), flips[flag].mean() if flag.any() else float("nan")


def main():
    L, add = [], None
    add = L.append
    add("Structural stratification of Taylor-2 / DDM-1 sign errors")
    add("CITY-ANCHORED n=16 audit, all 65535 coalitions, flip floor 1% of max|truth|")
    add("Stratified rates are standardised over coalition size k (see header)")
    add("lift = standardised rate with the feature - standardised rate without it")
    add("")

    for label in LABELS:
        frac = load_frac(label)
        raw, base_abs = load_truth(label)
        if len(frac) != 512 or len(raw) != 1 << NPLAY:
            add("!! %s incomplete -- skipped" % label)
            continue
        masks = np.arange(1, 1 << NPLAY)
        M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)
        ks = M.sum(axis=1).astype(int)
        nopt = (masks < 256) | (masks >= 256)          # placeholder, unused
        nox = M[:, :8].sum(axis=1) > 0
        voc = M[:, 8:].sum(axis=1) > 0

        for col, mname in METRICS:
            j = 1 if col == "v_mda8" else 0
            b = base_abs[j]
            tp = np.array([raw[m][j] - b for m in masks])
            d1, d2 = derivatives(frac, b, j)
            p1 = -(M @ d1)
            p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
            scale = np.abs(tp).max() or 1.0
            e1 = np.abs(p1 - tp) / scale
            e2 = np.abs(p2 - tp) / scale
            big = np.abs(tp) > FLIP_FLOOR * scale
            fl1 = (np.sign(p1) != np.sign(tp)) & big & (np.abs(p1) > FLIP_FLOOR * scale)
            fl2 = (np.sign(p2) != np.sign(tp)) & big & (np.abs(p2) > FLIP_FLOOR * scale)

            # dominant single player, from the truth itself
            single = np.array([abs(raw[1 << i][j] - b) for i in range(NPLAY)])
            dom = int(np.argmax(single))

            feats = [
                ("contains NOx 6-9 (bit 2)", M[:, 2] > 0),
                ("contains VOC 6-9 (bit 10)", M[:, 10] > 0),
                ("contains either 6-9 window", (M[:, 2] + M[:, 10]) > 0),
                ("contains dominant player bit %d (%s %s)" % (dom, "NOx" if dom < 8 else "VOC", WINS[dom % 8]), M[:, dom] > 0),
                ("contains any 6-15 h window", (M[:, [2, 3, 4, 5]] + M[:, [10, 11, 12, 13]]).sum(axis=1) > 0),
                ("mixed species (NOx and VOC)", nox & voc),
                ("NOx only", nox & ~voc),
                ("VOC only", voc & ~nox),
                ("contains all four daytime NOx windows", M[:, [2, 3, 4, 5]].sum(axis=1) == 4),
            ]

            add("=" * 96)
            add("%s   METRIC = %s   max|truth| = %.3f   dominant player = %s %s"
                % (label, mname, scale, "NOx" if dom < 8 else "VOC", WINS[dom % 8]))
            add("overall: flip1 = %.2f%%  flip2 = %.2f%%  median err1 = %.1f%%  median err2 = %.1f%%"
                % (100 * fl1.mean(), 100 * fl2.mean(),
                   100 * np.median(e1), 100 * np.median(e2)))
            add("=" * 96)
            add("%-38s %8s %8s %8s | %8s %8s %8s" %
                ("structure", "n_coal", "f1_std", "lift1", "n_coal", "f2_std", "lift2"))
            # baselines
            b1, _ = standardise(np.ones(len(ks), bool), fl1, ks)
            b2, _ = standardise(np.ones(len(ks), bool), fl2, ks)
            for name, flag in feats:
                s1, _ = standardise(flag, fl1, ks)
                s2, _ = standardise(flag, fl2, ks)
                add("%-38s %8d %7.2f%% %+7.2f | %8d %7.2f%% %+7.2f"
                    % (name, int(flag.sum()), 100 * s1, 100 * (s1 - b1),
                       int(flag.sum()), 100 * s2, 100 * (s2 - b2)))
            add("  (baseline standardised: flip1 = %.2f%%, flip2 = %.2f%%)"
                % (100 * b1, 100 * b2))

            # worst coalitions by HDDM-2 error
            order = np.argsort(-e2)[:8]
            add("")
            add("worst 8 coalitions by HDDM-2 relative error:")
            for m in order:
                on = [("NOx " if i < 8 else "VOC ") + WINS[i % 8]
                      for i in range(NPLAY) if M[m, i] > 0]
                add("  err=%5.1f%%  k=%2d  truth=%+7.3f  hd2=%+7.3f  flip=%s   {%s}"
                    % (100 * e2[m], ks[m], tp[m], p2[m], "Y" if fl2[m] else "-",
                       ", ".join(on)))
            add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
