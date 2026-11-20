#!/usr/bin/env python3
"""Taylor-2 (HDDM-2) vs DDM-1 scored on EVERY coalition of the CITY-ANCHORED
n=16 games.  Identical protocol to _hddm_coalitions.py (synthetic) -- same
delta, same coefficient set, same flip floor -- so the two tables are directly
comparable.

Model
  Switching off S means x_i: 1 -> 0, so with D1_i = df/dx_i and D2_ij =
  d2f/dx_i dx_j both at the all-on state,
      DDM-1  : Delta_hat(S) = -sum_{i in S} D1_i
      HDDM-2 : Delta_hat(S) = -sum_{i in S} D1_i + 1/2 sum_{i,j in S} D2_ij
  Truth is v(S) - v(0) from the exact enumeration.  152 coefficients against
  65535 coalitions, no combinatorial sampling on either side.

Inputs
  hddm_v_<label>_*.csv     the 512 derivative runs (fractional multipliers)
  n16_v_<label>_m8_*.csv   the exact truth table, 12 shards (all 2^16 masks)

Writes results_hddm_cities.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_hddm_cities.txt")
DELTA = 0.05
NPLAY = 16
LABELS = ("dallas_2016w", "deerpark_2017s")
METRICS = (("v_mda8", "MDA8 (primary)"), ("v_ppb", "d3_mean (secondary)"))
FLIP_FLOOR = 0.01   # ignore signs below 1% of max|truth| when counting flips


def load_frac(label):
    """the 512 derivative runs -> {name: (v48, mda8)}"""
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_h_v_%s_*.csv" % label))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def load_truth(label):
    """(raw, base): raw[mask] = absolute (v48, mda8); base = value at mask 0.

    NOTE the derivatives are taken w.r.t. ABSOLUTE concentrations, so the
    central-difference for D2_ii needs the absolute base value.
    """
    raw = {}
    for K in range(1, 13):
        p = os.path.join(HERE, "n16_vh_%s_m8_%d.csv" % (label, K))
        if not os.path.exists(p):
            continue
        for row in csv.DictReader(open(p)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            raw[int(row["name"][1:])] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return raw, raw[0]


def derivatives(frac, cen, j):
    """D1 (16,), D2 (16,16) for metric index j; cen = absolute baseline value."""
    d1 = np.zeros(NPLAY)
    d2 = np.zeros((NPLAY, NPLAY))
    for i in range(NPLAY):
        vp = frac["p%02d" % i][j]
        vm = frac["m%02d" % i][j]
        d1[i] = (vp - vm) / (2 * DELTA)
        d2[i, i] = (vp - 2 * cen + vm) / (DELTA ** 2)
    for i in range(NPLAY):
        for k in range(i + 1, NPLAY):
            vpp = frac["c%02d_%02d_pp" % (i, k)][j]
            vpm = frac["c%02d_%02d_pm" % (i, k)][j]
            vmp = frac["c%02d_%02d_mp" % (i, k)][j]
            vmm = frac["c%02d_%02d_mm" % (i, k)][j]
            v = (vpp - vpm - vmp + vmm) / (4 * DELTA ** 2)
            d2[i, k] = d2[k, i] = v
    return d1, d2


def main():
    L = []
    add = L.append
    add("Taylor-2 (HDDM-2) vs DDM-1 on ALL 2^16-1 = 65535 coalitions")
    add("CITY-ANCHORED slices (Dallas 2016 winter, Deer Park 2017 summer)")
    add("coefficients: n + n(n+1)/2 = 152   (from 512 derivative runs per slice)")
    add("error normalised by max|truth| of that slice+metric")
    add("sign flips counted only when both |pred| and |truth| > 1% of max|truth|")
    add("")

    summary = []
    stats = []
    for label in LABELS:
        frac = load_frac(label)
        raw, base_abs = load_truth(label)
        truth = raw
        if len(frac) != 512 or len(raw) != 1 << NPLAY:
            add("!! %s incomplete (frac %d/512, truth %d/65536) -- skipped"
                % (label, len(frac), len(raw)))
            continue
        masks = np.arange(1, 1 << NPLAY)
        M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)
        ks = M.sum(axis=1).astype(int)

        for col, mname in METRICS:
            j = 1 if col == "v_mda8" else 0
            b = base_abs[j]
            tp = np.array([truth[m][j] - b for m in masks])
            cen = b
            d1, d2 = derivatives(frac, cen, j)
            p1 = -(M @ d1)
            p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
            scale = np.abs(tp).max() or 1.0
            e1 = np.abs(p1 - tp) / scale
            e2 = np.abs(p2 - tp) / scale
            fl1 = (np.sign(p1) != np.sign(tp)) & (np.abs(p1) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)
            fl2 = (np.sign(p2) != np.sign(tp)) & (np.abs(p2) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)

            add("=" * 92)
            add("%s   METRIC = %s   max|truth| = %.3f" % (label, mname, scale))
            add("=" * 92)
            add("%4s %8s | %9s %9s | %9s %9s | %8s %8s" %
                ("k", "n_coal", "DDM1 med", "DDM1 max", "HD2 med", "HD2 max",
                 "flip1%", "flip2%"))
            for k in range(1, NPLAY + 1):
                s = ks == k
                if not s.any():
                    continue
                add("%4d %8d | %8.1f%% %8.1f%% | %8.1f%% %8.1f%% | %7.2f%% %7.2f%%" %
                    (k, int(s.sum()),
                     100 * np.median(e1[s]), 100 * e1[s].max(),
                     100 * np.median(e2[s]), 100 * e2[s].max(),
                     100 * fl1[s].mean(), 100 * fl2[s].mean()))
            add("%4s %8d | %8.1f%% %8.1f%% | %8.1f%% %8.1f%% | %7.2f%% %7.2f%%" %
                ("all", len(masks),
                 100 * np.median(e1), 100 * e1.max(),
                 100 * np.median(e2), 100 * e2.max(),
                 100 * fl1.mean(), 100 * fl2.mean()))
            add("")
            summary.append((label, mname, 100 * np.median(e1), 100 * np.median(e2),
                            100 * fl1.mean(), 100 * fl2.mean(),
                            int((e2 > 0.5).sum()),
                            int(fl1.sum()), int(fl2.sum())))
            stats.append((label, mname, "DDM-1", e1))
            stats.append((label, mname, "HDDM-2", e2))

    add("=" * 92)
    add("HEADLINE")
    add("=" * 92)
    add("%-18s %-18s %10s %10s %9s %9s %12s %11s %11s" %
        ("slice", "metric", "DDM1 med", "HD2 med", "flip1", "flip2", "HD2 >50%err",
         "flip1 n", "flip2 n"))
    for r in summary:
        add("%-18s %-18s %9.1f%% %9.1f%% %8.2f%% %8.2f%% %11d/65535 %8d %8d" % r)

    add("")
    add("=" * 92)
    add("DISTRIBUTION OF THE COALITION ERROR (percent of max|truth|)")
    add("=" * 92)
    add("%-18s %-18s %-8s %10s %10s %10s %10s" %
        ("slice", "metric", "method", "mean", "median", "p90", "max"))
    for label, mname, method, e in stats:
        add("%-18s %-18s %-8s %9.1f%% %9.1f%% %9.1f%% %9.1f%%" %
            (label, mname, method, 100 * e.mean(), 100 * np.median(e),
             100 * np.percentile(e, 90), 100 * e.max()))

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
