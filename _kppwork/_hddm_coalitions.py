#!/usr/bin/env python3
"""Taylor-2 (HDDM) vs DDM-1 scored on EVERY coalition of the n=16 game.

Inputs
  hddm_v_<scen>_*.csv   the 512 derivative runs (fractional multipliers)
  n16_v_<scen>_m8.csv   the exact truth table for all 2^16 coalitions

Model
  f(x) = peak-O3 response as a function of the retained fractions x_i.
  Switching off the set S means x_i: 1 -> 0, so with D1_i = df/dx_i at x=1 and
  D2_ij = d2f/dx_i dx_j at x=1,

      DDM-1  :  Delta_hat(S) = -sum_{i in S} D1_i
      HDDM-2 :  Delta_hat(S) = -sum_{i in S} D1_i + 1/2 sum_{i,j in S} D2_ij

  Both are extrapolations from the all-on base; the truth is v(S) - v(0).
  Coefficient count: n + n(n+1)/2 = 16 + 136 = 152, against 65535 coalitions.

Reports, stratified by coalition size k and by scenario (regime), for both
metrics: median and max |error| relative to max|truth|, and the sign-flip rate
(only counted when both prediction and truth exceed 1% of max|truth|, so that
numerical dust is not scored as a flip).

Writes results_hddm_coalitions.txt
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_hddm_coalitions.txt")
DELTA = 0.05
NPLAY = 16
SCENS = ("r_sat", "r_lim1")
METRICS = (("v_mda8", "MDA8 (primary)"), ("v_ppb", "d3_mean (secondary)"))
FLIP_FLOOR = 0.01   # ignore signs below 1% of max|truth| when counting flips


def load_frac(scen):
    """the 512 derivative runs -> {name: (v48, mda8)}"""
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_v_%s_*.csv" % scen))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def load_truth(scen):
    """(raw, base): raw[mask] = absolute (v48, mda8); base = value at mask 0.

    NOTE the derivatives are taken w.r.t. ABSOLUTE concentrations, so the
    central-difference for D2_ii needs the absolute base value -- using the
    baseline-subtracted zero there inflates D2 by a factor of ~2*base/delta^2.
    """
    raw = {}
    with open(os.path.join(HERE, "n16_v_%s_m8.csv" % scen)) as fh:
        for row in csv.DictReader(fh):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            raw[int(row["name"][1:])] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return raw, raw[0]


def derivatives(frac, cen, j):
    """D1 (16,), D2 (16,16) for metric index j; cen = baseline value."""
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
    add("coefficients: n + n(n+1)/2 = 152   (fitted from 512 runs per scenario)")
    add("error normalised by max|truth| of that scenario+metric")
    add("sign flips counted only when both |pred| and |truth| > 1% of max|truth|")
    add("")

    summary = []
    for scen in SCENS:
        frac = load_frac(scen)
        raw, base_abs = load_truth(scen)
        truth = raw
        if len(frac) != 512 or len(raw) != 1 << NPLAY:
            add("!! %s incomplete (frac %d/512, truth %d/65536) -- skipped"
                % (scen, len(frac), len(raw)))
            continue
        masks = np.arange(1, 1 << NPLAY)
        M = ((masks[:, None] >> np.arange(NPLAY)) & 1).astype(float)  # (65535, 16)
        ks = M.sum(axis=1).astype(int)

        for col, mname in METRICS:
            j = 1 if col == "v_mda8" else 0
            b = base_abs[j]
            tp = np.array([truth[m][j] - b for m in masks])   # truth relative to all-on
            cen = b                                          # absolute, for the derivatives
            d1, d2 = derivatives(frac, cen, j)
            p1 = -(M @ d1)
            p2 = p1 + 0.5 * np.einsum("mi,ik,mk->m", M, d2, M)
            scale = np.abs(tp).max() or 1.0
            e1 = np.abs(p1 - tp) / scale
            e2 = np.abs(p2 - tp) / scale
            fl1 = (np.sign(p1) != np.sign(tp)) & (np.abs(p1) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)
            fl2 = (np.sign(p2) != np.sign(tp)) & (np.abs(p2) > FLIP_FLOOR * scale) & (np.abs(tp) > FLIP_FLOOR * scale)

            add("=" * 92)
            add("%s   METRIC = %s   max|truth| = %.3f" % (scen, mname, scale))
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
            summary.append((scen, mname, 100 * np.median(e1), 100 * np.median(e2),
                            100 * fl1.mean(), 100 * fl2.mean(),
                            int((e2 > 0.5).sum())))

    add("=" * 92)
    add("HEADLINE")
    add("=" * 92)
    add("%-8s %-18s %10s %10s %9s %9s %12s" %
        ("scenario", "metric", "DDM1 med", "HD2 med", "flip1", "flip2", "HD2 >50%err"))
    for r in summary:
        add("%-8s %-18s %9.1f%% %9.1f%% %8.2f%% %8.2f%% %11d/65535" % r)

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
