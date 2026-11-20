#!/usr/bin/env python3
"""Shapley-Taylor k-truncated values: exact error vs the full-spectrum Shapley.

ST_k (Sundararajan et al. 2020) restores efficiency at order k by dividing the
residual r(k) = sum_{|A|>k} m(A) equally among players:

  phi_ST_i(k) = sum_{A ni i, |A|<=k} m(A)/|A| + r(k)/n

We can price ST exactly because we hold the enumerated spectra.  The point of
the test: does efficiency restoration rescue low-order attribution?  (Prediction
from the sign-alternating high-order family: no -- r(k) nearly cancels, the
uniform share is negligible, and the error is the same tail.)

Writes results_shapley_taylor.txt
"""
import csv
import os
import platform
import numpy as np

if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"

N = 16
LABELS = ["dallas_2016w", "deerpark_2017s", "r_sat", "r_lim1"]


def mobius_np(vals):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << N)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def shapley(m):
    idx = np.arange(1, 1 << N)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(N)
    for i in range(N):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    return phi


def read_vals(label, col):
    if label in ("dallas_2016w", "deerpark_2017s"):
        v = np.zeros(65536)
        for K in range(1, 13):
            p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (label, K))
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                for row in csv.DictReader(fh):
                    val = row.get(col, "")
                    if val not in ("", "NaN"):
                        v[int(row["name"][1:])] = float(val)
        return v
    v = np.zeros(65536)
    with open(os.path.join(W, "n16_v_%s_m8.csv" % label)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


L = []
add = L.append
add("Shapley-Taylor k-truncated attribution: exact error vs full-spectrum Shapley")
add("phi_ST_i(k) = sum_{A ni i, |A|<=k} m(A)/|A| + r(k)/n,  r(k) = sum_{|A|>k} m(A)")
add("err_ST(k) = max_i |phi_i - phi_ST_i(k)| / max_i |phi_i|")
add("")

for label in LABELS:
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3_mean")):
        v = read_vals(label, col)
        m = mobius_np(v)
        phi = shapley(m)
        scale = np.abs(phi).max() or 1.0
        idx = np.arange(1, 1 << N)
        k_all = np.array([bin(x).count("1") for x in idx])
        w = m[idx] / k_all
        add("=== %s  %s ===   max|phi| = %.4f" % (label, tag, scale))
        add("  %5s %9s %12s %9s" % ("k", "err_ST(k)", "r(k)/n / max|phi|", "sum phi_ST"))
        for q in range(1, N + 1):
            sel_lo = k_all <= q
            sel_hi = k_all > q
            r = m[idx[sel_hi]].sum()
            phi_st = np.zeros(N)
            for i in range(N):
                mem = ((idx[sel_lo] >> i) & 1) > 0
                phi_st[i] = w[sel_lo][mem].sum() + r / N
            err = np.abs(phi - phi_st).max() / scale
            share = (r / N) / scale
            add("  %5d %8.1f%% %11.3f%% %12.5f" %
                (q, 100 * err, 100 * share, phi_st.sum()))
        add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_shapley_taylor.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_shapley_taylor.txt")
