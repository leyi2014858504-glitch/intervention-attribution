#!/usr/bin/env python3
"""Truncation error curve for the city slices (post-processing, zero runs).

trunc(q) = max_i |phi_i - phi_i^{<=q}| / max_i |phi_i|   -- the error of the
ORDER-q TRUNCATED attribution, exactly what low-order attribution implicitly
does.  Computed from the exact n=16 spectra (both metrics, both cities), so it
is the honest "how bad is low-order" number to sit next to the mass share R_q.

Writes results_city_trunc.txt
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
LABELS = ["dallas_2016w", "deerpark_2017s"]


def mobius_np(vals, n):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def shapley_from_mobius(m, n):
    idx = np.arange(1, 1 << n)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(n)
    for i in range(n):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    return phi


def trunc_curve(m, n):
    idx = np.arange(1, 1 << n)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(n)
    phi_q = np.zeros(n)
    for i in range(n):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    scale = np.abs(phi).max() or 1.0
    errs = []
    for q in range(1, n + 1):
        phi_q = np.zeros(n)
        sel = k <= q
        for i in range(n):
            phi_q[i] = w[sel & (((idx >> i) & 1) > 0)].sum()
        errs.append(np.abs(phi - phi_q).max() / scale)
    return errs


def read_n16(label, col):
    v = np.zeros(65536)
    for K in range(1, 13):
        p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (label, K))
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh):
                mm = int(row["name"][1:])
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[mm] = float(val)
    return v


L = []
add = L.append
add("Truncation error of the ORDER-q attribution (exact spectra, city slices)")
add("trunc(q) = max_i |phi_i - phi_i^{<=q}| / max_i |phi_i|, in percent")
add("")

for label in LABELS:
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3_mean")):
        v = read_n16(label, col)
        m = mobius_np(v, 16)
        errs = trunc_curve(m, 16)
        add("=== %s  %s ===" % (label, tag))
        add("  " + "  ".join("q<=%-2d %6.1f%%" % (q, 100 * e) for q, e in enumerate(errs, 1)))
        add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_city_trunc.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_city_trunc.txt")
