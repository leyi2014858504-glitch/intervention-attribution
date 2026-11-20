#!/usr/bin/env python3
"""Unweighted (kernel-free) diagnostics on the exact n=16 tables, q = 1..16.

For each of the four games and both metrics:
  * bare order mass M(k)/M(1) = sum_{|A|=k}|m(A)| / sum_{|A|=1}|m(A)|
  * cumulative low-order mass share R2, R5 (sum_{k<=q} M(k) / sum_k M(k))
  * err_v(q) = max_S |v(S) - sum_{A subset S, |A|<=q} m(A)| / max_S |v(S)|

Reads the released oracle tables only (zero chemistry runs).
Writes results_errv_full.txt
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
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
METRICS = (("v_mda8", "MDA8"), ("v_ppb", "d3_mean"))
IDX = np.arange(1 << N)
K = np.array([bin(x).count("1") for x in IDX])


def mobius(vals):
    f = np.array(vals, float)
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        f[sel] -= f[IDX[sel] ^ bit]
    return f


def zeta(f):
    g = f.copy()
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        g[sel] += g[IDX[sel] ^ bit]
    return g


def read_vals(game, col):
    v = np.zeros(1 << N)
    if game in ("r_sat", "r_lim1"):
        files = [os.path.join(W, "n16_vh_%s_m8.csv" % game)]
    else:
        files = [os.path.join(W, "n16_vh_%s_m8_%d.csv" % (game, j)) for j in range(1, 13)]
    for p in files:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    return v


L = []
add = L.append
add("Unweighted diagnostics on the exact n=16 tables, q = 1..16")
add("M(k) = sum_{|A|=k} |m(A)|, normalised by M(1)  (no Shapley weighting)")
add("R_q  = cumulative share of M over orders <= q")
add("err_v(q) = max_S |v(S) - sum_{A subset S, |A|<=q} m(A)| / max_S |v(S)|")
add("")
for game in GAMES:
    for col, tag in METRICS:
        v = read_vals(game, col)
        v = v - v[0]
        m = mobius(v)
        mass = np.array([np.abs(m[K == k]).sum() for k in range(1, N + 1)])
        M1 = mass[0] or 1.0
        ref = np.abs(v).max() or 1.0
        r = [100.0 * mass[:q].sum() / mass.sum() for q in (2, 5)]
        add("=== %-16s %-8s  max|v| = %.3f ppb   M(1) = %.4f ppb ===" %
            (game, tag, ref, M1))
        add("  M(k)/M(1)  " + " ".join("k%02d:%7.2f" % (k, mass[k - 1] / M1)
                                       for k in range(1, N + 1)))
        add("  R2 = %.3f%%   R5 = %.3f%%" % (r[0], r[1]))
        e = []
        for q in range(1, N + 1):
            mq = np.where(K <= q, m, 0.0)
            e.append(100.0 * np.abs(v - zeta(mq)).max() / ref)
        add("  err_v(q)   " + " ".join("q%02d:%7.2f" % (q, x)
                                      for q, x in enumerate(e, 1)))
        add("  err_v peak %.2f%% at q=%d   err_v(q=10) = %.2f%%" %
            (max(e), e.index(max(e)) + 1, e[9]))
        add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_errv_full.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_errv_full.txt")
