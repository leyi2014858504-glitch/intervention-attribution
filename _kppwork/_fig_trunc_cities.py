#!/usr/bin/env python3
"""Truncation-error figure: synthetic + anchored slices (n=16).
trunc(q) = max_i |phi_i - phi_i^{<=q}| / max|phi_i|.
2x2: (a) r_sat  (b) r_lim1  (c) Dallas  (d) Deer Park.
Both curves on every panel: MDA8 bold, d3_mean thin.  The thin curve keeps
the same shape at a percent-level amplitude and shows that the high-order
structure is present in the chemistry itself, not introduced by the metric.
Writes figs/fig_trunc_cities.png
"""
import csv
import os
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"

N = 16
PANELS = [
    ("r_sat", "r_sat"),
    ("r_lim1", "r_lim1"),
    ("dallas_2016w", "Dallas 2016 winter"),
    ("deerpark_2017s", "Deer Park 2017 summer"),
]


def mobius_np(vals):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << N)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def trunc_curve(m):
    idx = np.arange(1, 1 << N)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(N)
    for i in range(N):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    scale = np.abs(phi).max() or 1.0
    errs = []
    for q in range(1, N + 1):
        sel = k <= q
        phi_q = np.zeros(N)
        for i in range(N):
            phi_q[i] = w[sel & (((idx >> i) & 1) > 0)].sum()
        errs.append(100 * np.abs(phi - phi_q).max() / scale)
    return errs


def read_vals(label, col):
    if label in ("dallas_2016w", "deerpark_2017s"):
        v = np.zeros(65536)
        for K in range(1, 13):
            p = os.path.join(W, "n16_vh_%s_m8_%d.csv" % (label, K))
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                for row in csv.DictReader(fh):
                    val = row.get(col, "")
                    if val not in ("", "NaN"):
                        v[int(row["name"][1:])] = float(val)
        return v
    v = np.zeros(65536)
    with open(os.path.join(W, "n16_vh_%s_m8.csv" % label)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.0), sharex=True)
qs = np.arange(1, N + 1)
for pi, (label, title) in enumerate(PANELS):
    ax = axes[pi // 2][pi % 2]
    for col, name, c, lw in (("v_mda8", "MDA8", "#c0392b", 2.0),
                             ("v_ppb", "d3_mean", "#7f8c8d", 1.3)):
        v = read_vals(label, col)
        m = mobius_np(v)
        ax.plot(qs, trunc_curve(m), color=c, lw=lw, marker="o", ms=3.2,
                label=name)
    ax.set_yscale("log")
    ax.set_title(title, fontsize=10)
    ax.grid(alpha=0.35, lw=0.5)
    ax.tick_params(labelsize=8.5)
    if pi % 2 == 0:
        ax.set_ylabel("trunc($q$)  (% of max|$\\phi$|)", fontsize=9.5)
    if pi // 2 == 1:
        ax.set_xlabel("truncation order $q$", fontsize=9.5)
    ax.legend(fontsize=8, frameon=False)
fig.suptitle("Truncation error of the order-$q$ attribution",
             fontsize=11, y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
out = os.path.join(ROOT, "figs", "fig_trunc_cities.png")
fig.savefig(out, dpi=200)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
