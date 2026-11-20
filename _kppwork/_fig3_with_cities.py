#!/usr/bin/env python3
"""Fig spectrum: add the two city-anchored curves to the synthetic pair.

4 curves per scenario panel is too busy; instead 2x2: rows = scenario
(r_sat, r_lim1 top) / (Dallas, Deer Park bottom) with the n=8/12/16 ladder
only where defined (cities are n=16 only).  Actually keep the original
2-panel ladder figure and add ONE overlay panel set: cities are n=16-only,
so plot r_sat/r_lim1 n=16 thin + Dallas/Deer Park n=16 bold in a single
second row.  Output: figs/fig3_spectrum.png (regenerated) and
figs/fig3b_cities.png (city panel).
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
WLAB = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24"]


def mobius_np(vals, n):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def mass_curve(m, n):
    idx = np.arange(1, 1 << n)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = np.abs(m[idx])
    mass = np.zeros(n + 1)
    for ki in range(1, n + 1):
        mass[ki] = w[k == ki].sum()
    return mass / mass[1]


def read_n16(label, col="v_mda8"):
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


# city panel: n=16 spectra, synthetic vs anchored
curves = {}
for lab in ["r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"]:
    v = read_n16(lab)
    m = mobius_np(v, N)
    curves[lab] = mass_curve(m, N)

fig, ax = plt.subplots(figsize=(6.2, 4.2))
ks = np.arange(1, N + 1)
styles = {
    "r_sat": ("#7f8c8d", "--", "r_sat"),
    "r_lim1": ("#bdc3c7", "--", "r_lim1"),
    "dallas_2016w": ("#c0392b", "-", "Dallas 2016 winter"),
    "deerpark_2017s": ("#2471a3", "-", "Deer Park 2017 summer"),
}
for lab, (c, ls, name) in styles.items():
    y = np.asarray(curves[lab][1:], dtype=float)
    y = np.where(y > 0, y, np.nan)          # log axis cannot show exact zeros
    ax.plot(ks, y, color=c, ls=ls, lw=2 if ls == "-" else 1.4,
            marker="o", ms=3.5, label=name)
ax.set_yscale("log")
ax.set_ylim(0.05, 500)
ax.set_xlabel("interaction order $k$", fontsize=11)
ax.set_ylabel("$M(k)\\,/\\,M(1)$  (MDA8, log scale)", fontsize=11)
ax.set_xticks(ks[::1])
ax.tick_params(labelsize=9)
ax.grid(alpha=0.35, lw=0.5)
ax.legend(fontsize=8.5, frameon=False)
ax.set_title("Spectral mass by order: synthetic vs anchored slices",
             fontsize=10.5)
fig.tight_layout()
out = os.path.join(ROOT, "figs", "fig3b_cities.png")
fig.savefig(out, dpi=200)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
for lab in curves:
    print(lab, "peak k =", int(np.argmax(curves[lab][1:])) + 1,
          "  R2 =", round(curves[lab][1:3].sum() / curves[lab][1:].sum(), 3))
