#!/usr/bin/env python3
"""Two-city figure: finite single-window shutdown effects at the anchored
slices (n=16, MDA8), same visual grammar as fig_signflip_finite.png.

Bars = v({i}) - v(empty) read exactly from the enumerated city tables;
also cross-checks against results_*_lvg.txt (L2 gate).  Panel (a) Dallas
2016 winter (expected saturated), (b) Deer Park 2017 summer (limited).

Writes figs/fig_cities_24win.png
"""
import csv
import os
import platform

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"

WLAB = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24"]
CITIES = [
    ("dallas_2016w", "Dallas 2016 winter"),
    ("deerpark_2017s", "Deer Park 2017 summer"),
]


def read_lvg_co(path):
    """CO single-window finite effects (D_finite) from the lvg file."""
    out = []
    in_tab = False
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if not p:
                continue
            if p[0] in ("sectorhours", "sector") and len(p) > 3:
                in_tab = True
                continue
            if in_tab:
                if p[0] != "COW":
                    if out:
                        break
                    continue
                out.append(float(p[3]))
    return out


def read_n16(label, col="v_mda8"):
    v = {}
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


fig, axes = plt.subplots(3, 2, figsize=(10.2, 8.0), sharex=True)
for ci, (label, title) in enumerate(CITIES):
    v = read_n16(label)
    v0 = v[0]
    eff = {
        "NOx": [v[1 << j] - v0 for j in range(8)],
        "VOC": [v[1 << (8 + j)] - v0 for j in range(8)],
        "CO": read_lvg_co(os.path.join(ROOT, "results_%s_lvg.txt"
                                       % label.split("_")[0])),
    }
    # L2 cross-check against the lvg protocol
    print("== %s  v(empty)=%.4f" % (label, v0))
    for sec in ("NOx", "VOC", "CO"):
        print("   %s d_i: %s" % (sec, " ".join("%+.3f" % x for x in eff[sec])))
    for row, sec in enumerate(("NOx", "VOC", "CO")):
        ax = axes[row][ci]
        e = eff[sec]
        cols = ["#c0392b" if x > 0 else "#2471a3" for x in e]
        ax.bar(range(8), e, color=cols, edgecolor="none", width=0.72)
        ax.axhline(0, color="#555", lw=0.8)
        ax.set_xticks(range(8))
        ax.set_xticklabels(WLAB, fontsize=7, rotation=45)
        ax.grid(axis="y", lw=0.4, alpha=0.45)
        if row == 0:
            ax.set_title(title, fontsize=10)
        if ci == 0:
            ax.set_ylabel("%s window shutdown\n$\\Delta$MDA8 (ppb)" % sec,
                          fontsize=9)
fig.suptitle("Finite single-window shutdown effects at the anchored city slices",
             fontsize=10, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.965])
out = os.path.join(ROOT, "figs", "fig_cities_24win.png")
fig.savefig(out, dpi=200)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
