#!/usr/bin/env python3
"""Local (DDM-1 / HDDM-2) vs exact finite shutdown, hourly MDA8, two cities.

Exact finite effects come from the hourly n=16 truth tables (single-player
coalitions).  DDM-1 and HDDM-2 come from the 512 central-difference runs with
delta = 0.05, evaluated on the same hourly metric.

Writes figs/fig_cities_lvg.png / .pdf
"""
import csv
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

W = r"D:\atmos\_kppwork"
ROOT = r"D:\atmos"
FIG = os.path.join(ROOT, "figs")
DELTA = 0.05
SECS = ("NOXW", "VOCW")
SECLAB = {"NOXW": r"NO$_x$", "VOCW": "VOC"}
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")
CITIES = (("dallas_2016w", "dallas", "Dallas 2016w"),
          ("deerpark_2017s", "deerpark", "Deer Park 2017s"))


def read_truth(label):
    v = np.zeros(1 << 16)
    p = os.path.join(W, "n16_vh_%s_m8.csv" % label)
    files = [p] if os.path.exists(p) else [
        os.path.join(W, "n16_vh_%s_m8_%d.csv" % (label, k)) for k in range(1, 13)]
    for f in files:
        with open(f) as fh:
            for row in csv.DictReader(fh):
                val = row.get("v_mda8", "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    return v


def read_frac(city):
    d = {}
    for f in sorted(glob.glob(os.path.join(W, "hddm_h_v_%s_*.csv" % city))):
        with open(f) as fh:
            for row in csv.DictReader(fh):
                if row["v_mda8"] in ("", "NaN"):
                    continue
                d[row["name"]] = float(row["v_mda8"])
    return d


fig, axes = plt.subplots(2, 2, figsize=(9.0, 5.6), sharex=True)
for ci, (label, city, title) in enumerate(CITIES):
    v = read_truth(label)
    base = v[0]
    frac = read_frac(city)
    for si, sec in enumerate(SECS):
        ax = axes[ci][si]
        x = np.arange(8)
        ddm, hd2, fin = [], [], []
        for w in range(8):
            idx = w if sec == "NOXW" else 8 + w
            vp, vm = frac["p%02d" % idx], frac["m%02d" % idx]
            d1 = (vp - vm) / (2 * DELTA)
            d2 = (vp - 2 * base + vm) / DELTA ** 2
            ddm.append(-d1)
            hd2.append(-d1 + 0.5 * d2)
            fin.append(v[1 << idx] - base)
        ddm, hd2, fin = map(np.array, (ddm, hd2, fin))
        bw = 0.26
        ax.bar(x - bw, ddm, width=bw, color="none", edgecolor="#7f8c8d",
               hatch="///", label=r"$T_1=-D_{1}$")
        ax.bar(x, hd2, width=bw, color="#f5b041", edgecolor="#b9770e",
               label=r"$T_2$")
        cols = ["#c0392b" if y > 0 else "#2471a3" for y in fin]
        ax.bar(x + bw, fin, width=bw, color=cols, label=r"$\Delta_{\mathrm{off}}$")
        ax.axhline(0, color="#555", lw=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels(WIN, fontsize=6, rotation=45)
        ax.grid(axis="y", lw=0.4, alpha=0.4)
        if ci == 0:
            ax.set_title(SECLAB[sec], fontsize=10)
        if si == 0:
            ax.set_ylabel("%s\nMDA8 change (ppb)" % title, fontsize=8.5)
        if label.startswith("dallas") and sec == "NOXW":
            for j in range(8):
                if abs(fin[j]) > 0.5 and abs(ddm[j]) > 1e-9:
                    r = fin[j] / ddm[j]
                    if abs(r - 1.0) > 0.5:
                        ax.annotate("x%.1f" % abs(r), (j - bw, ddm[j]),
                                    textcoords="offset points", xytext=(0, 4),
                                    fontsize=7, ha="center", color="#111",
                                    fontweight="bold")
fig.suptitle("Local predictions vs the exact finite shutdown at the anchored "
             "city slices", fontsize=11, y=0.99)
h, lab = axes[0][0].get_legend_handles_labels()
fig.legend(h, lab, loc="lower center", ncol=3, fontsize=9, frameon=False,
           bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=[0, 0.035, 1, 0.97])
out = os.path.join(FIG, "fig_cities_lvg.png")
fig.savefig(out, dpi=200, bbox_inches="tight")
fig.savefig(out.replace(".png", ".pdf"), bbox_inches="tight")
print("wrote", out)

# summary numbers for the text
for label, city, title in CITIES:
    v = read_truth(label)
    base = v[0]
    frac = read_frac(city)
    worst = 0.0
    print("\n%s  base %.4f" % (title, base))
    for sec in SECS:
        for w in range(8):
            idx = w if sec == "NOXW" else 8 + w
            vp, vm = frac["p%02d" % idx], frac["m%02d" % idx]
            d1 = (vp - vm) / (2 * DELTA)
            d2 = (vp - 2 * base + vm) / DELTA ** 2
            p1, p2 = -d1, -d1 + 0.5 * d2
            f = v[1 << idx] - base
            worst = max(worst, abs(p1 - f), abs(p2 - f))
            if abs(f) > 0.3:
                print("  %-5s %-6s DDM %+7.3f HD2 %+7.3f fin %+7.3f  fin/DDM %6.2f fin/HD2 %6.2f"
                      % (sec, WIN[w], p1, p2, f,
                         f / p1 if abs(p1) > 1e-9 else float("nan"),
                         f / p2 if abs(p2) > 1e-9 else float("nan")))
    print("  max |prediction - finite| = %.3f ppb" % worst)
