#!/usr/bin/env python3
"""Three-kernel allocation figure: Shapley bars + the two Owen values as markers.

Reads results_three_kernels_players.csv (16 players x 4 games, MDA8).
Writes figs/fig_alloc.png and figs/fig_alloc.pdf
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

D = os.environ.get("ATMOS_DIR", "/mnt/d/atmos")
FIG = os.path.join(D, "figs")
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
TITLE = {"r_sat": "r_sat", "r_lim1": "r_lim1",
         "dallas_2016w": "Dallas", "deerpark_2017s": "Deer Park"}
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")

rows = {}
with open(os.path.join(D, "results_three_kernels_players.csv"), newline="") as fh:
    for r in csv.DictReader(fh):
        if r["metric"] != "MDA8" or r["game"] not in GAMES:
            continue
        rows.setdefault(r["game"], []).append(
            (int(r["player_index"]),
             float(r["shapley_ppb"]), float(r["owen_time_ppb"]),
             float(r["owen_species_ppb"])))

fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.4))
for ax, g in zip(axes.ravel(), GAMES):
    data = sorted(rows[g])
    x = np.arange(16)
    phi = np.array([d[1] for d in data])
    psi_t = np.array([d[2] for d in data])
    psi_s = np.array([d[3] for d in data])
    cols = ["#4c72b0"] * 8 + ["#dd8452"] * 8
    ax.bar(x, phi, color=cols, alpha=0.60, width=0.72, zorder=2)
    lo = np.minimum(np.minimum(phi, psi_t), psi_s)
    hi = np.maximum(np.maximum(phi, psi_t), psi_s)
    ax.vlines(x, lo, hi, color="0.45", lw=1.0, zorder=2.5)
    ax.plot(x, psi_t, "o", ms=5.0, mfc="white", mec="k", mew=1.0, zorder=3)
    ax.plot(x, psi_s, "^", ms=5.4, mfc="none", mec="#c44e52", mew=1.0, zorder=3)
    ax.axhline(0, color="0.35", lw=0.7, zorder=1)
    ax.axvline(7.5, color="0.75", lw=0.8, ls="--", zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels(WIN * 2, fontsize=6.3, rotation=55, ha="right")
    lim = max(np.abs(phi).max(), np.abs(psi_t).max(), np.abs(psi_s).max()) * 1.18
    ax.set_ylim(-lim, lim)
    ax.text(3.5, -lim * 0.96, "NO$_x$ windows", ha="center", fontsize=7.5)
    ax.text(11.5, -lim * 0.96, "VOC windows", ha="center", fontsize=7.5)
    joint = phi.sum()
    ax.set_title("%s   joint effect %+.2f ppb" % (TITLE[g], joint),
                 loc="left", fontsize=8.5)
    ax.tick_params(axis="y", labelsize=6.5)
    ax.tick_params(axis="x", length=2)

for ax in axes[:, 0]:
    ax.set_ylabel("allocated share (ppb)", fontsize=8)
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
handles = [Patch(facecolor="#4c72b0", alpha=0.6, label="Shapley, NO$_x$ block"),
           Patch(facecolor="#dd8452", alpha=0.6, label="Shapley, VOC block"),
           Line2D([], [], marker="o", ls="", mfc="white", mec="k", label="time-first Owen"),
           Line2D([], [], marker="^", ls="", mfc="none", mec="#c44e52", label="species-first Owen")]
fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.0,
           frameon=False, bbox_to_anchor=(0.5, -0.005))
fig.tight_layout(rect=(0, 0.04, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, "fig_alloc." + ext), dpi=300, bbox_inches="tight")
print("wrote figs/fig_alloc.png/.pdf")
