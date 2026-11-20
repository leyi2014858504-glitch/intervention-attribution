#!/usr/bin/env python3
"""Exact budget-constrained control plan from the full n=16 table.

Budget currency: fraction of the daily anthropogenic emission of each precursor,
i.e. the two EKMA axes (R_N, R_V).  A coalition is feasible if the NOx windows
switched off remove at most R_N and the VOC windows at most R_V.

For every feasible coalition the MDA8 is known exactly, so the optimal plan is
an exhaustive scan, not a heuristic.  A greedy rule ranked by single-window
effect per unit emission is computed for comparison.

Writes results_optimal_plan.txt, results_optimal_plan.csv, figs/fig_optimal.*
"""
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

D = os.environ.get("ATMOS_DIR", "/mnt/d/atmos")
FIG = os.path.join(D, "figs")
N = 16
IDX = np.arange(1 << N)
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
TITLE = {"r_sat": "r_sat", "r_lim1": "r_lim1",
         "dallas_2016w": "Dallas", "deerpark_2017s": "Deer Park"}
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")
BUDGETS = ((0.10, 0.35), (0.20, 0.50))      # (R_N, R_V)


def gp(mu, sd, h):
    return math.exp(-(((h - mu) / sd) ** 2))


def profile(h):
    return 0.45 + 1.30 * (gp(8.0, 1.6, h) + gp(18.5, 2.0, h))


SHARE = np.array([sum(profile(w * 3 + m / 60.0) for m in range(180)) / 60.0
                  for w in range(8)])
SHARE = SHARE / SHARE.sum()


def load(tag):
    suf = ""
    v = np.full(1 << N, np.nan)
    for k in range(1, 13):
        p = os.path.join(D, "_kppwork", "n16_vh_%s_m8%s_%d.csv" % (tag, suf, k))
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get("v_mda8", "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    if np.isnan(v).any():
        raise SystemExit("%s missing rows" % tag)
    return v


# per-mask removed fractions
RN = np.zeros(1 << N)
RV = np.zeros(1 << N)
for w in range(8):
    RN[IDX & (1 << w) > 0] += SHARE[w]
    RV[IDX & (1 << (8 + w)) > 0] += SHARE[w]


def label_set(mask):
    off = [("NOx " + WIN[w]) for w in range(8) if mask >> w & 1]
    off += [("VOC " + WIN[w]) for w in range(8) if mask >> (8 + w) & 1]
    return ", ".join(off) if off else "none"


L = []
add = L.append
add("Exact budget-constrained MDA8 minimisation, n=16 table")
add("budget = fraction of daily anthropogenic emission removed, per precursor")
add("window emission shares (weekday profile): " +
    "  ".join("%s %.1f%%" % (WIN[w], 100 * SHARE[w]) for w in range(8)))
add("")

rows = []
for g in GAMES:
    v = load(g)
    base = v[0]
    for rn_max, rv_max in BUDGETS:
        feas = (RN <= rn_max + 1e-12) & (RV <= rv_max + 1e-12)
        costs = np.where(feas, v, np.inf)
        best = int(np.argmin(costs))
        d = v[best] - base
        rows.append((g, rn_max, rv_max, best, RN[best], RV[best], v[best], d))

        # greedy: single-window effect per unit emission, then budget
        eff = np.array([(v[1 << i] - base) / SHARE[i % 8] for i in range(N)])
        order = np.argsort(eff)
        gm = 0
        for i in order:
            trial = gm | (1 << i)
            if RN[trial] <= rn_max + 1e-12 and RV[trial] <= rv_max + 1e-12:
                gm = trial
        rows[-1] = rows[-1] + (gm, v[gm] - base)

for g in GAMES:
    add("=" * 74)
    add(TITLE[g])
    for (gg, rn_max, rv_max, best, un, uv, vb, db, gm, dg) in rows:
        if gg != g:
            continue
        add("  budget NOx<=%.0f%%, VOC<=%.0f%%" % (100 * rn_max, 100 * rv_max))
        add("    exact : %-52s  R_N %.1f%%  R_V %.1f%%  dMDA8 %+.3f ppb"
            % (label_set(best), 100 * un, 100 * uv, db))
        add("    greedy: %-52s                dMDA8 %+.3f ppb  (gap %+.3f)"
            % (label_set(gm), dg, dg - db))

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_optimal_plan.txt"), "w") as fh:
    fh.write(txt + "\n")

with open(os.path.join(D, "results_optimal_plan.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["game", "R_N", "R_V", "mask", "used_RN", "used_RV",
                 "mda8", "dMDA8", "greedy_mask", "greedy_dMDA8", "off"])
    for r in rows:
        wr.writerow([r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7],
                     r[8], r[9], label_set(r[3])])

# ---------------- figure: 2 budgets x 4 games, 2x8 window grid
fig, axes = plt.subplots(2, 4, figsize=(8.4, 4.4))
for j, g in enumerate(GAMES):
    for i, (rn_max, rv_max) in enumerate(BUDGETS):
        ax = axes[i, j]
        rec = [r for r in rows if r[0] == g and r[1] == rn_max][0]
        mask = rec[3]
        for sec in range(2):
            for w in range(8):
                bit = sec * 8 + w
                off = (mask >> bit) & 1
                ax.add_patch(plt.Rectangle((w, 1 - sec), 1, 1,
                                           facecolor="#c44e52" if off else "#e8e8e8",
                                           edgecolor="white", lw=1.2))
        ax.set_xlim(0, 8)
        ax.set_ylim(0, 2)
        ax.set_yticks([0.5, 1.5])
        ax.set_yticklabels(["VOC", "NO$_x$"], fontsize=6.5)
        ax.set_xticks(np.arange(8) + 0.5)
        ax.set_xticklabels(WIN, fontsize=5.6, rotation=55, ha="right")
        ax.set_title("%s$\\quad$NO$_x\\leq$%.0f%%, VOC$\\leq$%.0f%%\n"
                     "$\\Delta$MDA8 %+.2f ppb"
                     % (TITLE[g], 100 * rn_max, 100 * rv_max, rec[7]),
                     fontsize=6.8, loc="left")
        ax.tick_params(length=0)
fig.legend(handles=[Patch(facecolor="#c44e52", label="switched off"),
                    Patch(facecolor="#e8e8e8", label="left on")],
           loc="lower center", ncol=2, fontsize=7, frameon=False,
           bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.05, 1, 1))
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, "fig_optimal." + ext), dpi=300, bbox_inches="tight")
print("wrote results_optimal_plan.txt/.csv and figs/fig_optimal.png/.pdf")
