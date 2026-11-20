#!/usr/bin/env python3
"""Attribution-side sign flip figure -- NO local method involved.

Bars = finite single-window shutdown effect  d_i = v({i}) - v(empty),
read EXACTLY from the enumerated n=16 MDA8 oracle (all-on baseline).
Two sectors x two scenarios.  The adjacent-window sign flip in NOx
(12-15 h lowers MDA8 while the neighbouring 15-18 h window raises it)
is annotated.  Local derivatives enter only in Fig. 6 (later section).
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
WLAB = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24"]


def read_lvg_co(path):
    """CO single-window finite effects (D_finite column) from the lvg file."""
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
                out.append(float(p[3]))  # D_finite
    return out


def load(scen):
    v = {}
    with open(os.path.join(HERE, "n16_vh_%s_m8.csv" % scen)) as fh:
        for row in csv.DictReader(fh):
            val = row.get("v_mda8", "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


data = {}
for scen in ("r_sat", "r_lim1"):
    v = load(scen)
    data[scen] = {
        "NOx": [v[1 << j] - v[0] for j in range(8)],
        "VOC": [v[1 << (8 + j)] - v[0] for j in range(8)],
    }
    # verification print
    print(scen, "v(empty)=%.3f" % v[0])
    print(scen, "NOx d_i:", " ".join("%+.3f" % x for x in data[scen]["NOx"]))
    print(scen, "VOC d_i:", " ".join("%+.3f" % x for x in data[scen]["VOC"]))

# CO from the lvg perturbation protocol (same quantity, different source)
co_lvg = read_lvg_co(os.path.join(ROOT, "results_local_vs_global.txt"))
print("CO d_i (lvg, r_sat):", " ".join("%+.3f" % x for x in co_lvg))
data["r_sat"]["CO"] = co_lvg
# r_lim1 CO from the dedicated 9-run fill (_co_lvg_lim1.py)
co_l1 = []
with open(os.path.join(ROOT, "results_co_lim1.txt")) as fh:
    for line in fh:
        if line.startswith("CO d_i:"):
            co_l1 = [float(x) for x in line.split(":")[1].split()]
print("CO d_i (lvg, r_lim1):", " ".join("%+.3f" % x for x in co_l1))
data["r_lim1"]["CO"] = co_l1 if co_l1 else None

SECS3 = ("NOx", "VOC", "CO")
fig, axes = plt.subplots(3, 2, figsize=(10.2, 8.0), sharex=True)
for col, scen in enumerate(("r_sat", "r_lim1")):
    for row, sec in enumerate(SECS3):
        ax = axes[row][col]
        eff = data[scen][sec]
        if eff is None:
            ax.set_xlim(-0.5, 7.5)
            ax.text(0.5, 0.5, "CO not evaluated\nfor this scenario",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color="#999")
        else:
            cols = ["#c0392b" if x > 0 else "#2471a3" for x in eff]
            ax.bar(range(8), eff, color=cols, edgecolor="none", width=0.72)
            ax.axhline(0, color="#555", lw=0.8)
        ax.set_xticks(range(8))
        ax.set_xticklabels(WLAB, fontsize=7, rotation=45)
        ax.grid(axis="y", lw=0.4, alpha=0.45)
        if row == 0:
            ax.set_title(("r_sat" if scen == "r_sat"
                          else "r_lim1"), fontsize=10)
        if col == 0:
            ax.set_ylabel("%s shutdown\n$\\Delta$MDA8 (ppb)" % sec, fontsize=9)
        # no sign-flip annotation: the caption states trends only
fig.suptitle("Finite single-window shutdown effects $v(\\{i\\})-v(\\emptyset)$",
             fontsize=10, y=0.985)
fig.tight_layout(rect=[0, 0, 1, 0.965])
out = os.path.join(ROOT, "figs", "fig_signflip_finite.png")
fig.savefig(out, dpi=200)
fig.savefig(out.replace(".png", ".pdf"))
print("wrote", out)
