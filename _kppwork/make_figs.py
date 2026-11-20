#!/usr/bin/env python3
"""Generate the working-paper figures Fig 2-7 from released result files.

Outputs PNG (+PDF) into <root>/figs/.  All values are parsed from released
artifacts -- this script recomputes nothing chemical:
  Fig 2  _kppwork/spectrum_orders{,_m8}.txt        (truncation + order mass)
  Fig 3  _kppwork/co_v.csv, n16_v_r_*_m8.csv       (order masses, fast Mobius)
  Fig 4  _kppwork/scaling_dual_metrics.txt         (cost law, both metrics)
  Fig 5  results_local_vs_global.txt               (finite chords, MDA8 3h)
  Fig 6  results_local_vs_global.txt               (local vs finite, 24 windows)
  Fig 7  results_quota.txt                         (budget recipe per regime)
Run with the Windows Python (matplotlib):  python _kppwork/make_figs.py
"""
import csv
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIGS = os.path.join(ROOT, "figs")
os.makedirs(FIGS, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 110, "savefig.dpi": 220, "font.size": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "pdf.fonttype": 42,
})
CN = {8: "#1f77b4", 12: "#ff7f0e", 16: "#d62728"}
METC = {"d3": "#4c72b0", "m8": "#c44e52"}


def save(fig, name):
    for ext in ("png", "pdf"):
        p = os.path.join(FIGS, "%s.%s" % (name, ext))
        fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


# --------------------------------------------------------------- Fig 2
def parse_spectrum_orders(path):
    """-> {'r_sat': {'trunc': {n: {q: pct}}, 'mass': {n: {k: val}}}}, ...}"""
    out = {}
    scen = None
    for line in open(path):
        m = re.match(r"scenario (\S+)", line)
        if m:
            scen = m.group(1)
            out.setdefault(scen, {"trunc": {}, "mass": {}})
            mode = None
            continue
        if scen is None:
            continue
        if line.startswith("Irreducible floor"):
            mode = "trunc"
            continue
        if line.startswith("Spectral mass"):
            mode = "mass"
            continue
        m = re.match(r"\s*n=(\d+)", line)
        if m and mode:
            n = int(m.group(1))
            if mode == "trunc":
                out[scen]["trunc"][n] = {int(q): float(v) for q, v in
                                         re.findall(r"q<=(\d+):\s*([\d.]+)%", line)}
            else:
                out[scen]["mass"][n] = {int(k): float(v) for k, v in
                                        re.findall(r"k=(\d+):([\d.]+)", line)}
    return out


def _fig2_panels(dat, metric_label, prefix, suptitle):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4))
    ax = axes[0]
    for n in (8, 12, 16):
        t = dat["r_sat"]["trunc"][n]
        ax.plot(list(t), [t[q] for q in t], "o-", color=CN[n],
                label="n=%d" % n)
    ax.set_yscale("log")
    ax.set_xlabel("truncation order q")
    ax.set_ylabel("truncation error (% of max|phi|)")
    ax.set_title("(a) low-order truncation (r_sat, %s)" % metric_label,
                 loc="left")
    ax.legend(fontsize=7)
    ax = axes[1]
    for n in (8, 12, 16):
        t = dat["r_sat"]["mass"][n]
        ax.plot(list(t), [t[k] for k in t], "o-", color=CN[n], label="n=%d" % n)
    ax.set_xlabel("interaction order k")
    ax.set_ylabel("order mass $\Sigma|m(A)|$ (k = 1 normalized)")
    ax.set_title("(b) spectrum by order (r_sat, %s)" % metric_label, loc="left")
    fig.tight_layout()
    save(fig, prefix)


def fig2():
    """main text: MDA8 only (the regulatory metric); d3_mean goes to SI."""
    m8 = parse_spectrum_orders(os.path.join(HERE, "spectrum_orders_m8.txt"))
    _fig2_panels(m8, "MDA8", "fig2_truncation", None)


def fig2_supp():
    """supplementary: same panels on the daily-mean metric."""
    d3 = parse_spectrum_orders(os.path.join(HERE, "spectrum_orders.txt"))
    _fig2_panels(d3, "d3_mean", "figS1_truncation_d3", None)


# --------------------------------------------------------------- Fig 3
def mobius_fast(vals, n):
    f = np.array(vals, dtype=float)
    for i in range(n):
        bit = 1 << i
        idx = np.arange(1 << n)
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def load_oracle(path, n, col="v_ppb", scen=None):
    v = np.zeros(1 << n)
    with open(os.path.join(HERE, path)) as fh:
        for row in csv.DictReader(fh):
            name, s = row["name"], None
            if scen and "_m" in name:
                s, ms = name.rsplit("_m", 1)
                if s != scen:
                    continue
                mask = int(ms)
            else:
                mask = int(name[1:])
            if row[col] in ("", "NaN"):
                continue
            v[mask] = float(row[col])
    return v


def order_mass(v, n):
    m = mobius_fast(v, n)
    mass = np.zeros(n + 1)
    for mask in range(1, 1 << n):
        mass[bin(mask).count("1")] += abs(m[mask])
    return mass / mass[1]


def embed8(m8):
    """4-window (6 h) mask -> 8-window (3 h) mask; each 6 h block = its two
    adjacent 3 h blocks switched together (sub-lattice of the n=16 game)."""
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def load_shards(scen, n, col="v_mda8"):
    """co_v_m8_k.csv shards: name '<scen>_m<mask>' -> {mask: value}."""
    v = np.zeros(1 << n)
    for fn in sorted(os.listdir(HERE)):
        if not (fn.startswith("co_vh_m8_") and fn.endswith(".csv")):
            continue
        with open(os.path.join(HERE, fn)) as fh:
            for row in csv.DictReader(fh):
                name = row["name"]
                if "_m" not in name:
                    continue
                s, ms = name.rsplit("_m", 1)
                if s != scen or row[col] in ("", "NaN"):
                    continue
                v[int(ms)] = float(row[col])
    return v


def fig3():
    """MDA8 only (metric policy: d3_mean goes to SI); three curves that
    isolate the two clean ladders: 8->12 adds CO (weak player), 8->16
    refines the time axis (lifts the spectrum)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=True)
    for ax, scen in zip(axes, ("r_sat", "r_lim1")):
        v16 = load_oracle("n16_vh_%s_m8.csv" % scen, 16, col="v_mda8")
        v12 = load_shards(scen, 12)
        v8 = np.zeros(256)
        for m in range(256):
            v8[m] = v16[embed8(m)]
        for v, n, lab in ((v8, 8, "n = 8"),
                          (v12, 12, "n = 12"),
                          (v16, 16, "n = 16")):
            mass = order_mass(v, n)
            ax.plot(range(1, n + 1), mass[1:], "o-", color=CN[n], label=lab)
        ax.set_yscale("log")
        ax.set_ylim(0.05, 500)
        ax.set_xlabel("interaction order k")
        ax.set_title("scenario %s" % scen, loc="left")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("order mass $\\Sigma|m(A)|$ (k = 1 normalized, log scale)")
    fig.tight_layout()
    save(fig, "fig3_spectrum")


def fig3_supp_mass():
    """SI companion to Fig 3: the Shapley-weighted order mass
    W(k) = sum_{|A|=k}|m(A)|/|A|, the weighting whose tail bounds the Shapley
    truncation error.  The main-text figures use the unweighted M(k)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4), sharey=True)
    for ax, scen in zip(axes, ("r_sat", "r_lim1")):
        v16 = load_oracle("n16_vh_%s_m8.csv" % scen, 16, col="v_mda8")
        v12 = load_shards(scen, 12)
        v8 = np.zeros(256)
        for m in range(256):
            v8[m] = v16[embed8(m)]
        for v, n, lab in ((v8, 8, "n = 8 (2 sectors, 4 windows)"),
                          (v12, 12, "n = 12 (+ CO, 3 sectors)"),
                          (v16, 16, "n = 16 (2 sectors, 8 windows)")):
            m = mobius_fast(v, n)
            mass = np.zeros(n + 1)
            for mask in range(1, 1 << n):
                mass[bin(mask).count("1")] += abs(m[mask]) / bin(mask).count("1")
            mass /= mass[1]
            ax.plot(range(1, n + 1), mass[1:], "o-", color=CN[n], label=lab)
        ax.set_xlabel("interaction order k")
        ax.set_title("scenario %s" % scen, loc="left")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("Shapley-weighted mass (k = 1 normalized)")
    fig.tight_layout()
    save(fig, "figS2_mass_weighted")


# --------------------------------------------------------------- Fig 4
def parse_scaling(path):
    """-> {metric: {scen: {n: [(B, err%, q1%, q3%)]}}}, {metric: {scen: {n: runs or None}}}
    q1/q3 are the seed quartiles (may be nan for pre-2026-10 files)."""
    res, runs = {}, {}
    metric, scen = None, None
    mode = None
    for line in open(path):
        if "PART A" in line:
            metric = "d3" if "d3_mean" in line else "m8"
            res.setdefault(metric, {})
            runs.setdefault(metric, {})
            continue
        if metric is None:
            continue
        m = re.match(r"\s*scenario (\S+)", line)
        if m:
            scen = m.group(1)
            res[metric].setdefault(scen, {})
            runs[metric].setdefault(scen, {})
            mode = "err"
            continue
        if "runs needed" in line:
            mode = "runs"
            continue
        if mode == "err":
            m = re.match(r"\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)%"
                         r"(?:\s+([\d.]+)%\s+([\d.]+)%)?", line)
            if m:
                res[metric][scen].setdefault(int(m.group(1)), []).append(
                    (int(m.group(2)), float(m.group(4)),
                     float(m.group(5)) if m.group(5) else float("nan"),
                     float(m.group(6)) if m.group(6) else float("nan")))
        elif mode == "runs":
            m = re.match(r"\s+(\S+)\s+(\d+)\s+(\S+)", line)
            if m and m.group(1) in ("r_sat", "r_lim1"):
                runs[metric][m.group(1)][int(m.group(2))] = (
                    int(m.group(3)) if m.group(3).isdigit() else None)
    return res, runs


def fig4():
    res, runs = parse_scaling(os.path.join(HERE, "scaling_dual_metrics.txt"))
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.4))
    ax = axes[0]
    for metric, ls, lab in (("m8", "-", "MDA8"), ("d3", "--", "d3_mean")):
        for n in (8, 12, 16):
            pts = sorted(res[metric]["r_sat"][n])
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.plot(xs, ys, ls + "o", color=CN[n],
                    alpha=1.0 if metric == "m8" else 0.55,
                    label="%s, n=%d" % (lab, n))
            q1 = [p[2] for p in pts]
            q3 = [p[3] for p in pts]
            if not any(math.isnan(q) for q in q1 + q3):
                ax.fill_between(xs, q1, q3, color=CN[n], alpha=0.12,
                                linewidth=0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("coalition evaluations (budget B)")
    ax.set_ylabel("single-run error (% of max|phi|)")
    ax.set_title("(a) sampling error vs cost (r_sat; band = IQR over seeds)",
                 loc="left", fontsize=9)
    ax.legend(fontsize=6.5, ncol=2)
    ax = axes[1]
    ns = np.array([8, 12, 16])
    ax.plot(ns, 2.0 ** ns, "k-", label="full enumeration 2$^n$")
    for scen, ls, lab in (("r_lim1", "-", "d3_mean, r_lim1"),
                          ("r_sat", "--", "d3_mean, r_sat")):
        rr = [runs["d3"][scen].get(n) for n in ns]
        if any(rr):
            ax.plot(ns, [x if x else np.nan for x in rr], ls + "o",
                    color=METC["d3"],
                    label="perm_anti, %s (1.5%%, interpolated)" % lab)
    ax.set_yscale("log")
    ax.set_xticks(ns)
    ax.set_xlabel("players n")
    ax.set_ylabel("runs for 1.5% error")
    ax.set_title("(b) cost of equal precision", loc="left")
    ax.annotate("MDA8: 1.5% not reached\nwithin tested budgets (both scenarios)",
                xy=(0.98, 0.03), xycoords="axes fraction", ha="right",
                fontsize=7.5, color=METC["m8"])
    ax.legend(fontsize=7)
    fig.tight_layout()
    save(fig, "fig4_cost")


# --------------------------------------------------------------- Fig 5/6
def parse_lvg(path):
    """-> {sector: [(h0, D_local, D_finite)]} for the MDA8 3-hour block."""
    rows, in_m8, in_3h = [], False, False
    for line in open(path):
        if "MDA8 (primary metric)" in line:
            in_m8 = True
            continue
        if in_m8 and "d3_mean" in line and line.startswith("#"):
            break
        if in_m8 and "3-hour windows" in line:
            in_3h = True
            continue
        if in_3h and "6-hour blocks" in line:
            break
        if in_3h:
            m = re.match(r"(NOXW|VOCW|COW)\s+(\d+)-(\d+)\s+(-?[\d.]+)\s+(-?[\d.]+)",
                         line)
            if m:
                rows.append((m.group(1), int(m.group(2)),
                             float(m.group(4)), float(m.group(5))))
    out = {}
    for sec, h0, dloc, dfin in rows:
        out.setdefault(sec, []).append((h0, dloc, dfin))
    return out


def fig5():
    lvg = parse_lvg(os.path.join(ROOT, "results_local_vs_global.txt"))
    nox = lvg["NOXW"]
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    xs = [h for h, _, _ in nox]
    ys = [df for _, _, df in nox]
    cols = ["#c44e52" if y < 0 else "#4c72b0" for y in ys]
    ax.bar(xs, ys, width=2.4, color=cols)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels(["%d-%d" % (h, h + 3) for h in xs], fontsize=7)
    ax.set_xlabel("window (hour of day)")
    ax.set_ylabel("finite switch-off effect on MDA8 (ppb)")
    ax.set_title("NOx windows: one sign reversal at 12-15 h", loc="left")
    fig.tight_layout()
    save(fig, "fig5_signflip")


def fig6():
    lvg = parse_lvg(os.path.join(ROOT, "results_local_vs_global.txt"))
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.4), sharex=True)
    for ax, sec, title in zip(axes, ("NOXW", "VOCW", "COW"),
                              ("NOx", "VOC", "CO")):
        rows = lvg[sec]
        xs = np.arange(len(rows))
        dloc = [d for _, d, _ in rows]
        dfin = [f for _, _, f in rows]
        ax.bar(xs - 0.22, dloc, width=0.22, color="#9ecae1",
               label=r"$D_{local}$")
        ax.bar(xs, [-x for x in dloc], width=0.22, color="#4c72b0",
               label=r"-$D_{local}$ (linear)")
        ax.bar(xs + 0.22, dfin, width=0.22, color="#c44e52",
               label=r"$\Delta_{off}$ (finite)")
        ax.axhline(0, color="k", lw=0.7)
        ax.set_ylabel("MDA8 (ppb)")
        ax.set_title(title, loc="left", fontsize=9)
        if sec == "NOXW":
            ax.annotate("sign flip", xy=(4, dfin[4]), xytext=(4.4, dfin[4] + 4),
                        fontsize=8, color="#c44e52",
                        arrowprops=dict(arrowstyle="->", color="#c44e52"))
            ax.legend(fontsize=7, ncol=3)
    axes[-1].set_xticks(range(8))
    axes[-1].set_xticklabels(["%d-%d" % (3 * i, 3 * i + 3) for i in range(8)],
                             fontsize=7)
    axes[-1].set_xlabel("window (hour of day)")
    fig.tight_layout()
    save(fig, "fig6_lvg")


# --------------------------------------------------------------- Fig 7
def fig7():
    rows = []
    for line in open(os.path.join(ROOT, "results_quota.txt")):
        m = re.match(r"\s+(r_\S+|u0_base)\s+([\d.]+)\s+(\S+)\s+(\S+)\s+([+-][\d.]+)",
                     line)
        if m:
            rows.append((m.group(1), m.group(3), m.group(4), float(m.group(5))))
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    labels = []
    for i, (scen, voc, nox, delta) in enumerate(rows):
        y = len(rows) - i
        labels.append(scen)
        for base, wnames, col in ((0.13, voc, "#4c72b0"), (-0.13, nox, "#ff7f0e")):
            wins = [int(w) - 1 for w in re.findall(r"W(\d)", wnames)]
            if "all" in wnames:
                wins = list(range(4))
            for w in range(4):
                cut = w in wins
                ax.barh(y + base, 1.0, left=w, height=0.22,
                        color=col if cut else "#e0e0e0",
                        edgecolor="white", linewidth=0.6)
    ax.set_yticks(range(1, len(rows) + 1))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xticks([0.5, 1.5, 2.5, 3.5])
    ax.set_xticklabels(["W1 00-06", "W2 06-12", "W3 12-18", "W4 18-24"], fontsize=8)
    ax.set_xlim(0, 4)
    ax.set_xlabel("coloured = window switched off in the optimal recipe")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#ff7f0e", label="NOx cut"),
                       Patch(color="#4c72b0", label="VOC cut"),
                       Patch(color="#e0e0e0", label="left on")],
              fontsize=8, loc="lower right")
    ax.set_title(r"optimal recipe at budget $(R_N,R_V)=(1.0,0.35)$ flips across regimes",
                 loc="left", fontsize=9)
    fig.tight_layout()
    save(fig, "fig7_quota")


if __name__ == "__main__":
    fig2()
    fig2_supp()
    fig3()
    fig3_supp_mass()
    fig4()
    fig5()
    fig6()
    fig7()
    print("all figures ->", FIGS)
