#!/usr/bin/env python3
"""Calibration study: can oracle-free flags PREDICT the realized error, not
merely rank it?  Pure post-processing on the existing exact oracles.

Three questions (answers go in selfdiag_calib.txt):
  Q1 slope     : log e_q = a + b log h_q -- is b stable across corners
                 (scenario x n x metric), and does a leave-one-corner-out
                 fit predict the held-out corner?
  Q2 multiplier: c = median(e/h) per corner -- does ONE c per metric
                 transfer?  (spread of c across corners = transfer error)
  Q3 coverage  : what safety factor k makes k * c * h >= e_q in 90% / 95%
                 of cases -- the practical "error bar without truth" recipe.
Same for the bootstrap flag sp (expected worse) and the efficiency flag g.

Protocol matches exp_selfdiag.py exactly (same samplers, fits, splits).
"""
import importlib.util
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("sd", os.path.join(HERE, "exp_selfdiag.py"))
sd = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sd)

SEEDS = {8: 4, 12: 4, 16: 2}     # trimmed: the n=16 lstsq (6500x2516) is heavy
DELTA = 0.05


def collect(corner, metric, Bs, qcap):
    scen, n = corner
    v = sd.load_table(
        {("v_ppb", 8): "scen_v.csv", ("v_ppb", 12): "co_v.csv",
         ("v_ppb", 16): "n16_v_%s.csv" % scen,
         ("v_mda8", 12): "co_v_m8_*.csv", ("v_mda8", 16): "n16_v_%s_m8.csv" % scen
         }[(metric, n)],
        None if n == 16 else scen, col=metric, shards=(metric == "v_mda8" and n == 12))
    if len(v) != 1 << n:
        return None
    phi_true, t_q, scale = sd.true_decomposition(v, n)
    Nmask = (1 << n) - 1
    vN = v[Nmask]
    universe = np.arange(1, Nmask)
    rows = []
    for B in Bs:
        for seed in range(SEEDS[n]):
            rng = np.random.default_rng(1000 * seed + B + n)
            masks = rng.choice(universe, size=min(B, len(universe)), replace=False).tolist()
            if Nmask not in masks:
                masks.append(Nmask)
            rng.shuffle(masks)
            cut = int(0.8 * len(masks))
            fit, hold = masks[:cut], masks[cut:]
            yf = np.array([v[S] for S in fit])
            for q in range(1, qcap + 1):
                from math import comb
                P = sum(comb(n, k) for k in range(1, q + 1))
                if P > len(fit) / 2 or P > 3000:
                    break
                cols, col_of = sd.columns_upto(n, q)
                m_hat, *_ = np.linalg.lstsq(sd.design(fit, col_of), yf, rcond=None)
                phi_hat = sd.phi_from_cols(m_hat, cols, n)
                e_q = float(np.max(np.abs(phi_hat - phi_true)) / scale)
                g_q = abs(vN - float(m_hat.sum())) / (abs(vN) or 1.0)
                res = np.array([v[S] for S in hold]) - sd.design(hold, col_of) @ m_hat
                h_q = float(np.sqrt(np.mean(res ** 2)) / (float(np.std(yf)) or 1.0))
                if P <= 1000:   # bootstrap refit only while cheap (same guard
                    bs = rng.choice(np.array(fit), size=len(fit), replace=True).tolist()
                    mb, *_ = np.linalg.lstsq(sd.design(bs, col_of),
                                             np.array([v[S] for S in bs]), rcond=None)
                    sp_q = float(np.max(np.abs(sd.phi_from_cols(mb, cols, n) - phi_hat)) / scale)
                else:
                    sp_q = float("nan")
                rows.append(dict(scen=scen, n=n, B=B, q=q, e=e_q, g=g_q, h=h_q, sp=sp_q))
    return rows


def ols(x, y):
    A = np.vstack([np.ones_like(x), x]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coef  # a, b for log-log


def analyze(allrows, metric, L, add):
    rows = [r for rs in allrows for r in (rs or [])]
    if not rows:
        return
    e = np.array([r["e"] for r in rows])
    add("")
    add("=" * 68)
    add("METRIC %s   (%d fit points, corners: %s)" %
        (metric, len(rows), sorted({(r["scen"], r["n"]) for r in rows})))
    add("=" * 68)
    for flag in ("h", "sp", "g"):
        f = np.array([r[flag] for r in rows])
        m = (e > 1e-6) & (f > 1e-6)
        c = np.median(e[m] / f[m])
        ratio = e[m] / (c * f[m])
        k90, k95 = np.quantile(ratio, 0.90), np.quantile(ratio, 0.95)
        # per-corner c for transfer spread
        corners = sorted({(r["scen"], r["n"]) for r in rows})
        fv = np.array([r[flag] for r in rows])
        per_c = []
        for cn in corners:
            sel = np.array([(r["scen"], r["n"]) == cn for r in rows]) & m
            if sel.sum() >= 8:
                per_c.append(np.median(e[sel] / fv[sel]))
        spread = (max(per_c) / min(per_c)) if len(per_c) > 1 and min(per_c) > 0 else float("nan")
        # log-log slope
        b = ols(np.log(f[m]), np.log(e[m]))
        add("  flag=%-3s c=median(e/f)=%7.3f  per-corner c spread=%5.1fx  "
            "slope b=%.2f  k(90%%)=%5.1f  k(95%%)=%5.1f" %
            (flag, c, spread, b[1], k90, k95))
        # leave-one-corner-out prediction quality (median ratio)
        lloo = []
        for cn in corners:
            tr = np.array([(r["scen"], r["n"]) != cn for r in rows]) & m
            te = np.array([(r["scen"], r["n"]) == cn for r in rows]) & m
            if tr.sum() >= 20 and te.sum() >= 5:
                coef = ols(np.log(f[tr]), np.log(e[tr]))
                pred = np.exp(coef[0] + coef[1] * np.log(f[te]))
                lloo.append(np.median(e[te] / pred))
        if lloo:
            add("        LOO-corner e/pred ratio: median %.2f  range [%.2f, %.2f]" %
                (np.median(lloo), min(lloo), max(lloo)))


def main():
    L = []
    add = L.append
    add("Calibration of oracle-free flags: can they PREDICT (not just rank) the")
    add("realized error e_q?  Fit e = c*flag; transfer via leave-one-corner-out;")
    add("coverage factor k for 90/95 pct.  Protocol identical to exp_selfdiag.")
    for metric in ("v_ppb", "v_mda8"):
        allrows = []
        grid = [(s, n) for s in ("r_sat", "r_lim1")
                for n in ((8, 12, 16) if metric == "v_ppb" else (12, 16))]
        Bs = {8: [64, 128, 250], 12: [512, 1024, 2048], 16: [2048, 8192]}
        qcap = {8: 3, 12: 4, 16: 4}
        for scen, n in grid:
            add("")
            add("collecting %s n=%d %s ..." % (scen, n, metric))
            print("collect", scen, n, metric, flush=True)
            allrows.append(collect((scen, n), metric, Bs[n], qcap[n]))
        analyze(allrows, metric, L, add)
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "selfdiag_calib.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("wrote selfdiag_calib.txt")


if __name__ == "__main__":
    main()
