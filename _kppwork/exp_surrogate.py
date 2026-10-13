#!/usr/bin/env python3
"""Surrogate attribution (fit-then-explain) vs direct sampling, oracle-scored.

Design: at equal budget B of distinct coalition evaluations, two arms
  1. GBM surrogate: sample B distinct non-empty coalitions uniformly, fit a
     gradient-boosted tree ensemble on the coalition bit-vector, then take the
     EXACT Shapley of the fitted model (evaluate it on all 2^n masks, Mobius,
     Shapley).  For a tree ensemble this is the quantity TreeSHAP computes in
     polynomial time; we obtain it by exhaustive evaluation of the surrogate,
     which is free of chemistry and mathematically identical.  This is the
     "reduced representation" the paper's exact games are meant to test.
  2. perm_anti: antithetic permutation sampling of v directly (the working
     sampler of the paper; chains complete past the budget, as in Methods).
Error: max_i |phi_hat - phi| / max_i |phi| against the exact enumerated
oracle, per seed; median [q1, q3] over NSEED seeds.  Also the surrogate's own
value error err_v = max|v_hat - v| / max|v| over all masks (diagnostic).

Hyper-parameters fixed a priori: HistGradientBoostingRegressor(100 iterations,
max_depth=3, lr=0.1).  Budget currency = distinct non-empty coalition
evaluations; v(empty)=0 is free for both arms.  Zero chemistry runs.

Run with the Windows Python (sklearn):  python _kppwork/exp_surrogate.py
Writes results_surrogate.txt
"""
import csv
import math
import os
import random

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCEN = ("r_sat", "r_lim1")
NSEED = 24
BUDGETS = {8: [8, 16, 32, 64, 96, 128, 192],
           12: [64, 128, 256, 512, 1024],
           16: [64, 128, 256, 512, 1024]}


# --------------------------------------------------------------- oracle load
def read_tagged(path, col, keep):
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            name = row.get("name", "")
            if "_m" not in name:
                continue
            s, _, ms = name.rpartition("_m")
            if s not in keep:
                continue
            v = row.get(col, "")
            if v in ("", "NaN"):
                continue
            out.setdefault(s, {})[int(ms)] = float(v)
    return out


def read_n16(scen, col):
    v = {}
    with open(os.path.join(HERE, "n16_v_%s_m8.csv" % scen)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


def embed8(m8):
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def load_oracle(metric):
    """-> {scen: {n: np.array of v over all 2^n masks}}"""
    out = {s: {} for s in SCEN}
    if metric == "d3_mean":
        v8 = read_tagged(os.path.join(HERE, "scen_v.csv"), "v_ppb", SCEN)
        v12 = read_tagged(os.path.join(HERE, "co_v.csv"), "v_ppb", SCEN)
        for s in SCEN:
            v16 = read_n16(s, "v_ppb")          # read ONCE per scenario
            out[s][8] = np.array([v8[s][m] for m in range(256)])
            out[s][12] = np.array([v12[s][m] for m in range(4096)])
            out[s][16] = np.array([v16[m] for m in range(65536)])
    else:  # MDA8
        v12 = {}
        for K in range(1, 7):
            p = os.path.join(HERE, "co_v_m8_%d.csv" % K)
            if os.path.exists(p):
                part = read_tagged(p, "v_mda8", SCEN)
                for s, d in part.items():
                    v12.setdefault(s, {}).update(d)
        for s in SCEN:
            v16 = read_n16(s, "v_mda8")
            out[s][16] = np.array([v16[m] for m in range(65536)])
            out[s][8] = np.array([v16[embed8(m)] for m in range(256)])
            out[s][12] = np.array([v12[s][m] for m in range(4096)])
    return out


# --------------------------------------------------------- spectrum / Shapley
def mobius_np(vals, n):
    f = vals.astype(float).copy()
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def shapley_np(vals, n):
    m = mobius_np(vals, n)
    idx = np.arange(1, 1 << n)
    k = order_k(idx)
    w = m[idx] / k
    phi = np.zeros(n)
    for i in range(n):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    return phi


# -------------------------------------------------------------------- arms
_ALL_X = {}
_ALL_K = {}


def masks_matrix(masks, n):
    idx = np.asarray(masks, dtype=np.int64)
    return ((idx[:, None] >> np.arange(n, dtype=np.int64)) & 1).astype(np.float64)


def all_x(n):
    if n not in _ALL_X:
        _ALL_X[n] = masks_matrix(range(1 << n), n)
    return _ALL_X[n]


def order_k(idx):
    """popcount vector for the positive mask indices (cached)."""
    key = (idx[0], idx[-1], len(idx))
    if key not in _ALL_K:
        _ALL_K[key] = np.array([bin(x).count("1") for x in idx], dtype=float)
    return _ALL_K[key]


def arm_gbm(vals, n, B, seed):
    rng = random.Random(1000 * B + seed)
    train = rng.sample(range(1, 1 << n), B)
    Xtr = all_x(n)[train]
    ytr = vals[train]
    mdl = HistGradientBoostingRegressor(max_iter=100, max_depth=3,
                                        learning_rate=0.1, random_state=seed)
    mdl.fit(Xtr, ytr)
    vhat = mdl.predict(all_x(n))
    phi = shapley_np(vhat, n)
    scale = np.abs(vals).max() or 1.0
    errv = np.abs(vhat - vals).max() / scale
    return phi, B, errv


def arm_perm_anti(vals, n, B, seed):
    rng = random.Random(1000 * B + seed)
    s = np.zeros(n)
    c = np.zeros(n)
    seen = set()
    while len(seen) < B:
        order = list(range(n))
        rng.shuffle(order)
        for seq in (order, order[::-1]):
            mask = 0
            for i in seq:
                prev = mask
                mask |= 1 << i
                seen.add(mask)
                s[i] += vals[mask] - vals[prev]
                c[i] += 1
    return s / c, len(seen), float("nan")


def quant(xs, p):
    ys = sorted(xs)
    k = (len(ys) - 1) * p
    f, c = int(math.floor(k)), int(math.ceil(k))
    return ys[f] if f == c else ys[f] + (ys[c] - ys[f]) * (k - f)


def main():
    L = []
    add = L.append
    add("Surrogate attribution vs direct sampling, scored on the exact oracle")
    add("arm 1 = HistGBM(100 iters, depth 3, lr 0.1) on B uniform coalitions,")
    add("        exact Shapley of the fitted model (the TreeSHAP quantity;")
    add("        obtained here by exhaustive evaluation of the surrogate)")
    add("arm 2 = perm_anti (antithetic permutation sampling of v directly)")
    add("error = max_i |phi_hat - phi| / max_i |phi|, %d seeds; median [q1,q3]"
        % NSEED)
    add("budget = distinct non-empty coalition evaluations (v(empty) free)")
    for metric in ("MDA8", "d3_mean"):
        orc = load_oracle(metric)
        add("")
        add("=" * 78)
        add("metric = %s" % metric)
        add("=" * 78)
        for s in SCEN:
            for n in (8, 12, 16):
                vals = orc[s][n]
                phi = shapley_np(vals, n)
                scale = np.abs(phi).max() or 1.0
                add("")
                add("scenario %s  n=%d   max|phi| = %.4f" % (s, n, scale))
                add("  %6s | %22s %22s %8s | %22s"
                    % ("B", "GBM+exact-SHAP", "err_v(GBM)", "", "perm_anti"))
                for B in BUDGETS[n]:
                    if B > (1 << n) - 1:
                        continue
                    eg, ep, ev = [], [], []
                    for sd in range(NSEED):
                        phig, used, errv = arm_gbm(vals, n, B, sd)
                        eg.append(np.abs(phig - phi).max() / scale)
                        ev.append(errv)
                        phip, usedp, _ = arm_perm_anti(vals, n, B, sd)
                        ep.append(np.abs(phip - phi).max() / scale)
                    add("  %6d | %8.2f%% [%6.2f,%6.2f] %8.2f%% [%6.2f,%6.2f] %7s | %8.2f%% [%6.2f,%6.2f]"
                        % (B,
                           100 * quant(eg, .5), 100 * quant(eg, .25), 100 * quant(eg, .75),
                           100 * quant(ev, .5), 100 * quant(ev, .25), 100 * quant(ev, .75),
                           "",
                           100 * quant(ep, .5), 100 * quant(ep, .25), 100 * quant(ep, .75)))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(ROOT, "results_surrogate.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote results_surrogate.txt")


if __name__ == "__main__":
    main()