#!/usr/bin/env python3
"""Pilot: can the truncation tail t_i = sum_{A ni i, |A|>k} m(A)/|A| be imputed
from LOW-ORDER information, instead of spread uniformly (Shapley-Taylor)?

Shapley-Taylor: t_hat_i = r(k)/n with r(k) = sum_{|A|>k} m(A) -- computable
cheaply as r(k) = v(N) - v(empty) - sum_{|A|<=k} m(A), so an order-k pipeline
costs sum_{j<=k} C(n,j) + 1 runs.  We test learned redistribution
t_hat_i = r(k) * w_i(phi^{<=k}, m({i}), pairs, sector, window), sum w_i = 1.

Methods (all renormalised to sum t_hat = r(k)):
  uniform   r(k)/n                        (Shapley-Taylor)
  prop_low  r(k) * phi_i^{<=k} / sum phi^{<=k}
  prop_1    r(k) * m({i}) / sum m({j})
  ridge     leave-one-scenario-out ridge on low-order features
Error: max_i |t_i - t_hat_i| / max_i |phi_i|.  Zero chemistry runs.
"""
import csv
import math
import os
import platform

import numpy as np
from sklearn.linear_model import Ridge

if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"

N = 16
LABELS = ["r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"]
KS = [2, 5]


def mobius_np(vals):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << N)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def read_vals(label, col):
    if label in ("dallas_2016w", "deerpark_2017s"):
        v = np.zeros(65536)
        for K in range(1, 13):
            p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (label, K))
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                for row in csv.DictReader(fh):
                    val = row.get(col, "")
                    if val not in ("", "NaN"):
                        v[int(row["name"][1:])] = float(val)
        return v
    v = np.zeros(65536)
    with open(os.path.join(W, "n16_v_%s_m8.csv" % label)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


def shapley_and_weights(m):
    idx = np.arange(1, 1 << N)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(N)
    for i in range(N):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    return phi, idx, k, w


def build(label, ktrunc, col="v_mda8"):
    v = read_vals(label, col)
    m = mobius_np(v)
    phi, idx, k_all, w = shapley_and_weights(m)
    scale = np.abs(phi).max() or 1.0
    lo = k_all <= ktrunc
    phi_lo = np.zeros(N)
    for i in range(N):
        phi_lo[i] = w[lo & (((idx >> i) & 1) > 0)].sum()
    t_true = phi - phi_lo
    r_k = m[idx[k_all > ktrunc]].sum()
    # low-order features
    single = np.array([m[1 << i] for i in range(N)])
    pairs = np.zeros(N)
    for i in range(N):
        for j in range(N):
            if j != i:
                pairs[i] += m[(1 << i) | (1 << j)]
    feat = np.column_stack([
        phi_lo / scale,
        single / scale,
        pairs / scale,
        (np.arange(N) >= 8).astype(float),          # sector: VOC=1
        (np.arange(N) % 8) / 7.0,                   # window position
        np.full(N, r_k / scale),                    # scenario-level residual
    ])
    return t_true / scale, r_k / scale, feat, scale


L = []
add = L.append
add("Tail imputation pilot: t_hat_i = r(k) * w_i, learned vs uniform (Shapley-Taylor)")
add("error = max_i |t_i - t_hat_i| / max_i |phi_i|, MDA8, n=16")
add("budget of order-k pipeline = sum_{j<=k} C(16,j) + 1 runs")
add("")

data = {lab: {k: build(lab, k) for k in KS} for lab in LABELS}

for k in KS:
    budget = sum(math.comb(16, j) for j in range(k + 1)) + 1
    add("=== k=%d   budget = %d runs ===" % (k, budget))
    add("%-18s %9s %9s %9s %9s" % ("scenario", "uniform", "prop_low", "prop_1", "ridge(LOSO)"))
    err_all = {m_: [] for m_ in ("uniform", "prop_low", "prop_1", "ridge")}
    for lab in LABELS:
        t_true, r_k, feat, scale = data[lab][k]
        # common denominator for the reported error: max|phi|
        denom = 1.0          # t and r already normalised by max|phi|
        res = {}
        # uniform
        res["uniform"] = np.abs(t_true - r_k / N).max() / denom
        # proportional to low-order phi
        lo_sum = feat[:, 0].sum()
        w_lo = feat[:, 0] / lo_sum if abs(lo_sum) > 1e-12 else np.full(N, 1.0 / N)
        res["prop_low"] = np.abs(t_true - r_k * w_lo).max() / denom
        # proportional to singleton
        s_sum = feat[:, 1].sum()
        w_1 = feat[:, 1] / s_sum if abs(s_sum) > 1e-12 else np.full(N, 1.0 / N)
        res["prop_1"] = np.abs(t_true - r_k * w_1).max() / denom
        # ridge LOSO
        Xtr, ytr = [], []
        for lab2 in LABELS:
            if lab2 == lab:
                continue
            t2, r2, f2, _ = data[lab2][k]
            Xtr.append(f2)
            ytr.append(t2 - r2 / N)          # learn the DEVIATION from uniform
        Xtr = np.vstack(Xtr)
        ytr = np.concatenate(ytr)
        mdl = Ridge(alpha=1.0).fit(Xtr, ytr)
        dev = mdl.predict(feat)
        t_hat = r_k / N + dev
        t_hat = t_hat * (r_k / t_hat.sum()) if abs(t_hat.sum()) > 1e-12 else t_hat
        res["ridge"] = np.abs(t_true - t_hat).max() / denom
        for key in err_all:
            err_all[key].append(res[key])
        add("%-18s %8.1f%% %8.1f%% %8.1f%% %13.1f%%" %
            (lab, 100 * res["uniform"], 100 * res["prop_low"],
             100 * res["prop_1"], 100 * res["ridge"]))
    add("%-18s %8.1f%% %8.1f%% %8.1f%% %13.1f%%" %
        ("mean",
         100 * np.mean(err_all["uniform"]), 100 * np.mean(err_all["prop_low"]),
         100 * np.mean(err_all["prop_1"]), 100 * np.mean(err_all["ridge"])))
    add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_tail_impute.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_tail_impute.txt")
