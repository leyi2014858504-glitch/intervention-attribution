#!/usr/bin/env python3
"""Idea C: is a surrogate useful as a CONTROL VARIATE (not as a predictor)?

Estimator:
    phi_hat = phi_surrogate  +  mean_pi[ Delta_true - Delta_surrogate ]
The surrogate term is computed on the surrogate itself (free); the residual is
sampled.  The estimator is unbiased for ANY surrogate, so the surrogate's own
error floor is harmless -- what decides the gain is the correlation rho between
Delta_true and Delta_surrogate across sampled contexts:

    speed-up  ~  Var(Delta_true) / Var(Delta_true - Delta_surrogate) = 1/(1-rho^2)

Two surrogates:
  gbm    gradient boosting on the 16 switch bits (the earlier failed predictor)
  mob2   least-squares fit of the order<=2 Moebius model (linear + pairwise)

Budget accounting is honest: B_train true runs to fit the surrogate + B_sample
true runs to sample the residual, compared against plain perm_anti run at the
SAME total budget.

Writes results_surrogate_cv.txt
"""
import os
import random

import numpy as np
from sklearn.ensemble import GradientBoostingRegressor

import _sampler_city as sc

OUT = os.path.join(sc.ROOT, "results_surrogate_cv.txt")
B_TRAIN = 512
B_SAMPLES = (256, 512, 1024)
N = 16
ALL = np.arange(1 << N)
BITS = ((ALL[:, None] >> np.arange(N)) & 1).astype(float)


def fit_gbm(v, masks):
    g = GradientBoostingRegressor(n_estimators=200, max_depth=3,
                                  learning_rate=0.08, subsample=0.8,
                                  random_state=0)
    g.fit(BITS[masks], v[masks] - v[0])
    return g.predict(BITS) + v[0]


def fit_mob2(v, masks):
    """order<=2 Moebius model by least squares on the sampled coalitions."""
    cols = []
    cols.append(np.ones(len(masks)))
    for i in range(N):
        cols.append(BITS[masks, i])
    for i in range(N):
        for j in range(i + 1, N):
            cols.append(BITS[masks, i] * BITS[masks, j])
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, v[masks] - v[0], rcond=None)
    B = []
    B.append(np.ones(len(ALL)))
    for i in range(N):
        B.append(BITS[:, i])
    for i in range(N):
        for j in range(i + 1, N):
            B.append(BITS[:, i] * BITS[:, j])
    return np.column_stack(B) @ coef + v[0]


def residual_sample(v, vs, n, budget, seed, track=None):
    """perm+antithetic over v, accumulating Delta_true and Delta_surrogate."""
    rng = random.Random(seed)
    seen = np.zeros(1 << n, dtype=bool)
    seen[0] = True
    used = 0
    dt = [[] for _ in range(n)]
    ds = [[] for _ in range(n)]
    stall = 0
    while used < budget and stall < 3000:
        before = used
        order = list(range(n))
        rng.shuffle(order)
        for o in (order, order[::-1]):
            S = 0
            for i in o:
                prev = S
                S |= 1 << i
                if not seen[S]:
                    seen[S] = True
                    used += 1
                dt[i].append(v[S] - v[prev])
                ds[i].append(vs[S] - vs[prev])
        stall = stall + 1 if used == before else 0
    return dt, ds


def plain_sample(v, n, budget, seed):
    dt, _ = residual_sample(v, np.zeros_like(v), n, budget, seed)
    return np.array([np.mean(x) for x in dt])


def main():
    L, add = [], None
    add = L.append
    add("Surrogate as a CONTROL VARIATE (unbiased for any surrogate)")
    add("total budget = %d training runs + residual-sampling runs" % B_TRAIN)
    add("gain = plain error / control-variate error, same total cost")
    add("rho = corr(Delta_true, Delta_surrogate) at the dominant player")
    add("")

    rng = np.random.default_rng(7)
    for label in sc.LABELS:
        for col, mname in sc.METRICS:
            v = sc.load_table(label, col)
            v = v - v[0]
            true = sc.exact_shapley(v, N)
            scale = np.abs(true).max()
            dom = int(np.argmax(np.abs(true)))

            tr = rng.choice(ALL[1:], size=B_TRAIN, replace=False)
            tr = np.unique(np.concatenate([[0], tr]))
            surr = {"gbm": fit_gbm(v, tr), "mob2": fit_mob2(v, tr)}

            add("=" * 90)
            add("%s   metric = %s   max|phi| = %.4f   dominant %d"
                % (label, mname, scale, dom))
            add("=" * 90)
            for sname, vs in surr.items():
                phis = sc.exact_shapley(vs - vs[0] if sname == "gbm" else vs, N)
                if sname == "gbm":
                    phis = sc.exact_shapley(vs - vs[0], N)
                # surrogate's own shapley error (as a predictor)
                perr = 100 * np.abs(phis - true).max() / scale
                add("  surrogate %s: own Shapley error as predictor = %.1f%%" % (sname, perr))
                add("    %6s %10s %10s %10s %10s %8s"
                    % ("B_samp", "total", "plain", "cv", "gain", "rho"))
                for Bs in B_SAMPLES:
                    total = B_TRAIN + Bs
                    ep, ec, rhos = [], [], []
                    for s in range(sc.NSEED):
                        dt, ds = residual_sample(v, vs, N, Bs, 2000 + s)
                        est_p = np.array([np.mean(x) for x in dt])
                        est_c = phis + np.array([np.mean([a - b for a, b in zip(dt[i], ds[i])])
                                                 for i in range(N)])
                        ep.append(100 * np.abs(est_p - true).max() / scale)
                        ec.append(100 * np.abs(est_c - true).max() / scale)
                        a = np.array(dt[dom]); b = np.array(ds[dom])
                        if a.std() > 0 and b.std() > 0:
                            rhos.append(float(np.corrcoef(a, b)[0, 1]))
                    # plain at the SAME total budget
                    efull = []
                    for s in range(sc.NSEED):
                        est = plain_sample(v, N, total, 3000 + s)
                        efull.append(100 * np.abs(est - true).max() / scale)
                    mfull = float(np.median(efull))
                    mc = float(np.median(ec))
                    add("    %6d %10d %9.2f%% %9.2f%% %9.2fx %8.3f"
                        % (Bs, total, mfull, mc,
                           (mfull / mc if mc > 0 else float("nan")),
                           float(np.median(rhos))))
                add("")
            add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
