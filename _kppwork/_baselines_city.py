#!/usr/bin/env python3
"""Baseline table on the CITY n=16 games: KernelSHAP, plain permutation,
antithetic permutation.  Same currency (distinct coalition evaluations), same
oracle (the exact table), same seeds.

KernelSHAP here is the standard estimator: uniform coalition draws, Shapley-kernel
weighted least squares on v(S) ~ c + sum_{i in S} phi_i, with the efficiency
constraint sum phi_i = v(N) - v(0) imposed by eliminating the last coefficient.

Writes results_baselines_city.txt
"""
import csv
import glob
import os
import random

import numpy as np

import _sampler_city as sc

OUT = os.path.join(sc.ROOT, "results_baselines_city.txt")
N = 16
BUDGETS = (128, 256, 512, 1024, 2048)
FULL = (1 << N) - 1


def kernel(n, s):
    if s == 0 or s == n:
        return 1e6
    from math import comb
    return (n - 1) / (comb(n, s) * s * (n - s))


def kernelshap(v, budget, seed):
    rng = np.random.default_rng(seed)
    masks = rng.choice(np.arange(1, FULL), size=budget, replace=False)
    X = ((masks[:, None] >> np.arange(N)) & 1).astype(float)
    y = v[masks] - v[0]
    w = np.array([kernel(N, int(m)) for m in X.sum(axis=1)])
    vN = v[FULL] - v[0]
    # eliminate phi_{n-1}: x_i - x_{n-1}, rhs y - x_{n-1} * vN
    A = X[:, :-1] - X[:, -1:]
    b = y - X[:, -1] * vN
    sw = np.sqrt(w)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], b * sw, rcond=None)
    return np.append(coef, vN - coef.sum())


def plain_perm(v, budget, seed):
    rng = random.Random(seed)
    seen = np.zeros(1 << N, dtype=bool)
    seen[0] = True
    used = 0
    acc = [[] for _ in range(N)]
    while used < budget:
        o = list(range(N))
        rng.shuffle(o)
        S = 0
        for i in o:
            prev = S
            S |= 1 << i
            if not seen[S]:
                seen[S] = True
                used += 1
            acc[i].append(v[S] - v[prev])
    return np.array([np.mean(a) for a in acc])


def main():
    L, add = [], None
    add = L.append
    add("Baselines on the city n=16 games, same exact oracle and currency")
    add("error = max_i |est - phi| / max_i |phi|, median of 24 seeds")
    add("")
    for label in sc.LABELS:
        for col, mname in sc.METRICS:
            v = sc.load_table(label, col)
            v = v - v[0]
            true = sc.exact_shapley(v, N)
            ref = np.abs(true).max()
            add("=" * 76)
            add("%s   metric = %s   max|phi| = %.4f" % (label, mname, ref))
            add("=" * 76)
            add("  %7s %12s %12s %12s" % ("budget", "KernelSHAP", "perm", "perm+anti"))
            for B in BUDGETS:
                k = np.median([100 * np.abs(kernelshap(v, B, 500 + s) - true).max() / ref
                               for s in range(sc.NSEED)])
                p = np.median([100 * np.abs(plain_perm(v, B, 700 + s) - true).max() / ref
                               for s in range(sc.NSEED)])
                a = np.median([100 * np.abs(sc.perm_anti(v, N, B, 1000 + s)[0] - true).max() / ref
                               for s in range(sc.NSEED)])
                add("  %7d %11.2f%% %11.2f%% %11.2f%%" % (B, k, p, a))
            add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
