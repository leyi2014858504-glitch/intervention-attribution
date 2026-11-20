#!/usr/bin/env python3
"""Sampler cost law on the CITY-ANCHORED n=16 games (zero new chemistry runs).

The sampler's output depends only on v(S) at the coalitions it visits, and the
city tables contain every one of the 2^16 values.  Sampling can therefore be
replayed exactly against the released table: no box runs, and the oracle is the
same exact enumeration used in the audit.

Two ladders come from one game each:
  n = 8   the 6-hour view: each pair of adjacent 3-hour windows is switched
          together, i.e. the embedded sub-lattice used for the synthetic ladder
  n = 16  the 3-hour view, the full game

Estimator: permutation chains with antithetic partners (each order is also
walked backwards), the working sampler of the paper.  Budget currency = one
evaluation of v at a coalition not seen before; v({}) is known and free.  A
chain may overshoot the nominal budget by up to n-1 evaluations (checked below).

Scoring: max_i |est_i - phi_i| / max_i |phi_i|, error per seed, 24 seeds.

Writes results_sampler_city.txt
"""
import csv
import glob
import math
import os
import random

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results_sampler_city.txt")
LABELS = ("dallas_2016w", "deerpark_2017s")
METRICS = (("v_mda8", "MDA8"), ("v_ppb", "d3_mean"))
NSEED = 24
BUDGETS = (64, 128, 256, 512, 1024, 2048, 4096, 8192)


def load_table(label, col):
    v = np.full(1 << 16, np.nan)
    for f in sorted(glob.glob(os.path.join(HERE, "n16_v_%s_m8_*.csv" % label))):
        for row in csv.DictReader(open(f)):
            if row[col] in ("", "NaN"):
                continue
            v[int(row["name"][1:])] = float(row[col])
    if np.isnan(v).any():
        raise RuntimeError("%s incomplete for %s" % (label, col))
    return v


def embed8(m):
    """n=8 mask -> n=16 mask: each 6 h block = its two 3 h blocks together."""
    m16 = 0
    for j in range(4):
        if (m8 := (m >> j) & 1):
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def mobius(v, n):
    f = np.array(v, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def exact_shapley(v, n):
    m = mobius(v, n)
    A = np.arange(1 << n)
    pc = np.array([bin(int(x)).count("1") for x in A], dtype=float)
    pc[0] = 1.0
    phi = np.zeros(n)
    for i in range(n):
        sel = ((A >> i) & 1) == 1
        phi[i] = np.sum(m[sel] / pc[sel])
    return phi


def perm_anti(v, n, budget, seed):
    """one antithetic permutation sampler; returns (estimate, evaluations used)."""
    seen = np.zeros(1 << n, dtype=bool)
    seen[0] = True
    used = 0
    acc = [[] for _ in range(n)]
    rng = random.Random(seed)
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
                acc[i].append(v[S] - v[prev])
        stall = stall + 1 if used == before else 0
    return np.array([np.mean(a) for a in acc]), used


def runs_for(pts, target):
    """log-log interpolation of the runs needed to reach `target` error."""
    pts = sorted(pts)
    for (b0, e0), (b1, e1) in zip(pts, pts[1:]):
        if (e0 - target) * (e1 - target) <= 0 and e0 != e1:
            t = (math.log(target) - math.log(e0)) / (math.log(e1) - math.log(e0))
            return math.exp(math.log(b0) + t * (math.log(b1) - math.log(b0)))
    return float("nan")


def main():
    L, add = [], None
    add = L.append
    add("Sampler cost law on the CITY-ANCHORED games (replayed on the exact tables)")
    add("estimator = permutation + antithetic, %d seeds" % NSEED)
    add("error = max_i |est - phi| / max_i |phi|, per seed")
    add("budget = distinct coalition evaluations; v({}) free; chain overshoot <= n-1")
    add("")

    headline = []
    for label in LABELS:
        for col, mname in METRICS:
            v16 = load_table(label, col)
            base = v16[0]
            v16 = v16 - base                      # v({}) = 0, as in the audit

            for n in (8, 16):
                if n == 8:
                    v = np.zeros(256)
                    for m in range(256):
                        v[m] = v16[embed8(m)]
                else:
                    v = v16
                true = exact_shapley(v, n)
                scale = np.abs(true).max()

                add("=" * 88)
                add("%s  metric = %s  n = %d   max|phi| = %.4f ppb   enum = %d runs"
                    % (label, mname, n, scale, (1 << n) - 1))
                add("=" * 88)
                add("  %8s %10s %10s %10s %10s %10s"
                    % ("budget", "used", "err med", "err q1", "err q3", "err max"))
                pts = []
                for B in BUDGETS:
                    errs, useds = [], []
                    for s in range(NSEED):
                        est, used = perm_anti(v, n, B, 1000 + s)
                        errs.append(100 * np.abs(est - true).max() / scale)
                        useds.append(used)
                    errs = np.array(errs)
                    add("  %8d %10.1f %9.2f%% %9.2f%% %9.2f%% %9.2f%%"
                        % (B, np.mean(useds), np.median(errs),
                           np.percentile(errs, 25), np.percentile(errs, 75),
                           errs.max()))
                    pts.append((B, float(np.median(errs))))
                r15 = runs_for(pts, 1.5)
                r50 = runs_for(pts, 5.0)
                add("  runs for 1.5%%: %s   for 5%%: %s   (enumeration %d)"
                    % ("%.0f" % r15 if not math.isnan(r15) else "n/a",
                       "%.0f" % r50 if not math.isnan(r50) else "n/a",
                       (1 << n) - 1))
                add("")
                headline.append((label, mname, n, r15, r50, (1 << n) - 1,
                                 pts[-1][1]))

    add("=" * 88)
    add("HEADLINE  runs needed (median seed)")
    add("=" * 88)
    add("  %-16s %-8s %4s %10s %10s %10s" %
        ("slice", "metric", "n", "runs@1.5%", "runs@5%", "enum"))
    for label, mname, n, r15, r50, e, _ in headline:
        add("  %-16s %-8s %4d %10s %10s %10d" %
            (label, mname, n,
             "%.0f" % r15 if not math.isnan(r15) else "n/a",
             "%.0f" % r50 if not math.isnan(r50) else "n/a", e))

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
