#!/usr/bin/env python3
"""How far do the Owen values sit from Shapley, and does the distance grow with
partition coarseness?

Partitions: singletons (Shapley), 8 groups of 2 (time-first), 4 groups of 4
(sector x half-day), 2 groups of 8 (species-first).

Writes results_kernel_spread.txt
"""
import csv
import glob
import os
from math import comb, factorial

import numpy as np

D = os.environ.get("ATMOS_DIR", "/mnt/d/atmos")
W = os.path.join(D, "_kppwork")
N = 16
IDX = np.arange(1 << N)
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")
LAB = [("NOx", w) for w in WIN] + [("VOC", w) for w in WIN]

PARTS = {
    "time-first 8x2": [frozenset([w, 8 + w]) for w in range(8)],
    "sector-halfday 4x4": [frozenset([0, 1, 2, 3]), frozenset([4, 5, 6, 7]),
                           frozenset([8, 9, 10, 11]), frozenset([12, 13, 14, 15])],
    "species-first 2x8": [frozenset(range(8)), frozenset(range(8, 16))],
}


def load(tag):
    pat = ("n16_vh_%s_m8_*.csv" % tag) if tag in ("r_sat", "r_lim1") \
        else ("n16_vh_%s_m8_[0-9]*.csv" % tag)
    v = np.full(1 << N, np.nan)
    for p in sorted(glob.glob(os.path.join(W, pat))):
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get("v_mda8", "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    return v


def mobius(v):
    f = v.copy()
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        f[sel] -= f[IDX[sel] ^ bit]
    return f


def shapley(m):
    k = np.array([bin(x).count("1") for x in IDX])
    w = np.where(k > 0, m / np.maximum(k, 1), 0.0)
    return np.array([w[(IDX & (1 << i)) > 0].sum() for i in range(N)])


def S_coef(k, p):
    s = 0.0
    for t in range(k - 1, p):
        s += comb(p - k, t - k + 1) * factorial(t) * factorial(p - t - 1) / factorial(p)
    return s


def G_coef(g, q):
    s = 0.0
    for r in range(g - 1, q):
        s += comb(q - g, r - g + 1) * factorial(r) * factorial(q - r - 1) / factorial(q)
    return s


def owen(m, groups):
    q = len(groups)
    grp, p_of = {}, {}
    for gi, g in enumerate(groups):
        for i in g:
            grp[i] = gi
            p_of[i] = len(g)
    Sc = {(k, p): S_coef(k, p) for p in (2, 4, 8) for k in range(1, p + 1)}
    Gc = {g: G_coef(g, q) for g in range(1, q + 1)}
    psi = np.zeros(N)
    for mask in IDX[1:]:
        mv = m[mask]
        if mv == 0.0:
            continue
        mem = [i for i in range(N) if mask >> i & 1]
        gc = Gc[len({grp[i] for i in mem})]
        per = {}
        for i in mem:
            per[grp[i]] = per.get(grp[i], 0) + 1
        for i in mem:
            psi[i] += mv * Sc[(per[grp[i]], p_of[i])] * gc
    return psi


L = []
add = L.append
add("Distance of each Owen partition from Shapley (MDA8)")
add("d = max_i |psi_i - phi_i| / max_i |phi_i|")
add("%-16s %10s %10s %10s %10s" %
    ("game", "time 8x2", "4x4", "species 2x8", "max|phi|"))
for g in GAMES:
    v = load(g)
    m = mobius(v)
    phi = shapley(m)
    sc = np.abs(phi).max()
    ds = []
    for name in ("time-first 8x2", "sector-halfday 4x4", "species-first 2x8"):
        psi = owen(m, PARTS[name])
        ds.append(100 * np.abs(psi - phi).max() / sc)
    add("%-16s %9.1f%% %9.1f%% %9.1f%% %10.4f" % (g, ds[0], ds[1], ds[2], sc))

add("")
add("Largest player movements between the two reported partitions (ppb)")
for g in GAMES:
    v = load(g)
    m = mobius(v)
    pt = owen(m, PARTS["time-first 8x2"])
    ps = owen(m, PARTS["species-first 2x8"])
    d = ps - pt
    order = np.argsort(-np.abs(d))[:4]
    add("  %-16s " % g + "  ".join("%s@%s %+.2f" % (LAB[i][0], LAB[i][1], d[i])
                                  for i in order))

add("")
add("Sector totals (ppb): time-first vs species-first")
for g in GAMES:
    v = load(g)
    m = mobius(v)
    pt = owen(m, PARTS["time-first 8x2"])
    ps = owen(m, PARTS["species-first 2x8"])
    add("  %-16s NOx %+7.3f / %+7.3f   VOC %+7.3f / %+7.3f   shift %5.2f"
        % (g, pt[:8].sum(), ps[:8].sum(), pt[8:].sum(), ps[8:].sum(),
           abs(pt[:8].sum() - ps[:8].sum())))

add("")
add("What moves: sector split vs window totals (ppb), time-first -> species-first")
add("%-16s %12s %14s" % ("game", "max sector", "max window total"))
for g in GAMES:
    v = load(g)
    m = mobius(v)
    pt = owen(m, PARTS["time-first 8x2"])
    ps = owen(m, PARTS["species-first 2x8"])
    sec = max(abs(pt[:8].sum() - ps[:8].sum()), abs(pt[8:].sum() - ps[8:].sum()))
    win = max(abs((pt[w] + pt[8 + w]) - (ps[w] + ps[8 + w])) for w in range(8))
    add("%-16s %+12.3f %+14.3f" % (g, sec, win))

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_kernel_spread.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_kernel_spread.txt")
