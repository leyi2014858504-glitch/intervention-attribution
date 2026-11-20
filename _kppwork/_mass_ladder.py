#!/usr/bin/env python3
"""Bare-mass ladder n=8 (6-h sub-lattice of n=16), n=12 (co table), n=16.

Prints M(k)/M(1), R_q (q=1,2,3,5) and the k=3..5 share of total mass, matching
the definitions in the manuscript (M(k) = sum_{|A|=k}|m(A)|).
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def mobius(v, n):
    f = v.copy()
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def mass(v, n):
    m = mobius(v, n)
    out = np.zeros(n + 1)
    for mask in range(1, 1 << n):
        out[bin(mask).count("1")] += abs(m[mask])
    return out / out[1]


def embed8(m8):
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def load_n16(scen):
    v = np.zeros(1 << 16)
    with open(os.path.join(HERE, "n16_v_%s_m8.csv" % scen)) as fh:
        for row in csv.DictReader(fh):
            if row["v_mda8"] not in ("", "NaN"):
                v[int(row["name"][1:])] = float(row["v_mda8"])
    return v


def load_n12(scen):
    v = np.zeros(1 << 12)
    for K in range(1, 7):
        with open(os.path.join(HERE, "co_v_m8_%d.csv" % K)) as fh:
            for row in csv.DictReader(fh):
                name = row["name"]
                if "_m" not in name or "v_mda8" not in row:
                    continue
                s, ms = name.rsplit("_m", 1)
                if s == scen and row["v_mda8"] not in ("", "NaN"):
                    v[int(ms)] = float(row["v_mda8"])
    return v


for scen in ("r_sat", "r_lim1"):
    v16 = load_n16(scen)
    v12 = load_n12(scen)
    v8 = np.array([v16[embed8(m)] for m in range(256)])
    print("=" * 70)
    print(scen)
    for v, n in ((v8, 8), (v12, 12), (v16, 16)):
        mk = mass(v, n)
        tot = mk[1:].sum()
        r = [mk[1:q + 1].sum() / tot for q in (1, 2, 3, 5)]
        print("  n=%2d  M(k)/M(1): %s" %
              (n, " ".join("%.3f" % x for x in mk[1:])))
        print("        R1=%.3f R2=%.3f R3=%.3f R5=%.3f   k3-5 share=%.1f%%   peak k=%d (%.2fx)" %
              (r[0], r[1], r[2], r[3], 100 * mk[3:6].sum() / tot,
               1 + int(np.argmax(mk[1:])), mk[1:].max()))
