#!/usr/bin/env python3
"""Why is Deer Park's spectrum flat?  Sector-resolved diagnostics on the v3 table.

Prints, for each city and MDA8:
  * the 16 first-order coefficients m({i}) (ppb), NOx windows then VOC windows
  * the single-switch effect v({i}) - v(0)
  * order mass split by whether A contains only NOx players, only VOC players,
    or both
"""
import csv
import os
import sys

import numpy as np

W = "/mnt/d/atmos/_kppwork"
N = 16
IDX = np.arange(1 << N)
K = np.array([bin(x).count("1") for x in IDX])


def load(label):
    v = np.zeros(1 << N)
    for j in range(1, 13):
        p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (label, j))
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


for label in ("dallas_2016w", "deerpark_2017s"):
    v = load(label)
    m = mobius(v)
    print("=" * 78)
    print("%s   v(all-on)=%.3f  v(all-off)=%.3f  min=%.3f  max=%.3f"
          % (label, v[0], v[(1 << N) - 1], v.min(), v.max()))
    first = np.array([m[1 << i] for i in range(N)])
    print("  first order  NOx: " + " ".join("%6.2f" % x for x in first[:8]))
    print("  first order  VOC: " + " ".join("%6.2f" % x for x in first[8:]))
    print("  |m1| sum NOx %.3f   VOC %.3f" % (np.abs(first[:8]).sum(), np.abs(first[8:]).sum()))
    nox = IDX & 0x00FF == IDX & 0x00FF
    only_nox = (IDX >> 8) == 0
    only_voc = (IDX & 0xFF) == 0
    print("  %4s %10s %10s %10s" % ("k", "NOx-only", "VOC-only", "mixed"))
    for k in range(1, N + 1):
        s = K == k
        a = np.abs(m[s & only_nox]).sum()
        b = np.abs(m[s & only_voc]).sum()
        c = np.abs(m[s & ~only_nox & ~only_voc]).sum()
        print("  %4d %10.3f %10.3f %10.3f" % (k, a, b, c))
