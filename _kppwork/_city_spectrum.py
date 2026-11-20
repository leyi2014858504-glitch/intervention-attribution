#!/usr/bin/env python3
"""Spectral mass by order for city-anchored slices (post-processing, zero runs).

Reads the n=16 MDA8 oracle tables for dallas_2016w and deerpark_2017s,
computes the Mobius spectrum, and reports mass(k) = sum_{|A|=k} |m(A)|/|A|
normalised to mass(1), plus the cumulative share R_q = sum_{k<=q} mass(k) /
sum_{k<=n} mass(k).  Also reports the top-5 dividends per scenario.

Writes results_city_spectrum.txt
"""
import csv
import os
import numpy as np

import platform
if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"
LABELS = ["dallas_2016w", "deerpark_2017s"]


def mobius_np(vals, n):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def mass_by_order(m, n):
    mass = np.zeros(n + 1)
    idx = np.arange(1, 1 << n)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = np.abs(m[idx]) / k
    for ki in range(1, n + 1):
        mass[ki] = w[k == ki].sum()
    return mass


def read_n16(label, col):
    v = np.zeros(65536)
    for K in range(1, 13):
        p = os.path.join(W, "n16_v_%s_m8_%d.csv" % (label, K))
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh):
                m = int(row["name"][1:])
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[m] = float(val)
    return v


def top_dividends(m, n, kmax=5):
    idx = np.arange(1, 1 << n)
    k = np.array([bin(x).count("1") for x in idx])
    order = np.argsort(-np.abs(m[idx]))
    out = []
    for j in order[:kmax]:
        mask = idx[j]
        out.append((int(mask), int(k[j]), float(m[mask])))
    return out


def mask_to_str(mask, nwin=8):
    nox = [j for j in range(nwin) if (mask >> j) & 1]
    voc = [j for j in range(nwin) if (mask >> (nwin + j)) & 1]
    wl = ["0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24"]
    ns = "+".join(wl[j] for j in nox) if nox else "-"
    vs = "+".join(wl[j] for j in voc) if voc else "-"
    return "NOx[%s] VOC[%s]" % (ns, vs)


L = []
add = L.append
add("City-anchored spectral mass by order (n=16, MDA8)")
add("mass(k) = sum_{|A|=k} |m(A)|/|A|  normalised to mass(1)")
add("R_q = cumulative share = sum_{k<=q} mass(k) / sum_{k<=n} mass(k)")
add("")

for label in LABELS:
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3_mean")):
        v = read_n16(label, col)
        m = mobius_np(v, 16)
        mass = mass_by_order(m, 16)
        r = mass / mass[1]
        cum = np.cumsum(mass[1:])
        total = cum[-1]
        Rq = cum / total
        add("=== %s  %s ===" % (label, tag))
        add("  mass(k)/mass(1): " + "  ".join("k%d=%.3f" % (k, r[k]) for k in range(1, 17)))
        add("  R_q cumulative:  " + "  ".join("q%d=%.3f" % (q, Rq[q - 1]) for q in range(1, 17)))
        add("  peak k = %d  (mass/mass1 = %.3f)" % (int(np.argmax(r[1:])) + 1, r[1:].max()))
        add("")
        add("  top-5 |m(A)|:")
        for mask, k, mv in top_dividends(m, 16, 5):
            add("    k=%-3d m=%+9.3f  %s" % (k, mv, mask_to_str(mask)))
        add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_city_spectrum.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_city_spectrum.txt")
