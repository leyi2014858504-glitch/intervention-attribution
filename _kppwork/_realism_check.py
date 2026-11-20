#!/usr/bin/env python3
"""Is the synthetic pair inside the range spanned by the two anchored cities?

Compares, for each n=16 game: emission amplitudes, the implied loading scale
E*tau (mean profile 0.8), the resulting all-on MDA8 / d3_mean, the response
fraction, and the sector split of the first-order Möbius mass.

Writes results_realism_check.txt
"""
import csv
import os

import numpy as np

W = "/mnt/d/atmos/_kppwork"
N = 16
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
CFG = {
    "r_sat":      dict(NOX=16.0, VOC=12.0, ISO=1.0, CO=40.0, KDIL=1500.0,
                       BGO3=30.0, BGNO2=5.0, BGCO=120.0),
    "r_lim1":     dict(NOX=8.0, VOC=18.0, ISO=1.0, CO=40.0, KDIL=1500.0,
                       BGO3=30.0, BGNO2=5.0, BGCO=120.0),
    "dallas_2016w": dict(NOX=2.7312, VOC=0.50, ISO=0.0, CO=1.0, KDIL=2288.0,
                         BGO3=7.0625, BGNO2=2.0, BGCO=222.1875),
    "deerpark_2017s": dict(NOX=0.50, VOC=4.2406, ISO=10.0, CO=1.0, KDIL=2317.09,
                           BGO3=5.0, BGNO2=2.0, BGCO=143.4375),
}
IDX = np.arange(1 << N)
K = np.array([bin(x).count("1") for x in IDX])


def load(game, col):
    v = np.zeros(1 << N)
    if game in ("r_sat", "r_lim1"):
        files = [os.path.join(W, "n16_v_%s_m8.csv" % game)]
    else:
        files = [os.path.join(W, "n16_v_%s_m8_%d.csv" % (game, j)) for j in range(1, 13)]
    for p in files:
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get(col, "")
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


L = []
add = L.append
add("Synthetic vs anchored: what actually differs, and where each sits")
add("L_NOx = 0.796 * E_NOx * tau_day  [ppb]  (source-exchange indicator;")
add("not ambient NOx and not a full-day steady state)")
add("hdr weights: O3 1.0, NOx 1.0, CO 0.25, PAR 0.20, each normalised by its daily mean")
add("")
hdr = ("%-16s %6s %6s %6s %7s %7s %7s %6s %7s %8s %8s %7s" %
       ("game", "NOx", "VOC", "ISO", "CO", "KDIL", "BGO3", "V/N", "L_NOx", "MDA8", "d3", "resp%"))
add(hdr)
for g in GAMES:
    c = CFG[g]
    ratio = c["VOC"] / c["NOX"]
    loadscale = 0.8 * c["NOX"] * c["KDIL"] / 3600.0
    v8 = load(g, "v_mda8")
    vd = load(g, "v_ppb")
    resp = 100 * (v8[(1 << N) - 1] - v8[0]) / np.abs(v8).max()
    add("%-16s %6.2f %6.2f %6.2f %7.1f %7.0f %7.1f %6.2f %6.2f %8.2f %8.2f %7.1f" %
        (g, c["NOX"], c["VOC"], c["ISO"], c["CO"], c["KDIL"], c["BGO3"],
         ratio, loadscale, v8[0], vd[0], resp))

add("")
add("Sector composition of the first-order mass (MDA8, ppb) and mixed share")
add("%-16s %10s %10s %12s %10s" % ("game", "|m1|NOx", "|m1|VOC", "VOC share", "mixed share"))
for g in GAMES:
    v = load(g, "v_mda8")
    m = mobius(v)
    m1n = np.abs([m[1 << i] for i in range(8)]).sum()
    m1v = np.abs([m[1 << (8 + i)] for i in range(8)]).sum()
    mixed = np.abs(m[(K >= 2) & ((IDX & 0xFF) > 0) & ((IDX >> 8) > 0)]).sum()
    total = np.abs(m[1:]).sum()
    add("%-16s %10.2f %10.3f %11.1f%% %9.1f%%" %
        (g, m1n, m1v, 100 * m1v / (m1n + m1v), 100 * mixed / total))

txt = "\n".join(L)
print(txt)
with open("/mnt/d/atmos/results_realism_check.txt", "w") as fh:
    fh.write(txt + "\n")
print("wrote results_realism_check.txt")
