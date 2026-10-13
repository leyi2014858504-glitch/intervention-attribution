#!/usr/bin/env python3
"""Cumulative spectral-mass share R_q = sum_{k<=q} W(k) / sum_k W(k).

W(k) = sum_{|A|=k} |m(A)| / |A|   (Shapley-weighted spectral mass, k=1 normalised)
R_q  = fraction of the absolute-mass total captured by orders <= q.

Zero chemistry runs: released oracle tables only.
  n=16  _kppwork/n16_v_<scen>_m8.csv    (v_ppb = d3_mean, v_mda8)
  n=12  _kppwork/co_v.csv (d3_mean), _kppwork/co_v_m8_*.csv (MDA8)
  n=8   _kppwork/scen_v.csv (d3_mean), sub-lattice of the n=16 MDA8 table
Also prints the exact W(k)/W(1) at k=3,4,5 and their sum-share (the
"0.4-0.6 of k=1" phrase check).  Writes results_cum_share.txt.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCEN = ("r_sat", "r_lim1")
NP = {8: 256, 12: 4096, 16: 65536}


def mobius_fast(vals, n):
    f = list(vals)
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                f[mask] -= f[mask ^ bit]
    return f


def embed8(m8):
    """4-window (6 h) mask -> 8-window (3 h) mask (each 6 h block = its two
    3 h halves switched together)."""
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def read_scen(path, col):
    """-> {scen: {mask: v}} from a CSV with name '<scen>_m<mask>'."""
    out = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh)
        if col not in (rdr.fieldnames or []):
            raise SystemExit("!! %s has no column %s (has %s)"
                             % (path, col, rdr.fieldnames))
        for row in rdr:
            name = row.get("name", "")
            if "_m" not in name:
                continue
            s, _, ms = name.rpartition("_m")
            if s not in SCEN:
                continue
            try:
                mask = int(ms)
            except ValueError:
                continue
            v = row.get(col, "")
            if v in ("", "NaN"):
                continue
            out.setdefault(s, {})[mask] = float(v)
    return out


def read_n16(scen, col):
    v = {}
    path = os.path.join(HERE, "n16_v_%s_m8.csv" % scen)
    with open(path) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val in ("", "NaN"):
                continue
            v[int(row["name"][1:])] = float(val)
    return v


def mass_and_share(vals, n):
    m = mobius_fast(vals, n)
    w = [0.0] * (n + 1)
    for mask in range(1, 1 << n):
        k = bin(mask).count("1")
        w[k] += abs(m[mask]) / k
    ww = [x / w[1] for x in w[1:]]          # normalised to k=1
    tot = sum(ww)
    R, c = [], 0.0
    for x in ww:
        c += x
        R.append(c / tot)
    return ww, R


def main():
    L = []
    add = L.append
    add("Cumulative spectral-mass share  R_q = sum_{k<=q} W(k) / sum_k W(k)")
    add("W(k) = sum_{|A|=k} |m(A)|/|A|  (Shapley-weighted spectral mass, k=1 normalised)")
    add("Zero chemistry runs: released oracle tables only.")
    for metric, col, src in (("MDA8", "v_mda8", "m8"), ("d3_mean", "v_ppb", "d3")):
        add("")
        add("=" * 78)
        add("metric = %s" % metric)
        add("=" * 78)
        for scen in SCEN:
            v16 = read_n16(scen, col)
            if src == "m8":
                v12 = {}
                for K in range(1, 7):
                    p = os.path.join(HERE, "co_v_m8_%d.csv" % K)
                    if os.path.exists(p):
                        v12.update(read_scen(p, col).get(scen, {}))
                v8 = {m: v16[embed8(m)] for m in range(NP[8])}
            else:
                v12 = read_scen(os.path.join(HERE, "co_v.csv"), col).get(scen, {})
                v8 = read_scen(os.path.join(HERE, "scen_v.csv"), col).get(scen, {})
                dmax = max(abs(v16[embed8(m)] - v8[m]) for m in range(NP[8]))
                add("  [%s] n=8 sub-lattice crosscheck vs scen_v.csv: max|d| = %.2e"
                    % (scen, dmax))
            add("")
            add("scenario %s" % scen)
            for n in (8, 12, 16):
                vd = {8: v8, 12: v12, 16: v16}[n]
                assert len(vd) == NP[n], "n=%d incomplete: %d of %d" % (n, len(vd), NP[n])
                vals = [vd[m] for m in range(NP[n])]
                ww, R = mass_and_share(vals, n)
                share35 = sum(ww[2:5]) / sum(ww)   # k=3..5 as fraction of total
                add("  n=%2d  W(k)/W(1):" % n)
                add("        " + " ".join("k%02d=%.3f" % (k + 1, x) for k, x in enumerate(ww)))
                add("      R_q: " + " ".join("q%d=%.3f" % (q + 1, r) for q, r in enumerate(R)))
                add("      k=3..5 of k=1: %.3f %.3f %.3f   (range %.3f-%.3f; %.1f%% of total)"
                    % (ww[2], ww[3], ww[4], min(ww[2:5]), max(ww[2:5]), 100 * share35))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(ROOT, "results_cum_share.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote results_cum_share.txt")


if __name__ == "__main__":
    main()