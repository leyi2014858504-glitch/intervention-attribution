#!/usr/bin/env python3
"""Truncation divergence and spectral mass under the MDA8 metric.

Same order-decomposition analysis as _n16_spectrum.py, run on the v_mda8
column of the MDA8 campaign tables.  Purpose: show that claim (1) -- the
truncation curve is non-monotone and can diverge, and the spectral mass
peaks at high order -- does not depend on the choice of metric.

n = 8 is NOT re-run.  The 6 h game is an exact sub-lattice of the 16-player
game: turning off 6 h window j == turning off 3 h windows 2j and 2j+1, i.e.
the SAME forcing history, hence the SAME trajectory, hence ANY functional of
it agrees.  So the n=8 MDA8 values are read off the n=16 table at
union-of-pairs masks.  This derivation is licensed by first checking the
identity numerically in v_ppb against scen_v.csv (the gate7 sub-lattice
check verified it to ~3.6e-14 across backends; we re-verify per run).

n = 12 comes from the co_v_m8 shards (3 sectors x 4 windows, kpp6).

Writes spectrum_orders_m8.txt.
"""
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = HERE
OUT = os.path.join(W, "spectrum_orders_m8.txt")

_spec = importlib.util.spec_from_file_location("expspec", os.path.join(HERE, "_n16_spectrum.py"))
_ex = importlib.util.module_from_spec(_spec)
sys.argv = [sys.argv[0]]  # _n16_spectrum reads argv? it does not; keep guard
_spec.loader.exec_module(_ex)
analyse = _ex.analyse

SCENS = ("r_sat", "r_lim1")
TOL = 1e-9  # sub-lattice cross-check tolerance on v_ppb


def load_n16_m8(scen):
    """merged n16_v_<scen>_m8.csv -> {mask: (v48, v8)}"""
    fn = os.path.join(W, "n16_v_%s_m8.csv" % scen)
    d = {}
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[int(row["name"][1:])] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def embed8(m8):
    """6 h mask (NOx bits 0-3, VOC bits 4-7) -> 3 h mask (NOx 0-7, VOC 8-15)."""
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:            # NOx 6 h window j off
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:      # VOC 6 h window j off
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def derive_n8(v16, scen):
    """n=8 table from the n=16 MDA8 table + sub-lattice cross-check in v48."""
    ref = {}
    with open(os.path.join(W, "scen_v.csv")) as fh:
        for row in csv.DictReader(fh):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if s == scen and row["v_ppb"] not in ("", "NaN"):
                ref[int(ms)] = float(row["v_ppb"])
    v48, v8 = {}, {}
    worst = 0.0
    for m8 in range(256):
        p48, p8 = v16[embed8(m8)]
        v48[m8], v8[m8] = p48, p8
        if m8 in ref:
            worst = max(worst, abs(p48 - ref[m8]))
    flag = "OK" if worst < TOL else "!! EXCEEDS TOLERANCE -- derivation unsafe"
    print("[%s] sub-lattice check (n=8 derived vs scen_v.csv): max|dv48| = %.3e  %s"
          % (scen, worst, flag))
    return v8


def load_co_m8(scen):
    """co_v_m8_*.csv shards -> {mask: v_mda8} + v48 sub-lattice check vs co_v.csv"""
    v48, v8 = {}, {}
    for K in range(1, 7):
        fn = os.path.join(W, "co_v_m8_%d.csv" % K)
        if not os.path.exists(fn):
            continue
        with open(fn) as fh:
            for row in csv.DictReader(fh):
                nm = row["name"]
                if not nm.startswith(scen + "_m"):
                    continue
                m = int(nm.rsplit("_m", 1)[1])
                if row["v_ppb"] not in ("", "NaN"):
                    v48[m] = float(row["v_ppb"])
                if row["v_mda8"] not in ("", "NaN"):
                    v8[m] = float(row["v_mda8"])
    ref = {}
    with open(os.path.join(W, "co_v.csv")) as fh:
        for row in csv.DictReader(fh):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if s == scen and row["v_ppb"] not in ("", "NaN"):
                ref[int(ms)] = float(row["v_ppb"])
    worst = max((abs(v48[m] - ref[m]) for m in ref if m in v48), default=0.0)
    flag = "OK" if worst < TOL else "!! EXCEEDS TOLERANCE"
    print("[%s] co12 v48 regression vs co_v.csv: max|dv48| = %.3e  %s (%d of %d)"
          % (scen, worst, flag, len(v8), 1 << 12))
    return v8


def main():
    L = []
    add = L.append
    add("Order decomposition under the MDA8 metric (v_mda8 column)")
    add("Same definitions as spectrum_orders.txt: trunc(q) normalised by max|phi|;")
    add("mass by order normalised to k=1.")
    add("n=8 derived from the n=16 table via the 6h sub-lattice (checked in v48).")
    add("")
    for scen in SCENS:
        add("=" * 78)
        add("scenario %s   (metric = MDA8)" % scen)
        add("=" * 78)
        v16 = load_n16_m8(scen)
        if len(v16) != 1 << 16:
            add("n=16 incomplete (%d of %d) -- abort this scenario" % (len(v16), 1 << 16))
            continue
        tabs = {}
        v8n = derive_n8(v16, scen)
        if len(v8n) == 256:
            tabs[8] = analyse(v8n, 8, 8)
        vc = load_co_m8(scen)
        if len(vc) == 1 << 12:
            tabs[12] = analyse(vc, 12, 12)
        tabs[16] = analyse({m: p[1] for m, p in v16.items()}, 16, 16)

        add("%-6s %9s %9s %9s %9s %9s %9s" %
            ("n", "max|phi|", "trunc(1)", "trunc(2)", "trunc(3)", "trunc(4)", "trunc(6)"))
        for n in sorted(tabs):
            true, trunc, mass, scale = tabs[n]
            add("%-6d %9.2f %8.1f%% %8.1f%% %8.1f%% %8.1f%% %8.1f%%" % (
                n, scale,
                100 * trunc.get(1, float("nan")),
                100 * trunc.get(2, float("nan")),
                100 * trunc.get(3, float("nan")),
                100 * trunc.get(4, float("nan")),
                100 * trunc.get(6, float("nan"))))
        add("")
        add("Irreducible floor per order:")
        for n in sorted(tabs):
            true, trunc, mass, scale = tabs[n]
            vals = ", ".join("q<=%d: %.1f%%" % (q, 100 * trunc[q])
                             for q in sorted(trunc)[:8])
            add("  n=%-3d %s" % (n, vals))
        add("")
        add("Spectral mass by order, normalised to k=1:")
        for n in sorted(tabs):
            true, trunc, mass, scale = tabs[n]
            base = mass[1] if mass[1] else 1.0
            add("  n=%-3d " % n + " ".join("k=%d:%.3f" % (k, mass[k] / base)
                                          for k in range(1, min(n, 12) + 1)))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
