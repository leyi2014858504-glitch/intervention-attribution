#!/usr/bin/env python3
"""
Order decomposition of phi at n = 8, 12, 16 from the FULLY ENUMERATED oracles.

WHY
---
Claim (1) of the paper is that low-order truncation of the interaction spectrum
cannot work, so the spectrum must be measured rather than assumed.  That claim
was previously supported at n=8 and n=12 only.  The n=16 enumeration (65536
coalitions, exact) lets us put a third, much larger point on it -- and this is a
statement about the ATTRIBUTION ARCHITECTURE, not about any sampler.

For a coalition game with Moebius coefficients m(A),

    phi_i = sum_{A containing i} m(A)/|A|

so splitting the sum by |A| gives an exact order decomposition of every Shapley
value, and the truncation error at order <= q is measurable with no sampling at
all.  Two numbers come out of it:

  * truncation curve  : max_i |phi_i - phi_i^{<=q}| / max_i |phi_i|
                        the best a perfect q-th-order model could ever do
  * spectral mass     : sum_{|A|=k} |m(A)| / |A|   by order
                        shows how the interaction weight is distributed

Writes spectrum_orders.txt
"""
import csv
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W = "/mnt/d/atmos/_kppwork"
OUT = "/mnt/d/atmos/_kppwork/spectrum_orders.txt"

_spec = importlib.util.spec_from_file_location("expn16", os.path.join(HERE, "exp_n16.py"))
_ex16 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ex16)
mobius_fast = _ex16.mobius_fast

SCEN = "r_sat"


def load_simple(fn, scen):
    v = {}
    with open(os.path.join(HERE, fn)) as fh:
        for row in csv.DictReader(fh):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if s != scen or row["v_ppb"] in ("", "NaN"):
                continue
            v[int(ms)] = float(row["v_ppb"])
    return v


def load_n16(scen):
    v = {}
    for K in range(1, 13):
        fn = os.path.join(W, "n16_v_%s_%d.csv" % (scen, K))
        if not os.path.exists(fn):
            continue
        with open(fn) as fh:
            for row in csv.DictReader(fh):
                if row.get("name", "").startswith("m") and row["v_ppb"] not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(row["v_ppb"])
    return v


def analyse(v, n, tag):
    """returns (true, trunc table, mass by order)"""
    vals = [v[m] for m in range(1 << n)]
    mo = mobius_fast(vals, n)
    pc = [bin(A).count("1") for A in range(1 << n)]

    # exact phi, and the same split by order
    true = [0.0] * n
    part = [{q: 0.0 for q in range(n + 1)} for _ in range(n)]
    mass = [0.0] * (n + 1)
    for A in range(1, 1 << n):
        k = pc[A]
        c = mo[A] / k
        mass[k] += abs(c)
        for i in range(n):
            if (A >> i) & 1:
                true[i] += c
                part[i][k] += c
    scale = max(abs(x) for x in true)

    trunc = {}
    for q in range(1, n):
        err = 0.0
        for i in range(n):
            approx = sum(part[i][k] for k in range(1, q + 1))
            err = max(err, abs(true[i] - approx))
        trunc[q] = err / scale
    return true, trunc, mass, scale


def main():
    L = []
    add = L.append
    add("Order decomposition of the Shapley values from fully enumerated oracles")
    add("error normalised by max_i |phi_i|")
    add("")
    add("phi_i = sum_{A containing i} m(A)/|A|, split by |A| -- exact, no sampling")
    add("trunc(q) = max_i |phi_i - phi_i^{<=q}| / max_i |phi_i|")
    add("  = the error a PERFECT q-th-order model would still carry")
    add("")

    for SCEN in ("r_sat", "r_lim1"):
        add("=" * 78)
        add("scenario %s" % SCEN)
        add("=" * 78)
        tabs = {}
        for n, v in ((8, load_simple("scen_v.csv", SCEN)),
                     (12, load_simple("co_v.csv", SCEN)),
                     (16, load_n16(SCEN))):
            if len(v) != (1 << n):
                add("n=%d incomplete (%d of %d) -- skipped" % (n, len(v), 1 << n))
                continue
            tabs[n] = analyse(v, n, n)

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
        add("Irreducible floor per order (the number the paper needs):")
        for n in sorted(tabs):
            true, trunc, mass, scale = tabs[n]
            vals = ", ".join("q<=%d: %.1f%%" % (q, 100 * trunc[q])
                             for q in sorted(trunc)[:8])
            add("  n=%-3d %s" % (n, vals))
        add("")
        add("Spectral mass sum_{|A|=k} |m(A)|/|A| by order, normalised to order k=1:")
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
