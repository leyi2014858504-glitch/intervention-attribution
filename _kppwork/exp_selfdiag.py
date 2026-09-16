#!/usr/bin/env python3
"""Oracle-free error flags for truncated-spectrum attribution -- validated
against the EXACT oracles we already enumerated.

MOTIVATION
  Every headline number in the paper (truncation divergence 1033%, cost law
  1.6-1.8x) is measured against a full-enumeration oracle.  On a 3D model no
  oracle will ever exist, so the architecture transfers only if the ERROR can
  be estimated WITHOUT ground truth.  The 0-D tables give the true errors, so
  candidate flags can be scored against reality.  Zero new box runs.

CANDIDATE FLAGS (none may look at the truth)
  g_q  efficiency violation : LS-fit order-<=q Mobius model m_hat on sampled
       coalitions, then g_q = |v(N) - sum m_hat| / |v(N)|.  A truncated model
       must miss the grand-coalition total by the tail mass it dropped.
       (Permutation estimators cannot use this flag -- telescoping makes them
       exactly efficient; the flag belongs to MODEL fits.)
  h_q  holdout misfit       : RMS residual on held-out coalitions / sd(v_fit).
  sp_q bootstrap spread     : max|phi_hat - phi_hat_boot| / max|phi| (weak by
       expectation: heavy tails make spread a poor statistic, cf. n=24 note).

TARGETS (oracle side, scoring only)
  t_q  true truncation error of phi (the Table-6 quantity)
  e_q  realized error of the sampled fit, max_i|phi_hat_i - phi_i|/max|phi|

PRE-REGISTERED VERDICTS (fixed before looking)
  V1 tracking   : Spearman(g_q, t_q) over q >= +0.8 at n=12 AND n=16, both
                  scenarios (the flag must rise where the true error rises).
  V2 safety gate: "largest q with g_q <= 0.05" must have t_q <= 0.20; in
                  divergent cases must not certify q >= 3.
  V3 hierarchy  : |Spearman(flag, e_q)| of g_q or h_q exceeds that of sp_q.

FEASIBILITY CAPS: P_q = sum_{k<=q} C(n,k) <= B_fit/2 and <= 3000.
Writes selfdiag.txt.  Run: python3 exp_selfdiag.py [scen]
"""
import csv
import importlib.util
import itertools
import os
import sys
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SEEDS = 8

_spec = importlib.util.spec_from_file_location("expn16", os.path.join(HERE, "exp_n16.py"))
_ex16 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ex16)
mobius_fast = _ex16.mobius_fast


def load_table(fn, scen=None):
    raw = {}
    with open(os.path.join(HERE, fn)) as fh:
        for row in csv.DictReader(fh):
            nm = row["name"]
            if row["v_ppb"] in ("", "NaN"):
                continue
            if scen is not None:
                if "_m" not in nm:
                    continue
                s, ms = nm.rsplit("_m", 1)
                if s != scen:
                    continue
                m = int(ms)
            else:
                if not nm.startswith("m"):
                    continue
                m = int(nm[1:])
            raw[m] = float(row["v_ppb"])
    base = raw[0]
    return {m: v - base for m, v in raw.items()}


def true_decomposition(v, n):
    vals = [v[m] for m in range(1 << n)]
    mo = mobius_fast(vals, n)
    pc = [bin(A).count("1") for A in range(1 << n)]
    true = [0.0] * n
    part = [{q: 0.0 for q in range(n + 1)} for _ in range(n)]
    for A in range(1, 1 << n):
        k = pc[A]
        c = mo[A] / k
        for i in range(n):
            if (A >> i) & 1:
                true[i] += c
                part[i][k] += c
    scale = max(abs(x) for x in true) or 1.0
    t = {}
    for q in range(1, n):
        t[q] = max(abs(true[i] - sum(part[i][k] for k in range(1, q + 1)))
                   for i in range(n)) / scale
    return np.array(true), t, scale


def columns_upto(n, q):
    cols = [sum(1 << i for i in c)
            for k in range(1, q + 1) for c in itertools.combinations(range(n), k)]
    return cols, {A: j for j, A in enumerate(cols)}


def design(masks, col_of):
    D = np.zeros((len(masks), len(col_of)))
    for r, S in enumerate(masks):
        A = S
        while A:
            j = col_of.get(A)
            if j is not None:
                D[r, j] = 1.0
            A = (A - 1) & S
    return D


def phi_from_cols(mv, cols, n):
    phi = np.zeros(n)
    for A, c in zip(cols, mv):
        if c == 0.0:
            continue
        k = bin(A).count("1")
        for i in range(n):
            if (A >> i) & 1:
                phi[i] += c / k
    return phi


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        return float("nan")
    rx = x.argsort().argsort().astype(float)
    ry = y.argsort().argsort().astype(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def run_case(v, n, Bs, qcap, phi_true, t_q, scale, add):
    Nmask = (1 << n) - 1
    vN = v[Nmask]
    universe = np.arange(1, Nmask)
    seeds = {8: 8, 12: 6, 16: 4}.get(n, SEEDS)
    rows = []
    for B in Bs:
        for seed in range(seeds):
            rng = np.random.default_rng(1000 * seed + B + n)
            masks = rng.choice(universe, size=min(B, len(universe)),
                               replace=False).tolist()
            if Nmask not in masks:
                masks.append(Nmask)
            rng.shuffle(masks)
            cut = int(0.8 * len(masks))
            fit, hold = masks[:cut], masks[cut:]
            yf = np.array([v[S] for S in fit])
            for q in range(1, qcap + 1):
                P = sum(comb(n, k) for k in range(1, q + 1))
                if P > len(fit) / 2 or P > 3000:
                    break
                cols, col_of = columns_upto(n, q)
                m_hat, *_ = np.linalg.lstsq(design(fit, col_of), yf, rcond=None)
                phi_hat = phi_from_cols(m_hat, cols, n)
                e_q = float(np.max(np.abs(phi_hat - phi_true)) / scale)
                g_q = abs(vN - float(m_hat.sum())) / (abs(vN) or 1.0)
                res = np.array([v[S] for S in hold]) - design(hold, col_of) @ m_hat
                h_q = float(np.sqrt(np.mean(res ** 2)) / (float(np.std(yf)) or 1.0))
                if P <= 1000:   # bootstrap refit only while cheap
                    bs = rng.choice(np.array(fit), size=len(fit), replace=True).tolist()
                    mb, *_ = np.linalg.lstsq(design(bs, col_of),
                                             np.array([v[S] for S in bs]), rcond=None)
                    sp_q = float(np.max(np.abs(phi_from_cols(mb, cols, n) - phi_hat)) / scale)
                else:
                    sp_q = float("nan")
                rows.append((B, q, e_q, g_q, h_q, sp_q))
    agg = {}
    for B, q, e_q, g_q, h_q, sp_q in rows:
        agg.setdefault((B, q), []).append((e_q, g_q, h_q, sp_q))
    add("  %-5s %-3s %8s %8s %8s %8s %8s" %
        ("B", "q", "t_q*", "e_q", "g_q", "h_q", "sp_q"))
    for (B, q) in sorted(agg):
        med = np.median(np.array(agg[(B, q)]), axis=0)
        add("  %-5d %-3d %8.1f %8.1f %8.1f %8.3f %8.1f" %
            (B, q, 100 * t_q[q], 100 * med[0], 100 * med[1], med[2], 100 * med[3]))
    return agg


def main():
    L = []
    add = L.append
    add("Oracle-free flags vs true errors (all % except h_q which is RMS/sd;")
    add("t_q* = Table-6 oracle truncation error; e_q = realized fit error)")
    scen_arg = sys.argv[1] if len(sys.argv) > 1 else None
    scen_list = [scen_arg] if scen_arg else ["r_sat", "r_lim1"]
    v1 = []
    v3 = []
    for scen in scen_list:
        for n, fn, Bs, qcap in (
                (8, "scen_v.csv", [64, 128, 250], 3),
                (12, "co_v.csv", [512, 1024, 2048], 4),
                (16, "n16_v_%s.csv" % scen, [1024, 2048, 4096, 8192], 4)):
            path = os.path.join(HERE, fn)
            if not os.path.exists(path):
                continue
            # n16_v_<scen>.csv: one file per scenario, names m<05d>
            # scen_v.csv / co_v.csv: mixed scenarios, names <scen>_m<ddd>
            v = load_table(fn, None if fn.startswith("n16_v_") else scen)
            if len(v) != 1 << n:
                add("\n%s n=%d incomplete (%d of %d) -- skipped" % (scen, n, len(v), 1 << n))
                continue
            add("\n" + "=" * 66)
            add("%s   n=%d" % (scen, n))
            add("=" * 66)
            phi_true, t_q, scale = true_decomposition(v, n)
            agg = run_case(v, n, Bs, qcap, phi_true, t_q, scale, add)
            Bbig = max(B for (B, q) in agg)
            qs = sorted(q for (B, q) in agg if B == Bbig)
            if len(qs) >= 3:
                med = np.array([np.median(np.array(agg[(Bbig, q)]), axis=0)
                                for q in qs])          # (nq, 4)
                rho_g = spearman(med[:, 1], [t_q[q] for q in qs])
                rho_h = spearman(med[:, 2], [t_q[q] for q in qs])
                v1.append((scen, n, rho_g, rho_h, qs))
                allr = {}
                for j, name in ((0, "e"), (1, "g"), (2, "h"), (3, "sp")):
                    flat = [(B, q) for (B, q) in agg]
                    x = np.array([np.median(np.array(agg[bq])[:, j]) for bq in flat])
                    e = np.array([np.median(np.array(agg[bq])[:, 0]) for bq in flat])
                    allr[name] = abs(spearman(x, e))
                v3.append((scen, n, allr))
    add("\n" + "=" * 66)
    add("V1  Spearman(flag_q, t_q) over q at max feasible B")
    add("=" * 66)
    for scen, n, rho_g, rho_h, qs in v1:
        add("  %-7s n=%-3d  rho(g)=%+.2f  rho(h)=%+.2f   q=%s" % (scen, n, rho_g, rho_h, qs))
    need = [r for s, n, r, rh, q in v1 if n in (12, 16)]
    add("  V1 (rho(g) >= +0.8 at n=12,16): %s" %
        ("PASS" if need and all(r >= 0.8 for r in need) else "FAIL/see rows"))
    add("\nV3  |Spearman(flag, e_q)| hierarchy at max B (want g or h > sp)")
    for scen, n, allr in v3:
        add("  %-7s n=%-3d  e:%.2f g:%.2f h:%.2f sp:%.2f" %
            (scen, n, allr["e"], allr["g"], allr["h"], allr["sp"]))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "selfdiag.txt"), "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
