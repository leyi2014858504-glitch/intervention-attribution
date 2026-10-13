#!/usr/bin/env python3
"""CO mass share by interaction order:  f_CO(k) = 1 - W_8(k)/W_12(k).

Premise (must hold bit-exactly): the n=8 game is the sub-lattice of the n=12
game obtained by leaving all four CO windows ON.  For a CO-free coalition A,
m_12(A) = m_8(A) exactly (Mobius inversion is local: every B subset of A is
also CO-free, and v_12(B) = v_8(B) there).  Hence

    RAW_W_12(k) = sum_{|A|=k} |m_12(A)|/k        (all 12 players)
    RAW_W_8(k)  = sum_{|A|=k, CO-free} |m_8(A)|/k
    W_12(k) - W_8(k) = mass at order k carried by CO-touching dividends

and f_CO(k) = 1 - W_8(k)/W_12(k) is the CO-touching share of order-k
Shapley-weighted mass (RAW ppb masses; no k=1 renormalisation, which would
corrupt the ratio).

Context number: the COUNT share of CO-touching terms at order k is
1 - C(8,k)/C(12,k); if CO were a generic-strength sector, f_CO(k) should
track that count share.  How far below it f_CO(k) sits measures CO weakness.

Zero chemistry runs -- released oracle tables only.
Writes results_co_share.txt
"""
import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCEN = ("r_sat", "r_lim1")


def mobius_fast(vals, n):
    f = list(vals)
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                f[mask] -= f[mask ^ bit]
    return f


def mass_by_order(vals, n):
    """raw Shapley-weighted mass W(k) in ppb, k=1..n"""
    m = mobius_fast(vals, n)
    w = [0.0] * (n + 1)
    for mask in range(1, 1 << n):
        k = bin(mask).count("1")
        w[k] += abs(m[mask]) / k
    return w


def read_simple(path, col, keep):
    """name '<scen>_m<mask>' -> {scen: {mask: v}}"""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            name = row.get("name", "")
            if "_m" not in name:
                continue
            s, _, ms = name.rpartition("_m")
            if s not in keep:
                continue
            v = row.get(col, "")
            if v in ("", "NaN"):
                continue
            try:
                out.setdefault(s, {})[int(ms)] = float(v)
            except ValueError:
                continue
    return out


def read_n16(scen, col):
    v = {}
    with open(os.path.join(HERE, "n16_v_%s_m8.csv" % scen)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


def embed8(m8):
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def load_pair(metric):
    """-> {scen: (v8, v12)} for the given metric column."""
    out = {}
    if metric == "d3_mean":
        v8all = read_simple(os.path.join(HERE, "scen_v.csv"), "v_ppb", SCEN)
        v12all = read_simple(os.path.join(HERE, "co_v.csv"), "v_ppb", SCEN)
        for s in SCEN:
            out[s] = (v8all[s], v12all[s])
    else:  # MDA8
        v12all = {}
        for K in range(1, 7):
            p = os.path.join(HERE, "co_v_m8_%d.csv" % K)
            part = read_simple(p, "v_mda8", SCEN)
            for s, d in part.items():
                v12all.setdefault(s, {}).update(d)
        for s in SCEN:
            v16 = read_n16(s, "v_mda8")
            v8 = {m: v16[embed8(m)] for m in range(256)}
            out[s] = (v8, v12all[s])
    return out


def main():
    L = []
    add = L.append
    add("CO-touching share of spectral mass by order:  f_CO(k) = 1 - W_8(k)/W_12(k)")
    add("W = RAW Shapley-weighted mass sum_{|A|=k}|m(A)|/|A| in ppb (no renormalisation)")
    add("n=8 = n=12 restricted to CO-on (sub-lattice); premise re-verified below.")
    for metric in ("MDA8", "d3_mean"):
        pairs = load_pair(metric)
        add("")
        add("=" * 78)
        add("metric = %s" % metric)
        add("=" * 78)
        for s in SCEN:
            v8, v12 = pairs[s]
            # ---- premise check: CO-free masks must reproduce v8 bit-exactly
            dmax, n_cf = 0.0, 0
            for m12, val in v12.items():
                if (m12 >> 8) != 0:
                    continue
                n_cf += 1
                dmax = max(dmax, abs(val - v8[m12 & 0xFF]))
            add("")
            add("scenario %s  [%s]  CO-free coalitions checked: %d"
                % (s, metric, n_cf))
            add("  sub-lattice gate: max|v12(CO-on) - v8| = %.3e ppb  %s"
                % (dmax, "OK" if dmax < 1e-9 else "!! PREMISE FAILS"))
            vals8 = [v8[m] for m in range(256)]
            vals12 = [v12[m] for m in range(4096)]
            w8 = mass_by_order(vals8, 8)
            w12 = mass_by_order(vals12, 12)
            add("  %5s %12s %12s %10s %12s" %
                ("k", "W_8(k)", "W_12(k)", "f_CO(k)", "count share"))
            tot_num = tot_den = 0.0
            for k in range(1, 13):
                if k <= 8:
                    f = 1.0 - w8[k] / w12[k] if w12[k] > 0 else float("nan")
                    cnt = 1.0 - math.comb(8, k) / math.comb(12, k)
                    add("  %5d %12.4f %12.4f %9.2f%% %11.2f%%"
                        % (k, w8[k], w12[k], 100 * f, 100 * cnt))
                    tot_num += w12[k] - w8[k]
                else:
                    add("  %5d %12s %12.4f %9s %11s"
                        % (k, "-", w12[k], "100.00%*", "100.00%"))
                tot_den += w12[k]
            add("  (* for k>8 every order-k divisor must touch CO: f_CO = 1 by counting)")
            add("  TOTAL mass share touching CO: %.2f%%  (sum_k (W12-W8) / sum_k W12)"
                % (100 * tot_num / tot_den))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(ROOT, "results_co_share.txt"), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote results_co_share.txt")


if __name__ == "__main__":
    main()