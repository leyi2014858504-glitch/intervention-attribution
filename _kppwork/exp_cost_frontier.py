"""Is enumeration the only way?  Cost-accuracy frontier from data we ALREADY
have (co_v.csv: 2 scenarios x 4096 exact coalitions).  Zero new chemistry runs.

Methods compared, all scored against the exact Shapley (4096 evaluations):

  exact          phi_i from the full Moebius spectrum
  trunc-q        phi_i using only |A| <= q  (runs = sum_{k<=q} C(n,k))
  two-stage      hierarchical allocation: 3-sector Shapley (2^3 = 8 runs) to
                 get sector totals, then within-sector Shapley on each 4-window
                 sub-game (3 x 2^4 = 48 runs) to split them  -> 56 runs
  owen           exact OWEN value with the coalition structure
                 {NOx}{VOC}{CO} -- needs the full table too, but it is the
                 axiomatically right object when the block structure is known
                 (its block sums equal the 3-sector Shapley exactly)

Reported: runs used, and max |est - true| / max |true| for individual phi_i,
plus the same for the three sector totals.
"""
import csv
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
N = 12
M = 1 << N
SCENS = ["r_sat", "r_lim1"]
BLOCKS = [[0, 1, 2, 3], [4, 5, 6, 7], [8, 9, 10, 11]]   # sectors
BID = {i: k for k, b in enumerate(BLOCKS) for i in b}


def pc(m):
    return bin(m).count("1")


def role(i):
    return "NOx" if i < 4 else ("VOC" if i < 8 else "CO")


def mobius(v):
    out = [0.0] * M
    for mask in range(M):
        s, sub = 0.0, mask
        while True:
            s += (1.0 if pc(mask ^ sub) % 2 == 0 else -1.0) * v[sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out[mask] = s
    return out


def shapley_from(m, n=N):
    return [sum(m[A] / pc(A) for A in range(1, 1 << n) if (A >> i) & 1)
            for i in range(n)]


def owen(v, blocks):
    """Exact Owen value for a known coalition structure."""
    m = len(blocks)
    out = [0.0] * N
    for k, B in enumerate(blocks):
        others = [j for j in range(m) if j != k]
        for i in B:
            rest = [x for x in B if x != i]
            tot = 0.0
            for csub in range(1 << len(others)):
                C = [others[t] for t in range(len(others)) if (csub >> t) & 1]
                pre = sum(1 << x for j in C for x in blocks[j])
                wC = (math.factorial(len(C)) * math.factorial(m - len(C) - 1)
                      / math.factorial(m))
                for dsub in range(1 << len(rest)):
                    D = [rest[t] for t in range(len(rest)) if (dsub >> t) & 1]
                    s = pre | sum(1 << x for x in D) | (1 << i)
                    wD = (math.factorial(len(D))
                          * math.factorial(len(B) - len(D) - 1)
                          / math.factorial(len(B)))
                    tot += wC * wD * (v[s] - v[s & ~(1 << i)])
            out[i] = tot
    return out


def two_stage(v, blocks):
    """8 + 48 evaluations: sector Shapley, then split each sector's total."""
    m = len(blocks)
    # --- top level: one binary per sector (whole sector on/off)
    def V(sel):
        mask = 0
        for k, B in enumerate(blocks):
            if (sel >> k) & 1:
                for i in B:
                    mask |= 1 << i
        return v[mask]
    sec = {}
    for k in range(m):
        others = [j for j in range(m) if j != k]
        tot = 0.0
        for csub in range(1 << len(others)):
            C = [others[t] for t in range(len(others)) if (csub >> t) & 1]
            w = (math.factorial(len(C)) * math.factorial(m - len(C) - 1)
                 / math.factorial(m))
            sel = sum(1 << j for j in C)
            tot += w * (V(sel | (1 << k)) - V(sel))
        sec[k] = tot
    # --- within sector: 4-window sub-game, other sectors held ON
    out = [0.0] * N
    for k, B in enumerate(blocks):
        sub = {}
        for amask in range(1 << 4):
            mask = 0
            for i in range(N):
                if i in B:
                    if not (amask >> B.index(i)) & 1:
                        mask |= 1 << i        # ON unless the sub-mask turns it off
                else:
                    mask |= 1 << i            # other sectors ON
            sub[amask] = v[mask]
        base = sub[0]
        w4 = [0.0] * 4
        for t in range(4):
            tot = 0.0
            for A in range(1, 16):
                if (A >> t) & 1:
                    tot += (sub[A] - base) / pc(A)
            w4[t] = tot
        tot4 = sum(w4)
        for t, i in enumerate(B):
            out[i] = sec[k] * (w4[t] / tot4) if abs(tot4) > 1e-12 else 0.0
    return out


def main():
    tab = {}
    with open(os.path.join(HERE, "co_v.csv")) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if s in SCENS and row["v_ppb"] not in ("", "NaN"):
                tab.setdefault(s, {})[int(ms)] = float(row["v_ppb"])

    L = []
    for s in SCENS:
        if len(tab.get(s, {})) != M:
            L.append("%s incomplete" % s)
            continue
        v = [tab[s][m] for m in range(M)]
        vr = [x - v[0] for x in v]
        mo = mobius(vr)
        true = shapley_from(mo)
        true_sec = [sum(true[i] for i in b) for b in BLOCKS]
        scale = max(abs(x) for x in true)

        L.append("")
        L.append("=" * 78)
        L.append("scenario %s   n=12, exact enumeration = %d runs, max|phi|=%.2f ppb"
                 % (s, M, scale))
        L.append("=" * 78)
        L.append("  %-14s %8s %10s %10s" %
                 ("method", "runs", "err(phi_i)", "err(sector)"))

        # truncation q = 0..5
        for q in range(0, 6):
            runs = sum(math.comb(N, k) for k in range(q + 1))
            est = [sum(mo[A] / pc(A) for A in range(1, M)
                       if (A >> i) & 1 and pc(A) <= q) for i in range(N)]
            est_sec = [sum(est[i] for i in b) for b in BLOCKS]
            L.append("  %-14s %8d %9.2f%% %9.2f%%"
                     % ("trunc q<=%d" % q, runs,
                        100 * max(abs(est[i] - true[i]) for i in range(N)) / scale,
                        100 * max(abs(est_sec[k] - true_sec[k])
                                  for k in range(3)) / scale))

        # two-stage hierarchical (8 + 48)
        ts = two_stage(vr, BLOCKS)
        ts_sec = [sum(ts[i] for i in b) for b in BLOCKS]
        L.append("  %-14s %8d %9.2f%% %9.2f%%"
                 % ("two-stage", 8 + 3 * 16,
                    100 * max(abs(ts[i] - true[i]) for i in range(N)) / scale,
                    100 * max(abs(ts_sec[k] - true_sec[k])
                              for k in range(3)) / scale))

        # exact Owen with the sector structure
        ow = owen(vr, BLOCKS)
        ow_sec = [sum(ow[i] for i in b) for b in BLOCKS]
        L.append("  %-14s %8d %9.2f%% %9.2f%%"
                 % ("owen(exact)", M,
                    100 * max(abs(ow[i] - true[i]) for i in range(N)) / scale,
                    100 * max(abs(ow_sec[k] - true_sec[k])
                              for k in range(3)) / scale))

        # how far is the two-stage from Owen (its real target)?
        L.append("")
        L.append("  two-stage vs OWEN (its true target): err(phi_i) = %.2f%%,"
                 " err(sector) = %.2f%%"
                 % (100 * max(abs(ts[i] - ow[i]) for i in range(N)) / scale,
                    100 * max(abs(ts_sec[k] - ow_sec[k])
                              for k in range(3)) / scale))
        L.append("  sector Shapley from 8 runs: NOx %+7.2f  VOC %+7.2f  CO %+7.2f"
                 % tuple(ow_sec))
        L.append("  exact  phi sector totals  : NOx %+7.2f  VOC %+7.2f  CO %+7.2f"
                 % tuple(true_sec))

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "cost_frontier.txt"), "w") as f:
        f.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
