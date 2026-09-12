"""Adaptivity arbiter: does stochastic planning have value in this chemistry?

Reads scen_v.csv (scenario x coalition mean-O3 metric from /root/kpp5).

U-SET (10 weather/parameter perturbations of ONE planning problem) decides
route A. PRE-REGISTERED rule, fixed before any scenario was run:

    gap(U) = mean_s v_s(committed optimum) - mean_s v_s(per-scenario optimum)
    gap <= 0.5 ppb  -> planning trivial even under uncertainty: route A CLOSED
    gap >  0.5 ppb  -> genuine stochastic value of information: route A LIVE

plus a trap census: from the committed recipe m_c, can constrained greedy
(repair rule) reach each scenario's own optimum by monotone/repair steps, or
does it need lookahead depth > 1?

R-SET (6 EKMA-regime walk) is route-B material: per-regime Shapley vector,
spectral levels R_2/R_3 vs the null baseline, optimum location, and the
sign pattern of "reduce-NOx-only" -- the regime dependence of attribution.
"""
import csv
import itertools
import os
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY, MAX1 = 8, 256
SETS = {}


def load():
    tab = {}
    with open(os.path.join(HERE, "scen_v.csv")) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_" not in nm:
                continue
            scen, ms = nm.rsplit("_m", 1)
            assert int(row["nfail"]) == 0 and int(row["rc"]) == 0, row
            tab.setdefault(scen, {})[int(ms)] = float(row["v_ppb"])
    with open(os.path.join(HERE, "scen_list.txt")) as f:
        for ln in f:
            s, grp = ln.split()
            SETS[s] = grp
    return tab


def pc(m):
    return bin(m).count("1")


def mask_name(m):
    return "+".join("%s@W%d" % ("NOx" if i < 4 else "VOC", (i % 4) + 1)
                    for i in range(NPLAY) if m & (1 << i)) or "none"


def mobius(v):
    out = [0.0] * MAX1
    for mask in range(MAX1):
        s, sub = 0.0, mask
        while True:
            s += (1 if (pc(mask ^ sub) % 2) == 0 else -1) * v[sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out[mask] = s
    return out


def shapley(m):
    return [sum(m[A] / pc(A) for A in range(1, MAX1) if A & (1 << i))
            for i in range(NPLAY)]


def spectral_levels(m):
    byk = {}
    for mask in range(1, MAX1):
        byk.setdefault(pc(mask), []).append(abs(m[mask]))
    mx1 = max(abs(m[1 << i]) for i in range(NPLAY))
    return {k: max(v) / mx1 for k, v in byk.items()}, mx1


def greedy_repair_reachable(v, start, target, feas):
    """constrained greedy from `start`: descend where possible, else take the
    least-harmful feasible-extension step; can it reach the target optimum?"""
    cur, seen = start, set()
    for _ in range(64):
        if cur == target:
            return True
        if cur in seen:
            return False
        seen.add(cur)
        cands = [cur | (1 << i) for i in range(NPLAY) if not (cur >> i) & 1]
        if not cands:
            return False
        nxt = min(cands, key=lambda c: v[c])
        if v[nxt] < v[cur] or not feas(cur):
            cur = nxt
        else:
            return False
    return False


def main():
    tab = load()
    us = [s for s, g in SETS.items() if g == "U"]
    rs = [s for s, g in SETS.items() if g == "R"]
    L = []
    L.append("Adaptivity arbiter (pre-registered: gap(U)<=0.5 ppb closes route A)")
    L.append("")

    # ---------- U set: value of adaptivity ------------------------------
    mean_v = [sum(tab[s][m] for s in us) / len(us) for m in range(MAX1)]
    m_c = min(range(1, MAX1), key=lambda m: mean_v[m])
    per_opt = {s: min(range(1, MAX1), key=lambda m: tab[s][m]) for s in us}
    committed = sum(tab[s][m_c] for s in us) / len(us)
    adaptive = sum(tab[s][per_opt[s]] for s in us) / len(us)
    gap = committed - adaptive
    L.append("U-set (%d scenarios x 256):" % len(us))
    L.append("  committed optimum   = %-24s  mean v = %.3f ppb"
             % (mask_name(m_c), committed))
    opts = {}
    for s in us:
        opts.setdefault(per_opt[s], []).append(s)
    L.append("  per-scenario optima : %d distinct recipes over %d scenarios"
             % (len(opts), len(us)))
    for m, ss in sorted(opts.items(), key=lambda kv: -len(kv[1])):
        L.append("     %-24s x%-3d  %s" % (mask_name(m), len(ss), ",".join(ss)))
    L.append("  ADAPTIVITY GAP(U)   = %.4f ppb   -> %s"
             % (gap, "ROUTE A CLOSED (trivial under uncertainty)" if gap <= 0.5
                else "ROUTE A LIVE (real stochastic value)"))
    # trap census: repair-greedy from m_c vs each scenario optimum
    feas = lambda m: True
    trapped = sum(1 for s in us
                  if not greedy_repair_reachable(tab[s], m_c, per_opt[s], feas))
    L.append("  greedy-repair from m_c reaches each scenario optimum: %d/%d"
             % (len(us) - trapped, len(us)))
    # regret of the best SINGLE weekday plan under each scenario (spread)
    sp = [tab[s][m_c] - tab[s][per_opt[s]] for s in us]
    L.append("  per-scenario regret at m_c: max %.3f  median %.3f ppb"
             % (max(sp), st.median(sp)))
    L.append("")

    # ---------- R set: regime dependence of attribution -----------------
    L.append("R-set: attribution structure across EKMA regimes")
    L.append("  regime      max|phi| R_2    R_3    argmin recipe            reduce-NOx-only?")
    for s in rs:
        base = tab[s][0]
        v = [tab[s][m] - base for m in range(MAX1)]
        msp = mobius(v)
        lv, mx1 = spectral_levels(msp)
        phi = shapley(msp)
        arg = min(range(MAX1), key=lambda m: tab[s][m])
        nox_only = tab[s][0x0F] < tab[s][0]        # all 4 NOx off vs baseline
        L.append("  %-9s %7.2f %5.2f %6.2f  %-24s   %s (NOx-all-off %+0.2f ppb)"
                 % (s, mx1, lv.get(2, 0), lv.get(3, 0), mask_name(arg),
                    "worsens" if nox_only else "improves", tab[s][0x0F] - base))
        # top-2 Shapley players per regime
        top2 = sorted(range(NPLAY), key=lambda i: -abs(phi[i]))[:2]
        L.append("            phi: %s"
                 % "  ".join("%s %+0.2f" % (mask_name(1 << i), phi[i]) for i in top2))
    L.append("")
    L.append("VERDICT: gap(U)=%.4f ppb %s" % (
        gap,
        "<= 0.5 -> route A is CLOSED even under weather uncertainty; "
        "paper = route B (framework + regime-dependent attribution, "
        "planning boundary proven)" if gap <= 0.5 else
        "> 0.5 -> stochastic planner question is LIVE; build the "
        "closed-loop instance next"))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "stoch_gap.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    sys.exit(main())
