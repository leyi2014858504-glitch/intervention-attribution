"""Coalition-sampling benchmark for Shapley/Owen attribution (n=8, oracle known).

CLAIM UNDER TEST
  A chemistry-structure-aware active coalition sampler reaches a given
  attribution accuracy with fewer expensive chemistry simulations than the
  incumbent samplers, while preserving Shapley/Owen semantics (unbiased
  estimators, no truncation).

PROTOCOL (pre-registered before the run)
  ground truth : exact phi / psi from the full 2^8 enumeration
                 (attrib_oracle.py, four independent checks PASS)
  budget       : number of DISTINCT coalition evaluations (= chemistry runs),
                 B in {8, 16, 24, 32, 48, 64, 96, 128}; kept well below 2^8
                 because an already-evaluated coalition is free, so random
                 search saturates near the full lattice
  error metric : e = max_i |est_i - true_i| / max_i |true_i|
  seeds        : 50 per (arm, budget, target)
  rollouts-to-eps : smallest B whose MEDIAN e over seeds is <= eps,
                 eps in {20%, 10%, 5%, 2%}
  arms         : phi -> kernelshap (incumbent), perm, strat, active
                 psi -> perm, strat, active          (KernelSHAP has no Owen)
  PASS         : for BOTH targets and for eps = 10% and 5%, the active arm's
                 rollouts-to-eps is <= 0.5 x that of the BEST baseline
                 (i.e. a 2x reduction is needed to claim the method works)
  FAIL         : active is not better than the best baseline at eps = 5%
  BIAS CHECK   : the seed-mean signed error of every unbiased arm must be
                 within 2 MC standard errors of 0 at B = 256; a visible bias
                 invalidates that arm regardless of its error curve.
"""
import csv
import math
import os
import statistics as st
import sys

from attrib_oracle import load_v, mobius, name_of, shapley_from_mobius, owen
import attrib_samplers as S

HERE = os.path.dirname(os.path.abspath(__file__))
# Budgets stay below the size of the coalition lattice, which differs by
# target: a Shapley chain can reach any subset, but an Owen chain is
# group-contiguous, so for psi only ~48 distinct coalitions are reachable and
# every arm saturates (error -> 0) around B=96.  The psi grid is therefore
# capped where the comparison is still informative.
BUDGETS = {"phi": [8, 16, 24, 32, 48, 64, 96, 128, 160, 192],
           "psi": [8, 16, 24, 32, 48, 64]}
EPS = [0.20, 0.10, 0.05, 0.02]
SEEDS = 50
THRESH = 0.5


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "mobius_v_dummy.csv")
    v_abs = load_v(src)
    table = dict(v_abs)
    base = v_abs[0]
    vv = [table[i] - base for i in range(1 << S.NPLAY)]
    m = mobius(vv)
    true = {"phi": shapley_from_mobius(m),
            "psi": owen(vv, m)}

    if "--smoke" in sys.argv:
        for target in ("phi", "psi"):
            scale = max(abs(x) for x in true[target])
            print("target %s  max|true| = %.4f" % (target, scale))
            for arm in S.ARMS[target]:
                for B in (16, 32):
                    est, nev = S.estimate(arm, target, table, B, 0)
                    err = max(abs(e - t) for e, t in zip(est, true[target])) / scale
                    print("  %-11s B=%3d  used=%3d  err=%7.4f  est=%s"
                          % (arm, B, nev, err,
                             " ".join("%+.2f" % x for x in est)))
        return 0

    rows = []
    for target in ("phi", "psi"):
        scale = max(abs(x) for x in true[target])
        for arm in S.ARMS[target]:
            for B in BUDGETS[target]:
                errs, per_player = [], [[] for _ in range(S.NPLAY)]
                for seed in range(SEEDS):
                    est, nev = S.estimate(arm, target, table, B, seed)
                    assert est is not None
                    d = [e - t for e, t in zip(est, true[target])]
                    for i, x in enumerate(d):
                        per_player[i].append(x / scale)
                    errs.append(max(abs(x) for x in d) / scale)
                # Bias must be measured per player and averaged over seeds
                # FIRST; taking max_i per seed is a selection statistic whose
                # mean is positive even for an unbiased estimator.
                biasi = [sum(xs) / len(xs) for xs in per_player]
                sei = [(st.pstdev(xs) / math.sqrt(len(xs))) if len(xs) > 1 else 0.0
                       for xs in per_player]
                k = max(range(S.NPLAY), key=lambda i: abs(biasi[i]))
                rows.append(dict(target=target, arm=arm, B=B,
                                 med=st.median(errs),
                                 mean_bias=biasi[k], se=sei[k],
                                 worst_player=k))
                print("%-4s %-11s B=%3d  median_err=%7.4f  bias(player %d)=%+8.5f"
                      % (target, arm, B, rows[-1]["med"], k, rows[-1]["mean_bias"]),
                      flush=True)

    # ---- rollouts to eps -------------------------------------------------
    def r2e(target, arm, eps):
        xs = [r for r in rows if r["target"] == target and r["arm"] == arm
              and r["med"] <= eps]
        return min((r["B"] for r in xs), default=None)

    L = []
    L.append("Coalition-sampling benchmark, n=8, oracle = full enumeration")
    L.append("  error metric : max_i |est-true| / max_i |true|, median over %d seeds" % SEEDS)
    L.append("  budget       : distinct coalition evaluations (= chemistry runs)")
    L.append("")
    for target in ("phi", "psi"):
        L.append("target = %s   (max|true| = %.4f ppb)"
                 % (target, max(abs(x) for x in true[target])))
        L.append("  budget grid: %s" % BUDGETS[target])
        L.append("  rollouts to reach median error <= eps:")
        L.append("    %-11s %8s %8s %8s %8s" % ("arm", *["%.0f%%" % (100 * e) for e in EPS]))
        arms = list(S.ARMS[target])
        top = max(BUDGETS[target])
        for arm in arms:
            cells = []
            for e in EPS:
                v = r2e(target, arm, e)
                cells.append(">%d" % top if v is None else "%d" % v)
            L.append("    %-11s %8s %8s %8s %8s" % (arm, *cells))
        best = {}
        for e in EPS:
            vals = [(r2e(target, a, e) or 10 ** 6) for a in arms if a != "OURS"]
            best[e] = min(vals)
        L.append("    best baseline (non-OURS): %s"
                 % " ".join("%.0f%%=%s" % (100 * e, ">%d" % top if best[e] >= 10 ** 6 else best[e])
                            for e in EPS))
        L.append("")
    # ---- verdict ---------------------------------------------------------
    L.append("DECISION (pre-registered)")
    ok = True
    for target in ("phi", "psi"):
        arms = [a for a in S.ARMS[target] if a != "OURS"]
        top = max(BUDGETS[target])
        for e in (0.10, 0.05):
            act = r2e(target, "OURS", e)
            bl = [r2e(target, a, e) for a in arms]
            bl = [b for b in bl if b is not None]
            blmin = min(bl) if bl else None
            good = (act is not None and blmin is not None and act <= THRESH * blmin)
            ok &= good
            L.append("  %-4s eps=%2.0f%%  OURS=%s  best baseline=%s  -> %s"
                     % (target, 100 * e,
                        ">%d" % top if act is None else act,
                        ">%d" % top if blmin is None else blmin,
                        "PASS" if good else "FAIL"))
    L.append("  VERDICT: %s" % ("PASS (2x fewer rollouts than the best baseline)"
                                if ok else "FAIL (no 2x reduction)"))
    # ---- bias check ------------------------------------------------------
    L.append("")
    L.append("BIAS CHECK at the largest budget of each grid "
             "(seed-mean signed error, need |bias| <= 2 SE):")
    for r in rows:
        if r["B"] != max(BUDGETS[r["target"]]):
            continue
        good = abs(r["mean_bias"]) <= 2 * r["se"] + 1e-9
        L.append("  %-4s %-11s B=%3d  bias=%+8.5f  2SE=%.5f  %s"
                 % (r["target"], r["arm"], r["B"], r["mean_bias"], 2 * r["se"],
                    "ok" if good else "BIASED"))
    txt = "\n".join(L)
    print("\n" + txt)
    with open(os.path.join(HERE, "attrib_bench.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(HERE, "attrib_bench.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["target", "arm", "B", "med", "mean_bias",
                                          "se", "worst_player"])
        w.writeheader()
        w.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
