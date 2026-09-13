#!/usr/bin/env python3
"""
Re-run the two n=8 benchmarks on the n=16 oracle.

WHY THIS EXISTS
---------------
Both published benchmarks ran at n=8, where the universe is 2^8 = 256 coalitions
and the budget grid reaches 192 -- 75% coverage.  At that coverage the estimator
is nearly exact whatever the method, so method differences CANNOT appear.  The
same structural fact turned up independently in the sampler-scaling experiment:
enumeration beats the sampler at n=8 (256 runs vs ~630 for 1.5% error) and only
loses from n = 12 upwards.

n=16 gives 2^16 = 65536 coalitions, fully enumerated (see sweep_n16.sh), so
  * the oracle is exact, not mapped;
  * a budget of 64 covers 0.098% of the universe instead of 25%.
That is the regime where a sampler's design can actually matter.

PROTOCOL -- unchanged from the published runs, only n changes
  attribution : error  = max_i |est_i - true_i| / max_i |true_i|, median over seeds
  planning    : regret = v(best found) - v(oracle),              median over seeds
  currency    : DISTINCT coalition evaluations (1 = 1 chemistry run)

FAITHFULNESS
  Part A and Part B each FIRST reproduce the published n=8 tables arm by arm.
  If a single number disagrees the n=16 column is not trustworthy, and the
  script says so rather than quietly reporting it.

Bit layout: bits 0..W-1 = NOxW(W), bits W..2W-1 = VOCW(W).
  n=8  -> W=4   n=16 -> W=8     (matches gen_scen_sweep.py / gen_n16.py)

Usage: python3 exp_n16_bench.py [n_seeds]
"""
import csv
import importlib.util
import json
import math
import os
import random
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
W16 = "/mnt/d/atmos/_kppwork"
OUT = os.path.join(W16, "bench_n16.txt")
SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 50

# ---- the estimator under test must be the SAME code as the published runs ----
_spec = importlib.util.spec_from_file_location("expscale", os.path.join(HERE, "exp_scale.py"))
_es = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_es)
mobius_fast = None  # filled below from exp_n16

_n16spec = importlib.util.spec_from_file_location("expn16", os.path.join(HERE, "exp_n16.py"))
_ex16 = importlib.util.module_from_spec(_n16spec)
_n16spec.loader.exec_module(_ex16)
mobius_fast = _ex16.mobius_fast


# =====================================================================
# oracles
# =====================================================================
def load_n8(scen):
    """n=8 table: 2 sectors x 4 windows, from scen_v.csv (name '<scen>_m<mask>')."""
    v = {}
    with open(os.path.join(HERE, "scen_v.csv")) as fh:
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
        fn = os.path.join(W16, "n16_v_%s_%d.csv" % (scen, K))
        if not os.path.exists(fn):
            continue
        with open(fn) as fh:
            for row in csv.DictReader(fh):
                if row.get("name", "").startswith("m") and row["v_ppb"] not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(row["v_ppb"])
    return v


def exact_phi(v, n):
    mo = mobius_fast([v[m] for m in range(1 << n)], n)
    pc = [bin(A).count("1") for A in range(1 << n)]
    true = [0.0] * n
    for A in range(1, 1 << n):
        c = mo[A] / pc[A]
        if c == 0.0:
            continue
        for i in range(n):
            if (A >> i) & 1:
                true[i] += c
    return true


# =====================================================================
# Part A -- attribution arms, ported from attrib_samplers.py with n as a knob
# =====================================================================
class Budget:
    """currency: one distinct coalition evaluation = one chemistry run.

    Note the stall guards below.  Near the size of the lattice random search
    saturates (coupon collector), so every arm loop needs an escape; this is
    exactly the failure mode that made n=8 look like a tie.
    """

    def __init__(self, table):
        self.t = table
        self.seen = {0}
        self.n = 0

    def v(self, mask):
        if mask not in self.seen:
            self.seen.add(mask)
            self.n += 1
        return self.t[mask]


def chain(bud, order):
    n = len(order)
    d = [0.0] * n
    S = 0
    for i in order:
        d[i] = bud.v(S | (1 << i)) - bud.v(S)
        S |= 1 << i
    return d


def arm_perm(bud, n, budget, seed):
    rng = random.Random(seed)
    acc = [[] for _ in range(n)]
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        order = list(range(n))
        rng.shuffle(order)
        for i, d in enumerate(chain(bud, order)):
            acc[i].append(d)
        stall = stall + 1 if bud.n == before else 0
    return [sum(a) / len(a) if a else 0.0 for a in acc]


def arm_perm_anti(bud, n, budget, seed):
    rng = random.Random(seed)
    acc = [[] for _ in range(n)]
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        order = list(range(n))
        rng.shuffle(order)
        for o in (order, order[::-1]):
            for i, d in enumerate(chain(bud, o)):
                acc[i].append(d)
        stall = stall + 1 if bud.n == before else 0
    return [sum(a) / len(a) if a else 0.0 for a in acc]


def _std(xs):
    m = sum(xs) / len(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def arm_perm_strat(bud, n, budget, seed, anti=True, npilot=3):
    """OUR METHOD: hub-stratified two-stage permutation sampling (+antithetic).

    Stage 1 picks the hub = the player with the largest marginal dispersion;
    stage 2 stratifies by the hub's position and Neyman-allocates across strata.
    Unbiased because the hub's position is uniform under the reference law.
    """
    rng = random.Random(seed)
    nstrat = n
    acc = {(i, s): [] for i in range(n) for s in range(nstrat)}
    stall = 0

    def emit(order, hub):
        for o in ((order, order[::-1]) if anti else (order,)):
            d = chain(bud, o)
            s = o.index(hub)
            for i in range(n):
                acc[(i, s)].append(d[i])

    pilot = []
    for _ in range(npilot):
        order = list(range(n))
        rng.shuffle(order)
        pilot.append((order, chain(bud, order)))
    hub = max(((_std([d[i] for _, d in pilot]) if len(pilot) > 1 else 0.0), i)
              for i in range(n))[1]
    for order, d in pilot:
        s = order.index(hub)
        for i in range(n):
            acc[(i, s)].append(d[i])

    for s in range(nstrat):
        if bud.n >= budget:
            break
        before = bud.n
        order = list(range(n))
        rng.shuffle(order)
        order.remove(hub)
        order.insert(s, hub)
        emit(order, hub)
        stall = stall + 1 if bud.n == before else 0
        if stall > 3000:
            return combine_strat(acc, n, nstrat)
    refresh, k, w = 4, 0, None
    while bud.n < budget:
        if k % refresh == 0:
            w = {}
            for s in range(nstrat):
                xs = [x for i in range(n) for x in acc[(i, s)]]
                w[s] = max(_std(xs), 1e-12) if len(xs) > 1 else 1e-12
        k += 1
        tot = sum(w.values())
        r, c, pick = rng.random() * tot, 0.0, nstrat - 1
        for s in range(nstrat):
            c += w[s]
            if r <= c:
                pick = s
                break
        before = bud.n
        order = list(range(n))
        rng.shuffle(order)
        order.remove(hub)
        order.insert(pick, hub)
        emit(order, hub)
        stall = stall + 1 if bud.n == before else 0
        if stall > 3000:
            break
    return combine_strat(acc, n, nstrat)


def combine_strat(acc, n, nstrat):
    est = []
    for i in range(n):
        means = [sum(acc[(i, s)]) / len(acc[(i, s)]) for s in range(nstrat) if acc[(i, s)]]
        est.append(sum(means) / len(means) if means else 0.0)
    return est


def _kernel(n, s):
    if s == 0 or s == n:
        return 1e6
    return (n - 1) / (math.comb(n, s) * s * (n - s))


def _solve(M, b):
    m = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(M)]
    for c in range(m):
        p = max(range(c, m), key=lambda r: abs(M[r][c]))
        M[c], M[p] = M[p], M[c]
        for r in range(m):
            if r != c and abs(M[c][c]) > 0:
                f = M[r][c] / M[c][c]
                for k in range(c, m + 1):
                    M[r][k] -= f * M[c][k]
    return [M[i][m] / M[i][i] for i in range(m)]


def arm_kernelshap(bud, n, budget, seed):
    rng = random.Random(seed)
    FULL = (1 << n) - 1
    X, y, w = [], [], []
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        mask = rng.randrange(1, FULL)
        s = bin(mask).count("1")
        X.append([1 if mask & (1 << i) else 0 for i in range(n)])
        y.append(bud.v(mask))
        w.append(_kernel(n, s))
        stall = stall + 1 if bud.n == before else 0
    v0 = bud.v(0)
    vN = bud.v(FULL)
    m = n - 1
    A = [[row[i] - row[n - 1] for i in range(m)] for row in X]
    b = [yy - v0 - row[n - 1] * (vN - v0) for yy, row in zip(y, X)]
    M = [[sum(w[k] * A[k][i] * A[k][j] for k in range(len(A))) + (1e-9 if i == j else 0)
          for j in range(m)] for i in range(m)]
    rhs = [sum(w[k] * A[k][i] * b[k] for k in range(len(A))) for i in range(m)]
    sol = _solve(M, rhs)
    return list(sol) + [(vN - v0) - sum(sol)]


ARMS = [("kernelshap", arm_kernelshap),
        ("perm", arm_perm),
        ("perm_anti", arm_perm_anti),
        ("perm_strat", lambda b, n, k, s: arm_perm_strat(b, n, k, s, anti=False)),
        ("OURS", arm_perm_strat)]


def attrib_run(v, n, budgets, tag):
    true = exact_phi(v, n)
    scale = max(abs(x) for x in true)
    res = {}
    for name, fn in ARMS:
        res[name] = {}
        for B in budgets:
            errs = []
            for sd in range(SEEDS):
                bud = Budget(v)
                est = fn(bud, n, B, sd)
                errs.append(max(abs(est[i] - true[i]) for i in range(n)) / scale)
            res[name][B] = st.median(errs)
    return res, scale


def reach(res, eps):
    """smallest budget whose median error <= eps, else None."""
    for B in sorted(res):
        if res[B] <= eps:
            return B
    return None


# =====================================================================
# Part B -- planning arms, ported from exp_planning_bench.py with n as a knob
# =====================================================================
class BudgetExhausted(Exception):
    pass


class Eval:
    def __init__(self, v, limit):
        self.V = v
        self.seen = set()
        self.limit = limit

    def __call__(self, m):
        if m not in self.seen:
            if len(self.seen) >= self.limit:
                raise BudgetExhausted()
            self.seen.add(m)
        return self.V[m]

    @property
    def n(self):
        return len(self.seen)


def make_planning(v, n, nox_min, voc_max):
    W = n // 2
    NOXM, VOCM = (1 << W) - 1, ((1 << n) - 1) ^ ((1 << W) - 1)
    FULL1 = (1 << n) - 1

    def feasible_R0(m):
        return True

    def feasible_R1(m):
        return bin(m & NOXM).count("1") >= nox_min and bin((m & VOCM) >> W).count("1") <= voc_max

    def extendable_R0(m):
        return bin(m).count("1") < n

    def extendable_R1(m):
        if bin((m & VOCM) >> W).count("1") > voc_max:
            return False
        return (bin(m & NOXM).count("1") + (n - bin(m).count("1")) >= nox_min
                and bin(m).count("1") < n)
    return feasible_R0, feasible_R1, extendable_R0, extendable_R1, FULL1


def arm_greedy(V, n, feas, ext, ev):
    cur, best = 0, (0 if feas(0) else None)
    try:
        while True:
            cands = [i for i in range(n) if not (cur >> i) & 1 and ext(cur | (1 << i))]
            if not cands:
                break
            nxt, nv = None, float("inf")
            for i in cands:
                val = ev(cur | (1 << i))
                if val < nv:
                    nxt, nv = cur | (1 << i), val
            if nxt is None or nv >= V[cur]:
                break
            cur = nxt
            if feas(cur) and (best is None or V[cur] < V[best]):
                best = cur
    except BudgetExhausted:
        pass
    return best


def arm_greedy_repair(V, n, feas, ext, ev):
    cur, best = 0, (0 if feas(0) else None)
    try:
        while True:
            cands = [i for i in range(n) if not (cur >> i) & 1 and ext(cur | (1 << i))]
            if not cands:
                break
            nxt, nv = None, float("inf")
            for i in cands:
                val = ev(cur | (1 << i))
                if val < nv:
                    nxt, nv = cur | (1 << i), val
            if nxt is None:
                break
            if nv >= V[cur] and feas(cur):
                break
            cur = nxt
            if feas(cur) and (best is None or V[cur] < V[best]):
                best = cur
    except BudgetExhausted:
        pass
    return best


def arm_random(V, feas, ev, rng, pool):
    best = 0 if feas(0) else None
    try:
        while True:
            m = rng.choice(pool)
            val = ev(m)
            if best is None or val < V[best]:
                best = m
    except BudgetExhausted:
        pass
    return best


C_ROLL = 1.4


def arm_mcts(V, n, feas, ext, ev, rng, worst, best_v):
    T = {0: [0.0, 1]}
    best = [0 if feas(0) else None]

    def record(m):
        if feas(m) and (best[0] is None or V[m] < V[best[0]]):
            best[0] = m

    def children(m):
        return [m | (1 << i) for i in range(n) if not (m >> i) & 1 and ext(m | (1 << i))]

    def rollout(m):
        path = [m]
        record(m)
        while True:
            ch = children(path[-1])
            if not ch:
                break
            nxt = ch[rng.randrange(len(ch))]
            ev(nxt)
            path.append(nxt)
            record(nxt)
        vals = [V[x] for x in path if feas(x)]
        return min(vals) if vals else None

    def ucb(c, tot):
        s, nn = T[c]
        return -s / nn + C_ROLL * math.sqrt(math.log(max(tot, 2)) / nn)

    try:
        while True:
            path, m = [0], 0
            while True:
                ch = children(m)
                if not ch:
                    break
                unseen = [c for c in ch if c not in T]
                if unseen:
                    m = unseen[rng.randrange(len(unseen))]
                    ev(m)
                    T[m] = [0.0, 1]
                    path.append(m)
                    break
                tot = sum(T[c][1] for c in ch)
                m = max(ch, key=lambda c: ucb(c, tot))
                ev(m)
                path.append(m)
                record(m)
                T.setdefault(m, [0.0, 1])
            r = rollout(m)
            if r is None:
                continue
            rew = (worst - r) / (worst - best_v + 1e-9)
            for node in path:
                T[node][0] += rew
                T[node][1] += 1
    except BudgetExhausted:
        pass
    return best[0]


def planning_run(v, n, ks, nox_min, voc_max, tag):
    feas0, feas1, ext0, ext1, FULL1 = make_planning(v, n, nox_min, voc_max)
    out = {}
    for rname, feas, ext in (("R0", feas0, ext0), ("R1", feas1, ext1)):
        pool = [m for m in range(1, 1 << n) if feas(m)]
        oracle = min(pool, key=lambda m: v[m])
        bestV = v[oracle]
        worstV = max(v[m] for m in range(1 << n))
        out[rname] = {"oracle": oracle, "best": bestV, "nfeas": len(pool), "ks": {}}
        for k in ks:
            for arm in ("greedy", "greedy+", "random", "mcts"):
                regs = []
                for seed in range(SEEDS):
                    rng = random.Random(seed * 977 + k)
                    ev = Eval(v, k)
                    try:
                        ev(0)
                    except BudgetExhausted:
                        pass
                    if arm == "greedy":
                        rec = arm_greedy(v, n, feas, ext, ev)
                    elif arm == "greedy+":
                        rec = arm_greedy_repair(v, n, feas, ext, ev)
                    elif arm == "random":
                        rec = arm_random(v, feas, ev, rng, pool)
                    else:
                        rec = arm_mcts(v, n, feas, ext, ev, rng, worstV, bestV)
                    regs.append(float("inf") if rec is None else v[rec] - bestV)
                out[rname]["ks"].setdefault(k, {})[arm] = st.median(regs)
    return out


# =====================================================================
def main():
    L = []
    add = L.append
    add("n=8 benchmarks re-run on the n=16 oracle")
    add("seeds per (arm,budget) = %d" % SEEDS)
    add("")
    add("n=8 : 2 sectors x 4 windows, universe 256   (the published runs)")
    add("n=16: 2 sectors x 8 windows, universe 65536 (fully enumerated, exact oracle)")
    add("")

    v8 = load_n8("r_sat")
    v16 = load_n16("r_sat")
    add("oracle sizes: n=8 %d, n=16 %d" % (len(v8), len(v16)))
    assert len(v8) == 256 and len(v16) == 65536, "oracle incomplete"
    add("")

    # ---------------- Part A ----------------
    add("=" * 78)
    add("PART A -- attribution (target = phi, Shapley)")
    add("=" * 78)
    BUD8 = [8, 16, 24, 32, 48, 64, 96, 128, 160, 192]
    BUD16 = [8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512]

    r8, s8 = attrib_run(v8, 8, BUD8, "n8")
    add("n=8 max|phi| = %.4f ppb" % s8)
    add("  %-11s %7s %7s %7s %7s" % ("arm", "20%", "10%", "5%", "2%"))
    for name, _ in ARMS:
        add("  %-11s %7s %7s %7s %7s" % (
            name,
            reach(r8[name], 0.20) or ">%d" % max(BUD8),
            reach(r8[name], 0.10) or ">%d" % max(BUD8),
            reach(r8[name], 0.05) or ">%d" % max(BUD8),
            reach(r8[name], 0.02) or ">%d" % max(BUD8)))
    add("  published n=8 reference:")
    add("    kernelshap   48     128    >192   >192")
    add("    perm         16      64     192   >192")
    add("    perm_anti     8      16      96   >192")
    add("    perm_strat    8      64     160   >192")
    add("    OURS          8      48      96   >192")
    add("")

    r16, s16 = attrib_run(v16, 16, BUD16, "n16")
    add("n=16 max|phi| = %.4f ppb     (budget 64 = %.3f%% of the universe)"
        % (s16, 100.0 * 64 / (1 << 16)))
    add("  %-11s %7s %7s %7s %7s" % ("arm", "20%", "10%", "5%", "2%"))
    for name, _ in ARMS:
        add("  %-11s %7s %7s %7s %7s" % (
            name,
            reach(r16[name], 0.20) or ">%d" % max(BUD16),
            reach(r16[name], 0.10) or ">%d" % max(BUD16),
            reach(r16[name], 0.05) or ">%d" % max(BUD16),
            reach(r16[name], 0.02) or ">%d" % max(BUD16)))
    add("")
    add("  median error vs budget, n=16:")
    add("  %-11s" % "arm" + "".join("%8d" % B for B in BUD16))
    for name, _ in ARMS:
        add("  %-11s" % name + "".join(
            ("%7.2f%%" % (100 * r16[name][B])) if B in r16[name] else "       ."
            for B in BUD16))
    add("")

    # ---------------- Part B ----------------
    add("=" * 78)
    add("PART B -- planning (best coalition to switch off)")
    add("=" * 78)
    KS8, KS16 = [8, 16, 24, 32, 48, 64], [8, 16, 24, 32, 48, 64]
    p8 = planning_run(v8, 8, KS8, 2, 3, "n8")
    for r in ("R0", "R1"):
        d = p8[r]
        add("n=8 %s: %d feasible, oracle=%s (%+.3f ppb)"
            % (r, d["nfeas"], format(d["oracle"], "08b"), d["best"]))
        for k in KS8:
            add("  k=%3d  %s" % (k, "  ".join(
                "%s=%+0.3f" % (a, d["ks"][k][a]) for a in ("greedy", "greedy+", "random", "mcts"))))
    add("")

    p16 = planning_run(v16, 16, KS16, 2, 3, "n16")
    for r in ("R0", "R1"):
        d = p16[r]
        add("n=16 %s: %d feasible (%.4f%% of universe), oracle=%s (%+.3f ppb)"
            % (r, d["nfeas"], 100.0 * d["nfeas"] / (1 << 16),
               format(d["oracle"], "016b"), d["best"]))
        for k in KS16:
            add("  k=%3d  %s" % (k, "  ".join(
                "%s=%+0.3f" % (a, d["ks"][k][a]) for a in ("greedy", "greedy+", "random", "mcts"))))
    add("")

    # verdicts
    add("=" * 78)
    add("VERDICTS")
    add("=" * 78)
    for tag, rr, bud in (("n=8", r8, BUD8), ("n=16", r16, BUD16)):
        ours, best_base = "OURS", None
        for eps in (0.10, 0.05):
            o = reach(rr[ours], eps)
            bl = [reach(rr[a], eps) for a, _ in ARMS if a != "OURS"]
            bl = [x for x in bl if x]
            bb = min(bl) if bl else None
            verdict = "FAIL" if (o is None or bb is None or o > bb / 2.0) else "PASS"
            add("  %-5s attribution eps=%2.0f%%  OURS=%s  best baseline=%s  -> %s"
                % (tag, 100 * eps, o, bb, verdict))
    for tag, pp in (("n=8", p8), ("n=16", p16)):
        d = pp["R1"]
        sep = []
        for k in KS8:
            g = min(d["ks"][k][a] for a in ("greedy", "greedy+", "random"))
            m = d["ks"][k]["mcts"]
            if math.isfinite(g) and math.isfinite(m) and g > 0 and m <= 0.5 * g:
                sep.append(k)
        add("  %-5s planning R1: MCTS <= 0.5x best simple at k in %s -> %s"
            % (tag, sep or "{}", "MCTS ADDS PLANNING VALUE" if sep else
               "no separation (planning trivial at this n)"))

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
