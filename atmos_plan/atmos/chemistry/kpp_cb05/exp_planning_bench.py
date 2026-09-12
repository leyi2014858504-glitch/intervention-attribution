"""Planning-benchmark rerun: MCTS vs greedy vs random, on validated numerics.

WHY
  The original "MCTS loses to greedy" verdict was produced on the pykpp stack
  whose numbers are voided (solver accuracy failed gates G1/G2).  The rerun is
  therefore required for an honest record -- and it is FREE: the full 2^8
  value table from the Möbius pilot serves as an exact, noiseless oracle.

  Rerunning only the ORIGINAL (unconstrained) instance would be unfair to
  MCTS: we measured that its optimum is the corner "all VOC off", dominated by
  four main effects, so greedy hill-climbing is provably near-optimal there and
  no amount of lookahead can pay.  This script therefore adds a CONSTRAINED
  regime where the interaction spectrum is decision-relevant.

REGIMES
  R0 unconstrained : feasible = all 256 coalitions.  Control; the expected
                     result is a tie (reproducing the old negative verdict on
                     clean numerics).
  R1 treaty+cap    : feasible masks satisfy  >=2 NOx windows OFF (treaty)  and
                     <=3 VOC windows OFF (cost cap).  Because switching NOx off
                     RAISES O3 in this NOx-saturated chemistry, every NOx pick
                     looks locally bad, and the pair NOx@W3+W4 carries a +7.5
                     ppb synergy: a greedy, myopic sequence can be steered into
                     a bad NOx pair.  This is the classic lookahead trap.

ARMS (all score best-feasible-found; budget k = distinct coalition
      evaluations, all served from the exact table)
  greedy  : steepest-descent additions with feasibility propagation; may stop
            only on feasible states; keeps the best feasible visited.
  random  : uniform sample over feasible masks.
  mcts    : UCT over the addition sequence (actions = extendable single
            additions), random feasible playouts, terminal option at feasible
            nodes, anytime best-feasible reporting.

PRE-REGISTERED VERDICT (fixed before running)
  In R1, "MCTS adds planning value" iff its median regret at some k <= 48 is
  <= 0.5 x min(greedy, random) at that k.  Otherwise the planning problem on
  this instance is trivial for MCTS too, and the old negative verdict is
  UPGRADED from broken-numerics artifact to verified fact.
  R0 is a control, not a claim: any arm may win; ties are the expected outcome.

Usage: python exp_planning_bench.py [n_seeds]
"""
import csv
import math
import os
import random as _random
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8
MAX1 = 1 << NPLAY
SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
KS = [8, 16, 24, 32, 48, 64]
NOX_MIN = 2      # treaty: >= 2 NOx windows off
VOC_MAX = 3      # cap   : <= 3 VOC windows off


def pc(m):
    return bin(m).count("1")


def load():
    v = {}
    with open(os.path.join(HERE, "mobius_v_dummy.csv")) as f:
        for row in csv.DictReader(f):
            nm = str(row.get("name", ""))
            if nm.startswith("mob_"):
                v[int(nm.split("_")[1])] = float(row["v_ppb"])
    return v


V = load()


def feasible_R0(m):
    return True


def feasible_R1(m):
    return pc(m & 0x0F) >= NOX_MIN and pc(m >> 4) <= VOC_MAX


def extendable_R0(m):
    return pc(m) < NPLAY


def extendable_R1(m):
    """can a feasible superset still be reached from partial m?"""
    if pc(m >> 4) > VOC_MAX:
        return False
    # NOx: up to 4 slots remain for the treaty
    return pc(m & 0x0F) + (NPLAY - pc(m)) >= NOX_MIN and pc(m) < NPLAY


# ------------------------------------------------------------------ arms
class BudgetExhausted(Exception):
    pass


class Eval:
    """budget meter: every UNSEEN mask costs one evaluation once, after the
    budget is spent any new mask raises -- arms must stop cleanly and report
    their best-so-far."""

    def __init__(self, limit):
        self.seen = set()
        self.limit = limit

    def __call__(self, m):
        if m not in self.seen:
            if len(self.seen) >= self.limit:
                raise BudgetExhausted()
            self.seen.add(m)
        return V[m]

    @property
    def n(self):
        return len(self.seen)

    def full(self):
        return self.n >= self.limit


def arm_greedy(feas, ext, ev):
    """steepest-descent additions; returns the best FEASIBLE mask visited."""
    cur, best = 0, (0 if feas(0) else None)
    try:
        while True:
            cands = [i for i in range(NPLAY)
                     if not (cur >> i) & 1 and ext(cur | (1 << i))]
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


def arm_greedy_repair(feas, ext, ev):
    """myopic policy that must obey the constraints: normally steepest descent;
    while still infeasible it takes the least-harmful addition even if v rises
    (a repair step).  This is the honest 'simple strategy' baseline -- plain
    descent can never enter a region that requires a short-term sacrifice."""
    cur, best = 0, (0 if feas(0) else None)
    try:
        while True:
            cands = [i for i in range(NPLAY)
                     if not (cur >> i) & 1 and ext(cur | (1 << i))]
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


def arm_random(feas, ext, ev, rng):
    best = 0 if feas(0) else None
    pool = [m for m in range(1, MAX1) if feas(m)]
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


def arm_mcts(feas, ext, ev, rng, worst, best_v, pool):
    """UCT on the addition sequence; rollouts are random feasible playouts.

    EVERY coalition touched (tree node or playout step) is billed through ev;
    on BudgetExhausted the search stops and reports best-feasible-visited
    (anytime behaviour, same as the other arms).
    """
    import math
    T = {0: [0.0, 1]}
    best = [0 if feas(0) else None]

    def record(m):
        if feas(m) and (best[0] is None or V[m] < V[best[0]]):
            best[0] = m

    def children(m):
        return [m | (1 << i) for i in range(NPLAY)
                if not (m >> i) & 1 and ext(m | (1 << i))]

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
        s, n = T[c]
        return -s / n + C_ROLL * math.sqrt(math.log(max(tot, 2)) / n)

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


# ------------------------------------------------------------------ driver
def run():
    out = []
    for rname, feas, ext in (("R0", feasible_R0, extendable_R0),
                             ("R1", feasible_R1, extendable_R1)):
        pool = [m for m in range(MAX1) if feas(m) and m != 0] or [m for m in range(MAX1) if feas(m)]
        oracle = min(pool, key=lambda m: V[m])
        BEST, WORST = V[oracle], max(V[m] for m in range(MAX1))
        print("%s: %d feasible coalitions, oracle = %s (%+0.3f ppb)"
              % (rname, len(pool), format(oracle, "08b"), V[oracle]))
        for k in KS:
            for arm in ("greedy", "greedy+", "random", "mcts"):
                regs = []
                for seed in range(SEEDS):
                    rng = _random.Random(seed * 977 + k)
                    ev = Eval(k)
                    ev(0)
                    if arm == "greedy":
                        rec = arm_greedy(feas, ext, ev)
                    elif arm == "greedy+":
                        rec = arm_greedy_repair(feas, ext, ev)
                    elif arm == "random":
                        rec = arm_random(feas, ext, ev, rng)
                    else:
                        rec = arm_mcts(feas, ext, ev, rng, WORST, BEST, pool)
                    if rec is None:
                        regs.append(float("inf"))
                    else:
                        regs.append(V[rec] - V[oracle])
                med = st.median(regs)
                out.append(dict(r=rname, k=k, arm=arm, med=med,
                                fail=sum(1 for x in regs if x == float("inf")) / len(regs)))
                print("  %-6s k=%3d %-7s median_regret=%+7.3f  infeasible=%.2f"
                      % ("", k, arm, med, out[-1]["fail"]), flush=True)
    # verdict
    L = ["Planning rerun on the exact 2^8 oracle (validated numerics)", ""]
    for r in ("R0", "R1"):
        L.append("%s:" % r)
        for k in KS:
            cells = []
            for arm in ("greedy", "greedy+", "random", "mcts"):
                x = next(o for o in out if o["r"] == r and o["k"] == k and o["arm"] == arm)
                cells.append("%s=%+0.3f" % (arm, x["med"]))
            L.append("  k=%3d  %s" % (k, "  ".join(cells)))
    mcts_claim = False
    for k in [kk for kk in KS if kk <= 48]:
        g = min(o["med"] for o in out if o["r"] == "R1" and o["k"] == k
                and o["arm"] in ("greedy", "greedy+", "random"))
        m = next(o["med"] for o in out if o["r"] == "R1" and o["k"] == k and o["arm"] == "mcts")
        if math.isfinite(g) and math.isfinite(m) and g > 0 and m <= 0.5 * g:
            mcts_claim = True
    L.append("")
    L.append("PRE-REGISTERED VERDICT (R1, k<=48): %s"
             % ("MCTS adds planning value" if mcts_claim else
                "no arm separated -> planning trivial, negative verdict VERIFIED on clean numerics"))
    txt = "\n".join(L)
    print("\n" + txt)
    with open(os.path.join(HERE, "planning_bench.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    run()
