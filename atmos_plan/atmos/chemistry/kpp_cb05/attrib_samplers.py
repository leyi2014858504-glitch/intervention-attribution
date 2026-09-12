"""Coalition samplers for Shapley / Owen attribution of the n=8 emission game.

CURRENCY
  One unit of budget = one evaluation of v at a coalition that has not been
  evaluated before (= one expensive chemistry run).  v({}) = 0 is known and
  free.  Everything is scored against the exact oracle (attrib_oracle.py).

WHY NOT A TRUNCATED SPECTRUM MODEL
  phi_i = sum over A containing i of m(A)/|A|, and the order decomposition of
  this game is 79% (k=1), 27% (k=2), 16% (k=3), 12% (k=4) of max|phi|, so the
  irreducible truncation floor is 17.7% at order<=2 and 6.7% at order<=3.  No
  low-order model can reach 5% accuracy, and truncation also breaks the
  efficiency axiom.  The method must therefore be UNBIASED and win on SAMPLING
  DESIGN alone.

ARMS
  kernelshap : the incumbent.  Uniform random coalitions, Shapley-kernel
               weighted constrained least squares.  (No Owen variant exists.)
  strat      : unbiased marginal estimator.  phi_i is an expectation of the
               marginal contribution Delta_i over contexts drawn uniformly
               within each layer, so estimates are unbiased by construction.
               Allocation is uniform, i.e. not adapted.
  active     : OUR METHOD.  Same unbiased estimator, but the budget is
               allocated across (player, context-layer) cells by a two-stage
               Neyman rule: a pilot stage measures the marginal-contribution
               dispersion sigma_{i,c} of every cell, then the remaining budget
               is split proportionally to sigma.  This is the
               "chemistry-structure-aware" part: the chemistry decides which
               (player, time-context) cells carry the variance.
  perm       : random-order (permutation) sampling, the classic unbiased
               Shapley estimator; for Owen the same thing with the order drawn
               group-wise.

TARGETS
  phi : Shapley values (8).  Context for player i = subset S of the others,
        layer s = |S|; layer mass is uniform (1/n each), which reproduces the
        exact Shapley weights.
  psi : Owen values with the coalition structure given by the four time
        windows.  Context for player i in window k = (T over the other
        windows, Q over the other member of k), layer = (|T|, |Q|), uniform
        over the 4x2 layer combinations, which reproduces the exact Owen
        weights.
"""
import itertools
import math
import random

NPLAY = 8
WINDOW = 4
FULL = (1 << NPLAY) - 1


def pc(m):
    return bin(m).count("1")


class Budget:
    """counts distinct coalition evaluations against a fixed oracle table.

    Note the currency: a coalition already seen is FREE.  Near the size of the
    lattice that makes random search saturate (coupon collector), so every arm
    loop carries a stall guard and the benchmark keeps B well below 2^n.
    """

    def __init__(self, table):
        self.t = table
        self.seen = {0}            # v({}) = 0 is known analytically
        self.n = 0

    def v(self, mask):
        if mask not in self.seen:
            self.seen.add(mask)
            self.n += 1
        return self.t[mask]


# --------------------------------------------------------------------------
# run-1 design: independent coalition-marginal estimation (NOT in ARMS).
# Kept for provenance -- results_attrib_bench_run1.txt was produced with it,
# and it is the experiment that showed why this foundation loses in the
# "distinct coalition evaluations" currency.  Do not use for new results.
# --------------------------------------------------------------------------
def shapley_contexts(i, n=NPLAY):
    """(layer, subset) drawn uniformly over the n layers, uniform within."""
    others = [j for j in range(n) if j != i]
    out = []
    for s in range(n):                      # layers 0..n-1
        for combo in itertools.combinations(others, s):
            m = 0
            for j in combo:
                m |= 1 << j
            out.append((s, m))
    return out


def owen_contexts(i, groups):
    """(layer tuple, (T, Q)) with T over other groups, Q over i's teammates."""
    k = next(g for g, B in enumerate(groups) if i in B)
    others = [h for h in range(len(groups)) if h != k]
    mates = [j for j in groups[k] if j != i]
    out = []
    for t in range(len(others)):
        for tcomb in itertools.combinations(others, t):
            T = 0
            for h in tcomb:
                for j in groups[h]:
                    T |= 1 << j
            for q in range(len(mates) + 1):
                for qcomb in itertools.combinations(mates, q):
                    Q = 0
                    for j in qcomb:
                        Q |= 1 << j
                    out.append(((t, q), T, Q))
    return out


CONTEXTS = {
    "phi": {i: shapley_contexts(i) for i in range(NPLAY)},
    "psi": {i: owen_contexts(i, [[w, 4 + w] for w in range(WINDOW)])
            for i in range(NPLAY)},
}
CELLS = {"phi": sorted({c[0] for c in CONTEXTS["phi"][0]}),
         "psi": sorted({c[0] for c in CONTEXTS["psi"][0]})}


def delta(bud, target, i, ctx):
    """marginal contribution of i in the given context (cost 1-2 evaluations)."""
    if target == "phi":
        _, S = ctx
        return bud.v(S | (1 << i)) - bud.v(S)
    _, T, Q = ctx
    base = T | Q
    return bud.v(base | (1 << i)) - bud.v(base)


# --------------------------------------------------------------------------
# arms
# --------------------------------------------------------------------------
def arm_strat(bud, target, budget, seed, active=False, pilot_frac=0.35):
    """unbiased marginal estimator; active=True turns on the two-stage Neyman
    allocation across (player, layer) cells."""
    rng = random.Random(seed)
    by_cell = {}                      # (i, layer) -> [ctx, ...]
    for i in range(NPLAY):
        for ctx in CONTEXTS[target][i]:
            by_cell.setdefault((i, ctx[0]), []).append(ctx)
    cells = sorted(by_cell)
    acc = {(i, c): [] for (i, c) in cells}

    def draw(i, c):
        ctx = rng.choice(by_cell[(i, c)])
        acc[(i, c)].append(delta(bud, target, i, ctx))

    stall, MAXSTALL = 0, 600

    def guard(before):
        nonlocal stall
        stall = stall + 1 if bud.n == before else 0
        return stall > MAXSTALL

    if not active:
        # uniform over cells, round-robin until the budget is spent
        while bud.n < budget:
            for (i, c) in cells:
                if bud.n >= budget:
                    break
                before = bud.n
                draw(i, c)
                if guard(before):
                    return _combine(acc, target)
    else:
        # stage 1: pilot, one sample per cell (uniform coverage)
        for (i, c) in cells:
            if bud.n >= budget:
                break
            before = bud.n
            draw(i, c)
            if guard(before):
                return _combine(acc, target)
        # stage 2: Neyman allocation from the pilot dispersions (recomputed
        # every REFRESH draws to keep the sweep cheap)
        refresh, k, wcell = 16, 0, None
        while bud.n < budget:
            if k % refresh == 0:
                wcell = {}
                for key in cells:
                    xs = acc[key]
                    wcell[key] = max(_std(xs), 1e-12) if len(xs) > 1 else 1e-12
            k += 1
            tot = sum(wcell.values())
            r, acc_w, pick = rng.random() * tot, 0.0, cells[-1]
            for key in cells:
                acc_w += wcell[key]
                if r <= acc_w:
                    pick = key
                    break
            before = bud.n
            draw(*pick)
            if guard(before):
                break

    return _combine(acc, target)


def _combine(acc, target):
    """unbiased combination: every context layer carries equal mass.

    Below full coverage (budget < number of cells) some layers have no sample;
    those fall back to the player's pooled mean over its covered layers, which
    is a bias.  exp_attrib.py reports the coverage point so the reader knows
    from which budget the estimator is exactly unbiased.
    """
    est = []
    for i in range(NPLAY):
        vals = [sum(xs) / len(xs) for c in CELLS[target] if (xs := acc[(i, c)])]
        est.append(sum(vals) / len(vals) if vals else 0.0)
    return est


def n_cells(target):
    return NPLAY * len(CELLS[target])


def _std(xs):
    m = sum(xs) / len(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def _chain(bud, target, order):
    """evaluate one visit order; return the marginal contribution of each player.

    The chain of nested sets is what makes permutation sampling cheap in the
    "distinct coalition evaluations" currency: n evaluations yield a marginal
    for ALL n players (compare with the independent (i,S) estimator, which
    spends ~2 evaluations per single marginal).
    """
    d = [0.0] * NPLAY
    S = 0
    for i in order:
        d[i] = bud.v(S | (1 << i)) - bud.v(S)
        S |= 1 << i
    return d


def _rand_order(rng, target, hub=None, hub_pos=None):
    """uniform order; for psi the visit order is group-contiguous.

    If hub/hub_pos are given, return a uniform order conditioned on the hub
    (psi: the hub's whole window group) sitting at that position.
    """
    if target == "phi":
        players = list(range(NPLAY))
        rest = [i for i in players if i != hub] if hub is not None else players
        rng.shuffle(rest)
        if hub is None:
            return rest
        return rest[:hub_pos] + [hub] + rest[hub_pos:]
    groups = [[w, 4 + w] for w in range(WINDOW)]
    gid = list(range(WINDOW))
    if hub is not None:
        hg = next(g for g, B in enumerate(groups) if hub in B)
        rest = [g for g in gid if g != hg]
        rng.shuffle(rest)
        rest = rest[:hub_pos] + [hg] + rest[hub_pos:]
    else:
        rest = gid
        rng.shuffle(rest)
    out = []
    for g in rest:
        inside = list(groups[g])
        rng.shuffle(inside)
        out += inside
    return out


def arm_perm(bud, target, budget, seed):
    """random-order sampling (baseline): unbiased Shapley / Owen estimator."""
    rng = random.Random(seed)
    acc = [[] for _ in range(NPLAY)]
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        for i, d in enumerate(_chain(bud, target, _rand_order(rng, target))):
            acc[i].append(d)
        stall = stall + 1 if bud.n == before else 0
    return [sum(a) / len(a) for a in acc]


def arm_perm_anti(bud, target, budget, seed):
    """antithetic orders: each draw also walks the reversed order.

    Both orders are uniformly distributed, so the pair average stays unbiased;
    the two chains visit complementary contexts, which is where the variance
    reduction comes from.  No chemistry prior and no pilot cost, so it works at
    every budget.
    """
    rng = random.Random(seed)
    acc = [[] for _ in range(NPLAY)]
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        order = _rand_order(rng, target)
        for i, d in enumerate(_chain(bud, target, order)):
            acc[i].append(d)
        for i, d in enumerate(_chain(bud, target, list(reversed(order)))):
            acc[i].append(d)
        stall = stall + 1 if bud.n == before else 0
    return [sum(a) / len(a) for a in acc]


def arm_perm_strat(bud, target, budget, seed, anti=True):
    """OUR METHOD: hub-stratified, two-stage permutation sampling (+antithetic).

    Stage 1 draws a few uniform orders and identifies the HUB -- the player
    whose marginal contribution has the largest dispersion.  (In this chemistry
    it is the nocturnal-NOx switch, the player the interaction spectrum is
    organised around; that is the structure-awareness.)

    Stage 2 stratifies the order space by the POSITION of the hub, samples
    uniformly WITHIN a stratum, and allocates the remaining budget across
    strata by a Neyman rule proportional to the within-stratum dispersion.
    With anti=True each order is also walked backwards (a standard variance
    reduction that costs no extra chemistry once the two chains share
    coalitions).

    Unbiasedness: the hub's position is uniform under the reference (uniform
    order) distribution, so the stratification weights are exactly 1/n and the
    per-stratum means are unbiased for the conditional means; combining them
    with 1/n is therefore unbiased, and adaptive allocation does not break this
    because the combination weights stay fixed at the true ones.  Each order
    (forward or reversed) is filed under ITS OWN hub position.
    """
    rng = random.Random(seed)
    nstrat = NPLAY if target == "phi" else WINDOW
    acc = {(i, s): [] for i in range(NPLAY) for s in range(nstrat)}
    stall = 0

    def hubpos(order, hub):
        """stratum of this order = position of the hub (psi: of its group)."""
        k = order.index(hub)
        return k if target == "phi" else k // 2

    def emit(order, hub):
        """walk the order (and its reverse) and file each under its stratum."""
        for o in ((order, order[::-1]) if anti else (order,)):
            d = _chain(bud, target, o)
            s = hubpos(o, hub)
            for i in range(NPLAY):
                acc[(i, s)].append(d[i])

    # stage 1: uniform pilot, enough to pick the hub.  The pilot orders are
    # kept and re-filed once the hub is known, so nothing is wasted.
    npilot = 3 if target == "phi" else 2
    pilot = []
    for _ in range(npilot):
        order = _rand_order(rng, target)
        pilot.append((order, _chain(bud, target, order)))
    hub = max(((_std([d[i] for _, d in pilot]) if len(pilot) > 1 else 0.0), i)
              for i in range(NPLAY))[1]
    for order, d in pilot:
        s = hubpos(order, hub)
        for i in range(NPLAY):
            acc[(i, s)].append(d[i])

    # stage 2: cover every stratum once, then Neyman-allocate
    for s in range(nstrat):
        if bud.n >= budget:
            break
        before = bud.n
        emit(_rand_order(rng, target, hub, s), hub)
        if (stall := stall + 1 if bud.n == before else 0) > 3000:
            return _combine_strat(acc, nstrat)
    refresh, k, w = 4, 0, None
    while bud.n < budget:
        if k % refresh == 0:
            w = {}
            for s in range(nstrat):
                xs = [x for i in range(NPLAY) for x in acc[(i, s)]]
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
        emit(_rand_order(rng, target, hub, pick), hub)
        if (stall := stall + 1 if bud.n == before else 0) > 3000:
            break
    return _combine_strat(acc, nstrat)


def _combine_strat(acc, nstrat):
    """stratum means combined with the true (uniform) stratification weights."""
    est = []
    for i in range(NPLAY):
        means = [sum(acc[(i, s)]) / len(acc[(i, s)])
                 for s in range(nstrat) if acc[(i, s)]]
        est.append(sum(means) / len(means) if means else 0.0)
    return est


def arm_kernelshap(bud, target, budget, seed):
    """incumbent: uniform random coalitions + Shapley-kernel constrained WLS."""
    if target != "phi":
        return None
    rng = random.Random(seed)
    X, y, w = [], [], []
    n = NPLAY
    stall = 0
    while bud.n < budget and stall < 3000:
        before = bud.n
        mask = rng.randrange(1, FULL)          # exclude {} and N
        s = pc(mask)
        X.append([1 if mask & (1 << i) else 0 for i in range(n)])
        y.append(bud.v(mask))
        w.append(_kernel(n, s))
        stall = stall + 1 if bud.n == before else 0
    bud.v(FULL)                                # needed for the constraint
    buds = bud.v(0)
    return _kernel_wls(X, y, w, n, buds, bud.v(FULL))


def _kernel(n, s):
    if s == 0 or s == n:
        return 1e6
    return (n - 1) / (math.comb(n, s) * s * (n - s))


def _kernel_wls(X, y, w, n, v0, vN):
    """KernelSHAP: v(S) ~ c + sum_{i in S} phi_i with sum phi_i = vN - v0."""
    # eliminate phi_{n-1} using the efficiency constraint
    A, b, W = [], [], []
    for row, yy, ww in zip(X, y, w):
        a = [row[i] - row[n - 1] for i in range(n - 1)]
        A.append(a)
        b.append(yy - v0 - row[n - 1] * (vN - v0))
        W.append(ww)
    # weighted normal equations with a tiny ridge for stability
    m = n - 1
    M = [[sum(W[k] * A[k][i] * A[k][j] for k in range(len(A))) + (1e-9 if i == j else 0)
          for j in range(m)] for i in range(m)]
    rhs = [sum(W[k] * A[k][i] * b[k] for k in range(len(A))) for i in range(m)]
    sol = _solve(M, rhs)
    phi = list(sol) + [(vN - v0) - sum(sol)]
    return phi


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


ARMS = {
    # the run-1 coalition-marginal arms (strat / active) are not repeated here:
    # they lose structurally in this currency and are recorded in
    # results_attrib_bench_run1.txt.
    "phi": {"kernelshap": arm_kernelshap,
            "perm": arm_perm,
            "perm_anti": arm_perm_anti,
            "perm_strat": lambda b, t, k, s: arm_perm_strat(b, t, k, s, anti=False),
            "OURS": arm_perm_strat},
    "psi": {"perm": arm_perm,
            "perm_anti": arm_perm_anti,
            "perm_strat": lambda b, t, k, s: arm_perm_strat(b, t, k, s, anti=False),
            "OURS": arm_perm_strat},
}


def estimate(arm, target, table, budget, seed):
    bud = Budget(table)
    out = ARMS[target][arm](bud, target, budget, seed)
    return out, bud.n
