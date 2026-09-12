"""Surrogate-based optimization of a time-resolved intervention recipe.

OBJECTIVE (policy framing, minimisation)
  v(S) = mean O3 over the last 24 h minus the baseline, so the goal is to find
  the switch set S that REDUCES ozone the most:  minimise v.
  Ground truth is the full 2^8 enumeration (256 coalitions), therefore regret
  and feasible-set accuracy are exact and this experiment costs no chemistry.

WHY THIS IS NOT THE SAME QUESTION AS BEFORE
  The earlier STOP verdict (Moebius spectrum not k-sparse; truncating at order
  <=3 still leaves 6.7% attribution error) is about EXACT ATTRIBUTION.
  Optimisation of a recipe needs only the argmin and the feasible level set --
  a much weaker requirement.  So the negative result does NOT carry over, and
  it has to be measured separately.

ARMS (all scored at the same number of simulated coalitions k)
  random      : k uniformly random coalitions
  dopt        : greedy D-optimal design for the order-q model
  size_strat  : round robin over coalition size |S| (classical size strata)
  struct      : screening then structure-aware.  An order-1 D-optimal design of
                9 runs screens the factors; the pairs among the top-3 factors
                (standard effect-sparsity assumption, fixed in advance, not
                read off the answer) form a reduced model, and the remaining
                budget is a greedy D-optimal design for that reduced model.
  coordinate  : greedy coordinate ascent on the TRUE v -- the "simple strategy"
                baseline.  Repeatedly evaluate all single flips from the
                current point (cached), move to the best, until k is spent.
  full256     : fit on every coalition (reference upper bound)

METRICS (exact, against the 256-point ground truth)
  regret      v(recommended recipe) - v(true best)                        [ppb]
  hit         recommended recipe truly reaches the target level
  prec        fraction of the surrogate's feasible set that truly is feasible
  topk_prec   fraction of the surrogate's top-|F_true| ranking that is feasible
  TARGET: T = v_best + 0.25 (v_worst - v_best) -- 75% of the achievable gain.

Usage: python exp_surrogate_opt.py [n_seeds]
"""
import csv
import itertools
import math
import os
import statistics as st
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8
MAX1 = 1 << NPLAY
FULL = MAX1 - 1
SEEDS = int(sys.argv[1]) if len(sys.argv) > 1 else 40
KS = [8, 16, 32, 48, 64, 96, 128]
ORDERS = [1, 2, 3]
RIDGE = 1e-6


def pc(m):
    return bin(m).count("1")


def load_v():
    v = {}
    with open(os.path.join(HERE, "mobius_v_dummy.csv")) as f:
        for row in csv.DictReader(f):
            nm = str(row.get("name", ""))
            if nm.startswith("mob_"):
                v[int(nm.split("_")[1])] = float(row["v_ppb"])
    base = v[0]
    return {m: v[m] - base for m in v}


# --------------------------------------------------------------------------
# model: v(S) ~ c + sum_{1<=|A|<=q} m(A) * 1[A subset S]
# --------------------------------------------------------------------------
def cols_of(order):
    out = [0]
    for k in range(1, order + 1):
        for combo in itertools.combinations(range(NPLAY), k):
            m = 0
            for i in combo:
                m |= 1 << i
            out.append(m)
    return out


def design_matrix(masks, cols):
    X = np.zeros((len(masks), len(cols)))
    for r, S in enumerate(masks):
        for c, A in enumerate(cols):
            if A == 0:
                X[r, c] = 1.0
            elif (S & A) == A:
                X[r, c] = 1.0
    return X


def fit(vv, masks, order):
    cols = cols_of(order)
    X = design_matrix(masks, cols)
    y = np.array([vv[S] for S in masks])
    A = X.T @ X + RIDGE * np.trace(X.T @ X) / len(cols) * np.eye(len(cols))
    beta = np.linalg.solve(A, X.T @ y)
    return cols, beta


def predict(vv, mask, fitted):
    cols, beta = fitted
    x = np.array([1.0 if (A == 0 or (mask & A) == A) else 0.0 for A in cols])
    return float(x @ beta)


# --------------------------------------------------------------------------
# designs
# --------------------------------------------------------------------------
def greedy_dopt(order, k, banned=(), weights=None):
    """greedy D-optimal rows for the order-`order` model (leverage rule)."""
    cols = cols_of(order)
    X = design_matrix(range(MAX1), cols)
    if weights is not None:
        X = X * weights
    chosen = []
    if weights is None:
        w = np.ones(len(cols))
    else:
        w = weights
    for _ in range(k):
        cur = np.array(chosen, dtype=int)
        if len(cur):
            Xc = X[cur]
            A = Xc.T @ Xc + RIDGE * np.eye(len(cols))
        else:
            A = RIDGE * np.eye(len(cols))
        Ainv = np.linalg.inv(A)
        lev = np.einsum("ij,jk,ik->i", X, Ainv, X)
        best, bv = None, -1.0
        for r in range(MAX1):
            if r in chosen or r in banned:
                continue
            if lev[r] > bv:
                bv, best = lev[r], r
        if best is None:
            break
        chosen.append(best)
    return chosen


def design_random(k, rng):
    pool = list(range(1, MAX1))
    rng.shuffle(pool)
    return pool[:k]


def design_size_strat(k, rng):
    by_layer = {}
    for m in range(1, MAX1):
        by_layer.setdefault(pc(m), []).append(m)
    for v in by_layer.values():
        rng.shuffle(v)
    out, layer = [], 1
    while len(out) < k:
        for s in sorted(by_layer):
            if len(out) >= k:
                break
            if by_layer[s]:
                out.append(by_layer[s].pop())
        layer += 1
        if layer > k + 2:
            break
    return out


ACTIVE_T = 3          # effect-sparsity assumption, declared in advance


def design_struct(k, rng, vv):
    """screen (order-1 D-optimal, 9 runs) then a structure-aware reduced design."""
    if k <= 9:
        return greedy_dopt(1, k)
    screen = greedy_dopt(1, 9)
    fitted = fit(vv, screen, 1)
    cols, beta = fitted
    mains = sorted(((abs(beta[c]), (cols[c]).bit_length() - 1)
                    for c in range(1, len(cols))), reverse=True)
    act = [i for _, i in mains[:ACTIVE_T]]
    keep = [0] + [1 << i for i in act]
    for a, b in itertools.combinations(act, 2):
        keep.append((1 << a) | (1 << b))
    # greedy D-optimal for the reduced (structured) model, rest of the budget
    cols_r = keep
    X = design_matrix(range(MAX1), cols_r)
    chosen = list(screen)
    for _ in range(k - len(screen)):
        Xc = X[np.array(chosen, dtype=int)]
        A = Xc.T @ Xc + RIDGE * np.eye(len(cols_r))
        Ainv = np.linalg.inv(A)
        lev = np.einsum("ij,jk,ik->i", X, Ainv, X)
        best, bv = None, -1.0
        for r in range(MAX1):
            if r in chosen:
                continue
            if lev[r] > bv:
                bv, best = lev[r], r
        if best is None:
            break
        chosen.append(best)
    return chosen


def coordinate_ascent(vv, k):
    """simple strategy: steepest coordinate ascent from the empty set.

    Returns the best point actually visited within the budget (that is the
    recipe you would report), not the last point reached.
    """
    seen = {0}
    cur = 0
    while len(seen) < k:
        improved, bv = None, vv[cur]
        for i in range(NPLAY):
            if cur & (1 << i):
                continue
            nxt = cur | (1 << i)
            if nxt in seen or len(seen) >= k:
                continue
            seen.add(nxt)
            if vv[nxt] < bv:
                bv, improved = vv[nxt], nxt
        if improved is None:
            break
        cur = improved
    return min(seen, key=lambda m: vv[m]), len(seen)


# --------------------------------------------------------------------------
# scoring
# --------------------------------------------------------------------------
def score(vv, fitted, T, feas, k_feas):
    preds = np.array([predict(vv, m, fitted) for m in range(MAX1)])
    rec = int(np.argmin(preds))
    pf = set(int(m) for m in np.where(preds <= T)[0])
    prec = (len(pf & feas) / len(pf)) if pf else 1.0
    topk = set(int(m) for m in np.argsort(preds, kind="stable")[:max(1, k_feas)])
    topk_prec = len(topk & feas) / len(topk)
    return dict(regret=vv[rec] - vbest, hit=bool(vv[rec] <= T),
                prec=prec, topk_prec=topk_prec, recipe=rec)


vv = load_v()
vbest = min(vv.values())
vworst = max(vv.values())
T = vbest + 0.25 * (vworst - vbest)
feas = set(m for m in vv if vv[m] <= T)
k_feas = len(feas)
best_mask = min(vv, key=lambda m: vv[m])

print("ground truth over the 2^8 enumeration")
print("  v(emptyset)=0   best %+.3f ppb   worst %+.3f ppb" % (vbest, vworst))
print("  target T = %+.3f ppb ; %d of 256 coalitions are feasible" % (T, k_feas))
print("  best recipe: %s" % "+".join(
    "%s@W%d" % ("NOx" if i < 4 else "VOC", (i % 4) + 1)
    for i in range(NPLAY) if best_mask & (1 << i)))
print("")

rows = []
cache = {}


def get_design(arm, q, k, seed):
    """dopt/struct are deterministic -- build once and reuse across seeds."""
    key = (arm, q, k)
    if arm in ("random", "size_strat"):
        rng = np.random.default_rng(seed)
        return (design_random(k, rng) if arm == "random"
                else design_size_strat(k, rng))
    if key not in cache:
        cache[key] = (greedy_dopt(q, k) if arm == "dopt"
                      else design_struct(k, None, vv))
    return cache[key]


for arm in ("random", "dopt", "size_strat", "struct"):
    for q in ORDERS:
        for k in KS:
            res = []
            for seed in range(SEEDS):
                d = get_design(arm, q, k, seed)
                res.append(score(vv, fit(vv, d, q), T, feas, k_feas))
            rows.append(dict(arm=arm, q=q, k=k,
                             regret=st.median(r["regret"] for r in res),
                             hit=sum(r["hit"] for r in res) / len(res),
                             prec=st.median(r["prec"] for r in res),
                             topk=st.median(r["topk_prec"] for r in res)))
            print("%-11s q=%d k=%3d  regret=%+7.3f  hit=%.2f  prec=%.2f  topk=%.2f"
                  % (arm, q, k, rows[-1]["regret"], rows[-1]["hit"],
                     rows[-1]["prec"], rows[-1]["topk"]), flush=True)

for k in KS:
    regs = []
    for seed in range(SEEDS):
        rec, used = coordinate_ascent(vv, k)
        regs.append(vv[rec] - vbest)
    rows.append(dict(arm="coordinate", q=-1, k=k, regret=st.median(regs),
                     hit=float(np.mean([r <= T - vbest for r in regs])),
                     prec=float("nan"), topk=float("nan")))
    print("%-11s          k=%3d  regret=%+7.3f  hit=%.2f"
          % ("coordinate", k, rows[-1]["regret"], rows[-1]["hit"]), flush=True)

cols_all, beta_all = fit(vv, list(range(MAX1)), NPLAY)
ref = score(vv, (cols_all, beta_all), T, feas, k_feas)
print("\nreference (fit on all 256): regret=%+.3f hit=%s prec=%.2f topk=%.2f"
      % (ref["regret"], ref["hit"], ref["prec"], ref["topk_prec"]))

L = []
L.append("Surrogate-based optimisation of a time-resolved intervention recipe")
L.append("  objective: minimise v (reduce O3); ground truth = full 2^8 enumeration")
L.append("  best %+.3f ppb, worst %+.3f ppb, target T = %+.3f ppb, %d/256 feasible"
         % (vbest, vworst, T, k_feas))
L.append("  reference (all 256 points): regret %+.3f, prec %.2f"
         % (ref["regret"], ref["prec"]))
L.append("")
L.append("  arm            q    k   regret     hit   prec  topk_prec")
for r in rows:
    L.append("  %-11s %3s %4d  %+8.3f  %5.2f  %5.2f  %6.2f"
             % (r["arm"], r["q"] if r["q"] > 0 else "-", r["k"], r["regret"],
                r["hit"], r["prec"], r["topk"]))
txt = "\n".join(L)
print("\n" + txt)
with open(os.path.join(HERE, "surrogate_opt.txt"), "w") as f:
    f.write(txt + "\n")
with open(os.path.join(HERE, "surrogate_opt.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["arm", "q", "k", "regret", "hit", "prec", "topk"])
    w.writeheader()
    w.writerows(rows)
