#!/usr/bin/env python3
"""
Measured (not mapped) scaling of the permutation sampler in n.

The problem this file fixes
---------------------------
exp_scale.py reported n=24 by MAPPING the n=12 spread->error ratio, because 2**24
has no oracle.  That mapping is exactly what the heavy-tail objection attacks: if
the marginal-contribution distribution is heavy-tailed, the 8-seed spread is a
weak statistic and the mapped error is optimistic.

So this file does the honest thing.  n = 8, 12 and 16 can all be FULLY
ENUMERATED (2**8, 2**12, 2**16 = 256, 4096, 65536 box runs), so at all three
points the true error is MEASURED against a real oracle.  The estimator and the
sampler are imported verbatim from exp_scale.py, so the comparison is
like-for-like and differences cannot come from the code path.

The ladder, and what each step isolates
---------------------------------------
  n = 8   2 sectors x 4 windows      time refinement (sectors held at 2)
  n = 16  2 sectors x 8 windows      <-- the new, fully measured point
  n = 24  3 sectors x 8 windows      adding a genuinely new sector

  n = 8 -> 16 therefore isolates TIME REFINEMENT, with the sector count fixed.
  n = 16 -> 24 isolates SECTOR ADDITION, with the time grid fixed.

Reads
  _kppwork/scen_v.csv   n=8  oracle (2 sectors x 4 windows)
  _kppwork/co_v.csv     n=12 oracle (3 sectors x 4 windows)
  _kppwork/n16_v_*.csv  n=16 oracle (2 sectors x 8 windows)   <- new
  _kppwork/scale_v.csv  n=24 sampled (no oracle)

Writes scaling_n16.txt
"""
import csv
import importlib.util
import json
import math
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
W = "/mnt/d/atmos/_kppwork"
OUT = "/mnt/d/atmos/_kppwork/scaling_n16.txt"

SCEN = ["r_sat", "r_lim1"]
BUDGETS = [64, 128, 256, 512, 1024]
# 24 seeds: the headline metric is the per-run error, and with only 8 seeds that
# statistic is too noisy to read a scaling law off (it came out non-monotone).
NSEED = 24

# ---- reuse the EXACT sampler / estimator from exp_scale.py (same provenance) --
_spec = importlib.util.spec_from_file_location("exp_scale", os.path.join(HERE, "exp_scale.py"))
_es = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_es)
gen = _es.gen            # identical RNG call sequence -> budgets are comparable
phi_perm = _es.phi_perm
mobius = _es.mobius


def load_simple(path, keep, mdim=None):
    """scen_v.csv / co_v.csv format: name like '<scen>_m<mask>'."""
    tab = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if keep and s not in keep:
                continue
            if row["v_ppb"] in ("", "NaN"):
                continue
            m = int(ms)
            if mdim and m >= mdim:
                continue
            tab.setdefault(s, {})[m] = float(row["v_ppb"])
    return tab


def load_n16(tag):
    """Merge the 12 worker shards of the n=16 enumeration."""
    vd = {}
    for K in range(1, 13):
        fn = os.path.join(W, "n16_v_%s_%d.csv" % (tag, K))
        if not os.path.exists(fn):
            continue
        with open(fn) as fh:
            for row in csv.DictReader(fh):
                if row.get("name", "").startswith("m"):
                    if row["v_ppb"] in ("", "NaN"):
                        continue
                    vd[int(row["name"][1:])] = float(row["v_ppb"])
    return vd


def mobius_fast(v, n):
    """Möbius transform in O(n 2^n) by in-place subset inversion.

    exp_scale.py evaluates the defining sum directly, which is O(3^n).  That is
    fine at n=12 (531k terms) but at n=16 it is 43M string-popcount terms, so
    this replaces it with the equivalent in-place inversion, 16 x 65536 = 1M
    cheap ops.  exact_shapley() cross-checks it against exp_scale.mobius at n=8
    before it is trusted.
    """
    f = list(v)
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                f[mask] -= f[mask ^ bit]
    return f


def exact_shapley(vd, n, verify=False):
    """Möbius transform then the closed-form Shapley linear functional."""
    vals = [vd[m] for m in range(1 << n)]
    mo = mobius_fast(vals, n)
    if verify:
        ref = mobius(vals, n)
        d = max(abs(a - b) for a, b in zip(mo, ref))
        print("    Möbius cross-check n=%d: max|fast - exp_scale| = %.3e" % (n, d), flush=True)
    pc = [bin(A).count("1") for A in range(1 << n)]
    true = [0.0] * n
    for A in range(1, 1 << n):
        c = mo[A] / pc[A]
        if c == 0.0:
            continue
        for i in range(n):
            if (A >> i) & 1:
                true[i] += c
    return true, mo


# The Möbius transform is O(3^n) and is by far the slowest step (3^16 = 43M
# iterations), while the sampler probing below is cheap.  Cache it to disk so a
# reaped session does not have to recompute it -- background jobs here get killed
# about 45 s after the launching session ends.
CACHE = os.path.join(W, "n16_exact_shapley.json")


def exact_shapley_cached(scen, n, vd):
    blob = {}
    if os.path.exists(CACHE):
        try:
            with open(CACHE) as fh:
                blob = json.load(fh)
        except Exception:
            blob = {}
    key = "%s|%d" % (scen, n)
    if key in blob and blob[key].get("n") == n:
        return blob[key]["true"], blob[key]["scale"]
    print("    computing exact Shapley for %s n=%d ..." % (scen, n), flush=True)
    true, _ = exact_shapley(vd, n, verify=(n <= 8))
    scale = max(abs(x) for x in true)
    blob[key] = {"n": n, "true": true, "scale": scale}
    with open(CACHE, "w") as fh:
        json.dump(blob, fh)
    return true, scale


def probe(scen, vd, n, budgets=BUDGETS, nseed=NSEED):
    """Sampler precision vs a REAL oracle.

    Three quantities, and it matters which is which:

      err(mean)  = max_i |mean_over_seeds(est_i) - true_i| / scale
                   This measures BIAS only.  It shrinks toward zero as seeds are
                   added regardless of budget, so it is the WRONG precision
                   metric -- and with 8 seeds it is noisy enough to come out
                   non-monotone in the budget.

      err(run)   = mean over seeds of max_i |est_i - true_i| / scale
                   How far off ONE sampler run typically is.  This is what a user
                   actually experiences, and it is the metric the scaling law
                   should be read from.

      spread     = max_i std_over_seeds(est_i) / scale
                   The estimator's own dispersion; the only one available at
                   n=24, where no oracle exists.
    """
    true, scale = exact_shapley_cached(scen, n, vd)
    # A budget larger than 2^n is UNREACHABLE: gen() keeps drawing permutations
    # until it has that many distinct coalitions, so it would spin forever.  This
    # is what hung the earlier runs at n=8 with B=512, and took the WSL VM down
    # with them.  Cap the budget at the size of the universe.
    budgets = [b for b in budgets if b <= (1 << n)]
    rows = []
    for B in budgets:
        ests = []
        used = 0
        for sd in range(nseed):
            rng = random.Random(1000 * B + sd)
            orders, seen = gen(rng, B, n)
            if any(m not in vd for m in seen):
                break
            ests.append(phi_perm({m: vd[m] for m in seen}, orders, n))
            used = len(seen)
        if not ests:
            rows.append((B, None, None, None, None))
            continue
        nseed_used = len(ests)
        mu = [sum(e[i] for e in ests) / nseed_used for i in range(n)]
        sd_i = []
        for i in range(n):
            xs = [e[i] for e in ests]
            m = sum(xs) / len(xs)
            sd_i.append(math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) if len(xs) > 1 else 0.0)
        err_mean = max(abs(mu[i] - true[i]) for i in range(n)) / scale
        single = [max(abs(e[i] - true[i]) for i in range(n)) / scale for e in ests]
        err_run = sum(single) / len(single)
        sprd = max(sd_i) / scale
        rows.append((B, used, err_run, sprd, err_mean))
        print("  %-8s n=%2d B=%5d runs=%5d  err(run)=%7.3f%%  spread=%7.3f%%  err(mean)=%6.3f%%"
              % (scen, n, B, used, 100 * err_run, 100 * sprd, 100 * err_mean), flush=True)
    return rows, true, scale


def main():
    L = []
    add = L.append
    add("Measured scaling of the permutation+antithetic sampler in n")
    add("estimator and sampler imported verbatim from exp_scale.py")
    add("error and spread are both normalised by max_i |phi_i|")
    add("")
    add("oracle sizes: n=8 -> 256, n=12 -> 4096, n=16 -> 65536 (all enumerated)")
    add("              n=24 -> none (only spread is available)")
    add("")

    print("loading oracles ...", flush=True)
    t8 = load_simple(os.path.join(W, "scen_v.csv"), set(SCEN))
    t12 = load_simple(os.path.join(W, "co_v.csv"), set(SCEN), mdim=4096)
    t16 = {s: load_n16(s) for s in SCEN}
    for s in SCEN:
        add("n=8  %-8s %6d values (expect 256)" % (s, len(t8.get(s, {}))))
        add("n=12 %-8s %6d values (expect 4096)" % (s, len(t12.get(s, {}))))
        add("n=16 %-8s %6d values (expect 65536)" % (s, len(t16.get(s, {}))))
    add("")

    results = {}
    for s in SCEN:
        add("=" * 78)
        add("scenario %s" % s)
        add("=" * 78)
        print("scenario %s" % s, flush=True)
        for n, tab in ((8, t8), (12, t12), (16, t16)):
            vd = tab.get(s, {})
            if len(vd) < (1 << n):
                add("  n=%d: incomplete (%d of %d) -- skipped" % (n, len(vd), 1 << n))
                print("  n=%d incomplete (%d)" % (n, len(vd)), flush=True)
                continue
            rows, true, scale = probe(s, vd, n)
            results[(s, n)] = dict(rows=rows, scale=scale, true=true)
            add("")
            add("  n=%d  max|phi| = %.3f ppb   oracle = %d coalitions" % (n, scale, 1 << n))
            add("  %6s %7s %11s %11s %11s %9s"
                % ("B", "runs", "err(run)", "spread", "err(mean)", "err*run^.5"))
            for (B, used, err_run, sprd, err_mean) in rows:
                if err_run is None:
                    add("  %6d  (insufficient)" % B)
                    continue
                add("  %6d %7d %10.3f%% %10.3f%% %10.3f%% %9.4f"
                    % (B, used, 100 * err_run, 100 * sprd, 100 * err_mean,
                       err_run * math.sqrt(used)))
            add("")

    # ---------------- the headline: runs needed for a fixed target -----------
    TARGET = 0.015
    add("=" * 78)
    add("runs needed for a %.1f%% normalised error (log-log interpolation)" % (100 * TARGET))
    add("=" * 78)
    add("  %-8s %4s %10s %10s" % ("scen", "n", "runs", "vs n=12"))
    need = {}
    for s in SCEN:
        for n in (8, 12, 16):
            rows = [r for r in results.get((s, n), {}).get("rows", [])
                    if r[2] is not None and r[2] > 0]
            if len(rows) < 2:
                continue
            # invert err ~ C * runs**-a using the two nearest points bracketing target
            lo = hi = None
            for r in rows:
                if r[2] <= TARGET:
                    hi = r
                if lo is None or r[2] < lo[2]:
                    pass
            pts = sorted(rows, key=lambda r: r[1])
            est = None
            for a, b in zip(pts, pts[1:]):
                if a[2] >= TARGET >= b[2]:
                    la, lb = math.log(a[2]), math.log(b[2])
                    ra, rb = math.log(a[1]), math.log(b[1])
                    t = (math.log(TARGET) - la) / (lb - la)
                    est = math.exp(ra + t * (rb - ra))
                    break
            if est is None:
                est = math.nan
            need[(s, n)] = est
            base = need.get((s, 12), math.nan)
            add("  %-8s %4d %10s %10s"
                % (s, n, ("%.0f" % est) if est == est else "n/a",
                   ("%.2fx" % (est / base)) if (est == est and base == base) else "n/a"))
    add("")
    add("enumeration for comparison: n=8 256, n=12 4096, n=16 65536, n=24 16777216")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
