#!/usr/bin/env python3
"""Dual-metric re-derivation of the sampling exhibits (audit fix).

WHY
  Every released sampling output -- attrib_bench.txt, bench_n16.txt,
  scaling_n16.txt -- was computed on the d3_mean column (v_ppb), while the
  paper's primary metric is MDA8.  The MDA8 oracles exist for every n:
    n = 16  enumerated, both columns in n16_v_<scen>_m8.csv
    n = 12  co_v_m8_*.csv shards (3 sectors x 4 windows)
    n =  8  read off the n=16 table at union-of-pairs masks (identical
            forcing history => identical trajectory => any functional agrees),
            licensed by a v_ppb cross-check against scen_v.csv each run
  So both metrics can be recomputed with zero chemistry runs.

WHAT
  PART A  sampler cost law (perm + antithetic), both metrics, n = 8/12/16:
          err(run) at matched budget + honest same-budget ratios + runs for a
          target error (reported n/a when the interpolation is not bracketed).
  PART B  arm comparison (kernelshap / perm / perm_anti / perm_strat / OURS)
          at n = 8 and n = 16, both metrics.

OUTPUTS
  scaling_dual_metrics.txt
  attrib_dual_metrics.txt

NOTE ON HONESTY
  Ratios are printed per (scenario, budget) pair, never pooled into a single
  headline number: the released data shows 1.6-1.8x on r_sat but a wider and
  budget-dependent range on r_lim1, and the earlier "same accuracy costs
  2.5-3x" statement has no artifact behind it (it is dropped).
"""
import csv
import importlib.util
import math
import os
import random
import statistics as st

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fn):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, fn))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_es = _load("es", "exp_scale.py")          # gen / phi_perm (identical RNG stream)
gen, phi_perm = _es.gen, _es.phi_perm
_b = _load("bn", "exp_n16_bench.py")       # arms + exact_phi + Budget
_b.SEEDS = 50

SCEN = ["r_sat", "r_lim1"]
BUDGETS = [64, 128, 256, 512, 1024]
NSEED = 24
TARGET = 0.015
METRICS = [("v_ppb", "d3_mean (secondary)"), ("v_mda8", "MDA8 (primary)")]
N8BUD = [8, 16, 24, 32, 48, 64, 96, 128, 160, 192]
N16BUD = [64, 128, 256, 512]


# ------------------------------------------------------------------ oracles
def load_tagged(path, col, keep=None, mdim=None):
    """scen_v.csv / co_v.csv / co_v_m8_k.csv format: name '<scen>_m<mask>'."""
    tab = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            nm = row.get("name", "")
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if not ms.isdigit() or (keep and s not in keep):
                continue
            if row.get(col, "") in ("", "NaN"):
                continue
            m = int(ms)
            if mdim is not None and m >= mdim:
                continue
            tab.setdefault(s, {})[m] = float(row[col])
    return tab


def load_n16(scen):
    """n16_v_<scen>_m8.csv -> {mask: (v48, v_mda8)}"""
    fn = os.path.join(HERE, "n16_v_%s_m8.csv" % scen)
    d = {}
    with open(fn) as fh:
        for row in csv.DictReader(fh):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            d[int(row["name"][1:])] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return d


def embed8(m8):
    """6 h mask (NOx 0-3, VOC 4-7) -> 3 h mask (NOx 0-7, VOC 8-15)."""
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16


def build_oracles(L):
    """returns oracles[metric][scen][n] = {mask: value}, with checks logged."""
    ref8 = load_tagged(os.path.join(HERE, "scen_v.csv"), "v_ppb",
                       keep=set(SCEN), mdim=256)
    o12 = {"v_ppb": load_tagged(os.path.join(HERE, "co_v.csv"), "v_ppb",
                                keep=set(SCEN), mdim=4096)}
    m12 = {}
    for K in range(1, 7):
        fn = os.path.join(HERE, "co_v_m8_%d.csv" % K)
        if os.path.exists(fn):
            for s, d in load_tagged(fn, "v_mda8", keep=set(SCEN), mdim=4096).items():
                m12.setdefault(s, {}).update(d)
    o12["v_mda8"] = m12

    orc = {k: {} for k, _ in METRICS}
    for s in SCEN:
        v16 = load_n16(s)
        d48_8, d8_8 = {}, {}
        worst = 0.0
        for m8 in range(256):
            p48, p8 = v16[embed8(m8)]
            d48_8[m8], d8_8[m8] = p48, p8
            if m8 in ref8.get(s, {}):
                worst = max(worst, abs(p48 - ref8[s][m8]))
        L.append("  [%s] n=8 sub-lattice check (derived v_ppb vs scen_v.csv): "
                 "max|dv48| = %.3e  %s"
                 % (s, worst, "OK" if worst < 1e-9 else "!! EXCEEDS TOL"))
        orc["v_ppb"][s] = {8: dict(d48_8), 12: dict(o12["v_ppb"].get(s, {})),
                           16: {m: p[0] for m, p in v16.items()}}
        orc["v_mda8"][s] = {8: dict(d8_8), 12: dict(o12["v_mda8"].get(s, {})),
                            16: {m: p[1] for m, p in v16.items()}}
    return orc


# ------------------------------------------------------------- exact Shapley
def mobius_fast(vals, n):
    f = list(vals)
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                f[mask] -= f[mask ^ bit]
    return f


def exact_shapley(vd, n):
    v0 = vd[0]
    vals = [vd[m] - v0 for m in range(1 << n)]
    mo = mobius_fast(vals, n)
    out = [0.0] * n
    for A in range(1, 1 << n):
        c = mo[A] / bin(A).count("1")
        for i in range(n):
            if (A >> i) & 1:
                out[i] += c
    return out


# ------------------------------------------------------------------- part A
def scaling_tables(orc, L):
    ratios = {}
    for col, mname in METRICS:
        L.append("")
        L.append("=" * 78)
        L.append("PART A -- sampler cost law, metric = %s" % mname)
        L.append("=" * 78)
        L.append("  estimator = permutation + antithetic, %d seeds, "
                 "error = max_i |est-true| / max_i |true| (single-run)" % NSEED)
        for s in SCEN:
            L.append("")
            L.append("  scenario %s" % s)
            L.append("  %6s %9s %11s %11s" % ("n", "B", "runs", "err(run)"))
            errs = {}
            for n in (8, 12, 16):
                vd = orc[col][s][n]
                if len(vd) < (1 << n):
                    L.append("  %6d  (incomplete: %d of %d)" % (n, len(vd), 1 << n))
                    continue
                true = exact_shapley(vd, n)
                scale = max(abs(x) for x in true)
                for B in BUDGETS:
                    if B > (1 << n):
                        continue
                    es = []
                    for sd in range(NSEED):
                        rng = random.Random(1000 * B + sd)
                        orders, seen = gen(rng, B, n)
                        est = phi_perm({m: vd[m] for m in seen}, orders, n)
                        es.append(max(abs(est[i] - true[i])
                                      for i in range(n)) / scale)
                    errs[(n, B)] = st.median(es)
                    L.append("  %6d %9d %11d %10.3f%%"
                             % (n, B, B, 100 * errs[(n, B)]))
            # honest ratio reporting: per-budget, no pooling
            L.append("    same-budget error ratios (clean ladders only):")
            for B in BUDGETS:
                parts = []
                for lo, hi, lab in ((8, 12, "8->12 sector+1"),
                                    (8, 16, "8->16 time-refined"),
                                    (12, 16, "12->16 NON-NESTED")):
                    a, b = errs.get((lo, B)), errs.get((hi, B))
                    if a and b:
                        parts.append("%s %.2fx" % (lab, b / a))
                        ratios.setdefault((col, s, lab), []).append(b / a)
                if parts:
                    L.append("      B=%4d  %s" % (B, "   ".join(parts)))
        # runs for target error
        L.append("")
        L.append("  runs needed for %.1f%% error (log-log interpolation; "
                 "'n/a' = not bracketed)" % (100 * TARGET))
        L.append("  %-8s %4s %10s" % ("scen", "n", "runs"))
        for s in SCEN:
            for n in (8, 12, 16):
                vd = orc[col][s][n]
                if len(vd) < (1 << n):
                    continue
                pts = []
                true = exact_shapley(vd, n)
                scale = max(abs(x) for x in true)
                for B in BUDGETS:
                    if B > (1 << n):
                        continue
                    es = []
                    for sd in range(NSEED):
                        rng = random.Random(1000 * B + sd)
                        orders, seen = gen(rng, B, n)
                        est = phi_perm({m: vd[m] for m in seen}, orders, n)
                        es.append(max(abs(est[i] - true[i])
                                      for i in range(n)) / scale)
                    pts.append((B, st.median(es)))
                pts.sort()
                run = None
                for (b1, e1), (b2, e2) in zip(pts, pts[1:]):
                    if e1 >= TARGET >= e2 and e1 > e2:
                        run = math.exp(math.log(b1)
                                       + (math.log(TARGET) - math.log(e1))
                                       * (math.log(b2) - math.log(b1))
                                       / (math.log(e2) - math.log(e1)))
                        break
                L.append("  %-8s %4d %10s" % (s, n, "%.0f" % run if run else "n/a"))
    return ratios


# ------------------------------------------------------------------- part B
def arms_table(orc, L):
    for col, mname in METRICS:
        L.append("")
        L.append("=" * 78)
        L.append("PART B -- attribution arms, metric = %s" % mname)
        L.append("=" * 78)
        for s in SCEN:
            for n, buds in ((8, N8BUD), (16, N16BUD)):
                vd = orc[col][s][n]
                if len(vd) < (1 << n):
                    L.append("  %s n=%d: oracle incomplete" % (s, n))
                    continue
                res, scale = _b.attrib_run(vd, n, buds, "%s_n%d" % (s, n))
                L.append("")
                L.append("  %s  n=%d  max|phi| = %.4f ppb  (%d arms x %d seeds)"
                         % (s, n, scale, len(_b.ARMS), _b.SEEDS))
                L.append("  %-12s %8s %8s %8s %8s"
                         % ("arm", "20%", "10%", "5%", "2%"))
                for name, _ in _b.ARMS:
                    row = [name]
                    for eps in (0.20, 0.10, 0.05, 0.02):
                        r = _b.reach(res[name], eps)
                        row.append(str(r) if r else ">%d" % max(buds))
                    L.append("  %-12s %8s %8s %8s %8s" % tuple(row))


def main():
    L = []
    add = L.append
    add("Dual-metric re-derivation of the released sampling exhibits")
    add("Rationale: attrib_bench.txt / bench_n16.txt / scaling_n16.txt were all")
    add("computed on the d3_mean column; this recomputes them on both metrics")
    add("with zero chemistry runs (all values read from released oracle tables).")
    add("")
    add("oracle provenance:")
    add("  n=16  n16_v_<scen>_m8.csv (enumerated, both columns)")
    add("  n=12  co_v.csv (d3_mean) / co_v_m8_*.csv shards (MDA8)")
    add("  n=8   scen_v.csv (d3_mean) / sub-lattice of the n=16 MDA8 table")
    orc = build_oracles(L)
    ratios = scaling_tables(orc, L)
    arms_table(orc, L)

    add("")
    add("=" * 78)
    add("SUMMARY -- same-budget error ratios, per (metric, scenario, ladder)")
    add("=" * 78)
    for (col, s, lab), xs in sorted(ratios.items()):
        add("  %-9s %-7s %-16s n=%d  median %.2fx  range %.2f-%.2f"
            % (col, s, lab, len(xs), st.median(xs), min(xs), max(xs)))
    add("")
    add("Reading: the 8->12 (sector addition) and 8->16 (time refinement) ladders")
    add("are clean nested comparisons; 12->16 is NOT nested (different games) and")
    add("is reported only for completeness.  Ratios are budget-dependent, so no")
    add("single pooled factor is quoted.")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "scaling_dual_metrics.txt"), "w") as fh:
        fh.write(txt + "\n")
    with open(os.path.join(HERE, "attrib_dual_metrics.txt"), "w") as fh:
        fh.write("See scaling_dual_metrics.txt PART B (single run writes both).\n")
    print("\nwrote scaling_dual_metrics.txt")


if __name__ == "__main__":
    main()