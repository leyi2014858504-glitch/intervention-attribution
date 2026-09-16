#!/usr/bin/env python3
"""Does HDDM (second-order Taylor) rescue the local method?

WHY
  A reviewer's strongest objection to "local methods fail" is: "we would use
  HDDM, not DDM -- second order captures curvature."  That must be answered
  with numbers, not argument.  The 109 runs behind results_local_vs_global.txt
  already contain everything needed: at the all-on base x=1 we have f(1) (base),
  f(1+delta) (plus) and f(1-delta) (minus), so BOTH derivatives are already
  measurable -- no new box runs.

  DDM-1  (what CMAQ-DDM / adjoint give)   :  delta_hat = -f'(1)
  HDDM-2 (second-order Taylor at x=1)     :  delta_hat = -f'(1) + f''(1)/2
  truth  (finite switch-off, x: 1 -> 0)   :  delta     = f(0) - f(1)

  Both are extrapolations from the base point across the full unit interval;
  the question is whether the parabola reaches the endpoint better than the
  tangent does.  Sign correctness and |error| are reported for both, under
  both metrics, for 24 three-hour windows and 12 six-hour blocks.

  Note on the metric: this is EXACTLY what the paper needs -- the reviewer
  claim is about the method, so the answer must hold under both metrics or be
  reported as metric-dependent.

Writes results_hddm_check.txt
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "g4runs", "scale")
OUT = os.path.join(HERE, "..", "results_hddm_check.txt")
DELTA = 0.05
SECTORS = ["NOXW", "VOCW", "COW"]
NW = 8


def metrics(path):
    """(v48, mda8) from a traj, same rules as _local_vs_global.py"""
    t, o3 = [], []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            p = ln.split()
            t.append(float(p[0]))
            o3.append(float(p[1]))
    sel = [i for i in range(len(t)) if t[i] >= 48.0]
    v48 = sum(o3[i] for i in sel) / len(sel)
    dt = t[1] - t[0]
    w = max(1, int(round(8.0 / dt)))
    best = float("-inf")
    for i in range(w - 1, len(o3)):
        if t[i] < 48.0:
            continue
        a = sum(o3[i - w + 1:i + 1]) / w
        if a > best:
            best = a
    return v48, best


def load(tag):
    return metrics(os.path.join(RUNS, "lvg_%s.traj" % tag))


def main():
    base = load("base")
    L = []
    add = L.append
    add("Does HDDM (2nd-order Taylor) rescue the local method?")
    add("DDM-1 : delta_hat = -f'(1)              [tangent extrapolation]")
    add("HDDM-2: delta_hat = -f'(1) + f''(1)/2   [parabola extrapolation]")
    add("truth : delta     = f(0) - f(1)         [finite switch-off]")
    add("derivatives by central difference at x=1, delta = %.2f" % DELTA)
    add("base v48 = %.6f   base MDA8 = %.6f" % base)
    add("")

    rows = []
    for label, sfx in (("3h windows", ""), ("6h blocks", "blk")):
        for sec in SECTORS:
            nw = NW if not sfx else 4
            for i in range(nw):
                try:
                    vp = load("%s_%sp%d" % (sec, sfx, i))
                    vm = load("%s_%sm%d" % (sec, sfx, i))
                    vo = load("%s_%soff%d" % (sec, sfx, i))
                except FileNotFoundError:
                    continue
                span = 3 if not sfx else 6
                lo = 3 * i if not sfx else 6 * i
                r = dict(kind=label, sec=sec, win="%d-%d" % (lo, lo + span))
                for tag, j in (("48", 0), ("8", 1)):
                    d1 = (vp[j] - vm[j]) / (2 * DELTA)
                    d2 = (vp[j] - 2 * base[j] + vm[j]) / (DELTA ** 2)
                    act = vo[j] - base[j]
                    p1 = -d1
                    p2 = -d1 + 0.5 * d2
                    r["d1" + tag], r["d2" + tag] = d1, d2
                    r["act" + tag], r["p1" + tag], r["p2" + tag] = act, p1, p2
                    r["e1" + tag] = p1 - act
                    r["e2" + tag] = p2 - act
                    r["s1" + tag] = (p1 * act) > 0
                    r["s2" + tag] = (p2 * act) > 0
                rows.append(r)

    for tag, mname in (("8", "MDA8 (primary)"), ("48", "d3_mean (secondary)")):
        add("=" * 94)
        add("METRIC = %s" % mname)
        add("=" * 94)
        add("%-11s %-8s %8s %8s | %8s %8s | %8s %8s | %s" %
            ("sector", "window", "D1", "D2", "DDM1pred", "HDDM2pr",
             "err1", "err2", "sign 1/2"))
        for r in rows:
            add("%-11s %-8s %8.3f %8.3f | %8.3f %8.3f | %+8.3f %+8.3f | %s/%s" %
                (r["sec"], r["win"], r["d1" + tag], r["d2" + tag],
                 r["p1" + tag], r["p2" + tag],
                 r["e1" + tag], r["e2" + tag],
                 "ok " if r["s1" + tag] else "BAD",
                 "ok " if r["s2" + tag] else "BAD"))
        add("")
        bad1 = [r for r in rows if not r["s1" + tag]]
        bad2 = [r for r in rows if not r["s2" + tag]]
        worse = [r for r in rows if abs(r["e2" + tag]) > abs(r["e1" + tag])]
        # normalised by the size of the true effect
        den = max(abs(r["act" + tag]) for r in rows) or 1.0
        med1 = sorted(abs(r["e1" + tag]) for r in rows)[len(rows) // 2] / den
        med2 = sorted(abs(r["e2" + tag]) for r in rows)[len(rows) // 2] / den
        mx1 = max(abs(r["e1" + tag]) for r in rows) / den
        mx2 = max(abs(r["e2" + tag]) for r in rows) / den
        add("  sign errors : DDM-1 %d/%d   HDDM-2 %d/%d" %
            (len(bad1), len(rows), len(bad2), len(rows)))
        add("  |err|/max|effect| : median  DDM-1 %.3f -> HDDM-2 %.3f" % (med1, med2))
        add("                       max     DDM-1 %.3f -> HDDM-2 %.3f" % (mx1, mx2))
        add("  cases where HDDM-2 is WORSE than DDM-1: %d/%d" % (len(worse), len(rows)))
        for r in bad2:
            add("    HDDM-2 still flips sign: %s %s  pred2=%+.3f truth=%+.3f" %
                (r["sec"], r["win"], r["p2" + tag], r["act" + tag]))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("wrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
