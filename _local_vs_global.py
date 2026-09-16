#!/usr/bin/env python3
"""
Local derivative (what DDM/adjoint can give) vs global finite attribution.

WHY THIS EXPERIMENT EXISTS
--------------------------
The first question a reviewer asks about a Shapley/Mobius attribution paper is
"why not just use CMAQ-DDM-3D or the adjoint?" Those tools are cheap (one run,
independent of the number of players) but they are LOCAL: they linearise at a
base state and report dC/dtheta there.

Our object is different in kind. We attribute FINITE interventions -- "shut this
sector off in this 3-hour window" -- which is a multiplier going from 1 to 0, not
an infinitesimal nudge.

The two coincide only if the response is linear in that multiplier. So we measure
both on the same model:

    D_i     local, semi-normalised first derivative at the all-on base state
            D_i = (v(1+delta) - v(1-delta)) / (2*delta)          [ppb per unit]
            this is exactly what DDM-1 reports (theta * dC/dtheta at theta=1)

    Delta_i finite intervention: turn the window fully off
            Delta_i = v(multiplier_i = 0) - v(all on)            [ppb]

If the window were linear, Delta_i = -D_i (a full unit decrease along the local
slope). The ratio Delta_i / (-D_i) is therefore a direct measure of how much of
the effect a local method cannot see, and a SIGN disagreement means a local
method gets the direction wrong, not merely the magnitude.

Run count: 1 base + 3 sectors x 8 windows x 3 (off/plus/minus) = 73,
plus the same aggregated to the four 6-hour blocks (+36). 109 box runs.

CO is included deliberately: its response turns out to be linear to within 2 %, which
makes it the control case -- a local method is exact on CO, but CO is also irrelevant
(<=0.14 ppb), while NOx, 40-200x larger, is where the local sign fails.

METRIC: every trajectory is scored under BOTH the d3_mean (v48) statistic and the
regulatory MDA8 statistic, using the same 109 runs.  MDA8 is the paper's primary
metric (matching the v_mda8 column of the enumeration sweep); v48 is kept so the
sign structure can be shown to be metric-independent and the magnitude amplification
visible.  The two are computed identically to exp_mda8.mean_and_mda8(t0=48).
"""
import sys, os, csv, json, math
sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box

D = "/mnt/d/atmos"
DELTA = 0.05                      # 5 % central difference
SECTORS = ["NOXW", "VOCW", "COW"]   # CO included so the comparison spans a weak, near-linear player too
NW = 8                            # eight 3-hour windows
ALL1 = [1.0] * NW
NAMES = [f"[{3*i:02d},{3*i+3:02d})" for i in range(NW)]


def run(overrides, tag):
    """One box run on the 8-window backend; returns (v48, mda8) from the traj."""
    cfg = dict(overrides)
    cfg.setdefault("NOXW", ALL1)
    cfg.setdefault("VOCW", ALL1)
    cfg.setdefault("COW", ALL1)
    m = kpp_box.run_scenario(cfg, "lvg_" + tag, "scale")
    return _metrics(m["traj"])


def _v48(traj):
    """Mean O3 over t >= 48 h. This is the frozen V-gate statistic."""
    xs = [r[1] for r in traj if r[0] >= 48.0]
    return sum(xs) / len(xs)


def _metrics(traj):
    """Return (v48, mda8): d3_mean and max rolling 8h mean over windows ending
    at t >= 48 h -- identical rule to the sweep's v_mda8 column."""
    tt = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    sel = [i for i in range(len(tt)) if tt[i] >= 48.0]
    v48 = sum(o3[i] for i in sel) / len(sel) if sel else float("nan")
    if len(tt) < 2:
        return v48, float("nan")
    dt = tt[1] - tt[0]
    w = max(1, int(round(8.0 / dt)))
    best = float("-inf")
    for i in range(w - 1, len(o3)):
        if tt[i] < 48.0:
            continue
        a = sum(o3[i - w + 1:i + 1]) / w
        if a > best:
            best = a
    return v48, (best if best > float("-inf") else float("nan"))


def _mk(sec, hours, base, v_off, v_p, v_m):
    """One comparison -> a row carrying BOTH metrics.  idx 0 = v48, 1 = mda8."""
    d = dict(sector=sec, hours=hours)
    for k, tag in ((0, "48"), (1, "8")):
        Dloc = (v_p[k] - v_m[k]) / (2 * DELTA)
        Dfin = v_off[k] - base[k]
        d["local" + tag] = Dloc
        d["finite" + tag] = Dfin
        d["ratio" + tag] = Dfin / (-Dloc) if abs(Dloc) > 1e-9 else float("nan")
        d["sign" + tag] = (Dloc * (-Dfin)) > 0
    return d


def _table(fh, rs, tag, label):
    """Write one metric's table; return (nbad, n)."""
    fh.write(label + "\n")
    fh.write(f"{'sector':<6}{'hours':<10}{'D_local':>11}{'D_finite':>11}"
             f"{'-D_local':>11}{'ratio':>9}{'sign_ok':>9}\n")
    for r in rs:
        fh.write(f"{r['sector']:<6}{r['hours']:<10}{r['local'+tag]:11.3f}"
                 f"{r['finite'+tag]:11.3f}{-r['local'+tag]:11.3f}"
                 f"{r['ratio'+tag]:9.2f}{str(r['sign'+tag]):>9}\n")
    nbad = sum(1 for r in rs if not r["sign" + tag])
    fh.write(f"  -> sign disagreements: {nbad} of {len(rs)}\n\n")
    return nbad, len(rs)


def main():
    rows = []

    base = run({}, "base")
    print(f"base v48 = {base[0]:.6f}  mda8 = {base[1]:.6f}   "
          f"(V-gate v48 reference 21.182293)")

    # ---------------- 3-hour windows ----------------
    for sec in SECTORS:
        for i in range(NW):
            off = list(ALL1); off[i] = 0.0
            plus = list(ALL1); plus[i] = 1.0 + DELTA
            minus = list(ALL1); minus[i] = 1.0 - DELTA
            v_off = run({sec: off}, f"{sec}_off{i}")
            v_p = run({sec: plus}, f"{sec}_p{i}")
            v_m = run({sec: minus}, f"{sec}_m{i}")
            r = _mk(sec, f"{3*i}-{3*i+3}", base, v_off, v_p, v_m)
            rows.append(r)
            print(f"  {sec:<5}{NAMES[i]:<10} "
                  f"[v48] loc={r['local48']:8.3f} fin={r['finite48']:8.3f} "
                  f"rat={r['ratio48']:6.2f} ok={r['sign48']:1}  "
                  f"[M8] loc={r['local8']:8.3f} fin={r['finite8']:8.3f} "
                  f"rat={r['ratio8']:6.2f} ok={r['sign8']}")

    # ---------------- the same, aggregated to four 6-hour blocks ----------------
    print("\n  --- aggregated to four 6-hour blocks (pairs of 3h windows) ---")
    blocks = []
    for sec in SECTORS:
        for b in range(4):
            idx = [2 * b, 2 * b + 1]
            off = list(ALL1)
            plus = list(ALL1)
            minus = list(ALL1)
            for i in idx:
                off[i] = 0.0
                plus[i] = 1.0 + DELTA
                minus[i] = 1.0 - DELTA
            v_off = run({sec: off}, f"{sec}_blkoff{b}")
            v_p = run({sec: plus}, f"{sec}_blkp{b}")
            v_m = run({sec: minus}, f"{sec}_blkm{b}")
            r = _mk(sec, f"{6*b}-{6*b+6}", base, v_off, v_p, v_m)
            blocks.append(r)
            print(f"  {sec:<5}[{6*b:02d},{6*b+6:02d})   "
                  f"[v48] loc={r['local48']:8.3f} fin={r['finite48']:8.3f} "
                  f"rat={r['ratio48']:6.2f} ok={r['sign48']:1}  "
                  f"[M8] loc={r['local8']:8.3f} fin={r['finite8']:8.3f} "
                  f"rat={r['ratio8']:6.2f} ok={r['sign8']}")

    # ---------------- report ----------------
    out = os.path.join(D, "results_local_vs_global.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("Local derivative vs global finite attribution -- BOTH metrics\n")
        fh.write(f"central-difference delta = {DELTA}\n")
        fh.write(f"base v48 = {base[0]:.6f}  (must equal 21.182293)   "
                 f"base MDA8 = {base[1]:.6f}\n\n")
        fh.write("D_local  = semi-normalised first derivative at the all-on base\n")
        fh.write("           (= what CMAQ-DDM-1 / the adjoint report at theta=1)\n")
        fh.write("D_finite = change in the metric from turning the window fully OFF\n")
        fh.write("           linear response would give D_finite = -D_local, ratio = 1\n\n")
        for tag, mname in (("8", "MDA8 (primary metric)"),
                           ("48", "d3_mean / v48 (secondary)")):
            fh.write("#" * 66 + "\n")
            fh.write("# " + mname + "\n")
            fh.write("#" * 66 + "\n")
            _table(fh, rows, tag, "=== 3-hour windows ===")
            _table(fh, blocks, tag, "=== 6-hour blocks ===")
    print()
    for tag, mname in (("8", "MDA8"), ("48", "v48")):
        nb3 = sum(1 for r in rows if not r["sign" + tag])
        nb6 = sum(1 for r in blocks if not r["sign" + tag])
        print(f"{mname}: sign disagreements 3h {nb3}/{len(rows)}, "
              f"6h {nb6}/{len(blocks)}")
    print("wrote", out)


if __name__ == "__main__":
    main()
