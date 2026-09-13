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
    """One box run on the 8-window backend; returns the V-gate metric v48."""
    cfg = dict(overrides)
    cfg.setdefault("NOXW", ALL1)
    cfg.setdefault("VOCW", ALL1)
    cfg.setdefault("COW", ALL1)
    m = kpp_box.run_scenario(cfg, "lvg_" + tag, "scale")
    return _v48(m["traj"])


def _v48(traj):
    """Mean O3 over t >= 48 h. This is the frozen V-gate statistic."""
    xs = [r[1] for r in traj if r[0] >= 48.0]
    return sum(xs) / len(xs)


def main():
    rows = []

    base = run({}, "base")
    print(f"base v48 = {base:.6f}   (V-gate reference 21.182293)")

    # ---------------- 3-hour windows ----------------
    for sec in SECTORS:
        for i in range(NW):
            off = list(ALL1); off[i] = 0.0
            plus = list(ALL1); plus[i] = 1.0 + DELTA
            minus = list(ALL1); minus[i] = 1.0 - DELTA
            v_off = run({sec: off}, f"{sec}_off{i}")
            v_p = run({sec: plus}, f"{sec}_p{i}")
            v_m = run({sec: minus}, f"{sec}_m{i}")
            Dloc = (v_p - v_m) / (2 * DELTA)
            Dfin = v_off - base
            rows.append(dict(sector=sec, win=NAMES[i], hours=f"{3*i}-{3*i+3}",
                             local=Dloc, finite=Dfin,
                             ratio=(Dfin / (-Dloc) if abs(Dloc) > 1e-9 else float("nan")),
                             sign_ok=(Dloc * (-Dfin) > 0)))
            print(f"  {sec:<5}{NAMES[i]:<10} local={Dloc:8.3f}  "
                  f"finite={Dfin:8.3f}  -local={-Dloc:8.3f}  "
                  f"ratio={rows[-1]['ratio']:7.2f}  sign_ok={rows[-1]['sign_ok']}")

    # ---------------- the same, aggregated to four 6-hour blocks ----------------
    print("\n  --- aggregated to four 6-hour blocks (pairs of 3h windows) ---")
    blocks = []
    for sec in SECTORS:
        for b in range(4):
            idx = [2 * b, 2 * b + 1]
            off = list(ALL1)
            for i in idx:
                off[i] = 0.0
            plus = list(ALL1)
            minus = list(ALL1)
            for i in idx:
                plus[i] = 1.0 + DELTA
                minus[i] = 1.0 - DELTA
            v_off = run({sec: off}, f"{sec}_blkoff{b}")
            v_p = run({sec: plus}, f"{sec}_blkp{b}")
            v_m = run({sec: minus}, f"{sec}_blkm{b}")
            Dloc = (v_p - v_m) / (2 * DELTA)
            Dfin = v_off - base
            blocks.append(dict(sector=sec, hours=f"{6*b}-{6*b+6}", local=Dloc,
                               finite=Dfin,
                               ratio=(Dfin / (-Dloc) if abs(Dloc) > 1e-9 else float("nan")),
                               sign_ok=(Dloc * (-Dfin) > 0)))
            print(f"  {sec:<5}[{6*b:02d},{6*b+6:02d})   local={Dloc:8.3f}  "
                  f"finite={Dfin:8.3f}  -local={-Dloc:8.3f}  "
                  f"ratio={blocks[-1]['ratio']:7.2f}  sign_ok={blocks[-1]['sign_ok']}")

    # ---------------- report ----------------
    out = os.path.join(D, "results_local_vs_global.txt")
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("Local derivative vs global finite attribution\n")
        fh.write(f"central-difference delta = {DELTA}\n")
        fh.write(f"base v48 = {base:.6f}  (must equal 21.182293)\n\n")
        fh.write("D_local  = semi-normalised first derivative at the all-on base\n")
        fh.write("           (= what CMAQ-DDM-1 / the adjoint report at theta=1)\n")
        fh.write("D_finite = change in v48 from turning the window fully OFF\n")
        fh.write("           linear response would give D_finite = -D_local, ratio = 1\n\n")
        for title, rs in (("=== 3-hour windows ===", rows),
                          ("=== 6-hour blocks ===", blocks)):
            fh.write(title + "\n")
            fh.write(f"{'sector':<6}{'hours':<10}{'D_local':>11}{'D_finite':>11}"
                     f"{'-D_local':>11}{'ratio':>9}{'sign_ok':>9}\n")
            for r in rs:
                hs = r.get("hours", r.get("win"))
                fh.write(f"{r['sector']:<6}{hs:<10}{r['local']:11.3f}{r['finite']:11.3f}"
                         f"{-r['local']:11.3f}{r['ratio']:9.2f}{str(r['sign_ok']):>9}\n")
            nbad = sum(1 for r in rs if not r["sign_ok"])
            fh.write(f"  -> sign disagreements: {nbad} of {len(rs)}\n\n")
    nbad3 = sum(1 for r in rows if not r["sign_ok"])
    nbad6 = sum(1 for r in blocks if not r["sign_ok"])
    print(f"\nsign disagreements: 3h windows {nbad3}/{len(rows)}, "
          f"6h blocks {nbad6}/{len(blocks)}")
    print("wrote", out)


if __name__ == "__main__":
    main()
