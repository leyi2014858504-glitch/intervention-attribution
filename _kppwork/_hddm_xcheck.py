#!/usr/bin/env python3
"""Cross-validate the hddm single-axis runs against the lvg runs.

WHY
  The first hddm sweep overlapped two scenarios on the same build dirs and
  wrote plausible-but-wrong values.  A NaN filter cannot catch that, but this
  can: the hddm run "p<i>"/"m<i>" for NOx window i uses EXACTLY the same
  namelist as lvg's "lvg_NOXW_p<i>"/"lvg_NOXW_m<i>" (same backend, same
  multiplier, all other windows at 1.0).  The trajectories therefore have to
  agree bit-for-bit.  Any disagreement means the table is contaminated.

Also checks value sanity: a +-5% change in one 3-hour window must move the
metric by a small amount, so |v - base| should be far below the metric itself
for the vast majority of runs.

usage: python3 _hddm_xcheck.py
"""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "g4runs", "scale")


def metrics(path):
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


def main():
    base = metrics(os.path.join(RUNS, "lvg_base.traj"))
    print("lvg base v48 = %.6f  mda8 = %.6f" % base)

    hddm = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_v_r_sat_*.csv"))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN"):
                continue
            hddm[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))

    print("hddm r_sat rows: %d" % len(hddm))
    SECT = {"p": "p", "m": "m"}
    worst = 0.0
    checked = 0
    for i in range(8):          # windows 0..7 = first 8 players = NOXW
        for k in ("p", "m"):
            tag = "%s%02d" % (k, i)
            if tag not in hddm:
                continue
            try:
                ref = metrics(os.path.join(RUNS, "lvg_NOXW_%s%d.traj" % (k, i)))
            except FileNotFoundError:
                continue
            d = abs(hddm[tag][0] - ref[0])
            worst = max(worst, d)
            checked += 1
    print("single-axis cross-check vs lvg: %d pairs, max|dv48| = %.6e" % (checked, worst))
    print("VERDICT:", "OK -- table consistent" if worst < 1e-6 else
          "!! CONTAMINATED -- rerun before use")

    # sanity: how far do +-5% runs move the metric?
    if hddm:
        devs = [abs(v[0] - base[0]) for v in hddm.values()]
        devs.sort()
        print("|v48 - base| over all rows: median %.3f  p90 %.3f  max %.3f"
              % (devs[len(devs) // 2], devs[int(0.9 * len(devs))], devs[-1]))
        lo = min(v[0] for v in hddm.values())
        hi = max(v[0] for v in hddm.values())
        print("v48 range across rows: %.3f .. %.3f   (base %.3f)" % (lo, hi, base[0]))


if __name__ == "__main__":
    main()
