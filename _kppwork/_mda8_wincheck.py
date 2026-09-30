#!/usr/bin/env python3
"""Do MDA8 maxima ever land in windows straddling the 48 h equilibration cut?

All campaigns compute MDA8 as the max over 8 h rolling windows ENDING at
t >= 48 h (10-min output).  A window ending at t < 56 h reaches back into
[40, 48) -- the tail of the equilibration period.  If, on spot-check
configurations spanning both diurnal shapes, the maximum always coincides
with the fully-contained variant (windows ending at t >= 56 h), the paper
can adopt the contained definition WITHOUT re-running anything.

Configs: Deer Park anchored slice (base, NOx/VOC 12-15 off, all-NOx off)
plus a high-NOx variant to cover the saturated-side diurnal shape.
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402

D = "/mnt/d/atmos"
T = []
with open(os.path.join(D, "met_deerpark_2017s.csv")) as fh:
    for r in fh:
        p = r.strip().split(",")
        if p[0] == "time":
            continue
        try:
            T.append(float(p[1]))
        except (ValueError, IndexError):
            pass
TMIN, TMAX = min(T) + 273.15, max(T) + 273.15

BASE = dict(EMIS_NOX=0.5, EMIS_VOC=4.2406, EMIS_ISO=10.0, EMIS_CO=1.0,
            KDIL_DAY=2233.674, KDIL_NIGHT=86400.0,
            BG_O3=5.0, BG_CO=136.875,
            BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
            TMIN_K=TMIN, TMAX_K=TMAX)

CFGS = [
    ("base", {}),
    ("nox12_15off", {"NOXW": [1, 1, 1, 1, 0, 1, 1, 1]}),
    ("voc12_15off", {"VOCW": [1, 1, 1, 1, 0, 1, 1, 1]}),
    ("allnoxoff", {"NOXW": [0] * 8}),
    ("highnox_base", {"EMIS_NOX": 5.0}),
]


def mda8_variants(o3, tt):
    dt = tt[1] - tt[0]
    w = int(round(8.0 / dt))
    rs = 0.0
    bs, bc = -1.0, -1.0
    ts, tc = None, None
    for i in range(len(o3)):
        rs += o3[i]
        if i >= w:
            rs -= o3[i - w]
        if i >= w - 1:
            a = rs / w
            if tt[i] >= 48.0 and a > bs:
                bs, ts = a, tt[i]
            if tt[i] >= 56.0 and a > bc:
                bc, tc = a, tt[i]
    return bs, ts, bc, tc, w


print("%-14s %4s %10s %7s | %10s %7s | %s" %
      ("cfg", "w", "mda8>=48", "t_max", "mda8>=56", "t_max", "verdict"))
for name, pert in CFGS:
    cfg = dict(BASE)
    cfg.update(pert)
    m = kpp_box.run_scenario(cfg, "winchk_" + name, "scale")
    traj = m["traj"]
    tt = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    bs, ts, bc, tc, w = mda8_variants(o3, tt)
    same = abs(bs - bc) < 1e-9
    print("%-14s %4d %10.4f %7.2f | %10.4f %7.2f | %s" %
          (name, w, bs, ts, bc, tc,
           "IDENTICAL" if same else "DIFFERS %+.4f" % (bs - bc)))
