"""Is the Dallas '27.7x' robust to the finite-difference step?

The derivative at NOx 06-09 h is small and sits next to an argmax switch, so it
may be step-dependent.  Recompute D1 at delta = 0.02 / 0.05 / 0.10 for that
player (and for NOx 12-15 h as a control) with the production build, using the
same hourly-MDA8 convention as the tables.
"""
import os
import sys
import numpy as np

sys.path.insert(0, r"D:\atmos\_kppwork")
os.environ["ATMOS_DIR"] = r"D:\atmos"
import kpp_box                                      # noqa: E402
from _gen16_city import build_base                  # noqa: E402

NW = 8


def hourly_mda8(traj):
    t = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    hm = {}
    for i in range(len(t)):
        hm.setdefault(int(t[i]), []).append(o3[i])
    hm = {h: sum(v) / len(v) for h, v in hm.items() if h >= 40}
    best = -1e9
    for e in range(48, 72):
        vals = [hm.get(e - 7 + j) for j in range(8)]
        if any(v is None for v in vals):
            continue
        best = max(best, sum(vals) / 8)
    return best


def cfg_frac(base, pert):
    nox = [1.0] * NW
    voc = [1.0] * NW
    for i, m in pert.items():
        if i < NW:
            nox[i] = m
        else:
            voc[i - NW] = m
    c = dict(base)
    c.update({"NOXW": nox, "VOCW": voc, "COW": [1.0] * NW})
    return c


base, _ = build_base("dallas_2016w")
v0 = hourly_mda8(kpp_box.run_scenario(cfg_frac(base, {}), "step_base", "scale")["traj"])
print("base hourly MDA8 = %.4f" % v0)

EXACT = {2: None, 4: None}
z = np.load(r"D:\atmos\_kppwork\aw_h_dallas_2016w.npz")
mh = z["mda8_h"]
for i in (2, 4):
    EXACT[i] = mh[1 << i] - mh[0]

for i in (2, 4):
    print("=" * 66)
    print("player %d (NOx %s), exact single-shutdown effect %+.4f ppb"
          % (i, ("06-09" if i == 2 else "12-15"), EXACT[i]))
    for delta in (0.02, 0.05, 0.10, 0.20):
        vp = hourly_mda8(kpp_box.run_scenario(
            cfg_frac(base, {i: 1.0 + delta}), "step_p%d_%02d" % (i, 100 * delta),
            "scale")["traj"])
        vm = hourly_mda8(kpp_box.run_scenario(
            cfg_frac(base, {i: 1.0 - delta}), "step_m%d_%02d" % (i, 100 * delta),
            "scale")["traj"])
        d1 = (vp - vm) / (2 * delta)
        print("   delta=%.2f  D1 %+.5f  -D1 %+.5f  exact/(-D1) = %s"
              % (delta, d1, -d1,
                 ("%.1f" % (EXACT[i] / -d1)) if abs(d1) > 1e-9 else "inf"))
