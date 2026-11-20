"""Close the DT=150 s endpoint of the '150-600 s, <0.06%' claim.

The archived MLMC probe only covers 300/600/900/1800 s.  This reruns the
all-on state at DT=150 s with the production build (kpp7) and recomputes the
hourly MDA8 exactly as the tables do (hourly means h=40..71, then 8-h rolling
windows ending 48..71).
"""
import os
import sys
import numpy as np

sys.path.insert(0, r"D:\atmos\_kppwork")
os.environ["ATMOS_DIR"] = r"D:\atmos"
import kpp_box                                      # noqa: E402
from _gen16_city import build_base                  # noqa: E402


def hourly_mda8(traj):
    t = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    dt = t[1] - t[0]
    per_h = max(1, int(round(1.0 / dt)))
    hm = {}
    for i in range(len(t)):
        h = int(t[i])          # floor: the convention used by the sweep
        hm.setdefault(h, []).append(o3[i])
    hm = {h: sum(v) / len(v) for h, v in hm.items() if h >= 40}
    best = -1e9
    for e in range(48, 72):
        vals = [hm.get(e - 7 + j) for j in range(8)]
        if any(v is None for v in vals):
            continue
        best = max(best, sum(vals) / 8)
    return best


L = []
for lab in ("dallas_2016w", "deerpark_2017s"):
    base, _ = build_base(lab)
    z = np.load(r"D:\atmos\_kppwork\aw_h_%s.npz" % lab)
    ref = float(z["mda8_h"][0])
    print("=" * 60)
    print("%-14s production (300 s) all-on hourly MDA8 = %.4f" % (lab, ref))
    L.append("%-14s production (300 s) all-on hourly MDA8 = %.4f" % (lab, ref))
    for dt in (150, 300, 600):
        m = kpp_box.run_scenario(dict(base, DT_S=float(dt)),
                                 "dt%d_%s" % (dt, lab), "scale")
        val = hourly_mda8(m["traj"])
        print("   DT=%4d s  MDA8 %.4f  dMDA8 %+.4f ppb (%+.3f%%)"
              % (dt, val, val - ref, 100 * (val - ref) / ref))
        L.append("   DT=%4d s  MDA8 %.4f  dMDA8 %+.4f ppb (%+.3f%%)"
                 % (dt, val, val - ref, 100 * (val - ref) / ref))

with open(r"D:\atmos\results_dt_check.txt", "w") as fh:
    fh.write("Outer-step check of the all-on hourly MDA8 (production build kpp7,\n"
             "hourly means h=40..71 then 8-h rolling windows ending 48..71).\n\n")
    fh.write("\n".join(L) + "\n")
print("wrote results_dt_check.txt")
