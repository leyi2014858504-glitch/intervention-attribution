"""Does a frozen-window DDM and a max-of-windows derivative agree?

Our T1 differentiates the max over the 24 rolling windows; a DDM-style
implementation differentiates a single frozen window (the baseline argmax).
Compare the two at the players that carry the manuscript's headline numbers.
"""
import os
import sys
import numpy as np

sys.path.insert(0, r"D:\atmos\_kppwork")
os.environ["ATMOS_DIR"] = r"D:\atmos"
import kpp_box                                      # noqa: E402
from _gen16_city import build_base                  # noqa: E402

NW = 8
DELTA = 0.05


def aw24(traj):
    t = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    hm = {}
    for i in range(len(t)):
        hm.setdefault(int(t[i]), []).append(o3[i])
    hm = {h: sum(v) / len(v) for h, v in hm.items() if h >= 40}
    A = []
    for e in range(48, 72):
        A.append(sum(hm[e - 7 + j] for j in range(8)) / 8)
    return np.array(A)


def cfg(base, pert):
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


L = []
for lab in ("dallas_2016w", "deerpark_2017s"):
    base, _ = build_base(lab)
    A0 = aw24(kpp_box.run_scenario(cfg(base, {}), "froz_base_" + lab, "scale")["traj"])
    jstar = int(A0.argmax())
    print("=" * 70)
    print("%-14s baseline max at window ending %d h (A_%d = %.4f, runner-up %.4f)"
          % (lab, jstar, jstar, A0[jstar], np.sort(A0)[-2]))
    L.append("%-14s baseline max at window ending %d h (A_%d = %.4f, runner-up %.4f)"
             % (lab, jstar, jstar, A0[jstar], np.sort(A0)[-2]))
    nsw = 0
    worst = 0.0
    for i in range(16):
        Ap = aw24(kpp_box.run_scenario(cfg(base, {i: 1 + DELTA}),
                                       "froz_p%d_%s" % (i, lab), "scale")["traj"])
        Am = aw24(kpp_box.run_scenario(cfg(base, {i: 1 - DELTA}),
                                       "froz_m%d_%s" % (i, lab), "scale")["traj"])
        d_max = (Ap.max() - Am.max()) / (2 * DELTA)
        d_froz = (Ap[jstar] - Am[jstar]) / (2 * DELTA)
        sw = (int(Ap.argmax()) != jstar) or (int(Am.argmax()) != jstar)
        nsw += sw
        if abs(d_max) > 1e-9:
            worst = max(worst, abs(d_froz / d_max - 1))
        print("  p%02d (%s %s): d(max) %+.5f | d(A_%d) %+.5f | ratio %.4f%s"
              % (i, "NOx" if i < 8 else "VOC",
                 ("%d-%d" % (3 * (i % 8), 3 * (i % 8) + 3)), d_max, jstar,
                 d_froz, d_froz / d_max if abs(d_max) > 1e-9 else float("nan"),
                 "  ARGMAX SWITCH" if sw else ""))
        L.append("  p%02d %s %s: d(max) %+.5f  d(A_%d) %+.5f  ratio %.4f%s"
                 % (i, "NOx" if i < 8 else "VOC",
                    ("%d-%d" % (3 * (i % 8), 3 * (i % 8) + 3)), d_max, jstar,
                    d_froz, d_froz / d_max if abs(d_max) > 1e-9 else float("nan"),
                    "  ARGMAX SWITCH" if sw else ""))
    print("  -> %d/16 players switch the winning window at +/-%.0f%%; "
          "worst |frozen/max - 1| = %.4f"
          % (nsw, 100 * DELTA, worst))
    L.append("  -> %d/16 players switch at +/-%.0f%%; worst |frozen/max-1| = %.4f\n"
             % (nsw, 100 * DELTA, worst))

with open(r"D:\atmos\results_frozen_window_check.txt", "w") as fh:
    fh.write("Frozen-window (DDM-style) derivative vs the derivative of the max,\n"
             "central difference delta = 0.05, production build kpp7, hourly MDA8.\n\n")
    fh.write("\n".join(L) + "\n")
print("wrote results_frozen_window_check.txt")
