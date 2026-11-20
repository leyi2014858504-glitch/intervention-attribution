#!/usr/bin/env python3
"""Find the BASELINE winning 8-h window, to define the frozen-window control.

The control keeps the metric as an 8-h mean over one fixed time span (the one
that wins at the all-on state) and removes only the argmax.  Comparing it with
MDA8 isolates the selector while holding the functional form fixed.

Writes _kppwork/frozen_<label>.txt with "i0 w" (index of the window's last
sample and the window length in samples).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kpp_box  # noqa: E402
import _gen16_city as g  # noqa: E402

g.D = r"D:\atmos"
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    for label in ("dallas_2016w", "deerpark_2017s"):
        base, _ = g.build_base(label)
        cfg = dict(base)
        cfg.update(g.cfg_of(0))
        m = kpp_box.run_scenario(cfg, "frozen_%s" % label, "scale")
        traj = m["traj"]
        t = [r[0] for r in traj]
        o = [r[1] for r in traj]
        dt = t[1] - t[0]
        w = max(1, int(round(8.0 / dt)))
        best, i0 = -1e9, 0
        rs = 0.0
        for i in range(len(o)):
            rs += o[i]
            if i >= w:
                rs -= o[i - w]
            if i >= w - 1 and t[i] >= 48.0:
                a = rs / w
                if a > best:
                    best, i0 = a, i
        with open(os.path.join(HERE, "frozen_%s.txt" % label), "w") as fh:
            fh.write("%d %d %.4f\n" % (i0, w, t[i0]))
        print("%-16s baseline MDA8 %.4f at index %d (t=%.2f h), w=%d samples"
              % (label, best, i0, t[i0], w))


if __name__ == "__main__":
    main()
