#!/usr/bin/env python3
"""Step 1 of the MLMC test: is a coarser integration step actually cheaper, and
how far does its MDA8 depart from the fine level?

Level idea: same chemistry, same coalition, larger DT_S.  The box cost scales
with the number of integration steps, so DT_S = 900 should cost about a third of
DT_S = 300.  MLMC's ceiling is exactly that cost ratio, so it has to be measured
before any estimator work.

Writes results_mlmc_levels.txt
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kpp_box  # noqa: E402
import _gen16_city as g  # noqa: E402

g.D = r"D:\atmos"          # the generator defaults to the WSL mount
build_base, cfg_of = g.build_base, g.cfg_of

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "results_mlmc_levels.txt")
LEVELS = (300.0, 600.0, 900.0, 1800.0)


def main():
    L, add = [], None
    add = L.append
    add("MLMC level probe: cost and MDA8 against the integration step")
    add("")
    for label in ("dallas_2016w", "deerpark_2017s"):
        base, _ = build_base(label)
        add("=" * 76)
        add(label)
        add("=" * 76)
        add("  %7s %9s %11s %9s %8s %8s"
            % ("DT_S", "wall_s", "MDA8", "v48", "nfail", "nrow"))
        ref = None
        for dt in LEVELS:
            cfg = dict(base)
            cfg.update(cfg_of(0))          # all players ON
            cfg["DT_S"] = dt
            t0 = time.time()
            m = kpp_box.run_scenario(cfg, "mlmc_%s_%d" % (label, int(dt)), "scale")
            el = time.time() - t0
            if ref is None:
                ref = m["mda8"]
            add("  %7.0f %9.2f %11.4f %9.4f %8d %8d   dMDA8=%+.4f"
                % (dt, el, m["mda8"], m["d3_mean"], m["nfail"], len(m["traj"]),
                   m["mda8"] - ref))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
