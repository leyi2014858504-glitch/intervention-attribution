#!/usr/bin/env python3
"""Is the near-zero tangent at the ridge points physical or truncation error?

WHY THIS REPLACES THE EARLIER "219x the floor" CLAIM
  _zero_effect.py compared |D1| against delta^2 * |D2|.  That comparison is
  DIMENSIONALLY INCONSISTENT: delta^2*|D2| has the units of the response (ppb),
  while D1 has units of ppb per unit multiplier, so their ratio is not a
  signal-to-noise number.  (The arithmetic did run: 0.0885/4.04e-04 = 219, but
  the quantity is meaningless.)  A third derivative would give the right units
  -- the central-difference truncation error is (delta^2/6) f''' -- but f'''
  is NOT measurable from runs at 1+-delta alone: a third-order stencil needs a
  wider point set.

  So we settle it empirically instead, with the test that needs no f''':
  DELTA REFINEMENT.  For a central difference,

      D1(delta) = f'(1) + (delta^2/6) f'''(1) + O(delta^4)

  so if D1 barely moves when delta is halved, the truncation term is negligible
  and a flat tangent is a property of the response surface, not of the stencil.
  Three deltas (0.20, 0.10, 0.05) also allow a Richardson extrapolation to
  f'(1) itself.

Writes results_ridge_delta.txt
"""
import glob
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_ridge_delta.txt")
NWIN = 8
DELTAS = (0.20, 0.10)          # 0.05 already exists in the flipscan tables
CASES = [
    # (tag, NOX_FACTOR, VOC_FACTOR, window index, why)
    ("n075v100", 0.75, 1.00, 4, "zero-effect (12-15)"),
    ("n075v100", 0.75, 1.00, 0, "control (0-3)"),
    ("n050v100", 0.50, 1.00, 3, "zero-effect (9-12)"),
    ("n050v100", 0.50, 1.00, 0, "control (0-3)"),
]


def run(noxf, vocf, i, mult, name):
    nox = [1.0] * NWIN
    nox[i] = mult
    cfg = {"NOX_FACTOR": noxf, "VOC_FACTOR": vocf,
           "NOXW": nox, "VOCW": [1.0] * NWIN, "COW": [1.0] * NWIN}
    m = kpp_box.run_scenario(cfg, "ridge_" + name, "scale")
    return _metrics(m["traj"])


def _metrics(traj):
    tt = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    sel = [i for i in range(len(tt)) if tt[i] >= 48.0]
    v48 = sum(o3[i] for i in sel) / len(sel)
    dt = tt[1] - tt[0]
    w = max(1, int(round(8.0 / dt)))
    best = float("-inf")
    for i in range(w - 1, len(o3)):
        if tt[i] < 48.0:
            continue
        a = sum(o3[i - w + 1:i + 1]) / w
        if a > best:
            best = a
    return v48, best


def main():
    L = []
    add = L.append
    add("Delta-refinement test at the zero-effect windows and a control")
    add("D1(delta) = [f(1+delta) - f(1-delta)] / (2 delta);  f'(1) = D1 + O(delta^2)")
    add("")

    for tag, noxf, vocf, i, why in CASES:
        add("%s  NOX_FACTOR=%.2f VOC_FACTOR=%.2f  window %d-%d   [%s]" %
            (tag, noxf, vocf, 3 * i, 3 * i + 3, why))
        for j, mname in ((1, "MDA8"), (0, "v48")):
            vals = {}
            for d in DELTAS:
                vp = run(noxf, vocf, i, 1.0 + d, "%s_w%d_p%03d" % (tag, i, int(d * 1000)))
                vm = run(noxf, vocf, i, 1.0 - d, "%s_w%d_m%03d" % (tag, i, int(d * 1000)))
                vals[d] = (vp[j] - vm[j]) / (2 * d)
            # delta = 0.05 from the existing flipscan table
            if not hasattr(main, "_flip"):
                d = {}
                for f in glob.glob(os.path.join(HERE, "flip_v_*.csv")):
                    for ln in open(f):
                        p = ln.strip().split(",")
                        if len(p) < 3 or p[1] in ("NaN", "v_ppb"):
                            continue
                        try:
                            d[p[0]] = (float(p[1]), float(p[2]))
                        except ValueError:
                            pass
                main._flip = d
            d005 = main._flip
            kp = "%s_NOXW_p%d" % (tag, i)
            km = "%s_NOXW_m%d" % (tag, i)
            v005 = (d005[kp][j] - d005[km][j]) / (2 * 0.05) if (kp in d005 and km in d005) else None
            dens = {}
            for d in (0.05, 0.10, 0.20):
                dens[d] = v005 if d == 0.05 else vals[d]
            add("    %-5s D1: delta=0.05 %+8.4f | 0.10 %+8.4f | 0.20 %+8.4f"
                % (mname, dens[0.05], dens[0.10], dens[0.20]))
            if v005 is not None:
                # Richardson: f' = D1(d) - c d^2, fit on the two smallest deltas
                c = (dens[0.05] - dens[0.10]) / (0.05 ** 2 - 0.10 ** 2)
                fprime = dens[0.05] - c * 0.05 ** 2
                chg = abs(dens[0.10] - dens[0.05])
                add("             |D1(0.10)-D1(0.05)| = %.4f  -> extrapolated f'(1) = %+.4f"
                    % (chg, fprime))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("wrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
