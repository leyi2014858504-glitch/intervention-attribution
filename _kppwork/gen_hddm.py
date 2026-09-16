#!/usr/bin/env python3
"""Taylor-2 (HDDM) coefficient runs on the n=16 game.

WHY
  The reviewer's strongest counter to "local methods fail" is "use HDDM, not
  DDM".  We already showed that second order REPAIRS the one-axis sign flip
  (1/36 -> 0/36), because on a single axis the entire departure from linearity
  is the on-axis curvature D2_ii.  The policy question is not one-axis though:
  it is a COMBINATION of switch-offs.  Second order then predicts

      Delta(S) ~ -sum_{i in S} D1_i + 1/2 sum_{i,j in S} D2_ij

  with n + n(n+1)/2 coefficients.  n=16 -> 152 coefficients, to be compared
  against all 2^16 - 1 = 65535 coalitions.

  KEY COST REDUCTION: the coefficients are few, the TRUTH TABLE already exists
  (n16_v_*.csv, exact).  So we need no sampling of combinations and no rerun of
  the lattice -- only the derivative runs below, after which Taylor-2 can be
  scored on EVERY coalition.

RUNS PER SCENARIO (delta = 0.05, same as _local_vs_global.py)
  base (already in the truth table, mask 0)
  16 players x (1+delta, 1-delta)                =  32
  C(16,2) = 120 pairs x 4 corners                = 480
  total new                                      = 512
  x 2 scenarios                                  = 1024

Player bit layout matches gen_n16.py: bit w = NOX window w OFF (w<8),
bit 8+w = VOC window w OFF.  Here "off" is fractional: multiplier 1 +- delta.

Output: /root/nmlhddm/<scen>/<name>.nml
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
from kpp_box import nml_text  # noqa: E402

ROOT = "/root/nmlhddm"
NWIN = 8
NPLAY = 16
DELTA = 0.05
SCEN = {"r_sat": {}, "r_lim1": {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5}}


def cfg_frac(pert):
    """pert: {player_index: multiplier} -> NOXW/VOCW lists (others stay 1.0)."""
    nox = [1.0] * NWIN
    voc = [1.0] * NWIN
    for i, mult in pert.items():
        if i < NWIN:
            nox[i] = mult
        else:
            voc[i - NWIN] = mult
    return {"NOXW": nox, "VOCW": voc, "COW": [1.0] * NWIN}


def main():
    n = 0
    for scen, base_cfg in SCEN.items():
        d = os.path.join(ROOT, scen)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.endswith(".nml"):
                os.remove(os.path.join(d, f))

        def w(name, pert):
            nonlocal n
            cfg = dict(base_cfg)
            cfg.update(cfg_frac(pert))
            with open(os.path.join(d, name + ".nml"), "w", newline="\n") as fh:
                fh.write(nml_text(cfg, "scale"))
            n += 1

        for i in range(NPLAY):
            w("p%02d" % i, {i: 1.0 + DELTA})
            w("m%02d" % i, {i: 1.0 - DELTA})
        for i in range(NPLAY):
            for j in range(i + 1, NPLAY):
                for si, sj, tag in ((1, 1, "pp"), (1, -1, "pm"),
                                    (-1, 1, "mp"), (-1, -1, "mm")):
                    w("c%02d_%02d_%s" % (i, j, tag),
                      {i: 1.0 + si * DELTA, j: 1.0 + sj * DELTA})
        print("%-8s -> %s" % (scen, d))
    print("total namelists: %d (expect 2 x 512 = 1024)" % n)


if __name__ == "__main__":
    main()
