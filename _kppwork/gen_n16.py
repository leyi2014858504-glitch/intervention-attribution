#!/usr/bin/env python3
"""
Generate all 2**16 namelists for the n=16 game: 2 sectors x 8 three-hour windows.

WHY n=16
--------
It is the largest game in this family that can be FULLY ENUMERATED
(2**16 = 65,536 box runs).  That yields a true oracle, which is exactly what the
n=24 point lacks: there the error had to be MAPPED from the n=12 spread/error
ratio, and that mapping is precisely the extrapolation the heavy-tail objection
attacks.  With a measured oracle at n=16 the scaling curve stops being an
extrapolation.

WHY THIS LADDER
---------------
Players are bits 0-7 = NOxW(8) and 8-15 = VOCW(8), with COW held at 1.0.  So the
sector count stays at 2 and the ONLY difference from n=8 (2 sectors x 4 windows,
whose oracle already exists in scen_v.csv) is the twofold refinement of the time
axis.  Comparing n=8 -> n=16 therefore isolates time refinement, while
n=16 -> n=24 (3 sectors x 8 windows) isolates the addition of a genuinely new
sector.

Bit convention matches gen_scale_masks.py: a SET bit means that player is
switched OFF, so mask 0 is the grand coalition (everything on).

Output: /root/nml16/<scen>/m<05d>.nml     plus /root/nml16/vgate.nml
"""
import os
import sys

from kpp_box import nml_text

NWIN = 8
N16 = 16
ROOT = "/root/nml16"
SCEN = {"r_sat": {}, "r_lim1": {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5}}


def cfg_of(mask):
    """mask -> RUNCFG overrides.  Set bit = player OFF."""
    nox = [0.0 if (mask >> w) & 1 else 1.0 for w in range(NWIN)]
    voc = [0.0 if (mask >> (NWIN + w)) & 1 else 1.0 for w in range(NWIN)]
    return {"NOXW": nox, "VOCW": voc, "COW": [1.0] * NWIN}


def main():
    total = 0
    for s, base in SCEN.items():
        d = os.path.join(ROOT, s)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.endswith(".nml"):
                os.remove(os.path.join(d, f))
        for mask in range(1 << N16):
            cfg = dict(base)
            cfg.update(cfg_of(mask))
            with open(os.path.join(d, "m%05d.nml" % mask), "w", newline="\n") as fh:
                fh.write(nml_text(cfg, "co"))
        print("%-8s %d namelists -> %s" % (s, 1 << N16, d))
        total += 1 << N16

    # smoke test target: all players ON must reproduce the frozen V-gate value
    cfg = dict(SCEN["r_sat"])
    cfg.update(cfg_of(0))
    with open(os.path.join(ROOT, "vgate.nml"), "w", newline="\n") as fh:
        fh.write(nml_text(cfg, "co"))
    print("total namelists: %d" % total)
    print("v-gate namelist written; expect v48 = 21.182293")
    return 0


if __name__ == "__main__":
    sys.exit(main())
