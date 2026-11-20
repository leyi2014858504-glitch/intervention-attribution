#!/usr/bin/env python3
"""Taylor-2 (HDDM) coefficient namelists for the CITY-ANCHORED n=16 games.

Same 512-run protocol as gen_hddm.py (16 players x {1+delta, 1-delta}, plus
C(16,2)=120 pairs x 4 corners, delta = 0.05), but the base config is the
observationally anchored slice from results_fit_<label>.txt -- byte-identical
to the base used by _gen16_city.py, so the derivatives are taken on the same
model that produced the truth table.

Player layout (same as _gen16_city.py): bit w = NOx window w OFF (w < 8),
bit 8+w = VOC window w OFF.  "Off" here is fractional: multiplier 1 +- delta.

Usage: python _gen_hddm_city.py dallas_2016w deerpark_2017s
Writes /root/nmlhddm/<label>/<name>.nml   (512 namelists per label)
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
from kpp_box import nml_text  # noqa: E402
from _gen16_city import build_base  # noqa: E402

ROOT = "/root/nmlhddm"
NWIN = 8
NPLAY = 16
DELTA = 0.05


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


def main(labels):
    n = 0
    for label in labels:
        base, anchor = build_base(label)
        d = os.path.join(ROOT, label)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.endswith(".nml"):
                os.remove(os.path.join(d, f))

        def w(name, pert):
            nonlocal n
            cfg = dict(base)
            cfg.update(cfg_frac(pert))
            with open(os.path.join(d, name + ".nml"), "w", newline="\n") as fh:
                fh.write(nml_text(cfg, "co"))
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
        print("%-18s 512 namelists -> %s   anchor=%s" % (label, d, anchor))
    print("total namelists: %d (expect 2 x 512 = 1024)" % n)


if __name__ == "__main__":
    main(sys.argv[1:] or ["dallas_2016w", "deerpark_2017s"])
