#!/usr/bin/env python3
"""Scan A generator: local-vs-global protocol on a 5x5 (NOX_FACTOR, VOC_FACTOR) grid.

Prevalence layer of claim (6): the baseline instance showed 1/24 sign flips at 3 h
under BOTH metrics; this grid asks WHERE on the regime map the flips live.

Per grid point, the same protocol as _local_vs_global.py:
    base, and for each of 3 sectors x 8 3-h windows: off / +5% / -5% runs
    -> D_local = (p - m) / (2*DELTA);  D_finite = off - base;  sign flip if
    sign(D_local) != sign(-D_finite).
73 namelists x 25 points = 1825 box runs (workers: sweep_flips.sh).

Backend "scale" (kpp7) with 8-element NOXW/VOCW/COW -- the same build that ran
the n=16 enumeration, so the 3 h windows match the oracle bit conventions.

Output: /root/flipscan/<ptag>[_<SEC>_<kind><i>].nml
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
from kpp_box import nml_text  # noqa: E402

ROOT = "/root/flipscan"
GRID = [0.5, 0.75, 1.0, 1.5, 2.0]
SECTORS = ["NOXW", "VOCW", "COW"]
NW = 8
ALL1 = [1.0] * NW
DELTA = 0.05


def ptag(noxf, vocf):
    return "n%03dv%03d" % (int(round(noxf * 100)), int(round(vocf * 100)))


def main():
    os.makedirs(ROOT, exist_ok=True)
    for f in os.listdir(ROOT):
        if f.endswith(".nml"):
            os.remove(os.path.join(ROOT, f))
    n = 0

    def w(name, cfg):
        nonlocal n
        full = {"NOXW": ALL1, "VOCW": ALL1, "COW": ALL1}
        full.update(cfg)
        with open(os.path.join(ROOT, name + ".nml"), "w", newline="\n") as fh:
            fh.write(nml_text(full, "scale"))
        n += 1

    for noxf in GRID:
        for vocf in GRID:
            pt = ptag(noxf, vocf)
            scen = {"NOX_FACTOR": noxf, "VOC_FACTOR": vocf}
            w(pt + "_base", dict(scen))
            for sec in SECTORS:
                for i in range(NW):
                    off = list(ALL1); off[i] = 0.0
                    plus = list(ALL1); plus[i] = 1.0 + DELTA
                    minus = list(ALL1); minus[i] = 1.0 - DELTA
                    w("%s_%s_off%d" % (pt, sec, i), {**scen, sec: off})
                    w("%s_%s_p%d" % (pt, sec, i), {**scen, sec: plus})
                    w("%s_%s_m%d" % (pt, sec, i), {**scen, sec: minus})
    print("wrote %d namelists to %s (expect 25*73 = 1825)" % (n, ROOT))


if __name__ == "__main__":
    main()
