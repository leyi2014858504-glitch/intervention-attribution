#!/usr/bin/env python3
"""Taylor-2 test at the DEER PARK anchored slice: coefficients + truth sample.

Extends the Table 7d experiment (r_sat/r_lim1 synthetic scenarios) to the real
anchored Deer Park 2017s slice.  The map predicts that on the LIMITED side
even DDM-1 should be nearly unbiased for combinations, and Taylor-2 with it;
this is the out-of-sample check of that prediction.

Generates into /root/nmlhddm/dp/:
  base, 16 players x (p,m), C(16,2) x 4 corners   = 512 derivative runs
  150 stratified random coalitions (t-prefixed)   = truth sample
    small k<=4: 40, mid 5<=k<=9: 40, large k>=10: 70
All runs carry the ANCHORED Deer Park cfg (results_fit_deerpark_2017s.txt +
the cost-function backgrounds BG_NO2=2.0/BG_NO=0.5, TMIN/TMAX from met).

Truth sample is STRATIFIED, not uniform: the r_sat failure tail lived at
large k, so large coalitions are oversampled (70 of 150) to make the
"limited side stays clean even at large k" prediction a sharp test.
"""
import os
import random
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
from kpp_box import nml_text  # noqa: E402

ROOT = "/root/nmlhddm/dp"
NWIN = 8
NPLAY = 16
DELTA = 0.05
TMIN, TMAX = 293.0, 304.4          # met_deerpark_2017s (from _fit_anchor read_met)

ANCHOR = dict(EMIS_NOX=0.5, EMIS_VOC=4.2406, EMIS_ISO=10.0, EMIS_CO=1.0,
              KDIL_DAY=2233.674, KDIL_NIGHT=86400.0,
              BG_O3=5.0, BG_CO=136.875,
              BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
              TMIN_K=TMIN, TMAX_K=TMAX)


def cfg_of(pert):
    """pert: {player: multiplier}; players 0-7 NOx 3h windows, 8-15 VOC."""
    nox = [1.0] * NWIN
    voc = [1.0] * NWIN
    for i, mult in pert.items():
        if i < NWIN:
            nox[i] = mult
        else:
            voc[i - NWIN] = mult
    return {**ANCHOR, "NOXW": nox, "VOCW": voc, "COW": [1.0] * NWIN}


def main():
    os.makedirs(ROOT, exist_ok=True)
    for f in os.listdir(ROOT):
        if f.endswith(".nml"):
            os.remove(os.path.join(ROOT, f))
    n = 0

    def w(name, pert):
        nonlocal n
        with open(os.path.join(ROOT, name + ".nml"), "w", newline="\n") as fh:
            fh.write(nml_text(cfg_of(pert), "scale"))
        n += 1

    w("base", {})
    for i in range(NPLAY):
        w("p%02d" % i, {i: 1.0 + DELTA})
        w("m%02d" % i, {i: 1.0 - DELTA})
    for i in range(NPLAY):
        for j in range(i + 1, NPLAY):
            for si, sj, tag in ((1, 1, "pp"), (1, -1, "pm"),
                                (-1, 1, "mp"), (-1, -1, "mm")):
                w("c%02d_%02d_%s" % (i, j, tag),
                  {i: 1.0 + si * DELTA, j: 1.0 + sj * DELTA})
    nderiv = n

    # stratified truth sample
    rng = random.Random(20260917)
    pools = {"s": [], "m": [], "l": []}
    for m in range(1, 1 << NPLAY):
        k = bin(m).count("1")
        if k <= 4:
            pools["s"].append(m)
        elif k <= 9:
            pools["m"].append(m)
        else:
            pools["l"].append(m)
    for grp, cnt in (("s", 40), ("m", 40), ("l", 70)):
        for m in rng.sample(pools[grp], cnt):
            pert = {i: 0.0 for i in range(NPLAY) if (m >> i) & 1}
            w("t%05d" % m, pert)
    print("wrote %d namelists to %s (deriv %d + truth %d)"
          % (n, ROOT, nderiv, n - nderiv))


if __name__ == "__main__":
    main()
