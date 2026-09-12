"""G4 qualitative gates for the KPP 3.5 CB05 scenario box (either backend).

Pre-registered decision rules (fixed before running; no post-hoc moving):
 G4a EKMA curvature: on a 6x6 (nox_factor x voc_factor) grid, metric MDA8,
     the finite-difference slope dMDA8/dnox must be POSITIVE at some voc
     level and NEGATIVE at another  -> NOx-limited and NOx-saturated
     regimes both exist (ridge interior to the grid counts either way).
 G4b reduce-NOx-only worsens: at voc_factor=1.5, MDA8(nox=0.7) must EXCEED
     MDA8(nox=1.0).
 G4c weekend effect: weekend scenario (WKEND=True, NOX_FACTOR=0.7,
     VOC_FACTOR=1.0) day-3 max O3 must exceed the weekday baseline
     (all factors 1.0) day-3 max O3.

PASS = G4a and G4b and G4c.

Usage: python exp_g4_kpp.py [backend]      ("dummy" default, or "split")
Writes g4_results_<backend>.txt and g4_grid_<backend>.csv (for the A/B diff).
"""
import sys
import time
from kpp_box import run_scenario

BACKEND = sys.argv[1] if len(sys.argv) > 1 else "dummy"

NOX = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]
VOC = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5]


def main():
    t0 = time.time()
    grid = {}
    for nf in NOX:
        for vf in VOC:
            name = "ekma_n%.2f_v%.2f" % (nf, vf)
            m = run_scenario({"NOX_FACTOR": nf, "VOC_FACTOR": vf}, name, BACKEND)
            grid[(nf, vf)] = m["mda8"]
            assert m["nfail"] == 0, name
            print("%s mda8=%.2f" % (name, m["mda8"]), flush=True)

    # G4a: slopes along NOx at each VOC level (adjacent pairs)
    slopes = {}
    for vf in VOC:
        slopes[vf] = [(grid[(NOX[i + 1], vf)] - grid[(NOX[i], vf)]) / (NOX[i + 1] - NOX[i])
                      for i in range(len(NOX) - 1)]
    g4a = any(s > 0 for sl in slopes.values() for s in sl) and \
          any(s < 0 for sl in slopes.values() for s in sl)

    # G4b: voc=1.5, nox 1.0 vs 0.7
    mb = run_scenario({"NOX_FACTOR": 0.7, "VOC_FACTOR": 1.5}, "g4b_n07_v15", BACKEND)
    g4b = mb["mda8"] > grid[(1.0, 1.5)]
    print("g4b: mda8(0.7,1.5)=%.2f vs base(1.0,1.5)=%.2f" %
          (mb["mda8"], grid[(1.0, 1.5)]))

    # G4c: weekend vs weekday baseline, metric day-3 max O3
    mw = run_scenario({"WKEND": True, "NOX_FACTOR": 0.7, "VOC_FACTOR": 1.0},
                      "g4c_weekend", BACKEND)
    base = run_scenario({}, "g4c_weekday", BACKEND)
    g4c = mw["d3_max"] > base["d3_max"]
    print("g4c: weekend d3max=%.2f vs weekday d3max=%.2f" %
          (mw["d3_max"], base["d3_max"]))

    L = []
    L.append("G4 results, backend=%s, wall=%.0fs" % (BACKEND, time.time() - t0))
    L.append("")
    L.append("EKMA grid: MDA8 ppb, rows=nox_factor, cols=voc_factor")
    L.append("        " + "".join("%8.1f" % v for v in VOC))
    for nf in NOX:
        L.append("%6.2f  " % nf + "".join("%8.1f" % grid[(nf, vf)] for vf in VOC))
    L.append("")
    L.append("slopes dMDA8/dNOx per VOC level:")
    for vf in VOC:
        L.append("  voc=%.2f: " % vf + " ".join("%+8.1f" % s for s in slopes[vf]))
    L.append("")
    L.append("G4a EKMA curvature (both signs present): %s" % ("PASS" if g4a else "FAIL"))
    L.append("G4b reduce-NOx-only worsens at voc=1.5 : %s" % ("PASS" if g4b else "FAIL"))
    L.append("G4c weekend positive anomaly (d3 max)  : %s" % ("PASS" if g4c else "FAIL"))
    L.append("VERDICT G4: %s" % ("PASS" if (g4a and g4b and g4c) else "FAIL"))
    txt = "\n".join(L)
    print("\n" + txt)
    with open("g4_results_%s.txt" % BACKEND, "w") as f:
        f.write(txt + "\n")
    with open("g4_grid_%s.csv" % BACKEND, "w") as f:
        f.write("# metric,nox_factor,voc_factor,value\n")
        for nf in NOX:
            for vf in VOC:
                f.write("ekma,%.2f,%.2f,%.6f\n" % (nf, vf, grid[(nf, vf)]))
        f.write("g4b,0.70,1.50,%.6f\n" % mb["mda8"])
        f.write("g4c_weekend,0.70,1.00,%.6f\n" % mw["d3_max"])
        f.write("g4c_weekday,1.00,1.00,%.6f\n" % base["d3_max"])


if __name__ == "__main__":
    main()
