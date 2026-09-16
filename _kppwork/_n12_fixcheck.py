"""Recompute the n=12 v48 spectrum from the CLEAN co_v_m8 table and compare
with the old spectrum_orders.txt numbers (which used 18 corrupt rows)."""
import csv
import glob
import importlib.util
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("expspec", os.path.join(HERE, "_n16_spectrum.py"))
_ex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ex)
analyse = _ex.analyse

for scen in ("r_sat", "r_lim1"):
    v = {}
    for f in glob.glob(os.path.join(HERE, "co_v_m8_*.csv")):
        for row in csv.DictReader(open(f)):
            nm = row["name"]
            if nm.startswith(scen + "_m") and row["v_ppb"] not in ("", "NaN"):
                v[int(nm.rsplit("_m", 1)[1])] = float(row["v_ppb"])
    assert len(v) == 4096, (scen, len(v))
    true, trunc, mass, scale = analyse(v, 12, 12)
    print("%s n=12 CLEAN v48:  max|phi|=%.3f" % (scen, scale))
    print("   trunc: " + ", ".join("q<=%d: %.1f%%" % (q, 100 * trunc[q])
                                   for q in sorted(trunc)[:8]))
    base = mass[1] or 1.0
    print("   mass : " + " ".join("k=%d:%.3f" % (k, mass[k] / base)
                                  for k in range(1, 13)))
