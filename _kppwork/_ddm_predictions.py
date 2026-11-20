#!/usr/bin/env python3
"""Per-window DDM-1 / HDDM-2 predictions from the 512 fractional runs.

Prints, for each city and sector, the tangent (-D1), the second-order
prediction (-D1 + D2_ii/2) and the ratio to the exact finite effect read from
results_<city>_lvg.txt.
"""
import csv
import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DELTA = 0.05
SECS = ("NOXW", "VOCW")
WIN = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")


def load_frac(city):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, "hddm_h_v_%s_*.csv" % city))):
        for row in csv.DictReader(open(f)):
            if row["v_mda8"] in ("", "NaN"):
                continue
            d[row["name"]] = float(row["v_mda8"])
    return d


def parse_lvg(path):
    rows, base = {}, None
    with open(path) as fh:
        for line in fh:
            if "base v48" in line and "mda8" in line:
                base = float(line.split("mda8 =")[1].split()[0])
                continue
            p = line.split()
            if len(p) >= 6 and p[0] in ("NOXW", "VOCW", "COW"):
                rows[(p[0], p[1])] = float(p[3])
    return rows, base


for city, lab in (("dallas_2016w", "dallas"), ("deerpark_2017s", "deerpark")):
    frac = load_frac(city)
    fin, base = parse_lvg(os.path.join(ROOT, "results_%s_lvg_h.txt" % lab))
    print("=" * 74)
    print("%s   base MDA8 %.4f" % (city, base))
    for sec in SECS:
        print("  %-5s %-6s %9s %9s %9s %8s %8s" %
              ("sec", "win", "DDM-1", "HDDM-2", "finite", "fin/DDM", "fin/HD2"))
        for w, wl in enumerate(WIN):
            idx = w if sec == "NOXW" else 8 + w
            vp, vm = frac["p%02d" % idx], frac["m%02d" % idx]
            d1 = (vp - vm) / (2 * DELTA)
            d2 = (vp - 2 * base + vm) / DELTA ** 2
            p1, p2 = -d1, -d1 + 0.5 * d2
            f = fin[(sec, wl)]
            r1 = f / p1 if abs(p1) > 1e-9 else float("nan")
            r2 = f / p2 if abs(p2) > 1e-9 else float("nan")
            print("  %-5s %-6s %9.3f %9.3f %9.3f %8.2f %8.2f" %
                  (sec, wl, p1, p2, f, r1, r2))
