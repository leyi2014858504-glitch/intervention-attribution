#!/usr/bin/env python3
"""Scan A report: sign-flip prevalence across the (NOX_FACTOR, VOC_FACTOR) grid.

Reads flip_v_*.csv shards (dual metric), applies the _local_vs_global protocol
per grid point, and writes results_scan_flips.txt:
  * per-point summary (flips under MDA8 and v48, ratio deviations by sector)
  * the 5x5 flip-count grid (the prevalence table)
  * a sentinel: point n100v100 MUST reproduce the baseline result (NOx 12-15h
    flips under both metrics), else the scan drifted and is not trustworthy.
"""
import csv
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_scan_flips.txt")
DELTA = 0.05
NW = 8
GRID = [0.5, 0.75, 1.0, 1.5, 2.0]
SECTORS = ["NOXW", "VOCW", "COW"]


def load():
    vals = {}
    for f in sorted(glob.glob(os.path.join(HERE, "flip_v_*.csv"))):
        for row in csv.DictReader(open(f)):
            if row["v_ppb"] in ("", "NaN") or row["v_mda8"] in ("", "NaN"):
                continue
            vals[row["name"]] = (float(row["v_ppb"]), float(row["v_mda8"]))
    return vals


def one_point(pt, vals):
    base = vals.get(pt + "_base")
    if base is None:
        return None
    out = {"pt": pt, "noxf": int(pt[1:4]) / 100.0, "vocf": int(pt[5:8]) / 100.0,
           "base48": base[0], "base8": base[1], "flips8": [], "flips48": [],
           "dev8": {"NOXW": [], "VOCW": [], "COW": []},
           "dev48": {"NOXW": [], "VOCW": [], "COW": []}}
    for sec in SECTORS:
        for i in range(NW):
            try:
                v_off = vals["%s_%s_off%d" % (pt, sec, i)]
                v_p = vals["%s_%s_p%d" % (pt, sec, i)]
                v_m = vals["%s_%s_m%d" % (pt, sec, i)]
            except KeyError:
                continue
            for tag, j in (("8", 1), ("48", 0)):
                Dloc = (v_p[j] - v_m[j]) / (2 * DELTA)
                Dfin = v_off[j] - base[j]
                flipped = (Dloc * (-Dfin)) <= 0 and abs(Dloc) > 1e-9 and abs(Dfin) > 1e-9
                if flipped:
                    out["flips" + tag].append("%s@%d-%d" % (sec, 3 * i, 3 * i + 3))
                if abs(Dloc) > 1e-9:
                    out["dev" + tag][sec].append(abs(Dfin / (-Dloc) - 1.0))
    return out


def med(xs):
    return sorted(xs)[len(xs) // 2] if xs else float("nan")


def main():
    vals = load()
    pts = sorted({n.split("_")[0] for n in vals})
    full = [r for r in (one_point(p, vals) for p in pts) if r]
    L = []
    add = L.append
    add("Scan A: sign-flip prevalence across the (NOX_FACTOR, VOC_FACTOR) grid")
    add("protocol = _local_vs_global (3h windows, DELTA=0.05, both metrics)")
    add("flip = sign(D_local) != sign(D_finite_off) with both above noise")
    add("")
    add("%-10s %5s %5s %6s %8s %8s  %6s %6s %6s | %s" %
        ("point", "NOxF", "VOCF", "base8", "flips/24", "(v48)",
         "NOxdev", "VOCdev", "COdev", "flipped windows (MDA8)"))
    for r in sorted(full, key=lambda r: (r["noxf"], r["vocf"])):
        add("%-10s %5.2f %5.2f %8.2f %6d %6d  %6.2f %6.2f %6.2f | %s" %
            (r["pt"], r["noxf"], r["vocf"], r["base8"],
             len(r["flips8"]), len(r["flips48"]),
             max(r["dev8"]["NOXW"], default=float("nan")),
             med(r["dev8"]["VOCW"]), med(r["dev8"]["COW"]),
             ", ".join(r["flips8"]) or "-"))
    add("")
    add("5x5 flip-count grid (MDA8), rows NOX_FACTOR, cols VOC_FACTOR:")
    add("      " + "".join("%7.2f" % v for v in GRID))
    for nf in GRID:
        row = []
        for vf in GRID:
            m = [r for r in full if r["noxf"] == nf and r["vocf"] == vf]
            row.append(len(m[0]["flips8"]) if m else -1)
        add("%5.2f " % nf + "".join("%7s" % x for x in row))
    add("")
    # sentinel: baseline point must reproduce NOx 12-15h flip under both metrics
    b = [r for r in full if r["pt"] == "n100v100"]
    if b:
        ok8 = any("NOXW@12-15" in x for x in b[0]["flips8"])
        ok48 = any("NOXW@12-15" in x for x in b[0]["flips48"])
        add("SENTINEL n100v100 (must match published baseline): "
            "MDA8 flip NOx12-15 = %s, v48 flip = %s -> %s" %
            (ok8, ok48, "OK" if (ok8 and ok48) else "!! DRIFTED, DO NOT TRUST"))
    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("\nwrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
