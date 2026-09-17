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
           "dev48": {"NOXW": [], "VOCW": [], "COW": []},
           "abs8": {"NOXW": [], "VOCW": [], "COW": []},
           "abs48": {"NOXW": [], "VOCW": [], "COW": []},
           "maxfin8": 0.0, "maxfin48": 0.0, "ridge8": [], "ridge48": []}
    # first pass: collect everything, track the largest finite effect per metric
    tmp = []
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
                tmp.append((sec, i, tag, Dloc, Dfin))
                out["maxfin" + tag] = max(out["maxfin" + tag], abs(Dfin))
    # second pass: deviations with the sane normalisation
    for sec, i, tag, Dloc, Dfin in tmp:
        scale = out["maxfin" + tag] or 1.0
        flipped = (Dloc * (-Dfin)) <= 0 and abs(Dloc) > 1e-9 and abs(Dfin) > 1e-9
        if flipped:
            out["flips" + tag].append("%s@%d-%d" % (sec, 3 * i, 3 * i + 3))
        if abs(Dloc) > 1e-9:
            out["dev" + tag][sec].append(abs(Dfin / (-Dloc) - 1.0))
        # absolute error of the DDM-1 prediction, normalised by max|effect|:
        # unlike |ratio-1| this does NOT blow up on the ridge (D1 -> 0)
        out["abs" + tag][sec].append(abs(Dfin - (-Dloc)) / scale)
        # ridge marker: the real "zero-effect" failure -- the tangent is at
        # least 10x smaller than the true finite effect, i.e. the local method
        # reports "no effect" for a lever that matters.  (An earlier, looser
        # test -- small |D1| relative to the point's max effect -- also flagged
        # small-but-consistent windows such as night-time VOC; those are not
        # ridges, they are simply weak players.)
        if abs(Dfin) > 0.5 and abs(Dloc) < 0.1 * abs(Dfin):
            out["ridge" + tag].append("%s@%d-%d" % (sec, 3 * i, 3 * i + 3))
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
    add("NOxdev/VOCdev/COdev = |ratio-1|, RATIO-BASED (blows up on ridge, D1~0)")
    add("NOxabs/VOCabs/COabs = |D_finite - (-D_local)| / max|D_finite|  (RIDGE-SAFE)")
    add("")
    add("%-10s %5s %5s %6s %8s %8s  %6s %6s %6s  %6s %6s %6s | %s" %
        ("point", "NOxF", "VOCF", "base8", "flips/24", "(v48)",
         "NOxabs", "VOCabs", "COabs", "NOxdev", "VOCdev", "COdev",
         "flipped windows (MDA8)"))
    for r in sorted(full, key=lambda r: (r["noxf"], r["vocf"])):
        flips = ", ".join(r["flips8"]) or "-"
        ridge = ", ".join(r["ridge8"]) or ""
        line = ("%-10s %5.2f %5.2f %8.2f %6d %6d  %6.2f %6.2f %6.2f  %6.2f %6.2f %6.2f | %s" %
                (r["pt"], r["noxf"], r["vocf"], r["base8"],
                 len(r["flips8"]), len(r["flips48"]),
                 max(r["abs8"]["NOXW"], default=float("nan")),
                 med(r["abs8"]["VOCW"]), med(r["abs8"]["COW"]),
                 max(r["dev8"]["NOXW"], default=float("nan")),
                 med(r["dev8"]["VOCW"]), med(r["dev8"]["COW"]),
                 flips))
        add(line)
        if ridge:
            add("%10s Zero-effect (tangent >=10x smaller than truth): %s" % ("", ridge))
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
