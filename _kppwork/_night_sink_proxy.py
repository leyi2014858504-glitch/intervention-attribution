#!/usr/bin/env python3
"""Cheap hedge for the omitted nocturnal dry deposition.

Deposition is negligible by day (<=8% of the sink) but comparable to the 24 h
dilution at night (43-56% for O3).  Halving KDIL_NIGHT doubles the nocturnal
first-order sink, which brackets what O3 deposition would add.  This script
re-runs the all-on state and the eight NOx single-window finite effects at that
perturbed sink and compares them with the published anchored values.

No refit: the other knobs stay at their fitted values, so the shift measured
here is an UPPER bound (re-optimisation would absorb part of it).

Usage: python3 _night_sink_proxy.py
Writes results_night_sink_proxy.txt
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402
from _gen16_city import build_base  # noqa: E402

D = "/mnt/d/atmos"
LABELS = ("dallas_2016w", "deerpark_2017s")
NW = 8
FIT_SUF = os.environ.get("FIT_SUF", "")     # e.g. "_night2" (pinned-sink refit)
KDIL_OVERRIDE = os.environ.get("KDIL_OVERRIDE", "")
KNOBS = ("EMIS_NOX", "EMIS_VOC", "EMIS_ISO", "EMIS_CO",
         "KDIL_DAY", "KDIL_NIGHT", "BG_O3", "BG_CO")


def read_fit(label):
    knobs = {}
    with open(os.path.join(D, "results_fit_%s%s.txt" % (label, FIT_SUF))) as fh:
        for line in fh:
            if "=" not in line:
                continue
            k, _, v = line.partition("=")
            if k.strip() in KNOBS:
                knobs[k.strip()] = float(v)
    return knobs


def tminmax(label):
    T = []
    with open(os.path.join(D, "met_%s.csv" % label)) as fh:
        for r in fh:
            p = r.strip().split(",")
            if p[0] == "time":
                continue
            try:
                T.append(float(p[1]))
            except (ValueError, IndexError):
                pass
    return min(T) + 273.15, max(T) + 273.15


def run(cfg, tag):
    m = kpp_box.run_scenario(cfg, "nightsink_" + tag, "scale")
    return metrics(m["traj"])


def metrics(traj):
    tt = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    dt = tt[1] - tt[0]
    w = max(1, int(round(8.0 / dt)))
    best = float("-inf")
    for i in range(w - 1, len(o3)):
        if tt[i] < 48.0:
            continue
        a = sum(o3[i - w + 1:i + 1]) / w
        if a > best:
            best = a
    sel = [i for i in range(len(tt)) if tt[i] >= 48.0]
    return best, sum(o3[i] for i in sel) / len(sel)


def published(label):
    path = os.path.join(D, "results_%s_lvg.txt" % label.split("_")[0])
    mda8, rows = None, {}
    with open(path) as fh:
        for line in fh:
            if "base v48" in line and "mda8" in line:
                mda8 = float(line.split("mda8 =")[1].split()[0])
                continue
            p = line.split()
            if len(p) >= 6 and p[0] == "NOXW":
                rows[p[1]] = float(p[3])
    return mda8, rows


L = []
add = L.append
add("Nocturnal-sink test: comparison = published anchor")
if FIT_SUF:
    add("anchor = results_fit_<label>%s.txt (re-optimised with the pinned sink)" % FIT_SUF)
else:
    add("anchor = production fit; KDIL_NIGHT -> %s (no refit, upper bound)"
        % (KDIL_OVERRIDE or "fitted"))
add("")
for label in LABELS:
    if FIT_SUF:
        anchor = read_fit(label)
        tmin, tmax = tminmax(label)
        base = dict(anchor, BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
                    TMIN_K=tmin, TMAX_K=tmax)
    else:
        base, anchor = build_base(label)
    cfg = dict(base)
    if KDIL_OVERRIDE:
        cfg["KDIL_NIGHT"] = float(KDIL_OVERRIDE)
    m0 = kpp_box.run_scenario(cfg, "nightsink_base_" + label, "scale")
    mda8, v48 = metrics(m0["traj"])
    ref80, refrows = published(label)
    add("=" * 72)
    add("%s   all-on MDA8 %.4f -> %.4f  (%+.4f ppb)" %
        (label, ref80, mda8, mda8 - ref80))
    add("  window     published   proxy     delta")
    pos = 0
    maxd = 0.0
    for i, wl in enumerate(("0-3", "3-6", "6-9", "9-12", "12-15", "15-18",
                            "18-21", "21-24")):
        nox = [1.0] * NW
        nox[i] = 0.0
        mi = kpp_box.run_scenario(dict(cfg, NOXW=nox), "nightsink_n%d_%s" % (i, label), "scale")
        mi8, _ = metrics(mi["traj"])
        d = mi8 - mda8
        ref = refrows[wl]
        pos += 1 if ref > 0 else 0
        maxd = max(maxd, abs(d - ref))
        add("  %-6s %10.3f %9.3f %+9.3f" % (wl, ref, d, d - ref))
    add("  NOx windows raising MDA8: %d/8 (published 6/8 Dallas, 0/8 Deer Park)" % pos)
    add("  max |finite effect shift| = %.3f ppb" % maxd)
    add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_night_sink_proxy.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_night_sink_proxy.txt")
