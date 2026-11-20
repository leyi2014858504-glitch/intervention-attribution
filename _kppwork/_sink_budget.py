#!/usr/bin/env python3
"""Sink budget: how large is the dry-deposition term the box currently omits?

The box's only sink is first-order dilution (KDIL, fitted).  If deposition were
added, the loss rate per species would be k_dep = v_d / H(t).  With the observed
hourly boundary-layer height we can compare the two terms directly, per city,
day and night, without running the model.

Outputs (per city): median H by period, tau_dil, tau_dep for literature v_d
ranges, and the deposition share of the total first-order sink.

Usage: python3 _sink_budget.py
Writes results_sink_budget.txt
"""
import csv
import os
import statistics
import sys

D = "/mnt/d/atmos"
LABELS = ("dallas_2016w", "deerpark_2017s")
# literature dry-deposition velocities, m/s (grass/urban, O3/NO2/CO)
VD = {"O3": (0.003, 0.005), "NO2": (0.005, 0.020), "CO": (0.0003, 0.001)}
KNOBS = ("KDIL_DAY", "KDIL_NIGHT")


def read_met(label):
    H = {h: [] for h in range(24)}
    with open(os.path.join(D, "met_%s.csv" % label), newline="") as fh:
        for r in csv.DictReader(fh):
            h = int(r["time"][11:13])
            try:
                H[h].append(float(r["boundary_layer_height"]))
            except (TypeError, ValueError):
                pass
    return H


def read_fit(label):
    d = {}
    with open(os.path.join(D, "results_fit_%s.txt" % label)) as fh:
        for line in fh:
            if "=" in line:
                k, _, v = line.partition("=")
                if k.strip() in KNOBS:
                    d[k.strip()] = float(v)
    return d["KDIL_DAY"], d["KDIL_NIGHT"]


def med(xs):
    ys = [x for x in xs if x == x]
    return statistics.median(ys) if ys else float("nan")


L = []
add = L.append
add("First-order sink budget with and without dry deposition")
add("k_dil = 1/tau_dil (fitted)   k_dep = v_d / H (literature v_d, observed H)")
add("share = k_dep / (k_dil + k_dep)")
add("")
for label in LABELS:
    H = read_met(label)
    tau_d, tau_n = read_fit(label)
    Hd = med([x for h in range(7, 19) for x in H[h]])       # 07-19 local
    Hn = med([x for h in list(range(0, 7)) + list(range(19, 24)) for x in H[h]])
    add("=" * 78)
    add("%s   KDIL_DAY = %.0f s (tau %.2f h)   KDIL_NIGHT = %.0f s (tau %.2f h)"
        % (label, tau_d, tau_d / 3600.0, tau_n, tau_n / 3600.0))
    add("   median H: day %.0f m   night %.0f m" % (Hd, Hn))
    add("   %-4s %8s %10s %10s %10s %10s" %
        ("sp", "v_d", "tau_dep_day", "share_day", "tau_dep_night", "share_night"))
    for sp, (vlo, vhi) in VD.items():
        for v in (vlo, vhi):
            td = Hd / v
            tn = Hn / v
            sd = (1.0 / td) / (1.0 / tau_d + 1.0 / td)
            sn = (1.0 / tn) / (1.0 / tau_n + 1.0 / tn)
            add("   %-4s %7.3f %9.1f h %9.1f%% %9.1f h %9.1f%%" %
                (sp, v, td / 3600.0, 100 * sd, tn / 3600.0, 100 * sn))
    add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_sink_budget.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_sink_budget.txt")
