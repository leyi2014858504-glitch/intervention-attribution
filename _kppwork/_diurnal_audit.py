#!/usr/bin/env python3
"""Diurnal audit of the four baselines.

Checks the items raised in review:
  1. unit chain: ppb/h -> molec/cm3/s, NO/NO2 and VOC splits, profile factor
  2. full diurnal cycle: O3, NOx, NO, NO2, PAR, ISOP for days 3-5
  3. day-to-day drift (is the 2-day spin-up enough?)
  4. nocturnal storage: how much of the night-time accumulation survives to 06h,
     and how much of it is set by the 86 400 s exchange (sealed-night control)
  5. city model vs observed hourly composite (amplitude and shape)

Writes results_diurnal_audit.txt and results_diurnal_composite.csv
"""
import csv
import math
import os
import statistics
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402
from _gen16_city import build_base  # noqa: E402

D = "/mnt/d/atmos"
COLS = {"o3": 1, "nox": 2, "co": 7, "no": 8, "no2": 9,
        "par": 10, "tol": 11, "xyl": 12, "isop": 13}
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
TEND = 432000.0          # five days; days 3-5 are diagnostic


def synth_cfg(name):
    cfg = dict(TEND_S=TEND)
    if name == "r_lim1":
        cfg.update(NOX_FACTOR=0.5, VOC_FACTOR=1.5)
    return cfg


def base_of(game):
    if game in ("r_sat", "r_lim1"):
        return synth_cfg(game)
    base, _ = build_base(game)
    cfg = dict(base)
    cfg["TEND_S"] = TEND
    return cfg


def hourly(traj, days, col):
    acc = {h: [] for h in range(24)}
    for r in traj:
        t = r[0]
        d = int(t // 24) + 1
        if d < days[0] or d > days[-1]:
            continue
        acc[int(round(t)) % 24].append(r[COLS[col]])
    return {h: (statistics.mean(v) if v else float("nan")) for h, v in acc.items()}


def day_mean(traj, day, col):
    """day 3 means t in [48,72) h."""
    lo = (day - 1) * 24.0
    v = [r[COLS[col]] for r in traj if lo <= r[0] < lo + 24.0]
    return statistics.mean(v) if v else float("nan")


def mda8(traj):
    o3 = [r[1] for r in traj]
    t = [r[0] for r in traj]
    dt = t[1] - t[0]
    w = max(1, int(round(8.0 / dt)))
    best = -1e9
    for i in range(w - 1, len(o3)):
        if t[i] < 48.0:
            continue
        best = max(best, sum(o3[i - w + 1:i + 1]) / w)
    return best


def gp(mu, sd, h):
    return math.exp(-(((h - mu) / sd) ** 2))


prof = [0.45 + 1.30 * (gp(8.0, 1.6, h / 60.0) + gp(18.5, 2.0, h / 60.0))
        for h in range(24 * 60)]
PR_MIN, PR_MAX, PR_MEAN = min(prof), max(prof), statistics.mean(prof)

L = []
add = L.append
add("Diurnal audit of the four baselines (production build, 5-day runs)")
add("profile factor (weekday): min %.3f  peak %.3f  daily mean %.3f" %
    (PR_MIN, PR_MAX, PR_MEAN))
add("unit chain: E[ppb/h] * profile * N_AIR/1e9 / 3600 -> molec cm-3 s-1")
add("splits: NO/NO2 0.75/0.25; VOC -> PAR/TOL/XYL 0.55/0.25/0.20 by mixing ratio")
add("")

runs = {}
for g in GAMES:
    cfg = base_of(g)
    runs[g] = kpp_box.run_scenario(cfg, "audit_" + g, "anchor")["traj"]

add("A. Day-to-day drift (daily means, days 3/4/5)")
add("%-16s %-5s %9s %9s %9s %8s" % ("game", "sp", "day3", "day4", "day5", "drift%"))
for g in GAMES:
    for sp in ("o3", "nox", "no2", "par", "isop"):
        d3 = day_mean(runs[g], 3, sp)
        d4 = day_mean(runs[g], 4, sp)
        d5 = day_mean(runs[g], 5, sp)
        dr = 100 * (d5 - d3) / d3 if d3 else float("nan")
        add("%-16s %-5s %9.3f %9.3f %9.3f %8.2f" % (g, sp, d3, d4, d5, dr))
    add("")
add("B. Nocturnal storage and a bounded five-fold slower night exchange")
add("exchange-only retention over 12 h at tau=86400 s = %.3f" % math.exp(-12 * 3600 / 86400))
add("%-16s %9s %9s %9s %10s %10s" %
    ("game", "NOx18", "NOx06", "ratio", "K5_NOx06", "K5_MDA8"))
night = {}
for g in GAMES:
    t = runs[g]
    h = hourly(t, (3, 4, 5), "nox")
    r18, r06 = h[18], h[6]
    cfg = base_of(g)
    cfg["KDIL_NIGHT"] = 5.0 * 86400.0
    slow = kpp_box.run_scenario(cfg, "audit_k5_" + g, "anchor")["traj"]
    hs = hourly(slow, (3, 4, 5), "nox")
    night[g] = (r18, r06, hs[6], mda8(t), mda8(slow))
    add("%-16s %9.3f %9.3f %9.3f %10.3f %10.3f" %
        (g, r18, r06, r06 / r18, hs[6], night[g][4]))
add("")
add("C. Night-time ozone loss and modelled MDA8 under the five-fold slower night")
add("%-16s %9s %9s %9s %9s %9s" % ("game", "O3_18", "O3_06", "ratio", "MDA8", "K5"))
for g in GAMES:
    hO = hourly(runs[g], (3, 4, 5), "o3")
    v8, v8s = night[g][3], night[g][4]
    add("%-16s %9.3f %9.3f %9.3f %9.3f %9.3f" %
        (g, hO[18], hO[6], hO[6] / hO[18], v8, v8s))
add("")

add("D. City composite against the observed hourly mean")
for g in ("dallas_2016w", "deerpark_2017s"):
    obs = os.path.join(D, "obs_%s_cb05.csv" % ("dallas_2016w" if g == "dallas_2016w" else "deerpark_2017s"))
    acc = {sp: {h: [] for h in range(24)} for sp in ("o3", "nox")}
    with open(obs) as fh:
        for r in csv.DictReader(fh):
            hh = int(r["datetime"][11:13])
            for sp, key in (("o3", "O3"), ("nox", "NOx")):
                try:
                    acc[sp][hh].append(float(r[key]))
                except (TypeError, ValueError):
                    pass
    add("%s" % g)
    add("  %-3s %10s %10s %10s %10s" % ("h", "O3mod", "O3obs", "NOxmod", "NOxobs"))
    for h in range(0, 24, 3):
        mo = hourly(runs[g], (3, 4, 5), "o3")[h]
        mn = hourly(runs[g], (3, 4, 5), "nox")[h]
        oo = statistics.mean(acc["o3"][h]) if acc["o3"][h] else float("nan")
        on = statistics.mean(acc["nox"][h]) if acc["nox"][h] else float("nan")
        add("  %-3d %10.3f %10.3f %10.3f %10.3f" % (h, mo, oo, mn, on))
    add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_diurnal_audit.txt"), "w") as fh:
    fh.write(txt + "\n")

with open(os.path.join(D, "results_diurnal_composite.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["game", "hour", "o3", "nox", "no", "no2", "co", "par", "tol", "xyl", "isop"])
    for g in GAMES:
        hs = {sp: hourly(runs[g], (3, 4, 5), sp) for sp in COLS}
        for h in range(24):
            wr.writerow([g, h] + ["%.4f" % hs[sp][h] for sp in
                                  ("o3", "nox", "no", "no2", "co", "par", "tol", "xyl", "isop")])
print("wrote results_diurnal_audit.txt / results_diurnal_composite.csv")
