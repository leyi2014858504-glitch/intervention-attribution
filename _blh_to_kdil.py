#!/usr/bin/env python3
"""
Bound the box's dilution timescales with observed boundary-layer height.

The model treats dilution as an exchange with background air:

    K_DIL = 1 / (KDIL_DAY + (KDIL_NIGHT - KDIL_DAY) * (1 - SDAY))     [1/s]
    dC/dt = E - K_DIL * (C - CBG)

so KDIL_DAY / KDIL_NIGHT are TIMESCALES IN SECONDS, not rates.

BLH gives us the daytime deepening rate, hence a physical UPPER bound:
    tau_BLH = H / (dH/dt)
If the air above the boundary layer is clean, entrainment dilutes the column
mean at exactly that rate; horizontal ventilation and entrainment of not-so-clean
air make the real exchange faster, i.e. tau smaller. So tau_BLH is an upper bound
and the fit is allowed to move inside [0.1*tau_BLH, tau_BLH].

Output: results_kdil.txt
"""
import csv, glob, os, re, statistics

D = "/mnt/d/atmos"
GROWTH_HOURS = range(8, 15)          # local hours when the BLH is deepening


def load(path):
    rows = []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            rows.append(r)
    return rows


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def analyse(path, label):
    rows = load(path)
    H = [f(r["boundary_layer_height"]) for r in rows]
    t = [r["time"] for r in rows]
    hr = [int(x[11:13]) for x in t]

    day = [h for h, i in zip(H, hr) if h is not None and 10 <= i <= 16]
    night = [h for h, i in zip(H, hr) if h is not None and (i <= 4 or i >= 22)]

    # tau from the morning deepening, hour by hour
    taus = []
    for i in range(1, len(H) - 1):
        if H[i] is None or H[i - 1] is None or H[i + 1] is None:
            continue
        if hr[i] not in GROWTH_HOURS:
            continue
        dHdt = (H[i + 1] - H[i - 1]) / (2 * 3600.0)      # m/s
        if dHdt > 0.005 and H[i] > 50:
            taus.append(H[i] / dHdt)

    print(f"--- {label} ({os.path.basename(path)}) ---")
    print(f"    days covered     : {len(set(x[:10] for x in t))}")
    print(f"    BLH day  (10-16h): mean {statistics.mean(day):7.0f} m   max {max(day):7.0f} m")
    print(f"    BLH night(22-04h): mean {statistics.mean(night):7.0f} m   min {min(night):7.0f} m")
    if taus:
        med = statistics.median(taus)
        print(f"    tau_BLH from morning deepening: median {med:8.0f} s "
              f"({med/3600:.2f} h), n={len(taus)}, p10={sorted(taus)[len(taus)//10]:.0f}s")
    else:
        med = None
        print("    tau_BLH: no usable deepening segments")
    return dict(label=label, day=statistics.mean(day), night=statistics.mean(night),
                tau=med, tau_lo=(med * 0.1 if med else None),
                tau_hi=med, ratio=(statistics.mean(day) / statistics.mean(night)))


if __name__ == "__main__":
    pats = ["met_dallas_2016w", "met_dallas_2017s",
            "met_deerpark_2016w", "met_deerpark_2017s",
            "met_houston_2016w", "met_houston_2017s",
            "met_elpaso_2016w", "met_elpaso_2017s"]
    out = []
    for p in pats:
        fp = os.path.join(D, p + ".csv")
        if os.path.exists(fp):
            out.append(analyse(fp, p))
    with open(os.path.join(D, "results_kdil.txt"), "w", encoding="utf-8") as fh:
        fh.write("label             BLH_day  BLH_night  day/night  tau_BLH(s)  "
                 "suggested KDIL_DAY range\n")
        for r in out:
            if r["tau"]:
                fh.write(f"{r['label']:<18}{r['day']:8.0f}{r['night']:11.0f}"
                         f"{r['ratio']:11.1f}{r['tau']:12.0f}   "
                         f"[{r['tau_lo']:.0f}, {r['tau_hi']:.0f}]\n")
            else:
                fh.write(f"{r['label']:<18}{r['day']:8.0f}{r['night']:11.0f}"
                         f"{r['ratio']:11.1f}{'n/a':>12}\n")
    print("wrote results_kdil.txt")
