#!/usr/bin/env python3
"""
Turn a pulled AQS CSV into an hourly matrix we can fit against.

Handles three things that will silently corrupt a fit if ignored:

1. MULTIPLE POCs. A parameter can come from several instruments (different POC).
   Mixing them creates artificial jumps, so for each parameter we keep the POC
   with the most usable hourly rows and discard the rest.

2. UNITS. Criteria gases come back in "Parts per billion" (ppb), but the VOC
   species come back in "Parts per billion Carbon" (ppbC). ppbC is NOT ppb:
   toluene at 5 ppbC is 5/7 ppb. The units of every species are recorded in
   results_obs_<label>_units.txt so the CB05 mapping step can divide by nC.

3. DERIVED SPECIES. NOx is not measured; it is NO + NO2. We build it here so the
   box model has something to match.

Usage:
  python3 _obs_build.py <label>
Reads  aqs_sites_data/site_<label>.csv
Writes obs_<label>_hourly.csv, results_obs_<label>.txt, results_obs_<label>_units.txt
"""
import sys, csv, os, collections

D = "/mnt/d/atmos"
SRC = os.path.join(D, "aqs_sites_data")
CORE = {"42101": "CO", "42601": "NO", "42602": "NO2", "44201": "O3"}


def main():
    label = sys.argv[1]
    src = os.path.join(SRC, f"site_{label}.csv")
    rows = list(csv.DictReader(open(src, newline="", encoding="utf-8", errors="replace")))
    print(f"{label}: {len(rows)} hourly rows in {src}")

    # --- what POCs exist per parameter, and how many rows each carries
    poc_count = collections.defaultdict(collections.Counter)
    units = {}
    names = {}
    for r in rows:
        pc = r["parameter_code"]
        poc_count[pc][r["poc"]] += 1
        units[pc] = r["units_of_measure"]
        names[pc] = r["parameter"]

    # --- keep the best POC per parameter
    best_poc = {pc: c.most_common(1)[0][0] for pc, c in poc_count.items()}

    # --- pivot to (datetime) -> {species: value}
    series = collections.defaultdict(dict)
    for r in rows:
        pc = r["parameter_code"]
        if r["poc"] != best_poc[pc]:
            continue
        try:
            v = float(r["sample_measurement"])
        except (TypeError, ValueError):
            continue
        dt = f"{r['date_local']} {r['time_local']}"
        series[dt][pc] = v

    times = sorted(series)
    codes = sorted(poc_count)
    out = os.path.join(D, f"obs_{label}_hourly.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["datetime"] + codes + ["NOx"])
        for t in times:
            d = series[t]
            nox = None
            if "42601" in d and "42602" in d:
                nox = d["42601"] + d["42602"]
            w.writerow([t] + [d.get(c, "") for c in codes] + ["" if nox is None else nox])

    # --- report
    rep = open(os.path.join(D, f"results_obs_{label}.txt"), "w", encoding="utf-8")
    rep.write(f"label: {label}\nrows in: {len(rows)}\n"
              f"distinct hourly timestamps: {len(times)}\n"
              f"time span: {times[0]} .. {times[-1]}\n\n")
    rep.write(f"{'code':<7}{'n':>5}  {'poc':>4}  {'unit':<38} name\n")
    rep.write("-" * 110 + "\n")
    for pc in codes:
        rep.write(f"{pc:<7}{poc_count[pc][best_poc[pc]]:>5}  {best_poc[pc]:>4}  "
                  f"{(units.get(pc) or ''):<38} {names.get(pc,'')}\n")

    # diurnal means for the species that matter most
    rep.write("\n=== diurnal means (ppb unless noted) ===\n")
    for pc, tag in [("44201", "O3"), ("42601", "NO"), ("42602", "NO2"),
                    ("42101", "CO"), ("45202", "toluene(ppbC)"),
                    ("43243", "isoprene(ppbC)"), ("45201", "benzene(ppbC)")]:
        if pc not in poc_count:
            rep.write(f"  {tag:<16} not present\n")
            continue
        byhr = collections.defaultdict(list)
        for t in times:
            v = series[t].get(pc)
            if v is None:
                continue
            byhr[int(t[11:13])].append(v)
        prof = []
        for h in range(24):
            xs = byhr.get(h, [])
            prof.append(f"{sum(xs)/len(xs):7.2f}" if xs else "      .")
        rep.write(f"  {tag:<16} " + " ".join(prof) + "\n")
    rep.write("  hour            " + " ".join(f"{h:>7d}" for h in range(24)) + "\n")
    rep.close()

    unitsf = os.path.join(D, f"results_obs_{label}_units.txt")
    with open(unitsf, "w", encoding="utf-8") as fh:
        fh.write(f"# units per parameter for {label}: ppb vs ppbC matters, ppbC must be divided by nC\n")
        for pc in codes:
            fh.write(f"{pc}\t{best_poc[pc]}\t{units.get(pc)}\t{names.get(pc,'')}\n")
    print("wrote", out, "and the two results_obs_* files")


if __name__ == "__main__":
    main()
