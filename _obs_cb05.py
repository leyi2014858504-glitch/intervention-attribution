#!/usr/bin/env python3
"""
Observations -> hourly CB05 lumped amounts, ready to fit against.

Merges the three earlier pieces:
  * _obs_build.py   best-POC selection, pivoting, NOx = NO + NO2
  * _voc_map.py     species -> CB05 group (with per-row confidence)
  * _voc_nc.py      carbon number, and ppbC -> carbon bonds

The unit chain, which is the whole reason this file exists:
  AQS reports VOC in "Parts per billion Carbon"
  -> molecules    = ppbC / nC
  -> carbon bonds = group-dependent (PAR gets nC bonds; OLE/IOLE get 1 bond plus
                    nC-2 PAR bonds; TOL/XYL/ISOP/FORM/ALD2/ALDX/ETHA/MEOH/ETOH
                    are per-molecule lumped species, so 1:1)

Outputs
  obs_<label>_cb05.csv    datetime + observed cores + CB05 lumped amounts
  results_obs_<label>_cb05.txt   diurnal composites and the species ledger
"""
import sys, csv, os, collections, importlib.util

D = "/mnt/d/atmos"
SRC = os.path.join(D, "aqs_sites_data")

CORE = {"42101": "CO", "42601": "NO", "42602": "NO2", "44201": "O3"}

# AQS mixes concentration units between species. O3 and CO come back in PPM
# (0.02 ppm = 20 ppb), while NO and NO2 come back in PPB. Reading 0.021 as
# "0.021 ppb of ozone" would silently destroy the fit, so scale explicitly.
UNIT_SCALE = {
    "parts per million": 1000.0,
    "parts per billion": 1.0,
    "parts per billion carbon": 1.0,     # VOC: carbon handled separately via nC
    "micrograms/cubic meter": None,      # not usable without a molar conversion
}
# CB05 groups we emit, in output order
GROUPS = ["PAR", "OLE", "IOLE", "ETH", "ETHA", "TOL", "XYL", "ISOP",
          "FORM", "ALD2", "ALDX", "APIN", "BPIN", "LIM", "MEOH", "ETOH", "MGLY"]


def _load(mod, path):
    spec = importlib.util.spec_from_file_location(mod, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main(label):
    vm = _load("vocmap", os.path.join(D, "_voc_map.py"))
    nc = _load("vocnc", os.path.join(D, "_voc_nc.py"))

    rows = list(csv.DictReader(open(os.path.join(SRC, f"site_{label}.csv"),
                                    newline="", encoding="utf-8", errors="replace")))
    if not rows:
        print(f"{label}: no rows"); return

    # --- best POC per parameter
    poc = collections.defaultdict(collections.Counter)
    for r in rows:
        poc[r["parameter_code"]][r["poc"]] += 1
    best = {pc: c.most_common(1)[0][0] for pc, c in poc.items()}

    # --- route every parameter
    routed, skipped = {}, []
    for r in rows:
        pc = r["parameter_code"]
        if pc in routed or pc in [s[0] for s in skipped]:
            continue
        name = r["parameter"]
        if pc in CORE:
            routed[pc] = ("CORE", CORE[pc], None, 0.0)
            continue
        grp, conf, note = vm.classify(name)
        if grp in ("-", "?"):
            skipped.append((pc, name, grp, conf))
            continue
        nC = nc.carbon_number(name)
        if not nC:
            skipped.append((pc, name, "no-nC", conf))
            continue
        routed[pc] = ("VOC", grp, nC, 0.0)

    # --- accumulate hourly
    series = collections.defaultdict(lambda: collections.defaultdict(float))
    units = {}
    unit_hits = collections.Counter()
    for r in rows:
        pc = r["parameter_code"]
        if pc not in routed or r["poc"] != best[pc]:
            continue
        u = (r["units_of_measure"] or "").strip()
        units[pc] = u
        scale = UNIT_SCALE.get(u.lower(), 1.0)
        unit_hits[u] += 1
        try:
            v = float(r["sample_measurement"])
        except (TypeError, ValueError):
            continue
        kind, grp, nC, _ = routed[pc]
        dt = f"{r['date_local']} {r['time_local']}"
        if kind == "CORE":
            if scale is None:
                continue
            series[dt][grp] += v * scale
        else:
            for g, amt in nc.bonds(v, nC, grp).items():
                series[dt][g] += amt

    # NOx from NO + NO2
    for dt, d in series.items():
        if "NO" in d and "NO2" in d:
            d["NOx"] = d["NO"] + d["NO2"]

    cols = ["O3", "NO", "NO2", "NOx", "CO"] + GROUPS
    times = sorted(series)
    out = os.path.join(D, f"obs_{label}_cb05.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["datetime"] + cols)
        for t in times:
            d = series[t]
            w.writerow([t] + [round(d.get(c), 4) if c in d else "" for c in cols])

    # --- report
    rep = open(os.path.join(D, f"results_obs_{label}_cb05.txt"), "w", encoding="utf-8")
    rep.write(f"label {label}: {len(rows)} rows -> {len(times)} hourly timestamps\n")
    rep.write(f"routed parameters: {len(routed)}   skipped: {len(skipped)}\n")
    rep.write("unit strings seen (raw rows): "
              + ", ".join(f"{k!r}={v}" for k, v in unit_hits.most_common()) + "\n\n")
    rep.write("=== diurnal composites (mean by local hour) ===\n")
    for c in cols:
        byhr = collections.defaultdict(list)
        for t in times:
            v = series[t].get(c)
            if v is not None:
                byhr[int(t[11:13])].append(v)
        if not byhr:
            rep.write(f"{c:<6} (empty)\n")
            continue
        prof = [f"{(sum(byhr[h])/len(byhr[h])):8.3f}" if byhr.get(h) else "       ."
                for h in range(24)]
        rep.write(f"{c:<6} " + " ".join(prof) + "\n")
    rep.write("hour   " + " ".join(f"{h:>8d}" for h in range(24)) + "\n\n")
    rep.write("=== skipped parameters ===\n")
    for pc, name, grp, conf in skipped:
        rep.write(f"  {pc}  {grp:<3} {conf}  {name}\n")
    rep.close()
    print(f"{label}: wrote {out} ({len(routed)} routed, {len(skipped)} skipped)")


if __name__ == "__main__":
    main(sys.argv[1])
