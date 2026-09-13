#!/usr/bin/env python3
"""
Inventory vetting: dump EVERY parameter the two chosen sites report, with the
monitor type / network / method strings that decide whether it is hourly.

Why: the first screen only looked at 43xxx, which silently dropped the aromatic
hydrocarbons -- benzene/toluene/ethylbenzene/xylene/styrene live in 452xx
(45201 Benzene, 45202 Toluene, 45203 Ethylbenzene, 45204 o-Xylene, 45220 Styrene).
Those are not optional: toluene and xylene are the dominant aromatic reactants in
a VOC-limited urban regime.

Output: results_site_params.txt
"""
import csv, os
from collections import defaultdict

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
OUT = "/mnt/d/atmos/results_site_params.txt"

SITES = {
    "LA_06-037-1103":        ("06", "037", "1103"),
    "Charlotte_37-119-0041": ("37", "119", "0041"),
}


def main():
    # site -> code -> {name, types, networks, methods, first, last, npoc}
    acc = {k: defaultdict(lambda: {"name": "", "types": set(), "net": set(),
                                   "meth": set(), "first": 9999, "last": "", "npoc": 0})
           for k in SITES}

    with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            key = (row["State Code"], row["County Code"], row["Site Number"])
            for lab, tgt in SITES.items():
                if key != tgt:
                    continue
                pc = row["Parameter Code"]
                a = acc[lab][pc]
                a["name"] = (row.get("Parameter Name") or "").strip()
                if row.get("Monitor Type"):
                    a["types"].add(row["Monitor Type"].strip())
                if row.get("Networks"):
                    a["net"].add(row["Networks"].strip())
                m = (row.get("Last Method") or "").strip()
                if m and m != " - ":
                    a["meth"].add(m)
                fy = (row.get("First Year of Data") or "").strip()
                if fy.isdigit():
                    a["first"] = min(a["first"], int(fy))
                ld = (row.get("Last Sample Date") or "").strip()
                if ld:
                    a["last"] = max(a["last"], ld)
                a["npoc"] += 1

    log = open(OUT, "w", encoding="utf-8")
    for lab in SITES:
        log.write(f"########## {lab}: {len(acc[lab])} distinct parameters ##########\n")
        pref = defaultdict(int)
        for pc in acc[lab]:
            pref[pc[:2] + "xxx"] += 1
        log.write("prefix histogram: " + ", ".join(f"{p}={n}" for p, n in sorted(pref.items())) + "\n\n")
        for pc in sorted(acc[lab]):
            a = acc[lab][pc]
            log.write(f"{pc}  {a['name']}\n")
            log.write(f"      type={sorted(a['types'])}  network={sorted(a['net'])}  "
                      f"poc={a['npoc']}  years={a['first']}..{a['last']}\n")
            if a["meth"]:
                log.write(f"      methods={sorted(a['meth'])[:3]}\n")
        log.write("\n")
    log.close()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
