#!/usr/bin/env python3
"""
Final site screen: everything we need must be LIVE, not merely present.

The two previous screens both failed the same way -- they accepted parameters
that exist in the inventory but whose instruments died years ago:
  * Beaumont's PAMS flag came from a monitor that ended in 2001 (LA) / 2007 (Charlotte)
  * Beaumont CO ended 1996-12-31, Fort Worth CO ended 2014-10-31

So this pass requires, per site, that ALL of
    44201 O3, 42601 NO, 42602 NO2, 42101 CO        (hourly criteria gases)
and an auto-GC VOC analyser                        (hourly speciated VOC)
have a last sample date at or after CUTOFF.

Output: results_sites_final.txt
"""
import csv, os
from collections import defaultdict

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
SIT = os.path.join(D, "aqs_sites.csv")
OUT = "/mnt/d/atmos/results_sites_final.txt"

CUTOFF = "2016-06-01"          # must be alive across both study windows
CORE = {"44201": "O3", "42601": "NO", "42602": "NO2", "42101": "CO"}
PAMS_HOURLY = ("AUTO GC", "DUAL FID", "PE8700", "AUX GC")
TOXICS_24H = ("CRYOGENIC PRECON", "GC/MS", "GC-ECD", "ANASORB", "TD/AGILENT")


def is_voc(pc):
    return (pc.startswith("43") and pc != "43000") or pc.startswith("45")


def main():
    last = defaultdict(dict)          # site -> param -> max last-sample-date
    pams = defaultdict(lambda: {"last": "", "n": 0})

    with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            pc = row["Parameter Code"]
            key = (row["State Code"], row["County Code"], row["Site Number"])
            ld = (row.get("Last Sample Date") or "").strip()
            if pc in CORE:
                prev = last[key].get(pc, "")
                if ld > prev:
                    last[key][pc] = ld
            if is_voc(pc):
                m = (row.get("Last Method") or "").strip().upper()
                if any(s in m for s in PAMS_HOURLY) and not any(t in m for t in TOXICS_24H):
                    if ld > pams[key]["last"]:
                        pams[key]["last"] = ld

    meta = {}
    with open(SIT, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            meta[(row["State Code"], row["County Code"], row["Site Number"])] = row

    ok = []
    for key in pams:
        if pams[key]["last"] < CUTOFF:
            continue
        core = {pc: last[key].get(pc, "") for pc in CORE}
        if all(v >= CUTOFF for v in core.values()):
            ok.append((key, pams[key]["last"], core))

    ok.sort(key=lambda t: t[1], reverse=True)
    log = open(OUT, "w", encoding="utf-8")
    log.write(f"cutoff: every core gas and the auto-GC VOC must reach {CUTOFF}\n")
    log.write(f"live auto-GC VOC sites: {sum(1 for k in pams if pams[k]['last'] >= CUTOFF)}\n")
    log.write(f"SITES WITH ALL FIVE, LIVE: {len(ok)}\n\n")
    for key, pl, core in ok:
        m = meta.get(key, {})
        log.write(f"{key}  pamsVOC_last={pl}\n")
        log.write(f"    core last-sample-dates: " +
                  ", ".join(f"{CORE[pc]}={core[pc]}" for pc in CORE) + "\n")
        log.write(f"    {(m.get('City Name') or '').strip()}, {(m.get('State Name') or '').strip()}"
                  f" | {(m.get('Local Site Name') or '').strip()}"
                  f" | {m.get('Latitude','')},{m.get('Longitude','')}\n")
        log.write(f"    CBSA={(m.get('CBSA Name') or '').strip()}\n\n")
    log.close()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
