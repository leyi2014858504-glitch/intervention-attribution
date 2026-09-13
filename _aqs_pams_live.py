#!/usr/bin/env python3
"""
Is ANY hourly PAMS auto-GC VOC still running?

What we just learned, the hard way:
  * AQS bulk publishes hourly files for only 8 parameters, all criteria pollutants.
    There is no hourly VOC bulk file at all.
  * Through the API, the aromatic/VOC parameters come back as "24 HOUR" with
    sample_frequency "EVERY 6TH DAY" -- that is the NATTS Air-Toxics canister,
    not PAMS.
  * Charlotte's toluene POC 1 is the real PAMS auto-GC
    ("...PE8700; AUTO GC; SUBAMBIENT-DUAL FID") and it DIED on 2007-12-26.

So the question this script answers: over all 20,856 sites and 370,090 monitor
records, which sites still have an auto-GC / dual-FID PAMS VOC monitor whose last
sample date reaches 2016 or later?

Output: results_pams_live.txt
"""
import csv, os
from collections import defaultdict

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
OUT = "/mnt/d/atmos/results_pams_live.txt"
SIT = os.path.join(D, "aqs_sites.csv")

# the auto-GC / dual-FID signatures of a real hourly PAMS VOC analyser
PAMS_HOURLY = ("AUTO GC", "DUAL FID", "PE8700", "AUX GC")
# 24-hour canister toxics: explicitly NOT what we want
TOXICS_24H = ("CRYOGENIC PRECON", "GC/MS", "GC-ECD", "ANASORB", "TD/AGILENT")


def is_voc(pc):
    return (pc.startswith("43") and pc != "43000") or pc.startswith("45")


def main():
    hits = defaultdict(lambda: {"n": 0, "params": set(), "last": "", "meth": set()})
    any_voc_recent = defaultdict(lambda: {"n": 0, "params": set(), "last": ""})

    with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            pc = row["Parameter Code"]
            if not is_voc(pc):
                continue
            key = (row["State Code"], row["County Code"], row["Site Number"])
            ld = (row.get("Last Sample Date") or "").strip()
            meth = (row.get("Last Method") or "").strip().upper()

            # any VOC parameter at all that is still alive
            if ld >= "2016":
                a = any_voc_recent[key]
                a["n"] += 1
                a["params"].add(pc)
                a["last"] = max(a["last"], ld)

            if any(sig in meth for sig in PAMS_HOURLY) and not any(t in meth for t in TOXICS_24H):
                h = hits[key]
                h["n"] += 1
                h["params"].add(pc)
                h["last"] = max(h["last"], ld)
                h["meth"].add(meth)

    meta = {}
    with open(SIT, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            meta[(row["State Code"], row["County Code"], row["Site Number"])] = row

    log = open(OUT, "w", encoding="utf-8")
    log.write("=== sites with an auto-GC / dual-FID PAMS VOC monitor (any date) ===\n")
    alive = [(k, v) for k, v in hits.items() if v["last"] >= "2016"]
    dead = [(k, v) for k, v in hits.items() if v["last"] < "2016"]
    log.write(f"  total such sites: {len(hits)}   still alive in 2016+: {len(alive)}   "
              f"dead before 2016: {len(dead)}\n\n")

    log.write("=== STILL ALIVE in 2016+ ===\n")
    for k, v in sorted(alive, key=lambda t: -t[1]["n"]):
        m = meta.get(k, {})
        log.write(f"  {k}  nVOC={v['n']:>3}  last={v['last']}  "
                  f"{(m.get('City Name') or '').strip()}, {(m.get('State Name') or '').strip()}"
                  f" | {(m.get('Local Site Name') or '').strip()}\n")
        for s in sorted(v["meth"])[:2]:
            log.write(f"       {s}\n")
    log.write("\n")

    log.write("=== last sample date distribution for those PAMS sites ===\n")
    yr = defaultdict(int)
    for k, v in hits.items():
        yr[v["last"][:4] or "?"] += 1
    for y in sorted(yr):
        log.write(f"  {y}: {yr[y]}\n")
    log.write("\n")

    log.write(f"=== sites with ANY 43xxx/45xxx VOC parameter alive in 2016+: {len(any_voc_recent)} ===\n")
    log.write("  (top 30 by parameter count)\n")
    for k, v in sorted(any_voc_recent.items(), key=lambda t: -t[1]["n"])[:30]:
        m = meta.get(k, {})
        log.write(f"  {k}  nVOC={v['n']:>3}  last={v['last']}  "
                  f"{(m.get('City Name') or '').strip()}, {(m.get('State Name') or '').strip()}\n")
    log.close()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
