#!/usr/bin/env python3
"""
Pick candidate AQS sites for the observation-anchoring study.

`aqs_monitors.zip` is the MONITOR INVENTORY (not the sample-level file), so it
has no "Sample Duration" column. What it does have is `Monitor Type` and
`Networks`, and the PAMS network is exactly the one that reports HOURLY
speciated VOC. So we screen on:

  * the four core species present        O3 44201 / NO 42601 / NO2 42602 / CO 42101
    (NO and NO2 must exist separately -- collapsing them destroys the very
     NOx-vs-VOC distinction the study is about)
  * a wide 43xxx VOC panel (PAMS is ~57 species)
  * PAMS membership, or a 43xxx panel wide enough to be PAMS-like

24-HOUR Air-Toxics VOC data never enters, because it would erase the diurnal
cycle that the box model is anchored to.

Outputs:
  results_aqs_pick.txt     ranked candidate sites
  results_aqs_voc_panel.txt  every 43xxx parameter seen (for the CB05 mapping table)
"""
import csv, os
from collections import defaultdict, Counter

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
SIT = os.path.join(D, "aqs_sites.csv")
OUT = "/mnt/d/atmos/results_aqs_pick.txt"
OUT2 = "/mnt/d/atmos/results_aqs_voc_panel.txt"

log = open(OUT, "w", encoding="utf-8")
def say(*a):
    print(*a); print(*a, file=log)

# VOC code families: 43xxx = hydrocarbons/oxygenates, 45xxx = AROMATICS
# (45109 m/p-xylene, 45201 benzene, 45202 toluene, 45203 ethylbenzene,
#  45207/45208/45225 trimethylbenzenes, 45220 styrene). Missing 45xxx would
# silently drop every aromatic, which is exactly the VOC-limited lever.
CORE = {"44201": "O3", "42601": "NO", "42602": "NO2", "42101": "CO"}


def is_voc(pc):
    return (pc.startswith("43") and pc != "43000") or pc.startswith("45")

# ---------------- read monitors ----------------
# key -> set(param codes);  key -> {flags}
params = defaultdict(set)
pams = defaultdict(bool)
networks = defaultdict(Counter)
mtype = defaultdict(Counter)
firsty = {}
lastd = {}
names = {}
voc_panel = Counter()

with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
    r = csv.DictReader(fh)
    say("monitor columns:", r.fieldnames)
    nrow = 0
    for row in r:
        nrow += 1
        key = (row["State Code"], row["County Code"], row["Site Number"])
        pc = row["Parameter Code"]
        params[key].add(pc)
        names[pc] = row.get("Parameter Name", "")
        nt = (row.get("Monitor Type") or "").strip()
        nw = (row.get("Networks") or "").strip()
        if nt:
            mtype[key][nt] += 1
        if nw:
            networks[key][nw] += 1
        if "PAMS" in nt.upper() or "PAMS" in nw.upper():
            pams[key] = True
        fy = (row.get("First Year of Data") or "").strip()
        if fy.isdigit():
            firsty[key] = min(firsty.get(key, 9999), int(fy))
        ld = (row.get("Last Sample Date") or "").strip()
        if ld:
            lastd[key] = max(lastd.get(key, ""), ld)
        if is_voc(pc):
            voc_panel[pc] += 1

say("monitor rows:", nrow, " distinct sites:", len(params))
say("distinct 43xxx VOC parameter codes:", len(voc_panel))

with open(OUT2, "w", encoding="utf-8") as fh:
    fh.write(f"# {len(voc_panel)} distinct 43xxx VOC parameter codes in aqs_monitors\n")
    fh.write(f"# {'code':<6} {'n_monitors':>10}  name\n")
    for pc, n in sorted(voc_panel.items()):
        fh.write(f"  {pc:<6} {n:>10}  {names.get(pc,'')}\n")

# ---------------- screen ----------------
cands = []
for key, ps in params.items():
    if not all(c in ps for c in CORE):
        continue
    voc = sorted(p for p in ps if is_voc(p))
    if len(voc) < 10:
        continue
    cands.append((key, len(voc), voc))

say("")
say("candidates (O3+NO+NO2+CO present AND >=10 VOC species):", len(cands))
say("   of which flagged PAMS:", sum(1 for k, _, _ in cands if pams[k]))

# ---------------- join sites ----------------
meta = {}
with open(SIT, newline="", encoding="utf-8", errors="replace") as fh:
    r = csv.DictReader(fh)
    for row in r:
        meta[(row["State Code"], row["County Code"], row["Site Number"])] = row

cands.sort(key=lambda t: (not pams[t[0]], -t[1], t[0]))
say("")
hdr = (f"{'st':>3}{'cty':>4}{'site':>5}{'VOC':>5}{'PAMS':>5}{'1stYr':>7} "
       f"{'last':<11} city / site / agency / lat,lon")
say(hdr); say("-" * 118)
for key, nvoc, voc in cands:
    m = meta.get(key, {})
    city = (m.get("City Name") or "").strip()
    nm = (m.get("Local Site Name") or "").strip()
    ag = (m.get("Owning Agency") or "").strip()
    st = (m.get("State Name") or "").strip()
    say(f"{key[0]:>3}{key[1]:>4}{key[2]:>5}{nvoc:>5}{('Y' if pams[key] else '-'):>5}"
        f"{(firsty.get(key,'?')):>7} {(lastd.get(key,'')):<11} "
        f"{city}, {st} | {nm} | {ag} | {m.get('Latitude','')},{m.get('Longitude','')}")

# PAMS candidates that also span our target years
say("")
say("=== PAMS candidates with data reaching 2016+ ===")
for key, nvoc, voc in cands:
    if not pams[key]:
        continue
    ld = lastd.get(key, "")
    if ld and ld >= "2016":
        m = meta.get(key, {})
        say(f"{key}  VOC={nvoc}  first={firsty.get(key,'?')} last={ld}  "
            f"{(m.get('City Name') or '').strip()}, {(m.get('State Name') or '').strip()}"
            f" | {(m.get('Local Site Name') or '').strip()}"
            f" | {m.get('Latitude','')},{m.get('Longitude','')}")

log.close()
print("wrote", OUT, "and", OUT2)
