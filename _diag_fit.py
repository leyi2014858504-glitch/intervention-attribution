#!/usr/bin/env python3
"""
Two questions the scan could not answer:

1. Is the 3-day run (TEND_S=259200) actually at a diurnal steady state? PAR
   starts near 453 ppb in the default initial condition while the observation is
   ~49, so if PAR is still decaying at t=48..72h the whole fit is fitting a
   transient, not a cycle.
2. Why is CO_err ~ 4.0 in every configuration? Print the levels.

Run the same configuration for 3 days and for 10 days and compare the LAST day.
"""
import sys, csv, statistics
sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box

D = "/mnt/d/atmos"
label = "dallas_2016w"

obs = {k: {h: [] for h in range(24)} for k in ("O3", "NOx", "CO", "PAR")}
with open(f"{D}/obs_{label}_cb05.csv", newline="") as fh:
    for r in csv.DictReader(fh):
        h = int(r["datetime"][11:13])
        for k in obs:
            v = r.get(k)
            if v not in (None, ""):
                obs[k][h].append(float(v))
obsm = {k: {h: (sum(v) / len(v) if v else float("nan")) for h, v in obs[k].items()} for k in obs}

base = dict(EMIS_NOX=6.0, EMIS_VOC=1.0, EMIS_ISO=0.3, EMIS_CO=60.0,
            KDIL_DAY=22880.0, KDIL_NIGHT=86400.0,
            BG_O3=14.0, BG_CO=200.0, BG_NO2=2.0, BG_NO=0.5,
            BG_VOC=20.0, BG_ISO=0.2, TMIN_K=271.2, TMAX_K=294.6)

runs = {}
for tag, tend in (("3d", 259200.0), ("10d", 864000.0)):
    cfg = dict(base); cfg["TEND_S"] = tend
    m = kpp_box.run_scenario(cfg, f"diag_{tag}", "anchor")
    traj = m["traj"]
    tmax = traj[-1][0]
    last24 = tmax - 24.0
    # last-day composite, computed here so the 10-day run is not truncated by
    # diurnal()'s fixed 48 h spin-up
    acc = {c: {h: [] for h in range(24)} for c in ("o3", "nox", "co", "par")}
    for r in traj:
        if r[0] < last24:
            continue
        h = int(round(r[0])) % 24
        for i, c in enumerate(("o3", "nox", "co", "par")):
            acc[c][h].append(r[1 if c == "o3" else (2 if c == "nox" else (7 if c == "co" else 10))])
    runs[tag] = {c: {h: statistics.mean(v) if v else float("nan") for h, v in acc[c].items()}
                 for c in acc}
    runs[tag]["_tmax"] = tmax
    print(f"{tag}: tmax={tmax:.1f}h")

print("\nhour    " + " ".join(f"{h:>8d}" for h in range(24)))
for c, o in (("o3", "O3"), ("nox", "NOx"), ("co", "CO"), ("par", "PAR")):
    for tag in ("3d", "10d"):
        p = runs[tag][c]
        print(f"{tag:>4} {c:<4}" + " ".join(f"{p[h]:8.2f}" for h in range(24)))
    print(f"     obs {o:<4}" + " ".join(f"{obsm[o][h]:8.2f}" for h in range(24)))
    print()

for c in ("o3", "nox", "co", "par"):
    m3 = statistics.mean([runs['3d'][c][h] for h in range(24)])
    m10 = statistics.mean([runs['10d'][c][h] for h in range(24)])
    print(f"{c:<4} daily mean: 3d={m3:8.2f}   10d={m10:8.2f}   drift={100*(m3-m10)/max(m10,1e-9):6.1f}%")
