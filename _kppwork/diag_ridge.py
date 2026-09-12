"""Does the R-set actually cross the EKMA ridge?  Diagnose the sign of the
single-window NOx switches in each scenario.

If a scenario is truly on the NOx-limited side, turning OFF one window's NOx
should DECREASE O3 (m({NOx@Ww}) < 0).  If every main effect is positive, the
whole R-set sits on the NOx-saturated side and the "regime walk" is a walk in
parameter space, not across the ridge.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8


def pc(m):
    return bin(m).count("1")


tab = {}
with open(os.path.join(HERE, "scen_v.csv")) as f:
    for row in csv.DictReader(f):
        nm = row["name"]
        if "_m" not in nm:
            continue
        scen, ms = nm.rsplit("_m", 1)
        tab.setdefault(scen, {})[int(ms)] = float(row["v_ppb"])

RS = ["r_sat", "r_sat_hi", "r_trans", "r_lim1", "r_lim2", "r_wkend"]
print("single-window main effects  m({i}) = v(off i) - v(none)   [ppb]")
print("  positive = turning that window OFF *raises* O3 (NOx-saturated sign)")
print("")
hdr = "  %-9s %4s | " % ("scenario", "v0") + " ".join(
    "%9s" % ("NOx@W%d" % (w + 1)) for w in range(4)) + " | " + \
    " ".join("%9s" % ("VOC@W%d" % (w + 1)) for w in range(4))
print(hdr)
for s in RS:
    v0 = tab[s][0]
    row = ["%9.2f" % (tab[s][1 << i] - v0) for i in range(8)]
    print("  %-9s %4.1f | %s | %s" % (s, v0, " ".join(row[:4]), " ".join(row[4:])))

print("")
print("most-NOx-off point within the NOx windows (NOXW = 0000):")
for s in RS:
    v0 = tab[s][0]
    print("  %-9s v(all 4 NOx off) - v0 = %+8.2f ppb" % (s, tab[s][0x0F] - v0))

# the strongest NOx-reduction point reachable while keeping VOC untouched
print("")
print("best (lowest-O3) NOx-only recipe vs baseline (VOC untouched):")
for s in RS:
    v0 = tab[s][0]
    best = min((tab[s][m], m) for m in range(1, 16))
    print("  %-9s best NOx-only = %s  delta = %+8.2f ppb"
          % (s, format(best[1], "04b"), best[0] - v0))
