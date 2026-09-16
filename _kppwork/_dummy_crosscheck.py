"""Validate mobius_v_dummy.csv (n=8 oracle behind both benchmarks) against
the fresh n16_m8 table via the 6 h sub-lattice embedding."""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))

v16 = {}
with open(os.path.join(HERE, "n16_v_r_sat_m8.csv")) as fh:
    for r in csv.DictReader(fh):
        v16[int(r["name"][1:])] = float(r["v_ppb"])

def embed(m8):
    m16 = 0
    for j in range(4):
        if (m8 >> j) & 1:
            m16 |= 1 << (2 * j) | 1 << (2 * j + 1)
        if (m8 >> (4 + j)) & 1:
            m16 |= 1 << (8 + 2 * j) | 1 << (8 + 2 * j + 1)
    return m16

old = {}
with open(os.path.join(HERE, "mobius_v_dummy.csv")) as fh:
    for r in csv.DictReader(fh):
        nm = r["name"]
        if nm.startswith("mob_"):
            old[int(nm[4:])] = float(r["v_ppb"])

base16 = v16[0]
based = old[0]          # dummy table stores RAW values too; compare game-relative
bad = 0
for m8, v in sorted(old.items()):
    ref = v16[embed(m8)]
    if abs((ref - base16) - (v - based)) > 1e-6:
        bad += 1
        if bad <= 8:
            print("  m8=%-4d dummy=%10.4f n16-derived=%10.4f" % (m8, v, ref - base16))
print("mobius_v_dummy: %d rows checked, %d disagree" % (len(old), bad))
