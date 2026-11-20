#!/usr/bin/env python3
"""One-knob variant: KDIL_DAY 1500 -> 2300 s (city-fitted value).

Usage: python3 _make_tau2300.py <src tag> <dst tag>   (default r_sat)
Copies every namelist and changes only that line.
"""
import os
import sys

src_tag = sys.argv[1] if len(sys.argv) > 1 else "r_sat"
dst_tag = sys.argv[2] if len(sys.argv) > 2 else src_tag + "_tau2300"
SRC = "/root/nml16/" + src_tag
DST = "/root/nml16/" + dst_tag
OLD = "KDIL_DAY = 1500.0"
NEW = "KDIL_DAY = 2300.0"

os.makedirs(DST, exist_ok=True)
n = 0
for f in sorted(os.listdir(SRC)):
    if not f.endswith(".nml"):
        continue
    with open(os.path.join(SRC, f)) as fh:
        s = fh.read()
    if OLD not in s:
        raise SystemExit("pattern not found in %s" % f)
    with open(os.path.join(DST, f), "w", newline="\n") as fh:
        fh.write(s.replace(OLD, NEW))
    n += 1
print("wrote %d namelists to %s" % (n, DST))
with open(os.path.join(DST, "m00000.nml")) as fh:
    for line in fh:
        if "KDIL_DAY" in line or "NOXW" in line:
            print("check:", line.strip())
