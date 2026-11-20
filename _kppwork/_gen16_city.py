#!/usr/bin/env python3
"""Generate all 2**16 namelists for a CITY-ANCHORED n=16 game.

Same player layout as gen_n16.py (bits 0-7 = NOxW(8), 8-15 = VOCW(8), COW all
1.0), but the config base is a real-anchored slice: the 8 knobs from
results_fit_<label>.txt (8-knob protocol), backgrounds and TMIN/TMAX exactly as
_dallas_lvg.py assembles them.  This lets the SAME sweep_n16.sh enumerate the
real slices with zero chemistry-layer changes -> real-city full spectra,
top-dividends and truncation curves replace the synthetic-only caveat.

Usage: python _gen16_city.py dallas_2016w deerpark_2017s
Writes /root/nml16/<label>/m<05d>.nml  (+ vgate_<label>.nml for the smoke gate)
"""
import os
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
from kpp_box import nml_text  # noqa: E402

D = os.environ.get("ATMOS_DIR", "/mnt/d/atmos")
NWIN = 8
N16 = 16
ROOT = "/root/nml16"
KNOBS = ("EMIS_NOX", "EMIS_VOC", "EMIS_ISO", "EMIS_CO",
         "KDIL_DAY", "KDIL_NIGHT", "BG_O3", "BG_CO")


def read_anchor(label):
    path = os.path.join(D, "results_fit_%s.txt" % label)
    knobs = {}
    with open(path) as fh:
        for line in fh:
            if "=" not in line:
                continue
            k, _, v = line.partition("=")
            k = k.strip()
            if k in KNOBS:
                knobs[k] = float(v)
    missing = [k for k in KNOBS if k not in knobs]
    if missing:
        raise SystemExit("%s is not an 8-knob anchor file; missing %s" % (path, missing))
    return knobs


def read_met_tminmax(label):
    path = os.path.join(D, "met_%s.csv" % label)
    T = []
    with open(path) as fh:
        for r in fh:
            p = r.strip().split(",")
            if p[0] == "time":
                continue
            try:
                T.append(float(p[1]))
            except (ValueError, IndexError):
                pass
    return min(T) + 273.15, max(T) + 273.15


def cfg_of(mask):
    """mask -> RUNCFG overrides.  Set bit = player OFF."""
    nox = [0.0 if (mask >> w) & 1 else 1.0 for w in range(NWIN)]
    voc = [0.0 if (mask >> (NWIN + w)) & 1 else 1.0 for w in range(NWIN)]
    return {"NOXW": nox, "VOCW": voc, "COW": [1.0] * NWIN}


def build_base(label):
    anchor = read_anchor(label)
    tmin, tmax = read_met_tminmax(label)
    base = dict(anchor)
    # backgrounds the fit's cost function sets (see _dallas_lvg.py / _deerpark_lvg.py)
    base.update(BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
                TMIN_K=tmin, TMAX_K=tmax)
    return base, anchor


def main(labels):
    for label in labels:
        base, anchor = build_base(label)
        d = os.path.join(ROOT, label)
        os.makedirs(d, exist_ok=True)
        for f in os.listdir(d):
            if f.endswith(".nml"):
                os.remove(os.path.join(d, f))
        for mask in range(1 << N16):
            cfg = dict(base)
            cfg.update(cfg_of(mask))
            with open(os.path.join(d, "m%05d.nml" % mask), "w", newline="\n") as fh:
                fh.write(nml_text(cfg, "co"))
        # smoke gate target: all players ON must reproduce the lvg base value
        cfg = dict(base)
        cfg.update(cfg_of(0))
        with open(os.path.join(ROOT, "vgate_%s.nml" % label), "w", newline="\n") as fh:
            fh.write(nml_text(cfg, "co"))
        print("%-18s %d namelists -> %s   anchor=%s" % (label, 1 << N16, d, anchor))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["dallas_2016w", "deerpark_2017s"]))
