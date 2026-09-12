"""Probe: how much leverage does biogenic isoprene actually have here?

Before spending 2**13 coalition runs on making ISOP a player, measure its size
in the two regimes we already care about.  The model already EMITS isoprene
(E_ISOP = EMIS_ISO * SDAY**1.5, daytime-peaked); it is simply not switchable.

Runs (each: r_sat and r_lim1, backend "co"):
    mask   0 (all windows on)            x EMIS_ISO in {0, 1, 3, 6}
    mask 240 (anthropogenic VOC all off) x EMIS_ISO in {0, 1}

The ISO sweep gives the leverage curve; the 2x2 gives the interaction between
anthropogenic and biogenic VOC.  Reported as last-24h mean and MDA8, the same
two metrics used everywhere else.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpp_box import run_scenario          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "g4runs", "co")
SCEN = {"r_sat": {}, "r_lim1": {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5}}
ISO_SWEEP = [0.0, 1.0, 3.0, 6.0]


def read_traj(path):
    rows = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            rows.append([float(x) for x in ln.split()])
    return rows


def mean_and_mda8(traj, t0=48.0):
    t = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    sel = [i for i in range(len(t)) if t[i] >= t0]
    mean = sum(o3[i] for i in sel) / len(sel)
    dt = t[1] - t[0]
    w = max(1, int(round(8.0 / dt)))
    mda8 = 0.0
    for i in range(w - 1, len(t)):
        if t[i] < t0:
            continue
        mda8 = max(mda8, sum(o3[i - w + 1:i + 1]) / w)
    return mean, mda8


def run(s, mask, iso):
    cfg = dict(SCEN[s])
    cfg["VOCW"] = [0.0 if (mask >> (4 + w)) & 1 else 1.0 for w in range(4)]
    cfg["NOXW"] = [0.0 if (mask >> w) & 1 else 1.0 for w in range(4)]
    cfg["COW"] = [1.0] * 4
    cfg["EMIS_ISO"] = iso
    name = "isop_%s_m%03d_i%02d" % (s, mask, int(iso * 10))
    run_scenario(cfg, name, backend="co")
    traj = read_traj(os.path.join(RUNS, name + ".traj"))
    return mean_and_mda8(traj)


def main():
    L = []
    L.append("Biogenic isoprene leverage (EMIS_ISO in ppb/h; baseline = 1.0)")
    L.append("")
    L.append("A) ISO sweep with ALL windows on (anthro emissions at baseline):")
    L.append("   %-8s %6s %9s %9s %9s %9s"
             % ("scenario", "ISO", "mean24", "MDA8", "dMean", "dMDA8"))
    for s in SCEN:
        base = None
        for iso in ISO_SWEEP:
            mean, mda8 = run(s, 0, iso)
            if base is None:
                base = (mean, mda8)
            L.append("   %-8s %6.1f %9.2f %9.2f %+9.2f %+9.2f"
                     % (s, iso, mean, mda8, mean - base[0], mda8 - base[1]))
        L.append("")
    L.append("B) 2x2 anthro-VOC(off) x biogenic-ISO(on/off):")
    L.append("   %-8s %10s %6s %9s %9s"
             % ("scenario", "anthroVOC", "ISO", "mean24", "MDA8"))
    for s in SCEN:
        for mask, lab in ((0, "on"), (240, "OFF")):
            for iso in (1.0, 0.0):
                mean, mda8 = run(s, mask, iso)
                L.append("   %-8s %10s %6s %9.2f %9.2f"
                         % (s, lab, "on" if iso else "OFF", mean, mda8))
        L.append("")
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "isop_probe.txt"), "w") as f:
        f.write(txt + "\n")


if __name__ == "__main__":
    main()
