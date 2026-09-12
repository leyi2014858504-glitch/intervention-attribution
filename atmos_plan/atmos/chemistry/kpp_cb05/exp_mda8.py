"""MDA8 re-verification of the budget-constrained optima.

The 4096-run sweep stored only the last-24h MEAN O3 (traj.dat was overwritten
each run), so the budget-constrained recipes cannot be re-scored from
scen_v.csv.  Policy statements are normally made in MDA8, not in the daily
mean, so re-run the handful of recipes that matter and score both metrics on
the same trajectories.

Recipes (bit w = NOx window w+1 OFF, bit 4+w = VOC window w+1 OFF):
    0    none
    3    NOx@W1+W2
    64   VOC@W3
    67   NOx@W1+W2 | VOC@W3      <- constrained optimum at (B_N,B_V)=(1.0,0.35)
    240  VOC@all                 <- unconstrained optimum
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kpp_box import run_scenario            # noqa: E402
from gen_scen_sweep import R_SET, cfg_of    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "g4runs", "dummy")

SCENS = ["r_sat", "r_lim1", "r_lim2"]
RECIPES = [(0, "none"), (3, "NOx@W1+W2"), (64, "VOC@W3"),
           (67, "NOx@W1+W2|VOC@W3"), (240, "VOC@all")]


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


def main():
    print("scenario  recipe                    mean24    MDA8   d_mean  d_MDA8")
    print("-" * 74)
    for s in SCENS:
        base = R_SET[s]
        ref = None
        for mask, lab in RECIPES:
            cfg = dict(base)
            cfg.update(cfg_of(mask))
            name = "mda8_%s_m%03d" % (s, mask)
            m = run_scenario(cfg, name, backend="dummy")
            traj = read_traj(os.path.join(RUNS, name + ".traj"))
            mean, mda8 = mean_and_mda8(traj)
            if ref is None:
                ref = (mean, mda8)
            print("%-9s %-24s %7.2f %7.2f  %+7.2f %+7.2f"
                  % (s, lab, mean, mda8, mean - ref[0], mda8 - ref[1]))
            if m.get("nfail", 0):
                print("   WARNING nfail=%d" % m["nfail"])
        print("-" * 74)


if __name__ == "__main__":
    main()
