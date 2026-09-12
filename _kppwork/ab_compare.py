"""A/B: dummy-reaction forcing vs operator-splitting forcing.

Both backends run the same scenario set (see exp_g4_kpp.py); this compares the
MDA8 grid cell by cell, so the difference measures exactly what the switch from
"forcing applied after each INTEGRATE step" to "forcing inside the ODE" costs.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(backend):
    d = {}
    with open(os.path.join(HERE, "g4_grid_%s.csv" % backend)) as f:
        for row in csv.reader(l for l in f if not l.startswith("#")):
            if not row:
                continue
            metric, nf, vf, val = row[0], float(row[1]), float(row[2]), float(row[3])
            d[(metric, nf, vf)] = val
    return d


def main():
    a, b = load("dummy"), load("split")
    keys = sorted(set(a) & set(b))
    grid = [k for k in keys if k[0] == "ekma"]
    deltas = [(k, a[k], b[k]) for k in grid]
    ds = [x[1] - x[2] for x in deltas]
    rel = [abs(x[1] - x[2]) / x[2] * 100.0 for x in deltas]
    absd = [abs(d) for d in ds]

    print("A/B over the %d EKMA cells (dummy - split), MDA8 ppb:" % len(grid))
    print("  max |delta| = %.3f ppb   mean |delta| = %.3f ppb"
          % (max(absd), sum(absd) / len(absd)))
    print("  max rel     = %.3f %%     mean rel    = %.3f %%"
          % (max(rel), sum(rel) / len(rel)))
    print("  signed delta range: %+.3f .. %+.3f ppb" % (min(ds), max(ds)))
    print("")
    for k in keys:
        if k[0] != "ekma":
            print("  %-11s dummy=%8.3f  split=%8.3f  delta=%+.3f"
                  % (k[0], a[k], b[k], a[k] - b[k]))


if __name__ == "__main__":
    main()
