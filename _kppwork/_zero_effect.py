#!/usr/bin/env python3
"""Is the near-zero local derivative in Scan A a RIDGE (physics) or noise?

Scan A reported |ratio - 1| = 42 at point n075v100, i.e. a window where the
local tangent is ~0 while the finite switch-off effect is not.  Two readings:

  (a) RIDGE: dO3/dNOx = 0 is the DEFINITION of the EKMA ridge line.  On it a
      local method reports "no effect" for a lever that actually matters --
      worse than a sign error, because "zero" reads as safe.
  (b) NOISE: if D1 is pure differencing dust, ratio is meaningless.

Discriminator: compare |D1| against the numerical noise floor of the runs.
The box is deterministic (the hddm cross-check matched lvg bit-for-bit at
4.7e-07), so the floor here is the truncation error of the central difference,
O(delta^2 * f'''), estimated from the spread of D2 across windows -- NOT
machine epsilon.  If |D1| is orders of magnitude above that floor the zero is
physical.

Writes results_zero_effect.txt
"""
import glob
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "results_zero_effect.txt")
DELTA = 0.05
SECTORS = ["NOXW", "VOCW", "COW"]
NW = 8
GRID = ["n050v050", "n050v075", "n050v100", "n050v150", "n050v200",
        "n075v050", "n075v075", "n075v100", "n075v150", "n075v200",
        "n100v050", "n100v075", "n100v100", "n100v150", "n100v200",
        "n150v050", "n150v075", "n150v100", "n150v150", "n150v200",
        "n200v050", "n200v075", "n200v100", "n200v150", "n200v200"]


def load():
    v = {}
    for f in sorted(glob.glob(os.path.join(HERE, "flip_v_*.csv"))):
        for ln in open(f):
            p = ln.strip().split(",")
            if len(p) < 3 or p[1] in ("NaN", "v_ppb"):
                continue
            try:
                v[p[0]] = (float(p[1]), float(p[2]))
            except ValueError:
                pass
    return v


def main():
    v = load()
    L = []
    add = L.append
    add("Is the near-zero local derivative a RIDGE (physics) or noise?")
    add("protocol: 3h windows, delta=0.05, same runs as results_scan_flips.txt")
    add("")

    rows = []
    for pt in GRID:
        if pt + "_base" not in v:
            continue
        base = v[pt + "_base"]
        for sec in SECTORS:
            for i in range(NW):
                try:
                    vp = v["%s_%s_p%d" % (pt, sec, i)]
                    vm = v["%s_%s_m%d" % (pt, sec, i)]
                    vo = v["%s_%s_off%d" % (pt, sec, i)]
                except KeyError:
                    continue
                for tag, j in (("8", 1), ("48", 0)):
                    d1 = (vp[j] - vm[j]) / (2 * DELTA)
                    d2 = (vp[j] - 2 * base[j] + vm[j]) / (DELTA ** 2)
                    dfin = vo[j] - base[j]
                    if abs(d1) < 1e-12:
                        continue
                    rows.append(dict(pt=pt, sec=sec, win="%d-%d" % (3 * i, 3 * i + 3),
                                     tag=tag, d1=d1, d2=d2, dfin=dfin,
                                     ratio=dfin / (-d1),
                                     scale=abs(base[j])))
    # noise floor: the runs are deterministic, so the floor is the central-
    # difference truncation error.  Estimate it as delta^2/6 * |D2| * anything
    # else is smooth -- a conservative proxy is delta^2 * |D2|.
    for tag, mname in (("8", "MDA8"), ("48", "d3_mean")):
        rs = [r for r in rows if r["tag"] == tag]
        floors = [DELTA ** 2 * abs(r["d2"]) for r in rs]
        floor = sorted(floors)[len(floors) // 2]      # median proxy
        add("=" * 88)
        add("METRIC = %s    central-difference truncation floor ~ delta^2*|D2| = %.2e"
            % (mname, floor))
        add("=" * 88)
        add("%-12s %-6s %-8s %10s %10s %10s %10s %9s" %
            ("point", "sector", "window", "D1", "D2", "D_finite", "ratio", "|D1|/floor"))
        hot = [r for r in rs if abs(r["ratio"]) > 5]
        for r in sorted(hot, key=lambda r: -abs(r["ratio"])):
            add("%-12s %-6s %-8s %10.4f %10.4f %10.4f %10.2f %9.0f" %
                (r["pt"], r["sec"], r["win"], r["d1"], r["d2"], r["dfin"],
                 r["ratio"], abs(r["d1"]) / floor if floor else float("inf")))
        add("")
        add("windows with |ratio| > 5 : %d of %d" % (len(hot), len(rs)))
        if hot:
            big = [r for r in hot if abs(r["dfin"]) > 0.5]
            add("  of those, |D_finite| > 0.5 ppb (i.e. the lever MATTERS): %d" % len(big))
            for r in big:
                add("    %s %s@%s : local says %.3f ppb (tangent), finite truth %.3f ppb"
                    % (r["pt"], r["sec"], r["win"], -r["d1"], r["dfin"]))
            mind1 = min(abs(r["d1"]) for r in hot)
            add("  smallest |D1| among them: %.4f ppb  (= %.0f x the truncation floor)"
                % (mind1, mind1 / floor if floor else float("inf")))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(os.path.abspath(OUT), "w") as fh:
        fh.write(txt + "\n")
    print("wrote", os.path.abspath(OUT))


if __name__ == "__main__":
    main()
