"""Budget-constrained (EKMA-style) time-resolved control allocation.

WHY THIS EXISTS.  Route B's unconstrained optimum is regime-INVARIANT
("switch all VOC off" in all six regimes), so the regime dependence shows up
only in the *explanation*, not in the *policy*.  This script adds the one
constraint that is standard practice, and checks whether the policy itself
then becomes time-resolved and regime-dependent.

BACKING (why a reduction budget, not an ad-hoc "at most k windows").
  1. The EKMA control-strategy diagram -- the canonical O3 control framework
     -- has *percent VOC reduction* and *percent NOx reduction* as its two
     axes (Dodge 1977; EPA EKMA guideline).  A control strategy is a point
     (R_V, R_N) in that plane, not a list of switches.
  2. Attainment planning / State Implementation Plans impose emission
     *budgets* per precursor; cost-effectiveness studies allocate reductions
     under such a budget.
  So the constraint currency is a fraction of the daily anthropogenic
  emission, and the constraint is per precursor -- exactly the two EKMA axes.
  A joint "total mass" budget would need an arbitrary NOx/VOC exchange rate,
  so it is NOT used.

CURRENCY IS EXACT, NO NEW RUNS.  The driver sets
    ENOX = EMIS_NOX*NOX_FACTOR*NOXW(IW)*tmpd(HH)
(mm/hh of cb05d_Main.f90), so the daily emission removed by switching window
w OFF is proportional to  int_w tmpd dh.  The constant factors cancel in the
fraction, so the per-window emission shares depend only on the weekday /
weekend diurnal profile -- computable here to machine precision.

THE POINT.  For any budget strictly below 100 %, the only time-UNIFORM VOC
policies are "off everywhere" (100 %) and "on everywhere" (0 %).  Neither is
feasible, so **every feasible policy with a reduction is necessarily
time-resolved**.  The open question is then *which* windows carry the cut --
and that is what this script measures across regimes.
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8
MASK_FULL = (1 << NPLAY) - 1

R_SET = ["r_sat", "r_sat_hi", "r_trans", "r_lim1", "r_lim2", "r_wkend"]
REF = ["u0_base"] + R_SET            # u0_base shown as the weekday reference
WKEND_OF = {"r_wkend": True}
BASELINE_V = {"u0_base": 10.70}      # filled from data below if present


# ---------------------------------------------------------------- profiles
def gposs(mu, sd, hh):
    x = hh - mu
    if x > 12.0:
        x -= 24.0
    if x < -12.0:
        x += 24.0
    return pow(2.718281828459045, -(x / sd) ** 2)


def tmpd(hh, weekend):
    if weekend:
        ga = gposs(11.0, 3.0, hh) + gposs(17.0, 3.0, hh)
        return 0.55 + 0.70 * ga
    ga = gposs(8.0, 1.6, hh) + gposs(18.5, 2.0, hh)
    return 0.45 + 1.30 * ga


def window_fracs(weekend, nsub=20000):
    """f[w] = share of the daily emission emitted inside window w (1..4)."""
    mass = [0.0] * 4
    dt = 24.0 / nsub
    for k in range(nsub):
        hh = (k + 0.5) * dt
        w = min(4, int(hh / 6.0) + 1)
        mass[w - 1] += tmpd(hh, weekend) * dt
    tot = sum(mass)
    return [m / tot for m in mass]


# ---------------------------------------------------------------- truth
def load():
    tab = {}
    with open(os.path.join(HERE, "scen_v.csv")) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            scen, ms = nm.rsplit("_m", 1)
            v = row["v_ppb"]
            if v in ("", "NaN"):
                continue
            tab.setdefault(scen, {})[int(ms)] = float(v)
    return tab


def redux(mask, f):
    """(R_NOx, R_VOC) = daily emission fraction switched off."""
    rn = sum(f[w] for w in range(4) if (mask >> w) & 1)
    rv = sum(f[w] for w in range(4) if (mask >> (4 + w)) & 1)
    return rn, rv


def uniform_species(mask):
    """True if, per species, all four windows are treated identically."""
    for base in (0, 4):
        bits = [(mask >> (base + w)) & 1 for w in range(4)]
        if len(set(bits)) != 1:
            return False
    return True


def wnames(mask, role):
    """e.g. 'NOx@W2+W4' for the windows switched OFF."""
    base = 4 if role == "VOC" else 0
    ws = [w + 1 for w in range(4) if (mask >> (base + w)) & 1]
    if not ws:
        return "%s@none" % role
    if len(ws) == 4:
        return "%s@all" % role
    return "%s@%s" % (role, "+".join("W%d" % w for w in ws))


def solve(tab_s, f, bn, bv):
    """min v over feasible masks; returns (mask, v, rn, rv) or None."""
    best = None
    for m, v in tab_s.items():
        rn, rv = redux(m, f)
        if rn <= bn + 1e-9 and rv <= bv + 1e-9:
            if best is None or v < best[1]:
                best = (m, v, rn, rv)
    return best


def main():
    tab = load()
    FW = {"wk": window_fracs(False), "we": window_fracs(True)}

    print("=" * 78)
    print("window emission shares  f_w  (fraction of the daily emission)")
    print("=" * 78)
    print("  profile   " + "  ".join("%8s" % ("W%d" % (w + 1)) for w in range(4)))
    for k, lab in (("wk", "weekday"), ("we", "weekend")):
        print("  %-9s " % lab + "  ".join("%8.4f" % x for x in FW[k]))

    for scen in REF:
        if scen not in tab:
            continue
        f = FW["we"] if WKEND_OF.get(scen) else FW["wk"]
        t = tab[scen]
        v0 = t[0]
        print("")
        print("=" * 78)
        print("scenario %s   baseline v(none) = %.2f ppb" % (scen, v0))
        print("=" * 78)
        # unconstrained reference
        um = min(t, key=lambda m: t[m])
        unr = redux(um, f)
        print("  UNCONSTRAINED best = %-28s v=%.2f (%.2f ppb)  R_N=%.2f R_V=%.2f"
              % (wnames(um, "NOx") + " | " + wnames(um, "VOC"), t[um],
                 t[um] - v0, unr[0], unr[1]))

        print("  %-18s %-30s %7s %8s %6s %6s %s"
              % ("budget(N,V)", "optimal recipe", "v", "delta", "R_N", "R_V",
                 "uniform?"))
        for bn, bv in [(1.0, 0.50), (1.0, 0.35), (1.0, 0.25), (1.0, 0.15),
                       (0.50, 1.0), (0.25, 1.0),
                       (0.35, 0.35), (0.25, 0.25)]:
            r = solve(t, f, bn, bv)
            if r is None:
                print("  %-18s %-30s %7s" % ("%.2f/%.2f" % (bn, bv),
                                             "(no feasible recipe)", "-"))
                continue
            m, v, rn, rv = r
            recipe = wnames(m, "NOx") + " | " + wnames(m, "VOC")
            print("  %-18s %-30s %7.2f %+8.2f %6.2f %6.2f %s"
                  % ("%.2f/%.2f" % (bn, bv), recipe, v, v - v0, rn, rv,
                     "yes" if uniform_species(m) else "NO -> time-resolved"))

        # is any fully time-uniform recipe feasible at a tight budget?
        feas = [m for m in t if uniform_species(m)
                and redux(m, f)[0] <= 0.50 + 1e-9
                and redux(m, f)[1] <= 0.35 + 1e-9]
        print("  time-uniform recipes feasible at (0.50,0.35): %d  %s"
              % (len(feas), sorted(wnames(m, "NOx") + "|" + wnames(m, "VOC")
                                   for m in feas)))

    # ---------------- regime comparison of the *chosen windows* ------------
    print("")
    print("=" * 78)
    print("REGIME COMPARISON: which windows carry the cut at (B_N,B_V)=(1.0,0.35)")
    print("(pure-VOC budget: NOx unconstrained, VOC limited to 35% of daily)")
    print("=" * 78)
    print("  %-9s %-8s %-26s %-26s %8s" %
          ("scenario", "v0", "VOC windows cut", "NOx windows cut", "delta"))
    for scen in REF:
        if scen not in tab:
            continue
        f = FW["we"] if WKEND_OF.get(scen) else FW["wk"]
        t = tab[scen]
        r = solve(t, f, 1.0, 0.35)
        if r is None:
            continue
        m, v, rn, rv = r
        print("  %-9s %8.2f %-26s %-26s %+8.2f" %
              (scen, t[0], wnames(m, "VOC"), wnames(m, "NOx"), v - t[0]))

    print("")
    print("VERDICT: a budget strictly below 100% admits no time-uniform policy,")
    print("so the policy is forced to be time-resolved; see whether the chosen")
    print("window set moves across regimes above.")


if __name__ == "__main__":
    main()
