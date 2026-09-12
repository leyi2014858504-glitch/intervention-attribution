"""Scale stress test, part 2: read the n=24 results and calibrate against the
exact n=12 table.

Two halves:
  A) CALIBRATION at n=12 (oracle = the exact 4096 table in co_v.csv).  For the
     same estimator and budgets we can see BOTH the true error and the
     seed-to-seed spread, so we get the spread -> error ratio.
  B) n=24 (3 sectors x 8 three-hour windows).  No oracle exists (2**24 = 16.7M
     coalitions), so we report the chemistRY cost actually spent (distinct
     evaluations) and the seed spread, then convert with the n=12 ratio.

Estimator = permutation + antithetic (the strongest arm at n=8).
Error/spread are normalised by max_i |phi_i|.
"""
import csv
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
N24, NWIN = 24, 8
N12 = 12
SCEN = ["r_sat", "r_lim1"]
BUDGETS = [64, 128, 256]
NSEED = 8


def gen(rng, budget, n):
    """Exactly the generator used by gen_scale_masks (same RNG calls)."""
    seen = set()
    orders = []
    while len(seen) < budget:
        order = list(range(n))
        rng.shuffle(order)
        for seq in (order, order[::-1]):
            mask = 0
            seen.add(0)
            for i in seq:
                mask |= 1 << i
                seen.add(mask)
            seen.add((1 << n) - 1)
            orders.append(seq)
            if len(seen) >= budget:
                break
    return orders, seen


def phi_perm(vd, orders, n):
    s = [0.0] * n
    c = [0] * n
    for seq in orders:
        mask = 0
        for i in seq:
            prev = mask
            mask |= 1 << i
            s[i] += vd[mask] - vd[prev]
            c[i] += 1
    return [s[i] / c[i] for i in range(n)]


def mobius(v, n):
    out = [0.0] * (1 << n)
    for mask in range(1 << n):
        s, sub = 0.0, mask
        while True:
            s += (1.0 if bin(mask ^ sub).count("1") % 2 == 0 else -1.0) * v[sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out[mask] = s
    return out


def load_csv(fn, keep=None, mdim=None):
    tab = {}
    with open(os.path.join(HERE, fn)) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if keep and s not in keep:
                continue
            if row["v_ppb"] in ("", "NaN"):
                continue
            if mdim and int(ms) >= mdim:
                continue
            tab.setdefault(s, {})[int(ms)] = float(row["v_ppb"])
    return tab


def mean_sd(xs):
    m = sum(xs) / len(xs)
    if len(xs) < 2:
        return m, 0.0
    sd = math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))
    return m, sd


def main():
    L = []
    add = L.append

    # ---------------- A) calibration at n=12 -------------------------
    t12 = load_csv("co_v.csv", keep=set(SCEN), mdim=4096)
    add("=" * 78)
    add("A) CALIBRATION at n=12   (oracle = exact 4096-point table)")
    add("=" * 78)
    add("  %-8s %6s %8s %10s %10s %9s"
        % ("scen", "B", "runs", "err(true)", "spread", "err/spread"))
    ratio = {}
    for s in SCEN:
        v12 = [t12[s][m] for m in range(4096)]
        v0 = v12[0]
        mo = mobius([x - v0 for x in v12], N12)
        true = [sum(mo[A] / bin(A).count("1")
                    for A in range(1, 4096) if (A >> i) & 1) for i in range(N12)]
        scale = max(abs(x) for x in true)
        for B in BUDGETS:
            ests = []
            for sd in range(NSEED):
                rng = random.Random(1000 * B + sd)
                orders, seen = gen(rng, B, N12)
                ests.append(phi_perm({m: v12[m] for m in seen}, orders, N12))
            mu = [mean_sd([e[i] for e in ests])[0] for i in range(N12)]
            sdi = [mean_sd([e[i] for e in ests])[1] for i in range(N12)]
            err = max(abs(mu[i] - true[i]) for i in range(N12)) / scale
            sprd = max(sdi) / scale
            ratio[(s, B)] = err / sprd if sprd > 0 else float("nan")
            add("  %-8s %6d %8d %9.2f%% %9.2f%% %9.2f"
                % (s, B, B, 100 * err, 100 * sprd, err / sprd))

    # ---------------- B) n=24 ----------------------------------------
    t24 = load_csv("scale_v.csv", keep=set(SCEN))
    add("")
    add("=" * 78)
    add("B) n=24   (3 sectors x 8 three-hour windows, 16.7M coalitions)")
    add("=" * 78)
    for s in SCEN:
        if s not in t24:
            add("  %s: no data" % s)
            continue
        vd = t24[s]
        add("  %s: %d coalition values collected" % (s, len(vd)))
        add("  %-8s %8s %10s %12s %14s"
            % ("B", "runs/seed", "spread", "est.err (map)", "sector(8-run)"))
        for B in BUDGETS:
            ests, nseen = [], []
            for sd in range(NSEED):
                rng = random.Random(1000 * B + sd)
                orders, seen = gen(rng, B, N24)
                missing = seen - set(vd)
                if missing:
                    add("    B=%d seed=%d: %d coalitions missing"
                        % (B, sd, len(missing)))
                    continue
                ests.append(phi_perm({m: vd[m] for m in seen}, orders, N24))
                nseen.append(len(seen))
            if not ests:
                continue
            mu = [mean_sd([e[i] for e in ests])[0] for i in range(N24)]
            sdi = [mean_sd([e[i] for e in ests])[1] for i in range(N24)]
            sc = max(abs(x) for x in mu)
            sprd = max(sdi) / sc
            r = ratio.get((s, B), float("nan"))
            add("  %-8d %8d %9.2f%% %11.2f%%"
                % (B, sum(nseen) // len(nseen), 100 * sprd, 100 * sprd * r))

            # window ranking at the largest budget
            if B == max(BUDGETS):
                add("    top-6 players at n=24 (mean over seeds, ppb):")
                top = sorted(range(N24), key=lambda i: -abs(mu[i]))[:6]
                for i in top:
                    role = "NOx" if i < 8 else ("VOC" if i < 16 else "CO")
                    w = i % 8
                    add("      %-4s 3h-window #%d [%02d-%02d h)  phi = %+7.2f"
                        % (role, w + 1, 3 * w, 3 * w + 3, mu[i]))
                # within-6h split: do the two 3h halves of the SAME 6h block
                # carry the same weight?  The 4-window game cannot see this.
                add("    within-6h split (the two 3h halves of each 6h block):")
                for role, base in (("NOx", 0), ("VOC", 8), ("CO", 16)):
                    for b in range(4):
                        p, q = base + 2 * b, base + 2 * b + 1
                        hi = "half1" if abs(mu[p]) >= abs(mu[q]) else "half2"
                        add("      %-4s block W%d [%02d-%02d h): half1 %+6.2f  "
                            "half2 %+6.2f   |ratio| = %4.2f  (%s dominates)"
                            % (role, b + 1, 6 * b, 6 * b + 6, mu[p], mu[q],
                               (max(abs(mu[p]), abs(mu[q]))
                                / max(1e-9, min(abs(mu[p]), abs(mu[q])))),
                               hi))

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "scale_result.txt"), "w") as f:
        f.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
