"""Exact Shapley / Owen attribution oracle for the n=8 emission-switch game.

Ground truth comes from the 256 fully enumerated configurations of
exp_mobius.py (mobius_v_dummy.csv), so no approximation and no extra chemistry
is involved: this is the reference against which every coalition sampler is
scored.

GAME
  8 players, index i<4 -> NOx emission off in window i+1,
             index i>=4 -> anthropogenic-VOC emission off in window i-3
  value  v(S) = mean O3 over the last 24 h of the 72 h episode, minus the same
  metric at S = {} (so v({}) = 0);  v is therefore a game in difference form.

TARGETS
  Shapley   phi_i = sum_{A ∋ i} m(A)/|A|                (m = Moebius spectrum)
            equivalently the textbook weighted-average-of-marginals formula;
            both are computed here and cross-checked.
  Owen      psi_i for the coalition structure B = {{W1 NOx, W1 VOC}, ..., }
            i.e. players are grouped BY TIME WINDOW (4 groups of 2), which is
            the "time-resolved" refinement: psi gives a per-window attribution
            plus the within-window sector split.

  Note: the Moebius transform on the Boolean cube is the fast zeta inversion;
  the Walsh (Fourier on {+-1}^n) spectrum of the same v differs only by a
  triangular basis change, so the two are informationally equivalent but their
  order-k decompositions are NOT the same partition.  We work in the Moebius
  basis; phi/psi are linear functionals of it -- this is what makes "sample
  the spectrum" and "estimate the attribution" the same problem with a known
  objective.

CHECKS (all must hold, else the oracle is wrong)
  * the two Shapley implementations agree to 1e-9
  * efficiency:      sum_i phi_i = v(N) - v({})
  * Owen efficiency: sum_i psi_i = v(N) - v({})
  * group consistency: sum_{i in B_w} psi_i == Shapley value of group w of the
    quotient game  v_B(T) = v(union of groups in T)
"""
import csv
import itertools
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8
WINDOW = 4


def load_v(path):
    """read the sweep csv: header name,v_ppb,nfail,rc,nrow (+ SWEEP_DONE)."""
    v = {}
    with open(path) as f:
        for row in csv.DictReader(f):
            nm = str(row.get("name", ""))
            if not nm.startswith("mob_"):
                continue
            assert int(row["nfail"]) == 0 and int(row["rc"]) == 0, row
            v[int(nm.split("_")[1])] = float(row["v_ppb"])
    assert len(v) == 1 << NPLAY, len(v)
    return v


def mobius(vv):
    out = [0.0] * (1 << NPLAY)
    for mask in range(1 << NPLAY):
        s, sub = 0.0, mask
        while True:
            if bin(mask ^ sub).count("1") % 2 == 0:
                s += vv[sub]
            else:
                s -= vv[sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out[mask] = s
    return out


def pc(m):
    return bin(m).count("1")


def name_of(mask):
    return "+".join("%s@W%d" % ("NOx" if i < 4 else "VOC", (i % 4) + 1)
                    for i in range(NPLAY) if mask & (1 << i)) or "base"


def shapley_from_mobius(m):
    return [sum(m[A] / pc(A) for A in range(1 << NPLAY) if A & (1 << i))
            for i in range(NPLAY)]


def shapley_direct(vv):
    """textbook: weighted average of marginal contributions over all orders."""
    phi = [0.0] * NPLAY
    rest = list(range(NPLAY))
    for i in range(NPLAY):
        others = [j for j in rest if j != i]
        for bits in range(1 << (NPLAY - 1)):
            S = 0
            for k, j in enumerate(others):
                if bits & (1 << k):
                    S |= 1 << j
            s = pc(S)
            w = math.factorial(s) * math.factorial(NPLAY - s - 1) / math.factorial(NPLAY)
            phi[i] += w * (vv[S | (1 << i)] - vv[S])
    return phi


def owen(vv, m):
    """Owen value, coalition structure = the WINDOW (each group = one window)."""
    groups = [[w, 4 + w] for w in range(WINDOW)]
    psi = [0.0] * NPLAY
    for g, B in enumerate(groups):
        others = [h for h in range(WINDOW) if h != g]
        for Tbits in range(1 << len(others)):
            ST = 0
            for k, h in enumerate(others):
                if Tbits & (1 << k):
                    for i in groups[h]:
                        ST |= 1 << i
            t = pc(Tbits)
            wT = (math.factorial(t) * math.factorial(WINDOW - t - 1)
                  / math.factorial(WINDOW))
            for i in B:
                Q = [j for j in B if j != i]
                for qbits in range(1 << len(Q)):
                    SQ = 0
                    for k, j in enumerate(Q):
                        if qbits & (1 << k):
                            SQ |= 1 << j
                    q = pc(SQ)
                    # random-order model: P(T before k) x P(Q before i in B_k)
                    wQ = (math.factorial(q) * math.factorial(len(B) - q - 1)
                          / math.factorial(len(B)))
                    psi[i] += wT * wQ * (vv[ST | SQ | (1 << i)] - vv[ST | SQ])
    return psi


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        HERE, "mobius_v_dummy.csv")
    v_abs = load_v(path)
    base = v_abs[0]
    vv = [v_abs[i] - base for i in range(1 << NPLAY)]
    m = mobius(vv)

    phi_mob = shapley_from_mobius(m)
    phi_dir = shapley_direct(vv)
    d1 = max(abs(a - b) for a, b in zip(phi_mob, phi_dir))
    tot = vv[(1 << NPLAY) - 1]
    d2 = abs(sum(phi_mob) - tot)

    psi = owen(vv, m)
    d3 = abs(sum(psi) - tot)

    # group consistency: sum of psi over a window == Shapley of that window in
    # the quotient game v_B(T) = v(union of groups in T)
    groups = [[w, 4 + w] for w in range(WINDOW)]
    qv = {}
    for T in range(1 << WINDOW):
        mask = 0
        for h in range(WINDOW):
            if T & (1 << h):
                for i in groups[h]:
                    mask |= 1 << i
        qv[T] = vv[mask]
    qphi = []
    for h in range(WINDOW):
        acc = 0.0
        for bits in range(1 << (WINDOW - 1)):
            others = [x for x in range(WINDOW) if x != h]
            T = 0
            for k, x in enumerate(others):
                if bits & (1 << k):
                    T |= 1 << x
            t = pc(T)
            w = (math.factorial(t) * math.factorial(WINDOW - t - 1)
                 / math.factorial(WINDOW))
            acc += w * (qv[T | (1 << h)] - qv[T])
        qphi.append(acc)
    d4 = max(abs(sum(psi[i] for i in g) - qphi[h])
             for h, g in enumerate(groups))

    L = []
    L.append("Exact attribution oracle from the full 2^8 enumeration")
    L.append("  source            : %s" % os.path.basename(path))
    L.append("  baseline metric   : %.4f ppb (mean O3, last 24 h)" % base)
    L.append("  v(N) - v({})      : %.4f ppb" % tot)
    L.append("")
    L.append("Shapley values [ppb]          (time-resolved: one number per sector x window)")
    for i in range(NPLAY):
        L.append("  %-9s %+9.4f   (%.1f%% of total)"
                 % (name_of(1 << i), phi_mob[i], 100 * phi_mob[i] / tot))
    L.append("  sum = %+.4f" % sum(phi_mob))
    L.append("")
    L.append("Owen values  [ppb]            (coalition structure = the 4 windows)")
    for w in range(WINDOW):
        lo = psi[w]
        hi = psi[w + 4]
        L.append("  W%d : NOx %+9.4f   VOC %+9.4f   window total %+9.4f"
                 % (w + 1, lo, hi, lo + hi))
    L.append("  sum = %+.4f" % sum(psi))
    L.append("")
    L.append("largest |Shapley| players:")
    for i in sorted(range(NPLAY), key=lambda i: -abs(phi_mob[i]))[:8]:
        L.append("  %-9s %+9.4f" % (name_of(1 << i), phi_mob[i]))

    # order decomposition of the Shapley value: how much of phi comes from the
    # |A|=k part of the spectrum.  phi_i = sum_{A ni i} m(A)/|A|, so this gives
    # the *truncation floor* of any low-order model: the irreducible bias of a
    # sampler that only resolves orders <= q.
    L.append("")
    L.append("order decomposition of Shapley  (contribution of |A|=k to phi):")
    L.append("     k   max_i |contrib|   as %% of max|phi|")
    maxphi = max(abs(x) for x in phi_mob)
    for k in range(1, NPLAY + 1):
        contrib = [sum(m[A] / k for A in range(1 << NPLAY)
                       if pc(A) == k and A & (1 << i)) for i in range(NPLAY)]
        cmax = max(abs(x) for x in contrib)
        L.append("  %4d   %12.4f   %14.1f%%" % (k, cmax, 100 * cmax / maxphi))
    for q in (1, 2, 3, 4):
        tr = []
        for i in range(NPLAY):
            low = sum(m[A] / pc(A) for A in range(1 << NPLAY)
                      if A & (1 << i) and pc(A) <= q)
            tr.append(abs(low - phi_mob[i]))
        L.append("  truncation floor at order <= %d : max_i |phi_i - phi_i^(<=%d)| "
                 "= %.4f ppb (%.2f%% of max|phi|)"
                 % (q, q, max(tr), 100 * max(tr) / maxphi))
    L.append("")
    L.append("CHECKS")
    L.append("  Moebius vs marginal-form Shapley : max diff %.2e" % d1)
    L.append("  efficiency  sum(phi) == v(N)     : %.2e" % d2)
    L.append("  efficiency  sum(psi) == v(N)     : %.2e" % d3)
    L.append("  window totals == quotient Shapley: %.2e" % d4)
    ok = max(d1, d2, d3, d4) < 1e-9
    L.append("  --> %s" % ("PASS" if ok else "FAIL"))
    txt = "\n".join(L)
    print(txt)
    if not ok:
        return 1
    with open(os.path.join(HERE, "attrib_oracle.txt"), "w") as f:
        f.write(txt + "\n")
    with open(os.path.join(HERE, "attrib_oracle.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["player", "kind", "set", "value_ppb", "shapley_ppb", "owen_ppb"])
        for i in range(NPLAY):
            w.writerow([i, "NOx" if i < 4 else "VOC", name_of(1 << i),
                        "", "%.9f" % phi_mob[i], "%.9f" % psi[i]])
    return 0


if __name__ == "__main__":
    sys.exit(main())
