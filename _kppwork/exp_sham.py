"""SHAM-SWITCH CONTROL, readout.

EMIS_CO = 0 makes the four CO windows NULL PLAYERS (E_CO is identically 0.0,
not merely small).  Every check below is therefore an EXACT equality, and a
failure cannot be dismissed as noise.

Note on the reference: setting EMIS_CO = 0 also removes the CO emission from
the baseline, so the frozen 8-player table (which has CO on) is NOT a valid
reference.  The checks are made self-contained instead -- the 8-player
sub-game is extracted from the same sham data -- plus one deliberately
approximate cross-check against the frozen table.

  A  nullity / bit mapping  v(mask) == v(mask & 0xFF) for all 4096 masks
  B  dummy axiom            phi_CO == 0
  C  null players are inert in the Shapley value
                            phi_live(12-player) == phi(8-player sub-game)
  D  efficiency             sum(phi) == v(full) - v(none)
  E  spectrum               m(A) == 0 for every A touching a CO bit
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
N12, M12 = 12, 1 << 12
N8, M8 = 8, 1 << 8
COMASK = 0xF00
PAIRS = [("sham_sat", "r_sat"), ("sham_lim1", "r_lim1")]


def pc(m):
    return bin(m).count("1")


def role(i):
    return "NOx" if i < 4 else ("VOC" if i < 8 else "CO")


def pname(i):
    return "%s@W%d" % (role(i), i % 4 + 1)


def mobius(v, n):
    out = [0.0] * (1 << n)
    for mask in range(1 << n):
        s, sub = 0.0, mask
        while True:
            s += (1.0 if pc(mask ^ sub) % 2 == 0 else -1.0) * v[sub]
            if sub == 0:
                break
            sub = (sub - 1) & mask
        out[mask] = s
    return out


def shapley(m, n):
    return [sum(m[A] / pc(A) for A in range(1, 1 << n) if (A >> i) & 1)
            for i in range(n)]


def load(fn, keep=None):
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
            assert int(row["nfail"]) == 0 and int(row["rc"]) == 0, row
            tab.setdefault(s, {})[int(ms)] = float(row["v_ppb"])
    return tab


def main():
    sham = load("sham_v.csv", keep={p[0] for p in PAIRS})
    frozen = load("scen_v.csv", keep={p[1] for p in PAIRS})
    L = []
    add = L.append
    add("=" * 78)
    add("SHAM-SWITCH CONTROL   EMIS_CO = 0  ->  4 CO windows are NULL PLAYERS")
    add("all checks are exact equalities; a failure is a wiring bug, not noise")
    add("=" * 78)

    allpass = True
    for s_sham, s_ref in PAIRS:
        if s_sham not in sham or len(sham[s_sham]) != M12:
            add("")
            add("scenario %s: only %d/%d coalitions -- sweep incomplete"
                % (s_sham, len(sham.get(s_sham, {})), M12))
            allpass = False
            continue
        v = [sham[s_sham][m] for m in range(M12)]
        v0 = v[0]
        vr = [x - v0 for x in v]
        mo = mobius(vr, N12)
        phi = shapley(mo, N12)

        # --- is the panel live?  (guards against a vacuous test) ----------
        nlive = len({v[a] for a in range(M8)})

        # A: the value must not depend on bits 8..11 at all
        chkA = max(abs(v[m] - v[m & 0xFF]) for m in range(M12))
        # E: no Moebius coefficient may touch a CO bit
        chkE = max(abs(mo[a]) for a in range(1, M12) if a & COMASK)
        # B: the dummy axiom
        phiCO = sum(phi[i] for i in range(8, 12))
        # C: extract the 8-player sub-game from the SAME data and compare
        v8 = [v[a] for a in range(M8)]
        mo8 = mobius([x - v8[0] for x in v8], N8)
        phi8 = shapley(mo8, N8)
        chkC = max(abs(phi[i] - phi8[i]) for i in range(8))
        # D: efficiency
        chkD = sum(phi) - (v[M12 - 1] - v0)

        # approximate cross-check vs the frozen same-regime table (CO ON)
        vref = [frozen[s_ref][a] for a in range(M8)]
        ref8 = shapley(mobius([x - vref[0] for x in vref], N8), N8)
        dcross = max(abs(phi8[i] - ref8[i]) for i in range(8))

        marks = []
        for name, val, tol in (("A", chkA, 1e-9), ("B", abs(phiCO), 1e-9),
                               ("C", chkC, 1e-9), ("D", abs(chkD), 1e-9),
                               ("E", chkE, 1e-9)):
            ok = val <= tol
            allpass &= ok
            marks.append("%s=%s" % (name, "PASS" if ok else "FAIL"))

        add("")
        add("-" * 78)
        add("scenario %s   (EMIS_CO=0, baseline %.3f ppb)   "
            "distinct values over the 256 live combos: %d"
            % (s_sham, v0, nlive))
        add("-" * 78)
        add("  A  max |v(mask) - v(mask & 0xFF)| over 4096 masks  = %.3e"
            % chkA)
        add("  E  max |m(A)| over A touching a CO bit             = %.3e"
            % chkE)
        add("  B  phi_CO (sum over the 4 CO windows)              = %.3e"
            % phiCO)
        add("  C  max |phi_live(12-ply) - phi_live(8-ply subgame)| = %.3e"
            % chkC)
        add("  D  sum(phi) - [v(full) - v(none)]                   = %.3e"
            % chkD)
        add("     -> %s" % "   ".join(marks))
        add("")
        add("  live players, 8-player sub-game extracted from the sham data:")
        for i in range(8):
            add("    %-8s phi = %+8.4f" % (pname(i), phi8[i]))
        add("  cross-check vs frozen %s (CO ON; expected to differ only by"
            % s_ref)
        add("  the CO sector's ~1%% share): max |dphi| = %.4f ppb" % dcross)

    add("")
    add("=" * 78)
    add("VERDICT: %s" % ("ALL CHECKS PASS -- the attribution pipeline is "
                         "numerically exact on a known zero"
                         if allpass else "FAILURE -- see above"))
    add("=" * 78)
    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "sham_result.txt"), "w") as f:
        f.write(txt + "\n")
    return 0 if allpass else 1


if __name__ == "__main__":
    sys.exit(main())
