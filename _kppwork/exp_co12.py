"""12-player Möbius/Shapley analysis: does the third sector (CO) change the
regime story?

Reads co_v.csv (2 scenarios x 4096 coalitions, /root/kpp6).

Player layout (bit -> switch):
    bit 0..3   NOx window W1..W4 OFF
    bit 4..7   VOC window W1..W4 OFF
    bit 8..11  CO  window W1..W4 OFF

Three questions:
  Q1  How much of the total effect does the CO sector carry at all?
  Q2  Does adding CO as a player move the NOx/VOC Shapley values?  The 8-player
      game is embedded exactly (CO bits held at 0), and a subset Möbius
      coefficient m(A) with A inside the 8-player set is IDENTICAL in both
      games -- but phi_i is not, because the averages run over coalitions that
      may now contain CO players with weights 1/|A|.  So the comparison is
      exactly "does the third sector re-weight the other two's attribution?".
  Q3  Does the regime flip (saturated vs limited) survive?

Self-checks (must hold to machine precision):
    sum_A m(A) == v(all off) - v(none)          (Möbius telescoping)
    sum_i phi_i == v(all off) - v(none)         (Shapley efficiency)
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
N12, M12 = 12, 1 << 12
N8, M8 = 8, 1 << 8
SCENS = ["r_sat", "r_lim1"]


def pc(m):
    return bin(m).count("1")


def role(i):
    return "NOx" if i < 4 else ("VOC" if i < 8 else "CO")


def pname(i):
    return "%s@W%d" % (role(i), i % 4 + 1)


def mask_name(m, n=N12):
    return "+".join(pname(i) for i in range(n) if (m >> i) & 1) or "none"


def mobius(v, n):
    """Exact subset Möbius transform (Harsanyi dividends), O(3**n).

    Not the Walsh spectrum: same information, but related to it by a
    triangular basis change, so the order-k splits differ.
    """
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
    """phi_i = sum_{A contains i} m(A)/|A|  -- the exact linear functional."""
    return [sum(m[A] / pc(A) for A in range(1, 1 << n) if (A >> i) & 1)
            for i in range(n)]


def levels(m, n):
    """R_k = max_{|A|=k} |m(A)| / max_i |m({i})|  (R_1 == 1 by construction)."""
    byk = {}
    for mask in range(1, 1 << n):
        byk.setdefault(pc(mask), []).append(abs(m[mask]))
    mx1 = max(abs(m[1 << i]) for i in range(n))
    return {k: max(v) / mx1 for k, v in byk.items()}, mx1


def load():
    tab = {}
    with open(os.path.join(HERE, "co_v.csv")) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            scen, ms = nm.rsplit("_m", 1)
            if row["v_ppb"] in ("", "NaN"):
                continue
            assert int(row["nfail"]) == 0 and int(row["rc"]) == 0, row
            tab.setdefault(scen, {})[int(ms)] = float(row["v_ppb"])
    return tab


def load8():
    """The frozen 8-player table, for the sub-game comparison."""
    tab = {}
    with open(os.path.join(HERE, "scen_v.csv")) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            scen, ms = nm.rsplit("_m", 1)
            if scen in SCENS:
                tab.setdefault(scen, {})[int(ms)] = float(row["v_ppb"])
    return tab


def main():
    tab = load()
    t8 = load8()
    L = []
    add = L.append

    for s in SCENS:
        if s not in tab or len(tab[s]) != M12:
            add("scenario %s: only %d/%d coalitions -- sweep incomplete"
                % (s, len(tab.get(s, {})), M12))
            continue
        v = [tab[s][m] for m in range(M12)]
        v0 = v[0]
        vr = [x - v0 for x in v]
        m = mobius(vr, N12)
        phi = shapley(m, N12)
        lv, mx1 = levels(m, N12)

        add("")
        add("=" * 78)
        add("scenario %s   12 players (3 sectors x 4 windows)   baseline %.2f ppb"
            % (s, v0))
        add("=" * 78)

        # ---- self-checks ------------------------------------------------
        chk_m = sum(m) - (v[M12 - 1] - v0)
        chk_p = sum(phi) - (v[M12 - 1] - v0)
        add("  SELF-CHECK  sum m(A) - [v(full)-v(none)] = %.3e" % chk_m)
        add("  SELF-CHECK  sum phi   - [v(full)-v(none)] = %.3e" % chk_p)
        # V-gate: with CO held ON the 12-player table IS the frozen 8-player
        # table, so the two must agree bit for bit on all 256 coalitions.
        vg = max(abs(tab[s][mm] - t8[s][mm]) for mm in range(M8))
        add("  V-GATE      max |v12(m) - v8(m)| over 256 CO-on coalitions"
            " = %.2e ppb" % vg)

        # ---- Q1: sector totals -----------------------------------------
        add("")
        add("  Q1  sector Shapley totals:      (share of sum|phi|)")
        tot = sum(abs(x) for x in phi)
        for r in ("NOx", "VOC", "CO"):
            idx = [i for i in range(N12) if role(i) == r]
            sr = sum(phi[i] for i in idx)
            add("      %-4s sum phi = %+8.2f ppb   (%5.1f%%)   windows: %s"
                % (r, sr, 100 * abs(sr) / tot,
                   "  ".join("%s %+0.2f" % (pname(i), phi[i]) for i in idx)))

        # ---- Q2: 8-player sub-game, same scenario ----------------------
        v8 = [t8[s][mm] for mm in range(M8)]
        r8 = [x - v8[0] for x in v8]
        m8 = mobius(r8, N8)
        phi8 = shapley(m8, N8)
        add("")
        add("  Q2  NOx/VOC attribution with and without the CO player:")
        add("      %-8s %12s %12s %8s" % ("player", "phi (8-ply)", "phi (12-ply)", "delta"))
        for i in range(8):
            add("      %-8s %12.2f %12.2f %+8.2f"
                % (pname(i), phi8[i], phi[i], phi[i] - phi8[i]))
        add("      %-8s %12.2f %12.2f %+8.2f"
            % ("NOx tot", sum(phi8[:4]), sum(phi[:4]),
               sum(phi[:4]) - sum(phi8[:4])))
        add("      %-8s %12.2f %12.2f %+8.2f"
            % ("VOC tot", sum(phi8[4:8]), sum(phi[4:8]),
               sum(phi[4:8]) - sum(phi8[4:])))

        # ---- Q3: regime flip -------------------------------------------
        add("")
        add("  Q3  regime signature")
        best = min(range(M12), key=lambda mm: v[mm])
        nox_all_off = v[0x00F] - v0
        add("      best recipe   = %s   v = %.2f  (%+.2f ppb vs none)"
            % (mask_name(best), v[best], v[best] - v0))
        add("      all-NOx off   = %+.2f ppb   -> NOx control %s"
            % (nox_all_off, "worsens O3" if nox_all_off > 0 else "improves O3"))
        add("      all-CO off    = %+.2f ppb" % (v[0xF00] - v0))
        add("      spectral levels R_k (max|m_k| / max|m_1|):")
        add("        " + "  ".join("R%d=%.3f" % (k, lv[k])
                                  for k in sorted(lv) if k <= 4))

        # ---- top interactions ------------------------------------------
        inter = sorted((abs(m[a]), a) for a in range(1, M12) if pc(a) >= 2)
        add("      top-5 interaction coefficients |m(A)|, |A|>=2:")
        for mag, a in inter[-5:][::-1]:
            add("        %-44s |m|=%7.3f  m=%+8.3f  |A|=%d"
                % (mask_name(a), mag, m[a], pc(a)))

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "co12_analysis.txt"), "w") as f:
        f.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
