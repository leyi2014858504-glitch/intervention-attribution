"""Why does adding the CO player barely move the NOx/VOC Shapley values?

Exact reason: in the game with players P u {i}, the Shapley value of j != i
changes only through Moebius coefficients m(A) with i in A.  So the size of
the perturbation is bounded by the size of the CO-touching coefficients.
This script measures exactly that, instead of arguing informally.

Also verifies the embedding identity: m(A) for A inside the 8-player set must
be IDENTICAL in the 12- and 8-player games (because v agrees on that subcube).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
N12, M12 = 12, 1 << 12
N8, M8 = 8, 1 << 8
COMASK = 0xF00
SCENS = ["r_sat", "r_lim1"]


def pc(m):
    return bin(m).count("1")


def role(i):
    return "NOx" if i < 4 else ("VOC" if i < 8 else "CO")


def pname(i):
    return "%s@W%d" % (role(i), i % 4 + 1)


def mask_name(m):
    return "+".join(pname(i) for i in range(N12) if (m >> i) & 1) or "none"


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


def load(fn, scens, mdim):
    tab = {}
    with open(os.path.join(HERE, fn)) as f:
        for row in csv.DictReader(f):
            nm = row["name"]
            if "_m" not in nm:
                continue
            s, ms = nm.rsplit("_m", 1)
            if s not in scens:
                continue
            if row["v_ppb"] in ("", "NaN"):
                continue
            tab.setdefault(s, {})[int(ms)] = float(row["v_ppb"])
    return tab


def main():
    t12 = load("co_v.csv", SCENS, M12)
    t8 = load("scen_v.csv", SCENS, M8)
    L = []
    for s in SCENS:
        if len(t12.get(s, {})) != M12 or len(t8.get(s, {})) != M8:
            L.append("%s incomplete" % s)
            continue
        v12 = [t12[s][m] for m in range(M12)]
        v8 = [t8[s][m] for m in range(M8)]
        m12 = mobius([x - v12[0] for x in v12], N12)
        m8 = mobius([x - v8[0] for x in v8], N8)

        # ---- embedding identity: coefficients inside the 8-player set ----
        emb = max(abs(m12[a] - m8[a]) for a in range(M8))

        # ---- CO-touching vs CO-free coefficient mass ---------------------
        free_mag = sum(abs(m12[a]) for a in range(1, M12) if not (a & COMASK))
        touch = [a for a in range(1, M12) if a & COMASK]
        touch_mag = sum(abs(m12[a]) for a in touch)
        free_max = max(abs(m12[a]) for a in range(1, M12) if not (a & COMASK))
        touch_max = max(abs(m12[a]) for a in touch)

        L.append("")
        L.append("=" * 78)
        L.append("scenario %s   (8-player m dropped into 12-player table)"
                 % s)
        L.append("=" * 78)
        L.append("  EMBEDDING CHECK  max |m12(A) - m8(A)|, A inside 8-player set"
                 " = %.3e" % emb)
        L.append("")
        L.append("  coefficient mass, split by whether the coalition touches CO:")
        L.append("    A does NOT touch CO :  max|m| = %7.3f   sum|m| = %8.2f"
                 % (free_max, free_mag))
        L.append("    A DOES touch CO     :  max|m| = %7.3f   sum|m| = %8.2f"
                 % (touch_max, touch_mag))
        L.append("    ratio sum|m| (touch/free) = %.3f" % (touch_mag / free_mag))

        # per-order breakdown
        L.append("")
        L.append("  by order k (sum |m_k|, split):")
        L.append("    %-4s %12s %12s" % ("k", "won't-touch", "touches-CO"))
        for k in range(1, N12 + 1):
            af = sum(abs(m12[a]) for a in range(1, M12)
                     if pc(a) == k and not (a & COMASK))
            at = sum(abs(m12[a]) for a in range(1, M12)
                     if pc(a) == k and (a & COMASK))
            if af or at:
                L.append("    %-4d %12.3f %12.3f" % (k, af, at))

        L.append("")
        L.append("  top-5 coefficients that TOUCH CO:")
        top = sorted(touch, key=lambda a: -abs(m12[a]))[:5]
        for a in top:
            L.append("    %-58s m=%+8.4f  |A|=%d" % (mask_name(a), m12[a], pc(a)))

        # ---- direct sensitivity: phi_j with and without CO --------------
        def shap(m, n, i):
            return sum(m[A] / pc(A) for A in range(1, 1 << n) if (A >> i) & 1)
        ph12 = [shap(m12, N12, i) for i in range(N12)]
        ph8 = [shap(m8, N8, i) for i in range(N8)]
        L.append("")
        L.append("  perturbation of NOx/VOC Shapley values: max |dphi| = %.4f ppb"
                 % max(abs(ph12[i] - ph8[i]) for i in range(8)))

    txt = "\n".join(L)
    print(txt)
    with open(os.path.join(HERE, "co12_coeff.txt"), "w") as f:
        f.write(txt + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
