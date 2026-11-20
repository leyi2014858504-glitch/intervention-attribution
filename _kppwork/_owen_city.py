#!/usr/bin/env python3
"""Owen values under both coalition structures, exact, from the spectra.

Owen is a linear functional of the dividends: psi_i = sum_{A ni i} m(A) w_i(A),
with the closed-form weight (Owen 1977; Winter 1989)

  w_i(A) = S(|A n P(i)|, |P(i)|) * G(|A/P|, q)
  S(k,p) = sum_{t=k-1}^{p-1} C(p-k, t-k+1) * t!(p-t-1)!/p!
  G(g,q) = sum_{r=g-1}^{q-1} C(q-g, r-g+1) * r!(q-r-1)!/q!

where q = number of groups, |A/P| = groups hit by A.  Two structures on the
16-player game (bits 0-7 = NOx windows, 8-15 = VOC windows):
  P_window  : 8 groups of 2   {NOx@w, VOC@w}          (time-first)
  P_species : 2 groups of 8   {all NOx}, {all VOC}    (species-first)
Shapley is the singleton-partition case.  Truncated Owen at order q drops
dividends with |A| > q.  Efficiency check: sum_i psi_i = v(N) - v(empty).

Writes results_owen_city.txt
"""
import csv
import os
import platform
from math import comb, factorial

import numpy as np

if platform.system() == "Windows":
    W = r"D:\atmos\_kppwork"
    ROOT = r"D:\atmos"
else:
    W = "/mnt/d/atmos/_kppwork"
    ROOT = "/mnt/d/atmos"

N = 16
GAMES = {
    "dallas_2016w": "city",
    "deerpark_2017s": "city",
    "r_sat": "synth",
    "r_lim1": "synth",
}
STRUCTS = {
    "window": [frozenset([w, 8 + w]) for w in range(8)],      # 8 groups of 2
    "species": [frozenset(range(8)), frozenset(range(8, 16))],  # 2 groups of 8
}
WINDOW_LABELS = ("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-21", "21-24")
PLAYER_LABELS = [("NOx", w) for w in WINDOW_LABELS] + [("VOC", w) for w in WINDOW_LABELS]


def quotient_shapley(v, groups):
    """Shapley values of the quotient game whose players are the groups."""
    q = len(groups)
    masks = [sum(1 << i for i in g) for g in groups]
    out = np.zeros(q)
    for gi in range(q):
        for T in range(1 << q):
            if T & (1 << gi):
                continue
            sm = 0
            for j in range(q):
                if T & (1 << j):
                    sm |= masks[j]
            coef = factorial(bin(T).count("1")) * factorial(q - bin(T).count("1") - 1) / factorial(q)
            out[gi] += coef * (v[sm | masks[gi]] - v[sm])
    return out


def sign_mismatch(a, b, floor):
    """Number of sign disagreements above the 1%-of-max|phi| floor."""
    return int(((np.sign(a) != np.sign(b)) & (np.abs(a) > floor) & (np.abs(b) > floor)).sum())


def S_coef(k, p):
    """Owen within-group weight: sum over supersets T of A_T."""
    s = 0.0
    for t in range(k - 1, p):
        s += comb(p - k, t - k + 1) * factorial(t) * factorial(p - t - 1) / factorial(p)
    return s


def G_coef(g, q):
    """Owen between-group weight."""
    s = 0.0
    for r in range(g - 1, q):
        s += comb(q - g, r - g + 1) * factorial(r) * factorial(q - r - 1) / factorial(q)
    return s


def owen_weights(m, groups):
    """psi_i = sum_{A ni i} m(A) w_i(A).  Returns (psi, W) with W[nz, i]."""
    q = len(groups)
    grp_of = {}
    for gi, g in enumerate(groups):
        for i in g:
            grp_of[i] = gi
    idx = np.arange(1, 1 << N)
    k_all = np.array([bin(x).count("1") for x in idx])
    # cache S and G per (k, p) and (g, q)
    p_of = {i: len(groups[grp_of[i]]) for i in range(N)}
    S_cache = {(k, p): S_coef(k, p) for k in range(1, 9) for p in (2, 8)}
    G_cache = {}
    for g in range(1, q + 1):
        G_cache[g] = G_coef(g, q)
    psi = np.zeros(N)
    W = np.zeros((len(idx), N))
    for nzi, mask in enumerate(idx):
        mv = m[mask]
        if mv == 0.0:
            continue
        mem = [i for i in range(N) if (mask >> i) & 1]
        # groups hit by A
        ghit = len({grp_of[i] for i in mem})
        gc = G_cache[ghit]
        # members per group for A
        per = {}
        for i in mem:
            per[grp_of[i]] = per.get(grp_of[i], 0) + 1
        for i in mem:
            k = per[grp_of[i]]
            w = S_cache[(k, p_of[i])] * gc
            W[nzi, i] = w
            psi[i] += mv * w
    return psi, W


def mobius_np(vals):
    f = np.array(vals, dtype=float)
    idx = np.arange(1 << N)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def shapley(m):
    idx = np.arange(1, 1 << N)
    k = np.array([bin(x).count("1") for x in idx], dtype=float)
    w = m[idx] / k
    phi = np.zeros(N)
    for i in range(N):
        phi[i] = w[((idx >> i) & 1) > 0].sum()
    return phi, w


def read_vals(label, col):
    if label in ("dallas_2016w", "deerpark_2017s"):
        v = np.zeros(65536)
        for K in range(1, 13):
            p = os.path.join(W, "n16_vh_%s_m8_%d.csv" % (label, K))
            if not os.path.exists(p):
                continue
            with open(p) as fh:
                for row in csv.DictReader(fh):
                    val = row.get(col, "")
                    if val not in ("", "NaN"):
                        v[int(row["name"][1:])] = float(val)
        return v
    v = np.zeros(65536)
    with open(os.path.join(W, "n16_vh_%s_m8.csv" % label)) as fh:
        for row in csv.DictReader(fh):
            val = row.get(col, "")
            if val not in ("", "NaN"):
                v[int(row["name"][1:])] = float(val)
    return v


L = []
add = L.append
player_rows = []
group_rows = []
summary_rows = []
add("Owen values under two coalition structures, exact, from the spectra")
add("P_window  = 8 groups of 2 {NOx@w, VOC@w}       (time-first)")
add("P_species = 2 groups of 8 {all NOx}, {all VOC}  (species-first)")
add("trunc(q) = max_i |psi_i - psi_i^{<=q}| / max_i |psi_i|")
add("")

for label in GAMES:
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3_mean")):
        v = read_vals(label, col)
        m = mobius_np(v)
        phi, w_sh = shapley(m)
        scale_sh = np.abs(phi).max() or 1.0
        add("=== %s  %s ===   max|phi| = %.4f" % (label, tag, scale_sh))
        vN = v[-1] - v[0]
        rows = {}
        for sname, groups in STRUCTS.items():
            psi, Wmat = owen_weights(m, groups)
            rows[sname] = (psi, Wmat)
            eff = psi.sum() - vN
            scale = np.abs(psi).max() or 1.0
            add("  %-8s psi: " % sname +
                " ".join("NOx%d=%+.3f" % (w, psi[w]) for w in range(4)) +
                " ...  max|psi|=%.4f  eff-err=%.2e" % (scale, eff))
            add("           " +
                " ".join("VOC%d=%+.3f" % (w, psi[8 + w]) for w in range(4)))
        # full player/group export and independent quotient-game checks
        group_err = {}
        for i, (sector, win) in enumerate(PLAYER_LABELS):
            player_rows.append([label, tag, i, sector, win,
                                "%.12g" % phi[i],
                                "%.12g" % rows["window"][0][i],
                                "%.12g" % rows["species"][0][i]])
        for sname, groups in STRUCTS.items():
            qsh = quotient_shapley(v, groups)
            totals = np.array([rows[sname][0][list(g)].sum() for g in groups])
            shapley_totals = np.array([phi[list(g)].sum() for g in groups])
            group_err[sname] = float(np.abs(totals - qsh).max())
            for gi, g in enumerate(groups):
                glabel = WINDOW_LABELS[gi] if sname == "window" else ("NOx" if gi == 0 else "VOC")
                group_rows.append([label, tag, sname, gi + 1, glabel,
                                   "%.12g" % shapley_totals[gi], "%.12g" % totals[gi],
                                   "%.12g" % qsh[gi]])
        # partition sensitivity: max_i |psi_win - psi_spc| / max|phi|
        d = np.abs(rows["window"][0] - rows["species"][0]).max() / scale_sh
        dsh_w = np.abs(rows["window"][0] - phi).max() / scale_sh
        dsh_s = np.abs(rows["species"][0] - phi).max() / scale_sh
        add("  partition spread |psi_win - psi_spc|_max / max|phi| = %.3f" % d)
        add("  vs Shapley: |psi_win - phi| = %.3f   |psi_spc - phi| = %.3f" % (dsh_w, dsh_s))
        # truncation curves
        idx = np.arange(1, 1 << N)
        k_all = np.array([bin(x).count("1") for x in idx])
        peaks = {}
        curves = {}
        for sname in ("window", "species"):
            psi, Wmat = rows[sname]
            scale = np.abs(psi).max() or 1.0
            errs = []
            for qq in range(1, N + 1):
                sel = k_all <= qq
                psi_q = (m[idx[sel]][:, None] * Wmat[sel]).sum(axis=0)
                errs.append(np.abs(psi - psi_q).max() / scale)
            add("  trunc(%-7s):" % sname +
                "  ".join("q<=%-2d %6.1f%%" % (qq, 100 * e) for qq, e in enumerate(errs, 1)))
            peaks[sname] = (max(errs), int(np.argmax(errs)) + 1)
            curves[sname] = errs
        # Shapley truncation for reference
        errs = []
        for qq in range(1, N + 1):
            sel = k_all <= qq
            phi_q = np.zeros(N)
            for i in range(N):
                phi_q[i] = (w_sh[sel & (((idx >> i) & 1) > 0)]).sum()
            errs.append(np.abs(phi - phi_q).max() / scale_sh)
        add("  trunc(shapley):" +
            "  ".join("q<=%-2d %6.1f%%" % (qq, 100 * e) for qq, e in enumerate(errs, 1)))
        peaks["shapley"] = (max(errs), int(np.argmax(errs)) + 1)
        curves["shapley"] = errs
        floor = 0.01 * scale_sh
        summary_rows.append({
            "game": label, "metric": tag, "total": vN, "max_phi": scale_sh,
            "spread": d, "win_shapley": dsh_w, "species_shapley": dsh_s,
            "sign_win": sign_mismatch(rows["window"][0], phi, floor),
            "sign_species": sign_mismatch(rows["species"][0], phi, floor),
            "sign_owen": sign_mismatch(rows["window"][0], rows["species"][0], floor),
            "eff_win": abs(rows["window"][0].sum() - vN),
            "eff_species": abs(rows["species"][0].sum() - vN),
            "group_err_win": group_err["window"],
            "group_err_species": group_err["species"],
            "top_shapley": PLAYER_LABELS[int(np.argmax(np.abs(phi)))],
            "top_window": PLAYER_LABELS[int(np.argmax(np.abs(rows["window"][0])))],
            "top_species": PLAYER_LABELS[int(np.argmax(np.abs(rows["species"][0])))],
            "peaks": peaks, "curves": curves,
        })
        add("")

txt = "\n".join(L)
print(txt)
with open(os.path.join(ROOT, "results_owen_city.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_owen_city.txt")

with open(os.path.join(ROOT, "results_three_kernels_players.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["game", "metric", "player_index", "sector", "window",
                 "shapley_ppb", "owen_time_ppb", "owen_species_ppb"])
    wr.writerows(player_rows)
with open(os.path.join(ROOT, "results_three_kernels_groups.csv"), "w", newline="") as fh:
    wr = csv.writer(fh)
    wr.writerow(["game", "metric", "partition", "group_index", "group_label",
                 "shapley_sum_ppb", "owen_sum_ppb", "quotient_shapley_ppb"])
    wr.writerows(group_rows)

S = ["Three allocation kernels, exact from the n=16 response tables",
     "shapley = singleton partition; window = time-first Owen (8 x 2); species = species-first Owen (2 x 8)",
     "group_err = max |sum_i in group psi_i - quotient Shapley|; expected ~0 for Owen",
     "sign mismatches use a 1%-of-max|phi| floor",
     ""]
for r in summary_rows:
    S.append("=== %s  %s ===" % (r["game"], r["metric"]))
    S.append("  total v(N)-v(0) = %.4f ppb   max|phi| = %.4f ppb" % (r["total"], r["max_phi"]))
    S.append("  partition spread = %.3f; vs Shapley: window %.3f, species %.3f" %
             (r["spread"], r["win_shapley"], r["species_shapley"]))
    S.append("  sign mismatches: window-vs-Shapley %d, species-vs-Shapley %d, window-vs-species %d" %
             (r["sign_win"], r["sign_species"], r["sign_owen"]))
    S.append("  efficiency errors: window %.2e, species %.2e; quotient errors: window %.2e, species %.2e" %
             (r["eff_win"], r["eff_species"], r["group_err_win"], r["group_err_species"]))
    S.append("  dominant |player|: Shapley %s, time-first %s, species-first %s" %
             ("%s@%s" % r["top_shapley"], "%s@%s" % r["top_window"], "%s@%s" % r["top_species"]))
    for kn in ("shapley", "window", "species"):
        c = r["curves"][kn]
        S.append("  %-8s trunc peak %.1f%% at q=%d; q1 %.1f%%, q2 %.1f%%, q5 %.1f%%, q10 %.1f%%" %
                 (kn, 100 * r["peaks"][kn][0], r["peaks"][kn][1],
                  100 * c[0], 100 * c[1], 100 * c[4], 100 * c[9]))
    S.append("")
summary_txt = "\n".join(S)
with open(os.path.join(ROOT, "results_three_kernels_summary.txt"), "w") as fh:
    fh.write(summary_txt + "\n")
print(summary_txt)
print("wrote results_three_kernels_players.csv / _groups.csv / _summary.txt")
