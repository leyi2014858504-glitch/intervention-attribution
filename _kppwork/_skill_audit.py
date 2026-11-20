#!/usr/bin/env python3
"""Audit the box-model setup against the transferable parts of the Cantera skill:
provenance hashes, numerical refinement (output spacing), and a conservation
budget on a nearly conserved species (CO).

Writes results_skill_audit.txt
"""
import hashlib
import math
import os
import statistics
import sys

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402
from _gen16_city import build_base  # noqa: E402

D = "/mnt/d/atmos"
GAMES = ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s")
COLS = {"o3": 1, "nox": 2, "co": 7, "sday": 4, "kdil": 5}


def base_of(game):
    if game in ("r_sat", "r_lim1"):
        cfg = {"TEND_S": 259200.0}
        if game == "r_lim1":
            cfg.update(NOX_FACTOR=0.5, VOC_FACTOR=1.5)
        return cfg
    cfg = dict(build_base(game)[0])
    cfg["TEND_S"] = 259200.0
    return cfg


def gp(mu, sd, h):
    return math.exp(-(((h - mu) / sd) ** 2))


def profile(hh):
    return 0.45 + 1.30 * (gp(8.0, 1.6, hh) + gp(18.5, 2.0, hh))


def mda8(traj):
    o3 = [r[1] for r in traj]
    t = [r[0] for r in traj]
    dt = t[1] - t[0]
    w = max(1, int(round(8.0 / dt)))
    best = -1e9
    for i in range(w - 1, len(o3)):
        if t[i] < 48.0:
            continue
        best = max(best, sum(o3[i - w + 1:i + 1]) / w)
    return best


def d3(traj):
    v = [r[1] for r in traj if r[0] >= 48.0]
    return statistics.mean(v)


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()[:16]


L = []
add = L.append
add("Audit against the Cantera skill's verification contract")
add("(provenance, numerical refinement, conservation, scope)")
add("")

# ---------- 1. provenance ----------
add("1. Provenance")
for p in ("/root/kpp7/cb05_Main.f90", "/root/kppa/cb05_Main.f90"):
    if os.path.exists(p):
        add("  %-34s sha256:%s  %d bytes" % (p, sha(p), os.path.getsize(p)))
mech = [os.path.join("/root/kpp7", f) for f in sorted(os.listdir("/root/kpp7"))
        if f.endswith((".eqn", ".def", ".spc", ".kpp"))]
for p in mech:
    add("  %-34s sha256:%s" % (os.path.basename(p), sha(p)))
add("")

# ---------- 2. output-spacing refinement ----------
add("2. Output-spacing refinement (DT_S; 3-day runs, all players ON)")
add("%-16s %6s %10s %10s %10s %10s" % ("game", "DT", "MDA8", "d3", "dMDA8%", "dd3%"))
for g in GAMES:
    got = {}
    for dt in (150.0, 300.0, 600.0):
        cfg = base_of(g)
        cfg["DT_S"] = dt
        m = kpp_box.run_scenario(cfg, "audit_dt%d_%s" % (int(dt), g), "anchor")
        got[dt] = (mda8(m["traj"]), d3(m["traj"]))
        if dt == 300.0:
            neg = min(min(r[1:]) for r in m["traj"])
    ref = got[300.0]
    for dt in (150.0, 300.0, 600.0):
        v8, v3 = got[dt]
        add("%-16s %6d %10.4f %10.4f %+10.3f %+10.3f" %
            (g, int(dt), v8, v3,
             100 * (v8 - ref[0]) / ref[0], 100 * (v3 - ref[1]) / ref[1]))
    add("  min over all columns of the DT=300 trajectory: %.3e" % neg)
    add("")

# ---------- 3. conservation budget on CO ----------
def co_budget(cfg, tag):
    m = kpp_box.run_scenario(cfg, tag, "anchor")
    tr = [r for r in m["traj"] if 48.0 <= r[0] <= 72.0]
    emis = cfg.get("EMIS_CO", kpp_box.DEFAULTS["EMIS_CO"])
    cbg = cfg.get("BG_CO", kpp_box.DEFAULTS["BG_CO"])
    e_tot = loss_tot = 0.0
    for a, b in zip(tr, tr[1:]):
        dth = b[0] - a[0]
        ha, hb = a[0] % 24.0, b[0] % 24.0
        e_tot += 0.5 * emis * (profile(ha) + profile(hb)) * dth
        fa = a[COLS["kdil"]] * (a[COLS["co"]] - cbg)
        fb = b[COLS["kdil"]] * (b[COLS["co"]] - cbg)
        loss_tot += 0.5 * (fa + fb) * dth
    dC = tr[-1][COLS["co"]] - tr[0][COLS["co"]]
    resid = dC - e_tot + loss_tot
    return dC, e_tot, loss_tot, resid


add("3. CO budget closure, day 3 (48-72 h), Dallas all-on")
dC, e_tot, loss_tot, resid = co_budget(base_of("dallas_2016w"), "audit_co_budget")
add("  full chemistry: dC %.3f; emission %.3f; dilution %.3f; residual %.3f (%.1f%% of emission)" %
    (dC, e_tot, loss_tot, resid, 100 * resid / e_tot))
add("  the residual here is net secondary CO production from VOC oxidation")
cfg = base_of("dallas_2016w")
cfg.update(EMIS_VOC=0.0, EMIS_ISO=0.0, BG_VOC=0.0, BG_ISO=0.0)
dC, e_tot, loss_tot, resid = co_budget(cfg, "audit_co_budget_novoc")
add("  no-VOC control:  dC %.3f; emission %.3f; dilution %.3f; residual %.3f (%.2f%% of emission)" %
    (dC, e_tot, loss_tot, resid, 100 * resid / e_tot))
add("  same control at three output spacings (is the residual numerical or chemical?)")
for dt in (150.0, 300.0, 600.0):
    c = dict(cfg)
    c["DT_S"] = dt
    dC, e_tot, loss_tot, resid = co_budget(c, "audit_co_novoc_dt%d" % int(dt))
    add("    DT=%3d s: residual %+7.3f ppb (%+6.2f%% of emission)" %
        (int(dt), resid, 100 * resid / e_tot))

txt = "\n".join(L)
print(txt)
with open(os.path.join(D, "results_skill_audit.txt"), "w") as fh:
    fh.write(txt + "\n")
print("wrote results_skill_audit.txt")
