#!/usr/bin/env python3
"""Regime positioning + failure-map test at the DALLAS 2016 winter anchored slice.

Second leg of the two-city bracket (Table 9).  Deer Park 2017 summer covers the
NOx-limited side (results_deerpark_lvg.txt); Dallas winter is expected to sit
on the NOx-saturated side, i.e. where the failure map says local derivatives
should break down.  Same protocol as the Deer Park script: at the anchored
slice, run local-vs-global (3 sectors x 8 three-hour windows, off/+delta/-delta,
73 runs) and read the regime position off the NOx response signs.

Map predictions, stated BEFORE the run (results_scan_flips.txt):
  * saturated side -> NOx finite effects should flip sign vs D_local somewhere;
  * limited side   -> NOx monotone, failure disappears (what Deer Park found);
  * near ridge     -> zero-effect / flip windows in the afternoon (12-15h).

The anchor knobs are read from results_fit_dallas_2016w.txt (8-knob protocol);
require them all, so a stale 6-knob file fails loudly instead of silently
running the wrong slice.

Writes results_dallas_lvg.txt
"""
import sys, os

sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box  # noqa: E402

D = "/mnt/d/atmos"
LABEL = "dallas_2016w"
OUT = os.path.join(D, "results_dallas_lvg.txt")
DELTA = 0.05
SECTORS = ["NOXW", "VOCW", "COW"]
NW = 8
ALL1 = [1.0] * NW
KNOBS = ("EMIS_NOX", "EMIS_VOC", "EMIS_ISO", "EMIS_CO",
         "KDIL_DAY", "KDIL_NIGHT", "BG_O3", "BG_CO")


def read_anchor(label):
    path = os.path.join(D, "results_fit_%s.txt" % label)
    knobs = {}
    with open(path) as fh:
        for line in fh:
            if "=" not in line:
                continue
            k, _, v = line.partition("=")
            k = k.strip()
            if k in KNOBS:
                knobs[k] = float(v)
    missing = [k for k in KNOBS if k not in knobs]
    if missing:
        raise SystemExit("%s is not an 8-knob anchor file; missing %s" % (path, missing))
    return knobs


def read_met_tminmax(label):
    path = os.path.join(D, "met_%s.csv" % label)
    T = []
    with open(path) as fh:
        for r in fh:
            p = r.strip().split(",")
            if p[0] == "time":
                continue
            try:
                T.append(float(p[1]))
            except (ValueError, IndexError):
                pass
    return min(T) + 273.15, max(T) + 273.15


def run(cfg, tag):
    m = kpp_box.run_scenario(cfg, "dal_" + tag, "scale")
    return _metrics(m["traj"])


def _metrics(traj):
    tt = [r[0] for r in traj]
    o3 = [r[1] for r in traj]
    sel = [i for i in range(len(tt)) if tt[i] >= 48.0]
    v48 = sum(o3[i] for i in sel) / len(sel)
    dt = tt[1] - tt[0]
    w = max(1, int(round(8.0 / dt)))
    best = float("-inf")
    for i in range(w - 1, len(o3)):
        if tt[i] < 48.0:
            continue
        a = sum(o3[i - w + 1:i + 1]) / w
        if a > best:
            best = a
    return v48, best


def main():
    anchor = read_anchor(LABEL)
    tmin, tmax = read_met_tminmax(LABEL)
    base_cfg = dict(anchor)
    # backgrounds the fit's cost function sets (make_cost in _fit_anchor.py);
    # omitting BG_NO2/BG_NO inflates the steady-state O3.
    base_cfg.update(BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
                    TMIN_K=tmin, TMAX_K=tmax)
    L = []
    add = L.append
    add("Dallas 2016 winter: regime positioning at the anchored slice")
    add("anchor: %s  (TMIN %.1f K, TMAX %.1f K)" % (anchor, tmin, tmax))
    add("")
    base = run(base_cfg, "base")
    add("base v48 = %.4f  mda8 = %.4f" % base)
    add("")
    add("%-6s %-8s %10s %10s %10s %9s %7s" %
        ("sector", "window", "D_local", "D_finite", "-D_local", "ratio", "sign"))
    flips8 = []
    flips48 = []
    zero8 = []
    for sec in SECTORS:
        for i in range(NW):
            off = list(ALL1); off[i] = 0.0
            plus = list(ALL1); plus[i] = 1.0 + DELTA
            minus = list(ALL1); minus[i] = 1.0 - DELTA
            vo = run({**base_cfg, sec: off}, "%s_off%d" % (sec, i))
            vp = run({**base_cfg, sec: plus}, "%s_p%d" % (sec, i))
            vm = run({**base_cfg, sec: minus}, "%s_m%d" % (sec, i))
            for tag, j in (("8", 1), ("48", 0)):
                Dloc = (vp[j] - vm[j]) / (2 * DELTA)
                Dfin = vo[j] - base[j]
                flip = (Dloc * (-Dfin)) <= 0 and abs(Dloc) > 1e-9 and abs(Dfin) > 1e-9
                zero = abs(Dfin) > 0.5 and abs(Dloc) < 0.1 * abs(Dfin)
                if tag == "8":
                    flips8.append((sec, i, Dloc, Dfin, flip, zero))
                else:
                    flips48.append((sec, i, Dloc, Dfin, flip, zero))
            Dloc = (vp[1] - vm[1]) / (2 * DELTA)
            Dfin = vo[1] - base[1]
            ratio = Dfin / (-Dloc) if abs(Dloc) > 1e-9 else float("nan")
            s = "FLIP" if (Dloc * (-Dfin)) <= 0 and abs(Dfin) > 1e-9 else "ok"
            add("%-6s %-8s %10.3f %10.3f %10.3f %9.2f %7s" %
                (sec, "%d-%d" % (3 * i, 3 * i + 3), Dloc, Dfin, -Dloc, ratio, s))
    add("")
    nf8 = sum(1 for f in flips8 if f[4])
    nf48 = sum(1 for f in flips48 if f[4])
    nz8 = sum(1 for f in flips8 if f[5])
    add("MDA8  sign flips: %d/24   zero-effect windows: %d" % (nf8, nz8))
    add("v48   sign flips: %d/24" % nf48)
    for sec, i, d1, df, fl, ze in flips8:
        if fl or ze:
            add("  %s@%d-%d  D1=%+.3f  Dfin=%+.3f  flip=%s zero=%s" %
                (sec, 3 * i, 3 * i + 3, d1, df, fl, ze))
    add("")
    # regime verdict from the NOx response signs (saturation test):
    # switching NOx OFF raises O3  <=> finite effect positive <=> saturated
    sat = sum(1 for sec, i, d1, df, fl, ze in flips8 if sec == "NOXW" and df > 0)
    add("NOx windows where switch-off RAISES MDA8 (saturated behaviour): %d of 8" % sat)
    add("=> regime read-off: %s" %
        ("SATURATED side" if sat >= 5 else
         "LIMITED side" if sat <= 2 else
         "NEAR RIDGE (mixed)"))
    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()