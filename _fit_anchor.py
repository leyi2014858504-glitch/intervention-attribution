#!/usr/bin/env python3
"""
Anchor the box to AQS observations by fitting its physical knobs.

Comparand: the 24-hour COMPOSITE diurnal cycle. The box is driven by a repeating
diurnal forcing, so day-by-day agreement is meaningless; what has to match is the
mean shape over the window.

Fitted (8):
    EMIS_NOX, EMIS_VOC, EMIS_ISO, EMIS_CO   [ppb/h]
    KDIL_DAY, KDIL_NIGHT                    [s]  (dilution TIMESCALES, see below)
    BG_O3, BG_CO                            [ppb]

BG_O3 has to be fitted, not defaulted: the stock value is 30 ppb while winter
Dallas ozone runs at 10-21 ppb, so a fixed 30 makes the observed minimum
literally unreachable and the optimiser just collapses every emission onto its
lower bound.

Bounded by observation/met:
  * KDIL_DAY is capped by tau_BLH = H/(dH/dt) from the observed boundary layer.
    For a well-mixed column, dC/dt = E/H - C*(dH/dt)/H, so entrainment alone
    dilutes on tau_BLH; any ventilation makes it faster, i.e. tau <= tau_BLH.
    The old default of 1500 s sits ~15x outside that bound.
  * TMIN/TMAX come from the observed 2 m temperature.
  * Backgrounds come from the observed nocturnal floors.

Usage: python3 _fit_anchor.py <label> [maxiter]
Writes results_fit_<label>.txt
"""
import sys, os, csv, math, statistics, json
sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box

D = "/mnt/d/atmos"
TARGETS = [("o3", "O3", 1.0), ("nox", "NOx", 1.0), ("co", "CO", 0.25), ("par", "PAR", 0.20)]


def read_obs(label):
    path = os.path.join(D, f"obs_{label}_cb05.csv")
    byhr = {col: {h: [] for h in range(24)} for col, _, _ in TARGETS}
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            h = int(r["datetime"][11:13])
            for col, key, _ in TARGETS:
                v = r.get(key)
                if v not in (None, ""):
                    byhr[col][h].append(float(v))
    obs = {}
    for col, key, _ in TARGETS:
        obs[col] = {h: (sum(v) / len(v) if v else float("nan")) for h, v in byhr[col].items()}
    # scale each target by its own daily mean so the weights mean the same thing
    scale = {col: max(statistics.mean([v for v in obs[col].values() if v == v]), 1.0)
             for col, _, _ in TARGETS}
    return obs, scale


def read_met(label):
    path = os.path.join(D, f"met_{label}.csv")
    T, H, hr = [], [], []
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh):
            T.append(float(r["temperature_2m"]))
            H.append(float(r["boundary_layer_height"]))
            hr.append(int(r["time"][11:13]))
    tmin, tmax = min(T) + 273.15, max(T) + 273.15
    taus = []
    for i in range(1, len(H) - 1):
        if hr[i] not in range(8, 15):
            continue
        dHdt = (H[i + 1] - H[i - 1]) / 7200.0
        if dHdt > 0.005 and H[i] > 50:
            taus.append(H[i] / dHdt)
    tau = statistics.median(taus) if taus else 20000.0
    return tmin, tmax, tau


def make_cost(label, obs, scale, tmin, tmax, tau, verbose=False):
    n = [0]

    def cost(x):
        (emis_nox, emis_voc, emis_iso, emis_co,
         kdil_day, kdil_night, bg_o3, bg_co) = x
        cfg = dict(EMIS_NOX=emis_nox, EMIS_VOC=emis_voc, EMIS_ISO=emis_iso,
                   EMIS_CO=emis_co, KDIL_DAY=kdil_day, KDIL_NIGHT=kdil_night,
                   BG_O3=bg_o3, BG_CO=bg_co,
                   BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
                   TMIN_K=tmin, TMAX_K=tmax)
        try:
            m = kpp_box.run_scenario(cfg, "fitrun", "anchor")
        except Exception:
            return 1e9
        traj = m["traj"]
        tot = 0.0
        for col, _, w in TARGETS:
            mod = kpp_box.diurnal(traj, col)
            s = 0.0
            cnt = 0
            for h in range(24):
                o = obs[col][h]
                mm = mod[h]
                if o != o or mm != mm:
                    continue
                s += ((mm - o) / scale[col]) ** 2
                cnt += 1
            if cnt:
                tot += w * s / cnt
        n[0] += 1
        if verbose and n[0] % 25 == 0:
            print(f"    eval {n[0]:4d}  cost={tot:.5f}  x={[round(v,1) for v in x]}", flush=True)
        return tot

    return cost


def main():
    label = sys.argv[1]
    maxiter = int(sys.argv[2]) if len(sys.argv) > 2 else 40

    obs, scale = read_obs(label)
    tmin, tmax, tau = read_met(label)
    print(f"{label}: T=[{tmin:.1f},{tmax:.1f}]K  tau_BLH={tau:.0f}s")
    for col, _, _ in TARGETS:
        prof = " ".join(f"{obs[col][h]:7.2f}" if obs[col][h] == obs[col][h] else "      ."
                        for h in range(24))
        print(f"  obs {col:<4} {prof}")

    cost = make_cost(label, obs, scale, tmin, tmax, tau, verbose=True)

    # ---- start from the BLH-consistent box, not from the over-ventilated default
    x = [12.0, 20.0, 0.5, 150.0, tau, 5.0 * 86400.0, 30.0, 150.0]
    lo = [0.5, 0.5, 0.0, 1.0, 0.1 * tau, 86400.0, 5.0, 50.0]
    hi = [60.0, 200.0, 10.0, 600.0, tau, 30.0 * 86400.0, 60.0, 400.0]
    NP = 8

    best = cost(x)
    print(f"  start cost = {best:.5f}")
    step = [(hi[i] - lo[i]) * 0.3 for i in range(NP)]

    for it in range(maxiter):
        improved = False
        for i in range(NP):
            for sgn in (+1, -1):
                y = list(x)
                y[i] = min(hi[i], max(lo[i], y[i] + sgn * step[i]))
                if y[i] == x[i]:
                    continue
                c = cost(y)
                if c < best:
                    best, x, improved = c, y, True
        if not improved:
            step = [s * 0.5 for s in step]
            if max(step[i] / max(hi[i] - lo[i], 1e-9) for i in range(NP)) < 0.01:
                print("  converged")
                break
        if it % 5 == 0:
            print(f"  iter {it:3d} cost={best:.5f} x={[round(v,2) for v in x]}", flush=True)

    print(f"FINAL cost={best:.5f}")
    print(f"FINAL x={x}")
    names = ["EMIS_NOX", "EMIS_VOC", "EMIS_ISO", "EMIS_CO",
             "KDIL_DAY", "KDIL_NIGHT", "BG_O3", "BG_CO"]
    with open(os.path.join(D, f"results_fit_{label}.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"label {label}\ncost {best:.6f}\ntau_BLH {tau:.0f} s\n\n")
        for n, v in zip(names, x):
            fh.write(f"  {n:<12} = {v:.4f}\n")
        cfg = dict(EMIS_NOX=x[0], EMIS_VOC=x[1], EMIS_ISO=x[2], EMIS_CO=x[3],
                   KDIL_DAY=x[4], KDIL_NIGHT=x[5], BG_O3=x[6], BG_CO=x[7],
                   BG_NO2=2.0, BG_NO=0.5, BG_VOC=20.0, BG_ISO=0.2,
                   TMIN_K=tmin, TMAX_K=tmax)
        m = kpp_box.run_scenario(cfg, "fitbest", "anchor")
        traj = m["traj"]
        fh.write("\n=== modelled vs observed diurnal composites ===\n")
        for col, _, _ in TARGETS:
            mod = kpp_box.diurnal(traj, col)
            fh.write(f"{col}  model " + " ".join(f"{mod[h]:8.3f}" for h in range(24)) + "\n")
            fh.write(f"{col}  obs   " + " ".join(f"{obs[col][h]:8.3f}" for h in range(24)) + "\n")
        fh.write("hour      " + " ".join(f"{h:>8d}" for h in range(24)) + "\n")
        fh.write(f"\nmetrics {json.dumps({k: v for k, v in m.items() if k.startswith('d3') or k == 'mda8'})}\n")
    print(f"wrote results_fit_{label}.txt")


if __name__ == "__main__":
    main()
