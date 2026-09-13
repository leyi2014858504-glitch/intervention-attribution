#!/usr/bin/env python3
"""
Can the box reach the observations at all? A 2-D scan, not an optimiser.

Before spending half an hour in an optimiser it is worth knowing whether the
attainable region even contains the observed diurnal cycle. This scans the two
emission scales that dominate ozone production and reports, for each point, the
modelled minus observed error on O3 / NOx / PAR / CO.

Usage: python3 _scan_anchor.py <label>
"""
import sys, os, csv, statistics
sys.path.insert(0, "/mnt/d/atmos/_kppwork")
import kpp_box

D = "/mnt/d/atmos"
TARGETS = [("o3", "O3", 1.0), ("nox", "NOx", 1.0), ("co", "CO", 0.25), ("par", "PAR", 0.20)]


def read_obs(label):
    byhr = {col: {h: [] for h in range(24)} for col, _, _ in TARGETS}
    with open(os.path.join(D, f"obs_{label}_cb05.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            h = int(r["datetime"][11:13])
            for col, key, _ in TARGETS:
                v = r.get(key)
                if v not in (None, ""):
                    byhr[col][h].append(float(v))
    obs = {c: {h: (sum(v) / len(v) if v else float("nan")) for h, v in byhr[c].items()}
           for c, _, _ in TARGETS}
    scale = {c: max(statistics.mean([v for v in obs[c].values() if v == v]), 1.0)
             for c, _, _ in TARGETS}
    return obs, scale


def read_met(label):
    T, H, hr = [], [], []
    with open(os.path.join(D, f"met_{label}.csv"), newline="") as fh:
        for r in csv.DictReader(fh):
            T.append(float(r["temperature_2m"]))
            H.append(float(r["boundary_layer_height"]))
            hr.append(int(r["time"][11:13]))
    taus = []
    for i in range(1, len(H) - 1):
        if hr[i] not in range(8, 15):
            continue
        dHdt = (H[i + 1] - H[i - 1]) / 7200.0
        if dHdt > 0.005 and H[i] > 50:
            taus.append(H[i] / dHdt)
    return min(T) + 273.15, max(T) + 273.15, (statistics.median(taus) if taus else 20000.0)


def main():
    label = sys.argv[1]
    obs, scale = read_obs(label)
    tmin, tmax, tau = read_met(label)
    print(f"{label}: tau_BLH={tau:.0f}s  obs daily means: "
          + ", ".join(f"{k}={scale[col]:.1f}" for col, k, _ in TARGETS))

    NOXS = [1.0, 3.0, 6.0, 10.0, 16.0]
    VOCS = [1.0, 3.0, 8.0, 20.0, 50.0]
    print(f"\n{'ENOx':>6}{'EVOC':>7}   " +
          "".join(f"{k+'_err':>10}" for _, k, _ in TARGETS) +
          f"{'NOx_ss':>9}{'O3max':>8}{'cost':>9}")
    for enox in NOXS:
        for evoc in VOCS:
            cfg = dict(EMIS_NOX=enox, EMIS_VOC=evoc, EMIS_ISO=0.3, EMIS_CO=60.0,
                       KDIL_DAY=tau, KDIL_NIGHT=86400.0,
                       BG_O3=14.0, BG_CO=200.0, BG_NO2=2.0, BG_NO=0.5,
                       BG_VOC=20.0, BG_ISO=0.2, TMIN_K=tmin, TMAX_K=tmax)
            try:
                m = kpp_box.run_scenario(cfg, "scanrun", "anchor")
            except Exception as e:
                print(f"{enox:>6}{evoc:>7}   FAILED {e}")
                continue
            traj = m["traj"]
            errs, cost = [], 0.0
            for col, _, w in TARGETS:
                mod = kpp_box.diurnal(traj, col)
                s, cnt = 0.0, 0
                for h in range(24):
                    o, mm = obs[col][h], mod[h]
                    if o != o or mm != mm:
                        continue
                    s += ((mm - o) / scale[col]) ** 2
                    cnt += 1
                e = (s / cnt) ** 0.5 if cnt else float("nan")
                errs.append(e)
                cost += w * s / cnt
            nox_ss = kpp_box.diurnal(traj, "nox")[5]
            o3max = max(kpp_box.diurnal(traj, "o3").values())
            print(f"{enox:>6.0f}{evoc:>7.0f}   " +
                  "".join(f"{e:>10.2f}" for e in errs) +
                  f"{nox_ss:>9.1f}{o3max:>8.1f}{cost:>9.3f}")


if __name__ == "__main__":
    main()
