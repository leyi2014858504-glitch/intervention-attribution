import csv, os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 12
FULL = 1 << N
IDX = np.arange(FULL)
K = np.array([bin(x).count("1") for x in IDX])


def mobius(vals):
    f = np.array(vals, float)
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        f[sel] -= f[IDX[sel] ^ bit]
    return f


def zeta(f):
    g = np.array(f, float)
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        g[sel] += g[IDX[sel] ^ bit]
    return g


def load_hm(scen):
    hm = np.full((32, FULL), np.nan)
    for k in range(1, 7):
        with open(os.path.join(W, "co_hm_v_%d.csv" % k)) as fh:
            for row in csv.DictReader(fh):
                nm = row["name"]
                if not nm.startswith(scen + "_m"):
                    continue
                mask = int(nm.split("_m")[1])
                for j, h in enumerate(range(40, 72)):
                    hm[j, mask] = float(row["h%02d" % h])
    assert not np.isnan(hm).any(), scen
    return hm


def load_old(scen):
    v = np.full(FULL, np.nan)
    for k in range(1, 7):
        with open(os.path.join(W, "co_v_m8_%d.csv" % k)) as fh:
            for row in csv.DictReader(fh):
                nm = row["name"]
                if not nm.startswith(scen + "_m"):
                    continue
                mask = int(nm.split("_m")[1])
                if row["v_mda8"] not in ("", "NaN"):
                    v[mask] = float(row["v_mda8"])
    assert not np.isnan(v).any(), scen
    return v


def diag(v, tag):
    vc = v - v[0]
    m = mobius(vc)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    order = np.argsort(-a)
    cum = np.cumsum(a[order]) / l1
    top100 = 100 * cum[99]
    n50 = int(np.searchsorted(cum, 0.50) + 1)
    Mk = np.array([a[K[1:] == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    R2 = float(np.cumsum(Mk)[1] / l1)
    errs = []
    for q in range(1, N + 1):
        mq = np.where(K <= q, m, 0.0)
        errs.append(np.abs(vc - zeta(mq)).max() / np.abs(vc).max())
    return pr, top100, n50, peak, R2, 100 * max(errs), int(np.argmax(errs) + 1)


for scen in ("r_sat", "r_lim1"):
    hm = load_hm(scen)
    A = np.full((24, FULL), np.nan)
    for j, h in enumerate(range(48, 72)):
        A[j] = hm[h - 40 - 7:h - 40 + 1].mean(axis=0)
    h = A.max(axis=0)
    v5 = load_old(scen)
    print("\n%s  hourly all-on %.3f | 5-min %.3f" % (scen, h[0], v5[0]))
    print("  hourly: PR %8.1f top100 %5.1f%% N50 %4d peak k%-2d R2 %.3f err %.0f%%@q%d"
          % diag(h, scen))
    print("  5-min : PR %8.1f top100 %5.1f%% N50 %4d peak k%-2d R2 %.3f err %.0f%%@q%d"
          % diag(v5, scen))
