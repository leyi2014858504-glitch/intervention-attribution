import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 16
FULL = 1 << N
POP = np.array([bin(i).count("1") for i in range(FULL)])
POP_NON = POP[1:]
F = np.array([0.071, 0.078, 0.220, 0.107, 0.072, 0.156, 0.216, 0.080])


def mobius(v):
    f = np.array(v, dtype=float)
    idx = np.arange(FULL)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def zeta(f):
    f = np.array(f, dtype=float)
    idx = np.arange(FULL)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] += f[idx[sel] ^ bit]
    return f


def load(tag, col, hourly=True):
    v = np.full(FULL, np.nan)
    p = os.path.join(W, ("n16_vh_%s_m8.csv" if hourly else "n16_v_%s_m8.csv") % tag)
    files = [p] if os.path.exists(p) else [
        os.path.join(W, ("n16_vh_%s_m8_%d.csv" if hourly else "n16_v_%s_m8_%d.csv") % (tag, k))
        for k in range(1, 13)]
    for p in files:
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    assert not np.isnan(v).any(), (tag, col)
    return v


def summary(v, tag):
    v = v - v[0]
    m = mobius(v)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    order = np.argsort(-a)
    cum = np.cumsum(a[order]) / l1
    top100 = 100 * cum[99]
    n50 = int(np.searchsorted(cum, 0.50) + 1)
    Mk = np.array([a[POP_NON == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    rk = int((np.linalg.svd(v.reshape(256, 256), compute_uv=False) > 1e-6 * np.linalg.svd(v.reshape(256, 256), compute_uv=False)[0]).sum())
    errs = []
    for q in range(1, N + 1):
        mt = m.copy()
        mt[POP > q] = 0.0
        errs.append(np.max(np.abs(zeta(mt) - v)) / np.max(np.abs(v)))
    e_peak = 100 * max(errs); e_q = int(np.argmax(errs) + 1)
    singles = np.array([m[1 << i] for i in range(N)])
    phi = np.zeros(N)
    for A in range(1, FULL):
        k = bin(A).count("1")
        for i in range(N):
            if (A >> i) & 1:
                phi[i] += m[A] / k
    rn = np.zeros(FULL); rv = np.zeros(FULL)
    idx = np.arange(FULL)
    for w in range(8):
        rn[(idx & (1 << w)) > 0] += F[w]
        rv[(idx & (1 << (8 + w))) > 0] += F[w]
    tight = (rn <= 0.10) & (rv <= 0.35)
    loose = (rn <= 0.20) & (rv <= 0.50)
    bt = int(np.argmin(np.where(tight, v, 1e9)))
    bl = int(np.argmin(np.where(loose, v, 1e9)))
    print("%-18s PR=%9.1f top100=%5.1f%% N50=%4d peak=k%-2d rank=%3d err=%.1f%%@q%d | NOx+ %d/8 VOC+ %d/8 | phiNOx %+.2f phiVOC %+.2f | tight %+.2f loose %+.2f"
          % (tag, pr, top100, n50, peak, rk, e_peak, e_q,
             int((singles[:8] > 0).sum()), int((singles[8:] > 0).sum()),
             phi[:8].sum(), phi[8:].sum(), v[bt], v[bl]))
    return singles


for tag in ("dallas_2016w", "r_sat", "r_lim1", "deerpark_2017s"):
    print("\n" + tag)
    print("  hourly:")
    summary(load(tag, "v_mda8", True), tag)
    print("  5-min :")
    summary(load(tag, "v_mda8", False), tag)
