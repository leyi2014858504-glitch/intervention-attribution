import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 16
FULL = 1 << N


def mobius(v):
    f = np.array(v, dtype=float)
    idx = np.arange(FULL)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


POP = np.array([bin(i).count("1") for i in range(1, FULL)])


def stats(v, label):
    v = v - v[0]
    m = mobius(v)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    order = np.argsort(-a)
    cum = np.cumsum(a[order]) / l1
    top100 = 100 * cum[99]
    n50 = int(np.searchsorted(cum, 0.50) + 1)
    n90 = int(np.searchsorted(cum, 0.90) + 1)
    Mk = np.array([a[POP == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    R2 = float(np.cumsum(Mk)[1] / l1)
    A = v.reshape(1 << 8, 1 << 8)
    s = np.linalg.svd(A, compute_uv=False)
    rk = int((s > 1e-6 * s[0]).sum())
    print("  %-28s PR=%9.1f  top100=%5.1f%%  N50=%5d  N90=%6d  peak=k%-2d  R2=%.3f  rank=%d/256"
          % (label, pr, top100, n50, n90, peak, R2, rk))
    return pr


A = np.full((24, FULL), np.nan)
for k in range(1, 13):
    with open(os.path.join(W, "aw_v_dallas_2016w_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            mask = int(row["name"][1:])
            for h in range(24):
                A[h, mask] = float(row["a%02d" % h])
assert not np.isnan(A).any()

groups = {
    "Morning sum (06-11)": list(range(6, 12)),
    "Noon sum (12-14)": list(range(12, 15)),
    "Afternoon sum (15-18)": list(range(15, 19)),
    "Evening/night sum (19-05)": list(range(19, 24)) + list(range(0, 6)),
}
print("Dallas n=16: period sums of the 24 hourly A_w tables")
print("=" * 100)
for label, hs in groups.items():
    stats(A[hs].sum(axis=0), label)
print("-" * 100)
stats(A.sum(axis=0), "24h sum (00-23)")
stats(A.max(axis=0), "MDA8 hourly (max)")

v5 = np.full(FULL, np.nan)
vd = np.full(FULL, np.nan)
for k in range(1, 13):
    with open(os.path.join(W, "n16_v_dallas_2016w_m8_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            mask = int(row["name"][1:])
            v5[mask] = float(row["v_mda8"])
            vd[mask] = float(row["v_ppb"])
print("-" * 100)
stats(v5, "MDA8 (5-min)")
stats(vd, "d3_mean")
