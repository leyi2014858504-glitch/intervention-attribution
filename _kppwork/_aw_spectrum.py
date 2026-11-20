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


def stats(v, label, rank=True):
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
    rk = ""
    if rank:
        A = v.reshape(1 << 8, 1 << 8)
        s = np.linalg.svd(A, compute_uv=False)
        rk = "%d/256" % int((s > 1e-6 * s[0]).sum())
    print("  %-10s PR=%9.1f  top100=%5.1f%%  N50=%5d  N90=%6d  peak=k%-2d  R2=%.3f  rank=%s"
          % (label, pr, top100, n50, n90, peak, R2, rk))
    return dict(label=label, pr=pr, top100=top100, n50=n50, n90=n90, peak=peak, r2=R2, rank=rk)


# ---- 24 A_w tables from the sweep ----
A = np.full((24, FULL), np.nan)
for k in range(1, 13):
    with open(os.path.join(W, "aw_v_dallas_2016w_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            mask = int(row["name"][1:])
            for h in range(24):
                A[h, mask] = float(row["a%02d" % h])
assert not np.isnan(A).any(), "missing A_w values"

print("Dallas n=16: 24 hourly A_w tables (8-h means ending at 48+h)")
print("=" * 100)
res = []
for h in range(24):
    res.append(stats(A[h], "A_%02d" % h, rank=True))

print("\nkey winning windows")
print("=" * 100)
for h in (15, 16, 17, 18, 19, 20, 21, 22, 23):
    stats(A[h], "A_%02d" % h)
print("\nmost/least dense A_w by PR")
print("  min PR: %s %.1f" % (min(res, key=lambda r: r["pr"])["label"], min(r["pr"] for r in res)))
print("  max PR: %s %.1f" % (max(res, key=lambda r: r["pr"])["label"], max(r["pr"] for r in res)))

# ---- hourly MDA8 = max over the 24 A_w ----
print("\nhourly MDA8 (max of the 24 A_w)")
print("=" * 100)
hourly = A.max(axis=0)
stats(hourly, "hourly", rank=True)

# ---- 5-min MDA8 and d3 from the existing table ----
v5 = np.full(FULL, np.nan)
vd = np.full(FULL, np.nan)
for k in range(1, 13):
    with open(os.path.join(W, "n16_v_dallas_2016w_m8_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            mask = int(row["name"][1:])
            v5[mask] = float(row["v_mda8"])
            vd[mask] = float(row["v_ppb"])
assert not np.isnan(v5).any()
print("\nreference metrics (existing table)")
print("=" * 100)
stats(v5, "MDA8-5m", rank=True)
stats(vd, "d3_mean", rank=True)
