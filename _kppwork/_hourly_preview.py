import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 16
FULL = 1 << N
POP = np.array([bin(i).count("1") for i in range(FULL)])
POP_NON = POP[1:]


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
    Mk = np.array([a[POP_NON == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    A = v.reshape(256, 256)
    s = np.linalg.svd(A, compute_uv=False)
    rk = int((s > 1e-6 * s[0]).sum())
    print("%-12s PR=%9.1f top100=%5.1f%% N50=%5d N90=%6d peak=k%-2d rank=%d/256"
          % (label, pr, top100, n50, n90, peak, rk))
    return m


A = np.full((24, FULL), np.nan)
for k in range(1, 13):
    with open(os.path.join(W, "aw_v_dallas_2016w_%d.csv" % k)) as fh:
        for row in csv.DictReader(fh):
            mask = int(row["name"][1:])
            for h in range(24):
                A[h, mask] = float(row["a%02d" % h])
assert not np.isnan(A).any()

hourly = A.max(axis=0)
print("Dallas hourly MDA8 preview")
print("=" * 80)
m = stats(hourly, "hourly")

# err_v(q)
qerr = []
for q in range(1, N + 1):
    mt = m.copy()
    mt[POP > q] = 0.0
    vq = zeta(mt)
    qerr.append(np.max(np.abs(vq - (hourly - hourly[0]))) / np.max(np.abs(hourly - hourly[0])))
print("err_v(q) peak %.1f%% at q=%d" % (100 * max(qerr), int(np.argmax(qerr) + 1)))

# Shapley and trunc(q)
phi = np.zeros(N)
for A_mask in range(1, FULL):
    k = bin(A_mask).count("1")
    for i in range(N):
        if (A_mask >> i) & 1:
            phi[i] += m[A_mask] / k
print("Shapley sum %.4f vs v(N) %.4f" % (phi.sum(), hourly[-1] - hourly[0]))
print("NOx Shapley:", " ".join("%+.3f" % x for x in phi[:8]))
print("VOC Shapley:", " ".join("%+.3f" % x for x in phi[8:]))

best = 0; bq = None
for q in range(1, N + 1):
    mt = m.copy()
    mt[POP > q] = 0.0
    phiq = np.zeros(N)
    for A_mask in range(1, FULL):
        if mt[A_mask] == 0:
            continue
        k = bin(A_mask).count("1")
        for i in range(N):
            if (A_mask >> i) & 1:
                phiq[i] += mt[A_mask] / k
    te = np.max(np.abs(phi - phiq)) / np.max(np.abs(phi))
    if te > best:
        best, bq = te, q
print("trunc(q) peak %.1f%% at q=%d" % (100 * best, bq))

# singles
singles = np.array([m[1 << i] for i in range(N)])
print("NOx singles:", " ".join("%+.3f" % x for x in singles[:8]))
print("VOC singles:", " ".join("%+.3f" % x for x in singles[8:]))
print("NOx positive %d/8, VOC positive %d/8" % ((singles[:8] > 0).sum(), (singles[8:] > 0).sum()))

# budget scan
f = np.array([0.071, 0.078, 0.220, 0.107, 0.072, 0.156, 0.216, 0.080])
v = hourly - hourly[0]
best_tight = (1e9, None, None, None)
best_loose = (1e9, None, None, None)
for mask in range(FULL):
    rn = f[[w for w in range(8) if (mask >> w) & 1]].sum() if mask & 0xFF else 0.0
    rv = f[[w for w in range(8) if (mask >> (8 + w)) & 1]].sum() if mask & 0xFF00 else 0.0
    if rn <= 0.10 and rv <= 0.35 and v[mask] < best_tight[0]:
        best_tight = (v[mask], rn, rv, mask)
    if rn <= 0.20 and rv <= 0.50 and v[mask] < best_loose[0]:
        best_loose = (v[mask], rn, rv, mask)
print("tight : dMDA8 %.3f  R_N %.1f%% R_V %.1f%%  mask %d" % (best_tight[0], 100 * best_tight[1], 100 * best_tight[2], best_tight[3]))
print("loose : dMDA8 %.3f  R_N %.1f%% R_V %.1f%%  mask %d" % (best_loose[0], 100 * best_loose[1], 100 * best_loose[2], best_loose[3]))
