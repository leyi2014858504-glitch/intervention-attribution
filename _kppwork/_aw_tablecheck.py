"""Verify every number in Table tab:aw against the npz."""
import numpy as np

W = r"D:\atmos\_kppwork"
d = np.load(W + r"\aw_h_dallas_2016w.npz")
A = d["A"].astype(float)
EXACT = d["mda8_h"].astype(float)
N = 16
FULL = 1 << N
POP = np.array([bin(i).count("1") for i in range(1, FULL)])

PERIODS = {
    "Morning 06-11": list(range(6, 12)),
    "Noon 12-14": list(range(12, 15)),
    "Afternoon 15-18": list(range(15, 19)),
    "Evening/night 19-05": list(range(19, 24)) + list(range(0, 6)),
    "24h sum": list(range(24)),
}


def mobius(v):
    f = np.array(v, dtype=float)
    idx = np.arange(FULL)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def stats(v):
    v = v - v[0]
    m = mobius(v)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    order = np.argsort(-a)
    cum = np.cumsum(a[order]) / l1
    Mk = np.array([a[POP == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    sv = np.linalg.svd(v.reshape(256, 256), compute_uv=False)
    return pr, 100 * cum[99], int(np.searchsorted(cum, 0.50) + 1), \
        int(np.searchsorted(cum, 0.90) + 1), peak, int((sv > 1e-6 * sv[0]).sum())


print("%-22s %10s %8s %7s %8s %5s %8s" % ("quantity", "PR", "top100", "N50", "N90", "peak", "rank"))
for name, hs in PERIODS.items():
    v = A[hs].sum(axis=0)
    print("%-22s %10.1f %7.1f%% %7d %8d %5d %8d" % ((name,) + stats(v)))
print("%-22s %10.1f %7.1f%% %7d %8d %5d %8d" % (("MDA8",) + stats(EXACT)))
print("%-22s %10.1f %7.1f%% %7d %8d %5d %8d"
      % (("d3_mean (24h mean)",) + stats(A.mean(axis=0))))
