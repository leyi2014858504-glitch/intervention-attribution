"""Best-case low-rank: truncate every A_w table to rank r using the FULL
sweep (no sampling error), then recompute MDA8 and its Mobius spectrum.
If this already moves the spectral diagnostics, no sampling scheme can
recover them from a low-rank surrogate."""
import numpy as np

W = r"D:\atmos\_kppwork"
d = np.load(W + r"\aw_h_dallas_2016w.npz")
A = d["A"].astype(float)
EXACT = d["mda8_h"].astype(float)
N = 16
FULL = 1 << N
POP = np.concatenate(([0], np.array([bin(i).count("1") for i in range(1, FULL)])))


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
    top100 = 100 * cum[99]
    n50 = int(np.searchsorted(cum, 0.50) + 1)
    n90 = int(np.searchsorted(cum, 0.90) + 1)
    Mk = np.array([a[POP[1:] == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    rk = int((np.linalg.svd(v.reshape(256, 256), compute_uv=False) > 1e-6 *
              np.linalg.svd(v.reshape(256, 256), compute_uv=False)[0]).sum())
    return dict(pr=pr, top100=top100, n50=n50, n90=n90, peak=peak, rank=rk, m=m)


def shapley(m):
    phi = np.zeros(N)
    for i in range(N):
        bit = 1 << i
        idx = np.arange(FULL)
        sel = (idx & bit) > 0
        phi[i] = (m[idx[sel]] / POP[idx[sel]]).sum()
    return phi


se = stats(EXACT)
phi_e = shapley(se["m"])
print("exact: PR=%.1f top100=%.1f%% N50=%d N90=%d peak=%d rank=%d"
      % (se["pr"], se["top100"], se["n50"], se["n90"], se["peak"], se["rank"]))
print("\nrank-r SVD truncation of each A_w (full data); downstream hourly MDA8")
print(" r  max|dA|  MDA8max  PR       rank  peak  top100  max|dphi|  relL2  signflip")
for r in (1, 2, 3, 4, 6, 8, 12, 16, 24):
    rec = np.empty_like(A)
    for h in range(24):
        M = A[h].reshape(256, 256)
        U, s, Vt = np.linalg.svd(M, full_matrices=False)
        rec[h] = ((U[:, :r] * s[:r]) @ Vt[:r]).ravel()
    mda8 = rec.max(0)
    sr = stats(mda8)
    phi_r = shapley(sr["m"])
    print(" %2d %8.4f %8.4f %9.1f %5d %5d %6.1f %9.4f %6.3f %5d"
          % (r, np.abs(rec - A).max(), np.abs(mda8 - EXACT).max(), sr["pr"],
             sr["rank"], sr["peak"], sr["top100"],
             np.abs(phi_e - phi_r).max(),
             np.linalg.norm(phi_e - phi_r) / np.linalg.norm(phi_e),
             int((np.sign(phi_e) != np.sign(phi_r)).sum())))
