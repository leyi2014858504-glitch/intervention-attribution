"""Sparse (low-order Mobius) reconstruction of each A_w, then max:
what is the MDA8 error?  This is the 'smooth the windows first' route."""
import numpy as np

W = r"D:\atmos\_kppwork"
d = np.load(W + r"\aw_h_dallas_2016w.npz")
A = d["A"].astype(float)
EXACT = d["mda8_h"].astype(float)
N = 16
FULL = 1 << N
POP = np.array([bin(i).count("1") for i in range(1, FULL)])
IDX = np.arange(FULL)


def mobius(v):
    f = np.array(v - v[0], dtype=float)
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        f[sel] -= f[IDX[sel] ^ bit]
    return f


def zeta(m):
    f = m.copy()
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        f[sel] += f[IDX[sel] ^ bit]
    return f


def stats(v):
    v = v - v[0]
    m = mobius(v)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    sv = np.linalg.svd(v.reshape(256, 256), compute_uv=False)
    return pr, int((sv > 1e-6 * sv[0]).sum())


print("direct MDA8 truncation (same convention)")
for q in range(1, 9):
    m = mobius(EXACT)
    m[1:][POP > q] = 0.0
    rec = zeta(m) + EXACT[0]
    print("%d  %10.4f  %10.4f"
          % (q, np.abs(rec - EXACT).max(),
             np.sqrt(((rec - EXACT) ** 2).mean())))

def shapley(m):
    phi = np.zeros(N)
    for i in range(N):
        bit = 1 << i
        sel = (IDX & bit) > 0
        phi[i] = (m[sel] / POP[IDX[sel] - 1]).sum()
    return phi


phi_e = shapley(mobius(EXACT))

print("\nper-window truncation then max")
print("q   MDA8 max err  err/max|v|  RMSE     PR        rank  "
      "dphi_max  relL2  flip")
for q in range(1, 9):
    rec = np.empty_like(A)
    for h in range(24):
        m = mobius(A[h])
        keep = POP <= q
        m[1:][~keep] = 0.0
        rec[h] = zeta(m) + A[h, 0]
    mda8 = rec.max(0)
    pr, rk = stats(mda8)
    phi_r = shapley(mm := mobius(mda8))
    print("%d  %10.4f  %14.2f%% %8.4f  %9.1f  %6d  %8.4f  %6.3f  %4d"
          % (q, np.abs(mda8 - EXACT).max(),
             100 * np.abs(mda8 - EXACT).max() / np.abs(EXACT).max(),
             np.sqrt(((mda8 - EXACT) ** 2).mean()), pr, rk,
             np.abs(phi_e - phi_r).max(),
             np.linalg.norm(phi_e - phi_r) / np.linalg.norm(phi_e),
             int((np.sign(phi_e) != np.sign(phi_r)).sum())))
