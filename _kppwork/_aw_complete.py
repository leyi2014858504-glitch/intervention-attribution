"""Minimal validation: can a low-rank / additive reconstruction of the 24
hourly A_w tables replace the full 2^16 sweep for the MDA8 spectrum?

Data: aw_h_dallas_2016w.npz  (A: 24 x 65536, mda8_h = max_h A)
Game: 16 players = 8 NOx windows (bits 0..7) x 8 VOC windows (bits 8..15).
Each A_w(h) is a 256x256 table  M[voc, nox].

Designs compared (all priced in box-model runs; one run returns all 24 h):
  cross  : all 256 NOx-only + all 256 VOC-only coalitions (511 runs),
           additive (rank-2) fit  A = c + f(nox) + g(voc)
  random : m uniformly sampled coalitions + empty set, same additive fit
           (least squares, 511 unknowns)
Reference: exact 65,536-run sweep.
"""
import os
import numpy as np

W = r"D:\atmos\_kppwork"
OUT = os.path.join(W, "results_aw_complete.txt")
N = 16
FULL = 1 << N
POP = np.array([bin(i).count("1") for i in range(1, FULL)])
POPF = np.concatenate(([0], POP))

d = np.load(os.path.join(W, "aw_h_dallas_2016w.npz"))
A = d["A"].astype(float)          # 24 x 65536
EXACT = d["mda8_h"].astype(float)  # 65536


def mobius(v):
    f = np.array(v, dtype=float)
    idx = np.arange(FULL)
    for i in range(N):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
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
    Mk = np.array([a[POP == k].sum() for k in range(1, N + 1)])
    peak = int(np.argmax(Mk) + 1)
    R2 = float(np.cumsum(Mk)[1] / l1)
    A2 = v.reshape(256, 256)
    s = np.linalg.svd(A2, compute_uv=False)
    rk = int((s > 1e-6 * s[0]).sum())
    return dict(label=label, pr=pr, top100=top100, n50=n50, n90=n90,
                peak=peak, r2=R2, rank=rk, m=m, l1=l1,
                m2=a[POP == 2].sum() / l1)


def shapley(m):
    phi = np.zeros(N)
    for i in range(N):
        bit = 1 << i
        idx = np.arange(FULL)
        sel = (idx & bit) > 0
        phi[i] = (m[idx[sel]] / POPF[idx[sel]]).sum()
    return phi


def fit_additive(obs_mask, vals):
    """Least squares for A[y,x] = c + f_x + g_y, f_0 = g_0 = 0."""
    m = len(obs_mask)
    X = np.zeros((m, 1 + 255 + 255))
    X[:, 0] = 1.0
    x = obs_mask & 255
    y = obs_mask >> 8
    rows = np.arange(m)
    fx = x > 0
    gy = y > 0
    X[rows[fx], x[fx]] = 1.0
    X[rows[gy], 255 + y[gy]] = 1.0
    coef, *_ = np.linalg.lstsq(X, vals, rcond=None)
    c = coef[0]
    f = np.concatenate(([0.0], coef[1:256]))
    g = np.concatenate(([0.0], coef[256:]))
    return c + f[np.arange(FULL) & 255] + g[np.arange(FULL) >> 8]


def fit_additive_cross(Ah):
    M = Ah.reshape(256, 256)            # [voc, nox]
    c = M[0, 0]
    f = M[0, :] - c                     # nox profile
    g = M[:, 0] - c                     # voc profile
    return (c + f[None, :] + g[:, None]).ravel()


def evaluate(Arec, tag, lines):
    err = np.abs(Arec - A)
    mda8 = Arec.max(axis=0)
    dmax = np.abs(mda8 - EXACT).max()
    drms = np.sqrt(((mda8 - EXACT) ** 2).mean())
    se = stats(EXACT, "exact")
    sr = stats(mda8, "recon")
    dm = sr["m"] - se["m"]
    phi_e = shapley(se["m"])
    phi_r = shapley(sr["m"])
    signflip = int((np.sign(phi_e) != np.sign(phi_r)).sum())
    line = ("  %-10s max|dA|=%7.4f  MDAdmax=%7.4f  dRMSE=%7.4f  "
            "PR %9.1f -> %9.1f  rank %2d -> %2d  peak %2d -> %2d  "
            "top100 %5.1f -> %5.1f  m2share %5.1f%% -> %5.1f%%  "
            "max|dphi|=%6.3f  relL2=%5.3f  signflip=%d"
            % (tag, err.max(), dmax, drms,
               se["pr"], sr["pr"], se["rank"], sr["rank"],
               se["peak"], sr["peak"], se["top100"], sr["top100"],
               se["m2"] * 100, sr["m2"] * 100,
               np.abs(phi_e - phi_r).max(),
               np.linalg.norm(phi_e - phi_r) / np.linalg.norm(phi_e),
               signflip))
    lines.append(line)
    print(line)


def main():
    lines = []
    print("A_w intrinsic compressibility (max abs residual, ppb)")
    print("hour   r=1    r=2    r=3    r=4    r=8")
    for h in range(24):
        M = A[h].reshape(256, 256)
        U, s, Vt = np.linalg.svd(M, full_matrices=False)
        row = "%02d  " % h
        for r in (1, 2, 3, 4, 8):
            Xr = (U[:, :r] * s[:r]) @ Vt[:r]
            row += "%7.4f" % np.abs(Xr - M).max()
        print(row)
        lines.append(row)

    se = stats(EXACT, "exact")
    print("\nexact hourly MDA8: PR=%.1f top100=%.1f%% N50=%d N90=%d peak=%d "
          "rank=%d/256 m2share=%.1f%%"
          % (se["pr"], se["top100"], se["n50"], se["n90"], se["peak"],
             se["rank"], se["m2"] * 100))

    print("\nstructured design: 511 cross coalitions (all NOx-only + VOC-only)")
    Arec = np.array([fit_additive_cross(A[h]) for h in range(24)])
    evaluate(Arec, "cross", lines)

    print("\nrandom designs + additive least squares")
    for m in (128, 256, 512, 1024, 2048, 4096):
        for seed in range(3):
            rng = np.random.default_rng(seed)
            pick = np.concatenate(([0], rng.choice(FULL, m - 1, replace=False)))
            Arec = np.array([fit_additive(pick, A[h][pick]) for h in range(24)])
            evaluate(Arec, "m=%d/s%d" % (m, seed), lines)

    with open(OUT, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
