import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
N = 16
FULL = 1 << N
IDX = np.arange(FULL)
K = np.array([bin(x).count("1") for x in IDX])


def mobius(v, n):
    f = np.array(v, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        f[sel] -= f[idx[sel] ^ bit]
    return f


def zeta(f, n):
    g = np.array(f, dtype=float)
    idx = np.arange(1 << n)
    for i in range(n):
        bit = 1 << i
        sel = (idx & bit) > 0
        g[sel] += g[idx[sel] ^ bit]
    return g


def load16(tag, col):
    v = np.zeros(FULL)
    p = os.path.join(W, "n16_vh_%s_m8.csv" % tag)
    files = [p] if os.path.exists(p) else [
        os.path.join(W, "n16_vh_%s_m8_%d.csv" % (tag, k)) for k in range(1, 13)]
    for p in files:
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    assert not np.isnan(v).any(), (tag, col)
    return v


def sub8(v16):
    v8 = np.zeros(256)
    for m8 in range(256):
        m16 = 0
        for j in range(8):
            if (m8 >> j) & 1:
                m16 |= (1 << (2 * j)) | (1 << (2 * j + 1))
        v8[m8] = v16[m16]
    return v8


def load12(scen, col=None):
    hm = np.full((32, 1 << 12), np.nan)
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


def diag(v, n, label):
    vc = v - v[0]
    m = mobius(vc, n)
    a = np.abs(m[1:])
    l1 = a.sum()
    pr = l1 ** 2 / (a ** 2).sum()
    order = np.argsort(-a)
    cum = np.cumsum(a[order]) / l1
    top100 = 100 * cum[99]
    n50 = int(np.searchsorted(cum, 0.50) + 1)
    n90 = int(np.searchsorted(cum, 0.90) + 1)
    n5 = int((a >= 0.05 * a.max()).sum())
    pc = np.array([bin(i).count("1") for i in range(1, 1 << n)])
    Mk = np.array([a[pc == k].sum() for k in range(1, n + 1)])
    R = np.cumsum(Mk) / l1
    peak = int(np.argmax(Mk) + 1)
    A = vc.reshape(1 << (n // 2), 1 << (n // 2))
    s = np.linalg.svd(A, compute_uv=False)
    rk = int((s > 1e-6 * s[0]).sum())
    errs = []
    for q in range(1, n + 1):
        mq = np.where(np.array([bin(x).count("1") for x in range(1 << n)]) <= q, m, 0.0)
        errs.append(np.abs(vc - zeta(mq, n)).max() / np.abs(vc).max())
    return dict(pr=pr, top100=top100, n50=n50, n90=n90, n5=n5, Mk=Mk, R=R,
                peak=peak, rank=rk, errs=np.array(errs), m=m, vc=vc)


print("=" * 90)
print("TABLE 1: cumulative R1/R2/R3/R5 for synthetic games (hourly MDA8)")
for scen in ("r_sat", "r_lim1"):
    v16 = load16(scen, "v_mda8")
    for n, v in ((8, sub8(v16)), (12, load12(scen, "h40")), (16, v16)):
        if n == 12:
            hm = v
            A = np.full((24, 1 << 12), np.nan)
            for j, h in enumerate(range(48, 72)):
                A[j] = hm[h - 40 - 7:h - 40 + 1].mean(axis=0)
            v = A.max(axis=0)
        d = diag(v, n, scen)
        print("  %-8s n=%2d  R1 %.3f R2 %.3f R3 %.3f R5 %.3f" % (scen, n, d["R"][0], d["R"][1], d["R"][2], d["R"][4]))

print("\nSPECTRUM PEAKS and mass shares")
for scen in ("r_sat", "r_lim1"):
    v16 = load16(scen, "v_mda8")
    for n, v in ((8, sub8(v16)), (16, v16)):
        d = diag(v, n, scen)
        print("  %-8s n=%2d  peak k=%d  M(k)/M1=%.2f  M3-5 share=%.3f  R2=%.3f" %
              (scen, n, d["peak"], d["Mk"][d["peak"] - 1] / d["Mk"][0], (d["R"][4] - d["R"][1]), d["R"][1]))

print("\nCITY SPECTRA")
for lab in ("dallas_2016w", "deerpark_2017s"):
    v = load16(lab, "v_mda8")
    d = diag(v, 16, lab)
    k = d["peak"]
    m = d["m"]
    a = np.abs(m[1:])
    pc = np.array([bin(i).count("1") for i in range(1, FULL)])
    nox = np.zeros(16); voc = np.zeros(16); mixed = 0.0
    for A_mask in range(1, FULL):
        if pc[A_mask - 1] != k:
            continue
        has_n = any((A_mask >> i) & 1 for i in range(8))
        has_v = any((A_mask >> (8 + i)) & 1 for i in range(8))
        if has_n and has_v:
            mixed += abs(m[A_mask])
        elif has_n:
            nox[k - 1] += abs(m[A_mask])
        else:
            voc[k - 1] += abs(m[A_mask])
    print("  %s  M1=%.3f  peak k=%d  M(k)/M1=%.2f  R2=%.4f%%" %
          (lab, d["Mk"][0], k, d["Mk"][k - 1] / d["Mk"][0], 100 * d["R"][1]))
    print("    at k=%d: NOx-only %.1f  VOC-only %.1f  mixed %.1f" % (k, nox[k - 1], voc[k - 1], mixed))
    singles = np.array([m[1 << i] for i in range(16)])
    print("    NOx singles:", " ".join("%+.3f" % x for x in singles[:8]))
    print("    VOC singles:", " ".join("%+.3f" % x for x in singles[8:]))

print("\nTRUNCATION (centered err_v, from regenerated results_errv_full.txt)")
for lab in ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"):
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3")):
        v = load16(lab, col)
        d = diag(v, 16, lab)
        e = d["errs"]
        pk = int(np.argmax(e) + 1)
        print("  %-16s %-5s peak %.1f%%@q%d  q10 %.1f%%  min(q4-12) %.1f%%" %
              (lab, tag, 100 * e[pk - 1], pk, 100 * e[9], 100 * e[3:12].min()))

print("\nD3 vs MDA8 at the error peak order")
for lab in ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"):
    vm = load16(lab, "v_mda8")
    dm = load16(lab, "v_ppb")
    em = diag(vm, 16, lab)["errs"]
    q = int(np.argmax(em) + 1)
    ed = diag(dm, 16, lab)["errs"]
    ratio = em[q - 1] / max(ed[q - 1], 1e-12)
    print("  %-16s q=%d  err_MDA8 %.1f%%  err_d3 %.1f%%  factor %.1f" %
          (lab, q, 100 * em[q - 1], 100 * ed[q - 1], ratio))

print("\nMASS RATIOS at spectrum peaks (MDA8/d3)")
for lab in ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"):
    vm = load16(lab, "v_mda8")
    dm = load16(lab, "v_ppb")
    dm_diag = diag(vm, 16, lab)
    dd_diag = diag(dm, 16, lab)
    k = dm_diag["peak"]
    rm = dm_diag["Mk"][k - 1] / dm_diag["Mk"][0]
    rd = dd_diag["Mk"][k - 1] / dd_diag["Mk"][0]
    print("  %-16s k=%d  MDA8/M1 %.2f  d3/M1 %.2f  ratio %.1f" % (lab, k, rm, rd, rm / max(rd, 1e-12)))

print("\nCOMPRESS TABLE (hourly MDA8 and d3)")
for lab in ("r_sat", "r_lim1", "dallas_2016w", "deerpark_2017s"):
    for col, tag in (("v_mda8", "MDA8"), ("v_ppb", "d3")):
        d = diag(load16(lab, col), 16, lab)
        print("  %-16s %-5s PR=%9.1f top100=%5.1f%% N>5%%=%4d N50=%5d N90=%6d rank=%d/256" %
              (lab, tag, d["pr"], d["top100"], d["n5"], d["n50"], d["n90"], d["rank"]))

print("\nA_w PERIOD SUMS (Dallas)")
z = np.load(os.path.join(W, "aw_h_dallas_2016w.npz"))
A = z["A"]
groups = {"Morning 06-11": list(range(6, 12)), "Noon 12-14": list(range(12, 15)),
          "Afternoon 15-18": list(range(15, 19)), "Evening/night 19-05": list(range(19, 24)) + list(range(0, 6))}
for label, hs in groups.items():
    d = diag(A[hs].sum(axis=0), 16, label)
    print("  %-24s PR=%8.1f top100=%5.1f%% N50=%4d N90=%5d peak=k%-2d rank=%d/256" %
          (label, d["pr"], d["top100"], d["n50"], d["n90"], d["peak"], d["rank"]))
d = diag(A.sum(axis=0), 16, "24h")
print("  %-24s PR=%8.1f top100=%5.1f%% N50=%4d N90=%5d peak=k%-2d rank=%d/256" %
      ("24h sum", d["pr"], d["top100"], d["n50"], d["n90"], d["peak"], d["rank"]))
d = diag(A.max(axis=0), 16, "hourly")
print("  %-24s PR=%8.1f top100=%5.1f%% N50=%4d N90=%5d peak=k%-2d rank=%d/256" %
      ("MDA8 hourly", d["pr"], d["top100"], d["n50"], d["n90"], d["peak"], d["rank"]))
