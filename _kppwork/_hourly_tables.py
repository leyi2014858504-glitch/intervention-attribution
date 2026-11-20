"""Build regulatory hourly MDA8 tables from the hm_v_<tag>_K.csv shards.

hm_v stores hourly mean O3 for hours 40..71.  The 24 MDA8 windows are 8-hour
means of hourly values ending at h = 48..71:
    A_h = mean(hm[h-7 .. h])
    hourly MDA8 = max_h A_h
"""
import csv
import os
import numpy as np

W = r"D:\atmos\_kppwork"
FULL = 1 << 16
H0, H1 = 40, 71


def load_hm(tag):
    hm = np.full((H1 - H0 + 1, FULL), np.nan)
    for k in range(1, 13):
        p = os.path.join(W, "hm_v_%s_%d.csv" % (tag, k))
        with open(p) as fh:
            for row in csv.DictReader(fh):
                mask = int(row["name"][1:])
                for j, h in enumerate(range(H0, H1 + 1)):
                    hm[j, mask] = float(row["h%02d" % h])
    assert not np.isnan(hm).any(), "missing hourly means for %s" % tag
    return hm


def aw_from_hm(hm):
    A = np.full((24, FULL), np.nan)
    for j, h in enumerate(range(48, 72)):
        A[j] = hm[h - H0 - 7:h - H0 + 1].mean(axis=0)
    assert not np.isnan(A).any()
    return A


def load_old(tag, col):
    v = np.full(FULL, np.nan)
    p = os.path.join(W, "n16_v_%s_m8.csv" % tag)
    files = [p] if os.path.exists(p) else [
        os.path.join(W, "n16_v_%s_m8_%d.csv" % (tag, k)) for k in range(1, 13)]
    for p in files:
        with open(p) as fh:
            for row in csv.DictReader(fh):
                val = row.get(col, "")
                if val not in ("", "NaN"):
                    v[int(row["name"][1:])] = float(val)
    assert not np.isnan(v).any(), "missing old column %s for %s" % (col, tag)
    return v


for tag in ("dallas_2016w", "r_sat", "r_lim1", "deerpark_2017s"):
    hm = load_hm(tag)
    A = aw_from_hm(hm)
    h = A.max(axis=0)
    v5 = load_old(tag, "v_mda8")
    vd = load_old(tag, "v_ppb")
    np.savez_compressed(os.path.join(W, "aw_h_%s.npz" % tag), hm=hm, A=A, mda8_h=h)

    with open(os.path.join(W, "n16_h_%s.csv" % tag), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "v_mda8_h", "v_mda8_5", "v_ppb", "nfail", "rc", "nrow"])
        for mask in range(FULL):
            w.writerow(["m%05d" % mask, "%.6f" % h[mask],
                        "%.6f" % v5[mask], "%.6f" % vd[mask], 0, 0, 864])

    # same header as the 5-min tables, so existing loaders only need a path swap
    shards = [[] for _ in range(12)]
    with open(os.path.join(W, "n16_vh_%s_m8.csv" % tag), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "v_ppb", "v_mda8", "nfail", "rc", "nrow"])
        for mask in range(FULL):
            w.writerow(["m%05d" % mask, "%.6f" % vd[mask],
                        "%.6f" % h[mask], 0, 0, 864])
            shards[mask % 12].append(mask)
    for k, masks in enumerate(shards, 1):
        with open(os.path.join(W, "n16_vh_%s_m8_%d.csv" % (tag, k)), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["name", "v_ppb", "v_mda8", "nfail", "rc", "nrow"])
            for mask in masks:
                w.writerow(["m%05d" % mask, "%.6f" % vd[mask],
                            "%.6f" % h[mask], 0, 0, 864])

    d = h - v5
    print("%-18s hourly all-on %.4f | 5-min %.4f | max diff %.3f | mean diff %.3f | >0.1 ppb: %d"
          % (tag, h[0], v5[0], d.max(), d.mean(), int((np.abs(d) > 0.1).sum())))
