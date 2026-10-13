#!/bin/bash
# Top dividends at n=16 (MDA8), with member names -- what is behind the k-peak.
cd /mnt/d/atmos
python3 - <<'EOF'
import csv

def load(path, col):
    v = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            if row[col] in ("", "NaN"):
                continue
            v[int(row["name"][1:])] = float(row[col])
    return v

def mobius_fast(vals, n):
    f = dict(vals)
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                f[mask] -= f[mask ^ bit]
    return f

def name_of(mask):
    nox = [w for w in range(8) if (mask >> w) & 1]
    voc = [w for w in range(8) if (mask >> (8 + w)) & 1]
    def ws(lst):
        return "+".join("%02d-%02d" % (3 * w, 3 * w + 3) for w in lst) if lst else "-"
    return "NOx[%s] VOC[%s]" % (ws(nox), ws(voc))

for scen in ("r_sat", "r_lim1"):
    v = load("_kppwork/n16_v_%s_m8.csv" % scen, "v_mda8")
    m = mobius_fast(v, 16)
    items = sorted(((abs(x), A, x) for A, x in m.items() if A), reverse=True)[:10]
    print("=== %s (MDA8): top-10 |m(A)| dividends" % scen)
    for a, A, x in items:
        print("  k=%2d  m=%+9.3f  %s" % (bin(A).count("1"), x, name_of(A)))
    # where does the k=10 mass live: top members of the k=10 class
    for K in (2, 3, 10):
        topk = sorted(((abs(m[A]), A) for A in m if bin(A).count("1") == K), reverse=True)[:5]
        print("  -- top-5 within k=%d --" % K)
        for a, A in topk:
            print("  m=%+9.3f  %s" % (m[A], name_of(A)))
EOF