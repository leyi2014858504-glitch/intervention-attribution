#!/usr/bin/env python3
"""G3: machine-derived conservation check for the KPP-3.5 CB05 box.

Usage:  python3 g3_check.py [root_dir] [n_chemistry_reactions]

Strategy (no manual atom tables):
  1. Parse the stoichiometric matrix S from the KPP-generated
     cb05_Function.f90 (aggregate form  Vdot(i) = c*A(j) + ...): this is the
     *merged* S the solver actually integrates.
     In the dummy-forcing build only the first N_CHEM reactions are chemistry;
     the appended forcing reactions (emissions / ventilation / background
     inflow) are dropped here, so the analysis is about the closed chemical
     sub-mechanism.
  2. Basis Z of the left null space of S -> every conservation law.
  3. Verify each law has zero drift between cb05_t0.dat and cb05_tend.dat
     (the pure-chemistry driver runs with every forcing global = 0, i.e. a
     genuinely closed box).
  4. Name C-count and N-count by fitting the null space to chemical anchors
     that are identified BY NAME (species order changes between builds).
"""
import re
import sys
import numpy as np

ROOT = (sys.argv[1] if len(sys.argv) > 1 else "/root/kpp2") + "/"
N_CHEM = int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 9

# ---------- species names ----------
mtxt = open(ROOT + "cb05_Monitor.f90", errors="ignore").read()
names = [n.strip() for n in
         re.findall(r"'([^']*)'", mtxt[mtxt.index("SPC_NAMES"):mtxt.index("EQN_NAMES")])]
idx = {n: i for i, n in enumerate(names)}
print("species: %d (%s ...)" % (len(names), ", ".join(names[:3])))

# ---------- 1. parse S ----------
txt = open(ROOT + "cb05_Function.f90", errors="ignore").read()
nA = max(int(m) for m in re.findall(r"\bA\((\d+)\)\s*=", txt))
joined = re.sub(r"&\s*\n\s*", " ", txt)            # glue continuations

entries = {}
nv = 0
for m in re.finditer(r"Vdot\((\d+)\)\s*=\s*([^&\n]+)", joined):
    i = int(m.group(1))
    nv = max(nv, i)
    for sign, coef, j in re.findall(r"([+-]?)\s*([0-9][0-9.eE+-]*)?\*?A\((\d+)\)", m.group(2)):
        j = int(j)
        if j > N_CHEM:
            continue                                # drop the dummy forcing
        s = -1.0 if sign == "-" else 1.0
        if coef:
            s *= float(coef)
        entries[(j, i)] = entries.get((j, i), 0.0) + s

S = np.zeros((min(nA, N_CHEM), nv))
for (j, i), v in entries.items():
    S[j - 1, i - 1] = v
print("S: %d chemistry reactions x %d variable species, nnz=%d"
      % (S.shape[0], S.shape[1], np.count_nonzero(S)))

# the dummy pools appear in no chemistry reaction -> their rows are all zero and
# would show up as trivial "conservation laws"; drop them.
live = np.where(np.any(S != 0.0, axis=0))[0]
print("all-zero (non-chemical) species columns dropped: %d" % (nv - len(live)))
S = S[:, live]

# ---------- 2. conservation laws ----------
rank = int(np.linalg.matrix_rank(S))
U, sv, Vt = np.linalg.svd(S, full_matrices=False)
tol = 1e-9
rank = int(np.sum(sv > tol * sv[0]))
Z = Vt[rank:]                                       # right null space of S
print("rank(S)=%d  dim null(S)=%d (species-space conservation laws)"
      % (rank, Z.shape[0]))


def integrality(z):
    m = np.abs(z) > 1e-6 * np.max(np.abs(z))
    z0 = z[np.argmax(np.abs(z))]
    v = z / z0
    return float(np.max(np.abs(v[m] - np.round(v[m]))))


ints = sum(1 for k in range(Z.shape[0]) if integrality(Z[k]) < 1e-6)
print("integer-valued basis vectors: %d/%d" % (ints, Z.shape[0]))

# ---------- 3. trajectory residual ----------
def load(fn):
    ls = [l.split() for l in open(fn) if l.strip()]
    n = int(ls[0][0])
    return [l[0] for l in ls[1:1 + n]], np.array([float(l[1]) for l in ls[1:1 + n]])


_, c0 = load(ROOT + "cb05_t0.dat")
_, ce = load(ROOT + "cb05_tend.dat")
assert len(c0) == len(names), (len(c0), len(names))
v0 = c0[live]
ve = ce[live]

worst, worst_k = 0.0, -1
for k in range(Z.shape[0]):
    z = Z[k]
    a, b = float(z @ v0), float(z @ ve)
    if abs(a) < 1e12:
        continue                                    # invariant ~0 at start
    rel = abs(b - a) / abs(a)
    if rel > worst:
        worst, worst_k = rel, k
print("max relative drift over ALL conservation laws: %.3e (basis #%d)"
      % (worst, worst_k))


# ---------- 4. named element counts ----------
def make_element(target):
    rows = np.array([Z[:, idx[s]] for s in target])   # (n_anchor, n_law)
    bb = np.array(list(target.values()))
    sol, *_ = np.linalg.lstsq(rows, bb, rcond=None)
    z = sol @ Z
    res = float(np.max(np.abs(rows @ sol - bb)))
    print("  anchor fit residual: %.3e" % res)
    return z


elem = {
    "C-count": {"co": 1.0, "ch4": 1.0, "no2": 0.0, "no": 0.0, "o3": 0.0},
    "N-count": {"n2o5": 2.0, "no2": 1.0, "no": 1.0, "o3": 0.0},
}
named_ok = True
for nm, target in elem.items():
    z = make_element(target)
    a, b = float(z @ v0), float(z @ ve)
    rel = abs(b - a) / abs(a)
    named_ok &= rel < 1e-6
    print("%s: Q0=%.6e  Qend=%.6e  rel_drift=%.3e  int_err=%.3e"
          % (nm, a, b, rel, integrality(z)))

ok = worst < 1e-6 and named_ok
print("VERDICT:", "PASS (G3 < 1e-6)" if ok else "CHECK (G3 >= 1e-6)")
sys.exit(0 if ok else 1)
