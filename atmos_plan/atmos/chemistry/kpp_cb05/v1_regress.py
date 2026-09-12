"""V1 structural regression: the chemistry must be untouched by the refactor.

The dummy-forcing mechanism appends reactions R192..R304 to the original 191.
If the refactor is clean, then in the newly generated cb05_Function.f90
  (a) every rate-constant definition A(1..191) is identical, and
  (b) every Vdot line, after dropping the terms A(j) with j > 191, is identical.

Species indices are NOT comparable between the two builds: KPP reorders the
species (sparse-LU fill-in) as soon as the reaction list grows.  Both sides are
therefore rewritten from index space into NAME space first, using the SPC_NAMES
table of the corresponding build.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "src")
NCHEM = 191


def species_names(path):
    txt = open(path).read()
    i = txt.index("SPC_NAMES")
    j = txt.index("EQN_NAMES")
    return [n.strip() for n in re.findall(r"'([^']*)'", txt[i:j])]


def params(path):
    t = open(path).read()
    return (int(re.search(r"NVAR\s*=\s*(\d+)", t).group(1)),
            int(re.search(r"NSPEC\s*=\s*(\d+)", t).group(1)),
            int(re.search(r"NREACT\s*=\s*(\d+)", t).group(1)))


def build(tag):
    """-> (nvar, nspec, nreact, A{j: expr}, Vdot{name: expr}) in NAME space."""
    nvar, nspec, nreact = params(os.path.join(SRC, tag + "_Parameters.f90"))
    names = species_names(os.path.join(SRC, tag + "_Monitor.f90"))
    assert len(names) == nspec, (tag, len(names), nspec)
    txt = re.sub(r"&\s*\n\s*", " ",
                 open(os.path.join(SRC, tag + "_Function.f90")).read())

    def sub(m):
        body = m.group(2)
        body = re.sub(r"\bV\((\d+)\)",
                      lambda k: "V<%s>" % names[int(k.group(1)) - 1], body)
        body = re.sub(r"\bF\((\d+)\)",
                      lambda k: "F<%s>" % names[nvar + int(k.group(1)) - 1], body)
        return body.strip()

    rates, vdot = {}, {}
    for m in re.finditer(r"^\s*A\((\d+)\)\s*=\s*(.*?)\s*$", txt, re.M):
        rates[int(m.group(1))] = sub(m)
    for m in re.finditer(r"^\s*Vdot\((\d+)\)\s*=\s*(.*?)\s*$", txt, re.M):
        vdot[names[int(m.group(1)) - 1]] = sub(m)
    return nvar, nspec, nreact, rates, vdot


def terms(e):
    """canonical term tuple, dummy reactions (A(j), j > NCHEM) removed.

    KPP emits the factors of a product in species-index order, and the species
    order changes between the two builds, so factors are sorted here: the
    comparison must be about the mathematics, not about KPP's emission order.
    """
    out = []
    for p in e.replace("-", "+-").split("+"):
        p = p.strip()
        if not p:
            continue
        m = re.search(r"A\((\d+)\)", p)
        if m and int(m.group(1)) > NCHEM:
            continue
        if "**" not in p:
            p = "*".join(sorted(f.strip() for f in p.split("*")))
        out.append(p)
    return tuple(sorted(out))


def raw_terms(e):
    return tuple(sorted(p.strip() for p in e.replace("-", "+-").split("+") if p.strip()))


def main():
    nv_o, ns_o, nr_o, ro, vo = build("old")
    nv_n, ns_n, nr_n, rn, vn = build("new")

    bad_r = [(j, ro.get(j), rn.get(j)) for j in range(1, NCHEM + 1)
             if terms(ro.get(j, "")) != terms(rn.get(j, ""))]
    bad_v = [(s, vo[s], vn.get(s, "")) for s in vo
             if terms(vo[s]) != terms(vn.get(s, ""))]

    print("old: NVAR=%d NSPEC=%d NREACT=%d" % (nv_o, ns_o, nr_o))
    print("new: NVAR=%d NSPEC=%d NREACT=%d" % (nv_n, ns_n, nr_n))
    print("chemistry rate definitions A(1..%d) mismatches : %d" % (NCHEM, len(bad_r)))
    for j, a, b in bad_r[:8]:
        print("   A(%d)\n     old: %s\n     new: %s" % (j, a, b))
    print("Vdot chemistry-part mismatches : %d (over %d species)"
          % (len(bad_v), len(vo)))
    for s, a, b in bad_v[:8]:
        print("   %s\n     old: %s\n     new: %s" % (s, a, b))
    gained = {s: len(raw_terms(vn[s])) - len(raw_terms(vo.get(s, ""))) for s in vn}
    gained = {s: k for s, k in gained.items() if k}
    print("species that gained dummy terms : %d %s"
          % (len(gained), sorted(gained)))
    ok = not bad_r and not bad_v
    print("VERDICT V1: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
