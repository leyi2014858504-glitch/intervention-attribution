#!/usr/bin/env python3
"""
ppbC -> CB05 lumped amounts.

The unit trap: AQS reports every VOC species in "Parts per billion Carbon"
(ppbC), while CB05 works in carbon bonds. A species measured at 5 ppbC is
5/nC ppb of molecules, and the mechanism's lumping then decides how many carbon
bonds that is:

  alkane Cn        -> n PAR bonds          => PAR += ppbC     (n cancels out)
  terminal alkene  -> 1 OLE + (n-2) PAR    => OLE += ppbC/n ; PAR += ppbC*(n-2)/n
  internal alkene  -> 1 IOLE + (n-2) PAR
  ethylene -> ETH ; isoprene -> ISOP ; terpene C10 -> APIN/BPIN/LIM
  benzene/ethylbenzene -> TOL ; xylene & trimethylbenzene -> XYL
  formaldehyde -> FORM ; acetaldehyde -> ALD2 ; C3+ aldehyde -> ALDX
  methanol -> MEOH ; ethanol -> ETOH ; ethane -> ETHA

So we need one thing AQS does not give: the carbon number nC.

--selftest checks the resolver against 78 hand-verified species. It caught two
real bugs, which is why it exists:
  * "methyl" contains the substring "ethyl" -> every methyl was double counted
  * chain matching took the SHORTEST key ("meth") before the longest ("hex")
"""
import re, sys

KNOWN = {
    "methane": 1, "ethane": 2, "ethylene": 2, "acetylene": 2, "propane": 3,
    "propylene": 3, "isoprene": 5, "1,3-butadiene": 4, "benzene": 6, "toluene": 7,
    "ethylbenzene": 8, "styrene": 8, "formaldehyde": 1, "acetaldehyde": 2,
    "acetone": 3, "methanol": 1, "ethanol": 2, "ethyl alcohol": 2,
    "methyl ethyl ketone": 4, "2-butanone": 4, "methyl vinyl ketone": 4,
    "methacrolein": 4, "crotonaldehyde": 4, "benzaldehyde": 7, "phenol": 6,
    "methylglyoxal": 3, "glyoxal": 2, "methyl tert-butyl ether": 5,
    "tert-amyl methyl ether": 6, "tert-butyl ethyl ether": 6,
    "formic acid": 1, "acetic acid": 2, "acetonitrile": 2, "acrylonitrile": 3,
    "2-nitropropane": 3, "propyne": 3, "isobutene": 4, "isobutylene": 4,
    "allene": 3, "propadiene": 3, "chloroprene": 4, "vinyl chloride": 2,
    "ethylene oxide": 2, "propylene oxide": 3, "cyclopentene": 5,
    "cyclohexene": 6, "methylcyclopentene": 6,
    # terpenes are C10 (sesquiterpenes C15): never guess these
    "alpha.-pinene": 10, "a-pinene": 10, "beta.-pinene": 10, "b-pinene": 10,
    "limonene": 10, "myrcene": 10, "ocimene": 10, "terpinene": 10,
    "carene": 10, "phellandrene": 10, "terpineol": 10, "terpene": 10,
    "humulene": 15, "caryophyllene": 15,
}

CHAIN = sorted([("meth", 1), ("eth", 2), ("prop", 3), ("but", 4), ("pent", 5),
                ("hex", 6), ("hept", 7), ("oct", 8), ("non", 9), ("dec", 10),
                ("undec", 11), ("dodec", 12), ("tridec", 13), ("tetradec", 14),
                ("pentadec", 15), ("hexadec", 16)],
               key=lambda t: -len(t[0]))            # longest key first: undec > dec > eth
AROMATIC = {"benzene": 6, "toluene": 7, "xylene": 8, "styrene": 8, "phenol": 6}
RING = {"cyclopent": 5, "cyclohex": 6, "cyclohept": 7, "cyclobut": 4}
# substituent -> carbons; matched longest-first so "isopropyl" wins over "propyl"
SUBST = [("isopropyl", 3), ("methyl", 1), ("ethyl", 2), ("propyl", 3),
         ("butyl", 4), ("pentyl", 5), ("amyl", 5)]
MULT = [("tetra", 4), ("tri", 3), ("di", 2), ("", 1)]
STEREO = {"cis", "trans", "n", "iso", "sec", "tert", "o", "m", "p",
          "alpha", "beta", "and", "total", "unverified", "verified", "s", "r"}


def _strip_substituents(token):
    """Pull leading multipliers+substituents off a token, longest-first."""
    extra, rest, changed = 0, token, True
    while changed:
        changed = False
        for mult, k in MULT:
            for sub, cs in SUBST:                 # SUBST is longest-first
                pat = mult + sub
                if pat and pat in rest:
                    extra += k * cs
                    rest = rest.replace(pat, "", 1)
                    changed = True
                    break
            if changed:
                break
    return extra, rest


def carbon_number(name):
    n = name.lower().strip()
    if n in KNOWN:
        return KNOWN[n]

    # co-eluting pairs "A & B": take the larger half (usually the dominant one)
    if "&" in n:
        vals = [carbon_number(p.strip(" .,")) for p in n.split("&")]
        vals = [v for v in vals if v]
        return max(vals) if vals else None

    # "m/p Xylene" and friends
    n = re.sub(r"^(m/p|m\s*\(and\s*p\)|o/m/p|m&p)\s*", "", n)
    n = re.sub(r"[(),;.]", " ", n)
    toks = [t for t in n.replace("-", " ").split() if t and not t.isdigit() and t not in STEREO]

    # Tokens are scanned in order and substituent carbons ACCUMULATE across them.
    # "2-Methyl-2-butene" splits into ["methyl", "butene"]: the +1 from methyl
    # lives in one token and the C4 parent in the next, so dropping the running
    # total between tokens loses the substituent entirely.
    extra_total = 0
    for tok in toks:
        extra, rest = _strip_substituents(tok)
        if not rest:                      # this token was pure substituent(s)
            extra_total += extra
            continue
        extra += extra_total
        for base, c in AROMATIC.items():
            if base in rest:
                return c + extra
        for base, c in RING.items():
            if base in rest:
                return c + extra
        if "acetate" in rest or "acrylate" in rest or "methacrylate" in rest:
            # esters must be caught BEFORE the chain test: "methacrylate"
            # contains "meth", which the chain matcher would happily read as C1.
            # methacrylic acid is C4, acrylic acid C3, acetic acid C2.
            if "methacrylate" in rest:
                acid = 4
            elif "acrylate" in rest:
                acid = 3
            else:
                acid = 2
            return extra + acid
        for key, c in CHAIN:
            if key in rest:
                return c + extra
        if "alcohol" in rest or rest.endswith("ol"):
            # reached only when there is no chain root in this token, e.g. the
            # bare word in "n-Propyl alcohol" where the alkyl carries the carbons
            return extra if extra else None

    # nothing matched: try the raw string one more time for chain roots
    for key, c in CHAIN:
        if key in n.replace(" ", ""):
            return c
    return None


def bonds(ppbC, nC, group):
    """Split a ppbC amount into CB05 groups, accounting for carbon bonds."""
    if not nC or nC <= 0:
        return {"PAR": ppbC}
    ppb = ppbC / nC
    if group == "PAR":
        return {"PAR": ppbC}                       # n carbons -> n PAR bonds
    if group in ("OLE", "IOLE"):
        if nC > 2:
            return {group: ppb, "PAR": ppb * (nC - 2)}
        return {group: ppb}
    return {group: ppb}


SELF_TEST = [
    ("n-Butane", 4), ("Isobutane", 4), ("Ethane", 2), ("Propane", 3),
    ("n-Pentane", 5), ("Isopentane", 5), ("n-Hexane", 6), ("n-Heptane", 7),
    ("n-Octane", 8), ("n-Nonane", 9), ("n-Decane", 10), ("n-Undecane", 11),
    ("n-Dodecane", 12), ("n-Tridecane", 13), ("Cyclohexane", 6),
    ("Methylcyclohexane", 7), ("Methylcyclopentane", 6),
    ("2,2,4-Trimethylpentane", 8), ("2,3,4-Trimethylpentane", 8),
    ("2,3-Dimethylpentane", 7), ("2,4-Dimethylpentane", 7),
    ("2,2-Dimethylpentane", 7), ("3-Methylheptane", 8), ("2-Methylhexane", 7),
    ("3-Ethylhexane", 8), ("Ethylcyclohexane", 8), ("2,2-Dimethylbutane", 6),
    ("2,3-Dimethylbutane", 6), ("2-Methylpentane", 6), ("3-Methylpentane", 6),
    ("2,5-Dimethylhexane", 8), ("2,4-Dimethylhexane", 8),
    ("cis-1,2-Dimethylcyclopentane", 7), ("trans-1,3-Dimethylcyclopentane", 7),
    ("Propylene", 3), ("1-Butene", 4), ("trans-2-Butene", 4), ("cis-2-Butene", 4),
    ("trans-2-Pentene", 5), ("cis-2-Pentene", 5), ("2-Methyl-2-butene", 5),
    ("1-Pentene", 5), ("3-Methyl-1-butene", 5), ("1-Hexene", 6),
    ("trans-2-Hexene", 6), ("Isoprene", 5), ("1,3-Butadiene", 4),
    ("Benzene", 6), ("Toluene", 7), ("Ethylbenzene", 8), ("o-Xylene", 8),
    ("m/p Xylene", 8), ("Styrene", 8), ("1,3,5-Trimethylbenzene", 9),
    ("1,2,4-Trimethylbenzene", 9), ("1,2,3-Trimethylbenzene", 9),
    ("n-Propylbenzene", 9), ("Isopropylbenzene", 9), ("o-Ethyltoluene", 9),
    ("m-Diethylbenzene", 10), ("Formaldehyde", 1), ("Acetaldehyde", 2),
    ("Propionaldehyde", 3), ("Butyraldehyde", 4), ("Hexanaldehyde", 6),
    ("Benzaldehyde", 7), ("Acetone", 3), ("Methyl ethyl ketone", 4),
    ("3-Pentanone", 5), ("2-Pentanone", 5), ("Methanol", 1), ("Ethanol", 2),
    ("1-Hexanol", 6), ("n-Propyl alcohol", 3), ("n-Butyl alcohol", 4),
    ("alpha.-Pinene", 10), ("beta.-Pinene", 10), ("Limonene", 10),
    ("Methyl tert-butyl ether", 5), ("Ethyl acetate", 4),
    ("n-Propyl acetate", 5), ("Methyl methacrylate", 5), ("Ethyl acrylate", 5),
    ("Methyl Vinyl Ketone", 4), ("Acetylene", 2), ("Ethyne & isobutane", 4),
]

if __name__ == "__main__":
    if "--selftest" in sys.argv:
        bad = []
        for nm, want in SELF_TEST:
            got = carbon_number(nm)
            if got != want:
                bad.append((nm, want, got))
        print(f"self-test: {len(SELF_TEST)-len(bad)}/{len(SELF_TEST)} correct")
        for b in bad:
            print("  MISMATCH name=%-32s expected=%s got=%s" % b)
        sys.exit(1 if bad else 0)
    print(__doc__)
