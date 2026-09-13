#!/usr/bin/env python3
"""
Observation species -> CB05 lumped group.

Only species that actually exist at the chosen sites are emitted, taken from
aqs_monitors.csv, so the table stays short and checkable.

Confidence column:
  H  standard CB05 lumping, no judgement needed
  M  conventional but approximate (a lumped surrogate is being used)
  L  genuinely uncertain -- pin to a citation before it goes in the paper
  X  has no CB05 species at all; trace at ambient levels, dropped
  -  not a reactive VOC: diagnostic only

Known limitations of this table (worth stating in the paper):
  * CO-ELUTING PAIRS. AQS names such as "1-Hexene & 2-Methyl-1-pentene" are two
    compounds that the chromatographic method cannot separate. The pair gets the
    group of whichever half the regex hits first; both halves matter here.
  * KETONES. CB05 has no ketone group, so acetone/MEK/3-pentanone are pushed onto
    PAR, which OVERSTATES their reactivity (acetone's lifetime is ~weeks).
  * BENZENE. No CB05 benzene; folding it to TOL overstates it, since benzene is
    markedly less reactive than toluene.
  * ACETYLENE. No CB05 species; normally dropped, which is defensible.
"""
import csv, os, re
from collections import defaultdict

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
OUT = "/mnt/d/atmos/results_voc_map.txt"

SITES = {
    "LA_06-037-1103":        ("06", "037", "1103"),
    "Charlotte_37-119-0041": ("37", "119", "0041"),
}


def is_voc(pc):
    """43xxx hydrocarbons/oxygenates, plus 45xxx -- the AROMATICS.

    452xx/451xx must be included: benzene 45201, toluene 45202, ethylbenzene
    45203, o-xylene 45204, m/p-xylene 45109, styrene 45220, the tri- and
    tetra-methylbenzenes 45207/45208/45225, ethyltoluenes 45211-45213.
    """
    return (pc.startswith("43") and pc != "43000") or pc.startswith("45")

# --- rule table: (regex on the parameter name, CB05 group, confidence, note).
# First match wins, so the order below IS the logic: exact names, then things
# that are not reactive VOCs at all, then specific functional groups, then the
# generic "-ene"/alkane catch-alls last.
RULES = [
    # ---- 1. exact small molecules
    (r"^methane$",                          "CH4",  "H", ""),
    (r"^ethane$",                           "ETHA", "H", "CB05 keeps ethane separate"),
    (r"^ethylene$",                         "ETH",  "H", ""),
    (r"^acetylene$",                        "-",    "L", "alkyne: no CB05 species, low reactivity"),
    (r"^propane$",                          "PAR",  "H", "reactivity-weighted into PAR"),
    (r"^propylene$",                        "OLE",  "H", ""),
    (r"^isoprene$",                         "ISOP", "H", ""),
    (r"^1,3-butadiene$",                    "OLE",  "M", "no diene in CB05; folded to OLE"),
    (r"^benzene$",                          "TOL",  "L", "no benzene in CB05; TOL overstates it"),
    (r"^toluene$",                          "TOL",  "H", ""),
    (r"^ethylbenzene$",                     "TOL",  "M", "CB05 lumps C2-benzenes with TOL"),
    (r"^styrene$",                          "OLE",  "L", "no CB05 species; vinyl group -> OLE"),
    (r"^formaldehyde$",                     "FORM", "H", ""),
    (r"^acetaldehyde$",                     "ALD2", "H", ""),
    (r"^acetone$",                          "PAR",  "L", "CB05 has no ketone group; PAR overstates it"),
    (r"^methyl ethyl ketone$|^2-butanone$", "PAR",  "L", "CB05 has no ketone group"),
    (r"^methanol$",                         "MEOH", "H", ""),
    (r"^ethanol$|^ethyl alcohol$",          "ETOH", "H", ""),
    (r"^formic acid$",                      "FACD", "H", ""),
    (r"^acetic acid$",                      "AACD", "H", ""),
    (r"^methylglyoxal$",                    "MGLY", "H", ""),
    (r"^glyoxal$",                          "FORM", "L", "no glyoxal in CB05; FORM is the only option"),

    # ---- 2. not VOC, or no CB05 species at all
    (r"^radon$",                            "-",    "-", "not a VOC"),
    (r"mercaptan|reduced sulfur|hydrogen sulfide|carbonyl sulfide",
                                            "-",    "L", "no CB05 species (CB05 sulfur = SO2/SULF only)"),
    (r"nitrile|nitropropane|nitro",         "-",    "L", "no CB05 species; trace"),
    (r"total chlorinated|total hydrocarbon|total nmoc|^volatile organic compounds$|"
     r"non-methane organic|unidentified|^sum of|isomers of|duplicate use|^aldehydes$",
                                            "-",    "-", "aggregate/diagnostic, not a single species"),
    (r"heptachlor|endrin|endosulfan|aldrin|dieldrin|methoxychlor|guthion|atrazine|"
     r"dacthal|dursban|parathion|malathion|lindane|hexachloro|chlordane|toxaphene|"
     r"dichlorophenoxy|bis *\(|pyrene|fluoranthene|anthracene|phenanthrene|chrysene|"
     r"benzo|naphthalene|biphenyl|pcb|phosphorotrithioic",
                                            "-",    "-", "PAH/pesticide/PCB: not in the O3 mechanism"),

    # ---- 3. aromatics (before the bare "-ene" catch-all, else these fall through).
    #         Placed ahead of terpenes too, because the co-eluting pairs here are
    #         named "trimethylbenzene & beta-pinene" and the aromatic half dominates.
    (r"xylene",                             "XYL",  "H", ""),
    (r"trimethylbenzene|ethyltoluene|ethytoluene|trimethyltoluene|tetramethylbenzene|"
     r"butylbenzene|cumene|isopropylbenzene|propylbenzene|diethylbenzene",
                                            "XYL",  "M", "CB05 XYL covers di-/tri-methylbenzenes"),
    (r"cresol|phenol|benzyl",               "CRES", "M", "CB05 CRES covers cresols/higher aromatics"),

    # ---- 4. terpenes / sesquiterpenes (what the SORG section exists for)
    (r"alpha\.?-?pinene|^a-pinene",         "APIN", "H", "SORG species"),
    (r"beta\.?-?pinene|^b-pinene",          "BPIN", "H", "SORG species"),
    (r"limonene",                           "LIM",  "H", "SORG species"),
    (r"humulene|caryophyllene",             "HUM",  "H", "SORG species"),
    (r"myrcene|ocimene|terpinene|carene|phellandrene|terpineol|terpene",
                                            "TERP", "M", "SORG lumped terpene"),

    # ---- 5. C3+ aldehydes
    (r"propionaldehyde|butyraldehyde|valeraldehyde|isovaleraldehyde|isobutyraldehyde|"
     r"hexanal|heptanal|octanal|nonanal|decanal|undecanal|dodecanal|tridecanal|"
     r"tetradecanal|benzaldehyde|tolualdehyde|crotonaldehyde|methacrolein|acrolein|"
     r"aldehyde|dialdehyde",
                                            "ALDX", "M", "C3+ aldehyde lumped surrogate"),

    # ---- 6. alkynes
    (r"yne$|propyne|ethyne|butyne",         "-",    "L", "alkyne: no CB05 species"),

    # ---- 7. halogenated: trace, and no CB05 species for any of them
    (r"chlor|brom|iod|fluor|freon|halocarbon|halogen|carbon tetrachloride",
                                            "-",    "X", "halogenated: no CB05 species; trace at ambient levels"),

    # ---- 7b. remaining aromatics, caught last so that the pure halocarbons above
    #          (chlorobenzene, trichlorobenzene) do NOT get pulled in here
    (r"benzene",                            "XYL",  "M", "aromatic; lumped to XYL"),

    # ---- 8. C3+ alcohols / esters / ethers -- absent from CB05, surrogate is a choice
    (r"propanol|butanol|pentanol|hexanol|heptanol|octanol|amyl alcohol|"
     r"propyl alcohol|butyl alcohol|glycol|glycerol",
                                            "PAR",  "L", "alcohol not in CB05; PAR surrogate is a guess"),
    (r"acetate|acrylate|methacrylate|ether|mtbe|tame|etbe|tert-butyl|tert-amyl",
                                            "PAR",  "L", "ester/ether not in CB05; PAR surrogate is a guess"),

    # ---- 9. ketones (CB05 has no ketone group at all)
    (r"ketone|one$",                        "PAR",  "L", "CB05 has no ketone group"),

    # ---- 10. epoxides
    (r"oxide$|epoxide",                     "PAR",  "L", "no CB05 species; PAR surrogate is a guess"),

    # ---- 11. internal olefins (double bond between two substituted carbons)
    (r"(cis|trans)-[2-9]-|(^|[^0-9])[2-9]-[a-z]*ene|"
     r"methyl-2-butene|^2-pentene|cyclohexene|cyclopentene",
                                            "IOLE", "H", "internal double bond"),

    # ---- 12. terminal olefins (including the bare "-ene" names)
    (r"(^|[^0-9])1-[a-z]*ene|diene|ene$|ylene|vinyl",
                                            "OLE",  "H", "terminal double bond"),

    # ---- 13. everything else aliphatic = paraffin carbon bonds
    (r"butane|pentane|hexane|heptane|octane|nonane|decane|undecane|dodecane|tridecane|"
     r"tetradecane|pentadecane|hexadecane|cyclohexane|cyclopentane|methylcyclo|"
     r"dimethyl|trimethyl|methylpent|methylhex|methylhept|methylbutane|methylpropane|"
     r"paraffin|alkane",
                                            "PAR",  "H", "alkane/cycloalkane carbon bonds"),
]

COMP = [(re.compile(p, re.I), g, c, n) for p, g, c, n in RULES]


def classify(name):
    for rx, grp, conf, note in COMP:      # first match wins, order matters
        if rx.search(name):
            return grp, conf, note
    return "?", "?", "unclassified -- needs a manual decision"


def main():
    want = {k: {} for k in SITES}
    with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            key = (row["State Code"], row["County Code"], row["Site Number"])
            for lab, tgt in SITES.items():
                if key == tgt:
                    pc = row["Parameter Code"]
                    if is_voc(pc):
                        want[lab][pc] = row.get("Parameter Name", "").strip()

    log = open(OUT, "w", encoding="utf-8")
    for lab in SITES:
        log.write(f"===== {lab} : {len(want[lab])} VOC parameters in the inventory =====\n")
        log.write(f"{'code':<7}{'cb05':<6}{'conf':<5}name  |  note\n")
        log.write("-" * 104 + "\n")
        counts = defaultdict(int)
        for pc in sorted(want[lab]):
            nm = want[lab][pc]
            grp, conf, note = classify(nm)
            counts[grp] += 1
            log.write(f"{pc:<7}{grp:<6}{conf:<5}{nm}  |  {note}\n")
        real = sum(n for g, n in counts.items() if g not in ("-", "?"))
        log.write(f"\n  groups: " + ", ".join(f"{g}={n}" for g, n in sorted(counts.items())) + "\n")
        log.write(f"  reactive species actually mapped to a CB05 group: {real}\n\n")
    log.close()
    print("wrote", OUT)


if __name__ == "__main__":
    main()
