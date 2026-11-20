#!/usr/bin/env python3
"""Does freezing the winning window restore the antithetic correlation?

Control metric: the 8-h mean over ONE fixed time span (the span that wins at the
all-on state).  The functional form is unchanged; only the argmax is removed.

rho is computed PER PLAYER over the 40 antithetic chain pairs, then reported for
the dominant player and as the median across players.  Pooling players would mix
between-player variance into the correlation and is not informative.

Writes results_frozen_rho.txt
"""
import csv
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "results_frozen_rho.txt")
LABELS = ("dallas_2016w", "deerpark_2017s")


def load(pattern, col):
    d = {}
    for f in sorted(glob.glob(os.path.join(HERE, pattern))):
        for row in csv.DictReader(open(f)):
            try:
                d[int(row["name"][1:])] = float(row[col])
            except (ValueError, KeyError):
                continue
    return d


def per_player(label, tab):
    """{player: (deltas_forward, deltas_reverse)} over antithetic chain pairs."""
    ch = json.load(open(os.path.join(HERE, "probe_chains_%s.json" % label)))["chains"]
    acc = {i: ([], []) for i in range(16)}
    for p in range(0, len(ch) - 1, 2):
        a, b = ch[p], ch[p + 1]
        for k, i in enumerate(a["order"]):
            prev = a["masks"][k - 1] if k else 0
            if prev not in tab or a["masks"][k] not in tab:
                continue
            df = tab[a["masks"][k]] - tab[prev]
            kr = b["order"].index(i)
            prevr = b["masks"][kr - 1] if kr else 0
            if prevr not in tab or b["masks"][kr] not in tab:
                continue
            dr = tab[b["masks"][kr]] - tab[prevr]
            acc[i][0].append(df)
            acc[i][1].append(dr)
    return acc


def main():
    L, add = [], None
    add = L.append
    add("Frozen-window control: does removing the argmax restore rho_anti?")
    add("control = 8-h mean over the span that wins at the all-on state")
    add("rho per player over 40 antithetic chain pairs; factor = (1+rho)/2")
    add("")
    for label in LABELS:
        mad = load("probe_probe_%s_*.csv" % label, "v_mda8")
        fro = load("probe_frozen_%s_*.csv" % label, "v_frozen")
        pm, pf = per_player(label, mad), per_player(label, fro)
        add("=" * 80)
        add("%s   coalitions %d (MDA8) / %d (frozen)" % (label, len(mad), len(fro)))
        add("=" * 80)
        rows = []
        for tag, pp in (("MDA8", pm), ("frozen", pf)):
            rhos, facs, vs = [], [], []
            for i in range(16):
                df = np.array(pp[i][0], dtype=float)
                dr = np.array(pp[i][1], dtype=float)
                if len(df) < 5 or df.std() == 0 or dr.std() == 0:
                    continue
                r = float(np.corrcoef(df, dr)[0, 1])
                rhos.append(r)
                facs.append((1 + r) / 2)
                vs.append(float(np.var(df)))
            add("  %-7s rho: median %+.3f   min %+.3f   max %+.3f | variance factor median %.3f"
                % (tag, np.median(rhos), np.min(rhos), np.max(rhos), np.median(facs)))
            rows.append((np.median(rhos), vs))
        # per-player selector contribution
        share = []
        for i in range(16):
            dm = np.array(pm[i][0], dtype=float)
            dfz = np.array(pf[i][0], dtype=float)
            if dm.std() > 0:
                share.append(1 - float(np.var(dfz)) / float(np.var(dm)))
        add("  selector share of marginal variance: median %.1f%%   min %.1f%%   max %.1f%%"
            % (100 * np.median(share), 100 * np.min(share), 100 * np.max(share)))
        add("")

    txt = "\n".join(L)
    print(txt)
    with open(OUT, "w") as fh:
        fh.write(txt + "\n")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
