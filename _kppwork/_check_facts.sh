#!/bin/bash
# Owen-chain reachability for the other configurations (pure combinatorics, no chemistry):
#   n=12 CO extension: 4 groups x 3 players  -> ?
#   n=16 3-hour game : 8 groups x 2 players  -> ?
cd /mnt/d/atmos
python3 - <<'EOF'
import random

def reach_count(nsec, nwin, trials=200000, seed=0):
    groups = [[w * nsec + s for s in range(nsec)] for w in range(nwin)]
    n = nsec * nwin
    rng = random.Random(seed)
    reach = {0}
    for _ in range(trials):
        gid = list(range(nwin)); rng.shuffle(gid)
        S = 0
        for g in gid:
            inside = list(groups[g]); rng.shuffle(inside)
            for i in inside:
                S |= 1 << i
                reach.add(S)
    return len(reach), 1 << n

for nsec, nwin, tag in ((2, 4, "n=8  (4 groups x 2)"),
                        (3, 4, "n=12 (4 groups x 3, CO ext)"),
                        (2, 8, "n=16 (8 groups x 2)"),
                        (3, 8, "n=24 (8 groups x 3, stress only)")):
    r, tot = reach_count(nsec, nwin)
    print("%s: Owen chains reach %5d of %8d  (%.1f%%)" % (tag, r, tot, 100.0 * r / tot))
EOF