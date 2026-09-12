"""Scale stress test, part 1: decide WHICH coalitions the sampler will need.

Design (no emulator, real chemistry):
  * n=24 = 3 sectors x 8 three-hour windows.  2**24 = 16.7M coalitions: exact
    enumeration is out, so the sampler visits only what it needs and a cache
    makes repeats free.  The total chemistry cost is the size of the UNION of
    visited coalitions, which is far smaller than the product of seeds x budget.
  * At n=12 we already own the exact table, so the sampler can be calibrated
    there: for the SAME estimator and budget we can see both the seed-to-seed
    spread AND the true error.  That gives the spread->error map used to read
    the n=24 numbers, where no oracle exists.

Sampler = permutation + antithetic (the strongest arm in the n=8 benchmark).
A "pair" evaluates all prefixes of a random permutation and of its reverse.
Budget B = target number of DISTINCT coalition evaluations.

Output: nmldir_scale/<scen>_m<08d>.nml  (union over all seeds/budgets, deduped)
"""
import itertools
import json
import os
import random

from kpp_box import nml_text

HERE = os.path.dirname(os.path.abspath(__file__))
N = 24                      # 3 sectors x 8 three-hour windows
NWIN = 8
SCEN = {"r_sat": {}, "r_lim1": {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5}}
BUDGETS = [64, 128, 256]
NSEED = 8


def cfg_of(mask):
    def win(base):
        return [0.0 if (mask >> (base + w)) & 1 else 1.0 for w in range(NWIN)]
    return {"NOXW": win(0), "VOCW": win(8), "COW": win(16)}


def perm_coalitions(rng, budget, n=N):
    """Coalitions visited by permutation+antithetic until `budget` distinct."""
    seen = set()
    while len(seen) < budget:
        order = list(range(n))
        rng.shuffle(order)
        for seq in (order, order[::-1]):
            mask = 0
            seen.add(0)
            for i in seq:
                mask |= 1 << i
                seen.add(mask)
            seen.add((1 << n) - 1)
            if len(seen) >= budget:
                break
    return seen


def main():
    out = os.path.join(HERE, "nmldir_scale")
    os.makedirs(out, exist_ok=True)
    for f in os.listdir(out):
        if f.endswith(".nml"):
            os.remove(os.path.join(out, f))

    manifest = {}
    total = 0
    for scen in SCEN:
        need = {0, (1 << N) - 1}          # always include the two anchors
        for B in BUDGETS:
            for s in range(NSEED):
                rng = random.Random(1000 * B + s)
                need |= perm_coalitions(rng, B)
        manifest[scen] = {"n_seeds": NSEED, "budgets": BUDGETS,
                          "coalitions": len(need)}
        for mask in sorted(need):
            cfg = dict(SCEN[scen])
            cfg.update(cfg_of(mask))
            with open(os.path.join(out, "%s_m%08d.nml" % (scen, mask)), "w",
                      newline="\n") as f:
                f.write(nml_text(cfg, "co"))
        total += len(need)
        print("%-8s coalitions needed = %d" % (scen, len(need)))
    print("TOTAL chemistry runs = %d" % total)
    print("(exact enumeration would be 2 * %d = %d)"
          % (1 << N, 2 * (1 << N)))
    with open(os.path.join(HERE, "scale_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=1)


if __name__ == "__main__":
    main()
