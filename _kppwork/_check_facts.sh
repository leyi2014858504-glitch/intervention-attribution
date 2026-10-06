#!/bin/bash
# Verify: (a) reachable coalition set of group-contiguous Owen chains at n=8;
#         (b) budget overshoot of exp_scale.gen (chain completed past B).
cd /mnt/d/atmos
python3 - <<'EOF'
import random, itertools

# (a) Owen chain reachability, n=8, groups = [[w, 4+w] for w in range(4)]
groups = [[w, 4 + w] for w in range(4)]
reach = set()
rng = random.Random(0)
for _ in range(200000):
    gid = list(range(4)); rng.shuffle(gid)
    S = 0; reach.add(0)
    for g in gid:
        inside = list(groups[g]); rng.shuffle(inside)
        for i in inside:
            S |= 1 << i
            reach.add(S)
print("psi reachable coalitions (n=8):", len(reach), "of 256")
unions = {m for m in reach if all((m >> j) & 1 == (m >> (4 + j)) & 1
                                   for j in range(4)) or True}
full_unions = set()
for t in itertools.product([0, 1], repeat=4):
    m = 0
    for w in range(4):
        if t[w]:
            m |= 1 << w | 1 << (4 + w)
    full_unions.add(m)
print("  of which whole-group unions:", len(full_unions & reach))

# (b) budget overshoot of exp_scale.gen
def gen(rng, budget, n):
    seen = set(); orders = []
    while len(seen) < budget:
        order = list(range(n)); rng.shuffle(order)
        for seq in (order, order[::-1]):
            mask = 0; seen.add(0)
            for i in seq:
                mask |= 1 << i; seen.add(mask)
            seen.add((1 << n) - 1); orders.append(seq)
            if len(seen) >= budget:
                break
    return orders, seen

for n in (8, 12, 16):
    for B in (64, 128):
        over = []
        for sd in range(50):
            _, seen = gen(random.Random(1000 * B + sd), B, n)
            over.append(len(seen) - 1 - B)   # minus the free empty set
        print("gen: n=%2d B=%3d  realised runs: min %d  max %d  (nominal %d, overshoot up to %d)"
              % (n, B, B + min(over), B + max(over), B, max(over)))
EOF