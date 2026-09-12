"""Generate the 12-player (3 sectors x 4 windows) namelist sweep.

Player layout (mask bit -> switch):
    bit 0..3   NOx window W1..W4 OFF
    bit 4..7   VOC window W1..W4 OFF
    bit 8..11  CO  window W1..W4 OFF
so mask 0 = everything on (baseline) and 4095 = everything off.
2**12 = 4096 coalitions per scenario -> the Moebius spectrum is still exactly
enumerable, no surrogate and no sampling.

Only two scenarios for now (the pilot for the third sector):
    r_sat   NOx-saturated side
    r_lim1  NOx-limited side  (where the NOx lever actually works)

Output: nmldir_co/<scen>_m<04d>.nml  (+ co_list.txt)
"""
import os

from kpp_box import nml_text

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 12
MAX1 = 1 << NPLAY
WIDTH = 4

SCEN = {
    "r_sat": {},
    "r_lim1": {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5},
}


def cfg_of(mask):
    def win(base):
        return [0.0 if (mask >> (base + w)) & 1 else 1.0 for w in range(4)]
    return {"NOXW": win(0), "VOCW": win(4), "COW": win(8)}


def main():
    out = os.path.join(HERE, "nmldir_co")
    os.makedirs(out, exist_ok=True)
    for f in os.listdir(out):
        if f.endswith(".nml"):
            os.remove(os.path.join(out, f))
    n = 0
    for name, cfg0 in SCEN.items():
        for mask in range(MAX1):
            cfg = dict(cfg0)
            cfg.update(cfg_of(mask))
            with open(os.path.join(out, "%s_m%0*d.nml" % (name, WIDTH, mask)),
                      "w", newline="\n") as f:
                f.write(nml_text(cfg, "co"))
            n += 1
    with open(os.path.join(HERE, "co_list.txt"), "w") as f:
        for name in SCEN:
            f.write("%s CO12\n" % name)
    print("wrote %d namelists (%d scenarios x %d)"
          % (n, len(SCEN), MAX1))


if __name__ == "__main__":
    main()
