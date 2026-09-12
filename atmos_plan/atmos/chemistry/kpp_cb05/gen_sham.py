"""SHAM-SWITCH CONTROL.

Set EMIS_CO = 0 so that the four CO windows become NULL PLAYERS while the
eight NOx/VOC windows stay live.  Because E_CO = EMIS_CO*COW(IW)*... is then
identically 0.0 (not merely small), the value function must be EXACTLY
independent of bits 8..11, which makes every one of the checks below exact
rather than approximate:

  A  bit <-> window mapping   v12(mask) == v8(mask & 0xFF) for all 4096 masks
  B  dummy axiom              phi_CO  == 0
  C  null players change nothing  phi_NOx/VOC == the frozen 8-player values
  D  efficiency               sum(phi) == v(full) - v(none)
  E  spectrum                 m(A) == 0 for every A touching a CO bit

Checks A and C are the point: a null player contributes zero marginal value in
every permutation, so the Shapley values of the live players must be bit-for-
bit unchanged.  Any deviation is a wiring bug, not noise.

Output: /root/sham_nml/<scen>_m<04d>.nml   (2 scenarios x 4096 = 8192)

NOTE: written into the WSL ext4 filesystem, not the Windows mount -- creating
8k small files over 9p is slow enough that the sandbox kills the writer
half way through.
"""
import os

from kpp_box import nml_text

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = "/root/sham_nml"
NPLAY, WIDTH = 12, 4
MAX1 = 1 << NPLAY

SCEN = {
    "sham_sat":  {"EMIS_CO": 0.0},
    "sham_lim1": {"EMIS_CO": 0.0, "NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5},
}


def cfg_of(mask):
    def win(base):
        return [0.0 if (mask >> (base + w)) & 1 else 1.0 for w in range(4)]
    return {"NOXW": win(0), "VOCW": win(4), "COW": win(8)}


def main():
    out = OUT
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
    print("wrote %d namelists (%d scenarios x %d)"
          % (n, len(SCEN), MAX1))


if __name__ == "__main__":
    main()
