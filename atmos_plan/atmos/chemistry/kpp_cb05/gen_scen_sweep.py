"""Generate the scenario x coalition namelist sweep for the adaptivity test.

TWO SETS OF SCENARIOS (all reuse the /root/kpp5 dummy-forcing driver; only
RUNCFG knobs are touched, so every scenario is a legitimate weather/chem state
of the same box -- no code changes):

  U-set (10) -- uncertainty inside ONE planning problem: perturbations of the
      baseline saturated-city episode (temperature, background O3/NO2, mixing,
      CO, biogenics). The adaptivity gap over this set is the ARBITER for
      route A (does stochastic planning have value at all?).

  R-set (6)  -- EKMA regime walk (nox/voc factor combos incl. a weekend
      profile). Used for route B material: how the Möbiius spectrum, Shapley
      signs and the optimal recipe move between NOx-saturated and NOx-limited
      regimes. NOT pooled into the U-set gap.

Outputs:
  nmldir_scen/scN_mMM.nml   one namelist per (scenario, coalition)
  scen_list.txt             scenario names + which set they belong to
"""
import os

from kpp_box import nml_text

HERE = os.path.dirname(os.path.abspath(__file__))
NPLAY = 8
MAX1 = 1 << NPLAY

U_SET = {
    "u0_base":   {},
    "u1_hot":    {"TMAX_K": 304.0, "TMIN_K": 291.0},
    "u2_cool":   {"TMAX_K": 298.0, "TMIN_K": 285.0},
    "u3_bgO3hi": {"BG_O3": 40.0},
    "u4_bgNOxhi": {"BG_O3": 20.0, "BG_NO2": 8.0},
    "u5_mixdeep": {"KDIL_DAY": 900.0, "KDIL_NIGHT": 43200.0},
    "u6_mixshal": {"KDIL_DAY": 2700.0},
    "u7_cohigh": {"EMIS_CO": 80.0},
    "u8_bgvochi": {"BG_VOC": 30.0},
    "u9_isophi": {"EMIS_ISO": 2.0},
}

R_SET = {
    "r_sat":     {},
    "r_sat_hi":  {"NOX_FACTOR": 1.3, "VOC_FACTOR": 1.3},
    "r_trans":   {"VOC_FACTOR": 0.6},
    "r_lim1":    {"NOX_FACTOR": 0.5, "VOC_FACTOR": 1.5},
    "r_lim2":    {"NOX_FACTOR": 0.3},
    "r_wkend":   {"WKEND": True, "NOX_FACTOR": 0.7},
}


def cfg_of(mask):
    return {"NOXW": [0.0 if (mask >> w) & 1 else 1.0 for w in range(4)],
            "VOCW": [0.0 if (mask >> (4 + w)) & 1 else 1.0 for w in range(4)]}


def main():
    out = os.path.join(HERE, "nmldir_scen")
    os.makedirs(out, exist_ok=True)
    for f in os.listdir(out):
        if f.endswith(".nml"):
            os.remove(os.path.join(out, f))
    n = 0
    for setname, scen in list(U_SET.items()) + list(R_SET.items()):
        for mask in range(MAX1):
            cfg = dict(scen)
            cfg.update(cfg_of(mask))
            with open(os.path.join(out, "%s_m%03d.nml" % (setname, mask)), "w",
                      newline="\n") as f:
                f.write(nml_text(cfg, "dummy"))
            n += 1
    with open(os.path.join(HERE, "scen_list.txt"), "w") as f:
        for s in list(U_SET) + list(R_SET):
            f.write("%s %s\n" % (s, "U" if s in U_SET else "R"))
    print("wrote %d namelists (%d scenarios x 256)" % (n, len(U_SET) + len(R_SET)))


if __name__ == "__main__":
    main()
