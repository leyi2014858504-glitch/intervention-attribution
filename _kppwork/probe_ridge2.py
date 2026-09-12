"""Second ridge probe: absolute pollution level at fixed NOx:VOC=8:6.

Dilution x emissions (scaled together) so that the voc=1.5 column along
NOx peaks at or below factor 1.0 (saturated side) while weekday baseline
(nox=1.0, voc=1.0) MDA8 stays realistic (60-140 ppb).
"""
from kpp_box import run_scenario

CONFS = {
    "E_e10_kd25": dict(EMIS_NOX=10.0, EMIS_VOC=7.5, BG_VOC=20.0, KDIL_DAY=9000.0,  KDIL_NIGHT=86400.0),
    "F_e12_kd2":  dict(EMIS_NOX=12.0, EMIS_VOC=9.0,  BG_VOC=20.0, KDIL_DAY=7200.0,  KDIL_NIGHT=86400.0),
    "G_e16_kd15": dict(EMIS_NOX=16.0, EMIS_VOC=12.0, BG_VOC=20.0, KDIL_DAY=1500.0,  KDIL_NIGHT=86400.0),
    "H_e12_kd3":  dict(EMIS_NOX=12.0, EMIS_VOC=9.0,  BG_VOC=20.0, KDIL_DAY=10800.0, KDIL_NIGHT=43200.0),
}
NOXS = [0.7, 1.0, 1.5]

for name, cfg in CONFS.items():
    row = []
    for nf in NOXS:
        c = dict(cfg)
        c["NOX_FACTOR"] = nf
        c["VOC_FACTOR"] = 1.5
        c["DT_S"] = 600.0
        m = run_scenario(c, "probe2_%s_n%.1f" % (name, nf))
        row.append(m["mda8"])
    b = run_scenario({**cfg, "VOC_FACTOR": 1.0, "DT_S": 600.0}, "probe2_%s_base" % name)
    print("%-11s voc1.5 nox[0.7,1.0,1.5]= %s | baseline mda8=%.1f d3max=%.1f" %
          (name, " ".join("%.1f" % v for v in row), b["mda8"], b["d3_max"]), flush=True)
