"""Locate the EKMA ridge under candidate scenario configs.

For each config (emission strength x dilution), run the nox-column
{0.4,0.7,1.0,1.5,2.0} at voc=1.5 and print MDA8, to see which settings
place the baseline (nox=1.0) on the NOX-SATURATED side (peak at or below
nox ~0.7-1.0) with realistic absolute O3 levels. dt=600 for speed.
"""
from kpp_box import run_scenario

CONFS = {
    # name: overrides
    "A_emis8_kd3h":  dict(EMIS_NOX=8.0,  EMIS_VOC=6.0, BG_VOC=20.0, KDIL_DAY=10800.0, KDIL_NIGHT=86400.0),
    "B_emis15_kd3h": dict(EMIS_NOX=15.0, EMIS_VOC=5.0, BG_VOC=15.0, KDIL_DAY=10800.0, KDIL_NIGHT=86400.0),
    "C_emis15_kd6h": dict(EMIS_NOX=15.0, EMIS_VOC=5.0, BG_VOC=15.0, KDIL_DAY=21600.0, KDIL_NIGHT=129600.0),
    "D_emis25_kd6h": dict(EMIS_NOX=25.0, EMIS_VOC=5.0, BG_VOC=12.0, KDIL_DAY=21600.0, KDIL_NIGHT=129600.0),
}
NOXS = [0.4, 0.7, 1.0, 1.5, 2.0]

for name, cfg in CONFS.items():
    row = []
    for nf in NOXS:
        c = dict(cfg)
        c["NOX_FACTOR"] = nf
        c["VOC_FACTOR"] = 1.5
        c["DT_S"] = 600.0
        m = run_scenario(c, "probe_%s_n%.1f" % (name, nf))
        row.append(m["mda8"])
    print("%-14s voc=1.5 MDA8 over nox%s: %s" %
          (name, NOXS, " ".join("%.1f" % v for v in row)), flush=True)
