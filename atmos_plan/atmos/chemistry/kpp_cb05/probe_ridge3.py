"""Third probe: find a config with peak(nox) <= 0.7 at voc=1.5 AND a
realistic baseline (60-110 ppb MDA8)."""
from kpp_box import run_scenario

CONFS = {
    "G_full":  dict(EMIS_NOX=16.0, EMIS_VOC=12.0, BG_VOC=20.0, KDIL_DAY=1500.0,  KDIL_NIGHT=86400.0),
    "N_e14_15h": dict(EMIS_NOX=14.0, EMIS_VOC=10.5, BG_VOC=20.0, KDIL_DAY=5400.0,  KDIL_NIGHT=86400.0),
    "O_e20_2h":  dict(EMIS_NOX=20.0, EMIS_VOC=15.0, BG_VOC=20.0, KDIL_DAY=7200.0,  KDIL_NIGHT=43200.0),
}
for name, cfg in CONFS.items():
    row = []
    for nf in [0.4, 0.7, 1.0]:
        c = dict(cfg)
        c["NOX_FACTOR"] = nf
        c["VOC_FACTOR"] = 1.5
        c["DT_S"] = 600.0
        m = run_scenario(c, "probe3_%s_n%.1f" % (name, nf))
        row.append(m["mda8"])
    b = run_scenario({**cfg, "VOC_FACTOR": 1.0, "DT_S": 600.0}, "probe3_%s_base" % name)
    w = run_scenario({**cfg, "VOC_FACTOR": 1.0, "WKEND": True, "NOX_FACTOR": 0.7,
                      "DT_S": 600.0}, "probe3_%s_wknd" % name)
    print("%-10s voc1.5 nox[0.4,0.7,1.0]=%s | base mda8=%.1f d3max=%.1f | wknd d3max=%.1f (want > %.1f)" %
          (name, " ".join("%.1f" % v for v in row), b["mda8"], b["d3_max"], w["d3_max"], b["d3_max"]),
          flush=True)
