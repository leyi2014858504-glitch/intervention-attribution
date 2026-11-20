"""Python driver layer for the KPP 3.5 CB05 scenario box (/root/kpp4).

run_scenario(cfg, name) writes a Fortran namelist, copies it into the WSL
build dir, executes cb05.exe, pulls traj.dat + state_out.dat back to
d:/atmos/_kppwork/g4runs/, and returns parsed metrics.

All values map 1:1 to NAMELIST /RUNCFG/ in cb05s_Main.f90.
"""
import os
import subprocess
import math

HERE = os.path.dirname(os.path.abspath(__file__))

# Two scenario backends, same namelist interface:
#   split : /root/kpp4, cb05s_Main.f90 - emissions/dilution as an operator-
#           splitting step after each INTEGRATE (G4 reference build)
#   dummy : /root/kpp5, cb05d_Main.f90 - the same forcing expressed as dummy
#           reactions inside the mechanism (Möbius build)
BACKENDS = {"split": "/root/kpp4", "dummy": "/root/kpp5", "co": "/root/kpp6",
            "anchor": "/root/kppa", "scale": "/root/kpp7",
            "tight": "/root/kppt"}
DEFAULT_BACKEND = os.environ.get("KPP_BACKEND", "split")

# RUNCFG keys understood by every backend ...
BASE_KEYS = ["TEND_S", "DT_S", "STATE_IN", "EMIS_NOX", "EMIS_VOC", "EMIS_ISO",
             "EMIS_CO", "NOX_FACTOR", "VOC_FACTOR", "NOXW", "VOCW", "WKEND",
             "KDIL_DAY", "KDIL_NIGHT", "BG_O3", "BG_NO2", "BG_NO", "BG_CO",
             "BG_VOC", "BG_ISO", "TMIN_K", "TMAX_K"]
# ... and the ones only the dummy build knows about.  "co" is the same dummy
# build plus the third intervention sector (COW), i.e. the 12-player game;
# "dummy" deliberately omits COW so its namelists stay 8-player.
EXTRA_KEYS = {"dummy": ["DUMMYF", "CONSTJ", "C_H2O"],
              "co": ["DUMMYF", "CONSTJ", "C_H2O", "COW"],
              "anchor": ["DUMMYF", "CONSTJ", "C_H2O"],
              "scale": ["DUMMYF", "CONSTJ", "C_H2O", "COW"]}

DEFAULTS = dict(
    TEND_S=259200.0, DT_S=300.0, STATE_IN="none",
    EMIS_NOX=16.0, EMIS_VOC=12.0, EMIS_ISO=1.0, EMIS_CO=40.0,
    NOX_FACTOR=1.0, VOC_FACTOR=1.0, WKEND=False,
    # per-window emission multipliers: [00-06) [06-12) [12-18) [18-24)
    NOXW=[1.0, 1.0, 1.0, 1.0], VOCW=[1.0, 1.0, 1.0, 1.0],
    COW=[1.0, 1.0, 1.0, 1.0],
    KDIL_DAY=1500.0, KDIL_NIGHT=86400.0,
    BG_O3=30.0, BG_NO2=5.0, BG_NO=1.0, BG_CO=120.0, BG_VOC=20.0, BG_ISO=0.2,
    TMIN_K=288.0, TMAX_K=301.0,
    DUMMYF=True, CONSTJ=False, C_H2O=5.0e18,
)

def _fmt(v):
    if isinstance(v, bool):
        return ".TRUE." if v else ".FALSE."
    if isinstance(v, str):
        return v if v.startswith("'") else "'%s'" % v
    if isinstance(v, (list, tuple)):
        return ",".join(_fmt(x) for x in v)
    if isinstance(v, float):
        return repr(v)
    return str(v)

def nml_text(cfg, backend=None):
    """RUNCFG namelist text for a scenario config, with the backend's keys."""
    backend = backend or DEFAULT_BACKEND
    p = dict(DEFAULTS)
    p.update(cfg)
    keys = BASE_KEYS + EXTRA_KEYS.get(backend, [])
    return "&RUNCFG\n" + "".join("  %s = %s,\n" % (k, _fmt(p[k])) for k in keys) + "/\n"


def run_scenario(cfg, name, backend=None, attempts=3):
    """cfg: dict overriding DEFAULTS (RUNCFG names). Returns metrics dict.

    backend: "split" (/root/kpp4) or "dummy" (/root/kpp5); defaults to
    $KPP_BACKEND or "split".
    Retries a few times: killing a run mid-flight can leave a half-written
    traj.dat on the Windows-side mount, which is worth one more attempt in a
    256-run sweep rather than losing the whole sweep.
    """
    last = None
    for _ in range(attempts):
        try:
            return _run_once(cfg, name, backend)
        except (RuntimeError, ValueError, OSError) as e:
            last = e
    raise RuntimeError("run %s failed %d times: %s" % (name, attempts, last))


def _run_once(cfg, name, backend=None):
    backend = backend or DEFAULT_BACKEND
    wsl_dir = BACKENDS[backend]
    runs = os.path.join(HERE, "g4runs", backend)
    mnt_runs = "/mnt/d/atmos/_kppwork/g4runs/" + backend
    os.makedirs(runs, exist_ok=True)
    for stale in (name + ".traj", name + ".state", name + ".log"):
        try:
            os.remove(os.path.join(runs, stale))
        except OSError:
            pass

    p = dict(DEFAULTS)
    p.update(cfg)
    nml = os.path.join(HERE, "cb05.nml")
    with open(nml, "w") as f:
        f.write(nml_text(cfg, backend))

    script = (
        "set -e; cp /mnt/d/atmos/_kppwork/cb05.nml {d}/cb05.nml; "
        "cd {d}; rm -f traj.dat state_out.dat; timeout 300 ./cb05.exe > run.log 2>&1 || {{ "
        "echo RUNFAIL; tail -3 run.log; exit 1; }}; "
        "mkdir -p {m}; "
        "cp traj.dat {m}/{n}.traj; cp state_out.dat {m}/{n}.state; "
        "cp run.log {m}/{n}.log"
    ).format(d=wsl_dir, m=mnt_runs, n=name)
    r = subprocess.run(["wsl.exe", "bash", "-lc", script],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("run %s failed rc=%d: out=%s err=%s" %
                           (name, r.returncode, r.stdout[-300:], r.stderr[-300:]))

    traj = load_traj(os.path.join(runs, name + ".traj"))
    state = load_state(os.path.join(runs, name + ".state"))
    log = open(os.path.join(runs, name + ".log"), errors="replace").read()
    if "DONE" not in log:
        raise RuntimeError("run %s incomplete: %s" % (name, log[-300:]))
    nfail = int([ln for ln in log.splitlines() if "NFAIL=" in ln][0].split("=")[1])
    if len(traj) < 100:
        raise RuntimeError("run %s: traj only %d rows" % (name, len(traj)))
    m = metrics(traj)
    m["name"] = name
    m["nfail"] = nfail
    m["cfg"] = p
    m["state"] = state
    m["traj"] = traj          # diurnal() needs the full series, not just scalars
    return m

def load_traj(path):
    rows = []
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            rows.append([float(x) for x in ln.split()])
    return rows  # cols: t_h o3 nox temp sday kdil emis  [anchor adds co no no2 par tol xyl isop]


# The "anchor" backend appends the species we need to compare against AQS.
# Columns beyond index 6 only exist there, so this is guarded.
ANCHOR_COLS = ["co", "no", "no2", "par", "tol", "xyl", "isop"]


def diurnal(traj, col, spinup_h=48.0):
    """Mean of column `col` by local hour, using only t >= spinup_h.

    The box is driven by a repeating diurnal cycle, so the meaningful comparison
    with a station is the 24-h composite, not a day-by-day series.
    """
    acc = {h: [] for h in range(24)}
    idx = 1 if col == "o3" else (2 if col == "nox" else 7 + ANCHOR_COLS.index(col))
    for r in traj:
        if len(r) <= idx:
            continue
        t = r[0]
        if t < spinup_h:
            continue
        h = int(round(t)) % 24
        acc[h].append(r[idx])
    return {h: (sum(v) / len(v) if v else float("nan")) for h, v in acc.items()}

def _fort_float(s):
    # Fortran Ew.d may print exponents without 'E': 0.300167756076091-314
    s = s.strip()
    i = max(s.rfind("+"), s.rfind("-"))
    if i > 0 and "e" not in s.lower() and s[:i][-1:].isdigit():
        s = s[:i] + "e" + s[i:]
    return float(s)

def load_state(path):
    d = {}
    with open(path) as f:
        n = int(f.readline().split()[0])
        for _ in range(n):
            k, v = f.readline().split()
            d[k] = _fort_float(v)
    return d

def metrics(traj, spinup_h=24.0, dtol=0):
    o3 = [r[1] for r in traj]
    t = [r[0] for r in traj]
    n = len(o3)
    sel = [i for i in range(n) if t[i] >= spinup_h]
    day3 = [i for i in range(n) if t[i] >= 48.0]
    m = {
        "mean_all": _mean([o3[i] for i in sel]),
        "max_all": max(o3[i] for i in sel),
        "d3_mean": _mean([o3[i] for i in day3]),
        "d3_max": max(o3[i] for i in day3),
    }
    # MDA8: rolling 8h mean, max over t>=24h. assumes uniform dt
    if n > 2:
        dt = t[1] - t[0]
        w = max(1, int(round(8.0 / dt)))
        best = 0.0
        for i in range(w - 1, n):
            if t[i] < spinup_h:
                continue
            a = _mean(o3[i - w + 1 : i + 1])
            best = max(best, a)
        m["mda8"] = best
    return m

def _mean(xs):
    return sum(xs) / len(xs) if xs else float("nan")

if __name__ == "__main__":
    import sys
    be = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_BACKEND
    m = run_scenario({}, "smoke_default", be)
    print(be, {k: v for k, v in m.items() if k not in ("cfg", "state")})
