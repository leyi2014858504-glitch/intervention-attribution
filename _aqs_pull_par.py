#!/usr/bin/env python3
"""
Parallel AQS puller.

Three things learned the hard way, all baked in here:

1. Bandwidth. aqs.epa.gov serves us at roughly 10 KB/s per connection, and one
   species for a whole month is ~587 KB -> 30-100 s. Sequentially that is hours.
   Six concurrent connections gave a 4.7x speedup (121 s -> 26 s), so we use a
   thread pool.

2. Species count. Beaumont's inventory lists 342 VOC parameters, but most are
   Air-Toxics halocarbons, pesticides and aggregates with no CB05 counterpart.
   We reuse the classifier from _voc_map.py and keep only parameters that map to
   a real CB05 group, plus the four criteria species.

3. Window length. Cost scales with bytes = species x days. 14 days per season
   still gives 336 hourly points per species, which is plenty for inverse
   fitting, and cuts transfer by more than half versus a full month.

Usage:
  python3 _aqs_pull_par.py <state> <county> <site> <bdate> <edate> <label> [workers]
"""
import sys, os, csv, json, time, threading, importlib.util
import urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
CACHE = "/mnt/d/atmos/_site_param_cache.json"
OUTDIR = "/mnt/d/atmos/aqs_sites_data"
API = "https://aqs.epa.gov/data/api/sampleData/bySite"
EMAIL = "leyi2014858504@gmail.com"
KEY = "ochrehare24"

CORE = ["42101", "42601", "42602", "44201"]   # CO, NO, NO2, O3

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
_print_lock = threading.Lock()


def load_classifier():
    spec = importlib.util.spec_from_file_location("vocmap", "/mnt/d/atmos/_voc_map.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.classify


def site_params(state, county, site):
    """code -> name for one site, cached so the 160 MB inventory is read once."""
    if os.path.exists(CACHE):
        blob = json.load(open(CACHE))
    else:
        blob = {}
        with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
            for row in csv.DictReader(fh):
                k = f"{row['State Code']}-{row['County Code']}-{row['Site Number']}"
                blob.setdefault(k, {})[row["Parameter Code"]] = (row.get("Parameter Name") or "").strip()
        json.dump(blob, open(CACHE, "w"))
    return blob.get(f"{state}-{county}-{site}", {})


def fetch(param, state, county, site, bdate, edate, tries=4):
    q = (f"{API}?email={EMAIL}&key={KEY}&param={param}&bdate={bdate}&edate={edate}"
         f"&state={state}&county={county}&site={site}")
    last = None
    for k in range(tries):
        try:
            with OPENER.open(q, timeout=180) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
        except Exception as e:
            last = repr(e)
        time.sleep(2 + 3 * k)
    return {"Header": [{"status": "FAILED", "error": [last]}], "Data": []}


def main():
    state, county, site, bdate, edate, label = sys.argv[1:7]
    workers = int(sys.argv[7]) if len(sys.argv) > 7 else 10
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"site_{label}.csv")
    logf = os.path.join(OUTDIR, f"site_{label}.log")
    log = open(logf, "w", encoding="utf-8", buffering=1)

    classify = load_classifier()
    names = site_params(state, county, site)

    keep, dropped = [], []
    for pc, nm in sorted(names.items()):
        if pc in CORE:
            keep.append(pc)
            continue
        grp, conf, note = classify(nm)
        if grp in ("-", "?"):
            dropped.append((pc, nm, grp))
        else:
            keep.append(pc)

    msg = (f"{label}: {len(names)} parameters at {state}-{county}-{site} -> "
           f"keeping {len(keep)} (criteria + real CB05 species), dropping {len(dropped)} "
           f"(halocarbons/PAH/aggregates). workers={workers}, window {bdate}..{edate}")
    print(msg, flush=True); log.write(msg + "\n")

    cols = ["parameter_code", "parameter", "poc", "date_local", "time_local",
            "sample_measurement", "units_of_measure", "sample_duration",
            "method", "latitude", "longitude"]
    rows = []
    t0 = time.time()
    done = 0

    def work(pc):
        return pc, fetch(pc, state, county, site, bdate, edate)

    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(work, pc): pc for pc in keep}
        for f in as_completed(futs):
            pc, j = f.result()
            done += 1
            hdr = (j.get("Header") or [{}])[0]
            data = j.get("Data") or []
            hourly = [d for d in data if (d.get("sample_duration") or "").upper() == "1 HOUR"]
            for d in hourly:
                rows.append(d)
            with _print_lock:
                line = (f"  [{done}/{len(keep)}] {pc} {hdr.get('status','?'):<32} "
                        f"rows={len(data):<5} hourly={len(hourly):<5} "
                        f"elapsed={time.time()-t0:.0f}s")
                print(line, flush=True); log.write(line + "\n")

    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    tail = f"{label}: wrote {len(rows)} hourly rows -> {out}  ({time.time()-t0:.0f}s)"
    print(tail, flush=True); log.write(tail + "\n")
    log.close()


if __name__ == "__main__":
    main()
