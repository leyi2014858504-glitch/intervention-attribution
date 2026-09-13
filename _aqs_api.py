#!/usr/bin/env python3
"""
Pull one AQS site's hourly data through the official API.

Why the API rather than the bulk files: the bulk "hourly_<param>_<year>.zip"
set simply does not publish nitric oxide (42601 -> HTTP 404), and NO vs NO2 as
SEPARATE hourly series is the whole point of this study. The API also returns
one site instead of the whole country (70 MB for O3 alone in bulk).

Which parameters get requested is not guessed: we intersect the site's own
parameter list from aqs_monitors.csv with (core species + the 43xxx VOC panel),
then keep only rows whose sample_duration is "1 HOUR".

Usage:
  python3 _aqs_api.py <email> <key> <state> <county> <site> <bdate> <edate> <label>

Example:
  python3 _aqs_api.py me@x.com ABC123 06 037 1103 20161201 20161231 la_2016w
"""
import sys, os, csv, json, time, urllib.request, urllib.error

D = "/mnt/d/atmos/data_aqs"
MON = os.path.join(D, "aqs_monitors.csv")
OUTDIR = "/mnt/d/atmos/aqs_sites_data"
API = "https://aqs.epa.gov/data/api/sampleData/bySite"

CORE = ["44201", "42601", "42602", "42101"]   # O3, NO, NO2, CO


def is_voc(pc):
    """43xxx (hydrocarbons/oxygenates) plus 45xxx (the AROMATICS).

    452xx and 451xx are easy to miss: benzene is 45201, toluene 45202,
    45109 is m/p-xylene, 45220 styrene. They are not optional -- toluene and
    xylene dominate the aromatic reactivity in a VOC-limited urban regime.
    """
    return (pc.startswith("43") and pc != "43000") or pc.startswith("45")

# Build the opener once, with proxies disabled: the local 127.0.0.1:7890 proxy
# breaks TLS to aqs.epa.gov (see _aqs_probe2.txt).
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def site_params(state, county, site):
    """Parameter codes this site actually reports (from the monitor inventory)."""
    want = set()
    with open(MON, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.DictReader(fh):
            if (row["State Code"], row["County Code"], row["Site Number"]) != (state, county, site):
                continue
            pc = row["Parameter Code"]
            if pc in CORE or is_voc(pc):
                want.add(pc)
    return sorted(want)


def fetch(email, key, state, county, site, param, bdate, edate, tries=4):
    q = (f"{API}?email={email}&key={key}&param={param}"
         f"&bdate={bdate}&edate={edate}"
         f"&state={state}&county={county}&site={site}")
    for k in range(tries):
        try:
            with OPENER.open(q, timeout=90) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            # 429 = throttled; back off
            time.sleep(3 * (k + 1))
        except Exception:
            time.sleep(3 * (k + 1))
    return None


def main():
    email, key, state, county, site, bdate, edate, label = sys.argv[1:9]
    os.makedirs(OUTDIR, exist_ok=True)
    out = os.path.join(OUTDIR, f"site_{label}.csv")
    logf = os.path.join(OUTDIR, f"site_{label}.log")
    log = open(logf, "w", encoding="utf-8")

    params = site_params(state, county, site)
    print(f"{label}: {len(params)} parameters at {state}-{county}-{site}")
    log.write(f"parameters requested: {params}\n")

    cols = None
    nrow = 0
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = None
        for pc in params:
            j = fetch(email, key, state, county, site, pc, bdate, edate)
            if j is None:
                log.write(f"{pc} FETCH FAILED\n")
                continue
            hdr = j.get("Header", [{}])[0]
            if hdr.get("status") == "No data matched your selection":
                log.write(f"{pc} no data\n")
                continue
            data = j.get("Data") or []
            if not data:
                log.write(f"{pc} empty\n")
                continue
            # keep hourly only: anything else has already lost the diurnal cycle
            hourly = [d for d in data if (d.get("sample_duration") or "").upper() == "1 HOUR"]
            log.write(f"{pc} {data[0].get('parameter')}: {len(data)} rows, {len(hourly)} hourly\n")
            for d in hourly:
                if w is None:
                    cols = ["parameter_code", "parameter", "poc", "date_local", "time_local",
                            "sample_measurement", "units_of_measure", "sample_duration",
                            "method", "latitude", "longitude"]
                    w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
                    w.writeheader()
                w.writerow(d)
                nrow += 1
            time.sleep(0.6)
    log.write(f"TOTAL hourly rows written: {nrow}\n")
    log.close()
    print(f"{label}: wrote {nrow} hourly rows -> {out}")


if __name__ == "__main__":
    main()
