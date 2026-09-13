#!/usr/bin/env python3
"""
Pull hourly meteorology from Open-Meteo's Historical Weather API (no key).

We need the boundary-layer height above all: the box model's dilution term
(and therefore the NOx lifetime) is driven by BLH, so anchoring the model to
observations without BLH is not possible. T2m / shortwave / wind / RH come
along for the photolysis and transport terms.

Usage:
    python3 _met_pull.py <lat> <lon> <start> <end> <label>

Writes /mnt/d/atmos/met_<label>.csv
"""
import sys, csv, json, urllib.request, os

ARCHIVE = "https://archive-api.open-meteo.com/v1/archive"
HOURLY = [
    "temperature_2m",
    "relative_humidity_2m",
    "surface_pressure",
    "wind_speed_10m",
    "wind_direction_10m",
    "shortwave_radiation",
    "boundary_layer_height",
]

def pull(lat, lon, start, end, label, tz="America/New_York", outdir="/mnt/d/atmos"):
    q = (
        f"{ARCHIVE}?latitude={lat}&longitude={lon}"
        f"&start_date={start}&end_date={end}"
        f"&hourly={','.join(HOURLY)}&timezone={tz}"
    )
    req = urllib.request.Request(q, headers={"User-Agent": "curl/8"})
    # bypass any local proxy: the sandbox proxy breaks TLS to these hosts too
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    with opener.open(req, timeout=60) as r:
        j = json.load(r)
    h = j["hourly"]
    out = os.path.join(outdir, f"met_{label}.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["time"] + HOURLY)
        for i, t in enumerate(h["time"]):
            w.writerow([t] + [h[k][i] for k in HOURLY])
    print(f"{label}: {len(h['time'])} rows -> {out}")
    return out

if __name__ == "__main__":
    if len(sys.argv) not in (6, 7):
        print(__doc__)
        sys.exit(1)
    tz = sys.argv[6] if len(sys.argv) == 7 else "America/New_York"
    pull(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], tz)
