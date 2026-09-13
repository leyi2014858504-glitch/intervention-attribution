#!/usr/bin/env python3
"""
Fetch the AQS email-confirmation link from the command line.

The confirmation link is a plain HTTP GET, so it can be opened with code and we
never touch the browser. Proxies are disabled because the local 127.0.0.1:7890
proxy is exactly what breaks TLS to aqs.epa.gov.

Usage:
  python3 _aqs_confirm.py "<paste the confirmation link here>"
"""
import sys, urllib.request, urllib.error, re

if len(sys.argv) < 2:
    print(__doc__)
    sys.exit(1)

url = sys.argv[1].strip().strip('"').strip("'")

# paste from a mail client often brings angle brackets or line breaks
url = url.strip("<>").replace("\n", "").replace("\r", "")

OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))
print("GET", url)
try:
    with OPENER.open(url, timeout=90) as r:
        body = r.read().decode("utf-8", "replace")
    print("HTTP", r.status, "bytes", len(body))
    text = re.sub(r"<script.*?</script>|<style.*?</style>", " ", body, flags=re.S | re.I)
    text = re.sub(r"<[^>]+>", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    print("---- readable text ----")
    print(text[:2000])
    m = re.findall(r"[A-Za-z][A-Za-z]+[0-9]{1,4}", body)
    if m:
        print("---- key-shaped tokens found on the page ----")
        print(sorted(set(m))[:20])
except urllib.error.HTTPError as e:
    print("HTTPError", e.code)
    print(e.read().decode("utf-8", "replace")[:2000])
except Exception as e:
    print("FAILED:", repr(e))
