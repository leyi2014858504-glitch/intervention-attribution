#!/bin/bash
# Audit verification, part 3: the cost-law claim's actual source, and what
# quantity the anchoring cost function fits.
set -u
cd /mnt/d/atmos

echo "=== A. everywhere the 1.6-1.8x / 2.5-3x numbers appear"
grep -rn "1\.6–1\.8\|1\.6-1\.8\|2\.5–3\|2\.5-3" paper_outline.md methods_zh.md README.md 2>/dev/null | head -12

echo
echo "=== B. outline Fig.4 / Table 5 source notes"
sed -n '149,240p' paper_outline.md

echo
echo "=== C. anchoring cost function body (_fit_anchor.py lines 86-140)"
sed -n '86,140p' _fit_anchor.py