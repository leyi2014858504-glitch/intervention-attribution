"""M1 实验:CB05 NOx×VOC 响应面(EKMA)+ 等值线图 + 脊线检测。

跑:  python -m atmos.experiments.exp_rsm
产出: results/rsm/rsm_cb05.json  +  results/rsm/ekma_cb05.png
"""
from __future__ import annotations

import dataclasses
import json
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from ..chemistry.cb05 import CB05Backend
from ..configs import CB05Config, EnvConfig
from ..scenarios.rsm import nox_voc_grid, build_response_surface

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "rsm")


def detect_ridge(nox_axis, voc_axis, peak):
    """对每个 VOC 水平,找使峰值 O3 最大的 NOx → 脊线 NOx*(VOC)。"""
    peak = np.asarray(peak)
    na = np.asarray(nox_axis)
    ridge = []
    for j in range(peak.shape[0]):
        i = int(np.argmax(peak[j]))
        ridge.append(float(na[i]))
    return ridge


def main(n: int = 7, n_steps: int = 6):
    os.makedirs(OUT_DIR, exist_ok=True)
    cfg = dataclasses.replace(CB05Config(), rtol=1e-3, atol=1e6,
                              substeps_per_control=2)
    backend = CB05Backend(cfg)
    nox_axis, voc_axis = nox_voc_grid(n=n)

    t0 = time.perf_counter()
    surf = build_response_surface(backend, nox_axis, voc_axis,
                                  dt=6.0, n_steps=n_steps)
    elapsed = time.perf_counter() - t0

    na = np.asarray(surf["nox_axis"])
    va = np.asarray(surf["voc_axis"])
    peak = np.asarray(surf["o3_daily_peak"])
    ridge_nox = detect_ridge(na, va, surf["o3_daily_peak"])

    payload = {
        "config": {"nox_base_ppbh": cfg.nox_base, "voc_base_ppbh": cfg.voc_base,
                   "vent_tau_h": cfg.vent_tau_h, "h_min": cfg.h_min,
                   "h_max": cfg.h_max, "so2_base_ppbh": cfg.so2_base},
        "env": dataclasses.asdict(EnvConfig()),
        "surface": surf,
        "ridge_nox_star_over_voc": ridge_nox,
        "elapsed_s": elapsed,
        "note": "CB05 气相,常数稀释;PM2.5 不在本机制内(见 use_pm_penalty)。",
    }
    jpath = os.path.join(OUT_DIR, "rsm_cb05.json")
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)

    # 等值线图(EKMA):横 NOx 比、纵 VOC 比、色 = 日峰 O3
    NA, VA = np.meshgrid(na, va)
    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    cf = ax.contourf(NA, VA, peak, levels=18, cmap="viridis")
    cs = ax.contour(NA, VA, peak, levels=8, colors="k", linewidths=0.6)
    ax.clabel(cs, inline=True, fontsize=7, fmt="%.0f")
    ax.plot(ridge_nox, va, "r-o", ms=4, lw=1.5, label="EKMA ridge (argmax NOx)")
    ax.axhline(1.0, color="w", ls="--", lw=0.7)
    ax.axvline(1.0, color="w", ls="--", lw=0.7)
    fig.colorbar(cf, ax=ax, label="O3 daily-peak (ppb)")
    ax.set_xlabel("NOx emission ratio"); ax.set_ylabel("VOC emission ratio")
    ax.set_title("CB05 EKMA response surface (baseline = 1.0/1.0)")
    ax.legend(loc="upper right", fontsize=8)
    ppath = os.path.join(OUT_DIR, "ekma_cb05.png")
    fig.tight_layout(); fig.savefig(ppath, dpi=140); plt.close(fig)

    print(f"cells={len(surf['cells'])} elapsed={elapsed:.1f}s")
    print(f"base(1,1) peak O3 = {peak[np.argmin(abs(na-1)), np.argmin(abs(va-1))]:.1f} ppb")
    print(f"ridge NOx*(VOC) over voc_axis={va.tolist()}:\n   {ridge_nox}")
    print(f"-> {jpath}\n-> {ppath}")


if __name__ == "__main__":
    main()
