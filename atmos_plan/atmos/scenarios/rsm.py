"""响应面方法(RSM)情景:NOx × VOC 排放比网格 → O3 响应面(EKMA)。

因已锁纯气相 CB05(无 NH3/气溶胶),四变量 RSM 自然收敛为经典 **EKMA**:
横轴 NOx 比、纵轴 VOC 比(相对标定基准的比值),响应 = 72h 日峰 / 日均 O3。
这张面也是后续归因实验的"答案参照"。
"""
from __future__ import annotations

import dataclasses
from typing import Sequence

import numpy as np

from ..chemistry.cb05 import CB05Backend
from ..configs import CB05Config

PYKPP_T0 = 9 * 3600.0          # 积分从当日 09:00 起(有日照)


def diurnal_pbl(t_hr: float, h_min: float = 250.0, h_max: float = 1800.0) -> float:
    """边界层高度日循环(米);CB05 现用常数稀释,此函数供 met 传入与后续接入。"""
    hod = (PYKPP_T0 / 3600.0 + t_hr) % 24.0
    s = max(0.0, np.sin(np.pi * (hod - 6.0) / 12.0))
    return h_min + (h_max - h_min) * np.sqrt(s)


def nox_voc_grid(n: int = 9, lo: float = 0.1, hi: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    """NOx 比 × VOC 比网格坐标(返回一维轴,调用方 meshgrid)。"""
    axis = np.linspace(lo, hi, n)
    return axis, axis


def run_scenario(backend: CB05Backend, nox_r: float, voc_r: float,
                 *, dt: float = 6.0, n_steps: int = 12,
                 action: np.ndarray | None = None) -> dict:
    """给定 NOx/VOC 比值,常动作推演 72h,返回 O3 指标(ppb)。

    复用同一 backend(仅 set_scenario 改 world 排放标量),零重建成本。
    """
    backend.set_scenario(nox=nox_r, voc=voc_r)
    a = np.zeros(4) if action is None else np.asarray(action, float)
    y = backend.initial_state().copy()
    met = {"h_pbl_fn": lambda t: diurnal_pbl(t), "temp": backend.cfg.temp_K}
    o3_peak = 0.0
    o3_step_means = []
    for i in range(n_steps):
        res = backend.step(y, i * dt, dt, action=a, met=met)
        y = res.conc_out
        traj_o3 = [backend.project_pollutants(row)["O3"] for row in res.conc_traj]
        o3_peak = max(o3_peak, max(traj_o3))
        o3_step_means.append(float(np.mean(traj_o3)))
    return {
        "NOx_ratio": float(nox_r), "VOC_ratio": float(voc_r),
        "o3_daily_peak_ppb": float(o3_peak),
        "o3_mean_ppb": float(np.mean(o3_step_means)),
        "final_NOx_ppb": backend.project_pollutants(y)["NOx"],
    }


def build_response_surface(backend: CB05Backend, nox_axis: Sequence[float],
                           voc_axis: Sequence[float], *, verbose: bool = True,
                           **kw) -> dict:
    """扫 NOx×VOC 网格,返回可 json 化的响应面(含 meshgrid 轴与峰/均 O3 矩阵)。"""
    na = np.asarray(nox_axis, float)
    va = np.asarray(voc_axis, float)
    peak = np.zeros((len(va), len(na)))
    mean = np.zeros((len(va), len(na)))
    cells = []
    for j, vr in enumerate(va):
        for i, nr in enumerate(na):
            r = run_scenario(backend, nr, vr, **kw)
            peak[j, i] = r["o3_daily_peak_ppb"]
            mean[j, i] = r["o3_mean_ppb"]
            cells.append(r)
            if verbose:
                print(f"  [{i+1}/{len(na)},{j+1}/{len(va)}] NOx={nr:.2f} VOC={vr:.2f} "
                      f"-> peakO3={r['o3_daily_peak_ppb']:.1f} ppb", flush=True)
    return {
        "nox_axis": na.tolist(), "voc_axis": va.tolist(),
        "o3_daily_peak": peak.tolist(), "o3_mean": mean.tolist(),
        "cells": cells,
    }
