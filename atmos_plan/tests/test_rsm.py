"""RSM 响应面模块测试(快速、小网格)。

只锁**稳健**性质(角点对比 + 可序列化),不锁噪声敏感的 argmax 脊线。
真实 7×7 图由 `python -m atmos.experiments.exp_rsm` 生成。
"""
import json

import numpy as np
import pytest

from atmos.configs import CB05Config
from atmos.chemistry.cb05 import CB05Backend
from atmos.scenarios.rsm import nox_voc_grid, build_response_surface, run_scenario

pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def backend():
    cfg = CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2)
    return CB05Backend(cfg)


def test_scenario_returns_finite_metrics(backend):
    r = run_scenario(backend, 1.0, 1.0, dt=6.0, n_steps=3)
    assert np.isfinite(r["o3_daily_peak_ppb"]) and r["o3_daily_peak_ppb"] > 0
    assert set(r) >= {"NOx_ratio", "VOC_ratio", "o3_daily_peak_ppb", "o3_mean_ppb"}


def test_ekma_corner_contrast(backend):
    """EKMA 稳健特征:低NOx/高VOC 的 O3 峰 明显高于 高NOx/低VOC。"""
    hi_voc = run_scenario(backend, 0.1, 2.0, dt=6.0, n_steps=4)["o3_daily_peak_ppb"]
    hi_nox = run_scenario(backend, 2.0, 0.1, dt=6.0, n_steps=4)["o3_daily_peak_ppb"]
    assert hi_voc > hi_nox, f"EKMA 角点对比失败: {hi_voc:.1f} vs {hi_nox:.1f}"


def test_response_surface_shape_and_serializable(backend):
    na, va = nox_voc_grid(n=3)
    surf = build_response_surface(backend, na, va, verbose=False,
                                  dt=6.0, n_steps=3)
    peak = np.asarray(surf["o3_daily_peak"])
    assert peak.shape == (len(va), len(na))
    assert np.all(np.isfinite(peak)) and np.all(peak > 0)
    json.dumps(surf)   # 可 json 化(实验落盘前提)
