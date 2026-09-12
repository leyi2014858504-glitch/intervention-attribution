"""CB05 真实气相机制后端:VOC-limited / EKMA 脊线确证(慢测试)。

`pytest -m "not slow"` 可跳过。断言的是**符号结构**,绝对量级已标定到现实区
(基准 NOx 平衡 ≈79 ppb、O3 日峰 ≈115 ppb)。
"""
import numpy as np
import pytest

from atmos.configs import CB05Config
from atmos.chemistry.cb05 import CB05Backend

pytestmark = pytest.mark.slow

PYKPP_T0 = 9 * 3600.0


def _pbl(t_hr, h_min=250.0, h_max=1800.0):
    hod = (PYKPP_T0 / 3600.0 + t_hr) % 24.0
    s = max(0.0, np.sin(np.pi * (hod - 6.0) / 12.0))
    return h_min + (h_max - h_min) * np.sqrt(s)


def _cfg():
    return CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2)


@pytest.fixture(scope="module")
def backend():
    return CB05Backend(_cfg())


def _peak_o3(backend, nox_s, voc_s, n=12, dt=6.0):
    """给定排放比值下的 72h 日峰 O3(ppb)。直调 backend.step,不经 EnvModel 缓存。"""
    backend.set_scenario(nox=nox_s, voc=voc_s)
    y = backend.initial_state().copy()
    met = {"h_pbl_fn": _pbl, "temp": 298.15}
    peak = 0.0
    for i in range(n):
        res = backend.step(y, i * dt, dt, action=np.zeros(4), met=met)
        y = res.conc_out
        peak = max(peak, backend.project_pollutants(y)["O3"])
    return peak


def test_backend_constructs_and_projects():
    b = CB05Backend(_cfg())
    assert b.n_species == len(b.species) >= 50
    assert {"NO", "NO2", "O3", "ISOP", "TOL", "SO2"} <= set(b.species)
    y0 = b.initial_state()
    assert y0.shape == (b.n_species,) and np.all(y0 >= 0)
    pol = b.project_pollutants(y0)
    assert pol["PM25"] == 0.0          # cb05cl 无气溶胶:如实返回 0
    assert pol["NOx"] > 0


def test_background_inflow_prevents_o3_collapse(backend):
    """背景入流必须存在 —— 否则高排放下物种无限累加、O3 被过量 NO 滴定归零。

    这是本项目真实发生过的 bug:稀释曾实现为单向漏(`X = DUMMY`)而无回补,
    导致基准 NOx 达 450 ppb(= 排放×驻留时间)、O3 在 0↔100 ppb 间振荡、奖励失效。
    """
    backend.set_scenario(nox=1.0, voc=1.0)
    y = backend.initial_state().copy()
    met = {"h_pbl_fn": _pbl, "temp": 298.15}
    o3 = []
    for i in range(12):
        res = backend.step(y, i * 6.0, 6.0, action=np.zeros(4), met=met)
        y = res.conc_out
        o3.append(backend.project_pollutants(y)["O3"])
    assert min(o3) > 5.0, f"O3 塌零(背景入流失效?): min={min(o3):.2f}"


def test_nox_reduction_raises_o3_on_voc_limited_side(backend):
    """基准(VOC-limited):减 40% NOx、VOC 不变 ⇒ 日峰 O3 上升。"""
    base = _peak_o3(backend, 1.0, 1.0)
    cut = _peak_o3(backend, 0.6, 1.0)
    assert cut > base, f"期望 VOC-limited:base={base:.1f} cut={cut:.1f}"


def test_cb05_step_handoff_carries_state(backend):
    """真 CB05 的跨步接力:第 0 步减排必须改变步末状态(否则环境构造性解耦)。

    实测相对差 ≈ 0.034;若掉到 1e-3 以下,先去查 rollout/后端的状态交接,
    不要做任何环境改造。
    """
    c0 = backend.initial_state()
    r_base = backend.step(c0, 0.0, 6.0, action=np.zeros(4)).conc_out
    r_pert = backend.step(c0, 0.0, 6.0, action=np.ones(4)).conc_out
    rel = np.linalg.norm(r_pert - r_base) / (np.linalg.norm(r_base) or 1.0)
    assert rel > 1e-3, f"接力断裂:动作未影响步末状态(rel={rel:.2e})"


def test_ekma_ridge_sign_flips_with_voc(backend):
    """固定 NOx:低 VOC 侧 VOC-limited(减NOx升O3);高 VOC 侧翻转为 NOx-limited。"""
    lo = _peak_o3(backend, 0.6, 0.7) - _peak_o3(backend, 1.0, 0.7)
    hi = _peak_o3(backend, 0.6, 2.5) - _peak_o3(backend, 1.0, 2.5)
    assert lo > 0, f"低VOC侧应 VOC-limited(Δ>0),实得 {lo:.1f}"
    assert hi < 0, f"高VOC侧应 NOx-limited(Δ<0),实得 {hi:.1f}"
