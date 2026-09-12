"""化学层单测(新契约):MockBackend 的 VOC-limited 非单调 + 数值 sanity。

VOC-limited 断言是全项目化学可信度第一锚:固定 VOC、降 NOx ⇒ 72h O3 上升。
"""
import numpy as np

from atmos.configs import EnvConfig
from atmos.chemistry.mock import MockBackend, MOCK_SPECIES

_I = {s: i for i, s in enumerate(MOCK_SPECIES)}


def _run(backend, action, y0=None, n_steps=12, dt=6.0):
    """以常动作推演 72h,返回步内轨迹时间均值 O3 与末态。"""
    h_min, h_max = backend.cfg.h_min, backend.cfg.h_max
    hod = lambda t: h_min + (h_max - h_min) * max(0.0, np.sin(np.pi * (t % 24 - 6) / 12)) ** 0.5
    met = {"h_pbl_fn": hod, "temp": 298.15}
    conc = backend.initial_state() if y0 is None else y0.copy()
    o3_acc = []
    for i in range(n_steps):
        res = backend.step(conc, i * dt, dt, action=action, met=met)
        conc = res.conc_out
        assert np.all(np.isfinite(res.conc_traj))
        assert np.all(res.conc_traj >= 0.0), "浓度出现负值"
        o3_acc.append(res.conc_traj[:, _I["O3"]].mean())
    return float(np.mean(o3_acc)), conc


def test_voc_limited_nox_reduction_raises_o3():
    """电力(偏 NOx、VOC≈0)减排 ⇒ 72h 均值 O3 上升(VOC-limited)。"""
    b = MockBackend(EnvConfig())
    o3_base, _ = _run(b, np.array([0.0, 0.0, 0.0, 0.0]))
    o3_cut, _ = _run(b, np.array([0.0, 0.0, 0.9, 0.0]))   # 电力 NOx 大减,VOC 不动
    assert o3_cut > o3_base, f"VOC-limited 不成立: base={o3_base:.1f} cut={o3_cut:.1f}"


def test_positive_bounded_no_nan():
    """零动作下各物种保持非负、有限、量级合理。"""
    b = MockBackend(EnvConfig())
    _, concf = _run(b, np.zeros(4))
    assert np.all(concf >= 0.0) and np.all(concf < 1e6)


def test_project_pollutants_keys():
    b = MockBackend(EnvConfig())
    pol = b.project_pollutants(b.initial_state())
    assert set(pol) == set(MOCK_SPECIES)
    assert pol["O3"] > 0
