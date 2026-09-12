"""CB05 后端的**确定性/可复现性**回归测试。

背景(已实测的严重问题):pykpp 的 `stdfuncs.update_func_world` 会执行
`globals().update(world)`,把 world 的内容写进**模块级全局**。于是同一个进程内:
  - 同 seed 的结果依赖**调用历史**;
  - 甚至两个全新 backend 实例的结果也不一致。

后果:后端 `step()` **不是纯函数**,与 `ChemistryBackend` 契约("无副作用 + 确定")不符,
记忆化、精确 oracle、跨实验可复现性都会被污染。

修复方向:进程隔离(每个后端/worker 独占进程),或在 pykpp 层阻断全局写入。
修复前本测试标 xfail(修复后会自动变为 XPASS 提醒更新)。
"""
import numpy as np
import pytest

from atmos.configs import CB05Config, EnvConfig
from atmos.chemistry.cb05 import CB05Backend
from atmos.search.model import EnvModel

pytestmark = pytest.mark.slow

_W = dict(rtol=1e-3, atol=1e6, substeps_per_control=2,
          vent_tau_logsigma=0.6, pbl_hmax_logsigma=0.45)


def _peak(backend, n=12):
    m = EnvModel(backend, EnvConfig(horizon_steps=n))
    c = backend.initial_state()
    pk = 0.0
    for i in range(n):
        c, _r, _d = m.step(c, i, np.zeros(4))
        pk = max(pk, backend.project_pollutants(c)["O3"])
    return pk


_XFAIL = pytest.mark.xfail(
    reason="已知缺陷:pykpp 把 world 写入模块级全局 ⇒ step 非纯函数(需进程隔离修复)",
    strict=False)


@_XFAIL
def test_backend_step_is_pure_wrt_history():
    """同一 backend、同一 seed:在中间插入别的 episode 后,结果必须不变。"""
    b = CB05Backend(CB05Config(**_W))
    b.set_episode(1)
    first = _peak(b)
    b.set_episode(2)
    _peak(b)
    b.set_episode(1)
    again = _peak(b)
    assert first == pytest.approx(again), (
        f"step 依赖调用历史: {first:.4f} vs {again:.4f}(pykpp 全局污染?)")


@_XFAIL
def test_two_fresh_backends_agree():
    """同 seed、同配置的两个全新实例必须给出相同结果。"""
    a = CB05Backend(CB05Config(**_W)); a.set_episode(1)
    b = CB05Backend(CB05Config(**_W)); b.set_episode(1)
    assert _peak(a) == pytest.approx(_peak(b)), "两个干净实例结果不一致(进程级共享状态)"
