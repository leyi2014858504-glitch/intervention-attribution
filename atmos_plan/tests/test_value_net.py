"""M3-1 value-net 截断测试(MockBackend,horizon=12 贴近真实 CB05)。

核心闸门断言:同预算下开启 value-net 后**箱模型调用数显著下降**,且计划质量不塌。
另含一条实测得到的结构性结论:0<w<1 的"混合"**不省算力**(仍要跑完整 rollout),
只有 w=1(彻底砍掉 rollout)才省 —— 中间段只对 M6 鲁棒性实验有意义。
"""
import numpy as np
import pytest

from atmos.configs import EnvConfig
from atmos.chemistry.mock import MockBackend
from atmos.search.model import EnvModel
from atmos.search.mcts import MCTS, MCTSConfig
from atmos.search.value_net import ValueNet


def _model():
    # horizon=12 与真实 CB05 一致;低 O3 阈值使超标有信号
    cfg = EnvConfig(o3_limit=50.0, horizon_steps=12, dt_hours=6.0)
    return EnvModel(MockBackend(cfg), cfg)


def _run(*, use_vnet: bool, budget: int, seed: int = 1, **over):
    model = _model()
    root = model.backend.initial_state()
    m = MCTS(model, MCTSConfig(budget=budget, seed=seed,
                               use_value_net=use_vnet, **over))
    res = m.search(root)
    plan = m.best_path(res.root)
    fresh = _model()                       # 独立评估,不受搜索缓存影响
    return fresh.rollout_sequence(root, plan), model.box_calls, res


# ---- 单元:ValueNet 本体 ----

def test_value_net_needs_min_samples():
    vn = ValueNet(MockBackend(), min_samples=5, horizon=12)
    assert not vn.ready
    vn.add_sample(np.ones(5), 0, -1.0)
    assert not vn.fit()                    # 样本不足 ⇒ 不就绪
    assert not vn.ready


def test_value_net_fits_and_predicts_finite():
    model = _model()
    vn = ValueNet(model.backend, min_samples=20, horizon=model.horizon)
    rng = np.random.default_rng(0)
    c = model.backend.initial_state()
    for _ in range(40):
        c2, _, _ = model.step(c, 0, rng.random(4))
        vn.add_sample(c2, 1, model.rollout(c2, 1, rng=rng))
    assert vn.fit() and vn.ready
    assert np.isfinite(vn.predict(c, 0))


# ---- 闸门:截断确实省算力且质量不塌 ----

def test_value_net_cuts_box_calls():
    base_val, base_calls, _ = _run(use_vnet=False, budget=200)
    vn_cfg = dict(value_w=1.0, warmup_iters=60, retrain_every=40,
                  label_fraction=0.2, vnet_min_samples=24)
    vn_val, vn_calls, res = _run(use_vnet=True, budget=200, **vn_cfg)

    assert res.vnet_ready and res.vnet_samples >= 24
    assert vn_calls < 0.7 * base_calls, (
        f"截断未显著省算力: base={base_calls} vnet={vn_calls}")
    # 回报 ≤0,越接近 0 越好 ⇒ 用 |value| 比较(不劣于基线 15%)
    assert abs(vn_val) <= 1.15 * abs(base_val), (
        f"质量塌了: base={base_val:.4f} vnet={vn_val:.4f}")


def test_partial_w_saves_no_calls():
    """0<w<1 仍执行完整 rollout ⇒ 省不到算力(实测结论,防后人误用)。"""
    _, base_calls, _ = _run(use_vnet=False, budget=150)
    _, w05_calls, _ = _run(use_vnet=True, budget=150, value_w=0.5,
                           warmup_iters=40, retrain_every=40)
    assert w05_calls >= 0.85 * base_calls


def test_value_net_deterministic_same_seed():
    vn_cfg = dict(value_w=1.0, warmup_iters=40, retrain_every=20, label_fraction=0.2)
    v1, c1, _ = _run(use_vnet=True, budget=120, seed=3, **vn_cfg)
    v2, c2, _ = _run(use_vnet=True, budget=120, seed=3, **vn_cfg)
    assert v1 == pytest.approx(v2)
    assert c1 == c2


def test_disabled_vnet_matches_pure_rollout():
    """use_value_net=False 必须与 value_w=0 行为一致(消融可对照)。"""
    a_val, a_calls, a_res = _run(use_vnet=False, budget=100, seed=5)
    b_val, b_calls, _ = _run(use_vnet=True, budget=100, seed=5, value_w=0.0)
    assert a_res.vnet_samples == 0 and a_res.vnet_fits == 0
    assert a_calls == b_calls
    assert a_val == pytest.approx(b_val)
