"""M2 MCTS 测试(在 MockBackend 上,毫秒级)。

核心验收:MCTS 贪心计划的开环回报 ≥ 随机策略 95 分位;确定性;记忆化生效。
"""
import numpy as np
import pytest

from atmos.configs import EnvConfig
from atmos.chemistry.mock import MockBackend
from atmos.search.model import EnvModel
from atmos.search.mcts import MCTS, MCTSConfig


def _model():
    # 低 O3 阈值 ⇒ 超标有信号,便于检验 MCTS 是否优于随机
    cfg = EnvConfig(o3_limit=50.0, horizon_steps=6, dt_hours=6.0)
    return EnvModel(MockBackend(cfg), cfg)


def _rand_seq(rng, n, d):
    return [rng.random(d) for _ in range(n)]


def test_mcts_beats_random_95th():
    model = _model()
    root = model.backend.initial_state()
    m = MCTS(model, MCTSConfig(budget=400, seed=1))
    res = m.search(root)
    plan = m.best_path(res.root)
    val_mcts = model.rollout_sequence(root, plan)

    rng = np.random.default_rng(0)
    rand_vals = [model.rollout_sequence(root, _rand_seq(rng, model.horizon, 4))
                 for _ in range(30)]
    p95 = np.percentile(rand_vals, 95)
    assert val_mcts >= p95, f"MCTS={val_mcts:.3f} < 随机95分位={p95:.3f}"


def test_mcts_root_visit_count():
    model = _model()
    root = model.backend.initial_state()
    m = MCTS(model, MCTSConfig(budget=100, seed=2))
    res = m.search(root)
    assert res.root.N == 100
    assert res.best_action.shape == (4,)
    assert np.all((res.best_action >= 0) & (res.best_action <= 1))


def test_mcts_determinism_same_seed():
    def run(seed):
        model = _model()
        root = model.backend.initial_state()
        m = MCTS(model, MCTSConfig(budget=80, seed=seed))
        res = m.search(root)
        return np.asarray(m.best_path(res.root)), model.rollout_sequence(root, m.best_path(res.root))
    s1, v1 = run(7)
    s2, v2 = run(7)
    assert v1 == pytest.approx(v2)
    assert all(np.allclose(a, b) for a, b in zip(s1, s2))


def test_model_memoization_hit():
    model = _model()
    c = model.backend.initial_state()
    a = np.full(4, 0.3)
    model.step(c, 0, a)
    model.step(c, 0, a)          # 相同 (step, state, action) → 命中缓存
    st = model.cache_stats()
    assert st["hits"] >= 1 and st["misses"] >= 1
