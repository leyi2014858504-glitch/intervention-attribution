"""M3-2 策略先验 + PUCT 测试(MockBackend,毫秒级)。

核心断言:同预算下 PUCT + 策略先验的计划质量**优于** UCT(探索项衰减更快 + 先验聚焦)。
"""
import numpy as np
import pytest

from atmos.configs import EnvConfig
from atmos.chemistry.mock import MockBackend, MOCK_SPECIES
from atmos.search.model import EnvModel
from atmos.search.mcts import MCTS, MCTSConfig
from atmos.search.policy import PolicyPrior

H = 12


def _model():
    cfg = EnvConfig(o3_limit=50.0, horizon_steps=H, dt_hours=6.0)
    return EnvModel(MockBackend(cfg), cfg)


def _backend():
    return MockBackend(EnvConfig(o3_limit=50.0, horizon_steps=H))


def _run(**over):
    model = _model()
    root = model.backend.initial_state()
    m = MCTS(model, MCTSConfig(budget=200, seed=1, **over))
    res = m.search(root)
    plan = m.best_path(res.root)
    val = _model().rollout_sequence(root, plan)
    return val, model.box_calls, len(plan), res


# ---- PolicyPrior 本体 ----

def test_uniform_policy_samples_in_unit_range():
    p = PolicyPrior(_backend(), mode="uniform", horizon=H)
    rng = np.random.default_rng(0)
    s = _backend().initial_state()
    for _ in range(20):
        a = p.sample(s, 0, rng)
        assert a.shape == (4,) and np.all((a >= 0) & (a <= 1))
    assert p.weight(s, 0, np.full(4, 0.9)) == pytest.approx(1.0)   # 均匀 ⇒ 常数权重


def test_heuristic_prior_biases_by_regime():
    """NOx 相对富(VOC-limited)→ 偏向 VOC 减排(工业/燃料VOC 维);反之偏向 NOx 减排。"""
    p = PolicyPrior(_backend(), mode="heuristic", horizon=H)
    idx = {s: i for i, s in enumerate(MOCK_SPECIES)}
    nox_rich = np.zeros(len(MOCK_SPECIES)); nox_rich[idx["NOx"]] = 300; nox_rich[idx["VOC"]] = 30
    voc_rich = np.zeros(len(MOCK_SPECIES)); voc_rich[idx["NOx"]] = 20; voc_rich[idx["VOC"]] = 400
    # 动作维 = [工业(VOC), 交通(NOx), 电力(NOx), 燃料VOC]
    m_noxrich = p.mean_action(nox_rich, 0)
    m_voc_RICH = p.mean_action(voc_rich, 0)
    assert m_noxrich[0] > m_noxrich[1], f"VOC-limited 应偏 VOC 减排: {m_noxrich}"
    assert m_voc_RICH[1] > m_voc_RICH[0], f"NOx-limited 应偏 NOx 减排: {m_voc_RICH}"


def test_learned_policy_fits_and_predicts():
    p = PolicyPrior(_backend(), mode="learned", horizon=H, min_samples=12)
    assert not p.ready
    model = _model()
    rng = np.random.default_rng(0)
    c = model.backend.initial_state()
    for _ in range(30):
        a = rng.random(4)
        c2, r, _ = model.step(c, 0, a)
        v = model.rollout(c2, 1, rng=rng)
        p.add_sample(c, 0, a, v)
    assert p.fit() and p.ready
    mu = p.mean_action(c, 0)
    assert mu.shape == (4,) and np.all((mu >= 0.02) & (mu <= 0.98))


# ---- PUCT 集成 ----

def test_puct_runs_and_records_priors():
    val, calls, plen, res = _run(use_puct=True, policy_mode="heuristic")
    assert np.isfinite(val)
    assert res.root.children, "根节点应有子"
    assert len(res.root.priors) == len(res.root.children), "每个候选都应带 π 权重"


def test_puct_with_prior_beats_uct():
    """同预算(200)下,PUCT+启发式先验的计划质量应优于 UCT。"""
    val_uct, _, _, _ = _run()
    val_puct, calls_puct, _, _ = _run(use_puct=True, policy_mode="heuristic")
    assert val_puct > val_uct, f"PUCT 未改善: UCT={val_uct:.4f} PUCT={val_puct:.4f}"


def test_puct_determinism_same_seed():
    v1, c1, _, _ = _run(use_puct=True, policy_mode="heuristic")
    v2, c2, _, _ = _run(use_puct=True, policy_mode="heuristic")
    assert v1 == pytest.approx(v2) and c1 == c2


def test_disabled_puct_matches_uct_path():
    """use_puct=False 时不得实例化策略(消融可对照)。"""
    _, _, _, res = _run()
    assert res.policy_fits == 0 and not res.policy_ready and res.policy_samples == 0
