"""MDP 环境层单测:接口可跑 + 随机策略 + 确定性(记忆化前置)。"""
import numpy as np

from atmos.configs import EnvConfig
from atmos.chemistry.mock import MockBackend
from atmos.env import AtmosEnv
from atmos.search.model import EnvModel


def _rollout(env, actions):
    env.reset(seed=0)
    tot, steps = 0.0, 0
    done = False
    for a in actions:
        obs, r, done, _ = env.step(a)
        tot += r
        steps += 1
    return obs, tot, done, steps


def test_random_policy_episode():
    env = AtmosEnv(MockBackend(EnvConfig()), EnvConfig())
    obs0 = env.reset(seed=1)
    assert obs0.shape == (env.n_state,)
    rng = np.random.default_rng(7)
    acts = rng.random((env.horizon, env.action_dim))
    obs, tot, done, steps = _rollout(env, acts)
    assert steps == env.horizon and done
    assert np.isfinite(tot)
    assert obs.shape == (env.n_state,)
    # 状态末两维为 h_pbl / temp,量级合理
    assert 100 <= obs[-2] <= 5000 and 250 <= obs[-1] <= 340


def test_action_clipped_to_unit():
    env = AtmosEnv(MockBackend(EnvConfig()), EnvConfig())
    env.reset(seed=0)
    # 越界动作应被裁剪进 [0,1],不抛错
    obs, r, done, _ = env.step(np.array([2.0, -1.0, 0.5, 5.0]))
    assert np.isfinite(r)


def test_determinism_for_memoization():
    """相同初始+相同动作序列 → 相同总回报(rollout 记忆化的正确性前提)。"""
    acts = np.full((12, 4), 0.3)
    env1 = AtmosEnv(MockBackend(EnvConfig()), EnvConfig())
    env2 = AtmosEnv(MockBackend(EnvConfig()), EnvConfig())
    _, t1, _, _ = _rollout(env1, acts)
    _, t2, _, _ = _rollout(env2, acts)
    assert abs(t1 - t2) < 1e-6


def test_step_handoff_carries_state():
    """跨步接力:第 0 步的减排必须改变第 0 步**末态**(即下一步的起始浓度)。

    若此断言失败,说明环境构造上解耦(动作只改排放、不改起始状态),
    那么任何规划都不可能有用 —— 后续环境改造都是在给坏环境调参。
    """
    model = EnvModel(MockBackend(EnvConfig()), EnvConfig())
    c0 = model.backend.initial_state()
    c_base, _, _ = model.step(c0, 0, np.zeros(4))
    c_pert, _, _ = model.step(c0, 0, np.ones(4))
    rel = np.linalg.norm(c_pert - c_base) / (np.linalg.norm(c_base) or 1.0)
    assert rel > 1e-3, f"接力断裂:动作未影响步末状态(rel={rel:.2e})"


def test_zero_action_vs_full_reduction_differ():
    """全零动作与全减排的累积回报应不同(动作确实影响化学/奖励)。

    用 Mock 的 PM 通道验证动作→奖励耦合(Mock 是双污染物测试替身;
    CB05 科学配置默认关 PM,见 configs.EnvConfig.use_pm_penalty)。
    """
    cfg = EnvConfig(use_pm_penalty=True)
    env = AtmosEnv(MockBackend(cfg), cfg)
    _, t_zero, _, _ = _rollout(env, np.zeros((12, 4)))
    _, t_red, _, _ = _rollout(env, np.full((12, 4), 0.8))
    assert t_zero != t_red
