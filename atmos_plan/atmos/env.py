"""MDP 环境层:把化学后端封装成分段决策过程(物种无关)。

环境**不知道**具体物种:状态 = 后端完整状态 + [h_pbl, temp];
动作 = 4 维部门减排率 ∈[0,1](语义由后端翻译);奖励由后端污染物投影算。
故 Mock → CB05 → SAPRC 可换,env 不改。
"""
from __future__ import annotations

import numpy as np

from .chemistry.base import ChemistryBackend
from .configs import EnvConfig
from .physics import pbl_height, temperature   # 共享物理单一真源(向后兼容再导出)

__all__ = ["AtmosEnv", "pbl_height", "temperature"]


class AtmosEnv:
    """72h 分部门减排 MDP。状态维度由后端物种数决定。"""

    action_dim = 4                      # [工业, 交通, 电力, VOC控制]

    def __init__(self, backend: ChemistryBackend, config: EnvConfig | None = None):
        self.backend = backend
        self.cfg = config or EnvConfig()
        self.dt = self.cfg.dt_hours
        self.horizon = self.cfg.horizon_steps
        self.n_state = backend.n_species + 2      # + h_pbl + temp
        self._conc: np.ndarray | None = None
        self._step_idx = 0

    def _h_pbl_fn(self):
        # 经后端取 PBL —— 保证"观测里的 PBL"与"化学里的 PBL"一致(单一真源)
        return lambda t: self.backend.pbl_height_at(t)

    def reset(self, *, seed: int | None = None, options: dict | None = None):
        rng = np.random.default_rng(seed)
        base = self.backend.initial_state()
        jit = 1.0 + self.cfg.init_jitter * rng.standard_normal(base.shape)
        self._conc = np.clip(base * jit, 0.0, None)
        self._step_idx = 0
        return self._obs(0.0)

    def _obs(self, t: float) -> np.ndarray:
        return np.concatenate([self._conc, [self.backend.pbl_height_at(t),
                                            self.backend.temperature_at(t)]])

    def step(self, action: np.ndarray):
        assert self._conc is not None, "先调用 reset()"
        t0 = self._step_idx * self.dt
        act = np.clip(np.asarray(action, dtype=float), 0.0, 1.0)
        res = self.backend.step(
            self._conc, t0, self.dt,
            action=act,
            met={"h_pbl_fn": self._h_pbl_fn(), "temp": temperature(t0)},
        )
        self._conc = res.conc_out
        self._step_idx += 1
        reward = self._reward(res.conc_traj, act)
        done = self._step_idx >= self.horizon
        return self._obs(self._step_idx * self.dt), reward, done, {"t_end": self._step_idx * self.dt}

    def _reward(self, conc_traj: np.ndarray, action: np.ndarray | None = None) -> float:
        """步内回报:空气质量代价 + **减排努力成本**。

        默认单目标 O3(use_pm_penalty=False);PM 项需含气溶胶机制才启用。
        努力成本 `λ·mean(action)` 使行动有代价 ⇒ 必须择时(λ=0 时退化为免费减排)。
        """
        mean = conc_traj.mean(axis=0)
        pol = self.backend.project_pollutants(mean)
        if self.cfg.reward_mode == "exposure":
            r = -float(pol.get("O3", 0.0))          # 连续暴露:无阈值、无饱和
        else:
            exceed_o3 = max(0.0, pol.get("O3", 0.0) - self.cfg.o3_limit)
            r = -self.cfg.lam_o3 * exceed_o3
            if self.cfg.use_pm_penalty:
                exceed_pm = max(0.0, pol.get("PM25", 0.0) - self.cfg.pm25_limit)
                r -= exceed_pm
        if self.cfg.effort_lambda and action is not None:
            a = np.clip(np.asarray(action, float), 0.0, 1.0)
            r -= self.cfg.effort_lambda * float(a.mean())
        return r
