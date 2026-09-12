"""函数式环境模型:给 MCTS/oracle 用,可直接从任意状态分叉。

与 gym 风格 `AtmosEnv` 的区别:`AtmosEnv` 有可变内部状态(适合顺序交互),
而本类**无内部状态**——`step(conc, step_idx, action)` 是纯函数,天然支持树搜索
的分支/回溯。底层复用 `ChemistryBackend.step`(契约要求无副作用+确定性)。

记忆化:键 = (step_idx, conc 位哈希, action 位哈希)。同参数二次调用直接命中缓存,
不重复积分。跨搜索/oracle 复用同一实例即可省掉大量重复 rollout。
"""
from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from ..chemistry.base import ChemistryBackend
from ..configs import EnvConfig


class EnvModel:
    def __init__(self, backend: ChemistryBackend, config: EnvConfig | None = None):
        self.backend = backend
        self.cfg = config or EnvConfig()
        self.dt = self.cfg.dt_hours
        self.horizon = self.cfg.horizon_steps
        self._cache: dict[tuple, tuple] = {}
        self.hits = 0
        self.misses = 0

    # ---- 纯函数单步 ----
    def _reward(self, conc_traj: np.ndarray, action: np.ndarray | None = None) -> float:
        """步内回报:空气质量代价 + **减排努力成本**。

        `r = −(暴露 + λ·mean(action))`。λ=0 时退化为纯空气质量项(免费减排);
        λ>0 时行动有代价 ⇒ 必须**择时**,这正是制造时序结构的开关。
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

    def _met(self, t0: float) -> dict:
        # 经后端取环境量(单一真源);CB05 自持内部 PBL,met 仅为 Mock 一类沿用
        b = self.backend
        return {"h_pbl_fn": lambda t: b.pbl_height_at(t), "temp": b.temperature_at(t0)}

    def step(self, conc: np.ndarray, step_idx: int, action: np.ndarray):
        """返回 (next_conc, reward, done)。确定性 + 记忆化。"""
        conc = np.asarray(conc, float)
        action = np.clip(np.asarray(action, float), 0.0, 1.0)
        key = (int(step_idx), conc.tobytes(), action.tobytes())
        cached = self._cache.get(key)
        if cached is not None:
            self.hits += 1
            return cached
        self.misses += 1
        t0 = step_idx * self.dt
        res = self.backend.step(conc, t0, self.dt, action=action, met=self._met(t0))
        out = (res.conc_out.copy(), self._reward(res.conc_traj, action),
               step_idx + 1 >= self.horizon)
        self._cache[key] = out
        return out

    # ---- rollout ----
    def rollout(self, conc: np.ndarray, step_idx: int,
                policy: Optional[Callable[[np.ndarray, int], np.ndarray]] = None,
                rng: Optional[np.random.Generator] = None,
                max_steps: Optional[int] = None) -> float:
        """从 (conc, step_idx) 用 policy(默认随机均匀)推演到终局,返回累积回报。"""
        rng = rng or np.random.default_rng()
        policy = policy or (lambda c, i: rng.random(self.backend.action_dim))
        c, r_total, i = np.asarray(conc, float), 0.0, step_idx
        limit = self.horizon if max_steps is None else min(self.horizon, step_idx + max_steps)
        while i < limit:
            c, r, done = self.step(c, i, policy(c, i))
            r_total += r
            i += 1
            if done:
                break
        return r_total

    def rollout_sequence(self, conc: np.ndarray, actions: list[np.ndarray]) -> float:
        """按给定动作序列(开环)从 conc 推演,返回累积回报。用于评估规划质量。"""
        c, total = np.asarray(conc, float), 0.0
        for i, a in enumerate(actions):
            if i >= self.horizon:
                break
            c, r, done = self.step(c, i, a)
            total += r
            if done:
                break
        return total

    @property
    def box_calls(self) -> int:
        """已发生的箱模型积分次数(= 缓存未命中数),搜索算力成本的主指标。"""
        return self.misses

    def clear_cache(self) -> None:
        """清空 rollout 记忆化缓存。

        注意:缓存键 = (step_idx, 状态, 动作),**不含后端情景**。若更换了
        `backend.set_scenario(...)` 或后端配置,必须调用本方法,否则会命中旧结果。
        """
        self._cache.clear()
        self.hits = 0
        self.misses = 0

    def cache_stats(self) -> dict:
        total = self.hits + self.misses
        return {"hits": self.hits, "misses": self.misses,
                "hit_rate": (self.hits / total) if total else 0.0,
                "size": len(self._cache)}
