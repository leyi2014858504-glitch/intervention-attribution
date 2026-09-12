"""M3-2 策略先验 π(a|s):为 PUCT 提供"候选动作该往哪采"。

三种模式(可插拔 —— 消融必需):
- `uniform`  :高斯中心在动作空间正中,σ 大 ⇒ 近似均匀。**默认**(不注入任何知识)。
- `heuristic`:物理启发式。用当前 `NOx/VOC` 与**本配置自身基准**(构造时测定)之比
  判断偏向哪一侧限制区,据此偏置"该减 NOx 还是该减 VOC"。
- `learned`  :从搜索数据**自监督**学习(αZero 式自举:用搜索自身找到的高价值动作
  作回归目标),不引入外部答案。

**纪律**:禁止用 RSM/EKMA 拟合先验 —— 那会把静态响应面答案注入规划器,
破坏"规划器 vs 无规划基线"的独立性,使归因结论不可证伪。
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from sklearn.neural_network import MLPRegressor

from .features import state_features

MODES = ("uniform", "heuristic", "learned")


class PolicyPrior:
    def __init__(self, backend, *, mode: str = "uniform", sigma: float = 0.30,
                 sigma_max: float = 0.90, hidden: Sequence[int] = (32, 32), seed: int = 0,
                 min_samples: int = 24, horizon: int = 12, ratio_scale: float = 0.7):
        if mode not in MODES:
            raise ValueError(f"未知 policy mode: {mode}")
        self.backend = backend
        self.mode = mode
        self.sigma = float(sigma)              # 有明确方向时的窄 σ(聚焦)
        self.sigma_max = float(sigma_max)      # 方向不明时的宽 σ(≈ 均匀,避免无谓收窄)
        self.hidden = tuple(hidden)
        self.seed = seed
        self.min_samples = min_samples
        self.horizon = max(int(horizon), 1)
        self.ratio_scale = ratio_scale
        self.n_action = int(getattr(backend, "action_dim", 4))
        # 基准 NOx/VOC(用于启发式的相对判据,免去魔数阈值)
        pol0 = backend.project_pollutants(backend.initial_state())
        self._r_ref = max(pol0.get("NOx", 1.0), 1e-6) / max(pol0.get("VOC", 1.0), 1e-6)
        # 学习器
        self._X: list[np.ndarray] = []
        self._A: list[np.ndarray] = []
        self._V: list[float] = []
        self._net: Optional[MLPRegressor] = None
        self._mean = self._std = None
        self.n_fits = 0

    # ---- 限制区置信度 s∈[-1,1](>0:NOx 相对富 ⇒ 偏 VOC-limited) ----
    def regime_signal(self, state) -> float:
        pol = self.backend.project_pollutants(np.asarray(state, dtype=float))
        r = max(pol.get("NOx", 0.0), 1e-6) / max(pol.get("VOC", 0.0), 1e-6)
        return float(np.tanh(np.log(r / max(self._r_ref, 1e-6)) / self.ratio_scale))

    def _sigma_eff(self, state) -> float:
        """**只在有方向信息时才收窄**:|s|→1 用 σ;|s|→0 放宽到 σ_max(≈均匀)。

        不这样做时,方向不明(如基准恰好处于限制区分界)仍会把候选挤在 0.5 附近,
        变成纯粹的探索损失 —— 实测在真 CB05 上会退化约 3 %。
        """
        if self.mode == "uniform":
            return self.sigma_max
        if self.mode == "learned":
            return self.sigma
        s = abs(self.regime_signal(state))
        return float(self.sigma_max - (self.sigma_max - self.sigma) * s)

    # ---- 先验中心 μ(s) ----
    def mean_action(self, state, step_idx: int) -> np.ndarray:
        if self.mode == "heuristic":
            return self._heuristic_mean(state)
        if self.mode == "learned":
            if self._net is not None:
                x = self._raw(state, step_idx)
                return np.clip(self._net.predict(self._standardize(x[None, :]))[0], 0.02, 0.98)
        return np.full(self.n_action, 0.5)

    def _heuristic_mean(self, state) -> np.ndarray:
        s = self.regime_signal(state)
        # 动作 = [工业, 交通, 电力, 燃料VOC]。工业/燃料VOC 冲 VOC;交通/电力 冲 NOx
        direction = np.array([+0.30, -0.30, -0.30, +0.30][: self.n_action])
        return np.clip(0.5 + s * direction, 0.02, 0.98)

    # ---- 采样与权重 ----
    def sample(self, state, step_idx: int, rng: np.random.Generator) -> np.ndarray:
        if self.mode == "uniform":
            return rng.random(self.n_action)
        mu = self.mean_action(state, step_idx)
        sig = self._sigma_eff(state)
        return np.clip(mu + sig * rng.standard_normal(self.n_action), 0.0, 1.0)

    def weight(self, state, step_idx: int, action: np.ndarray) -> float:
        """π 的未归一化权重(PUCT 用;选择时按同节点子项归一)。

        uniform 模式返回常数 ⇒ 探索项对所有候选一致(标准 PUCT + 均匀先验)。
        """
        if self.mode == "uniform":
            return 1.0
        mu = self.mean_action(state, step_idx)
        d = (np.asarray(action, dtype=float) - mu) / max(self._sigma_eff(state), 1e-6)
        return float(np.exp(-0.5 * float(d @ d)))

    # ---- 自监督学习 ----
    def _raw(self, state, step_idx: int) -> np.ndarray:
        return state_features(self.backend, state, step_idx, self.horizon)

    def _standardize(self, X: np.ndarray) -> np.ndarray:
        return (X - self._mean) / self._std

    @property
    def ready(self) -> bool:
        return self._net is not None

    def add_sample(self, state, step_idx: int, action: np.ndarray, value: float) -> None:
        self._X.append(self._raw(state, step_idx))
        self._A.append(np.asarray(action, dtype=float))
        self._V.append(float(value))

    def fit(self, beta: float = 2.0) -> bool:
        """以"回报软最大"为权重回归动作(策略改进步;标签全部来自搜索自身)。"""
        if len(self._V) < self.min_samples:
            return False
        X = np.asarray(self._X, dtype=float)
        A = np.asarray(self._A, dtype=float)
        v = np.asarray(self._V, dtype=float)
        w = np.exp(beta * (v - v.max()))          # 回报越高权重越大
        w = np.clip(w, 1e-6, None)
        self._mean = X.mean(axis=0)
        self._std = X.std(axis=0)
        self._std[self._std < 1e-9] = 1.0
        net = MLPRegressor(hidden_layer_sizes=self.hidden, activation="tanh",
                           solver="lbfgs", max_iter=400, tol=1e-6,
                           random_state=self.seed)
        net.fit(self._standardize(X), A, sample_weight=w)
        self._net = net
        self.n_fits += 1
        return True
