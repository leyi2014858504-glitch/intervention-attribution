"""M3-1:价值网络 = rollout 截断器(把 CB05 搜索从"跑不动"救到"跑得动")。

角色(规格 M3 主力项 1):叶节点评估用混合
    value = w · value_net(s) + (1 − w) · rollout(s)
`w=1` 时**完全砍掉 rollout**,单次迭代成本 ≈ 1 步;`w=0` 退化为纯 rollout。
需先 warmup 收集 (s, 真实回报) 标签,网络 ready 后才允许 w>0。

硬约束
------
- 只压"单次迭代成本",**不改变被验证对象**;oracle 侧永不借用本网络。
- 确定性:固定 seed,同数据同预测(记忆化/可复现前提)。
- 可开关:上游 `use_value_net=False` 时本类完全不参与(消融必需)。

特征
----
默认用后端**污染物投影**(NOx/VOC/O3/SO2/PM25)取 log10 + 归一化深度,
低维、数据高效、后端无关。`feature_mode="full"` 改用全物种 log10(更完整但更吃数据)。
"""
from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
from sklearn.neural_network import MLPRegressor

from .features import state_features


class ValueNet:
    """从搜索中见过的 (状态, 真实回报) 增量训练的叶节点价值估计器。"""

    def __init__(self, backend, *, feature_mode: str = "projected",
                 hidden: Sequence[int] = (32, 32), seed: int = 0,
                 min_samples: int = 32, horizon: int = 12):
        self.backend = backend
        self.feature_mode = feature_mode
        self.hidden = tuple(hidden)
        self.seed = seed
        self.min_samples = min_samples
        self.horizon = max(int(horizon), 1)
        self._X: list[np.ndarray] = []
        self._y: list[float] = []
        self._net: Optional[MLPRegressor] = None
        self._mean: Optional[np.ndarray] = None
        self._std: Optional[np.ndarray] = None
        self.n_fits = 0

    # ---- 特征 ----
    def _raw(self, state, step_idx: int) -> np.ndarray:
        if self.feature_mode == "projected":
            return state_features(self.backend, state, step_idx, self.horizon)
        if self.feature_mode == "full":
            return np.asarray(state, dtype=float)
        raise ValueError(f"未知 feature_mode: {self.feature_mode}")

    def _standardize(self, X: np.ndarray) -> np.ndarray:
        return (X - self._mean) / self._std

    # ---- 数据 ----
    @property
    def n_samples(self) -> int:
        return len(self._y)

    @property
    def ready(self) -> bool:
        return self._net is not None

    def add_sample(self, state, step_idx: int, value: float) -> None:
        self._X.append(self._raw(state, step_idx))
        self._y.append(float(value))

    def fit(self) -> bool:
        """重训;样本不足返回 False(保持未就绪,上游 w 会保持 0)。"""
        if self.n_samples < self.min_samples:
            return False
        X = np.asarray(self._X, dtype=float)
        y = np.asarray(self._y, dtype=float)
        self._mean = X.mean(axis=0)
        self._std = X.std(axis=0)
        self._std[self._std < 1e-9] = 1.0
        net = MLPRegressor(hidden_layer_sizes=self.hidden, activation="tanh",
                           solver="lbfgs", max_iter=500, tol=1e-6,
                           random_state=self.seed)
        net.fit(self._standardize(X), y)
        self._net = net
        self.n_fits += 1
        return True

    def predict(self, state, step_idx: int) -> float:
        if self._net is None:
            raise RuntimeError("ValueNet 尚未训练(ready=False)")
        x = self._raw(state, step_idx)
        return float(self._net.predict(self._standardize(x[None, :]))[0])
