"""搜索层共享状态特征(value-net 与策略先验共用)。

设计要点:低维、后端无关、含**物理派生比**——
`NOx/VOC`(限制区指示)与 `O3/NOx`(氧化态)能让小样本下的网络更快抓住结构,
这是"把 value-net 喂饱"的一部分(仅 5 个绝对浓度 + 深度时样本效率低)。
"""
from __future__ import annotations

import numpy as np

_POLL_KEYS = ("NOx", "VOC", "O3", "SO2", "PM25")
N_FEATURES = len(_POLL_KEYS) + 3          # 5 污染物 + 2 比值 + 深度


def state_features(backend, state, step_idx: int, horizon: int) -> np.ndarray:
    """状态 → 特征向量(log10 尺度 + 物理比值 + 归一化深度)。"""
    pol = backend.project_pollutants(np.asarray(state, dtype=float))
    vals = np.array([pol.get(k, 0.0) for k in _POLL_KEYS], dtype=float)
    x = np.log10(np.clip(vals, 1e-6, None))
    nox, voc, o3 = vals[0], vals[1], vals[2]
    ratios = np.array([np.log10(max(nox, 1e-6) / max(voc, 1e-6)),
                       np.log10(max(o3, 1e-6) / max(nox, 1e-6))])
    return np.concatenate([x, ratios, [step_idx / max(horizon, 1)]])
