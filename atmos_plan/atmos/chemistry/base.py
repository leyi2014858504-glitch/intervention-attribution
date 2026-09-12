"""化学后端抽象契约(通用化:状态 = 后端自身物种空间,不再写死 5 维)。

上层(MDP / 搜索 / 归因 / oracle)依赖:
- `species`        :后端状态布局(如 CB05 的 79 物种);决定 MDP 状态维度。
- `initial_state`  :背景初始完整状态(供 reset 采样)。
- `step`           :把化学从 t0 推演 dt;动作(4 维部门减排)由后端翻译进自身机制。
- `project_pollutants` :完整状态 → 政策污染物 {NOx,VOC,O3,SO2,PM25}(供奖励/AQI)。

`step` 必须**无副作用 + 确定性**:同一 (conc_in,t0,dt,action,met) → 同一结果。
这是 rollout 记忆化与精确 oracle 的前提(记忆化键 = 完整状态哈希 + 动作序列)。
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class StepResult:
    """单个控制步(如 6h)的化学推演结果。conc_* 为**完整后端状态**(len=species)。"""
    conc_out: np.ndarray      # 步末完整状态
    conc_traj: np.ndarray     # 步内轨迹 (n_times, n_species),含步初
    times: np.ndarray         # 对应时间(小时)


class ChemistryBackend(ABC):
    """0-D 化学箱模型后端契约。"""

    #: 状态变量名(顺序即状态向量布局);各后端自定义
    species: Sequence[str]

    @property
    def n_species(self) -> int:
        return len(self.species)

    #: 动作维度:固定为 [工业, 交通, 电力, 燃料VOC控制] 4 维部门减排率 ∈[0,1]
    action_dim: int = 4

    #: 边界层日循环参数(后端是**环境的唯一真源**;环境层必须经 pbl_height_at 取值,
    #: 不得自行计算,否则观测里的 PBL 与化学里的 PBL 会不一致)
    h_min: float = 250.0
    h_max: float = 1800.0
    hod0: float = 9.0          # t=0 时的当地时(小时)

    def pbl_height_at(self, t_hours: float) -> float:
        """(单步开始时)该时刻的边界层高度 [m]。"""
        from ..physics import pbl_height
        return pbl_height(t_hours, self.h_min, self.h_max, self.hod0)

    def temperature_at(self, t_hours: float) -> float:
        from ..physics import temperature
        return temperature(t_hours, hod0=self.hod0)

    @abstractmethod
    def initial_state(self) -> np.ndarray:
        """背景完整初始状态(未加噪声),供环境 reset 采样。"""
        raise NotImplementedError

    @abstractmethod
    def step(self, conc_in, t0, dt, *, action, met) -> StepResult:
        """把化学从 t0 推演 dt 小时。

        Parameters
        ----------
        conc_in : (n_species,) 步初完整状态
        t0, dt  : 起始时间 / 步长(小时)
        action  : (action_dim,) 部门减排率 ∈[0,1],由后端翻译为其物种的排放注入
        met     : 至少含 {"h_pbl_fn": callable(t)->m, "temp": K}
        """
        raise NotImplementedError

    @abstractmethod
    def project_pollutants(self, conc) -> Mapping[str, float]:
        """完整状态 → {NOx, VOC, O3, SO2, PM25}(ppb;PM25 µg/m³)。缺项返回 0。"""
        raise NotImplementedError
