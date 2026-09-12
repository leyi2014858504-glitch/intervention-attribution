"""共享物理环境(PBL 日循环 / 温度 / 夹卷稀释)。

**单一真源**:环境层(`AtmosEnv`)与化学后端(`CB05Backend`)都从这里取 PBL,
避免"观测里的 PBL"和"化学里的 PBL"不一致(曾发生过:后端完全没用 PBL)。

对齐 Knote et al. (2015) AQMEII-2 理想化箱模型的范式:环境参数外生给定,
不从目标输出反推。
"""
from __future__ import annotations

import math


def pbl_height(t_hours: float, h_min: float = 250.0, h_max: float = 1800.0,
               hod0: float = 9.0) -> float:
    """边界层高度日循环(米)。`hod0` 为 t=0 时刻的当地时(小时)。"""
    hod = (t_hours + hod0) % 24.0
    s = max(0.0, math.sin(math.pi * (hod - 6.0) / 12.0))
    return h_min + (h_max - h_min) * math.sqrt(s)


def temperature(t_hours: float, t_ref: float = 298.15, amp: float = 6.0,
                hod0: float = 9.0) -> float:
    """近地面温度日循环(K)。"""
    hod = (t_hours + hod0) % 24.0
    return t_ref + amp * max(0.0, math.sin(math.pi * (hod - 6.0) / 12.0))


def entrainment_dilution(t_hours: float, h_min: float = 250.0,
                         h_max: float = 1800.0, hod0: float = 9.0) -> float:
    """边界层抬升夹卷自由对流层空气导致的稀释速率 [1/s]。

    k = (1/h)·(dh/dt),且**仅在边界层增长时为正**(收缩不产生稀释)。
    这是 Knote et al. (2015) 一类的理想化箱模型处理"通风/稀释"的标准做法,
    取代此前手调的常数 τ。
    """
    hod = (t_hours + hod0) % 24.0
    s = math.sin(math.pi * (hod - 6.0) / 12.0)
    if s <= 1e-6:
        return 0.0
    h = h_min + (h_max - h_min) * math.sqrt(s)
    dhdt = ((h_max - h_min)
            * (math.cos(math.pi * (hod - 6.0) / 12.0) * math.pi / 12.0)
            / (2.0 * math.sqrt(s)))          # m/h
    return max(0.0, dhdt / h) / 3600.0        # 1/s
