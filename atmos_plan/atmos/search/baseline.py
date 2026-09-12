"""搜索基线(对照):回答"这个确定性开环问题真的需要 MCTS 吗"。

- `random_search`  :采样 N 条**随机动作序列**,各跑完整 episode,留最佳。
- `constant_search`:采样 N 个**常量动作向量**(恒定减排),各跑完整 episode,留最佳。

两者都以**箱模型调用次数**计量算力,与 MCTS 在同一预算下对比。
若简单搜索不劣于 MCTS ⇒ 搜索层的树结构不产生价值(必须在论文里如实报告)。
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from .model import EnvModel


def _budget_plan(model: EnvModel, root, n_calls: int, make_plan, seed: int) -> dict:
    """通用:在 n_calls 次箱模型调用的预算内搜索,返回最佳计划/价值/实际用量。"""
    rng = np.random.default_rng(seed)
    start = model.box_calls
    best_plan, best_val = None, -np.inf
    while model.box_calls - start < n_calls:
        plan = make_plan(rng)
        v = model.rollout_sequence(root, plan)
        if v > best_val:
            best_val, best_plan = float(v), plan
    return {"plan": best_plan, "value": best_val,
            "box_calls": model.box_calls - start,
            "n_evaluated": None}


def random_search(model: EnvModel, root: np.ndarray, n_calls: int,
                  *, seed: int = 0, horizon: Optional[int] = None) -> dict:
    H = model.horizon if horizon is None else horizon
    d = model.backend.action_dim
    r = _budget_plan(model, root, n_calls,
                     lambda rng: [rng.random(d) for _ in range(H)], seed)
    r["baseline"] = "random"
    return r


def constant_search(model: EnvModel, root: np.ndarray, n_calls: int,
                    *, seed: int = 0, horizon: Optional[int] = None) -> dict:
    H = model.horizon if horizon is None else horizon
    d = model.backend.action_dim

    def make(rng):
        a = rng.random(d)
        return [a.copy() for _ in range(H)]

    r = _budget_plan(model, root, n_calls, make, seed)
    r["baseline"] = "constant"
    return r


def constant_grid_search(model: EnvModel, root: np.ndarray,
                         *, levels=(0.0, 0.5, 1.0),
                         horizon: Optional[int] = None) -> dict:
    """枚举**最优常量减排策略**(网格)。

    参数族 = [工业, 交通(=电力), 燃料VOC] 三参数,各取 `levels` ⇒ 3³=27 条恒定轨迹。
    这是"准静态环境下的理论最优"代理:若连它都赢不了,说明时序组合没有价值;
    若它被显著超越,说明时序结构存在。
    """
    H = model.horizon if horizon is None else horizon
    start = model.box_calls
    best_plan, best_val, best_cfg = None, -np.inf, None
    for a_ind in levels:
        for a_mob in levels:          # 交通与电力同档(都冲 NOx)
            for a_fvoc in levels:
                act = np.array([a_ind, a_mob, a_mob, a_fvoc], dtype=float)
                plan = [act.copy() for _ in range(H)]
                v = model.rollout_sequence(root, plan)
                if v > best_val:
                    best_val, best_plan, best_cfg = float(v), plan, act.copy()
    return {"plan": best_plan, "value": best_val, "box_calls": model.box_calls - start,
            "baseline": "constant_grid", "best_constant_action": best_cfg}


def greedy_search(model: EnvModel, root: np.ndarray, n_candidates: int,
                  *, seed: int = 0, horizon: Optional[int] = None) -> dict:
    """**逐步贪心**(1 步前瞻 + 闭环重规划)。

    每步在当前真实状态上试 `n_candidates` 个随机动作,提交"当步回报"最大者。
    环境确定且无观测噪声 ⇒ 其闭环轨迹 == 所选动作序列的开环轨迹,可与 MCTS 公平比较。
    它是"只看眼前一步"的对照:若 MCTS ≈ 贪心,则环境的时序组合没有增值。
    """
    rng = np.random.default_rng(seed)
    d = model.backend.action_dim
    H = model.horizon if horizon is None else horizon
    start = model.box_calls
    c = np.asarray(root, dtype=float)
    plan, total = [], 0.0
    for i in range(H):
        best_a, best_r, best_c = None, -np.inf, None
        for _ in range(max(1, n_candidates)):
            a = rng.random(d)
            c2, r, _done = model.step(c, i, a)
            if r > best_r:
                best_r, best_a, best_c = r, a, c2
        plan.append(best_a)
        total += best_r
        c = best_c
    return {"plan": plan, "value": total, "box_calls": model.box_calls - start,
            "baseline": "greedy"}
