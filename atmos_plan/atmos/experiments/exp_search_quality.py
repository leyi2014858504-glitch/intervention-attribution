"""搜索器有效性主要实验:同**箱模型调用预算**下,MCTS vs 简单基线。

回答三个问题:
1. 这个确定性开环问题真的需要 MCTS 吗?(vs 随机/常量搜索)
2. 浅树(plan_len≈2)是搜索质量瓶颈吗?
3. 加大分支宽度 width 是否换来质量?(已去除渐进加宽,候选数固定)

用法:
    python -m atmos.experiments.exp_search_quality --backend mock --budget 200 --widths 4,8
    python -m atmos.experiments.exp_search_quality --backend cb05 --budget 40 --widths 4,8,16

产出:results/search/search_quality.json
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from ..configs import CB05Config, EnvConfig
from ..search.baseline import constant_search, random_search
from ..search.mcts import MCTS, MCTSConfig
from ..search.model import EnvModel

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "search")


def _build(backend: str):
    if backend == "mock":
        from ..chemistry.mock import MockBackend
        cfg = EnvConfig(o3_limit=50.0, horizon_steps=12, dt_hours=6.0)
        return lambda: EnvModel(MockBackend(cfg), cfg), cfg, 50.0

    from ..chemistry.cb05 import CB05Backend
    bcfg = CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2)
    cfg = EnvConfig(horizon_steps=12)
    return lambda: EnvModel(CB05Backend(bcfg), cfg), cfg, None


def _probe_limit(mk, cfg):
    """关掉惩罚跑一遍零动作,取 0.7×峰值作阈值(制造可比信号)。"""
    m = mk()
    m.cfg = EnvConfig(o3_limit=1e9, horizon_steps=cfg.horizon_steps)
    c = m.backend.initial_state()
    peak = 0.0
    for i in range(cfg.horizon_steps):
        c, _, done = m.step(c, i, np.zeros(4))
        peak = max(peak, m.backend.project_pollutants(c)["O3"])
        if done:
            break
    return peak, round(0.7 * peak, 1)


def _eval_plan(mk, root, plan, cfg) -> float:
    """独立评估开环计划;**持续补足**到 horizon(保持最后动作)。

    注意:此处**必须**用持续补足而非零动作补足 —— 零补足等于"规划器没想清楚的部分
    直接关掉控制",会系统性惩罚浅计划(MCTS 是主要受害者),使对比失真。
    """
    m = mk()
    act = list(plan)
    if not act:
        act = [np.zeros(4)]
    tail = [np.asarray(act[-1], float)] * (cfg.horizon_steps - len(act))
    return m.rollout_sequence(root, (act + tail)[:cfg.horizon_steps])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=("mock", "cb05"), default="mock")
    ap.add_argument("--budget", type=int, default=200, help="MCTS 迭代数")
    ap.add_argument("--widths", type=str, default="4,8", help="各 MCTS 臂的固定分支宽度")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--puct", action="store_true", help="MCTS 用 PUCT+启发式先验")
    args = ap.parse_args()

    mk, cfg0, fixed_limit = _build(args.backend)
    if args.backend == "cb05":
        peak, limit = _probe_limit(mk, cfg0)
        cfg0 = EnvConfig(o3_limit=limit, horizon_steps=cfg0.horizon_steps)
        print(f"[probe] 基准日峰 O3 = {peak:.1f} ppb → o3_limit = {limit}", flush=True)
    else:
        limit = fixed_limit
        print(f"[cfg] mock, o3_limit = {limit}", flush=True)

    over = {"use_puct": True, "policy_mode": "heuristic"} if args.puct else {}
    rows = []

    # --- MCTS 各分支宽度 ---
    for w in [int(s) for s in args.widths.split(",") if s.strip()]:
        m = mk(); m.cfg = cfg0
        root = m.backend.initial_state()
        mc = MCTS(m, MCTSConfig(budget=args.budget, seed=args.seed, width=w, **over))
        t0 = time.perf_counter()
        res = mc.search(root)
        plan = mc.best_path(res.root)
        rows.append({"arm": f"mcts_w{w}", "box_calls": m.box_calls,
                     "value": _eval_plan(mk, root, plan, cfg0),
                     "plan_len": len(plan), "seconds": time.perf_counter() - t0})
        print(f"[run] mcts_w{w}  calls={rows[-1]['box_calls']} "
              f"value={rows[-1]['value']:.3f} plan_len={rows[-1]['plan_len']}", flush=True)

    # --- 基线:匹配"最高 MCTS 预算"(对基线不吃亏) ---
    ref_calls = max(r["box_calls"] for r in rows)
    for name, fn in (("random_search", random_search), ("constant_search", constant_search)):
        m = mk(); m.cfg = cfg0
        root = m.backend.initial_state()
        t0 = time.perf_counter()
        r = fn(m, root, ref_calls, seed=args.seed)
        rows.append({"arm": name, "box_calls": r["box_calls"],
                     "value": _eval_plan(mk, root, r["plan"], cfg0),
                     "plan_len": cfg0.horizon_steps,
                     "seconds": time.perf_counter() - t0})
        print(f"[run] {name}  calls={rows[-1]['box_calls']} "
              f"value={rows[-1]['value']:.3f}", flush=True)

    best = max(rows, key=lambda r: r["value"])
    for r in rows:
        r["value_vs_best"] = abs(r["value"]) / max(abs(best["value"]), 1e-9)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "search_quality.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"backend": args.backend, "budget": args.budget, "o3_limit": limit,
                   "seed": args.seed, "puct": args.puct, "arms": rows,
                   "note": "同箱模型调用预算;计划质量由全新模型独立评估(零动作补足)。"},
                  f, ensure_ascii=False, indent=2)

    print("\n arm                 calls   value    |v|/best  plan_len")
    for r in rows:
        print(f" {r['arm']:20s} {r['box_calls']:5d} {r['value']:9.3f} "
              f"{r['value_vs_best']:8.2f} {r['plan_len']:8d}")
    print(f"-> {path}", flush=True)


if __name__ == "__main__":
    main()
