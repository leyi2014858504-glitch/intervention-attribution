"""M3 闸门验证(真 CB05):比较 UCT / PUCT / 策略先验 / value-net 截断。

用法:
    python -m atmos.experiments.exp_gate_cb05 --smoke
    python -m atmos.experiments.exp_gate_cb05 --budget 40 --arms uct,puct_heuristic
    python -m atmos.experiments.exp_gate_cb05 --budget 40 --arms uct,puct_learned,puct_heur_vnet

指标:**箱模型调用数**(算力)与**独立评估的计划质量**(越接近 0 越好),并报计划深度。

- 奖励阈值:CB05 基准日峰 O3 通常 90–100 ppb,默认 160 会令奖励恒 0(无信号),
  故先探一次基准峰值,取 `o3_limit = 0.7 × peak` 制造可比信号(所有臂同一阈值)。
- 计划评估:贪心路径若不足 horizon,统一用**零动作补足**(各臂同规则,公平)。
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from ..chemistry.cb05 import CB05Backend
from ..configs import CB05Config, EnvConfig
from ..search.mcts import MCTS, MCTSConfig
from ..search.model import EnvModel

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "gate")

#: 各臂的 MCTSConfig 覆盖项
ARMS: dict[str, dict] = {
    "uct":          {},
    "puct_uniform": {"use_puct": True, "policy_mode": "uniform"},
    "puct_heuristic": {"use_puct": True, "policy_mode": "heuristic"},
    "puct_learned": {"use_puct": True, "policy_mode": "learned"},
    "vnet":         {"use_value_net": True, "value_w": 1.0, "warmup_iters": 20,
                     "retrain_every": 20, "label_fraction": 0.2, "vnet_min_samples": 16},
    "puct_heur_vnet": {"use_puct": True, "policy_mode": "heuristic",
                       "use_value_net": True, "value_w": 1.0, "warmup_iters": 20,
                       "retrain_every": 20, "label_fraction": 0.2, "vnet_min_samples": 16},
}


def _backend_cfg() -> CB05Config:
    # 速度优先(仅用于搜索动力学验证,不作策略结论)
    return CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2)


def probe_peak_o3(horizon: int) -> tuple[float, float]:
    cfg = EnvConfig(o3_limit=1e9, horizon_steps=horizon)
    model = EnvModel(CB05Backend(_backend_cfg()), cfg)
    c = model.backend.initial_state()
    peak, t0 = 0.0, time.perf_counter()
    for i in range(horizon):
        c, _, done = model.step(c, i, np.zeros(4))
        peak = max(peak, model.backend.project_pollutants(c)["O3"])
        if done:
            break
    return peak, time.perf_counter() - t0


def eval_plan(root: np.ndarray, plan, horizon: int, cfg: EnvConfig) -> float:
    model = EnvModel(CB05Backend(_backend_cfg()), cfg)
    acts = list(plan) + [np.zeros(4) for _ in range(horizon - len(plan))]
    return model.rollout_sequence(root, acts[:horizon])


def run_arm(name: str, *, o3_limit: float, horizon: int, budget: int,
            seed: int = 1) -> dict:
    over = ARMS[name]
    cfg = EnvConfig(o3_limit=o3_limit, horizon_steps=horizon)
    model = EnvModel(CB05Backend(_backend_cfg()), cfg)
    root = model.backend.initial_state()
    m = MCTS(model, MCTSConfig(budget=budget, seed=seed, **over))
    t0 = time.perf_counter()
    res = m.search(root)
    elapsed = time.perf_counter() - t0
    plan = m.best_path(res.root)
    val = eval_plan(root, plan, horizon, cfg)
    return {"arm": name, "value": val, "box_calls": model.box_calls,
            "seconds": elapsed, "plan_len": len(plan),
            "vnet_ready": res.vnet_ready, "vnet_samples": res.vnet_samples,
            "policy_ready": res.policy_ready, "policy_fits": res.policy_fits}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--budget", type=int, default=40)
    ap.add_argument("--horizon", type=int, default=12)
    ap.add_argument("--arms", type=str, default="uct,puct_heuristic")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    names = [s.strip() for s in args.arms.split(",") if s.strip()]
    bad = [n for n in names if n not in ARMS]
    if bad:
        raise SystemExit(f"未知臂 {bad};可选:{list(ARMS)}")
    print(f"[cfg] budget={args.budget} horizon={args.horizon} arms={names}", flush=True)

    peak, t_probe = probe_peak_o3(args.horizon)
    thr = round(0.7 * peak, 1)
    print(f"[probe] 基准日峰 O3 = {peak:.1f} ppb ({t_probe:.1f}s) → o3_limit = {thr}", flush=True)
    if args.smoke:
        print("[smoke] 通过:CB05 + EnvModel 可跑", flush=True)
        return

    rows = []
    for nm in names:
        print(f"[run] arm={nm} ...", flush=True)
        r = run_arm(nm, o3_limit=thr, horizon=args.horizon, budget=args.budget,
                    seed=args.seed)
        rows.append(r)
        print(f"  -> value={r['value']:9.3f} calls={r['box_calls']:5d} "
              f"{r['seconds']:6.1f}s plan_len={r['plan_len']} "
              f"pol={r['policy_ready']}/{r['policy_fits']} "
              f"vnet={r['vnet_ready']}/{r['vnet_samples']}", flush=True)

    base = rows[0]
    for r in rows:
        r["calls_vs_first"] = r["box_calls"] / max(base["box_calls"], 1)
        r["quality_ratio_abs"] = abs(r["value"]) / max(abs(base["value"]), 1e-9)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "gate_cb05.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"o3_limit": thr, "baseline_peak_o3": peak, "budget": args.budget,
                   "horizon": args.horizon, "seed": args.seed, "arms": rows,
                   "note": "真 CB05;阈值取基准日峰 0.7 以制造信号;计划零动作补足到 horizon。"},
                  f, ensure_ascii=False, indent=2)
    print(f"[done] -> {path}", flush=True)


if __name__ == "__main__":
    main()
