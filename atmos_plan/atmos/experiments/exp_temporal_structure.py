"""时序结构判别实验 —— 论文主张"MCTS 的价值在时序组合"的生死线。

**进程隔离要求(必读)**:`CB05Backend` 存在**跨 episode 的隐藏状态污染**(实测 ~10%,
见 `tests/test_reproducibility.py` 的 xfail 用例)。因此本脚本采用
**"一个臂 = 一个子进程"** 的调度:每个子进程只跑**一个气象 episode**,
从而绕开污染(同一 weather 下的重复计算是完全确定的)。

三臂(同箱模型调用预算):
    const_grid : 网格枚举的最优**常量**策略(准静态最优代理)
    greedy     : **逐步贪心**(1 步前瞻 + 闭环重规划)
    mcts       : 完整搜索

判决见 `verdict()`:分别回答"有没有时序结构"与"前瞻/树搜索是否增值"。

用法:
    python -m atmos.experiments.exp_temporal_structure --backend cb05 --budget 40
    python -m atmos.experiments.exp_temporal_structure --backend cb05 --budget 40 \
        --weather-sigma 0.6 --episode-seed 1
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time

import numpy as np

from ..configs import CB05Config, EnvConfig
from ..search.baseline import constant_grid_search, greedy_search
from ..search.mcts import MCTS, MCTSConfig
from ..search.model import EnvModel

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "search")
ARMS = ("mcts", "greedy", "const_grid")


# --------------------------------------------------------------------------
# 后端构造
# --------------------------------------------------------------------------
def _backend_factory(backend: str, cfg: EnvConfig, weather_sigma: float, seed: int):
    if backend == "mock":
        from ..chemistry.mock import MockBackend
        return lambda: EnvModel(MockBackend(cfg), cfg)
    from ..chemistry.cb05 import CB05Backend
    bcfg = CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2,
                      vent_tau_logsigma=weather_sigma, pbl_hmax_logsigma=0.75 * weather_sigma)
    b = CB05Backend(bcfg)
    b.set_episode(seed)
    return lambda: EnvModel(b, cfg)


def _eval_plan(model: EnvModel, root, plan, cfg) -> float:
    """持续补足(保持最后动作)到 horizon —— 零补足会系统性惩罚浅计划。

    复用**同一个 model**:同一 weather 下重复计算完全确定,且省一次重复积分。
    """
    act = list(plan) or [np.zeros(4)]
    tail = [np.asarray(act[-1], float)] * (cfg.horizon_steps - len(act))
    return model.rollout_sequence(root, (act + tail)[: cfg.horizon_steps])


# --------------------------------------------------------------------------
# 动态范围探针(所有 case 用同一 weather ⇒ 不受污染影响)
# --------------------------------------------------------------------------
def dynamic_range(backend: str, cfg: EnvConfig, weather_sigma: float) -> dict:
    mk = _backend_factory(backend, cfg, weather_sigma, 0)
    root = mk().backend.initial_state()
    cases = {"zero": np.zeros(4), "full": np.ones(4),
             "voc_only": np.array([1.0, 0.0, 0.0, 1.0]),
             "nox_only": np.array([0.0, 1.0, 1.0, 0.0])}
    out = {}
    for name, a in cases.items():
        m = mk()
        c, tot, peak = root.copy(), 0.0, 0.0
        for i in range(cfg.horizon_steps):
            c, r, _ = m.step(c, i, a)
            peak = max(peak, m.backend.project_pollutants(c)["O3"])
            tot += r
        out[name] = {"value": tot, "peak_o3": peak}
    vals = [v["value"] for v in out.values()]
    out["_span"] = {"min": min(vals), "max": max(vals),
                    "range": max(vals) - min(vals),
                    "range_frac_of_max": (max(vals) - min(vals)) / max(abs(max(vals)), 1e-9)}
    return out


# --------------------------------------------------------------------------
# 单臂(在子进程中运行)
# --------------------------------------------------------------------------
def run_arm(arm: str, backend: str, cfg: EnvConfig, weather_sigma: float,
            episode_seed: int, budget: int, seed: int, width: int = 4) -> dict:
    mk = _backend_factory(backend, cfg, weather_sigma, episode_seed)
    model = mk()
    root = model.backend.initial_state()
    t0 = time.perf_counter()

    if arm == "mcts":
        mc = MCTS(model, MCTSConfig(budget=budget, seed=seed, width=width))
        res = mc.search(root)
        plan = mc.best_path(res.root)
    elif arm == "greedy":
        n_cand = max(1, cfg.horizon_steps * 4)      # 与 MCTS 预算匹配(见驱动侧说明)
        g = greedy_search(model, root, n_cand, seed=seed, horizon=cfg.horizon_steps)
        plan = g["plan"]
    elif arm == "const_grid":
        cg = constant_grid_search(model, root, horizon=cfg.horizon_steps)
        plan = cg["plan"]
    else:
        raise ValueError(arm)

    return {"arm": arm, "box_calls": model.box_calls,
            "value": _eval_plan(model, root, plan, cfg),
            "plan_len": len(plan), "seconds": time.perf_counter() - t0,
            "weather": {"hmax": [round(x, 1) for x in model.backend._weather[0]],
                        "vent_tau": [round(x, 2) for x in model.backend._weather[1]]}
            if backend == "cb05" else None}


# --------------------------------------------------------------------------
# 判决
# --------------------------------------------------------------------------
def verdict(rows: list[dict], temporal_thresh: float = 0.05,
            planning_thresh: float = 0.02) -> dict:
    by = {r["arm"]: r["value"] for r in rows}
    c, g, m = by.get("const_grid", 0.0), by.get("greedy", 0.0), by.get("mcts", 0.0)
    g_t = (g - c) / max(abs(c), 1e-9)
    g_p = (m - g) / max(abs(g), 1e-9)
    if g_t > temporal_thresh:
        temporal = {"code": "TEMPORAL", "text": "时序结构存在:逐步决策显著优于恒定策略"}
    elif abs(g_t) <= temporal_thresh:
        temporal = {"code": "QUASI_STATIC",
                    "text": "时序结构不足:逐步决策相对恒定策略无显著增益 ⇒ 环境准静态"}
    else:
        temporal = {"code": "CONST_BEST", "text": "常量策略最优 ⇒ 环境设计仍有问题"}
    if g_p > planning_thresh:
        planning = {"code": "MCTS_ADDS", "text": "MCTS 优于逐步贪心 ⇒ 前瞻有价值"}
    else:
        planning = {"code": "MCTS_NO_GAIN",
                    "text": "MCTS 未优于逐步贪心 ⇒ 搜索层的树结构未产生价值"}
    return {"temporal": temporal, "planning": planning,
            "gain_temporal": g_t, "gain_planning": g_p,
            "const": c, "greedy": g, "mcts": m,
            "summary": f"{temporal['text']};{planning['text']}"}


# --------------------------------------------------------------------------
# 驱动:每个臂一个子进程(规避 CB05 的跨 episode 污染)
# --------------------------------------------------------------------------
def _spawn(arm: str, args, o3_limit: float, effort_lambda: float) -> dict:
    with tempfile.TemporaryDirectory() as td:
        out = os.path.join(td, f"{arm}.json")
        cmd = [sys.executable, "-m", "atmos.experiments.exp_temporal_structure",
               "--arm", arm, "--backend", args.backend, "--budget", str(args.budget),
               "--seed", str(args.seed), "--episode-seed", str(args.episode_seed),
               "--weather-sigma", str(args.weather_sigma),
               "--effort-lambda", str(effort_lambda),
               "--width", str(args.width),
               "--o3-limit", str(o3_limit), "--out", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"arm {arm} 子进程失败:\n{r.stdout[-2000:]}\n{r.stderr[-2000:]}")
        with open(out, encoding="utf-8") as f:
            return json.load(f)


def _run_case(args, effort_lambda: float, o3_limit: float) -> tuple[dict, list[dict], dict]:
    """跑一个 λ 案例:动态范围探针 + 三臂,返回 (span, rows, verdict)。"""
    cfg0 = EnvConfig(horizon_steps=12, dt_hours=6.0, reward_mode="exposure",
                     effort_lambda=effort_lambda)
    span = dynamic_range(args.backend, cfg0, args.weather_sigma)
    rows = [_spawn(arm, args, o3_limit, effort_lambda) for arm in ARMS]
    return span, rows, verdict(rows)


def _print_case(effort_lambda: float, span: dict, rows: list[dict], v: dict) -> None:
    s = span["_span"]
    print(f"\n=== λ={effort_lambda:g}  可达区间 [{s['min']:.2f}, {s['max']:.2f}] "
          f"跨度={s['range']:.2f} ({s['range_frac_of_max']:.0%}) ===", flush=True)
    for r in rows:
        print(f"  [run] {r['arm']:11s} calls={r['box_calls']:5d} value={r['value']:9.3f} "
              f"plan_len={r['plan_len']}", flush=True)
    print(f"  判决:{v['summary']}")
    print(f"    (时序增益 greedy−const={v['gain_temporal']:+.1%};"
          f" 前瞻增益 mcts−greedy={v['gain_planning']:+.1%})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=("mock", "cb05"), default="cb05")
    ap.add_argument("--budget", type=int, default=40, help="MCTS 迭代数")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--episode-seed", type=int, default=1, help="气象采样种子")
    ap.add_argument("--weather-sigma", type=float, default=0.0,
                    help="日际气象变率(0 = 固定 Knote 环境)")
    ap.add_argument("--effort-lambda", type=float, default=0.0,
                    help="减排努力成本系数 [ppb / 单位平均减排比例]")
    ap.add_argument("--width", type=int, default=16,
                    help="MCTS 固定分支宽度(已去除渐进加宽)")
    ap.add_argument("--sweep", type=str, default=None,
                    help="逗号分隔的 λ 列表(如 0,10,25,50,100);给定则做扫描,忽略 --effort-lambda")
    ap.add_argument("--o3-limit", type=float, default=None)
    ap.add_argument("--arm", choices=ARMS, default=None, help="(内部)只跑单臂")
    ap.add_argument("--out", type=str, default=None, help="(内部)单臂结果 JSON")
    args = ap.parse_args()

    lim = args.o3_limit if args.o3_limit is not None else 60.0

    # ---- worker 模式 ----
    if args.arm is not None:
        cfg = EnvConfig(horizon_steps=12, dt_hours=6.0,
                        reward_mode="exposure", o3_limit=lim,
                        effort_lambda=args.effort_lambda)
        row = run_arm(args.arm, args.backend, cfg, args.weather_sigma,
                      args.episode_seed, args.budget, args.seed, args.width)
        if args.out:
            with open(args.out, "w", encoding="utf-8") as f:
                json.dump(row, f)
        return

    # ---- driver 模式 ----
    if args.sweep:
        lambdas = [float(x) for x in args.sweep.split(",")]
    else:
        lambdas = [args.effort_lambda]
    print(f"[cfg] backend={args.backend} budget={args.budget} "
          f"weather_sigma={args.weather_sigma} episode_seed={args.episode_seed} "
          f"lambda={lambdas}", flush=True)

    cases = []
    for lam in lambdas:
        span, rows, v = _run_case(args, lam, lim)
        _print_case(lam, span, rows, v)
        cases.append({"effort_lambda": lam, "dynamic_range": span, "arms": rows, "verdict": v})

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "temporal_structure.json")
    payload = {"backend": args.backend, "budget": args.budget, "seed": args.seed,
               "episode_seed": args.episode_seed, "weather_sigma": args.weather_sigma,
               "reward_mode": "exposure", "o3_limit": lim, "cases": cases,
               "note": "每臂一个子进程(规避 CB05 跨 episode 污染);计划持续补足评估。"}
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"-> {path}", flush=True)


if __name__ == "__main__":
    main()
