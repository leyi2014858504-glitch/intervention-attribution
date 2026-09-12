"""跨步状态接力 + 记忆强度诊断(最便宜、最优先的嫌疑排查)。

`greedy − const = +2.6%` 有两种解释,化学上像、命运完全不同:
  (A) **物理短记忆**:O3 寿命短 ⇒ 各控制步近似独立 ⇒ 逐日贪心天然最优;
  (B) **实现未接力**:每步从基线重新锚定(动作只改排放、不改下一步起始浓度)
      ⇒ 环境在构造上解耦,任何规划都不可能有用。

两个探针:
1. **接力探针**:第 0 步给最大减排、其余零,比较第 0 步**末态**与基线末态。
   若两者几乎相同 ⇒ 接力断了(查代码);若显著不同 ⇒ 接力正常。
2. **记忆衰减**:第 0–3 步最大减排、其后零;比较第 4–11 步 O3 与基线的差。
   Δ≈0 ⇒ 记忆太弱(物理);Δ 显著 ⇒ 记忆存在,问题在别处。

用法:python -m atmos.experiments.exp_memory_check --backend cb05
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

from ..configs import CB05Config, EnvConfig
from ..search.model import EnvModel

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "search")


def _factory(backend: str, cfg: EnvConfig):
    if backend == "mock":
        from ..chemistry.mock import MockBackend
        return lambda: EnvModel(MockBackend(cfg), cfg)
    from ..chemistry.cb05 import CB05Backend
    return lambda: EnvModel(CB05Backend(
        CB05Config(rtol=1e-3, atol=1e6, substeps_per_control=2)), cfg)


def episode(mk, cfg, actions):
    """跑一条开环序列,返回逐步末态矩阵与逐步 O3。"""
    m = mk()
    c = m.backend.initial_state()
    states = [c.copy()]
    o3 = []
    for i, a in enumerate(actions[: cfg.horizon_steps]):
        c, _r, _done = m.step(c, i, np.asarray(a, float))
        states.append(c.copy())
        o3.append(m.backend.project_pollutants(c)["O3"])
    return np.asarray(states), np.asarray(o3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=("mock", "cb05"), default="cb05")
    args = ap.parse_args()

    cfg = EnvConfig(horizon_steps=12, reward_mode="exposure")
    mk = _factory(args.backend, cfg)
    H = cfg.horizon_steps
    Z, FULL = np.zeros(4), np.ones(4)

    # 基线:全程零动作
    st_base, o3_base = episode(mk, cfg, [Z] * H)

    # --- 探针 1:接力(仅第 0 步最大减排) ---
    st_pr, o3_pr = episode(mk, cfg, [FULL] + [Z] * (H - 1))
    denom = np.linalg.norm(st_base[1]) or 1.0
    rel_state = float(np.linalg.norm(st_pr[1] - st_base[1]) / denom)
    print("=== 探针 1:状态接力(仅第 0 步最大减排)===")
    print(f"  第 0 步末态相对差 = {rel_state:.4f}   (>1e-3 判定接力正常)")
    print(f"  第 0 步 O3: 基线={o3_base[0]:.2f}  扰动={o3_pr[0]:.2f}  "
          f"Δ={o3_pr[0]-o3_base[0]:+.2f} ppb", flush=True)

    # --- 探针 2:记忆衰减(第 0–3 步最大减排,其后零) ---
    acts = [FULL] * 4 + [Z] * (H - 4)
    _, o3_mem = episode(mk, cfg, acts)
    print("\n=== 探针 2:记忆衰减(第 0–3 步减排,第 4 步起零动作)===")
    print(f"  {'step':>4} {'hod':>5} {'base O3':>8} {'abated O3':>10} {'Δ':>8}")
    deltas = []
    for i in range(H):
        d = o3_mem[i] - o3_base[i]
        deltas.append(d)
        tag = "  <- 减排期" if i < 4 else ""
        print(f"  {i:4d} {(i*6+9)%24:5.0f} {o3_base[i]:8.2f} {o3_mem[i]:10.2f} "
              f"{d:+8.2f}{tag}", flush=True)
    post = np.array(deltas[4:])
    print(f"\n  减排期后(|Δ| 均值) = {np.abs(post).mean():.2f} ppb;"
          f" 最大 |Δ| = {np.abs(post).max():.2f} ppb", flush=True)

    # 衰减半程:扰动影响降到一半所需步数
    ref = abs(deltas[3]) or 1e-9
    half = next((i for i in range(4, H) if abs(deltas[i]) < 0.5 * ref), None)
    print(f"  第 3 步 |Δ|={abs(deltas[3]):.2f}; 衰减到一半的步数 = "
          f"{half if half is not None else '未衰减'}", flush=True)

    handoff_ok = rel_state > 1e-3
    memory_alive = np.abs(post).max() > 1.0
    if not handoff_ok:
        concl = "接力断裂 ⇒ 先去修 rollout 的状态交接,不要做环境改造"
    elif memory_alive:
        concl = "接力正常且记忆存在但偏弱 ⇒ 属物理问题,进入非平稳强迫改造"
    else:
        concl = "接力正常但记忆极弱(≈0)⇒ 物理短记忆,各步近乎独立"
    print(f"\n结论:{concl}", flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "memory_check.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"backend": args.backend, "handoff_rel_diff": rel_state,
                   "o3_base": o3_base.tolist(), "o3_abated": o3_mem.tolist(),
                   "delta": deltas, "post_abs_mean": float(np.abs(post).mean()),
                   "conclusion": concl}, f, ensure_ascii=False, indent=2)
    print(f"-> {path}", flush=True)


if __name__ == "__main__":
    main()
