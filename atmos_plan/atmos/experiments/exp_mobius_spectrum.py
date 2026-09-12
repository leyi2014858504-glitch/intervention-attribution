"""M4 前哨:Möbius(FS / ANOVA)交互谱枚举。

**值函数** v(S) = 干预后的**暴露均值**(= 每步 reward 的均值,越大越好;reward 口径
= `reward_mode="exposure"`,即 −mean(O₃))。S 为被打开的**干预单元**集合。

**单元**:T 个时段 × 2 个杠杆(NOx / VOC)⇒ n = 2T 个二元单元(mask 位 = t*2 + [0:NOx,1:VOC])。
枚举全部 2ⁿ 个联盟,做 **Möbius 变换** `m(S) = Σ_{T⊆S} (−1)^{|S\\T|} v(T)`
(等价于 FS 的"纯贡献 + 交互项"),回答三问:
1. 谱是否 **k-稀疏**(FS 文献 30 年先验:3 阶以上常可忽略)?
2. **跨时段交互** vs **同时段物种交互**:谁大?(预测:前者≈0、后者显著 —— 来自"记忆<1天"实测)
3. 谱若稀疏 ⇒ "MCTS 引导枚举"有戏;若稠密 ⇒ 应转做 block/Owen 聚合。

**动作杠杆**(与 `CB05Backend._emit_scalars` 对齐,互不污染):
- NOx 杠杆:置 `a_tra = level/nox_share_tra` ⇒ NOx 排放 ×(1−level);不动 VOC(a_ind=0)/SO2(a_pow=0);
- VOC 杠杆:置 `a_fvoc = level` ⇒ VOC 排放 ×(1−level);不动 NOx。

**成本**:2ⁿ 个联盟。共享**同一个 EnvModel**(键含 (step,state,action) ⇒ 前缀联盟自动复用)
后,实际积分次数 ≈ Σ_{k=1}^{T} 2^{2k}:
n=8(T=4) → 256 联盟 / ≈340 次积分;n=12(T=6) → 4096 / ≈5460 次。

用法:
    python -m atmos.experiments.exp_mobius_spectrum --blocks 4 --level 0.2
    python -m atmos.experiments.exp_mobius_spectrum --blocks 6 --level 0.2 --dt 6
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from ..configs import CB05Config, EnvConfig
from ..search.model import EnvModel

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "results", "mobius")


# --------------------------------------------------------------------------
# 值函数 v(S)
# --------------------------------------------------------------------------
def _action(level: float, nox_share: float, nox_on: bool, voc_on: bool) -> np.ndarray:
    """把一个 (NOx 开/关, VOC 开/关) 组合翻译成 CB05 的 4 维动作。

    只动 a_tra(NOx)与 a_fvoc(VOC),互不污染(不碰 a_ind/a_pow)。
    """
    a = np.zeros(4)
    if nox_on:
        a[1] = min(1.0, level / max(nox_share, 1e-9))
    if voc_on:
        a[3] = level
    return a


def _eval_subset(model: EnvModel, root: np.ndarray, mask: int, blocks: int,
                 level: float, nox_share: float) -> float:
    """v(S) = 暴露均值(每步 reward 均值)。mask 位 = t*2 + [0:NOx,1:VOC]。"""
    actions = []
    for t in range(blocks):
        nox_on = bool(mask & (1 << (t * 2)))
        voc_on = bool(mask & (1 << (t * 2 + 1)))
        actions.append(_action(level, nox_share, nox_on, voc_on))
    total = model.rollout_sequence(root, actions)
    return total / blocks


# --------------------------------------------------------------------------
# Möbius 变换(快速子集 DP)
# --------------------------------------------------------------------------
def mobius(v: np.ndarray, n: int) -> np.ndarray:
    m = np.asarray(v, float).copy()
    for i in range(n):
        bit = 1 << i
        for mask in range(1 << n):
            if mask & bit:
                m[mask] -= m[mask ^ bit]
    return m


# --------------------------------------------------------------------------
# 谱分析
# --------------------------------------------------------------------------
def _popcount(x: int) -> int:
    return bin(x).count("1")


def analyze(v: np.ndarray, m: np.ndarray, n: int, blocks: int) -> dict:
    """谱分析。**能量/稀疏度一律排除 m(∅)**(它只是基线常数,不是效应)。"""
    energy = m ** 2
    energy[0] = 0.0
    tot = float(energy.sum())          # 只统计 S≠∅ 的"效应能量"

    # ① 按阶能量占比(1..n 阶)
    order = {k: float(energy[[s for s in range(1, 1 << n) if _popcount(s) == k]].sum())
             for k in range(n + 1)}

    # ② 稀疏度:前 k 大系数占多少能量(排除 ∅)
    idx = np.argsort(-np.abs(m))[1:]   # 丢弃 m[0]
    cum = np.cumsum(energy[idx]) / max(tot, 1e-30)
    k90 = int(np.searchsorted(cum, 0.90) + 1)
    k95 = int(np.searchsorted(cum, 0.95) + 1)
    top = [{"mask": int(s), "order": _popcount(int(s)), "m": float(m[s])}
           for s in idx[: min(8, n)]]

    # ③ 一对交互:同时段物种对 vs 跨时段对
    same_t, cross_t, species_pairs = [], [], []
    for i in range(n):
        for j in range(i + 1, n):
            ti, tj = i // 2, j // 2
            mask = (1 << i) | (1 << j)
            if ti == tj:
                same_t.append(abs(float(m[mask])))
                species_pairs.append({"t": ti, "m": float(m[mask])})
            else:
                cross_t.append(abs(float(m[mask])))

    # ④ 一阶项(每个单元):按时段列出 NOx / VOC 纯贡献
    first = []
    for t in range(blocks):
        first.append({"t": t,
                      "nox": float(m[1 << (t * 2)]),
                      "voc": float(m[1 << (t * 2 + 1)])})

    def _stat(x):
        return {"n": len(x), "mean_abs": float(np.mean(x)) if x else 0.0,
                "max_abs": float(np.max(x)) if x else 0.0}

    return {
        "v_baseline": float(v[0]), "v_all_on": float(v[-1]),
        "improvement_all_on": float(v[-1] - v[0]),
        "m_empty": float(m[0]),
        "order_energy": order,
        "order_energy_frac": {k: (e / max(tot, 1e-30)) for k, e in order.items()},
        "top_coeffs": top,
        "k_for_90pct": k90, "k_for_95pct": k95, "n_effects": (1 << n) - 1,
        "pairs_same_time": _stat(same_t),
        "pairs_cross_time": _stat(cross_t),
        "species_pair_by_time": species_pairs,
        "first_order_by_time": first,
    }


# --------------------------------------------------------------------------
# 主流程
# --------------------------------------------------------------------------
def run(blocks: int, level: float, dt: float, seed: int,
        rtol: float = 1e-3, atol: float = 1e6,
        jac: bool = True, atol_rel: float = 0.0, atol_floor: float = 1e4,
        maxstep: int = 1000, substeps: int = 2, force_env: bool = False) -> dict:
    bcfg = CB05Config(rtol=rtol, atol=atol, jac=jac,
                      atol_rel=atol_rel, atol_floor=atol_floor,
                      maxstep=maxstep, substeps_per_control=substeps,
                      force_env_update=force_env)
    from ..chemistry.cb05 import CB05Backend
    b = CB05Backend(bcfg)
    b.set_episode(seed)
    cfg = EnvConfig(horizon_steps=blocks, dt_hours=dt, reward_mode="exposure")
    model = EnvModel(b, cfg)
    root = b.initial_state()
    nox_share = bcfg.nox_share_tra

    n = 2 * blocks
    v = np.zeros(1 << n)
    t0 = time.perf_counter()
    for mask in range(1 << n):
        v[mask] = _eval_subset(model, root, mask, blocks, level, nox_share)
    secs = time.perf_counter() - t0

    m = mobius(v, n)
    out = analyze(v, m, n, blocks)
    out.update({"blocks": blocks, "n_units": n, "n_coalitions": 1 << n,
                "level": level, "dt_hours": dt, "seed": seed,
                "nox_share_tra": nox_share,
                "rtol": rtol, "atol": atol, "jac": jac,
                "atol_rel": atol_rel, "atol_floor": atol_floor,
                "n_solver_fail": b.n_solver_fail,
                "last_solver_message": b.last_solver_message,
                "box_calls": model.box_calls, "seconds": secs,
                "v_spectrum": v.tolist(), "mobius_spectrum": m.tolist()})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", type=int, default=4, help="时段数 T(单元数 n=2T)")
    ap.add_argument("--level", type=float, default=0.2, help="单元减排比例(二元)")
    ap.add_argument("--dt", type=float, default=6.0, help="每个时段小时数")
    ap.add_argument("--seed", type=int, default=1, help="episode/气象种子")
    ap.add_argument("--rtol", type=float, default=1e-3, help="积分相对容差(数值噪声探针)")
    ap.add_argument("--atol", type=float, default=1e6, help="标量绝对容差(atol-rel<=0 时生效)")
    ap.add_argument("--atol-rel", type=float, default=0.0,
                    help=">0 启用逐物种分层 atol:max(atol_floor, atol-rel·|y|)")
    ap.add_argument("--atol-floor", type=float, default=1e4, help="分层 atol 的绝对下限(分子/cm³)")
    ap.add_argument("--no-jac", action="store_true", help="用数值 Jacobian(jac=False)对照解析 Dfun")
    ap.add_argument("--maxstep", type=int, default=1000, help="LSODA 内部步数上限(pykpp 默认 1000)")
    ap.add_argument("--substeps", type=int, default=2,
                    help="每个控制步的积分段数;段长=dt/substeps(断点对齐用 300s ⇒ substeps=dt*3600/300)")
    ap.add_argument("--force-update", action="store_true",
                    help="每段前强制刷新环境更新器(检验并规避 last_updated 跨调用泄漏)")
    args = ap.parse_args()

    print(f"[cfg] cb05 blocks={args.blocks} n_units={2*args.blocks} level={args.level} "
          f"dt={args.dt}h horizon={args.blocks*args.dt}h rtol={args.rtol:g} "
          f"atol={args.atol:g} atol_rel={args.atol_rel:g} floor={args.atol_floor:g} "
          f"jac={not args.no_jac} maxstep={args.maxstep} "
          f"substeps={args.substeps}(seg={args.dt*3600/args.substeps:.0f}s) "
          f"force_update={args.force_update}", flush=True)
    out = run(args.blocks, args.level, args.dt, args.seed, args.rtol, args.atol,
              not args.no_jac, args.atol_rel, args.atol_floor, args.maxstep,
              args.substeps, args.force_update)

    print(f"[run] 联盟={out['n_coalitions']} 箱模型积分={out['box_calls']} "
          f"耗时={out['seconds']:.1f}s 求解器失败={out['n_solver_fail']}"
          f"{('  (' + out['last_solver_message'][:60] + ')') if out['last_solver_message'] else ''}",
          flush=True)
    print(f"  暴露均值: 基线={out['v_baseline']:.2f}  全开={out['v_all_on']:.2f} ppb "
          f"(全干预改善={out['improvement_all_on']:+.2f} ppb)", flush=True)

    print("\n  ① 按阶能量占比(Möbius):")
    for k in range(out["n_units"] + 1):
        frac = out["order_energy_frac"][k]
        if frac > 1e-6:
            print(f"     阶 {k}: {frac:7.2%}", flush=True)
    print(f"  ② 稀疏度:{out['k_for_90pct']}/{out['n_effects']} 个效应系数覆盖 90% 能量,"
          f"{out['k_for_95pct']} 个覆盖 95%", flush=True)

    print("\n  ③ 一阶项(每个单元的纯贡献,ppb;>0 表示降暴露):")
    print("     时段    NOx       VOC")
    for f in out["first_order_by_time"]:
        print(f"     t={f['t']:<4d} {f['nox']:+8.3f}  {f['voc']:+8.3f}", flush=True)

    st, ct = out["pairs_same_time"], out["pairs_cross_time"]
    print(f"\n  ④ 二阶交互:同时段物种对(n={st['n']}) mean|m|={st['mean_abs']:.3f} "
          f"max|m|={st['max_abs']:.3f}  |  跨时段对(n={ct['n']}) "
          f"mean|m|={ct['mean_abs']:.3f} max|m|={ct['max_abs']:.3f}", flush=True)

    print("\n  ⑤ 前 8 大系数:")
    for c in out["top_coeffs"]:
        print(f"     m[{c['mask']:3d}]  阶={c['order']}  m={c['m']:+.3f}", flush=True)

    os.makedirs(OUT_DIR, exist_ok=True)
    path = os.path.join(OUT_DIR, "mobius_spectrum.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print(f"-> {path}", flush=True)


if __name__ == "__main__":
    main()
