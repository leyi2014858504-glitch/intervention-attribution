# 论文主线：摘要 + 图表清单

> 定位：GMD / Environmental Modelling & Software 首选，ACP（methods 框架）可冲。
> 主线一句话：**归因不能靠低阶截断 ⇒ 必须测全谱 ⇒ 给出代价对 n 平坦的采样器 ⇒ 谱揭示了局地方法看不见的结构。**

---

## 一、主张链（五个断言，图表必须各自承载一个）

| # | 断言 | 证据 | 承载图 |
|---|---|---|---|
| ① | 低阶截断不成立，谱必须实测 | 预注册 n=8 pilot 判 STOP；受限侧低阶只解释 R²=0.26 | Fig. 2 |
| ② | 全谱枚举是 2ⁿ，不可承受 | n=24 精确需 16,777,216 次 | Fig. 4 |
| ③ | 采样器代价对 n 几乎平坦 | n=12→24 仅 1.07×（1.5% 精度，274 次）；KernelSHAP 的 6–8× | Fig. 4 |
| ④ | 全谱才暴露的结构 | 块内符号翻转；配额最优随 regime 翻转；CO 弱玩家稳健性 | Fig. 5–7 + Table 3 |
| ⑤ | 整条链路可验（L1/L2 封死） | sham 已知零（1.07e-14）；V-门 6 个冻结值逐位复现 | Table 1 |
| ⑥ | 局地方法结构上看不见 ④，且**只在无关紧要的玩家上精确** | 相邻 3h 窗：局地同号、有限关停反号（−0.42 vs +6.98 ppb）；CO 精确到 2% 但总效应仅 0.14 ppb | Fig. 6 + Table 2 |

**注意 ⑥ 的位置**：它既是发现，又是定位武器。前五条说"我们造了什么、它是对的"，⑥ 说"不这么做会错"。

---

## 二、Abstract（英文稿）

**Title (draft):**
*Coalitional attribution of emission-control interventions: a cost-flat sampler and the interaction spectrum that low-order truncation misses*

**Abstract:**

Ozone control strategy asks which precursor reductions, in which sector and at which
time of day, actually lower ozone. The operational tools — adjoint models and
decoupled direct methods (DDM/HDDM) — answer this locally: they linearise at a base
state and return first (and, for HDDM, second) derivatives there. Real interventions
are not infinitesimal; a sector is switched off for a window, so the multiplier goes
from 1 to 0, not from 1 to 1+ε.

We recast the problem as a cooperative game over (sector × time-window) players and
compute the full Möbius/Walsh interaction spectrum, of which Shapley and Owen values
are known linear functionals. Three things follow.

First, low-order truncation fails, so the spectrum has to be measured rather than
assumed: a pre-registered eight-player pilot stops on its own criterion, and a
low-order truncation explains only R² = 0.26 of the response on the NOx-limited side.

Second, measuring the spectrum is a 2ⁿ problem, and we give a permutation-with-duality
sampler whose cost is essentially flat in the number of players — 274 runs reach ~1.5 %
precision at n = 24, 1.07× the n = 12 cost, against 16,777,216 runs for enumeration, and
6–8× more efficient than KernelSHAP on the same value function.

Third, the spectrum exposes structure that a local method cannot represent. Inside a
single six-hour block, NOx reductions in [12,15) and [15,18) are both predicted by the
local derivative to RAISE ozone; under finite intervention the first one LOWERS it
(−0.42 and +6.98 ppb respectively, with the linear extrapolation giving +3.03 and
+5.69). Aggregating to six hours hides the sign change entirely (0/12 sign disagreements
at 6 h versus 1/24 at 3 h). The divergence tracks whether the response is MONOTONE, not
how large the effect is: VOC windows stay within 11 % of the linear extrapolation and CO
windows, a third and weaker sector, within 2 %. But CO's entire effect is at most
0.14 ppb, so a local method is exact precisely where the player does not matter, while
the sector that dominates the outcome — NOx, 40–200× larger — is the one whose sign it
gets wrong. EKMA-style budget constraints additionally make the optimal allocation
regime-dependent — VOC-only on the NOx-saturated side, morning NOx plus afternoon VOC on
the NOx-limited side — and adding CO as a third sector moves the NOx/VOC attributions by
at most 0.18 ppb.

The pipeline is verified against physically forced zeros (a sham switch that makes four
CO windows null players; max |Möbius| = 1.1e-14) and a sub-lattice gate (the four-window
game is an exact sub-lattice of the eight-window game; six frozen values reproduced
bit-for-bit, including v = 21.182293).

We also report a data-driven validity diagnostic for the host box: two independent
nighttime constraints on its dilution timescale — O₃ titration requires ≈13 min, the
NOx budget requires ≈5.6 h — differ by a factor of 25, which proves that the
single-parcel assumption fails in the nocturnal boundary layer, consistent with the
known limitation of 0-D boxes and with the field's practice of excluding high-NOx
winter periods from them.

---

## 三、图表清单

### Fig. 1 — 问题的重述（示意图，无数据）
- **画什么**：左：领域现状——一个标量灵敏度 ∂O₃/∂NOx，在基点线性化。右：我们的对象——(部门 × 时间窗) 联盟格上的 2ⁿ 谱，Shapley/Owen 是它的线性泛函。
- **承载**：入口定义（不承载结果，但决定读者怎么读全文）
- **素材**：手绘/矢量图；代码依据 `_kppwork/exp_mobius.py`、`_kppwork/attrib_samplers.py`

### Fig. 2 — 低阶截断失效 ①
- **画什么**：截断阶数 k 对 R²（两条曲线：饱和侧 / 受限侧）。受限侧 0.26 是断点。
- **承载**：①；也是"必须实测全谱"的全部理由
- **素材**：`_kppwork/exp_mobius.py`、`mobius_results_dummy.txt`；pilot 的预注册判定在 spec 文档"⚠️ n=8 Möbius pilot"节
- **必须写进图注**：这是**预注册**判定 STOP，不是事后挑出来的

### Fig. 3 — 谱的形态（n=12 与 n=24 两格）
- **画什么**：Möbius 系数按阶分布（重尾、不稀疏）；两种 regime 分列，显示"饱和侧低阶主导 / 受限侧高阶主导"。
- **承载**：④ 的前半
- **素材**：`results_co12.txt`、`co_v.csv`、`results_scale.txt`、`scale_v.csv`

### Fig. 4 — 采样代价 vs n ③②
- **画什么**：横轴 n（8,12,24），纵轴运行次数（对数）。三条线：精确枚举（2ⁿ）、我们的采样器（达 ~1.5% 精度）、KernelSHAP。
- **关键数字**：枚举 4096 → 16,777,216（4096×）；我们 256 → 274（**1.07×**）
- **承载**：②③。这张图是方法贡献的主图
- **素材**：`_kppwork/exp_scale.py`、`results_scale.txt`、`attrib_bench.txt`

### Fig. 5 — 块内符号翻转（发现的主图）④
- **画什么**：8 个 3 小时窗的有限关停效应 Δ_off（柱状），按小时排列。12–15h 一根负柱紧邻 15–18h 一根正柱。
- **承载**：④；也是 ⑥ 的引子
- **素材**：`results_local_vs_global.txt`（3h 段）
- **⚠️ 口径一致性**：本轮用 V-门基线 + 度量 v48 得 −0.42 / +6.98 ppb；早先 n=24 扫描给出过 −1.60 / +11.32。**符号结构稳定、幅度依度量（v48 vs MDA8）而异**——全文必须统一到一个度量，并在图注写明是哪一个。

### Fig. 6 — 局地导数 vs 全局归因（定位武器）⑥
- **画什么**：三部门（NOx / VOC / CO）× 8 个 3h 窗。每个窗并排：(a) 局地导数 D_local（DDM-1 能给的），(b) −D_local（线性外推的关停效应），(c) 真实有限关停 Δ_off。NOx 12–15h 处 (b)(c) 反号，红框标出。
- **图注必须给的三条**：NOx \|D\|=1.0–5.7 ppb、ratio 0.99–2.21 且一处反号；VOC \|D\|=1.0–3.6、ratio ≤0.11；CO \|D\|=0.028–0.14、ratio ≤0.02。**结论句：局地方法只在无关紧要的玩家上精确——CO 精确到 2% 却只值 0.14 ppb，NOx 大 40–200× 却反号。**
- **承载**：⑥ —— 回答"为什么不用 DDM-3D"，并给出可判定规则（偏离跟随响应是否单调，而非效应大小）
- **素材**：`results_local_vs_global.txt`、`_local_vs_global.py`

### Fig. 7 — 配额约束下的最优配置随 regime 翻转 ④
- **画什么**：EKMA 式预算（0.50, 0.35）下，最优 (部门 × 时间窗) 配置在两侧 regime 的对比；预算 100% 时"平摊解"不存在的事实用一行注记说明。
- **承载**：④ 的可操作性那一块，也是唯一一句能被政策读者直接读的话
- **素材**：`_kppwork/exp_quota.py`、`results_quota.txt`

### Table 1 — 验证门（L1/L2 封死）⑤
- **列**：检验 | 期望 | 实测
- **行**：活性守卫（256/256）；A `max|v(mask)−v(mask&0xFF)|` = 0.000e+00；E `max|m(A)|`（触 CO 位）= 1.066e-14 / 0.000e+00；B φ_CO = 1.674e-14 / 0.000e+00；C 子博弈一致性 3.553e-14 / 0.000e+00；D `Σφ−[v(full)−v(none)]` = 9.948e-14 / 6.661e-16
- **另起一段**：V-门 6 个冻结值（21.182293 / 27.371391 / 34.342422 / 46.762071 / 10.247975 / 43.814426）逐位复现
- **素材**：`_kppwork/exp_sham.py`、`sham_result.txt`、`results_gate7.txt`

### Table 2 — Local vs global（审稿人正面回答表）⑥
- **列**：sector | hours | D_local | D_finite | −D_local | ratio | sign_ok
- **行**：24 个 3h 窗（NOx/VOC/CO 各 8）+ 12 个 6h 块，**全列不挑**
- **图注必须点出**：ratio = D_finite/(−D_local)，线性响应下恒为 1；符号不一致 1/24（3h）与 0/12（6h）
- **素材**：`results_local_vs_global.txt`

### Table 3 — 玩家集稳健性 ④
- **列**：φ（8 玩家）vs φ（12 玩家，含 CO 第三部门），逐玩家 |Δφ|
- **结论**：max |Δφ| ≤ 0.18 ppb；CO Σφ 仅占 ~1%
- **素材**：`_kppwork/exp_co12.py`、`results_co12.txt`

### Table 4 — 宿主模型有效域诊断（scope 声明，不是负结果）
- **列**：站点/季节 | BLH_day | BLH_night | τ_BLH | 夜间 O₃ 约束 τ | 夜间 NOx 约束 τ
- **结论**：两条夜间约束差 25×（≈13 min vs ≈5.6 h）⇒ 单气胞假设在夜间失效
- **引用**：Geyer & Stutz 2004（0-D 箱在 NBL 不足）；Wang et al. 2026 Sci Adv（0-D 箱排除 NOx > 50 ppb 的冬季时段）；Lu et al. 2017（2 天 spin-up + 风速缩放稀释率）
- **素材**：`results_kdil.txt`、`_blh_to_kdil.py`、`results_obs_dallas_2016w_cb05.txt`

---

## 四、Discussion 必须写的三段（否则会被拒）

1. **对位 adjoint / DDM / ISAM**（用 Table 2 撑）：局地 vs 全局；无穷小 vs 有限关停；拿不到高阶谱，且有反号反例。
2. **ISAM 的分摊任意性 = 我们的 L3 的工业级实例**：EPA 文档承认 "fair" 比例分摊产生 non-intuitive 结果（电厂烟羽里出现"生物源 NOx"），为此加了 5 种赋值选项，其中一种直接用 PH₂O₂/PHNO₃（Sillman 1995）作 regime 判据。这说明**归因依赖先验分摊选择**——正是 L3"公理是选择、不可验"。
3. **局限**：宿主为 0-D 箱、无干沉降、玩家集为手工设计而非清单派生、无真实大气锚定（Table 4 给出可计算的适用域判据，并把夜间明确划到范围外）。

---

## 五、局地方法可用性的判据（已完成的补充实验）

原计划的"把 CO 补进 Fig. 6"已经做完，36 个窗的结果给出一条可判定规则：

| sector | \|D_local\| 量级 | ratio 偏离 | 局地方法可靠吗 |
|---|---|---|---|
| NOx | 1.0–5.7 ppb | 0.99–2.21，一处**反号** | ✗ 有 EKMA 拐点 |
| VOC | 1.0–3.6 ppb | ≤0.11 | 基本可靠（单调） |
| CO | 0.028–0.14 ppb | ≤0.02 | ✓ 完全可靠，**但因此无关紧要** |

**要写进 Discussion 的一句**：偏离跟随的是**响应是否单调**，而不是效应大小；而"最需要精确的那个部门恰好就是失败的那个"。这句把 ⑥ 从"我们更好"升级成"我们知道什么时候可以不换"——后者审稿人更喜欢，因为它是可证伪的。

数据：`results_local_vs_global.txt`（24 个 3h 窗 + 12 个 6h 块），代码 `_local_vs_global.py`。
