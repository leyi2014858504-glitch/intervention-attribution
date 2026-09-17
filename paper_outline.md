# 论文主线：摘要 + 图表清单

> 定位：GMD / Environmental Modelling & Software 首选，ACP（methods 框架）可冲。
> 主线一句话：**归因不能靠低阶截断 ⇒ 必须测全谱 ⇒ 给出代价对 n 平坦的采样器 ⇒ 谱揭示了局地方法看不见的结构。**

---

## 一、主张链（五个断言，图表必须各自承载一个）

| # | 断言 | 证据 | 承载图 |
|---|---|---|---|
| ① | 低阶截断**不只是有地板，而是非单调、会发散** | 全枚举 n=8/12/16 × 双度量实测：v48 在 n=16 发散（8.9%→56.2%）；**MDA8 全级别发散**（n=16 q≤8 = 2722%，n=8 一阶截断已 118%）；谱质量峰在 n=16（k=8–10，17× 一阶），n=8/12 随阶衰减 | Fig. 2 + Table 6 |
| ② | 全谱枚举是 2ⁿ，不可承受 | n=24 精确需 16,777,216 次 | Fig. 4 |
| ③ | 采样器代价**远低于枚举，但并非与 n 无关**；且**小 n 时枚举反而更便宜** | 全枚举 n=8/12/16 实测：同等预算下误差 1.6–1.8×/倍 n ⇒ 同精度约 2.5–3×；枚举 16×/256×。交叉点在 n=8–12 之间 | Fig. 4 + Table 5 |
| ④ | 全谱才暴露的结构 | 块内符号翻转；配额最优随 regime 翻转；CO 弱玩家稳健性 | Fig. 5–7 + Table 3 |
| ⑤ | 整条链路可验（L1/L2 封死） | sham 已知零（1.07e-14）；V-门 6 个冻结值逐位复现 | Table 1 |
| ⑥ | 局地方法结构上看不见 ④，且**只在无关紧要的玩家上精确** | 相邻 3h 窗：局地同号、有限关停反号（MDA8 −0.91 vs +8.24 ppb；v48 同为反号）；CO 精确到 2% 但总效应仅 0.30 ppb；**普遍率**：EKMA 网格 8/25 点有反号（全 NOx）、2 点零效应失效（切线小 10×以上，`results_zero_effect.txt` 证实为物理脊线非噪声） | Fig. 6 + Table 2 + Table 8 |

**注意 ⑥ 的位置**：它既是发现，又是定位武器。前五条说"我们造了什么、它是对的"，⑥ 说"不这么做会错"。

**贡献分层（重要，别混）**

| 层 | 内容 | 状态 |
|---|---|---|
| **架构层 = 本文贡献** | ① 干预的**有限**博弈表述（部门 × 时间窗开关）② Möbius **全谱**（Harsanyi 红利）作为研究对象 ③ Shapley/Owen 是谱的**已知线性泛函** ④ 四层验证框架（L1/L2 可验、L3 公理选择、L4 外部效度）⑤ 诊断：截断发散、局地导数失效、regime 依赖 | 全部有实测支撑 |
| **实现层 ≠ 贡献** | 用哪个采样器估计谱 | 实测 n=16 上 `perm_anti`（标准、结构盲）**优于** hub-分层主动采样 4–6×（差距随 n 扩大），故**直接采用 `perm_anti`**，不再声称"结构感知主动采样"是方法贡献 |

采样器是可替换件：换掉它，架构一分钱不损。原"OURS"臂的失败应写入 Discussion，作为"分层粒度随 n 线性增长、收益不增长"的教训。

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
compute the full Möbius interaction spectrum (the Harsanyi dividends), of which Shapley and Owen values
are known linear functionals. Three things follow.

First, low-order truncation does not merely have a floor — it is NON-MONOTONE and can
diverge, so the spectrum has to be measured rather than assumed. Fully enumerating
n = 8, 12 and 16 lets the truncation error be computed exactly at every order, under two
metrics. On the NOx-saturated case with a daily-mean metric, truncating the spectrum at
order q — what low-order attribution implicitly does — leaves a floor of 6.7 % (q=3) at
n = 8 and 6.9 % at n = 12, and at n = 16 it bottoms at 8.9 % (q=2) before climbing to
56 % (q=7). Under the regulatory MDA8 metric — itself a max functional, which concentrates
rather than dilutes the afternoon response — the effect is far stronger: first-order
truncation already carries 118 % at n = 8, and at n = 16 the error climbs through 113 %
(q=4) to 2722 % (q=8). The reason is visible in the spectrum itself: the interaction mass
Σ_{|A|=k}|m(A)|/|A| decays with k at n = 8 and 12 but peaks at k = 10 — 17× the
first-order mass — at n = 16. High-order terms are individually large and cancel;
truncating away the cancellation is what destroys the estimate.

Second, measuring the spectrum is a 2ⁿ problem, and we use a STANDARD permutation-with-
antithetic estimator, characterising its cost rather than claiming it. Its cost grows far
more slowly than enumeration but is NOT independent of n.
Fully enumerating n = 8, 12 and 16 (256, 4096 and 65536 box runs) gives a measured,
oracle-backed scaling curve rather than an extrapolation: doubling the player count —
whether by refining the time axis (n = 8 → 16, sectors held at 2) or by adding a sector
(n = 8 → 12, windows held at 4) — raises the error at equal run count by 1.6–1.8×, so
equal precision costs roughly 2.5–3× more, against 16× and 256× for enumeration. The sampler
only overtakes enumeration between n = 8 and n = 12; above that its margin grows with n.
The residual growth is the combinatorial weight C(n−1, k−1) carried by order-k Möbius
terms, only partly cancelled by the decay of the spectrum.

Third, the spectrum exposes structure that a local method cannot represent. Inside a
single six-hour block, NOx reductions in [12,15) and [15,18) are both predicted by the
local derivative to RAISE ozone; under finite intervention the first one instead LOWERS
it (−0.91 and +8.24 ppb respectively, with the linear extrapolation giving +6.37 and
+8.20). Aggregating to six hours hides the sign change entirely (0/12 sign disagreements
at 6 h versus 1/24 at 3 h). The divergence tracks whether the response is MONOTONE, not
how large the effect is: VOC windows stay within 13 % of the linear extrapolation and CO
windows, a third and weaker sector, within 2 %. But CO's entire effect is at most
0.30 ppb, so a local method is exact precisely where the player does not matter, while
the sector that dominates the outcome — NOx, one to two orders of magnitude larger — is
the one whose sign it gets wrong. The sign flip is metric-independent: it appears
identically under the daily-mean statistic (same 1/24, 0/12 counts). EKMA-style budget
constraints additionally make the optimal allocation regime-dependent — VOC-only on the
NOx-saturated side, morning NOx plus afternoon VOC on the NOx-limited side — and adding
CO as a third sector moves the NOx/VOC attributions by at most 0.10 ppb per player
(sector totals by 0.31 ppb).

The pipeline is verified against physically forced zeros (a sham switch that makes four
CO windows null players; max |Möbius| = 1.1e-14) and a sub-lattice gate (the four-window
game is an exact sub-lattice of the eight-window game; six frozen values reproduced
bit-for-bit, including v = 21.182293). The same bit-for-bit regression gates later caught
— and localized to 18 rows — an interleaved-write corruption in an earlier 12-player
table, which the dual-metric re-run exposed and repaired: the verification layer protects
the data, not only the code.

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

### Fig. 2 — 低阶截断**非单调且会发散** ①（双度量）
- **画什么**：双面板。(a) 截断误差 vs 阶数 q，n=8/12/16 三条线，**每面板一个度量**（d3_mean | MDA8）——MDA8 的 n=16 那条**随 q 爆炸**是主视觉；(b) 谱质量 `Σ_{|A|=k}|m(A)|/|A|` 按 k 分布（以 k=1 归一）：n=8/12 随阶衰减 vs n=16 峰在 k≈10。
- **实测（r_sat，全枚举 oracle；2026-09-16 重建 co_v.csv 后的干净数字）**：

  d3_mean（v48）：
  | q | n=8 | n=12 | n=16 |
  |---|---|---|---|
  | ≤1 | 34.4% | 35.6% | 39.8% |
  | ≤2 | 17.7% | 18.6% | **8.9%** |
  | ≤3 | **6.7%** | 6.9% | 13.2% |
  | ≤4 | 8.7% | 9.8% | 16.1% |
  | ≤6 | 5.7% | 6.9% | 41.9% |
  | ≤7 | 0.3% | 0.5% | **56.2%** |

  MDA8（法定度量，主图）：
  | q | n=8 | n=12 | n=16 |
  |---|---|---|---|
  | ≤1 | **117.6%** | 119.2% | 81.6% |
  | ≤2 | 54.3% | 55.0% | 58.1% |
  | ≤3 | 37.2% | 39.2% | 56.2% |
  | ≤4 | 65.0% | 63.4% | **113.4%** |
  | ≤6 | 10.0% | 16.7% | **877.6%** |
  | ≤8 | 3.6%(q≤7) | 0.3% | **2722.0%** |

  r_lim1 同结构：v48 n=16 q≤3 最低 2.0% → q≤6 32.7%；MDA8 n=16 q≤3 9.6% → q≤7 170.9%。
- **新机制句（MDA8 版，必须写）**：MDA8 本身是 max 泛函（8h 滚动均值的日最大），它**集中**午后响应而非稀释——法定度量恰是交互结构最强的度量。一阶截断误差 >100% 意味着"只保留主效应"比"什么都不说"还差。
- **机制句（原版，保留）**：高阶 Möbius 项**各自很大但相互抵消**；截断掉的是抵消，不是小量。所以误差不降反升。
- **承载**：① ——"必须实测全谱、不能假设低阶"的全部理由；这是**架构**主张而非采样器主张
- **素材**：`_kppwork/_n16_spectrum.py`、`_kppwork/spectrum_orders.txt`、`_kppwork/_m8_spectrum.py`、`_kppwork/spectrum_orders_m8.txt`、`mobius_results_dummy.txt`
- **⚠️ 更正记录**：旧 n=12 列（26.7/79.7/1032.9%）**作废**——来自 co_v.csv 的 18 行交织写坏（WSL 崩溃恢复产物），双度量重跑被逐位回归门抓获并重建；干净 n=12 与 n=8 同形（地板 ~6.9%）。发散证据现在完全由 **n=16（两度量）+ MDA8 全级别**承载，更强而非更弱。
- **图注要点**：引 MDA8 面板做主视觉；v48 面板放补充材料，两度量并置本身就是"结论不依赖度量选择"的证据

### Fig. 3 — 谱的形态（n=12 与 n=24 两格）
- **画什么**：Möbius 系数按阶分布（重尾、不稀疏）；两种 regime 分列，显示"饱和侧低阶主导 / 受限侧高阶主导"。
- **承载**：④ 的前半
- **素材**：`results_co12.txt`、`co_v.csv`、`results_scale.txt`、`scale_v.csv`

### Fig. 4 — 采样代价 vs n ③②
- **画什么**：横轴 n（8, 12, **16**, 24），纵轴"达到归一化精度 X% 所需运行次数"（双对数）。两条线：精确枚举（2ⁿ）、我们的采样器。n=8/12/16 三点**全枚举、有真 oracle**；n=24 只有 spread。
- **实测（r_sat，`err(run)` = 单次运行误差，全部对真 oracle）**：

  | B | runs (n=8/16) | n=8 | n=16 | 比值 |
  |---|---|---|---|---|
  | 64 | 65 / 75 | 7.304% | 11.536% | 1.58× |
  | 128 | 132 / 137 | 4.463% | 7.642% | 1.71× |
  | 256 | 256 / 267 | 1.108%（**饱和**） | 4.894% | — |
  | 1024 | — / 1030 | — | 2.183% | — |

  部门增加（n=8→12）：B=64 7.304%→12.133%（1.66×），B=128 4.463%→8.186%（1.83×）。
- **两条干净阶梯**（**重要**：n=12 与 n=16 不嵌套，不可直接比）
  - 时间细化 n=8（2×4）→ n=16（2×8），部门固定
  - 部门增加 n=8（2×4）→ n=12（3×4），窗数固定
- **关键数字**：两种机制同等预算下误差都涨 **1.6–1.8×** ⇒ 同精度约 **3×**；枚举涨 16×（n=8→12）/ 256×（n=8→16）
- **交叉点（反直觉，值得单独一句）**：达到 1.5% 误差，n=8 需枚举 256 次 vs 采样 ~620 次 ⇒ **小 n 时枚举更便宜**；n=12 是 4096 vs ~1370，n=16 是 65536 vs ~1850 ⇒ 采样器优势随 n 增长
- **⚠️ 口径警告（必须写进图注）**
  1. n=8 的 B=256 是**全域饱和**（2⁸=256），误差被人为压低，**不可用于外推**
  2. n=24 的误差是 `spread × (err/spread 来自 n=12)` **映射**的，不是实测。实测的时间细化惩罚（1.6–1.8×）**大于**映射值（1.28×）⇒ 映射确实偏乐观
- **修正记录**：① 早期"274 次达 1.5%，仅 1.07×"——那是拿 1.44% 的价钱与 1.85% 的价钱相比，口径不等价，**作废**；② 中期"×1.28 ⇒ 1.8×"基于映射值，已被实测取代，**作废**
- **承载**：②③。这张图是方法贡献的主图
- **素材**：`_kppwork/exp_n16.py`（实测，主）、`_kppwork/scaling_n16.txt`、`_kppwork/gen_n16.py`、`_kppwork/sweep_n16.sh`；`_kppwork/exp_scale.py` + `results_scale.txt`（映射版，保留作对照）

### Fig. 5 — 块内符号翻转（发现的主图）④
- **画什么**：8 个 3 小时窗的有限关停效应 Δ_off（柱状），按小时排列。12–15h 一根负柱紧邻 15–18h 一根正柱。
- **承载**：④；也是 ⑥ 的引子
- **素材**：`results_local_vs_global.txt`（3h 段）
- **度量**：主口径 MDA8，得 −0.91 / +8.24 ppb；v48 为 −0.42 / +6.98。**符号结构两度量完全一致**（1/24、0/12 相同），幅度依度量而异——图注声明用 MDA8，v48 进补充材料做稳健性对照。

### Fig. 6 — 局地导数 vs 全局归因（定位武器）⑥
- **画什么**：三部门（NOx / VOC / CO）× 8 个 3h 窗。每个窗并排：(a) 局地导数 D_local（DDM-1 能给的），(b) −D_local（线性外推的关停效应），(c) 真实有限关停 Δ_off。NOx 12–15h 处 (b)(c) 反号，红框标出。
- **图注必须给的三条（MDA8 主口径）**：NOx \|D_local\|=2.1–8.2 ppb、ratio −0.14–1.34 且一处反号；VOC \|D_local\|=2.3–8.5、ratio ≤0.13；CO \|D_local\|=0.06–0.29、ratio ≤0.02。**结论句：局地方法只在无关紧要的玩家上精确——CO 精确到 2% 却只值 0.30 ppb，NOx 大一到两个数量级却反号。**
- **反号的度量不变性（必须点出）**：v48 下同一窗反号、计数同为 1/24（3h）与 0/12（6h）——核心定性结论不依赖度量选择。⚠️ 注意 v48 的幅度偏离更极端（ratio 达 2.21），MDA8 更温和（≤1.34）但绝对效应更大；两度量并置讲"符号翻转稳健、幅度依度量"。
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
- **列**：sector | hours | D_local | D_finite | −D_local | ratio | sign_ok（**主表 = MDA8**；v48 同格式表进补充材料）
- **行**：24 个 3h 窗（NOx/VOC/CO 各 8）+ 12 个 6h 块，**全列不挑**
- **图注必须点出**：ratio = D_finite/(−D_local)，线性响应下恒为 1；符号不一致 1/24（3h）与 0/12（6h），**两度量计数相同**
- **素材**：`results_local_vs_global.txt`（双度量已并置）

### Table 3 — 玩家集稳健性 ④
- **列**：φ（8 玩家）vs φ（12 玩家，含 CO 第三部门），逐玩家 |Δφ|
- **结论**（2026-09-16 干净表重跑）：逐玩家 max |Δφ| = **0.10 ppb**（旧值 0.18 出自被污染的 co_v.csv，作废）；部门总量移动 ≤ 0.31 ppb；CO Σφ 仅占 ~1%
- **素材**：`_kppwork/exp_co12.py`、`results_co12.txt`

### Table 4 — 宿主模型有效域诊断（scope 声明，不是负结果）
- **列**：站点/季节 | BLH_day | BLH_night | τ_BLH | 夜间 O₃ 约束 τ | 夜间 NOx 约束 τ
- **结论**：两条夜间约束差 25×（≈13 min vs ≈5.6 h）⇒ 单气胞假设在夜间失效
- **引用**：Geyer & Stutz 2004（0-D 箱在 NBL 不足）；Wang et al. 2026 Sci Adv（0-D 箱排除 NOx > 50 ppb 的冬季时段）；Lu et al. 2017（2 天 spin-up + 风速缩放稀释率）
- **素材**：`results_kdil.txt`、`_blh_to_kdil.py`、`results_obs_dallas_2016w_cb05.txt`

### Table 5 — 采样的实测缩放曲线（n=8/12/16 全枚举，真 oracle）③
- **列**：n | B（预算） | runs（实际化学次数） | err(run)（单次运行误差） | spread | 备注
- **行**：r_sat 的 n=8（3 点）、n=12（5 点）、n=16（5 点）；r_lim1 同结构
- **必带的两个注**
  1. `err(run) = mean over seeds of max_i |est_i − true_i| / max|φ|`，测的是"**跑一次**离真值多远"，不是"多 seed 均值的偏差"（后者测无偏性，不随预算单调，用它会读出错误的缩放律）
  2. n=8 的 B=256 是全域饱和，需在图/表中标出或剔除
- **另附**：Möbius 快速实现的交叉校验 `max|fast − exp_scale| = 1.6e-14`（n=8）；两个冻结锚点 `m00000 → 21.182293`、`m00012 → 27.371391` 逐位复现
- **素材**：`_kppwork/scaling_n16.txt`、`_kppwork/n16_exact_shapley.json`、`_kppwork/exp_n16.py`

### Table 6 — 截断地板与谱质量（全枚举 n=8/12/16 × 双度量）①
- **列**：n | max\|φ\| | trunc(q≤1..8) | 谱质量 `Σ_{|A|=k}|m(A)|/|A|` 按 k（k=1 归一）
- **行**：r_sat 与 r_lim1 各 3 个 n × 2 个度量（数字以 Fig. 2 双表为准，2026-09-16 重建后）
- **max|φ| 对照**：v48 16.66/16.61/13.89 vs MDA8 7.64/7.56/6.35 ppb（n=8/12/16）——MDA8 下主效应相互抵消更多、尾部更重，正是截断爆炸的来源
- **必带的机制句**：高阶项各自大、相互抵消；截断破坏抵消 ⇒ 误差非单调、可发散。MDA8 的 max 泛函性质放大之（法定度量恰是交互最强的度量）
- **与 Table 5 的区分**：Table 5 是**采样器**的代价（实现层），本表是**谱结构**的性质（架构层）——两者不可混引
- **素材**：`_kppwork/spectrum_orders.txt`、`_kppwork/spectrum_orders_m8.txt`、`_kppwork/_n16_spectrum.py`、`_kppwork/_m8_spectrum.py`

### Table 8 — Scan A：反号普遍率 × EKMA 网格（⑥ 的 prevalence 层，`results_scan_flips.txt`）
- **设计**：5×5 (NOX_FACTOR × VOC_FACTOR) ∈ [0.5, 2]²，每点跑 Fig. 6 同一协议（3h 窗、Δ=5%、双度量），25×73 = 1825 次运行；哨兵点 n100v100 必须复现基线反号——**已通过**。
- **反号计数网格（MDA8，/24）**：
  | NOx↓ VOC→ | 0.5 | 0.75 | 1.0 | 1.5 | 2.0 |
  |---|---|---|---|---|---|
  | 0.5 | 0 | 0 | 1 | 0 | 0 |
  | 0.75 | 1 | 1 | 1 | 0 | **5** |
  | 1.0 | 0 | 0 | **1** | 1 | 0 |
  | 1.5 | 0 | 0 | 0 | 0 | 1 |
  | 2.0 | 0 | 0 | 0 | 0 | 0 |
- **三条结论**：
  1. **8/25 格有反号，全部在 NOx 部门**；NOx_F = 2.0 整行为零（深度饱和侧反而安全），反号聚集在 NOx_F ≤ 1.5 带内、最密在 0.75 行（5 格中 4 格有），热点 (0.75, 2.0) 达 5/24——反号不是孤例，是**脊线可及区的结构性特征**；普遍性声明从"1/24"升级为这张表。
  2. **单调性判据全网格成立**：VOC 偏离 ≤0.21、CO 偏离 ≤0.01 在 25 个点无一例外——局地方法只在跨脊线的部门失效，与效应大小无关。
  3. **零效应失效（第三种、也是最危险的失效模式）**：切线比真值小 10 倍以上、局地方法把大效应报成"无影响"——网格上 **2 处，全在 NOx**（n050v100@9-12、n075v100@12-15）。n075v100 的 NOx@12-15：D1 = −0.089 ppb（切线，"无影响"），真值 **−3.63 ppb**，ratio = **−41**，且 |D1| 是中心差分截断地板的 **219 倍**——**是物理脊线（∂O₃/∂NOx = 0）而非数值噪声**。危害等级最高：报零的政策建议"这个部门不重要"，不会有人复查。
- **统计量更正（2026-09-17，必读）**：Table 8 原 `NOxdev` 列用 max\|ratio − 1\| 作指标，而 **ratio = D_finite/(−D1) 在 D1→0 时发散**——脊线点上该列爆炸（n075v100 达 42.09），是**指标缺陷而非物理量级**。改用**绝对误差** `|D_finite − (−D_local)| / max|D_finite|`（脊线安全）：**NOx 0.21–0.85、VOC ≤0.05、CO ≤0.01**。反号计数（8/25）是符号比较、不受影响；已作废的只有 ratio 幅度列。
- **三类失效模式的危险性排序（Discussion 用，全部有实测）**：
  | 模式 | 局地说 | 真相 | 危险度 |
  |---|---|---|---|
  | 幅度错 | −3.03 | −0.42 | 中（方向对） |
  | **符号反** | +3.03 | −0.42 | 高 |
  | **零效应** | −0.09 ≈ 0 | **−3.64** | **最高（"不重要"无人复查）** |
  一句话：**局地方法最危险的失效，正好发生在它最有信心的地方——报零**。
- **度量对比**：MDA8 反号严格多于 v48（n050v100: 1 vs 0；n075v200: 5 vs 1）——max 泛函对午后峰更敏感；主文用 MDA8。
- **素材**：`_kppwork/gen_flipscan.py`、`_kppwork/sweep_flips.sh`、`_kppwork/_flip_report.py`、`_kppwork/_zero_effect.py`、`results_scan_flips.txt`、`results_zero_effect.txt`
- **Scan B 已裁决不做（首篇）**：n=16 全枚举 × 更多 regime 点只是确认性稳健检查，无方法学增量；发散主张已由"脊线两侧 × 双度量"的 bracket 封住。**网格加密 + 真实城市坐标 + 观测弦（封锁/治理天然实验）打包为第二篇**（"divergence map + observational chords"），首篇在 future work 一句话指路。

### Table 7d — HDDM（二阶 Taylor）挡得住吗：单轴 + 全组合实测（`results_hddm_check.txt`、`results_hddm_coalitions.txt`）
- **动机**：审稿人最强的一击是"我们用 HDDM 不用 DDM，二阶能捕获曲率"。必须用数字回答，不能靠论证。**零新运行做单轴**（109 次已有运行的轨迹里有 ±δ 与 off），**512 次/情景新运行做组合**（测 152 个系数，然后对全 65535 个联盟打分）。
- **单轴结果（诚实，部分削弱原武器）**：二阶**修复了**一阶的符号错误——MDA8 与 d3_mean 双双 1/36 → 0/36；max 误差 0.712 → 0.154（MDA8）。**机制（关键）**：多线性函数沿单轴是直线 ⇒ **D2_ii ≡ 0 是多线性的判据**；实测 D2_ii ≠ 0（NOx 12–15h 达 −13.92，是 D1 的 2.2 倍）证明轴内曲率存在，而单轴的**全部**偏离源恰是这一项，二阶正好抓住它。
- **组合结果（新，决定性）**：全部 2¹⁶−1 个联盟，两个 regime，双度量：

  | 情景/度量 | 指标 | DDM-1 | **HDDM-2** |
  |---|---|---|---|
  | r_sat / MDA8 | 中位误差 | 15.4% | **5.7%** |
  | r_sat / MDA8 | 符号翻转率 | 17.88% | **7.14%** |
  | r_sat / MDA8 | 误差 >50% 的联盟 | — | **1731/65535 (2.6%)** |
  | r_lim1 / MDA8 | 中位误差 | 9.7% | **2.5%** |
  | r_lim1 / MDA8 | 符号翻转率 | 0% | **0%** |

- **⚠️ 结论修正（重要，不要过度声称）**：**"HDDM-2 在组合上失败"不成立**——二阶整体显著优于一阶且不爆炸。能站住的是三条更窄的陈述：
  1. **尾部仍失效**：r_sat 上 7.1% 的联盟符号翻转、2.6% 误差 >50%，且**随 k 单调恶化**（k≥10 翻转率 11–50%）；一阶的失败被缓解而非消除。
  2. **regime 依赖，且与全文一致**：NOx 受限侧（r_lim1）二阶近乎完美（0%），失效**只在 NOx 饱和/山脊可及侧（r_sat）**——与截断发散、局地反号、采样重尾**同一区域**，这是 Table 7d 最有价值的贡献：四条独立证据指向同一物理判据。
  3. **二阶无法自我诊断**：它给出的是单点预测，没有 holdout 式的误差棒；采样方法在同成本下既无该尾部又能自校准（Table 7a）。
- **成本对比（必须写准确，勿夸大）**：本实验用 512 次运行拟合 152 个系数（有限差分）。真实 CMAQ-HDDM 是**灵敏度方程与主 ODE 同时积分**，故非 152× 成本，但灵敏度变量数随 n² 增长（n=16→136 个二阶、n=24→300、n=48→1176），内存与稳定性是实际瓶颈，且文献中绝大多数应用只做一阶。我们的采样成本按 Table 5 实测为 1.6–1.8×/倍 n。**两条成本曲线的标度不同（n² vs n）**——这是可画成图的一手对比。
- **数据可信度（必读）**：首轮 sweep 把两个情景并发在**同一 kpp7 构建目录**上，worker 互相覆盖 `cb05.nml`/`traj.dat`，写出了**看似合理但错误**的值（NaN 过滤抓不住）。修法：情景串行 + 清表重跑，并加**逐位交叉验证**（单轴 ±δ 与 lvg 同配置运行，16 对，max|Δ|=4.7e-07）。另修正分析脚本一处真 bug：`∂²/∂x²` 的中心差分需**绝对基准值**，误用基线归零值会把 D2 放大 ~2·base/δ²（曾致 74 万% 假象）。
- **素材**：`_kppwork/gen_hddm.py`、`_kppwork/sweep_hddm.sh`、`_kppwork/_hddm_check.py`、`_kppwork/_hddm_coalitions.py`、`_kppwork/_hddm_xcheck.py`、`results_hddm_check.txt`、`results_hddm_coalitions.txt`

### Table 7a — 旗标校准：从"会排序"到"能报数"（`selfdiag_calib.txt`）
- **问题**：V3 只证了秩相关（会排序）；3D 上实践者要的是**数值**——"我这次误差多少 ppb"。旗标必须校准，不能只单调。零新运行，纯后处理（148 + 100 个配对样本）。
- **做法**：多角 × 多预算 × 多阶 × 多 seed 收集 (旗标, 真误差) 配对 → 拟合 log e = a + b log(flag) → **留一角交叉验证（LOO）**：在 5 个角拟合，预测第 6 个没见过的角。
- **实测（关键三行）**：

  | 度量 | 旗标 | c = median(e/f) | 跨角 c 散布 | 斜率 b | LOO 迁移 e/pred |
  |---|---|---|---|---|---|
  | **MDA8** | **holdout h** | **0.327** | **1.2×** | **1.03** | **1.00 [0.91, 1.10]** |
  | MDA8 | bootstrap sp | 2.02 | 2.5× | 0.21 | 1.48 [1.13, 1.65] |
  | MDA8 | 效率违反 g | 0.92 | 3.2× | 0.44 | 1.07 [0.59, 1.29] |
  | v48 | holdout h | 0.278 | 1.5× | 1.08 | 0.97 [0.80, 1.14] |
  | v48 | bootstrap sp | 2.30 | 3.3× | 0.51 | 1.39 [1.04, 1.91] |
  | v48 | 效率违反 g | 0.16 | **21.7×** | 0.28 | 1.47 [0.40, 2.89] |

- **三条结论**：
  1. **holdout 是唯一可校准的旗标**：斜率 b≈1（e 与 h 成正比，不是幂律怪形），跨角 c 散布仅 1.2×（MDA8），**LOO 在没见过的化学环境中预测真误差，中位误差比 1.00、范围 ±10%**。这直接把"旗标管用"从排序主张升级为**定量主张**。
  2. **交付配方（无真值误差棒）**：e_est ≈ 0.33·h（中位），e_95 ≈ **0.6·h**（95% 覆盖，k95=1.8）。3D 上跑一层 holdout 就能报"归因误差 < 0.6h，95% 置信"——**采样归因可以自带统计诚实的误差棒，这是本文的独立方法学贡献**。
  3. **另两个旗标被定量否决**：bootstrap spread 在 v48 上跨角漂移 3.3×（老"弱旗标"结论的定量版）；效率违反 21.7×（v48）——它们的 LOO 看起来还行纯属巧合（斜率 <1 的补偿），不能用于报数。
- **与 V1 的关系（重要）**：V1 说"截断模型 + 全局和旗标 = 看不见发散"；本表说"**换成拟合模型 + holdout 旗标，误差可定量预报**"。两句合起来才是完整建议：**别截断、别用全局和；用拟合 + holdout + 系数 0.33/0.6**。
- **诚实边界**：所有角共享同一 CB05 箱，"跨角"验证的是**化学环境间**迁移，不是**模型种类间**；CTM 迁移需在第二篇用少量 3D 全枚举点（n=8 级）做哨兵复核。最便宜的中间步是箱 vs 1D 柱。
- **素材**：`_kppwork/exp_calib.py`、`_kppwork/selfdiag_calib.txt`、`_kppwork/_run_calib.sh`

### Table 7b — 无真值自诊断旗标的 oracle 验证（3D 迁移的可行性实验，`selfdiag.txt`）
- **动机**：3D 上永远没有 oracle，"截断/采样误差多大"必须**不依赖真值**估计。0-D 表给了真值，是唯一能验证旗标本身的环境。零新运行，纯后处理。
- **旗标**：g_q 效率违反（截断模型对不上总账）、h_q holdout 失拟、sp_q bootstrap spread。预注册判据 V1/V2/V3 先于看数据写死。
- **V1 FAIL（有价值的负结果，双度量已闭环）**：全局和类旗标**看不见发散**——v48：n=16 发散段（t: 8.9→16.1%）g_q 反而降（12.1→1.0%）；**MDA8（头条度量）更决绝**：r_lim1 n=12 的 ρ(g_q)=**−1.00**（t: 51→97%，g_q: 3.0→1.1%），r_sat n=16 ρ(g_q)=0.00。机制：g_q 测总尾质量 Σ_{|A|>q}m(A)，高阶抵消使其恒小；发散住在**逐玩家尾部**。**结论：截断误差不可自我认证——"必须实测全谱"从方法偏好升级为逻辑必然。**
- **意外正发现（MDA8 下更夸张）**：同阶同预算，**拟合 q 阶模型远优于截断**——MDA8 差距达 **30–50×**（r_sat n=16 q≤4：e=3.3% vs t=113%；r_lim1 n=12 q≤4：1.9% vs 96.5%），v48 下为 10–25×。大联盟样本值携带尾部信息，LS 吸收之。**截断不仅发散，还是在 q 阶模型里选错了估计量。**
- **V3（MDA8 下收窄）**：holdout 是唯一双度量存活的旗标（0.59–0.92）；效率违反在 r_sat n=16 崩到 0.06（不可靠）；bootstrap spread 全面弱（0.05–0.47）。**3D 协议定稿：归因用 perm_anti + 其自带方差；模型拟合用 holdout 认证；截断归因 + 全局和旗标 = 双重错误。**
- **措辞联动**：摘要原句 "the best achievable q-th-order model" **作废**（拟合模型就是反例），已改为"truncating the spectrum ... what low-order attribution implicitly does"。
- **素材**：`_kppwork/exp_selfdiag.py`（支持 v_ppb/v_mda8 两列）、`_kppwork/selfdiag.txt`、`_kppwork/selfdiag_m8.txt`

### Table 7 — 两个 n=8 benchmark 在 n=16 oracle 上重跑（讨论用）
- **动机**：两个 benchmark 原都在 n=8，宇宙仅 256、预算达 192（覆盖 75%），**方法无法分离**——与 Table 5 的交叉点（枚举在 n=8 反而便宜）是同一事实。n=16 的预算 64 只覆盖 0.098%。
- **保真性**：n=8 侧**逐格复现**已发表数字（attribution 20/20 格全中；planning 的 greedy/greedy+ 逐位一致）
- **归因侧结果（n=16）**：方法确实分离了，但 `perm_anti`（标准、结构盲）在 20%/10% 容差上**胜过** hub-分层主动采样 6×/4×；差距随 n 扩大。机制：分层粒度 = n，覆盖成本随 n 线性增长而收益不增。
- **规划侧结果（n=16）**：MCTS 在**每一个 k** 上仍是最差臂（R1：k=64 时 greedy+ 0.000 / random 2.783 / **mcts 2.915**）。搜索空间从 256 涨到 65,536（256×）后负结果**依然成立** ⇒ 原判定从"n=8 太小"的嫌疑中解脱，**升级为已验证事实**
- **素材**：`_kppwork/bench_n16.txt`、`_kppwork/exp_n16_bench.py`、`_kppwork/sweep_n16.sh`
- **保留**：planning 的 random/mcts 用 50 seed，原文件默认 60，故二者略有偏差（greedy 系确定性、逐位一致）；要严格复现需 60 seed

---

## 四、Discussion 必须写的三段（否则会被拒）

1. **对位 adjoint / DDM / ISAM**（用 Table 2 撑）：局地 vs 全局；无穷小 vs 有限关停；拿不到高阶谱，且有反号反例。
2. **ISAM 的分摊任意性 = 我们的 L3 的工业级实例**：EPA 文档承认 "fair" 比例分摊产生 non-intuitive 结果（电厂烟羽里出现"生物源 NOx"），为此加了 5 种赋值选项，其中一种直接用 PH₂O₂/PHNO₃（Sillman 1995）作 regime 判据。这说明**归因依赖先验分摊选择**——正是 L3"公理是选择、不可验"。
3. **局限**：宿主为 0-D 箱、无干沉降、玩家集为手工设计而非清单派生、无真实大气锚定（Table 4 给出可计算的适用域判据，并把夜间明确划到范围外）。
4. **采样器的教训（诚实披露，用 Table 7）**：hub-位置分层的**开销随 n 线性增长**（层数 = n），而它买到的方差削减不随 n 增长，因此该构造在 n=16 上被**结构盲的**对偶排列采样反超 4–6 倍，且差距随 n 扩大。两个 n=8 benchmark 原本都落在"宇宙 256、覆盖 75%"的退化区，方法无法分离；把它们搬到 n=16 才第一次看清。**我们据此改用对偶排列采样**，不把采样器当贡献——这是架构与实现分离的直接体现。

---

## 五、局地方法可用性的判据（已完成的补充实验）

原计划的"把 CO 补进 Fig. 6"已经做完，36 个窗的结果给出一条可判定规则（**MDA8 主口径**；括号内为 v48，反号计数两度量相同）：

| sector | \|D_local\| 量级 | ratio 偏离 | 局地方法可靠吗 |
|---|---|---|---|
| NOx | 2.1–8.2 ppb (1.0–5.7) | −0.14–1.34 (0.99–2.21)，一处**反号** | ✗ 有 EKMA 拐点 |
| VOC | 2.3–8.5 ppb (1.0–3.6) | ≤0.13 (≤0.11) | 基本可靠（单调） |
| CO | 0.06–0.29 ppb (0.028–0.14) | ≤0.02 (≤0.02) | ✓ 完全可靠，**但因此无关紧要** |

**要写进 Discussion 的一句**：偏离跟随的是**响应是否单调**，而不是效应大小；而"最需要精确的那个部门恰好就是失败的那个"。这句把 ⑥ 从"我们更好"升级成"我们知道什么时候可以不换"——后者审稿人更喜欢，因为它是可证伪的。

数据：`results_local_vs_global.txt`（24 个 3h 窗 + 12 个 6h 块），代码 `_local_vs_global.py`。
