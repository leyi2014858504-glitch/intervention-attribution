"""数据类配置。所有可调参数集中在此,实验时快照进 results/。"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Dict


# 物种顺序(全项目统一索引)
SPECIES: tuple[str, ...] = ("NOx", "VOC", "O3", "SO2", "PM25")
_IDX = {s: i for i, s in enumerate(SPECIES)}


@dataclass(frozen=True)
class ChemParams:
    """Mock 机制的速率/系数(单位:浓度 ppb、时间 h;PM25 视作 µg/m³ 量纲)。

    刻意编码的**非单调性**:O3 生产对 NOx 饱和、损失随 NOx 线性增长(滴定),
    因此在高 NOx(VOC-limited)区,减 NOx 反而升 O3。
    """
    # 一级清除(/h)
    k_nox: float = 0.25
    k_voc: float = 0.15
    k_so2: float = 0.20
    k_pm: float = 0.30
    # O3 生产/损失
    prod_max: float = 40.0        # ppb/h,O3 最大生产速率
    k_nox_o3: float = 15.0        # 生产对 NOx 半饱和
    k_voc_o3: float = 250.0       # 生产对 VOC 半饱和
    k_dep_o3: float = 0.04        # O3 干沉降/基线损失 /h
    k_titr_o3: float = 0.0015     # NOx 滴定 O3 系数 /(h·ppb)
    # 边界层稀释:k_dil = v_dil * 3600 / h
    v_dil: float = 0.005          # m/s 沉降/通风速度
    # PM2.5 二次生成
    pm_from_so2: float = 0.15     # 每 ppb SO2 → µg/m³ PM2.5 的转化系数(/h 已并入一级)


@dataclass(frozen=True)
class EmissionConfig:
    """基准排放:三部门 × 物种(ppb/h 或 µg/m³·h),供 MockBackend 翻译动作。

    物理剖面(贴合城市源):交通=NOx高/VOC中/SO2≈0;电力=SO2主导/NOx中/VOC≈0;
    工业=VOC高/NOx中/SO2中。这样"减电力 NOx 且保持 VOC"是干净的 VOC-limited 探针。
    """
    # [industry, traffic, power] 各部门对每种物种的基准注入
    nox: tuple[float, float, float] = (12.0, 25.0, 15.0)
    voc: tuple[float, float, float] = (20.0, 15.0, 0.0)
    so2: tuple[float, float, float] = (6.0, 1.0, 30.0)
    pm: tuple[float, float, float] = (2.0, 1.0, 3.0)


@dataclass(frozen=True)
class CB05Config:
    """CB05 后端的框模型参数(ppb 尺度;经 CFACTOR 转分子数/秒积分)。

    说明:cb05cl 为纯气相机制,**无气溶胶/一次PM**,故 PM2.5 投影恒为 0。
    基准排放/背景以 ppb·h⁻¹ 给,后端内部换算为 pykpp 所需的分子数尺度。
    """
    # 环境
    temp_K: float = 298.15
    press_Pa: float = 101325.0
    rh: float = 50.0
    lat_deg: float = 40.0
    start_jday: int = 196          # 7/15 附近(盛夏)
    # 边界层日循环(Knote 范式:环境外生给定,后端自持)
    h_min: float = 250.0           # m,夜间边界层
    h_max: float = 1800.0          # m,午后边界层
    hod0: float = 9.0              # t=0 时的当地时(小时)
    h_ref: float = 1000.0          # 排放基准高度:排放按 h_ref/PBLH 缩放(ppb/h 定义在 1000 m)
    # 水平通风时间尺度 τ=L/u(L≈40 km 城区尺度,u≈3 m/s 风速)。**必须存在**:
    # 仅有夹卷(白天边界层增长)时夜间 18 h 无稀释 → NOx 累积到数百 ppb、O3 被滴定归零。
    vent_tau_h: float = 4.0
    # --- 逐日气象采样(打破"强迫完全周期"这一准静态根源)---
    # 说明:规格原本即要求"气象从配置分布采样"。此前气象是常数 ⇒ 每天相同 ⇒
    # 第 1/2/3 天可互换,"该在哪天行动"在构造上无收益。σ=0 时退化为固定 Knote 环境。
    vent_tau_logsigma: float = 0.0   # 通风 τ 的对数正态 σ(0 = 无日际变化)
    pbl_hmax_logsigma: float = 0.0   # 日最大 PBL 的对数正态 σ(0 = 无日际变化)
    episode_days: int = 3            # 一个 episode 覆盖的天数(72 h ⇒ 3 天)
    # 干沉降速度 vd [m/s](典型白天文献量级 1e-4–3e-2;与稀释分离,不再混为一谈)
    vd: Dict[str, float] = field(default_factory=lambda: {
        "O3": 0.005, "NO2": 0.0015, "NO": 0.0002, "HNO3": 0.030,
        "H2O2": 0.010, "SO2": 0.005, "CO": 0.0003, "FORM": 0.005,
        "ALD2": 0.004, "ALDX": 0.004, "PAN": 0.002, "HONO": 0.005,
        "ISOP": 0.003, "MGLY": 0.003, "NTR": 0.002,
    })
    # 基准排放(ppb/h),按部门拆分用于翻译动作;NOx 走交通/电力,VOC 走工业+交通
    nox_base: float = 4.0          # ppb/h@1000m(标定:NOx≈34 ppb,O3 日峰≈96 ppb,VOC-limited)
    nox_share_tra: float = 0.6     # 交通在 NOx 中的份额(用于 a_tra 减排)
    nox_share_pow: float = 0.15
    # VOC 基准(ppb/h)分散到几个活性烃;fvoc 控制额外乘 (1-a_fvoc)
    voc_base: float = 10.0
    voc_split: Dict[str, float] = field(default_factory=lambda: {
        "ISOP": 0.15, "TOL": 0.25, "XYL": 0.20, "OLE": 0.15, "PAR": 0.25,
    })
    so2_base: float = 4.0          # 电力为主 → a_pow 减排
    # 背景(ppb):既作初值,也作**入流回补**目标(BG_O3 是 O3 不塌零的保障)
    bg: Dict[str, float] = field(default_factory=lambda: {
        "O3": 30.0, "NO": 2.0, "NO2": 20.0, "CO": 150.0, "CH4": 1800.0, "SO2": 2.0,
    })
    # 积分控制(诊断得:odeint 快、lsoda 挂;分子尺度 atol)
    solver: str = "odeint"
    rtol: float = 1e-4
    atol: float = 1e6               # 标量 atol(atol_rel<=0 时生效;否则被分层向量取代)
    # --- 数值保真(收敛)---
    #: 解析 Jacobian 开关。True 用机制自带 Dfun;若 Dfun 定向/转置错误,LSODA 会报
    #: "Excess work done ... (perhaps wrong Dfun type)"。诊断时置 False 用数值 Jacobian。
    jac: bool = True
    #: **逐物种分层 atol**(解决"单一标量 atol 横跨 ~10 个量级"的问题):
    #: `atol_i = max(atol_floor, atol_rel · |y_i|)`。前体(ppb 级)由 rtol 主导;
    #: 痕量物种由 atol_rel/atol_floor 兜住,既不被放松到失控、也不因近零而被卡死。
    #: `atol_rel<=0` ⇒ 退回标量 `atol`(旧行为)。
    atol_rel: float = 0.0
    atol_floor: float = 1e4         # 分子/cm³,痕量物种绝对下限
    #: LSODA 内部步数上限(pykpp 默认 1000)。**若 "Excess work done" 频繁出现,
    #: 说明求解器在"放弃前"就打满了步数 ⇒ 结果并非容差收敛,只是尽力而为。**
    maxstep: int = 1000
    #: 诊断开关:每个积分段开始前**强制刷新环境更新器**(把 `mech.last_updated` 置为 −∞)。
    #: 缘由:pykpp 的 `resetworld()` **不重置 `last_updated`**,导致 `Update_World` 的触发
    #: (`t - last_updated > incr`)依赖**上一次调用**的历史 ⇒ 环境参数(PBLH/KDIL)不是
    #: (y,t) 的函数 ⇒ 结果随调用顺序漂移。置 True 可验证/规避该缺陷。
    force_env_update: bool = False
    substeps_per_control: int = 2  # 每个 6h 控制步的积分段数(轨迹分辨率)


@dataclass(frozen=True)
class EnvConfig:
    """MDP 环境配置。"""
    dt_hours: float = 6.0
    horizon_steps: int = 12                 # 12 × 6h = 72h
    pm25_limit: float = 35.0                # µg/m³(仅当 use_pm_penalty=True 生效)
    o3_limit: float = 160.0                 # µg/m³(以 ppb 量纲近似)
    lam_o3: float = 1.0                     # O3 超标权重
    #: 奖励模式。`exceedance` = AQI 超标惩罚(政策口径,但阈值设不好会**饱和**,
    #: 使所有方法同分、指标失去区分度);`exposure` = 连续暴露 −mean(O3),
    #: 无阈值、天然有动态范围,用于"搜索器/时序结构"这类需要精细排序的实验。
    reward_mode: str = "exceedance"          # exceedance | exposure
    #: **减排努力成本** [ppb 每单位平均减排比例]:`r = −(暴露 + λ·mean(action))`。
    #: λ=0 退化为"免费减排"(此时最优动作恒为"每天尽量减" ⇒ 环境准静态,
    #: 时序组合无价值)。λ>0 使"行动有代价" ⇒ 必须**择时**,才可能造出时序结构。
    effort_lambda: float = 0.0
    use_pm_penalty: bool = False            # 纯 CB05 气相无气溶胶 ⇒ PM2.5 恒 0,默认关闭
    # 注:补回 PM 项需换含气溶胶机制(CB6r5+ae7/SAPRC+SOA),当前锁定 CB05 单目标 O3。
    # 边界层日循环(后端自持;经 ChemistryBackend.pbl_height_at 取值)
    h_min: float = 300.0                    # m
    h_max: float = 2000.0                   # m
    hod0: float = 0.0                       # t=0 时的当地时(小时)
    # 初始浓度(ppb / µg·m⁻³)——reset 时加噪声
    init_conc: Dict[str, float] = field(default_factory=lambda: {
        "NOx": 40.0, "VOC": 150.0, "O3": 60.0, "SO2": 40.0, "PM25": 30.0,
    })
    init_jitter: float = 0.15               # 相对采样噪声
    chem: ChemParams = field(default_factory=ChemParams)
    emis: EmissionConfig = field(default_factory=EmissionConfig)

    def to_dict(self) -> dict:
        return asdict(self)
