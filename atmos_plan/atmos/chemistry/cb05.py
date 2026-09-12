"""CB05Backend:真实气相化学后端(经 pykpp 解析 CB05,74 物种)。

实现通用 `ChemistryBackend` 契约:状态 = CB05 全物种(分子数/cm³),
`step` 无状态且确定(每段先 resetworld→注入完整状态→设排放→积分)。

**物理环境对齐 Knote et al. (2015) AQMEII-2 理想化箱模型范式**
(环境参数外生给定,不从目标输出反推),三条关键做法:
1. **边界层日循环 PBLH(t)** 由后端自持(updater 逐时更新),并作为稀释/沉降/排放的共同尺度;
2. **稀释 = 边界层抬升的夹卷** `k_dil(t) = (1/h)·(dh/dt) ≥ 0`,夹卷带入自由对流层背景气
   (`EMISSION = X : k·BG_X`),即标准的 `dX/dt = k·(X_bg − X)` —— **不是**手调常数 τ;
3. **干沉降** `X = DUMMY : vd_X / PBLH`,vd 取文献典型白天值,与稀释分离;
   排放按 `h_ref/PBLH` 缩放(表面通量 → 混合比倾向)。

已验证事实(spike 实测):
- CB05 在 pykpp 中可解析积分;单步 6h ≈ 0.43 s。
- 求解器选 **odeint + jac=True**(Knote 式伪反应),`m.run()` 长积分与 `lsoda` 分段会崩/挂。
- **cb05cl 纯气相,无气溶胶/一次 PM** ⇒ PM2.5 投影恒 0(如实,不假造)。
"""
from __future__ import annotations

import contextlib
import io
import os
import tempfile
from typing import Sequence

import numpy as np
import pykpp
from pykpp.mech import Mech

from ..configs import CB05Config
from ..physics import entrainment_dilution, pbl_height, temperature
from .base import ChemistryBackend, StepResult

# CB05 里计入"反应性有机碳(VOC/ROG)"的物种(排除 CH4/CO 等惰性/弱活性)
_ROG = ["C2", "C3", "BTC2", "ETOH", "MEOH", "KETE", "PRPE", "13C4D", "C4D",
        "ISOP", "TERP", "SESQ", "OLE", "IOLE", "TOL", "XYL", "CRES", "BENZENE",
        "PAR", "ALD2", "ALDX", "FORM", "MGLY", "OPEN", "NTR"]

_TEMPLATE = """#INLINE PY_INIT
DT = 1800.
TSTART = 0.
TEND = 3600.
P = {press}
TEMP = {temp}
H2O = h2o_from_rh_and_temp({rh}, TEMP)
StartJday = {jday}
Latitude_Degrees = {lat}
Longitude_Degrees = 0.
HRS = 1.0/3600.
HMIN = {hmin}
HOD0 = {hod0}
HREF = {href}
NDAY = {nday}
HMAX_DAY = {hmax_day}
VENT_TAU_DAY = {vent_tau_day}
PBLH = {href}
KDIL = 0.0
{vds}
EMI_NO2 = {emi_no2}
EMI_NO  = {emi_no}
EMI_SO2 = {emi_so2}
{voc_inits}
{bg_inits}
#ENDINLINE

#INLINE PY_UTIL
add_world_updater(code_updater('''
_HOD = (t/3600. + HOD0) % 24.
_DAY = int(t/86400.)
if _DAY > NDAY-1: _DAY = NDAY-1
if _DAY < 0: _DAY = 0
_HMAX = HMAX_DAY[_DAY]
_VT = VENT_TAU_DAY[_DAY]
_S = max(0.0, sin(pi*(_HOD-6.)/12.))
PBLH = HMIN + (_HMAX-HMIN)*_S**0.5
if _S > 1e-6:
    _DHDT = (_HMAX-HMIN)*(cos(pi*(_HOD-6.)/12.)*pi/12.)/(2.*_S**0.5)
    KDIL = max(0.0, _DHDT/PBLH)/3600. + 1./(_VT*3600.)
else:
    KDIL = 1./(_VT*3600.)
''', incr=300., message='ENV'))
add_world_updater(func_updater(Update_M, incr=300.))
add_world_updater(func_updater(Update_THETA, incr=300.))
#ENDINLINE

#INITVALUES
CFACTOR = P * Avogadro / R / TEMP * centi **3 * nano
ALL_SPEC = 1e-32 * CFACTOR
{initvals}
#EQUATIONS
<EMNO2> EMISSION = NO2 : EMI_NO2 * HRS * CFACTOR * (HREF / PBLH) ;
<EMNO>  EMISSION = NO  : EMI_NO  * HRS * CFACTOR * (HREF / PBLH) ;
<EMSO2> EMISSION = SO2 : EMI_SO2 * HRS * CFACTOR * (HREF / PBLH) ;
{voc_emis}
{dilution}
{inflow}
{deposition}
{cb05_body}
"""


def _load_cb05_body() -> str:
    base = os.path.dirname(pykpp.__file__)
    txt = open(os.path.join(base, "models", "cb05cl.eqn"),
               encoding="utf-8", errors="ignore").read()
    return txt.replace("#EQUATIONS", "", 1)


class CB05Backend(ChemistryBackend):
    action_dim = 4
    species: Sequence[str] = ()

    def __init__(self, config: CB05Config | None = None):
        self.cfg = config or CB05Config()
        self._body = _load_cb05_body()
        # 环境由后端自持(base 契约:h_min/h_max/hod0 经 pbl_height_at 暴露)
        self.h_min = self.cfg.h_min
        self.h_max = self.cfg.h_max
        self.hod0 = self.cfg.hod0
        # 每个进程用独立的生成文件:多进程(根并行)下共用同一路径会互相覆写
        self._model_path = os.path.join(
            tempfile.gettempdir(), f"_cb05_gen_{os.getpid()}.kpp")
        # 情景缩放需在 _build_mech(内部调 _emit_scalars)之前就位
        self._ratio = {"NOx": 1.0, "VOC": 1.0, "SO2": 1.0}
        # 逐日气象(默认 = 常数,σ=0 时等价原固定环境);set_episode() 可重采样
        self._n_days = max(1, int(self.cfg.episode_days))
        self._weather = ([float(self.cfg.h_max)] * self._n_days,
                         [float(self.cfg.vent_tau_h)] * self._n_days)
        es0 = self._emit_scalars(np.zeros(4))
        # 两遍构建:稀释/沉降需按完整物种表展开,故先无环境构建一次换取物种表
        _probe = self._build_mech(es0, with_env=False)
        self.species = tuple(_probe.allspcs)
        self._idx = {s: i for i, s in enumerate(self.species)}
        self._dep_species = [s for s in self.cfg.vd if s in self._idx]
        self.mech = self._build_mech(es0, with_env=True)
        self.cfactor = float(self.mech.start_cfactor)
        self._bg = self.mech.get_y().copy()
        # 收敛诊断:pykpp 的失败信息经 stdout 打印,会被 step() 的 redirect_stdout 吞掉,
        # 故改从 mech.infodict['message'] 取,累计到这里。
        self.n_solver_fail = 0
        self.last_solver_message = ""

    def set_scenario(self, *, nox: float = 1.0, voc: float = 1.0, so2: float = 1.0):
        """设置相对基准的排放比值(RSM/EKMA 扫描);复用已解析机制,零重建成本。"""
        self._ratio = {"NOx": float(nox), "VOC": float(voc), "SO2": float(so2)}
        return self

    def set_episode(self, seed: int) -> "CB05Backend":
        """为当前 episode **采样逐日气象**(通风 τ、日最大 PBL),日循环形状保持不变。

        这是规格"气象从配置分布采样"的落实。日际变化是**跨日择时**有收益的前提:
        若每天相同,第 1/2/3 天可互换,"该在哪天减排"在构造上无收益 ⇒ 环境准静态。
        σ=0 时退化为固定环境(与旧行为一致,保证既有测试不变)。
        """
        c = self.cfg
        rng = np.random.default_rng(seed)
        n = self._n_days

        def _ln(median: float, sigma: float, lo: float, hi: float) -> list[float]:
            if sigma <= 0.0:
                return [float(median)] * n
            x = np.exp(np.log(median) + sigma * rng.standard_normal(n))
            return np.clip(x, lo, hi).tolist()

        hmax = _ln(c.h_max, c.pbl_hmax_logsigma, 500.0, 3000.0)
        vt = _ln(c.vent_tau_h, c.vent_tau_logsigma, 1.5, 12.0)
        self._weather = (hmax, vt)
        return self

    # ---- 契约:自有日循环(base.pbl_height_at 默认实现即用这三个属性) ----
    # h_min / h_max / hod0 已在 __init__ 赋值。

    # ---- 模板渲染 + Mech 构造(抑制 pykpp 的 stdout) ----
    def _emit_scalars(self, action: np.ndarray) -> dict:
        c = self.cfg
        a = np.clip(np.asarray(action, float), 0, 1)
        a_ind, a_tra, a_pow, a_fvoc = a
        # NOx:交通/电力/工业按份额减排,再乘情景比值
        nox = c.nox_base * self._ratio["NOx"] * (1 - a_tra * c.nox_share_tra
                                                 - a_pow * c.nox_share_pow - a_ind * 0.25)
        nox = max(nox, 0.0)
        voc_ef = max(1 - a_ind, 0.0) * max(1 - a_fvoc, 0.0) * self._ratio["VOC"]
        so2 = c.so2_base * self._ratio["SO2"] * max(1 - a_pow, 0.0)
        return dict(nox=nox, voc_ef=voc_ef, so2=so2)

    def _render(self, es: dict, with_env: bool = True) -> str:
        c = self.cfg
        no2 = es["nox"] * 0.8            # NOx 中 NO2 份额
        no = es["nox"] * 0.2
        voc_inits, voc_emis = [], []
        for spc, frac in c.voc_split.items():
            var = f"EMI_{spc}"
            voc_inits.append(f"{var} = {c.voc_base * frac * es['voc_ef']}")
            voc_emis.append(f"<EM{spc}> EMISSION = {spc} : {var} * HRS * CFACTOR * (HREF / PBLH) ;")
        # 稀释(夹卷)作用于全部物种;沉降作用于给定 vd 的物种(with_env=False 时留空)
        dilution = ([f"<DIL{s}> {s} = DUMMY : KDIL ;" for s in getattr(self, "species", ())]
                    if with_env else [])
        deposition = ([f"<DD{s}> {s} = DUMMY : vd_{s} / PBLH ;" for s in self._dep_species]
                      if with_env else [])
        vds = [f"vd_{s} = {v}" for s, v in c.vd.items() if s in getattr(self, "species", ())]
        bg_inits, inflow = [], []
        for spc, ppb in c.bg.items():
            if spc in self._idx_species() and ppb > 0:
                var = f"BGIN_{spc}"
                bg_inits.append(f"{var} = {ppb}")
                inflow.append(f"<IN{spc}> EMISSION = {spc} : KDIL * {var} * CFACTOR ;")
        initvals = [f"{spc} = {ppb}" for spc, ppb in c.bg.items()
                    if spc in self._idx_species()]
        return _TEMPLATE.format(
            press=c.press_Pa, temp=c.temp_K, rh=c.rh, jday=c.start_jday,
            lat=c.lat_deg, hmin=c.h_min, hod0=c.hod0, href=c.h_ref,
            nday=self._n_days,
            hmax_day=repr(list(self._weather[0])),
            vent_tau_day=repr(list(self._weather[1])),
            emi_no2=no2, emi_no=no, emi_so2=es["so2"],
            vds="\n".join(vds),
            voc_inits="\n".join(voc_inits), initvals="\n".join(initvals),
            bg_inits="\n".join(bg_inits), inflow="\n".join(inflow),
            voc_emis="\n".join(voc_emis), dilution="\n".join(dilution),
            deposition="\n".join(deposition),
            cb05_body=self._body,
        )

    def _idx_species(self):
        """已解析的物种集合(首遍构建前为空 ⇒ 背景初值全部保留)。"""
        sp = getattr(self, "species", ())
        return set(sp) if sp else set(self.cfg.bg)

    def _build_mech(self, es: dict, with_env: bool = True) -> Mech:
        open(self._model_path, "w", encoding="utf-8").write(self._render(es, with_env))
        with contextlib.redirect_stdout(io.StringIO()):
            return Mech(self._model_path, verbose=0)

    # ---- 契约实现 ----
    def initial_state(self) -> np.ndarray:
        return self._bg.copy()

    def _atol_for(self, y: np.ndarray):
        """逐物种分层 atol(见 CB05Config.atol_rel)。atol_rel<=0 时退回标量。"""
        if self.cfg.atol_rel <= 0.0:
            return self.cfg.atol
        ref = np.maximum(np.abs(np.asarray(y, float)), 1.0)
        return np.maximum(self.cfg.atol_rel * ref, self.cfg.atol_floor)

    def step(self, conc_in, t0, dt, *, action, met=None) -> StepResult:
        """推演一段。`met` 仅为接口兼容保留 —— CB05 的物理环境由后端自持
        (PBLH/KDIL 由内联 updater 按 t 逐时更新)。"""
        conc_in = np.asarray(conc_in, dtype=float)
        es = self._emit_scalars(np.clip(np.asarray(action, float), 0, 1))
        m = self.mech
        # 无状态:恢复干净 world(含 PY_INIT 默认),覆盖排放标量与逐日气象,注入完整状态
        m.resetworld()
        m.world["HMAX_DAY"] = list(self._weather[0])
        m.world["VENT_TAU_DAY"] = list(self._weather[1])
        for k, v in dict(EMI_NO2=es["nox"] * 0.8, EMI_NO=es["nox"] * 0.2,
                         EMI_SO2=es["so2"]).items():
            m.world[k] = v
        for spc, frac in self.cfg.voc_split.items():
            m.world[f"EMI_{spc}"] = self.cfg.voc_base * frac * es["voc_ef"]
        m.update_world_from_y(conc_in)

        atol = self._atol_for(conc_in)
        n_sub = max(1, self.cfg.substeps_per_control)
        seg = dt / n_sub
        traj = [conc_in.copy()]
        y = conc_in.copy()
        t_hr = t0
        for _ in range(n_sub):
            ts_s0, ts_s1 = t_hr * 3600.0, (t_hr + seg) * 3600.0
            m.world["t"] = ts_s0
            if self.cfg.force_env_update:
                # 强制刷新:PBLH/KDIL 必须是 (y,t) 的函数,不得依赖上一次调用的更新器历史
                m.last_updated = -1e18
            with contextlib.redirect_stdout(io.StringIO()):
                _, Y = m.integrate(ts_s0, ts_s1, y, solver=self.cfg.solver,
                                   jac=self.cfg.jac, rtol=self.cfg.rtol, atol=atol,
                                   mxstep=int(self.cfg.maxstep))
            # 收敛诊断:pykpp 的失败信息打到 stdout(被上面重定向吞掉),改从 infodict 取
            info = getattr(m, "infodict", None)
            msg = info.get("message") if isinstance(info, dict) else None
            if msg and msg != "Integration successful.":
                self.n_solver_fail += 1
                self.last_solver_message = str(msg)
            y = np.clip(Y[-1], 0.0, None)
            traj.append(y.copy())
            t_hr += seg
        arr = np.vstack(traj)
        times = np.linspace(t0, t0 + dt, n_sub + 1)
        return StepResult(conc_out=arr[-1].copy(), conc_traj=arr, times=times)

    def project_pollutants(self, conc) -> dict[str, float]:
        conc = np.asarray(conc, dtype=float)
        idx, cf = self._idx, self.cfactor

        def ppb(name):
            return float(conc[idx[name]] / cf) if name in idx else 0.0

        nox = ppb("NO") + ppb("NO2")
        voc = sum(float(conc[idx[s]] / cf) for s in _ROG if s in idx)
        return {"NOx": nox, "VOC": voc, "O3": ppb("O3"), "SO2": ppb("SO2"),
                "PM25": 0.0}   # cb05cl 无气溶胶模块,如实返回 0
