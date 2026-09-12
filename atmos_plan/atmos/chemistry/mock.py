"""MockBackend:教学用测试替身(实现通用 ChemistryBackend 契约)。

用途:①上层在无真实机制时毫秒级迭代;②显式编码 VOC-limited 非单调性。
科学结论一律在 CB05 / SAPRC-07 上得出。

非单调性来源:O3 生产对 NOx 饱和(米氏项),O3 损失含随 NOx 线性增长的滴定项
(k_titr_o3·NOx·O3)。高 NOx 区生产饱和、损失主导 ⇒ 减 NOx 反而升 O3。
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from ..configs import ChemParams, EnvConfig
from .base import ChemistryBackend, StepResult

# Mock 的 5 个状态物种(直接即政策污染物)
MOCK_SPECIES: tuple[str, ...] = ("NOx", "VOC", "O3", "SO2", "PM25")
_I = {s: i for i, s in enumerate(MOCK_SPECIES)}


class MockBackend(ChemistryBackend):
    """~5 物种简化机制,刚性 ODE(LSODA),单步 < 1 ms。"""

    species = MOCK_SPECIES

    def __init__(self, config: EnvConfig | None = None):
        self.cfg = config or EnvConfig()
        self.p: ChemParams = self.cfg.chem
        # 环境由后端自持(base 契约)
        self.h_min = self.cfg.h_min
        self.h_max = self.cfg.h_max
        self.hod0 = self.cfg.hod0

    # ---- 动作 → 5 维排放速率 ----
    def _emission_rate(self, action: np.ndarray) -> np.ndarray:
        a = np.clip(np.asarray(action, dtype=float), 0.0, 1.0)
        r_ind, r_tra, r_pow, r_fvoc = a
        e = self.cfg.emis

        def per_species(base):
            return (base[0] * (1 - r_ind) + base[1] * (1 - r_tra) + base[2] * (1 - r_pow))

        rate = np.zeros(self.n_species)
        rate[_I["NOx"]] = per_species(e.nox)
        rate[_I["VOC"]] = per_species(e.voc) * (1 - r_fvoc)
        rate[_I["SO2"]] = per_species(e.so2)
        rate[_I["PM25"]] = per_species(e.pm)
        return rate

    def _rhs(self, t, y, emis_rate, h_pbl_fn):
        p = self.p
        nox, voc, o3, so2, pm = y
        h = max(50.0, h_pbl_fn(t))
        k_dil = p.v_dil * 3600.0 / h
        prod = (p.prod_max * nox / (p.k_nox_o3 + nox) * voc / (p.k_voc_o3 + voc))
        loss_o3 = p.k_dep_o3 + p.k_titr_o3 * nox + k_dil
        return np.array([
            emis_rate[0] - (p.k_nox + k_dil) * nox,
            emis_rate[1] - (p.k_voc + k_dil) * voc,
            prod - loss_o3 * o3,
            emis_rate[3] - (p.k_so2 + k_dil) * so2,
            emis_rate[4] + p.pm_from_so2 * so2 - (p.k_pm + k_dil) * pm,
        ])

    # ---- 契约实现 ----
    def initial_state(self) -> np.ndarray:
        return np.array([self.cfg.init_conc[s] for s in MOCK_SPECIES], dtype=float)

    def step(self, conc_in, t0, dt, *, action, met) -> StepResult:
        conc_in = np.asarray(conc_in, dtype=float)
        emis_rate = self._emission_rate(action)
        h_pbl_fn = met["h_pbl_fn"]
        t_eval = np.linspace(t0, t0 + dt, 13)
        sol = solve_ivp(
            lambda t, y: self._rhs(t, y, emis_rate, h_pbl_fn),
            (t0, t0 + dt), conc_in, method="LSODA",
            t_eval=t_eval, rtol=1e-6, atol=1e-6,
        )
        if not sol.success:
            raise RuntimeError(f"MockBackend 积分失败: {sol.message}")
        traj = np.clip(sol.y.T, 0.0, None)
        return StepResult(conc_out=traj[-1].copy(), conc_traj=traj, times=sol.t.copy())

    def project_pollutants(self, conc) -> dict[str, float]:
        conc = np.asarray(conc, dtype=float)
        return {s: float(conc[_I[s]]) for s in MOCK_SPECIES}
