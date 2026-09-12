"""Chemistry-regime state abstraction for atmospheric MCTS.

Components
----------
RegimeClassifier : phi(s) via FNR = [HCHO]/[NO2]
ActionGrid       : shared coarse discretization of emission-reduction actions
                   (5^4 = 625 cells, matches the M1 scenario matrix)
RegimeStats      : per-regime aggregated (N, Q) over the coarse grid
RegimePrior      : soft action prior derived from the RSM fit (never hard-prunes)
RegimeAwareSelector : PUCT selection with annealed local/regime Q mixing
verification     : oracle-measured aggregation error + threshold sweep

Design rule (project-wide): all approximations live in priors and mixing
weights; Q values are always backed by true box-model rollouts.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import IntEnum
from typing import Callable, Dict, Mapping, Optional, Protocol, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# 1. Regime classification
# ---------------------------------------------------------------------------

class Regime(IntEnum):
    VOC_LIMITED = 0
    TRANSITIONAL = 1
    NOX_LIMITED = 2


N_REGIMES = 3


@dataclass(frozen=True)
class RegimeThresholds:
    """FNR thresholds. Defaults follow the satellite HCHO/NO2 convention;
    treat them as hyperparameters and re-fit via sweep_thresholds()."""
    theta1: float = 1.0   # FNR < theta1  -> VOC-limited
    theta2: float = 2.0   # FNR > theta2  -> NOx-limited


class RegimeClassifier:
    """phi(s): state vector -> regime, from FNR = [HCHO]/[NO2].

    Box-model states carry in-situ concentrations, so the ratio is exact —
    no column-retrieval error. Zero-cost query inside the search loop.
    """

    def __init__(
        self,
        species_index: Mapping[str, int],
        thresholds: RegimeThresholds = RegimeThresholds(),
        eps: float = 1e-12,
    ):
        self.i_hcho = species_index["HCHO"]
        self.i_no2 = species_index["NO2"]
        self.th = thresholds
        self.eps = eps

    def fnr(self, state: np.ndarray) -> float:
        return float(state[self.i_hcho] / max(float(state[self.i_no2]), self.eps))

    def __call__(self, state: np.ndarray) -> Regime:
        r = self.fnr(state)
        if r < self.th.theta1:
            return Regime.VOC_LIMITED
        if r > self.th.theta2:
            return Regime.NOX_LIMITED
        return Regime.TRANSITIONAL


# ---------------------------------------------------------------------------
# 2. Shared coarse action grid
# ---------------------------------------------------------------------------

class ActionGrid:
    """Discretizes continuous reduction ratios into a shared macro-action
    grid so that statistics are comparable across different states/nodes.

    Default: 4 precursors x 5 levels = 625 cells — identical to the M1
    scenario matrix, so M1 outputs can warm-start RegimeStats directly.
    """

    def __init__(self, n_precursors: int = 4,
                 levels: Sequence[float] = (0.0, 0.25, 0.5, 0.75, 1.0)):
        self.levels = np.asarray(levels, dtype=float)
        self.n_precursors = n_precursors
        self.size = len(self.levels) ** n_precursors

    def bin(self, action: np.ndarray) -> int:
        """Map a continuous action vector to its flat macro-cell index."""
        flat = 0
        for a in np.asarray(action, dtype=float):
            i = int(np.argmin(np.abs(self.levels - a)))
            flat = flat * len(self.levels) + i
        return flat

    def cell_center(self, flat: int) -> np.ndarray:
        """Inverse of bin(): representative action of a macro-cell."""
        idx = []
        for _ in range(self.n_precursors):
            flat, i = divmod(flat, len(self.levels))
            idx.append(i)
        return self.levels[list(reversed(idx))]


# ---------------------------------------------------------------------------
# 3. Aggregated per-regime statistics
# ---------------------------------------------------------------------------

class RegimeStats:
    """(N, Q) tables shared by all states within a regime.

    This is the state-abstraction mechanism: a newly expanded node in a
    previously visited regime starts with an informed value estimate
    instead of zero knowledge.
    """

    def __init__(self, grid: ActionGrid, default_q: float = 0.0):
        self.grid = grid
        self.n = np.zeros((N_REGIMES, grid.size), dtype=np.int64)
        self.q_sum = np.zeros((N_REGIMES, grid.size), dtype=float)
        self.default_q = default_q

    def update(self, regime: Regime, action: np.ndarray, value: float) -> None:
        b = self.grid.bin(action)
        self.n[regime, b] += 1
        self.q_sum[regime, b] += value

    def mean(self, regime: Regime, action_bin: int) -> float:
        n = self.n[regime, action_bin]
        if n == 0:
            return self.default_q
        return float(self.q_sum[regime, action_bin] / n)

    def visits(self, regime: Regime, action_bin: int) -> int:
        return int(self.n[regime, action_bin])


# ---------------------------------------------------------------------------
# 4. Regime-aware selection (annealed mixing + soft prior)
# ---------------------------------------------------------------------------

class ChildNode(Protocol):
    """Structural type matching the M2 MCTS node children."""
    action: np.ndarray
    prior: float
    visits: int
    value_sum: float


class RegimePrior:
    """Soft prior over candidate actions, shaped by regime response signs.

    Preferred construction: from the fitted RSM (M1 output) — gradient of
    predicted O3/PM2.5 response w.r.t. each precursor at the regime centroid.
    Fallback: manual sign table (must be validated against M1 data).

    NEVER assigns zero probability: ridge-line misclassification must stay
    recoverable. `floor` guarantees every candidate keeps >= floor prior mass.
    """

    def __init__(self, grid: ActionGrid, beta: float = 2.0, floor: float = 0.02):
        self.grid = grid
        self.beta = beta          # inverse temperature of the prior softmax
        self.floor = floor        # minimum prior mass per candidate
        # logit weight per (regime, precursor dim); + encourages reduction
        self.weights = np.zeros((N_REGIMES, grid.n_precursors))

    @classmethod
    def from_rsm(cls, grid: ActionGrid, rsm_gradient: Callable[[np.ndarray, int], np.ndarray],
                 regime_centroids: Mapping[Regime, np.ndarray],
                 beta: float = 2.0, floor: float = 0.02) -> "RegimePrior":
        """rsm_gradient(state, objective_id) -> dReward/dReduction per precursor."""
        prior = cls(grid, beta=beta, floor=floor)
        for r, s in regime_centroids.items():
            g = rsm_gradient(s, 0)          # objective 0 = primary pollutant
            w = -np.asarray(g, dtype=float)  # reward improves when gradient negative
            norm = np.abs(w).max()
            prior.weights[r] = w / norm if norm > 0 else 0.0
        return prior

    def logits(self, regime: Regime, actions: np.ndarray) -> np.ndarray:
        """Prior logit for each candidate action (rows) under `regime`."""
        return self.beta * actions @ self.weights[regime]

    def distribution(self, regime: Regime, actions: np.ndarray) -> np.ndarray:
        z = self.logits(regime, actions)
        z = z - z.max()
        p = np.exp(z)
        p /= p.sum()
        p = np.maximum(p, self.floor)
        return p / p.sum()


class RegimeAwareSelector:
    """PUCT with annealed local/regime Q mixing.

        Q_hat(s,a) = w * Q_local(s,a) + (1 - w) * Q_regime(phi(s), bin(a))
        w = n / (n + n0)

    Cold start (n << n0): search is guided by cross-state regime knowledge.
    Warm (n >> n0): converges to the node-local, rollout-backed estimate.
    Aggregation error therefore degrades gracefully — worst case is slower
    search, never wrong Q values.
    """

    def __init__(self, stats: RegimeStats, classifier: RegimeClassifier,
                 prior: Optional[RegimePrior] = None,
                 c_puct: float = 1.4, n0: float = 25.0):
        self.stats = stats
        self.classifier = classifier
        self.prior = prior
        self.c_puct = c_puct
        self.n0 = n0

    def mixed_q(self, child: ChildNode, regime: Regime, parent_value: float) -> float:
        n = child.visits
        q_local = child.value_sum / n if n > 0 else parent_value
        q_reg = self.stats.mean(regime, self.stats.grid.bin(child.action))
        w = n / (n + self.n0)
        return w * q_local + (1.0 - w) * q_reg

    def select(self, children: Sequence[ChildNode], state: np.ndarray,
               parent_visits: int, parent_value: float) -> ChildNode:
        regime = self.classifier(state)
        actions = np.stack([c.action for c in children])
        if self.prior is not None:
            priors = self.prior.distribution(regime, actions)
        else:
            priors = np.array([c.prior for c in children], dtype=float)
            priors = priors / priors.sum()
        sqrt_parent = math.sqrt(max(parent_visits, 1))
        best, best_score = None, -np.inf
        for child, p in zip(children, priors):
            q = self.mixed_q(child, regime, parent_value)
            u = self.c_puct * p * sqrt_parent / (1 + child.visits)
            score = q + u
            if score > best_score:
                best, best_score = child, score
        return best

    def backup(self, path: Sequence[Tuple[np.ndarray, ChildNode]], value: float) -> None:
        """Call after each rollout, in addition to the standard MCTS backup:
        feeds every (state, action) on the path into the regime tables."""
        for state, child in path:
            regime = self.classifier(state)
            self.stats.update(regime, child.action, value)


# ---------------------------------------------------------------------------
# 5. Verification: oracle-measured aggregation error + threshold sweep
# ---------------------------------------------------------------------------

def measure_aggregation_error(
    states: np.ndarray,
    q_tables: np.ndarray,
    fnr_values: np.ndarray,
    thresholds: RegimeThresholds,
) -> Dict[str, float]:
    """Empirical aggregation error e per regime.

    states    : (m, n_species) sampled states (from M1 scenario runs)
    q_tables  : (m, grid.size) oracle Q values on the coarse grid,
                obtained by exhaustive rollouts — reuse the M1 scenario matrix
    fnr_values: (m,) precomputed FNR per state (cache-friendly for sweeping)

    e_r = max over action bins of (max_s Q - min_s Q) within regime r.
    This is the quantity Hostetler / Bai & Russell could only bound in
    theory; with the box-model oracle it is directly measurable.
    """
    th = thresholds
    regimes = np.where(fnr_values < th.theta1, Regime.VOC_LIMITED,
                       np.where(fnr_values > th.theta2, Regime.NOX_LIMITED,
                                Regime.TRANSITIONAL))
    out: Dict[str, float] = {}
    total_e, total_w = 0.0, 0
    for r in Regime:
        mask = regimes == r
        if mask.sum() < 2:
            out[f"e_{r.name}"] = float("nan")
            continue
        qt = q_tables[mask]                      # (m_r, grid.size)
        e = float(np.max(qt.max(axis=0) - qt.min(axis=0)))
        out[f"e_{r.name}"] = e
        total_e += e * mask.sum()
        total_w += mask.sum()
    out["e_weighted"] = total_e / max(total_w, 1)
    return out


def sweep_thresholds(
    states: np.ndarray,
    q_tables: np.ndarray,
    species_index: Mapping[str, int],
    candidates_theta1: Sequence[float] = (0.5, 0.75, 1.0, 1.5),
    candidates_theta2: Sequence[float] = (1.5, 2.0, 3.0, 4.0),
) -> RegimeThresholds:
    """Pick (theta1, theta2) minimizing oracle-measured aggregation error.

    Rollouts are done ONCE; the sweep only re-classifies cached FNR values,
    so it costs zero additional simulations.
    """
    probe = RegimeClassifier(species_index)
    fnr = np.array([probe.fnr(s) for s in states])
    best_th, best_e = None, np.inf
    for t1 in candidates_theta1:
        for t2 in candidates_theta2:
            if t1 >= t2:
                continue
            th = RegimeThresholds(t1, t2)
            e = measure_aggregation_error(states, q_tables, fnr, th)["e_weighted"]
            if e < best_e:
                best_e, best_th = e, th
    assert best_th is not None
    return best_th


def regime_sign_consistency(q_tables: np.ndarray, regimes: np.ndarray) -> Dict[str, float]:
    """Assumption-A1 diagnostic: for each regime and action bin, the fraction
    of states whose response sign (vs. the no-control baseline column 0)
    agrees with the regime majority sign. Values near 1.0 validate the
    state abstraction; low values flag a bad threshold or a mixed regime.
    """
    signs = np.sign(q_tables - q_tables[:, [0]])     # (m, grid.size)
    out: Dict[str, float] = {}
    for r in Regime:
        mask = regimes == r
        if mask.sum() < 2:
            out[r.name] = float("nan")
            continue
        s = signs[mask]                              # (m_r, grid.size)
        majority = np.sign(s.sum(axis=0))
        agree = (s == majority).mean(axis=0)
        out[r.name] = float(agree[majority != 0].mean()) if (majority != 0).any() else 1.0
    return out
