# v1.0.0 — finite-intervention ozone attribution (paper submission snapshot)

Code and data companion for the manuscript *"Finite interventions, full
spectra: where local ozone attribution fails, and what measuring the whole
interaction spectrum costs"*.

This is the **frozen snapshot** referenced by the paper's code-availability
statement. Numbers in the manuscript trace to the `results_*.txt` files in
this tag; the mapping is tabulated in `README.md` and `paper_outline.md`.

## What is in this release

- **The box-model engine** (`_kppwork/kpp_box.py` + KPP/CB05 build via
  `build_kpp7.sh`): config -> namelist -> 72 h integration -> trajectory.
  Exact switch-off (window emission multipliers 1 -> 0), dummy-mechanism
  variant for the L1 gate.
- **Exact oracle tables** (the paper's ground truth, full enumeration):
  - `_kppwork/n16_v_r_sat.csv`, `_kppwork/n16_v_r_lim1.csv` — n=16, d3_mean
  - `_kppwork/n16_v_r_sat_m8.csv`, `_kppwork/n16_v_r_lim1_m8.csv` — n=16, MDA8
  - `_kppwork/co_v.csv` (n=12), `sham_v.csv`, `scen_v.csv`, `scale_v.csv`,
    `hddm_v_*.csv`, `dp_v_*.csv`
- **Campaign runners** (`sweep_*.sh`, `_sweep_dp.sh`): sharded, crash-resumable.
- **Analysis suite** (`exp_*.py`, `_*.py` under `_kppwork/`): Möbius spectrum
  and truncation curves, sampling benchmarks, local-vs-global protocol,
  Taylor-2 comparison, self-diagnosis and calibration, verification gates.
- **Observation-anchoring pipeline** (root): AQS/PAMS pullers, Open-Meteo
  meteorology (`met_*.csv`), PAMS->CB05 species routing (`_voc_map.py`),
  CB05 observation composites (`obs_*_cb05.csv`), the anchoring fit
  (`_fit_anchor.py`) and its results (`results_fit_*.txt`).
- **Documents**: `README.md` (repo map), `paper_outline.md` (claim/exhibit
  map with the numbers), `methods_zh.md` (methods draft),
  `vibecoding_prompt_atmos_mcts.md` (working log: decisions, bugs, verdicts,
  retractions).

## What is deliberately not in this release

- Raw AQS bulk downloads (~160 MB of site metadata plus per-site pulls).
  They are re-fetchable with the versioned pullers (`_aqs_api.py`,
  `_aqs_pull_par.py`); the compact derived form that *is* included is
  `obs_*_cb05.csv` + `met_*.csv`.
- Per-worker sweep shards (thousands of `_v_<k>.csv` partials). They are
  merged into the single oracle tables listed above.
- The LaTeX manuscript source (kept in a separate private working tree).
- Interactive-debugging scratch output (hundreds of one-shot probes); the
  reusable ones are listed in `README.md`.

## Scope notes (read before citing)

- The host is a **0-D single-cell box**; the framework is host-agnostic and
  the nocturnal boundary-layer regime is explicitly out of scope
  (`results_kdil.txt` gives the computable validity gate).
- The working sampler is **standard antithetic permutation sampling**. A
  structure-aware active sampler was developed as an arm, **failed its
  pre-registered 2x criterion at n=8** (48 vs 16 runs at 10% error) and was
  discarded; this is reported rather than hidden.
- The full spectrum is obtained **by enumeration only** (n <= 16). Sampling
  estimates the attribution functionals, not the spectrum; no claim of
  sub-2^n spectrum recovery is made.
- **Metric labelling (2026-09 audit).** The originally released sampling
  outputs (`scaling_n16.txt`, `attrib_bench.txt`, `bench_n16.txt`) were
  computed on the d3_mean column. `_dual_bench.py` recomputes every sampling
  exhibit on both metrics with zero chemistry runs
  (`scaling_dual_metrics.txt`); the pooled "1.6-1.8x per doubling" and
  "same accuracy costs 2.5-3x" statements are **retired**. In MDA8 the same
  estimator needs an order of magnitude more runs and a 1.5% single-run
  error is unreachable within the tested budgets.
- **Anchoring (2026-09 audit + 2026-09-30 re-run).** The anchoring cost
  function includes O3 with weight 1.0 (`TARGETS` in `_fit_anchor.py`), so
  it is a fit, not a held-out verification. Both city arms are now
  reproducible under the 8-knob protocol: Deer Park 2017s (cost 0.399;
  limited side: chord ratios 0.98-1.05, 0/24 sign flips) and Dallas 2016w
  (cost 0.441; saturated side: 0/24 sign flips but NOx chord ratios
  0.82-5.10 -- the city-side amplitude failure; VOC/CO stay ~1). The
  previously quoted Dallas "signature" (7/8, 1/24 flips) had no supporting
  run and is withdrawn; the pre-stated expectation that the saturated
  slice would flip signs somewhere did not materialise, and is disclosed
  as such.
- **Spin-up sufficiency (2026-09-30).** Resolved: at the Dallas anchored
  slice the 3-day and the 10-day run agree on the last-day composite to
  <=0.4% daily drift (`results_spinup_dallas.txt`), so the scored 3-day
  window is a diurnal cycle. The earlier caveat is retired.
- **Metric conventions (2026-09-30).** Measured, not assumed: the 289- vs
  288-point daily-mean choice moves Shapley values by at most 0.055 ppb
  (0.33% of max|phi|) and coalition differences by <0.12 ppb
  (`results_d3_conv.txt`, full n=8 game re-derived under both
  conventions); the MDA8 window-containment check
  (`results_mda8_wincheck.txt`, 5 configurations) is bit-identical
  between the t>=48 and t>=56 rules.

## Manuscript reference

Submitted to *Atmospheric Chemistry and Physics* (in preparation). A DOI
will be minted for this tag on acceptance via the journal's archive
recommendation.

## Changes since the previous tag

Initial release; no previous tag.