# intervention-attribution

Code and data companion for the manuscript *"Finite interventions, full
spectra: where local ozone attribution fails, and what measuring the
whole interaction spectrum costs"* (in preparation).

Sector-level ozone control is posed as a **finite cooperative game**:
players are (sector x time-window) binary emission switches, the value
function is the MDA8 response of a KPP-compiled CB05 box model, and the
object of study is the full **Möbius (Harsanyi dividend) interaction
spectrum**. Exact enumerations at n = 8 / 12 / 16 players provide the
ground truth for every claim; two observationally anchored city-seasons
(Dallas 2016 winter, Deer Park 2017 summer) provide the external checks.

Three things this repo is *not* claiming: sampling recovery of the full
spectrum from fewer than 2^n evaluations (the spectrum comes from
enumeration only); a new sampler (the working sampler is standard
antithetic permutation sampling — a structure-aware active sampler was
developed, failed its pre-registered criterion, and is reported as
such); anything beyond a 0-D host (the framework is host-agnostic; the
nocturnal boundary-layer regime is explicitly out of scope).

## Repository layout

| path | what it is |
|---|---|
| `_kppwork/` | the experiment engine: KPP/CB05 build, box-model driver, campaign sweeps, spectrum/benchmark analyses, oracle tables |
| `aqs_sites_data/`, `data_aqs/` | raw AQS/PAMS downloads (site metadata, hourly data) |
| `atmos_plan/` | legacy first-prototype tree, superseded by `_kppwork/` (kept for provenance) |
| `met_*.csv` | hourly meteorology per city-season (Open-Meteo): T2m, RH, pressure, wind, shortwave, boundary-layer height |
| `obs_*_cb05.csv` | observation composites routed to CB05 species (PAMS speciated VOC + NOx/O3) per city-season |
| `results_*.txt` | one file per paper claim (map below) |
| `paper_outline.md` | figure/table plan with the numbers behind every claim |
| `RELEASE_NOTES.md` | release body for the submission snapshot (what is / is not in the archive, scope caveats). No licence is attached yet; add one before publishing a release. |
| `methods_zh.md` | methods draft (Chinese working version) |
| `vibecoding_prompt_atmos_mcts.md` | the original project prompt (early SAPRC-07/MCTS planning notes; historical, superseded by the CB05 box) |
| everything `_*`-prefixed | scratch: probes, working logs, one-off diagnostics. Kept for provenance; not part of any pipeline. |

## The pipeline in one page

```
observations                game generation            campaigns                 analyses
-------------               -----------------           ----------------          ---------------
AQS/PAMS (_aqs_*.py)        gen_*.py                   sweep_*.sh (12 workers)  exp_*.py / _*.py
  -> obs_*_cb05.csv  ----->  namelists per   ------->   v_*.csv oracle  ------->  results_*.txt
Open-Meteo (_met_pull.py)     coalition mask             tables (2^n rows)         (paper claims)
  -> met_*.csv              gen via kpp_box.py          resumable, sharded        methods_zh.md
anchors (_fit_anchor.py)     cfg_of(mask)                                         paper_outline.md
  -> results_fit_*.txt
```

Every campaign follows the same pattern: a `gen_*.py` writes one
namelist per coalition (the switch-off is exact — window emission
multipliers 1 -> 0), a `sweep_*.sh` runs them on sharded KPP build
directories with crash-resume, and an analysis script reads the merged
value table. The currency of all benchmarks is **distinct coalition
evaluations** (= chemistry runs).

## Root scripts

| script | role |
|---|---|
| `_aqs_*.py` | AQS/PAMS acquisition suite: site screening, hourly pulls, PAMS live probes, parameter metadata (`results_aqs_pick.txt`, `results_sites_final.txt` document the site choice) |
| `_met_pull.py` | Open-Meteo hourly meteorology per city-season -> `met_*.csv` |
| `_blh_to_kdil.py` | boundary-layer height -> box dilution timescale (tau_BLH), the observation-bounded constraint used by anchoring |
| `_voc_map.py` | PAMS species -> CB05 model species routing (audit trail: `results_voc_map.txt`, `results_aqs_voc_panel.txt`) |
| `_voc_nc.py` | VOC netCDF ingest helper (early phase) |
| `_obs_build.py`, `_obs_cb05.py` | observation composites -> `obs_*_cb05.csv` (24 h composites per city-season) |
| `_fit_anchor.py` | city anchoring: 8 physical knobs fitted to observed composites with observation-bounded ranges (O3 and NOx carry weight 1.0 in the cost -- a fit, not a held-out check) -> `results_fit_*.txt` |
| `_scan_anchor.py`, `_diag_fit.py` | anchor diagnostics; `_diag_fit.py` is the spin-up sufficiency check (3-day vs 10-day last-day composites) -> `results_spinup_dallas.txt` |
| `_local_vs_global.py` | local-vs-global protocol driver: local central differences (delta=0.05) vs measured finite chords, both metrics -> `results_local_vs_global.txt` |
| `regime_mcts.py` | early planning prototype (MCTS over coalitions; the negative result that seeded the planning benchmark) |

## Core files in `_kppwork/`

| file | role |
|---|---|
| `kpp_box.py` | **the box-model driver**: config dict -> namelist -> KPP run -> trajectory; scenario defaults (r_sat / r_lim1); all switch mechanics |
| `build_kpp7.sh`, `cb05.kpp`, `cb05_Main.f90` | KPP/CB05 mechanism, build scripts, generated solver |
| `gen_n16.py`, `gen_co_sweep.py`, `gen_hddm.py`, `gen_hddm_dp.py`, `gen_flipscan.py`, `gen_sham.py`, `gen_scale_masks.py`, `gen_scen_sweep.py`, `mk_dummy.py` | namelist generators, one per experiment |
| `sweep_n16.sh`, `sweep_co*.sh`, `sweep_hddm.sh`, `sweep_flips.sh`, `sweep_sham.sh`, `sweep_scale.sh`, `sweep_mobius.sh`, `sweep_scen.sh`, `_sweep_dp.sh` | parallel campaign runners (sharded builds, crash-resume) |
| `attrib_oracle.py`, `attrib_samplers.py` | exact attribution from enumeration; the sampling arms (perm, perm_anti, perm_strat, kernelshap, OURS) with the budget currency |
| `exp_n16.py`, `exp_n16_bench.py` | n=16 game driver and benchmark rerun on the exact table |
| `_n16_spectrum.py`, `_m8_spectrum.py` | Möbius spectrum, order masses, truncation-error curves (both metrics) |
| `_local_vs_global` analyses: `_zero_effect.py`, `_ridge_delta.py`, `_flip_report.py` | three failure-mode statistics: near-zero chords on the ridge, delta-refinement Richardson check, ridge-safe flip report |
| `_hddm_check.py`, `_hddm_coalitions.py`, `_dp_taylor.py`, `_deerpark_lvg.py` | Taylor-2 (HDDM) single-axis and coalition scoring; Deer Park anchored-slice versions |
| `_dual_bench.py` | dual-metric (d3_mean + MDA8) re-derivation of every sampling exhibit, zero chemistry runs: builds the n=8/12/16 oracles for both metrics (n=8 MDA8 read off the n=16 table at union-of-pairs masks, licensed by a bit-exact v_ppb cross-check) and reruns the cost law and the five attribution arms |
| `exp_attrib.py`, `exp_calib.py`, `exp_selfdiag.py`, `exp_planning_bench.py`, `exp_quota.py`, `exp_sham.py`, `exp_co12*.py`, `exp_scale.py`, `exp_mda8.py`, `exp_stoch_gap.py`, `exp_cost_frontier.py`, `exp_surrogate_opt.py` | the analysis suite (one experiment per paper claim) |
| `v_*.csv` oracle tables | `n16_v_r_*.csv` / `n16_v_r_*_m8.csv` (n=16, both metrics), `co_v.csv` (n=12), `sham_v.csv`, `scen_v.csv`, `scale_v.csv`, `hddm_v_*.csv`, `dp_v_*.csv`, `flip_v_*.csv` |
| `_mda8_wincheck.py`, `_grid_check.sh` | metric-convention checks (window containment -> `results_mda8_wincheck.txt`, 5-min grid facts) |
| `_d3_conv.py` | 289- vs 288-point daily-mean convention, re-derived on the full n=8 game (256 runs, both conventions from the same trajectories) -> `results_d3_conv.txt` |
| early-phase numerics gates | `exp_g4_kpp.py`, `g3_check.py`, `ab_compare.py`, `v1_regress.py`, `probe_ridge*.py`, `probe_isop.py` |

## Results files -> paper claims

| results file | claim / exhibit |
|---|---|
| `_kppwork/sham_result.txt` | Table 1: L1 dummy-player gate (dividends ~1e-14) |
| `_kppwork/spectrum_orders*.txt` | Fig. 2 + Table 6: truncation divergence (MDA8, n=16) |
| `_kppwork/scaling_dual_metrics.txt` (from `_dual_bench.py`) | Fig. 4 + Table 5: sampling cost law, **both metrics** — the primary reference (the older `scaling_n16.txt` / `attrib_bench.txt` / `bench_n16.txt` are d3_mean-only; the pooled "1.6-1.8x per doubling" and "same accuracy 2.5-3x" statements are retired) |
| `results_local_vs_global.txt` | Fig. 5/6 + Table 2: local vs finite effects, three failure modes |
| `results_scan_flips.txt` | Table 8: failure map over the 25-point EKMA grid (8/25 cells flip) |
| `results_zero_effect.txt`, `results_ridge_delta.txt` | ridge zero-effect windows, delta-refinement verification |
| `results_co12.txt` (from `_kppwork/exp_co12.py`) | Table 3: player-set robustness (8 vs 12 players; max per-player shift 0.10 ppb) |
| `results_quota.txt` (from `_kppwork/exp_quota.py`) | Fig. 7: budget-constrained recipes flip across regimes |
| `figs/` (from `_kppwork/make_figs.py`) | working-draft figures Fig 2-7, regenerated from the released result files (index in `paper_outline.md`); Fig 1 is a schematic, pending |
| `results_fit_deerpark_2017s.txt`, `results_deerpark_lvg.txt`, `results_obs_*_cb05.csv`, `_kppwork/_deerpark_lvg.py` | Table 9 + L4: Deer Park anchored-slice out-of-sample test (limited side): all chord ratios 0.98-1.05, 0/24 sign flips -- "the failure disappears" reproduced |
| `results_fit_dallas_2016w.txt`, `results_dallas_lvg.txt`, `_kppwork/_dallas_lvg.py` | Table 9 + L4: Dallas anchored-slice test (saturated side, 8-knob re-anchor cost 0.441): 0/24 sign flips, but NOx chord ratios 0.82-5.10 (4/8 windows >= 1.9x; VOC/CO ~1) -- the city-side amplitude failure. The pre-stated "saturated -> flips" expectation did not materialise; disclosed, not hidden |
| `results_spinup_dallas.txt` | L4 + Methods: spin-up sufficiency (3-day vs 10-day last-day drift <= 0.4%) -- the 3-day window is a diurnal cycle |
| `results_surrogate.txt` (from `_kppwork/exp_surrogate.py`) | reduced-representation test on the exact oracles: GBM surrogate + exact (TreeSHAP-equivalent) Shapley vs perm_anti at equal budget -- the surrogate loses everywhere except a narrow n=16 MDA8 mid-budget band, and collapses to 100% error below B<=32 at n=8; fit-then-explain is not a substitute for unbiased sampling |
| `results_co_share.txt`, `results_cum_share.txt`, `results_top_dividends.txt` | post-hoc spectrum diagnostics: CO-touching mass share f_CO(k) (5.9% total vs 33% of slots), cumulative mass share R_q, named top dividends behind the order-k peak |
| `_kppwork/results_mda8_wincheck.txt` | Methods: MDA8 window containment (windows ending at t>=48 vs t>=56 bit-identical, 5 configurations) |
| `_kppwork/results_d3_conv.txt` (+ `d3_conv_n8.csv`) | Methods: 289- vs 288-point daily-mean convention measured end to end (Shapley shifts <= 0.055 ppb = 0.33% of max|phi|) |
| `results_hddm_check.txt`, `results_hddm_coalitions.txt`, `results_deerpark_taylor.txt` | Table 7d: Taylor-2 tail failure on the saturated side only |
| `_kppwork/selfdiag*.txt`, `selfdiag_calib.txt` | Table 7a/7b: self-diagnosis flags; holdout calibration (e ~ 0.33 h) |
| `results_kdil.txt` | Table 4: nocturnal validity diagnostics (scope gate) |

## Environment

Windows host + WSL2 (Ubuntu). The KPP-compiled box builds live under
`/root/kpp7*` inside WSL (built by `build_kpp7.sh`); campaign workers are
sharded copies of that build. Python 3 with numpy/scipy; no other
dependencies. Trajectories: 72 h integration, 5-minute output, first
48 h equilibration; MDA8 by right-endpoint rectangle rule (96 samples =
exactly 8.0 h); see `_kppwork/_grid_check.sh`.

## Reproducing

The repo is organised for *auditability* rather than one-command
reproduction: every paper number traces to a results file, every
results file to a script + oracle table, every oracle table to a
campaign sweep over the boxed CB05 host. Start from
`paper_outline.md` (the claim map), descend into the results file of
interest, then find its driver by name.

## Cite

Placeholder — a tagged release with DOI will accompany the submission.
