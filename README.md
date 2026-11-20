# intervention-attribution

Code and data companion for the manuscript **"Exact response tables reveal
high-order interactions in time-resolved ozone control"** (submitted to
*Advances in Atmospheric Sciences*).

Sector-level ozone control is posed as a **finite cooperative game**. Players
are (precursor sector x 3-h time-window) binary emission switches, the value
function is the MDA8 change of a KPP-compiled CB05 box model, and the object of
study is the full **Möbius (Harsanyi dividend) interaction spectrum**. Exhaustive
enumeration of all 2^16 coalitions in four n = 16 games — two idealized synthetic
scenarios (`r_sat`, `r_lim1`) and two observation-anchored city slices (Dallas
2016 winter, Deer Park 2017 summer) — provides the exact reference against which
low-order truncation, local derivatives (T1/T2, the DDM-1/HDDM-2 equivalents) and
greedy budget allocation are judged.

What this repository does **not** claim: sampling recovery of the spectrum from
fewer than 2^n runs (the spectrum comes from enumeration only), a new sampler, or
a host beyond the 0-D box.

## Repository layout

| path | what it is |
|---|---|
| `_kppwork/` | experiment engine: KPP/CB05 build, box-model driver, campaign sweeps, exact oracle tables, diagnostics and figure drivers |
| `figs/` | the eight manuscript figures, as PNG and PDF |
| `met_*.csv` | hourly meteorology per city-season (Open-Meteo): T2m, RH, pressure, wind, shortwave, boundary-layer height |
| `obs_*_cb05.csv` | observation composites routed to CB05 species (PAMS speciated VOC + NOx/O3) per city-season |
| `results_*.txt` | one file per reported diagnostic, written by the scripts below |
| `_*`-prefixed files | scratch probes and one-off diagnostics; not part of the pipeline |

## Exact tables included

| file | content |
|---|---|
| `_kppwork/n16_vh_r_sat_m8.csv`, `..._r_lim1_m8.csv`, `..._dallas_2016w_m8.csv`, `..._deerpark_2017s_m8.csv` | the four n = 16 MDA8 oracles: 65,536 rows each, one row per coalition. Columns `name` (mask `m00000`–`m65535`; bits 0–7 = NOx windows 00–03 … 21–24, bits 8–15 = VOC windows), `v_ppb` (daily-mean surface O3, ppb) and `v_mda8` (MDA8, ppb), both absolute — game values are obtained by subtracting the all-on row `m00000` — plus `nfail`, `rc` |
| `_kppwork/co_v_m8_*.csv` | n = 12 (NOx, VOC, CO at 6-h resolution) MDA8 shards behind the spectral ladder |
| `_kppwork/n16_v_*.csv` | earlier n = 8 / n = 16 synthetic tables (d3_mean and MDA8 variants) |
| `_kppwork/hddm_h_v_*.csv` | 512 central-difference runs per city slice (16 players × {1±δ}, plus 120 pairs × 4 corners, δ = 0.05), hourly MDA8 |
| `results_three_kernels_players.csv` | Shapley / time-first Owen / species-first Owen per player, four games |
| `results_optimal_plan.csv` | exact budget-constrained optimum per game and budget pair |

## Pipeline in one page

```
observations                anchoring                 campaigns                 analyses
-------------               ---------                 ---------                 --------
AQS/PAMS (_aqs_*.py)        _fit_anchor.py            _gen16_city.py            _errv_full.py
  -> obs_*_cb05.csv    --->  -> results_fit_*.txt ---> _gen_hddm_city.py   --->  _city_spectrum.py
Open-Meteo (_met_pull.py)   _scan_anchor.py             + sweep_n16.sh            _city_trunc.py
  -> met_*.csv              _diag_fit.py                (KPP workers)             _optimal_plan.py
                                                                                 _fig_*.py
```

Every campaign follows the same pattern: a `gen_*.py` writes one namelist per
coalition (switch-off is exact: window emission multipliers 1 → 0, applied on all
three simulated days), a `sweep_*.sh` runs the namelists on sharded KPP build
directories with crash-resume, and an analysis script reads the merged table. The
currency of all comparisons is **distinct coalition evaluations** (= chemistry
runs).

## Script map

### Model, generation and campaigns (`_kppwork/`)

| script | role |
|---|---|
| `kpp_box.py` | box-model driver: config dict → namelist → KPP run → trajectory; scenario defaults, switch mechanics, MDA8 / d3_mean metric |
| `build_kpp7.sh`, `cb05.kpp`, `cb05_sorg_aq.spc`, `cb05_sorg_aq.def`, `cb05_sorg_aq.eqn`, `cb05*_Main.f90` | CB05 mechanism and KPP build |
| `gen_n16.py`, `gen_co_sweep.py` | namelist generators for the synthetic n = 16 and n = 12 (CO) campaigns |
| `_gen16_city.py`, `_gen_hddm_city.py` | city-anchored n = 16 namelists and Taylor coefficients (8-knob anchored base) |
| `_hourly_tables.py`, `_hourly_core.py`, `_hourly_summary.py` | hourly tables and the `n16_vh_*_m8.csv` oracles written from them |
| `_make_co12_h.py`, `_n12_hourly.py`, `_make_tau2300.py` | n = 12 hourly shards; nocturnal-sink sensitivity build |
| `sweep_n16.sh`, `sweep_co_m8.sh`, `sweep_hddm_city_h.sh` | parallel campaign runners |

### Spectral and truncation diagnostics

| script | role |
|---|---|
| `_errv_full.py` | order mass M(k)/M(1), cumulative share R_q, response-level truncation err_v(q) on the exact n = 16 tables, both metrics |
| `_city_spectrum.py`, `_city_trunc.py` | city-slice spectra and attribution-level trunc(q) |
| `_n16_spectrum.py`, `_m8_spectrum.py` | Möbius spectrum / order masses for the synthetic games |
| `_aw_complete.py`, `_aw_periods.py`, `_aw_spectrum.py`, `_aw_rank_check.py`, `_aw_tablecheck.py`, `_aw_trunc_max.py` | rolling 8-h mean A_w diagnostics behind the SI tables (additivity of the 24 windows, period sums, ranks) |
| `_top_dividends.sh`, `_co_share.py`, `_cum_share.py` | named top dividends, CO-touching mass share, cumulative mass |

### Local vs exact, allocation, optimisation

| script | role |
|---|---|
| `_ddm_predictions.py`, `_ddm_step.py` | per-window T1 / T2 predictions from the 512 fractional runs; step-size check |
| `_hddm_cities.py`, `_hddm_amplitude.py`, `_hddm_strata.py` | T1/T2 scored on every coalition of the city games; error dependence on package composition |
| `_frozen_window_check.py`, `_frozen_rho.py`, `_frozen_base.py` | validity of the frozen maximizing-window approximation |
| `_error_structure.py`, `_single_vs_joint.py` | structure of the local error; single-window vs joint response |
| `_optimal_plan.py` | exact budget-constrained optimum, greedy ranking and order-q truncated selection (writes `results_optimal_plan.*`, `figs/fig_optimal.*`) |
| `attrib_oracle.py` | exact Shapley / Owen values from enumeration |
| `_shapley_taylor.py` | k-truncated Shapley–Taylor values and their error |

### Robustness and audits

| script | role |
|---|---|
| `_dt150_check.py` | outer-step convergence (150 / 300 / 600 s) of the all-on MDA8 |
| `_night_sink_proxy.py`, `_sink_budget.py` | nocturnal sink / deposition budget and the halved nocturnal timescale refit |
| `_realism_check.py` | synthetic scenarios inside the range spanned by the two anchored cities |
| `_kernel_spread.py`, `_tail_impute.py`, `_mlmc.py`, `_mlmc_levels.py` | kernel discretisation spread, tail treatment, multi-level sampling arithmetic |
| `_skill_audit.py`, `_surrogate_cv.py`, `_baselines_city.py`, `_sampler_city.py` | greedy/surrogate baselines and their cross-validated skill |
| `_diurnal_audit.py`, `_mda8_wincheck.py`, `_grid_check.sh`, `_d3_conv.py` | metric-convention and spin-up checks (window containment, 5-min grid, 288 vs 289-point daily mean) |

### Figures and manuscript build

| script | output |
|---|---|
| `make_figs.py` | draft figure set from released tables; in the manuscript it produces `fig3_spectrum.png` |
| `_fig_flip.py` | `fig_signflip_finite.png` (single-window shutdown effects, synthetic pair) |
| `_fig3_with_cities.py` | `fig3b_cities.png` (spectral mass: synthetic vs anchored) |
| `_fig_cities.py` | `fig_cities_24win.png` (single-window effects at the two cities) |
| `_fig_trunc_cities.py` | `fig_trunc_cities.png` (truncation error, four games, both metrics) |
| `_fig_cities_lvg_h.py` | `fig_cities_lvg.png` (T1 / T2 vs exact finite shutdown, hourly MDA8) |
| `_fig_alloc.py` | `fig_alloc.png` (Shapley bars + Owen markers; needs `ATMOS_DIR` set to the repo root on Windows) |
| `_build_aas.py`, `_build_aas_si.py` | assemble the AAS manuscript and SI from the npj-template source |

## Environment

- Chemistry: WSL2 Ubuntu, gfortran + KPP 3.5.0; builds live in `/root/kpp7*`
  (created by `build_kpp7.sh`), campaign workers are sharded copies.
- Analysis: Python 3 with `numpy`, `scipy`, `matplotlib`; no other dependencies.
- Manuscript: MiKTeX (pdflatex + bibtex) for the AAS template.
- Integration: 3 simulated days, 300-s outer step, first 48 h spin-up; MDA8 from
  hourly means with 8-h right-endpoint windows.

## Reproducing

1. Build the box model (`_kppwork/build_kpp7.sh`).
2. Synthetic oracles: `gen_n16.py` + `sweep_n16.sh`; city oracles:
   `_gen16_city.py` + `sweep_n16.sh` (65,536 runs per game). City anchoring
   itself is `_fit_anchor.py` → `results_fit_*.txt`; the 512-run Taylor sets are
   `_gen_hddm_city.py` + `sweep_hddm_city_h.sh` → `hddm_h_v_*.csv`.
3. Diagnostics: `_errv_full.py`, `_city_spectrum.py`, `_city_trunc.py`,
   `_optimal_plan.py`, `_hddm_cities.py` (each reads the released oracle tables,
   zero chemistry runs).
4. Figures: `make_figs.py` and the `_fig_*.py` drivers above; outputs go to
   `figs/`.
5. Manuscript: `_build_aas.py`, `_build_aas_si.py` (LaTeX sources are kept with
   the submission, not in this repository).

Every number in the manuscript traces to a `results_*.txt` written by one of
these scripts; every such script reads either an included oracle table or an
included observation product.

## Cite

Placeholder — a tagged release with DOI will accompany the submission.
