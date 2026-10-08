# Branch `stomata_fixes` (off `global_change_ecology` 40d1654, 8 Oct 2026)

Uncommitted on the branch; the user commits.

## Changes

1. **CGain (7) reference:** `k_unstressed` is now the plant's conductance at ψ = 0 without the soil-to-root link. This uses the new optional argument `k_path_zero_flow(..., l_plant_only=.TRUE.)`. Before, the link was included, so the reference shrank as the soil dried and the same embolism cost less in dry soil (Lu et al. 2020 use the plant's k_max).
2. **SOX_opt (5) comments:** model 5 is Eller et al. (2018) solved as Sabot et al. (2022) do. Eller's own definition is the halfway point ψ̄ = (ψ_root + ψ_leaf)/2 (Eqs 2.7–2.8). Sabot's Eq 16 (k at ψ_leaf) is not followed. The kcrit subtraction is Sabot's. Eller's cost is xylem-only; ours includes the soil link. No code change.
3. **CAP, new `stomata_model = 10`** (`cap_profit_model = 6`), after Dewar et al. (2018):
   - max A_n, where gross A is scaled by f = 1 − (ψ_leaf − ψ_nsl_onset)/(ψ_nsl0 − ψ_nsl_onset);
   - no stomatal cost; k > kcrit still applies;
   - reuses the `l_som_nsl` machinery, which is switched on with it;
   - parameters: `psi_nsl0_io` (Sabot's Ψφ,lim) and `psi_nsl_onset_io` (0 for Dewar's form);
   - needs `som_ci_search = 2`.
4. **`l_som_coupled_e` (jules_vegetation, default .false.):** the optimisation's transpiration for a trial g is
   - E = D·g′/(1 + r_ca(g′ + g′_other)), with g′ = g/(1 + g·r_bl);
   - previously E = g·D_c, with D_c held fixed from the last pass.
   - **Two-leaf with l_leaf_temp:**
     - D = qs(T_leaf) − q1;
     - r_bl = 1/gb from the leaf energy balance (`leaf_temp_update` now returns gbw per class and ra_m);
     - r_ca = the leaf_aero_model's resistance (0 for model 0);
     - g′_other = the other class.
   - **Without l_leaf_temp:** D = dq and r_ca = ra.
   - **Applies to:** big leaf and two-leaf, models 4–7 and 10. The supply cap (`apply_supply_limit`) and `edge_by_gl_cap` invert the same relation.
   - **State:** the first pass of a timestep takes last timestep's class gl and gb (module array `cpl_state`).
   - **Not supported:** multilayer and supply-loss (9) are rejected by the check.

## Tests (FR-Pue 2000–14, from the orcap3 raw_twb_ms best run's 2000 dump)

Runs are in `runs/roses/FR_Pue/FR_Pue_dev/stomata_fixes_tests/`; scores are in `score_2000_14.txt` (`score.py`).

- **Regression:** with the new switches off, `run.out.nc` is byte-identical to exe 029e020_or_rt0.
- **Coupled E is small at this site:**
  - leaf_aero_model 0 (r_bl only): annual GPP +0.8 %, TVeg −0.3 %; daytime TVeg changes by −2.2 % to +0.4 % (5–95 %).
  - leaf_aero_model 2 + coupling vs model 2: similar size.
  - No systematic run-time cost (229–271 s, noise under load).
- **CGain fix:** negligible here (TVeg +0.1 %). The soil link rarely dominates the path at FR-Pue.
- **CAP:**
  - ψφ,lim = −4 MPa: GPP 976 vs 1369 g C m⁻² yr⁻¹.
  - ψφ,lim = −6 MPa: GPP 1062.
  - Better Jul–Aug 2003 GPP (0.89–0.96 of obs vs 0.77).
  - Less negative midday ψ.
  - 2.3× slower.
  - Not calibrated; it uses the profit-max set.

## CMax: Wolf vs Anderegg, and options (8 Oct 2026)

Full report: `~/research/JULES/NOTE_cmax_theta_deep_dive_2026-10-08.md`.

**What the code implements:**
- Wolf et al. (2016) give the criterion, A_n − Θ(ψ_L) with Θ concave-up. They give no formula for Θ.
- Anderegg et al. (2018, Eqs 7–9) fit the slope, ∂Θ/∂ψ_L = aψ_L + b, so Θ = a/2 ψ² + bψ + c.
- Our Θ = a/2 p² + b p, with p = |ψ|, is the same model.

**Sign and scale conventions:** Anderegg's ψ is signed, so their b ≤ 0 is our b ≥ 0. Wang et al. (2020) write aP², so their a is half of ours. Convert fitted values before moving them between papers.

**b is a water price near ψ = 0:** Wolf's Eq 14 gives λ = −Θ′/K, which is about b/K at low tension. The FR-Pue a×b grid (`cmax_tests/g2003_a*_b*`) confirms that b matters there, against Sabot's "low influence".

**Options (from the deep dive; not implemented):**
- b from Medlyn's g1: b ≈ 1e3·K0·3Γ*·P_atm/(1.6·g1²), about 0.57 at FR-Pue. This is the worker's own derivation (wet-soil λ = b/K0 matched to g1's λ) and needs checking.
- a from a closure point: a = [(Ca − Γ)/(1.6·D_ref)·K(P_close) − b]/P_close.
- Optionally normalise both by K0.
- A viscosity (temperature) correction on K, for all hydraulic models.
- CMax with the hysteretic K of xylem_impairment.
- Test plan T1–T5 is in the report.
- **Caution:** February LE is +55 W m⁻² in every model, and CMax is low in May–Jun 2003. Sort out soil evaporation and interception before tuning stomatal costs to winter or spring LE.

## CGain: keep (8 Oct 2026)

In the calibrated comparison (`emulator/stomatal_models_comparison.csv`, gpp1_raw_tw), CGain has the best LE RMSE of models 4–7 (cal 12.7, val 13.7, against profit max 13.4/14.7). Its GPP is slightly worse (1.08/0.99 against 1.00/0.90).

Its cost is PLC-based, so it does not fix dim-light winter over-opening. That is not a reason to retire it.
