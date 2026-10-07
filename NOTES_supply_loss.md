# Supply–loss stomata (`stomata_model = 9`)

Branch `supply_loss` (off `global_change_ecology` 6fd40b3 (= the old 2e1e2ea/56fa2c9 after the 7 Oct history rewrite; identical tree)). Brief: `~/research/JULES/NOTE_supply_loss_handoff_2026-10-06.md`.
The option is **9** (first 7; moved to 9 on 7 Oct 2026 for the renumbering on branch `cmax`: 5 SOX_opt, 6 CMax, 7 CGain, 8 DESICA, 9 supply-loss). This branch still has the old 5 DESICA and 6 SOX profit.

## What it does

The regulation of Sperry et al. (2016, Eqn 5, with its saturation rule; the idea is from Sperry & Love 2015, Fig. 2), on the profit-max hydraulics. Sperry's fixed G_max is replaced by a carbon-coupled demand. There is no optimisation and no `som_gl_max`.

For each leaf path (`stom_opt_supply_loss` in `src/science/surface/stom_opt_jls_mod.f90`):

1. **Demand.**
   - ci_d = χ₀ ca, with χ₀ = `sl_cica_well_watered_io` (PFT, default 0.8). This is Medlyn et al. (2011) with the water price held fixed, so the D response comes from the hydraulics alone.
   - (A Medlyn √D demand, `l_sl_demand_vpd`, was tried in v1–v2 and **removed in v3**: it counts VPD twice, once in Medlyn's water cost and once in the hydraulics. Its switch-flip results are kept below for the record.)
   - A_d = A(ci_d), gl_d = 1.6 A_d R T / (ca − ci_d), and the unregulated E′ = gl_d D (the same gl↔E conventions as the profit max).
   - A_d ≤ 0 ⇒ closed.
2. **Supply function.** This is the same path as the profit max (`leaf_psi_jls`): single curve or root/stem/leaf segments, the soil-to-root link in series, and gravity in P₀.
   - P₀ = ψ_root_zone (zero flow), k₀ = k(P₀).
   - P′ = ψ_leaf(E′), k′ = −dE/dψ_leaf at E′.
3. **Loss rule (Sperry et al. 2016 Eqn 5).** ΔP = (P₀ − P′) k′/k₀.
   - As E′ rises, ΔP climbs to a maximum and then falls. Past the maximum it is held there.
   - "Past" means either:
     - beyond E_crit (k′ ≤ k_fail); or
     - ΔP falling at E′ (a finite difference at 1.01 E′).
   - The maximum is found by golden section on [0, E′], with states beyond E_crit scored −1.
4. **Regulated E.**
   - E solves ψ_leaf(E) = P₀ − ΔP (Illinois, bracketed on [0, E′] or [0, E_max]), and gl = E/D ≤ gl_d.
   - If gl < gl_d, ci solves A(ci)·1.6 R T = gl (ca − ci) on [ccp, ci_d] (Illinois), and A = gl (ca − ci)/(1.6 R T).
5. **Failure threshold.**
   - k_fail = kcrit (`kcrit_fractional_loss`, 5 % of kmax at 0.95), or `sl_kfail_frac` × kmax when > 0.
   - Sperry et al. (2016) use 0.0005 (0.05 %, "physiological zero").
   - k₀ ≤ k_fail ⇒ closed.

**Bracketed searches.** Every search in the new routine is bracketed:
- E_crit by bisection;
- an 8-point scan, then golden section, for the maximum of ΔP (ΔP need not be unimodal with segments and the soil link);
- Illinois for E and for ci.

So nothing in the routine can run away near E_crit, where dE/dP → 0. The supply function itself comes from `leaf_psi_jls`: the lookup table (`som_psi_solver = 3`, recommended) or Newton (2). The Taylor solver (1) is refused, because it floors k at kcrit. Sperry's model avoids the instability the same way, working on the monotone integral transform of each element and capping at E_crit.

**Soil supply limit** (`l_som_supply_limit`): it reaches the stomata through gl_max (`apply_supply_limit` in sf_stom). Model 7 caps the regulated gl at it, then recomputes E, ψ_leaf and k and re-solves ci. A cap of ~0 closes. This was missing in v1 (review B1, fixed in v2); at FR-Pue it changes the scores by < 0.01.

**Diagnostics** (the profit max's carbon gain / hydraulic cost outputs, reused):
- carbon_gain = A/A_d (the share of the demand's net A kept);
- hydraulic_cost = 1 − k(P′)/k₀ (the conductance the unregulated demand would lose).

## Switches and parameters

| name | where | default | meaning |
|---|---|---|---|
| `stomata_model = 9` | jules_vegetation | — | supply–loss; needs `som_base_parm = 1` |
| `sl_kfail_frac` | jules_vegetation | 0 (= kcrit) | failure threshold as a fraction of kmax |
| `sl_cica_well_watered_io` | pft_params | 0.8 | ci/ca of the demand, in (0, 1) |

- `som_gl_max` is ignored for 7: `gl_leaf_cap = 0` in `sf_stom`, so there is no cap in any canopy scheme and none in the `l_leaf_temp_gc_eq` gc.
- `l_som_rhizo_series` is already on for every `leaf_flux_stom_opt` option, so 7 included.
- The checks that allowed 4/6 now allow 7: can_rad_mod 7, segments, cuticular floor, gravity, supply limit, leaf T, vcmax_psi. `l_som_nsl` stays profit-max only.
- Model 7 needs `som_psi_solver` 2 or 3.

## Canopy schemes

The routine sits behind `stom_opt_mod`, so every scheme that calls it gets 7: big leaf (can_rad_mod 1), multilayer (6) and two-leaf (7).

## Review and versions

The independent review is in `REVIEW_supply_loss.md` (7 Oct 2026). Fixed in **v2** (`jules_sl_v2.exe`; model 4 still bit-for-bit):
- B1, the supply limit;
- B2, the g1 check;
- B3, the solver check;
- R1, the search scales with min(E′, E_crit), not E′ (no spurious closure in dry soil at high demand);
- R2, the pre-scan before the golden section.

Still open:
- R3: the canopy gc spikes with `l_leaf_temp_gc_eq`, which affect the reported gc only, not E;
- R4: full closure when k₀ ≤ 5 % kmax;
- R5: warnings for inputs model 9 ignores;
- rose-meta;
- new diagnostics (gl_demand, ΔP).

## Not done / notes

- **Irreversible cavitation.** Sperry et al. (2016) take ΔP from the original (uncavitated) curve and E from the damaged one. `supply_loss` has no impairment code. If it is merged with `xylem_impairment`, compute ΔP from kmax_ref/conductance_b_pft and E from the impaired curve.
- **Plant capacitance (non-steady state).** Like the profit max, this is steady state within a step. A later option could start the supply curve from DESICA's stem water potential.
- **The LUT path could invert E(P) directly** (E = kmax (S(P) − S(P₀))). At present it uses the bracketed root find as the segments do.
- **photo_al is duplicated** from `stom_opt_bounded_search` (a contained function), to leave the profit-max code untouched for the bit-for-bit test.

## Tests (FR-Pue, `runs/roses/FR_Pue/FR_Pue_dev/supply_loss_tests/`)

The base set is the calibrated profit max `calib_gce_rs_seg_lt_wave1_gpp1_ebc` run_0000: two-leaf, leaf T, root segment, soil link, gravity, supply limit, SoilGrids soils. Its namelists are copied by `make_run.py`; the emulator directory is never written.

1. **Spreadsheet check** (Sperry et al. 2016 sheet: segments 50/25/25, b = 2 MPa, c = 3, kmax 151, ψ_soil 0; driver `sheet_test.f90` linked against the build):
   - ψ_leaf at E = 146.627 is −1.001 MPa with both the Newton and the lookup-table solver (sheet −0.998).
   - E_crit at k = 0.05 % kmax is 270.3 (Newton) / 269.7 (table), against the sheet's 269.5.
   - P_crit is −3.93 MPa against the sheet's ≈ −3.49. The supply curve is flat there (E 267.6 → 269.7 between −2.9 and −3.9 MPa), so P_crit depends on the step and the cut-off; E_crit is the robust number.
2. **Regression (`stomata_model = 4`):**
   - The new exe gives a byte-identical 2003 `run.out.nc` to global_change_ecology 2e1e2ea (`pm2003_ref` vs `pm2003_new`).
   - The full 2000–14 run (`full_a_pm`) has every data variable identical to the emulator's `run.out.nc` (only the header differs).
3. **Smoke, 2003 (7 with ci/ca 0.8; 7 with `l_sl_demand_vpd`):**
   - No NaNs; night gs at min_gl as for the profit max; no night transpiration.
   - **Winter (DJF) daytime gs:** profit max 13–15 mm/s and *rising as light falls* (20 mm/s below 50 W m⁻² SWnet, 8.7 above 300); model 9 4–6 mm/s and flat in light. With √D: 9–16 mm/s.
   - **Annual 2003:** GPP 1106 (4) / 1115 (7) / 1149 (7 √D) g C m⁻²; TVeg 352 / 303 / 306 mm.
   - **Midday ψ_leaf:** 7 is 0.2–0.3 MPa less negative in winter and spring and similar in the July–August drought (−4.5 to −5.0 MPa).
   - **Rare gs spikes:** the reported canopy gs has a few spikes (28 steps a year > 0.1 m/s; with √D up to 51 m/s). They come from the `l_leaf_temp_gc_eq` conversion (gc that delivers the leaf-T transpiration at tstar when qs(tstar) − q1 is small), which `som_gl_max` used to cap. E stays as the stomata chose it: 0.11 mm a year in those steps. The same happens in the profit max with `som_gl_max = 0`. **Open:** cap that gc for model 9, e.g. at the demand conductance?
4. **Full 2000–14 switch flip on the profit-max calibration (favours the profit max).** Daily May–Sep vs ICOS (`score_daily.py`; LE vs LE_CORR):

   | run | GPP RMSE calib / hold | GPP bias | LE RMSE calib / hold | LE bias | TVeg mm/yr |
   |---|---|---|---|---|---|
   | (a) profit max (4) | 1.30 / 1.72 | 0.13 / 0.97 | 24.2 / 26.6 | −7.7 / −7.0 | 386 / 411 |
   | (b) 7, ci/ca 0.8 | 1.40 / 1.82 | 0.30 / 1.09 | 24.1 / 26.5 | −6.5 / −7.3 | 340 / 363 |
   | (c) 7, Medlyn √D demand | 1.61 / 1.98 | 0.49 / 1.27 | 22.1 / 27.9 | −6.9 / −7.9 | 338 / 360 |
   | ref: best profit max with Zhou floor (vfloor) | 1.22 / 1.90 | 0.20 / 1.11 | 22.1 / 24.8 | −8.2 / −7.8 | 371 / 395 |

   **VPD response** (`emulator_sl/plots/VPD_switchflip.png`, daily max VPD, rain-free days):
   - **Wet season:** model 9 gets the observed WUE at low VPD (0.21 vs obs 0.18 at < 0.5 kPa). The profit max is far too low (0.12): it over-transpires when water is cheap, which is the scale-free gain.
   - **Dry season:** model 9's LE is closer to the obs than the profit max's. GPP falls too little at high VPD (4.3 vs obs 3.6 at 3–3.5 kPa), so WUE is too high, before recalibration.
   - The D-free demand gives a VPD response from the hydraulics alone, close to the √D demand's.

   **Runtime:** full 2000–14 incl. spin-up 266 s (4) vs 141 s (7), so 1.9× faster for the whole model; ensemble members 1.5 vs 3.4 min.
   **Failure threshold** (`sl_kfail_frac = 0.0005`, Sperry's 0.05 %, vs kcrit 5 %; `full_d_sl_kf0005`):
   - It makes no difference at FR-Pue (GPP differs by < 1e-12, ψ_leaf by < 2e-6 MPa).
   - As expected: on these curves dP peaks while k is still 40–70 % of k(P0), long before either threshold. The threshold would only matter for full closure in very dry soil.

   **Winter** (Nov–Feb daytime, `emulator_sl/plots/winter_gs_switchflip.png`, same parameters):
   - Profit-max gs *rises* as it gets colder (median 20 mm/s at −2 °C vs 8 at 17 °C) and as light falls (28 mm/s below 50 W m⁻² vs 5 above 500): the scale-free gain.
   - Model 7's gs is 1–5 mm/s and rises with temperature and light, as expected.
   - Winter transpiration is about 40 % lower with model 9 (e.g. 2001: 28 vs 50 mm).

5. **Calibration:** model 9 is calibrated exactly as the best profit max (`gce_rs_seg_lt_vfloor_wave1`), plus `sl_cica_well_watered_io` free over 0.6–0.95 (7 free, 90 members). Harness copy in `supply_loss_tests/emulator_sl/` (ensemble `sl_seg_lt_vfloor_wave1`).
   - The members took 1.5 min, against 3.4 for the profit-max ensemble; none failed.
   - Scores are the harness's (`run_best.py`: daily, all months; calib 2000–07, hold-out 2008–12). ψ is from `evaluate.py`; predawn ψ is compared with JULES ψ_soil.

   | calibration | GPP RMSE calib / hold | LE_CORR RMSE calib / hold | predawn ψ RMSE 03–07 / 09–14 | midday ψ RMSE 03–07 / 09–14 |
   |---|---|---|---|---|
   | profit max, best (`gce_rs_seg_lt_vfloor_wave1` gpp1_ebc) | 0.94 / 0.88 | 18.5 / 20.2 | 0.73 / 0.47 | 0.85 / 0.83 |
   | 7, gpp1_ebc (as the profit max) | 1.01 / 0.99 | 20.6 / 19.3 | **1.84 / 1.57** | **1.79 / 1.81** |
   | 7, gpp1_ebc, root P50 fixed −2.41 (`_p50fix`) | 1.02 / 0.93 | 19.2 / 19.5 | 0.73 / 0.44 | 0.88 / 0.90 |
   | 7, gpp1_ebc + ψ_leaf 2003–07 in the fit (`_psi`) | 0.98 / 0.93 | 18.2 / 18.9 | 0.60 / 0.35 | 0.80 / 0.84 |

   **Best parameters:**

   | | vsl | depth (m) | catch0 | dcatch | root P50 (MPa) | sf_vcmax | ci/ca well watered |
   |---|---|---|---|---|---|---|---|
   | profit max | 22.7 | 1.2 | 1.0 | 0.4 | −2.41 | 1.85 | — |
   | 7 gpp1_ebc | 23.6 | 1.2 | 0.1 | 0.4 | **−8.05** | 1.63 | 0.76 |
   | 7 p50fix | 22.5 | 1.2 | 0.96 | 0.4 | −2.41 (fixed) | 3.94 | 0.79 |
   | 7 psi | 22.8 | 1.2 | 0.58 | 0.4 | −2.16 | 3.66 | 0.81 |

   - **Root P50 identifiability.** With GPP + LE alone, model 9 pushes root P50 to the end of the prior (−8 MPa: a root that hardly cavitates). GPP and LE are still fine, but the soil dries too far: predawn ψ is 1 MPa too negative, as is midday. The fluxes cannot see this (NROY 74 % vs 28 % for the profit max). It is the trade-off of `NOTE_sperry2016_followups` B.
   - **Fixing root P50 at the measured value** (Limousin −2.39; the profit max's own best −2.41), model 9 matches the profit max on GPP and LE and on ψ, out of sample too.
   - **With ψ in the fit,** model 9 finds root P50 −2.2 MPa by itself and is the best of all on LE and ψ. Note that the profit max was not calibrated with ψ, so that comparison is unequal.
   - **ci/ca when well watered** comes out at 0.76–0.81, i.e. Medlyn g1 3.2–4.3 at 1 kPa, against the site Medlyn g1 4.12. It is well constrained and physically sensible.
   - **VPD response, calibrated** (`emulator_sl/plots/VPD_calibrated.png`): model 9 and the profit max are almost the same in the dry season (GPP, LE, WUE vs VPD). At low VPD in the wet season model 9's WUE is above the obs (0.23 vs 0.18) and the profit max's below (0.15).
   - **Recommendation:** calibrate model 9 with root P50 fixed at the measured value (or ψ in the fit), never on fluxes alone.

6. **Unsegmented** (`l_som_plant_segments = .false.`: the whole plant on the stem curve, P50 −6.66 / P88 −7.75 MPa; soil link still in series). Calibrated like the best profit max, v2 exe; root P50 not used.

   | | GPP RMSE calib / hold | LE_CORR RMSE calib / hold | predawn ψ RMSE (bias) 03–07 / 09–14 | midday ψ RMSE (bias) 03–07 / 09–14 | NROY |
   |---|---|---|---|---|---|
   | profit max, segmented (best) | 0.94 / 0.88 | 18.5 / 20.2 | 0.73 (−0.19) / 0.47 | 0.85 (+0.19) / 0.83 | 0.28 |
   | profit max, unsegmented | 1.69 / 1.90 | 33.8 / 39.3 | 3.65 (−2.98) / 3.06 | 3.06 (−2.77) / 2.83 | **0.00** |
   | 7, unsegmented | 1.00 / 0.96 | 21.3 / 20.2 | 2.19 (−1.26) / 1.94 | 2.10 (−1.03) / 2.03 | 0.66 |

   - **Profit max cannot be calibrated without the vulnerable root.**
     - No member is plausible, and the best point sits on the prior bounds (vsl 35, catch0 0.1, dcatch 0.02, sf_vcmax 0.5).
     - LE is far too high in summer and ψ about 3 MPa too negative: with a resistant plant the hydraulic cost is too small, so it keeps transpiring.
   - **Model 7 still fits GPP and LE almost as well as the segmented profit max**, with χ0 0.77.
     - ψ is about 1–1.3 MPa too negative: the soil is dried further than observed, as in the segmented run with root P50 −8.
     - So the fluxes are robust to the missing root information, but ψ needs a vulnerable element (a root segment, or a less resistant whole-plant curve).


## Final status (7 Oct 2026): kept as a branch, not merged

**Decision.** `supply_loss` stays a branch. The hydraulic (Sperry) side is done and verified. The carbon side is defined and consistent, but it is not settled as science, so it is not going into `global_change_ecology` for now.

**Loss rule checked against the Sperry et al. (2016) spreadsheet** (T1, done):
- The test calls `stom_opt_supply_loss` directly on the sheet's plant: segments 50/25/25, b = 2 MPa, c = 3, kmax 151, ψ_soil 0. The demand is set so that E′ = Gmax·D/Patm, with Gmax = 12563.1, and D runs 0.25–6 kPa.
- The saturated state is E = 145.595 and ψ_leaf = −0.9936 MPa, with ΔP_max = 0.99363 MPa. That holds for solvers 2 and 3 and for k_fail at 5 % or 0.05 %.
- A separate Python calculation (exact Weibull integrals, continuous k) gives the same to 5 significant figures.
- The sheet gives 146.627 / −0.998. The 0.7 % gap is the sheet's coarse E steps: its saturated E is one of its own sample points, and estimating k by a forward difference over 1–5 E units moves ΔP_max by 0.3–1.5 %.
- E saturates between D = 1.5 and 2 kPa (sheet: 1.5). Below that the sheet cannot be matched row by row, because its D is the air VPD and it solves a leaf energy balance.
- With P₀ = −1.5 MPa, E saturates by D = 0.5 kPa at E = 24.9.
- The driver (`t1.f90`, plus a copy of the module with the routine made public) is not in the repo. Rebuild it from the description above if needed.

**Fidelity to Sperry, summary:**
- *Matches:* Eqn 5 with k₀ at zero flow on the current soil; the whole-path k′ (soil link and segments); the hold at maximum ΔP; E from ψ = P₀ − ΔP; gravity only through P₀.
- *Differs:*
  - no irreversible cavitation;
  - the rhizosphere is linear and fixed within the step, where Sperry uses van Genuchten per layer;
  - the failure threshold is applied to the whole path, not to each element;
  - sunlit and shaded leaves have separate parallel paths, each with its own canopy ψ (as in model 4).

**Carbon side, why it is not settled:**
- ci/ca is prescribed at χ₀ when water is plentiful (Medlyn with a fixed water price); nothing is optimised.
- χ₀ is a ceiling. The loss rule always cuts a little, so realised ci/ca < χ₀: 0.788 at D = 1 kPa on the sheet plant with χ₀ = 0.80. Mapping χ₀ = g1/(g1+1) from Medlyn tables therefore gives somewhat low gs.
- The default 0.8 is wrong for C4 (~0.4).
- χ₀ has no temperature dependence.
- Closure is on net A_d ≤ 0.
- Untested: a CO₂ ramp (280→800 ppm, wet and dry, 4 vs 9) and the Harwood winter case (RED, can_rad_mod 6, Jan–Mar 2002), the motivating failure.
- A possible way forward is χ₀ from least-cost theory (P-model ξ(T)), which adds temperature dependence.

**Still open in the code:**
- R3, the gc spikes with `l_leaf_temp_gc_eq`;
- R4, full closure when k₀ ≤ 5 % kmax;
- R5, warnings for ignored inputs;
- per-PFT/C4 defaults for `sl_cica_well_watered_io`;
- diagnostics gl_d, E′ and a saturated flag;
- rose-meta.

**Follow-up for the shared hydraulics (do it off `global_change_ecology`, not here).** The nonlinear rhizosphere, Kirchhoff transform with van Genuchten:
- E_j = G_j [Φ_j(ψ_soil,j) − Φ_j(ψ_root)];
- Φ is static for each soil layer, so it can be tabulated once.
- (A) a single lumped node replacing the linear `ksr_path`: ~1–2 days.
- (B) per-layer to a common root-crown ψ, as Sperry does: ~1 week, with a runtime cost.
- It sits in `leaf_psi_jls`, so profit max gets it too.
- First run an offline check of its size at FR-Pue: k′/k₀ at midday E in summer drought, linear vs nonlinear.
