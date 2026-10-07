# Review: supply–loss stomata (`stomata_model = 7`), branch `supply_loss`

Reviewed 2026-10-07 against the uncommitted `git diff` on base 2e1e2ea. I did not edit, build or run anything.
Line numbers refer to the working tree.
Main routine: `src/science/surface/stom_opt_jls_mod.f90` (SUBROUTINE `stom_opt_supply_loss`, L1858–2237; dispatch L617–644).

**Verdict.**
- The core algorithm is a faithful and numerically careful implementation of Sperry & Love (2015) Fig. 2 and Sperry et al. (2016) Eqn 5.
- The carbon-coupled demand is correctly derived.
- There is **one confirmed functional bug**: the soil-supply limit is silently dropped for model 7.
- There are a few robustness gaps, at tolerances, with the Taylor solver and with an unchecked g1.
- None of this invalidates the FR-Pue conclusions qualitatively. However, the 7-vs-4 comparisons ran with `l_som_supply_limit=.true.`, which acted on 4 but not on 7.

---

## 1. Fidelity to the papers

### What matches (verified against the PDFs)
- **Eqn 5:** ΔP = ΔP′ (dE′/dP_canopy)/(dE/dP_max). In Sperry 2016 (p. 580), dE/dP_max is "the derivative of the supply function … at P_canopy = P_soil".
  - The code has dP = (P0 − P′)·k′/k0 with k0 = supply(0) (L1976, `dp_of` L2148), which is the derivative at zero flow on the *current* soil.
  - **Correct.** It also matches SL15 ("maximum dE/dP_canopy is at the predawn start of the curve").
- **k = −dE/dψ_leaf of the whole path is the right quantity.** Sperry's supply function runs from bulk soil to canopy, so its derivative includes the rhizosphere and every element.
  - `leaf_psi_jls` returns exactly this. For the single curve with the soil link, k = k_p/(1 + k(ψ_in)/K_s) (`xylem_hydraulics_jls_mod.f90` L209–218). I re-derived it: dψ_l/dE = −(1/k_p)(1 + k_in/K_s) because dψ_l/dψ_in = k_in/k_p. ✔
  - For segments, the recursion dψ_out/dE = (k_in dψ_in/dE − 1)/k_out starts from −1/K_s (`…cumulative_weibull…` L930–985, L1061–1076). ✔
  - Gravity enters only through P0 (ψ_src), as Sperry has it.
- **Saturation rule.** Sperry 2016: "ΔP saturates at its maximum as E′ increases". SL15: past that point stomata "keep E and P_canopy constant".
  - The code holds dp at max_{E∈[0,E′]} dP(E) once E′ is past the peak (L1994–2045).
  - The regulated E then comes from ψ(E) = P0 − dP_max. This is < E_peak because dP_max ≤ P0 − ψ(E_peak), which is Sperry's step (4). ✔
  - It is continuous across the transition.
- **Unimodality.** Sperry states it ("ΔP rises to a maximum before decreasing back to zero").
  - For a single Weibull it is provable. Let x = P0 − P. Then d ln(x k)/dx = 1/x − (c/b^c)(x+|P0|)^(c−1). Its root is unique for any c > 0, because x(x+a)^(c−1) has derivative (x+a)^(c−2)(cx+a) > 0.
  - The linear soil link and the series segments are **not** proven unimodal. See R2.
- **Demand substitution.** Medlyn (g0 = 0): g_s = 1.6(1 + g1/√D)A/ca, so ci/ca = g1/(g1 + √D).
  - Holding D fixed gives ci = χ0·ca, gl_d = 1.6 A R T/(ca − ci_d) (L1970), and χ0 = g1/(g1+1) at 1 kPa (pftparm_mod comment).
  - The √D variant (L1960–1961) uses D in kPa, the same as `leaf_limits_mod.F90` L219–220. ✔
  - E′ = gl_d·VPD/(RT) is the same conversion as `hydraulic_state` (L1613). Sperry's E′ = G_max·D, with D as a mole fraction, is the same quantity. ✔

### Divergences and what is missing
1. **Failure threshold.**
   - Sperry fails each *element* at k < 0.05 % of continuum kmax.
   - Ours compares *whole-path* k (including the soil link) with `kcrit` (5 % of plant kmax), or with `sl_kfail_frac·kmax`.
   - It is used in only three places: k0 ≤ k_fail ⇒ closed (L1977); E′ beyond E_crit ⇒ saturated (L1994); and the −1 scoring (L2161).
   - Because dP peaks at k ≈ 40–70 % of k0, the threshold rarely matters. That agrees with `full_d_sl_kf0005`.
   - Two caveats:
     - With the soil link, a wet-soil k0 can fall below 5 % kmax purely through soil resistance. The point then **closes completely**, where Sperry would still regulate. 5 % is a harsher "physiological zero" for k0 than Sperry's.
     - With the Taylor solver it cannot work below kcrit (B3).
2. **Irreversible cavitation is not represented.** Sperry's default takes ΔP from the original curve and E from the damaged one. There is no impairment on this branch; this is documented in NOTES.
3. **Rhizosphere.**
   - Sperry has a van Genuchten k(P) per soil layer, solved jointly for a common root-crown P, and fits needed about 67 % rhizosphere resistance.
   - Ours has one bulk root node and a *linear* K_s that is fixed within the step. The soil term therefore lowers k0 but never makes k′/k0 fall with E.
   - Model 7 cannot express rhizosphere-driven regulation. That probably explains why calibration moves root P50 (follow-ups note B).
4. **No hydraulic redistribution** (`MAX(…,0)` uptake). This is the same as in profit max.
5. **The D response with the D-free demand comes only from the hydraulics.**
   - Sperry also has G = G_max until the loss bites ("D threshold causing G < G_max depended on…"). So with high kmax or resistant curves, gs barely responds to D.
   - At FR-Pue the vulnerable root segment provides the response. For global PFT parameters it may be too weak (S1).
6. **Lost temperature dependence of χ.** Medlyn's ξ ∝ √(Γ* λ) varies with T. Fixing χ removes that weak temperature dependence (S8).

---

## 2. Numerics and correctness

### Confirmed bugs

**B1. `l_som_supply_limit` is silently ignored for model 7.** (confirmed)
- `sf_stom` applies the soil-supply cap only through `gl_max_eff` (`apply_supply_limit`, `sf_stom_jls_mod.F90` L4582–4629). This happens in every scheme: multilayer L1782, big leaf L2385, two-leaf L2743/2765/2823/2845. The cap is then passed to `stom_opt_mod` as `gl_max`.
- The `CASE (supply_loss_model)` call (L623–631) **does not pass `gl_max`**, and `stom_opt_supply_loss` has no cap.
- Consequences:
  - The check message, which says "requires stomata_model=4, 6 or 7", and NOTES ("supply limit … now allow 7") are wrong in effect.
  - The FR-Pue namelists (`full_b_sl`, the ensembles) have `l_som_supply_limit=.true.`, so profit max ran with the cap and 7 without it. Only the `sf_evap` backstop remains.
- The size of the effect is unquantified. It is likely small, because it only binds when the extractable water per step is small, but it matters in dry spells and in global runs.
- **Fix:** pass `gl_max` in. After step 3, if `gl_max(l) > 0` and `gl_reg > gl_max(l)`, set `gl_reg = gl_max(l)`. Then set e_reg = conv_e·gl_reg and re-solve ψ, k with `supply(e_reg)` before the ci re-solve. Also make `gl_d ≤ gl_max` in the `gl_reg ≥ gl_d` branch. A TINY cap then closes, as intended.

**B2. `l_sl_demand_vpd` uses `g1_stomata` with no check.** (confirmed, minor)
- The default is `rmdi` (`pftparm_mod.F90` L621) and there is no validity check anywhere.
- With g1 unset, ci_d ≈ ca·(1+ε). Then ca − ci_d hits the 1e-2 Pa floor and gl_d becomes about 1e3× too large. The run proceeds with garbage, silently.
- **Fix:** in `pftparm_io_mod` add `IF (l_sl_demand_vpd .AND. stomata_model==7 .AND. ANY(g1_stomata<=0)) ereport 101`.

**B3. Model 7 is allowed with `som_psi_solver = 1` (Taylor), which is the namelist default (`jules_vegetation_mod.F90` L357).** (confirmed)
- Taylor floors k at kcrit (`…cumulative_weibull…` L477). So:
  - with `sl_kfail_frac` < kcrit/kmax, `k_d <= k_fail` (L1994) and the −1 scoring never fire;
  - the supply function is a midpoint approximation, so −dE/dψ is inconsistent with ψ(E).
- **Fix:** require `som_psi_solver = 3` (or 2/3) for model 7 in `check_jules_vegetation`. At minimum, give a warning.

### Concerns (not bugs)

**R1. The tolerances are relative to E′, not to the scale of the answer.**
- The golden search exits at `gb − ga ≤ 1e-3·E′` (L2016). Illinois exits at `xb − xa ≤ 1e-4·E′` (L2064).
- In dry soil with high demand, E_crit can be ≪ E′.
  - If E_crit < ~1e-3 E′, both golden probes stay in the −1 region until exit, which gives `dp_max ≤ 0` and a spurious **closure** (L2034). Whether it closes then depends on the demand, not on the hydraulics.
  - For E_peak ≈ 0.01 E′, E_max is resolved only to about 10 %. That is benign for dp_max, which is flat at the top, but the regulated E then inherits the coarse Illinois tolerance.
- **Fix:**
  1. When k_d ≤ k_fail, first bisect for E_crit on [0, E′] (k(E) = k_fail).
  2. Run the golden search on [0, min(E′, E_crit)], with tolerance relative to that bracket.
  3. In Illinois, use `e_rtol·xb_initial`, not `e_rtol·e_d`.

**R2. Unimodality is assumed, not guaranteed.** This applies to the segment chain plus the linear soil link, and to `psi_solver_newton`.
- `psi_solver_newton` runs four Newton iterations with a one-sided error near the ceiling (L370 comment: ψ not negative enough, k overstated). That can put a bump in dP(E).
- With two local maxima:
  - the 1 % slope test (L1998–1999) can flag "past the peak" at the first maximum;
  - golden can converge to either maximum;
  - in both cases dP is underestimated, so regulation is too strong.
- The LUT path is safe in practice. The grid is about 1 kPa, and a 1 % step in E′ near the peak moves ψ by tens of kPa, so the finite difference is not grid noise.
- Misclassification just *before* the peak is harmless: golden then returns ≈ E′ and dp ≈ dp_d (continuous).
- **Fix:** a cheap 6–8-point pre-scan of dP on [0, min(E′, E_crit)] to pick the bracket around the largest sample before the golden search. Also recommend (or require) `som_psi_solver = 3` for model 7.
- NOTES' "No Newton steps" is true of the new routine only. With solver 2 the underlying `leaf_psi_jls` is Newton.

**R3. Spikes in the reported canopy gc with `l_leaf_temp_gc_eq`** (`sf_stom` L2955–2969). The cap is `gl_max_sun_2l + gl_max_shd_2l`, which is 0 for model 7. This is already noted in NOTES as open.
- Suggested cap that needs no fixed g_max: `gc ≤ f × (gl_sun + gl_shd)` with f ≈ 3. Alternatively, return gl_d from `stom_opt_mod` and cap at Σ gl_d.
- Either way, E is unchanged; only the diagnostic or the energy-balance gc is bounded.

**R4. The 1e-2 Pa floor on ca − ci.** It is safe in the χ0 path because ca − ci_d = (1−χ0)ca ≥ 0.05 ca. In the √D path it is reached only through B2. At the dq_min floor (VPD ≈ 16 Pa), ci/ca ≈ 0.97 is still well above the floor.

### Verified OK
- **Units.** vpd = dq·p*/ε (Pa). conv_e = vpd/(R T) gives E [mol m⁻² s⁻¹] = conv_e·gl [m s⁻¹], identical to `hydraulic_state` L1612–1613. kmax/kcrit/k0 are all in the path's kmax units (mol m⁻² s⁻¹ Pa⁻¹); the ratio is 1.6.
- **ci re-solve bracket** (L2093–2117). h(ci) = gl(ca−ci) − 1.6RT·A(ci) is strictly decreasing.
  - h(max(ccp,0)) = gl(ca−ccp) + 1.6RT·rd > 0, because A(ccp) = −rd. This holds for Farquhar and Collatz C3; for C4, ccp = 0 and A(0) = −rd.
  - h(ci_d) = (gl−gl_d)(ca−ci_d) < 0 when gl_reg < gl_d.
  - `ci_d ≤ ccp` is closed earlier (L1966), so the bracket is always valid.
  - The returned A = gl(ca−ci)/(1.6RT) > 0 is consistent with gl. The A of a nearly closed open leaf therefore tends to 0, not −rd. That is fine physically, but it is a small discontinuity against `set_closed`.
- **E bracket.** f(0) = dp > 0 because ψ(0) = P0 for every solver. f(xb) ≤ 0 because dp ≤ P0 − ψ(xb) whenever k(xb) ≤ k0. If k(xb) > k0, which would be numerical, the loop is skipped and E = xb. This degrades gracefully.
- **No division by zero.** k0 > k_fail ≥ 0 before every `dp_of`. conv_e > 0 (dq ≥ dq_min). al_d > 0 before `/al_d`. The Illinois denominators have opposite-sign brackets.
- **No uninitialised outputs.** Every path ends in `set_state` or `set_closed`, and `hydraulic_cost_g` / `carbon_gain_g` are set on every path.
- **The closed state** (ci = ca, A = −rd, gl = 0, k = k(P0), ψ = P0, E = 0, CG = HC = 0) is identical to `stom_opt_bounded_search`'s `set_closed`. Points not in `open_index` keep `min_gl_pft` from `stom_opt_mod`, as for model 4.
- **`photo_al` copy:** byte-identical to the original (diffed L1639–1683 against L2191–2235).
- **Diagnostics.** They are reused with new meanings: CG = A/A_d, HC = 1 − k(E′)/k0. Both are bounded in [0,1]. Document that the HC definition differs from model 4's (k relative to kcrit) (D2).
- **OpenMP.** There are no SAVE or module variables in the routine; everything is local or host-associated within the call. `ksr_path` is read-only. The LUT builds are inside OMP CRITICAL. The new `gl_leaf_cap` assignments in `sf_stom` lie outside every `DEFAULT(NONE)` region (checked: L1477, L2383, L2575). ✔
- **Behaviour difference.** Model 7 closes on *net* A_d ≤ 0 at χ0·ca (L1966). Model 4 with `l_som_gain_gross=.true.` (the FR-Pue setting) closes only on gross A ≤ 0. So in very dim light, 7 closes earlier. This is intended, but `l_som_gain_gross` is silently ignored.
- **Efficiency.** `set_closed` repeats `supply(0)`, and k0 could come from `k_path_zero_flow`. E_max depends only on (P0, curves, K_s, kmax) and could be cached per point between the leaf-T iterations. Both are minor.

---

## 3. Plumbing

- **Namelists.**
  - `jules_vegetation`: n_real 18→19 and n_log 43→44 match the type. The type has 24 REAL lines, of which 5 are `n_photo_coef` arrays, giving 19 scalars; it has 45 LOGICAL lines, of which 1 is an npft_max array, giving 44 scalars.
  - The new members sit inside the contiguous REAL and LOGICAL blocks of the SEQUENCE type, and read, broadcast and print are all present. ✔
  - `jules_pftparm`: 136 REAL array members for n_real = 136·npft_max. ✔ It also has allocation, default 0.8, the rmdi-WHERE copy, the (0,1) check (only when 7) and print. ✔
- **Checks that include 7:**
  - can_rad_mod 7 (L917)
  - segments (L1384)
  - cuticular floor (L1394)
  - gravity
  - supply limit (L1410; ineffective, B1)
  - leaf T
  - vcmax_psi (L1471)
  - `som_base_parm = 1` (L1029)
  - `l_som_rhizo_series` set for all `leaf_flux_stom_opt` (L1428) ✔
- **Checks that correctly exclude 7:** `l_som_nsl` (L1480).
- **Missing checks:**
  - Taylor solver (B3).
  - g1 with `l_sl_demand_vpd` (B2).
  - A warning that `som_gl_max`, `som_ci_search`, `som_n_sample`, `l_som_skip_search_wellwatered` and `l_som_gain_gross` are ignored for 7.
  - `sl_kfail_frac > 0` with any model other than 7 (it is silently unused).
- **Old namelist form** (`leaf_flux_mod=2`, `som_profit_model=3`). It maps to model 4 and then fails the new som_profit_model check. That is an error, which is acceptable.
- **Canopy schemes.** Big leaf (1), multilayer (5/6) and two-leaf (7, with and without leaf T) all call `stom_opt_mod` with the same interface, so all of them get 7.
  - Multilayer: there is no per-leaf gl_max for 7 (`gl_max_lf = 0`), and segments and leaf T are not available there.
  - The cuticular floor is applied after the stomata for 7, as for 4.
- **Remaining `stomata_model` tests in `src`.** The full grep covers `smc_ext` L244 (gated on `leaf_flux_mod`), `sf_evap` / `sf_diags` (DESICA only) and the `dq_min` gate (updated). There is no other place that should include 7. ✔
- **rose-meta not updated.** `vn7.9/rose-meta.conf` L7995–8016 still lists `values=1,2,3`. This is a pre-existing gap for 4–6 as well, but 7 and the three new inputs should be added when the metadata is caught up.
- **Outputs.** There is no new diagnostic. It would help to have gl_d (demand), E′ and the saturated flag, so the regulation can be seen in output.

---

## 4. Viability against profit max

**Evidence (FR-Pue, NOTES and logs).**
- *Switch flip (favours 4):*
  - 7 is close on LE (24.1 vs 24.2 calib) and worse on GPP (1.40 vs 1.30).
- *Calibrated:*
  - With root P50 fixed or ψ in the fit, 7 matches or beats 4. Examples: ψ fit LE 18.2/18.9 vs 18.5/20.2; predawn ψ RMSE 0.60/0.35 vs 0.73/0.47.
  - It does so at 1.5 vs 3.4 min per member.
  - The VPD response in the dry season is almost identical (`plots/VPD_calibrated.png`, checked). At low VPD in the wet season, 7 overshoots WUE (≈0.23–0.24 vs obs 0.18) and 4 undershoots (0.15).
- *Fluxes alone (7):* root P50 is non-identifiable and goes to −8 MPa, and ψ is off by about 1 MPa.
- These are with B1 present.
- The unsegmented ensemble (`ensemble_sl_lt_vfloor_noseg_wave1`) was still running at review time (run_0017 started 09:24). It cannot be judged yet.

**Reasons to favour 7.**
- No `som_gl_max`, and no scale-free gain. gl ≤ gl_d ∝ A_d/ca, so dim, cold leaves cannot open wide. This removes the Harwood winter failure *by construction*, but it is untested there (see S5).
- About 1.9× faster per run, and deterministic (bracketed, no sampled grid).
- One interpretable PFT parameter. χ0 can be mapped from Medlyn/Lin g1 tables as χ0 = g1/(g1+1). That makes it attractive for **global runs**.
- A Medlyn-like CO2 response in the unregulated range (below).

**Reasons to favour 4, or to be cautious with 7.**
- 7 has no explicit carbon optimality: ci/ca is prescribed when unstressed. Model 4's ci/ca emerges from the A and k curves.
- In 7 the D response with the D-free demand relies entirely on the hydraulic vulnerability. That is fine with a vulnerable root segment, but probably too weak for resistant or high-kmax PFTs. Use `l_sl_demand_vpd` there.
- The isohydric/anisohydric spectrum is set entirely by curve shape and kmax in both models. 7 has no extra strategy parameter. Its saturation rule makes it strongly isohydric once saturated.
- Identifiability: 7 needs ψ data or a fixed root P50. The flux-only fit is degenerate.

**Expected CO2 response (to be tested).**
- **7 (χ0):**
  - Unstressed: ci = χ0 ca, so gs ∝ A/ca. A 1.5× rise in ca gives about −15 to −25 % gs (Medlyn-like). Lower E′ also means less regulation.
  - Saturated: regulated E = E(P0 − dP_max) is independent of ca, so **gs is CO2-insensitive in drought** and the gain is all in A, through higher ci.
- **4:**
  - The response emerges from the shape of A(ci) against the cost.
  - With this branch's normalised gain, (A+g_off)/max_al, a CO2 rise mainly rescales A. It shifts the optimum only through the curvature of A(ci), so gs probably decreases weakly and depends on whether photosynthesis is light- or Rubisco-limited.
  - At low cost it still goes to `som_gl_max`, which is not a CO2 response at all.
- A ca ramp (280→800 ppm) under wet and dry soil at FR-Pue for both models would settle this.

**Bottom line.** 7 is a viable alternative switch, especially for global or long runs. Make it production-ready after B1–B3 are fixed, with a Harwood test (S5) and a CO2 test (S4). Keep profit max as the "optimality" reference.

---

## 5. Prioritised to-do

| # | Type | Item | Suggested change |
|---|---|---|---|
| B1 | **Bug** | Soil supply limit ignored for 7 (L623–631: no `gl_max`) | Pass `gl_max`; cap `gl_reg`/`gl_d`; recompute E, ψ, k, ci. Re-run `full_b_sl` to quantify |
| B2 | Bug (input) | `g1_stomata` unchecked with `l_sl_demand_vpd` | ereport if g1 ≤ 0 / rmdi when 7 + `l_sl_demand_vpd` |
| B3 | Bug (config) | Taylor solver (default) allowed; k floored at kcrit | Require `som_psi_solver = 3` (or 2/3) for 7 |
| R1 | Robustness | Golden/Illinois tolerances ∝ E′ ⇒ spurious closure when E_crit ≪ E′ | Bisect E_crit first; golden on [0, min(E′,E_crit)]; Illinois tol ∝ upper bracket |
| R2 | Robustness | Unimodality assumed (segments + soil link, Newton noise) | 6–8-point pre-scan to bracket the global max; prefer LUT |
| R3 | Robustness | Canopy gc spikes with `l_leaf_temp_gc_eq` (no cap for 7) | Cap at f·(gl_sun+gl_shd) or at Σ gl_d (output gl_d from stom_opt) |
| R4 | Robustness | k0 ≤ 5 % kmax closes fully (soil link in wet soil) | Consider `sl_kfail_frac` default 0.0005 (Sperry), or apply k_fail to plant k only |
| R5 | Hygiene | Ignored inputs (`som_gl_max`, `l_som_gain_gross`, `som_ci_search`, …), `sl_kfail_frac` unused unless 7 | Warnings in `check_jules_vegetation` |
| S1 | Science | D response only via hydraulics with χ0 demand | Test resistant/high-kmax PFTs; recommend `l_sl_demand_vpd` for global |
| S2 | Science | No irreversible cavitation | On merge with `xylem_impairment`: dP from intact curve, E from impaired curve |
| S3 | Science | Linear, fixed K_s; one root node; no rhizosphere VG k(P) | Per-step nonlinear rhizosphere k(P) (Sperry Eqn 2) would let 7 regulate on the rhizosphere |
| S4 | Science/test | CO2 response untested | ca ramp, wet/dry, 4 vs 7 |
| S5 | Science/test | The motivating case (Harwood winter, can_rad_mod 6, no segments) is not run with 7 | Run RED Jan–Mar 2002 with 7 vs standard Jacobs vs 4 |
| S6 | Science | Root P50 non-identifiable from fluxes | Calibrate with ψ (as `_psi`) or fixed root P50; report conditional on kmax/soil |
| S7 | Science | Closes on net A_d ≤ 0 (gross for 4) | Document; optional gross-A variant for comparability |
| S8 | Science | χ0 has no temperature dependence (Medlyn ξ ∝ √Γ*) | Optional χ0(T) via ξ(T)/(ξ(T)+1) |
| D1 | Docs | rose-meta lacks 7 and the new inputs | Add `values=…,7`, `l_sl_demand_vpd`, `sl_kfail_frac`, `sl_cica_well_watered_io` |
| D2 | Docs | carbon_gain / hydraulic_cost change meaning under 7 | State in output docs/NOTES; add gl_d, E′, saturated-flag diagnostics |
| D3 | Docs | NOTES says the supply limit works for 7 and "No Newton steps" | Correct both (B1, R2) |
| T1 | Test | Spreadsheet check covers the supply function only, not the loss rule | Debug mode with fixed G_max demand; reproduce sheet G vs D and the dP saturation point |
| T2 | Test | Bit-for-bit for model 4 done ✔ | Repeat after B1 (B1 touches only the 7 branch) |
