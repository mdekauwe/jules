# CMax stomata (`stomata_model = 6`), branch `cmax`

**Stomatal model numbers on this branch** (renumbered 7 Oct 2026; 1–4 unchanged):

| value | constant | scheme |
|---|---|---|
| 1 | `stomata_jacobs` | original JULES (Jacobs closure; Best et al. 2011) |
| 2 | `stomata_medlyn` | Medlyn et al. (2011, Global Change Biology) |
| 3 | `stomata_sox_analytical` | SOX, semi-analytical (Eller et al. 2020); was `stomata_sox` |
| 4 | `stomata_profit_max` | ProfitMax (Sperry et al. 2017; De Kauwe et al. 2022, New Phytologist; Baguley et al. 2026, Biogeosciences) |
| 5 | `stomata_sox_opt` | SOX_opt (Eller et al. 2018); was 6, `stomata_sox_profit` |
| 6 | `stomata_cmax` | CMax (Wolf et al. 2016); was 8 |
| 7 | `stomata_cgain` | CGain (Lu et al. 2020); was 9 |
| 8 | `stomata_desica` | DESICA (De Kauwe et al. 2020); was 5 |
| 9 | `stomata_supply_loss` (branch `supply_loss`) | supply–loss (Sperry et al. 2016); was 7 |

Models 4–7 and supply–loss are checked bit-for-bit against their runs under the old numbers (`cmax_tests/v5_*`, `supply_loss_tests/sl2003_v4_m9`). DESICA appears in the code only as its constant.

Off `global_change_ecology` d80824c (7 Oct 2026). Carbon maximisation of Wolf et al. (2016, PNAS 113: E7222), in the form of Anderegg et al. (2018) and Sabot et al. (2022, JAMES, Eq. 12):

  maximise  A_n − Θ(ψ_leaf),   Θ = a/2·ψ² + b·|ψ|   (ψ in MPa; a, b per unit leaf area)

- **Absolute carbon units:** no normalisation by the instantaneous A_max (unlike profit max). There is no water price: in Wolf's argument, under competition water saved is taken by a neighbour, so the competitive optimum maximises carbon each instant.
- **Same hydraulics as profit max:** `leaf_psi_jls` (single curve or segments, soil-to-root link, gravity), feasibility k > kcrit, the gl_max cap (`som_gl_max` and the soil supply limit).
- **Search:** `stom_opt_bounded_search` with profit = A_n − Θ(ψ_leaf); A_n is net (`l_som_gain_gross` is not used). A − Θ is concave in Ci, so the golden section applies.
- **Θ basis:** Θ is scaled by kmax/kmax_pft, the leaf area of the call (canopy, sun/shade class or leaf), to match A.

## Switches and parameters

| name | where | default | meaning |
|---|---|---|---|
| `stomata_model = 6` | jules_vegetation | — | CMax; needs `som_base_parm = 1` and `som_ci_search = 2` |
| `cmax_a_io` | pft_params | none (required > 0) | a, µmol m⁻² s⁻¹ MPa⁻² |
| `cmax_b_io` | pft_params | 0 | b, µmol m⁻² s⁻¹ MPa⁻¹ (Sabot et al. 2022: low influence, poorly constrained) |

- `som_profit_model = 4` (`cmax_profit_model`). The values 3 and 7 are left free for the supply-loss model on branch `supply_loss`, so the two can be merged.
- **Checks extended to 8:** can_rad_mod 7, segments, cuticular floor, gravity, supply limit, leaf T, vcmax_psi. `l_som_nsl` stays profit-max only.
- **Diagnostics:** `carbon_gain` = A/A_max on the feasible range; `hydraulic_cost` = Θ/A_max.

## Tests (FR-Pue, `runs/roses/FR_Pue/FR_Pue_dev/cmax_tests/`)

The base set is the calibrated profit max `calib_gce_rs_seg_lt_wave1_gpp1_ebc` run_0000 with `knl_io = 0.5`. d80824c uses knl for the two-leaf N weights; 0.5 reproduces the previous kpar-based weights.

1. **Regression:** model 4 on this exe gives a 2003 `run.out.nc` byte-identical to `global_change_ecology` 2e1e2ea (`pm2003_ref` vs `pm2003_cmaxexe`).
2. **Smoke, 2003, a = 1 / 3 / 10 (b = 0):**
   - no NaNs; closed at night; same run time as the profit max.

   | | profit max | a = 1 | a = 3 | a = 10 |
   |---|---|---|---|---|
   | TVeg (mm) | 352 | 269 | 200 | 118 |
   | GPP (g C m⁻²) | 1106 | 1157 | 1095 | 934 |
   | min ψ_leaf (MPa) | −5.8 | −3.2 | −1.3 | −0.6 |
   | July GPP | 90 | 145 | 131 | 96 |

   - Much more isohydric than the profit max: it saves water in spring and uses it in the 2003 summer.
   - **Winter (DJF) daytime gs still rises as light falls** (a = 3: 7.4 → 3.0 mm/s from < 50 to > 300 W m⁻² SWnet; profit max 20 → 8.7). Lower than the profit max, but the same shape. The cost depends only on ψ_leaf, which barely moves when E is small, so CMax opens when cost is low (as Wolf et al. 2016 state), and the dim-light/cold problem of the profit max is not solved by absolute units alone.
3. **Calibration** (in progress): `supply_loss_tests/emulator_sl`, ensemble `cmax_seg_lt_vfloor_wave1`.
   - Set up as the best profit max (`gce_rs_seg_lt_vfloor_wave1`), plus cmax_a free (log 0.3–30), b = 0.
   - Calibrated three ways: gpp1_ebc; root P50 fixed −2.41; ψ in the fit.

## Context: Sabot et al. (2022)

Sabot et al. (2022, JAMES; 12 schemes in a LSM framework, leaf-level data):
- CMax's a is constrainable, but b is not.
- CMax closes about 45 % too early at high D, and gs(ψ_l) splits into branches.
- It reverses the ψ_l diurnal cycle in the afternoon, and was not among the best schemes.
- ProfitMax, SOXopt and CGain matched the observed closure.
- CGain (Lu et al. 2020), A_n − ϖ (kmax − k)/kmax, is the absolute-carbon version of the profit max (a carbon scale ϖ in place of the instantaneous A_max). That would be the natural next option if the absolute-units idea is pursued, with the same caveat on winter.

# CGain (`stomata_model = 7`) and SOX_opt fixes (7 Oct 2026, same branch)

## CGain

- **Model** (Lu et al. 2020; Sabot et al. 2022, Eq. 15): maximise A_n − ϖ·(k_max − k(ψ_leaf))/k_max.
  - k is the whole-path −dE/dψ_leaf, the same k as the profit max.
  - k_max is the unstressed whole-path conductance (`k_path_zero_flow` at ψ = 0, soil link included).
  - ϖ is `cgain_varpi_io` (µmol m⁻² s⁻¹ per leaf area, required), scaled to the call's leaf area like CMax.
- `som_profit_model = 5`; bounded search; net A.
- **Smoke, 2003** (`cmax_tests/cgain2003_w*`):

  | ϖ | TVeg (mm) | min ψ_leaf (MPa) | winter daytime gs, dim → bright (mm/s) |
  |---|---|---|---|
  | 3 | 374 | −7.1 | 18 → 9 |
  | 10 | 318 | −5.1 | 13 → 6 |
  | 30 | 264 | −2.7 | 9 → 4 |
  | profit max | 352 | −5.8 | 20 → 9 |

  - It is the profit max in absolute carbon units. As with CMax, winter gs still rises in dim light: the cost is purely hydraulic, so it vanishes when E is small.

## SOX_opt (`stomata_model = 5`; was "SOX profit", model 6) fixes

1. **Whole-path k at ψ̄** = (ψ_root + ψ_leaf)/2, from `k_path_zero_flow` (every segment and the soil link at ψ̄), instead of the single PFT curve (`xylem_conductance_jls`). The old cost priced risk on the resistant stem curve only, while ψ_leaf and the feasibility test used the whole path. Both the flat and the bounded paths use the fix.
2. **Bounded search** with `som_ci_search = 2` (profit A_n·(k(ψ̄) − kcrit)/(k_ref − kcrit), unimodal); the flat grid is kept for `som_ci_search = 1`.

**FR-Pue 2003** (same parameters, profit-max calibrated set):

| | old flat | new flat | new bounded | profit max |
|---|---|---|---|---|
| TVeg (mm) | **509** | 360 | 360 | 352 |
| GPP (g C m⁻²) | 946 | 1094 | 1094 | 1106 |
| July GPP | **39** | 83 | 83 | 90 |
| winter midday ψ (MPa) | **−2.4** | −0.66 | −0.67 | −0.66 |
| run time (s) | 28 | 31 | **16** | 19 |

- **The old SOX failed like the unsegmented profit max:** spring over-transpiration (May 141 vs 77 mm), then a summer crash.
- **The fixed SOX behaves like the profit max**, as Sabot et al. (2022) also found.
- **Bounded vs flat:** within grid resolution of each other, and the bounded search is 1.9× faster.
- **Model 6 results change** for any configuration with segments or the soil link (i.e. all of them now).

Model 4 is still bit-for-bit with 2e1e2ea (`pm2003_v3exe`), and CMax is unchanged.
