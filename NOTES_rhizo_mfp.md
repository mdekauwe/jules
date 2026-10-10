# rhizo_mfp: nonlinear (matric flux potential) rhizosphere — design check, 2026-10-10

Branch `rhizo_mfp` off gce 398a05c (worktree git/rhizo_mfp). No code changes yet.

## Pre-check: the -2.1 MPa root P50 came after l_som_rhizo_series
- l_som_rhizo_series: eec84b5 (10-05 21:01), always on for profit max in ccee572 (21:40).
- First root-P50-free ensemble with it: gce_rs_seg_wave1, started 10-05 22:38. Every calib_* with p50_root is later.
- Pre-series (gce_pm_root_wave1, 10-03, root supply + root_psi_crit): best -3.85 (NROY -7.5 to -2.2).
- Post-series best root P50 tracks the LE target: raw LE -2.0 to -2.4 (bound); EBC LE mostly -3.0 to -4.2.

## Units (verified in code)
- soil_k: kg m-2 s-1 per unit head gradient (hyd_con_ic on sthu; vG: Se = sthu, theta_r = 0,
  Mualem L = 0.5, linear K clamps for Se < 0.05 and > 0.95; b = 1/(n-1), sathh = 1/alpha).
- soil_to_root_k = B_i * soil_k, per m2 ground per m of head (kg m-3 s-1), with
  B_i = 4 pi / ln(rho_r dz / (f m)) * f m / (rho_r pi r^2)  [m-1], steady-state cylinder (Bonan A23).
- ksr = sum soil_to_root_k / (rho_w g M_h2o): mol m-2 s-1 Pa-1 ground = kmax_pft * LAI units.
- MFP: Phi(h) = int_h K dh' (kg m-1 s-1); E_i (kg m-2 s-1) = B_i [Phi(h_soil) - Phi(h_root)];
  dE/dpsi_root = sum B_i K_i(psi_root) / (rho_w g) -> / M_h2o for mol.
- Steady-rate geometry (Schroeder 2008 / de Jong van Lier 2008):
  B = 2 pi L / (rho^2 ln rho / (rho^2 - 1) - 1/2) vs 2 pi L / ln rho; at FR-Pue rho = r_s/r_r = 20-31, +17-20 %.
- Current code is internally inconsistent on HR: the solver's psi_root = psi_rz - E/ksr with
  psi_rz = sum(k psi)/sum(k) is the signed balance (implicit redistribution); physiol's weights clip
  q_i >= 0.

## Offline size check (rhizo_size_check_offline.py; best fit calib_gce_rs_seg_lt_vfloor_cca_..._k63, run_0000)
Offline psi_s* reproduces the model's psi_root_zone_pft output (e.g. -1.58 vs -1.58), so
geometry and soil functions match the code. Midday (11-14 UTC) medians, psi MPa,
k mmol m-2 s-1 MPa-1 ground:

| case | psi_s* | psi_r lin | psi_r MFP | soil k lin | soil k' MFP | kplant | E mm/h |
|---|---|---|---|---|---|---|---|
| 2003 JJA | -1.58 | -1.61 | -1.78 | 42 | 33 | 0.36 | 0.086 |
| 2003 JJA driest 5 % | -3.84 | -3.88 | -4.13 | 4.4 | 4.0 | 0.15 | 0.012 |
| 2006 JJA | -2.29 | -2.34 | -2.60 | 14 | 13 | 0.30 | 0.040 |
| 2006 JJA driest 5 % | -3.94 | -3.97 | -4.23 | 4.2 | 3.8 | 0.14 | 0.011 |
| wet spring | ~-0.03 | same | same | 1e5-1e6 | same | 0.56 | 0.12-0.15 |

- With Bonan geometry and root mass = lma * lai_bal (0.41 kg m-2, ~4400 m m-2 root length), the
  soil link is 30-100x the plant even in the driest middays; soil supply is 40-60x E. MFP moves
  psi_root by <= 0.3 MPa. It cannot stand in for a -2.1 MPa root P50.
- B x 0.1: MFP adds 0.3-0.6 MPa over linear and halves k'; soil k' ~ 2-5x plant.
- B x 0.03: soil supply/E ~ 1.3-2; psi_root collapses to -3.5 to -7.4 MPa (Carminati-type).
- So the nonlinearity only matters if the effective root-soil conductance is ~10-30x lower than
  Bonan's: absorptive fraction of root mass, stones (75-90 % at Puechabon; the fine-earth K is
  applied to the whole layer), root-soil contact/gaps.

## Implementation (2026-10-10, uncommitted)
- Switches (jules_vegetation): `l_som_rhizo_mfp` (default F), `l_som_rhizo_hr` (default F; only
  with mfp). mfp needs l_som_rhizo_series (trapped in check_jules_vegetation).
- New `src/science/surface/rhizo_mfp_mod.F90`: per-PFT layer state (c_i = B_i K_sat sathh / m_h2o,
  G at bulk soil, zero-flow root h0), G(x = h/sathh) tables per distinct b (JULES vG / BC K(Se) with
  hyd_con clamps, Simpson-integrated, cubic Hermite in ln x; K from the interpolant's slope),
  safeguarded Newton for the root node. DOUBLE precision inside (real_jlslsm is real_32: the
  difference G_s - G_root is 1e-3..1e-6 of G_s in wet soil, and single precision never converged).
- Call sites: set_ksr_path (rz_share), k_path_zero_flow (k' at h0), leaf_psi_segments_jls,
  leaf_psi_segments_lut_jls, single-curve wrapper in xylem_hydraulics_jls_mod (k' per sample),
  physiol (state after smc_ext; psi_root_zone_pft := zero-flow root psi; extraction weights =
  layer uptakes at el_pft), DESICA (chord conductance at the previous step's root uptake).
- root_psi_crit: unchanged (already unused under l_som_rhizo_series).
- Not done: steady-rate geometry option; the two-leaf paths still each solve their own root node
  with E_path/share (as the linear link).

## Tests (FR-Pue, cca k63 best fit namelists; runs/roses/FR_Pue/FR_Pue_dev/rhizo_mfp_tests)
1. Off vs gce 398a05c frozen exe: BIT-IDENTICAL, run.out.nc + all 23 dumps (v3 and v4 exes).
2. Wet MAM 2004/2010: GPP equal to 3 decimals, LE within 0.02 %, psi within 0.001 MPa.
3. Dry-downs, JJA means (off / mfp / mfp_hr):
   2003 GPP 3.539/3.555/3.495, LE 36.1/35.9/35.6, psi_md -3.42/-3.22/-3.43, psi_soil -2.03/-1.88/-2.07.
   Closure dates within 1 d; closed-daytime share unchanged (0.5-0.8 %).
   vs obs 2000-14: psi_pd RMSE 0.574/0.555/0.532, psi_md 0.594/0.613/0.570, psi_md bias
   +0.14/+0.24/+0.12; GPP 1.046/1.053/1.040; LE raw 12.98/13.02/12.93.
   mfp (no HR) effect is mostly the zero-flow root at the wettest layer (wetter psi by
   0.15-0.2 MPa in drought); the K-nonlinearity itself (mfp_hr) is ~0.02 MPa, GPP -1 %.
4. Root P50 (no refit; other parameters from the -2.30 fit), off / mfp:
   -4.9 MPa: GPP RMSE 1.38/1.38, LE raw 19.9/20.1, psi_pd bias -1.92/-1.50, psi_md -1.95/-1.50.
   -6.0 MPa: GPP 1.50/1.49, LE raw 22.6/22.8, psi_pd bias -2.34/-1.86.
   Resistant roots raise wet-spring LE 38.6 -> 55.5 (MAM 2004) and dry the soil: the vulnerable
   root mainly throttles E in WET conditions (midday psi ~ -1.5 MPa), where the rhizosphere has
   no leverage. mfp recovers only ~0.4-0.5 MPa of the ~2 MPa psi bias. A fair test needs a refit
   (kmax / vsl free) with mfp on: proposed as an emulator ensemble, not run.
5. Cost: 182 M root solves per 2000-14 run, mean 2.35 Newton iterations (max 18; HR max 46).
   Wall time mfp/off = 1.48 (v4, K from the table slope; v3 with vG powers was 2.8).

## 2026-10-10 later: redistribution on by default; recalibration
- `l_som_rhizo_hr` default .TRUE. (closest to the linear link's implicit assumption; best psi fit).
  The "ignored without mfp" warning was removed (it would fire on every run).
- Exe v5 frozen as emulator/exe/jules_rhizo_mfp_v5_hr.exe. Ensemble
  gce_rs_seg_lt_vfloor_cca_mfp_wave1 (= cca design and seed, l_som_rhizo_mfp on), run by
  emulator/run_cca_mfp.sh (3 workers; calibrate k63 and k73 as run_cca.sh). Comparison:
  plot_rhizo_mfp_compare.py, compare_sap_fits.py (mfp fits added).

## Recalibration result (cca_mfp vs cca, 2026-10-10 06:24)
- v5 switch-off: BIT-IDENTICAL to gce 398a05c (24 files).
- Root P50 unchanged: k63 best -2.29 (cca -2.30), NROY -4.28..-2.17 (cca -4.30..-2.18); k73 -2.00 both.
  Other parameters the same except sf_vcmax (poorly constrained: 6.4 vs 4.9 k63; 7.8 vs 6.0 k73).
- Skill (k63, cca_mfp vs cca): GPP RMSE cal 1.08 vs 1.06, val 1.05 vs 1.03; LE cal 12.55 vs 12.52,
  val 13.43 vs 13.50; MAM LE bias 5.2 vs 5.4; T/ET 0.760 both; psi_md bias cal +0.15 vs +0.14,
  r 0.761 vs 0.746; psi_pd bias +0.20 vs +0.18; sap decline MAE 0.121 vs 0.115 (val 0.110 both);
  PLC max 10.7 vs 11.8 %. I.e. a wash: no measurable gain or loss at FR-Pue.
- Figure: emulator/plots/rhizo_mfp_vs_cca.png. Table: emulator/compare_sap_fits.csv.
- Both fits share the same structural errors: spring LE too high (~+5 W m-2) and deep-drought LE
  (Aug 2003/2006) ~5 W m-2 vs ~20 observed.
