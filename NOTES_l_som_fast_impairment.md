# l_som_fast on the xylem_impairment branch: what to do

Working note, 2026-10-01. Written from `gs_opt_dev` commit `4c73d4b`
("Add l_som_fast stomatal-optimisation option; retire golden Ci search").
This branch (`xylem_impairment`, commit `823bbdf`) does not have it yet.

## What 4c73d4b adds (on gs_opt_dev)

- `l_som_fast` (`&jules_vegetation`, default `.false.`). With `.true.` and
  `leaf_flux_mod=2` it gives:
  - **Bounded Ci search** (`stom_opt_bounded_search` in `stom_opt_jls_mod.f90`).
    The feasible Ci range is always [ccp, ci_b], because gl, E and the loss of
    k all rise with Ci. The search finds ci_b by an Illinois root find, then
    runs golden-section on [ccp, ci_b] (`som_n_ci_golden_iter`). It works per
    point and exits early when the stomata are closed. With the table, the edge
    comes from a cap on gl, `g_cap = min(gl_max, E_crit*R*T/vpd)` with
    `E_crit = kmax*supply_lut_e_crit(...)`, so it needs photosynthesis only.
  - **Supply-function lookup table**, `som_psi_aprox_method=3`, which
    `l_som_fast` forces. It lives in
    `xylem_hydraulics_cumulative_weibull_jls_mod.f90`: `build_supply_lut`,
    `supply_lut_psi`, `supply_lut_f`, `supply_lut_e_crit`, `leaf_psi_lut_jls`.
    The table holds S(psi) = integral(f, psi, 0) once per PFT for f = k/kmax
    (Weibull or SOX), on 4001 uniform psi points. Solving
    S(psi_l) = S(psi_r) + E/kmax gives psi_leaf with no iteration.
- **Removed:** golden search, `som_ci_search_method`, `som_n_ci_prescan`.
- **Fixed:** the sign of the Newton update in `leaf_psi_SOX_jls`. It was
  `leaf_psi + (e - E)/k` and is now `leaf_psi - (e - E)/k`.
- **Validated at FR-Pue:** big-leaf, two-leaf (`can_rad_mod=7`) and
  multilayer (`can_rad_mod=6`); Farquhar and Collatz; Weibull and SOX.
  - Table vs Newton agree to 0.04% or better, for kmax = 0.2, 0.4 and 1.0.
  - `l_som_fast` vs a fine flat reference (8000 points, or 2000 for
    multilayer): gc within +0.3 to +2%. Flat-100 is out by -12 to -15%: its
    grid can't resolve the gl_max edge near ca at dawn and dusk.
  - Cost over profit max off (2000-2014 with spinup): +12%, against +28% for
    flat-100.

## Why it does NOT drop straight into this branch

The table is keyed by PFT and assumes the PFT's `conductance_b/c` are fixed.
On this branch:

1. **`pft_xylem_impairment_model` with kmax impairment**
   (`xylem_impairment_kmax_mod.f90`, `update_xylem_impairment_kmax`, around
   lines 414-455) recomputes **both b and c per land point every step**, so
   the curve still meets the impaired kmax at psi = 0:
   - CW: `b' = b (1 - ln(k'/kmax))^(1/c)`, `c' = c (b/b')^c`
   - SOX: `b' = b (2 kmax/k' - 1)^(1/c)`, `c' = c (kmax/k') (b'/b)^c`

   The hydraulics routines here take `conductance_b(land_pts)` and
   `conductance_c(land_pts)` as arrays, not `conductance_b(pft)`.
2. **Memory model** (`xylem_impairment_memory_mod.f90`,
   `pft_xylem_impairment_model = 3`) uses `k = MIN(k_intact(psi), k_cap)` and
   keeps b and c unchanged. **The existing table works here** with a
   two-piece integral (step 2 below).
3. **Whole-trunk model** (`xylem_impairment_whole_trunk_mod.f90`) has its own
   per-point Newton solvers: `leaf_psi_impaired_CW_jls` and
   `leaf_psi_impaired_SOX_jls`, using `incomplete_gamma_scalar` and
   `sox_2f1_scalar`. These are separate from `leaf_psi_jls`, so a table set
   there would not reach them. Their update sign is correct
   (`dpsi = -(E - e)/dE`, lines 650 and 853).
4. The branch's own `xylem_hydraulics_SOX_jls_mod.f90` (line 480) still
   has the **SOX sign bug**.

## Plan

1. **Merge or rebase.** Bring `gs_opt_dev` into this branch, or the reverse.
   Expect conflicts in `stom_opt_jls_mod.f90`, `sf_stom_jls_mod.F90` and the
   three `xylem_hydraulics_*` files: both sides change them heavily. Keep the
   SOX sign fix.
2. **Memory model (easy).** For each point, let
   psi_cap = `psi_at_k(k_cap/kmax)`. Above psi_cap (the wet end) k equals
   k_cap; below it the intact curve applies:
   `E(psi_l) = k_cap*(psi_r - max(psi_l, psi_cap))_+ + kmax*(S(min(psi_l,psi_cap)) - S(min(psi_r,psi_cap)))`
   (taking each part only where it applies). Invert piecewise: the linear part
   directly, the rest with `supply_lut_psi`. Also cap E_crit and `supply_lut_f`
   at k_cap.
3. **kmax impairment (b and c per point).** Make the table dimensionless.
   For Weibull, S(psi; b, c) = |b| * G_c(psi/b), with
   G_c(x) = integral(exp(-u^c), 0, x). SOX is the same with 1/(1+u^c).
   - b then needs no table: just rescale psi and S.
   - c needs either a 2-D table G(x, c) on a grid of c (about 30 values over
     the range impairment can produce, interpolated in c, since G is smooth
     in c), or a fall-back to the Newton solver for any point whose c differs
     from the PFT's c.
   - Pass per-point b and c into `supply_lut_psi`, `supply_lut_f` and
     `supply_lut_e_crit`, and into `stom_opt_bounded_search`'s `eval_ci`,
     which uses `supply_lut_f(pft, psi)` and the E_crit cap.
   - Check: impaired Newton vs impaired table at a few fixed impairment
     levels, and an FR-Pue 2003 drought run.
4. **Whole-trunk model.** Either route its psi solve through the table (steps
   2-3) or leave it on Newton with `l_som_fast`'s bounded search still in
   use. The search only needs `psi(E)`, `k(psi)` and, ideally, E_crit.
5. **Until this is done:** with `l_som_fast=.true.` and any impairment model
   on, use the Newton solver. One option: in `check_jules_vegetation`, keep
   `som_psi_aprox_method` unchanged (don't force 3) when impairment is
   active, or stop with an error, and make `stom_opt_bounded_search` use its
   Newton path (`l_lut = .false.`). The bounded search itself is still valid.

## How it was tested (to repeat after the merge)

- Base the test on the `FR_Pue_gs_opt` namelists:
  - Change `som_n_sample`.
  - Add `l_som_fast`.
  - Shorten the run to 2003-2004 with `max_spinup_cycles=0` for quick checks.
- Compare each method against flat with `som_n_sample=8000` and
  `som_psi_aprox_method=3` as the reference.
- For timing, build with `JULES_BUILD=fast` out of tree and set
  `output_main_run=.false.`. Serial `fcm make` only: `-j` deadlocks here.

## Status (2026-10-01, after merge f7cb920)

- Step 1 done: gs_opt_dev (0bf2e9d) merged into xylem_impairment. The SOX
  sign fix is kept (`- (e_leaf - e_leaf_current)`).
- Step 5 done as the error option: `check_jules_pftparm` stops with
  `l_som_fast` (and with DESICA, `stomata_model = 4`) when any
  `pft_xylem_impairment_model /= 0`. The bounded search passes the PFT's
  intact curve to `leaf_psi_jls`.
- Steps 2-4 (table for the impaired curves) are not started.
