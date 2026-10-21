! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_memory_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_MEMORY_MOD'

PUBLIC :: ximpair_memory_alloc,                                                &
          leaf_conductance_impaired_memory_jls,                                &
          leaf_psi_impaired_memory,                                            &
          update_xylem_impairment_memory,                                      &
          ximpair_store_npp

PRIVATE :: k_intact, antideriv, psi_at_k

! NPP (kg C m-2 s-1) of each PFT from the previous timestep, for the growth
! recovery term with ximpair_growth_basis = 2. NPP is computed in sf_stom
! after the impairment update, so the update uses last step's value.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_npp_prev(:,:)

! Locked-in loss of conductivity, 1 - k_cap/kmax, of each PFT at each land
! point, kept for TRIFFID phenology (l_ximpair_leaf_loss). Not allocated
! until the first update (none during spin-up), and then read as zero
! where not yet set.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_lock(:,:)

! Wood carbon of each PFT at the last update, for the growth recovery term
! (ximpair_growth_basis = 3); -1 = not yet seen.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_wood_prev(:,:)

! Slow recovery (ximpair_rec_years > 0): the running mean (s-1) of the
! renewed fraction with its weight (for the start-up bias correction).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_renew_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_renew_wt(:,:)

! Stem conductance rebuilt by xylem regrown from the TRIFFID store
! (l_ximpair_rec_nsc), as a fraction of kmax, queued by triffid_nsc_store
! at the TRIFFID step and applied (then cleared) at the next update.
! MDK, 17/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_repair_pend(:,:)

! Two-pool stem memory (nsc_xylem_replace > 0): the share of the stem's
! conducting xylem rebuilt from the store, and the worst loss of
! conductance of that rebuilt pool since it was built. The old pool's
! embolised share is the rest of the stem's loss.
! MDK, 20/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_rebuilt(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_lock_new(:,:)

! Time since the last TRIFFID step with the root zone wetter than
! ximpair_psi_growth, as a fraction of the TRIFFID period (l_ximpair_rec_nsc;
! cleared at the TRIFFID step).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_grow_acc(:,:)

! Lowest intact stem conductance fraction the water status imposed since the
! last TRIFFID step (l_ximpair_rec_nsc; 1 = none; reset at the TRIFFID step).
! Loss the current conditions re-impose is not repaired: only
! PLC - (1 - ximpair_fmin_acc) is, so a healthy stand with small midday
! losses does not pay again and again for the same conduits, and repair
! waits for a drought to ease.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_fmin_acc(:,:)
! ximpair_npp_prev, ximpair_lock, ximpair_wood_prev,
! ximpair_renew_mean, ximpair_renew_wt, ximpair_repair_pend,
! ximpair_grow_acc and ximpair_fmin_acc are written to and read from dumps
! (ximpair_memory_alloc allocates them).

CONTAINS

! ---------------------------------------------------------------------
! Allocate the impairment memory (n_land_pts, npft) with its start values,
! if not done yet (first update, or a dump / initial-condition read).
! ---------------------------------------------------------------------
SUBROUTINE ximpair_memory_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(ximpair_npp_prev)) THEN
  ALLOCATE(ximpair_npp_prev(n_land_pts, npft))
  ximpair_npp_prev(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_wood_prev)) THEN
  ALLOCATE(ximpair_wood_prev(n_land_pts, npft))
  ximpair_wood_prev(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(ximpair_renew_mean)) THEN
  ALLOCATE(ximpair_renew_mean(n_land_pts, npft))
  ximpair_renew_mean(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_renew_wt)) THEN
  ALLOCATE(ximpair_renew_wt(n_land_pts, npft))
  ximpair_renew_wt(:,:) = 0.0
END IF
! MDK, 17/10/26
IF (.NOT. ALLOCATED(ximpair_repair_pend)) THEN
  ALLOCATE(ximpair_repair_pend(n_land_pts, npft))
  ximpair_repair_pend(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_rebuilt)) THEN
  ALLOCATE(ximpair_rebuilt(n_land_pts, npft), ximpair_lock_new(n_land_pts, npft))
  ximpair_rebuilt(:,:)  = 0.0
  ximpair_lock_new(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_grow_acc)) THEN
  ALLOCATE(ximpair_grow_acc(n_land_pts, npft))
  ximpair_grow_acc(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_fmin_acc)) THEN
  ALLOCATE(ximpair_fmin_acc(n_land_pts, npft))
  ximpair_fmin_acc(:,:) = 1.0
END IF

END SUBROUTINE ximpair_memory_alloc

! *********************************************************************
! Embolism memory xylem impairment model (pft_xylem_impairment_model = 3).
!
! Embolism is irreversible over short timescales (no refilling while the
! xylem is under tension), so the conductance lost at the most negative
! water potential experienced is "remembered": the vulnerability curve is
! capped at
!
!   k(psi) = MIN( k_intact(psi), k_cap ),  k_cap = kmax * (1 - PLC_mem)
!
! i.e. the dry end of the intact curve is unchanged and only the wet end is
! lowered (the no-refilling hysteresis of e.g. Sperry et al. 2017, Venturas
! et al. 2018). The state is k_cap, held in the impaired kmax state array
! (k_max / kmax_impaired_pft, leaf basis); conductance_b/c are not changed.
!
! Each timestep (see update_xylem_impairment_memory):
!   - damage:   k_cap = MIN(k_cap, k_intact(psi_x)), psi_x the damage driver
!               (leaf, mean of leaf and root zone, or root zone water
!               potential)
!   - recovery: optional, from new xylem grown from carbon gain
!               (l_ximpair_rec_growth), plus
!               refilling (ximpair_tau_rec) and an annual reset
!               (ximpair_reset_mmdd).
!
! This is the target curve of Mackay et al. (2015, WRR, eqs 11-12; TREES):
! Kcav = Ksat f(psi_min) with the intact curve followed below psi_min. TREES
! then fits a new Weibull (b, c) to it (the retired model 4 did the same
! analytically), whereas this model uses the capped curve exactly.
!
! With the cap, transpiration from psi_r to psi_l is
!   E = k_cap * (psi_r - MAX(psi_l, psi_s))               [capped part]
!     + A(psi_l) - A(MIN(psi_r, psi_s))  if psi_l < psi_s [intact part]
! where psi_s is the water potential at which k_intact = k_cap and
! A(psi_1) - A(psi_2) is the intact-curve integral from psi_1 to psi_2.
! *********************************************************************

! ---------------------------------------------------------------------
! Intact curve conductance for a vector of water potentials at one point.
! ---------------------------------------------------------------------
FUNCTION k_intact( pft, n, psi, kmax, b, c ) RESULT( k )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance

INTEGER, INTENT(IN) :: pft, n
REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n), kmax, b, c
REAL(KIND=real_jlslsm) :: k(n)

IF (pft_conductance_model(pft) == CW_conductance) THEN
  k = kmax * EXP( -(psi/b)**c )
ELSE
  k = kmax / (1 + (psi/b)**c)
END IF

END FUNCTION k_intact

! ---------------------------------------------------------------------
! A(psi) such that the intact-curve integral of k from psi_1 to psi_2 is
! A(psi_1) - A(psi_2) (see the psi_aprox_NR branches of leaf_psi_CW_jls /
! leaf_psi_SOX_jls).
! ---------------------------------------------------------------------
FUNCTION antideriv( pft, n, psi, kmax, b, c ) RESULT( a )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance
USE xylem_hydraulics_CW_jls_mod, ONLY: incomplete_gamma
USE xylem_hydraulics_SOX_jls_mod, ONLY: SOX_2F1

INTEGER, INTENT(IN) :: pft, n
REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n), kmax, b, c
REAL(KIND=real_jlslsm) :: a(n)

IF (pft_conductance_model(pft) == CW_conductance) THEN
  ! Note: the negative sign is present because b is negative.
  a = kmax * (-b/c) * incomplete_gamma(n, 1/c, (psi/b)**c)
ELSE
  a = - kmax * psi * SOX_2F1(n, 1 + 1/c, -(psi/b)**c)
END IF

END FUNCTION antideriv

! ---------------------------------------------------------------------
! Water potential at which the intact curve falls to r * kmax (0 < r < 1).
! ---------------------------------------------------------------------
FUNCTION psi_at_k( pft, r, b, c ) RESULT( psi )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: r, b, c
REAL(KIND=real_jlslsm) :: psi

IF (pft_conductance_model(pft) == CW_conductance) THEN
  psi = b * (-LOG(r))**(1/c)
ELSE
  psi = b * (1/r - 1)**(1/c)
END IF

END FUNCTION psi_at_k

! ---------------------------------------------------------------------
! Leaf conductance on the capped curve for all land points.
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_impaired_memory_jls( pft,                          &
                                                 land_pnts,                    &
                                                 water_potential,              &
                                                 kmax_ref,                     &
                                                 kcap,                         &
                                                 kcrit,                        &
                                                 conductance_b,                &
                                                 conductance_c,                &
                                               ! INTENT OUT
                                                 leaf_k                        &
  )

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, land_pnts
                            ! Number of land points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(land_pnts)                                                   &
                            ! Water potential for each land point (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance (m/s)
, kcap(land_pnts)                                                              &
                            ! Conductance cap, kmax * (1 - PLC_mem) (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point (m/s)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_IMPAIRED_MEMORY_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL leaf_conductance_jls( pft, land_pnts, water_potential, kmax_ref, kcrit,   &
                           conductance_b, conductance_c, leaf_k )
leaf_k = MIN(leaf_k, MAX(kcap, kcrit))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_memory_jls


! ---------------------------------------------------------------------
! Leaf water potential from transpiration on the capped curve, for open
! points (stomatal optimisation layout). The transpiration integral is
! evaluated exactly (see the module header) and solved by Newton-Raphson,
! independently of som_psi_aprox_method.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_impaired_memory( pft,                                      &
                                     n_e_leaf,                                 &
                                     land_pnts,                                &
                                     open_pnts,                                &
                                     veg_index,                                &
                                     open_index,                               &
                                     e_leaf,                                   &
                                     root_zone_psi,                            &
                                     kmax_ref,                                 &
                                     kcap,                                     &
                                     kcrit,                                    &
                                     conductance_b,                            &
                                     conductance_c,                            &
                                   ! INTENT OUT
                                     leaf_psi,                                 &
                                     leaf_k                                    &
  )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut
USE pftparm, ONLY: conductance_b_pft, conductance_c_pft
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_capped_lut_jls
! MDK, 18/10/26

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pnts, open_pnts, veg_index(land_pnts),                  &
  open_index(land_pnts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
, root_zone_psi(land_pnts)                                                     &
                            ! Water potential in the root zone (Pa)
, kmax_ref(land_pnts), kcap(land_pnts), kcrit(land_pnts)                       &
, conductance_b(land_pnts), conductance_c(land_pnts)
                            ! See leaf_conductance_impaired_memory_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: i, j, l

REAL(KIND=real_jlslsm) :: kmx, kc, bb, cc, psi_r, psi_s, psi_top, psi_floor
                            ! Per-point curve, cap, root zone psi, psi where
                            ! the intact curve meets the cap, top of the
                            ! intact-curve part of the path, and lower
                            ! bound for the Newton iteration.
REAL(KIND=real_jlslsm) :: e_capped
                            ! Transpiration carried by the capped part of the
                            ! path when psi_l < psi_s.
REAL(KIND=real_jlslsm) :: a_top(1)
                            ! A(psi_top).
REAL(KIND=real_jlslsm) :: psi_vec(n_e_leaf), e_vec(n_e_leaf),               &
                          k_vec(n_e_leaf), dpsi_vec(n_e_leaf)
LOGICAL :: l_capped(n_e_leaf)
                            ! Samples solved within the capped part.
LOGICAL :: l_pft_curve
                            ! som_psi_solver = 3: every open point on the
                            ! PFT curve (the table's).
                            ! MDK, 18/10/26

INTEGER, PARAMETER :: max_nr_iter = 20
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Newton-Raphson convergence tolerance (Pa).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! som_psi_solver = 3: the same capped curve inverted exactly from the PFT's
! supply table (leaf_psi_capped_lut_jls), when every open point is on the
! PFT curve (b, c), which the memory model does not change.
! MDK, 18/10/26
IF ( som_psi_solver == psi_solver_lut ) THEN
  l_pft_curve = .TRUE.
  DO j = 1, open_pnts
    l = veg_index(open_index(j))
    IF ( conductance_b(l) /= conductance_b_pft(pft) .OR.                      &
         conductance_c(l) /= conductance_c_pft(pft) ) l_pft_curve = .FALSE.
  END DO
  IF ( l_pft_curve ) THEN
    CALL leaf_psi_capped_lut_jls( pft, n_e_leaf, land_pnts, open_pnts,         &
                                  veg_index, open_index, e_leaf,               &
                                  root_zone_psi, kmax_ref, kcap, kcrit,        &
                                  leaf_psi, leaf_k )
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF

DO j = 1, open_pnts
  l = veg_index(open_index(j))

  kmx   = kmax_ref(l)
  kc    = MAX(MIN(kcap(l), kmx), kcrit(l))
  bb    = conductance_b(l)
  cc    = conductance_c(l)
  psi_r = root_zone_psi(l)
  psi_floor = psi_r - 5.0 * ABS(bb)

  ! Water potential where the intact curve meets the cap (+HUGE when there
  ! is no cap, so the whole path is on the intact curve).
  IF (kc < kmx * (1.0 - 1.0e-6)) THEN
    psi_s = psi_at_k(pft, kc / kmx, bb, cc)
  ELSE
    psi_s = HUGE(1.0_real_jlslsm)
  END IF

  !-------------------------------------------------------------------------
  ! Capped part: constant conductance kc from psi_r down to psi_s.
  !-------------------------------------------------------------------------
  IF (psi_r > psi_s) THEN
    psi_vec(:) = psi_r - e_leaf(:,j) / kc
    l_capped(:) = psi_vec(:) >= psi_s
    e_capped = kc * (psi_r - psi_s)
  ELSE
    l_capped(:) = .FALSE.
    e_capped = 0.0
  END IF

  !-------------------------------------------------------------------------
  ! Intact part: solve e = e_capped + A(psi_l) - A(psi_top) by Newton-
  ! Raphson (dE/dpsi_l = -k_intact(psi_l)). E is concave and decreasing in
  ! psi_l, so starting from psi_top the iterates approach the root
  ! monotonically from above.
  !-------------------------------------------------------------------------
  IF (.NOT. ALL(l_capped)) THEN
    psi_top = MIN(psi_r, psi_s)
    a_top = antideriv(pft, 1, (/ psi_top /), kmx, bb, cc)

    WHERE (.NOT. l_capped) psi_vec(:) = psi_top

    DO i = 1, max_nr_iter
      e_vec(:) = e_capped + antideriv(pft, n_e_leaf, psi_vec, kmx, bb, cc)    &
                 - a_top(1)
      k_vec(:) = MAX(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc),           &
                     TINY(1.0_real_jlslsm))
      dpsi_vec(:) = (e_vec(:) - e_leaf(:,j)) / k_vec(:)
      WHERE (l_capped) dpsi_vec(:) = 0.0

      ! Only the samples on the intact part of the path are iterated (and
      ! bounded above by psi_top): the capped samples were solved exactly
      ! above and lie above psi_top.
      WHERE (.NOT. l_capped)                                                   &
        psi_vec(:) = MIN(MAX(psi_vec(:) + dpsi_vec(:), psi_floor), psi_top)

      IF (MAXVAL(ABS(dpsi_vec(:))) < psi_tol) EXIT
    END DO
  END IF

  leaf_psi(:,j) = psi_vec(:)
  leaf_k(:,j) = MAX(MIN(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc), kc),   &
                    kcrit(l))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_impaired_memory

! ---------------------------------------------------------------------
! Update the conductance cap (leaf basis): recovery, then damage.
!
! Damage: k_cap = MIN(k_cap, k_intact(psi_x)) at points with open stomata,
! where psi_x is the damage-driving water potential (ximpair_psi_driver):
!   1: leaf, 2: mean of leaf and root zone, 3: root zone (~predawn, as used
!   for Kcav by Mackay et al. 2015).
!
! Recovery (switch in jules_vegetation). A leaf-area term (a rise in LAI
! renewing the cap) was retired on 2026-10-08: the net (prescribed) LAI rise
! is not leaf production, and new leaves renew leaf and twig xylem, not the
! stem (see runs/roses/FR_Pue/recovery_review_2026-10-08/REPORT.md).
!   l_ximpair_rec_growth - new conducting xylem grown from carbon gain
!                          replaces damaged xylem: the fraction renewed per
!                          timestep is
!                            C_new / C_sapwood
!                          with C_sapwood = eta_sl * canht * LAI the live
!                          stemwood carbon (kg C m-2) from the JULES
!                          allometry and C_new the new wood carbon this
!                          timestep, set by ximpair_growth_basis:
!                            1: ximpair_wood_alloc * A_net * dt (canopy net
!                               assimilation)
!                            2: ximpair_wood_alloc * NPP * dt (NPP of the
!                               previous timestep)
!                            3: TRIFFID gross wood production, the rise in
!                               allometric wood carbon plus g_wood turnover
!                               (needs l_triffid; no free parameter).
!   l_ximpair_rec_nsc    - carbon-limited repair (Trugman et al. 2018): the
!                          stem cap rises by the xylem rebuilt from the
!                          TRIFFID store, queued at the TRIFFID step in
!                          ximpair_repair_pend (triffid_nsc_mod). It replaces
!                          the stem's other recovery (l_ximpair_rec_growth,
!                          ximpair_tau_stem); the leaf segment keeps
!                          ximpair_tau_leaf, in calendar years.
!
! TODO (with TRIFFID): carbon-limited recovery. Recovery is now a fixed
! timescale (ximpair_rec_years) on a growth clock that only sets its timing.
! With l_triffid the repair could instead be limited by carbon: the renewed
! fraction taken from the new sapwood TRIFFID actually allocates (basis 3),
! without the running-mean normalisation, and/or repair given an explicit
! carbon cost drawn from the wood allocation, so that recovery slows when
! drought leaves little carbon to grow (cf. Paschalis et al. 2023, who
! attribute their overestimated FR-Pue legacy to the missing extra carbon
! for xylem repair). See NOTES_recovery_review.md.
!
! Slow recovery (ximpair_rec_years > 0). Applied directly, the renewed
! fractions above recover the damage within ~1 year, as new leaf area and
! wood are compared only with the current (small) sapwood. Instead, new
! xylem replaces the damaged conduits over ximpair_rec_years years of growth
! (sapwood turnover; growth recovers over 3-5 years after drought, e.g.
! Anderegg et al. 2015, Kannenberg et al. 2019). The fraction renewed this
! timestep, f = f_growth, is scaled by its running mean, <f>
! (e-folding time ximpair_rec_years, bias-corrected at the start), to the
! growth g = f / (<f> * 1 year) in units of a typical year's growth, and the
! loss of conductivity falls exponentially with growth,
!   PLC' = PLC * EXP(-2 g / ximpair_rec_years),
! e-folding over ximpair_rec_years / 2 typical years of growth: the same
! integrated loss after an isolated event as a linear recovery over
! ximpair_rec_years (the one-pool limit of renewing the sapwood with new,
! intact conduits). Recovery has the seasonal timing of growth, and is
! slower after years of low growth. Read as sapwood turnover: the conducting sapwood holds
! about ximpair_rec_years years of typical growth, and each year's new
! (intact) wood replaces an equal amount of old wood.
!
! Also (per PFT, off by default): refilling with timescale ximpair_tau_rec
! while psi_x > ximpair_psi_refill, and an annual reset on
! ximpair_reset_mmdd.
!
! Per-segment memory (l_ximpair_seg_memory, with l_som_plant_segments):
! kcap is the stem segment's cap and ximpair_kcap_leaf the leaf segment's
! (as a fraction of the intact leaf segment). Each is damaged by its own
! outlet water potential on its own curve (conductance_b/c_seg): the stem
! by psi_stem, the leaf by psi_leaf (ximpair_psi_driver is not used). Each
! recovers with its own e-folding time, in typical years of growth with
! the growth clock (ximpair_rec_years then only sets the window of the
! running mean) or in calendar years without it,
!   PLC_stem' = PLC_stem EXP(-g / ximpair_tau_stem)   (sapwood turnover)
!   PLC_leaf' = PLC_leaf EXP(-g / ximpair_tau_leaf)   (leaf/twig turnover),
! so a drought's damage to the distal, fast-renewed tissue fades within a
! few years while the stem keeps a long, small legacy. Refilling acts on
! each with its own water potential; the reset restores both.
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_memory ( n_land_pts                         &
,                                           n_open_pts                         &
,                                           open_index                         &
,                                           pft                                &
,                                           psi_leaf                           &
,                                           psi_root                           &
,                                           psi_stem                           &
,                                           lai                                &
,                                           canht                              &
,                                           anetc                              &
,                                           kcap                               &
                                          )

USE calc_c_comps_triffid_mod, ONLY: lai_bal_from_ht, live_stemwood   ! MDK, 20/10/26
USE pftparm, ONLY: kmax_pft, kcrit, conductance_b_pft, conductance_c_pft,      &
                   eta_sl, ximpair_psi_driver, ximpair_reset_mmdd,             &
                   ximpair_tau_rec, ximpair_psi_refill, ximpair_wood_alloc,    &
                   ximpair_growth_basis, ximpair_rec_years, a_wl, a_ws, b_wl, &
                   ximpair_psi_growth, ximpair_tau_stem, ximpair_tau_leaf,     &
                   conductance_b_seg, conductance_c_seg, nsc_xylem_replace
USE jules_vegetation_mod, ONLY: ximpair_driver_leaf, ximpair_driver_mean,      &
                                ximpair_driver_root, ximpair_driver_stem,      &
                                l_ximpair_rec_growth, l_triffid,               &
                                l_ximpair_seg_memory, l_ximpair_rec_nsc,       &
                                l_ximpair_rec_leaf_turn,                       &
                                triffid_period
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                   &
                                          ximpair_psi_stem,                    &
                                          ximpair_seg_state_alloc
USE trif, ONLY: g_wood
USE jules_surface_types_mod, ONLY: npft
USE model_time_mod, ONLY: current_time, timestep_len
USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  n_land_pts                                                                   &
                            ! Number of land points
, n_open_pts                                                                   &
                            ! Number of points with open stomata
, open_index(n_open_pts)                                                       &
                            ! Land point of each point with open stomata
, pft
                            ! Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  psi_leaf(n_land_pts)                                                         &
                            ! Leaf water potential (Pa)
, psi_root(n_land_pts)                                                         &
                            ! Root zone water potential (Pa)
, psi_stem(n_land_pts)                                                         &
                            ! Stem water potential (Pa; outlet of the stem
                            ! segment, = psi_root without segments)
, lai(n_land_pts)                                                              &
                            ! Leaf area index
, canht(n_land_pts)                                                            &
                            ! Canopy height (m)
, anetc(n_land_pts)
                            ! Canopy net photosynthesis (mol CO2 m-2 s-1)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  kcap(n_land_pts)
                            ! Conductance cap, kmax_pft * (1 - PLC_mem)
                            ! (leaf basis).

! Local variables
INTEGER :: j, l, errcode
REAL(KIND=real_jlslsm) :: rep_l, r_n, m_n, e_new, e_old, d_old, d_new
LOGICAL :: l_two
                            ! Two-pool stem memory (nsc_xylem_replace > 0).

REAL(KIND=real_jlslsm) ::                                                      &
  psi_x(n_land_pts), k_x(n_land_pts), kmax_pts(n_land_pts),                   &
  kcrit_pts(n_land_pts), b_pts(n_land_pts), c_pts(n_land_pts),                &
  f_renew(n_land_pts), c_new(n_land_pts), wood(n_land_pts),                  &
  f_new(n_land_pts), growth(n_land_pts),                                      &
  plc(n_land_pts)
REAL(KIND=real_jlslsm) :: f_rec, a_mean, frac_crit, store_days
LOGICAL :: l_slow
                            ! Slow recovery (ximpair_rec_years > 0).

REAL(KIND=real_jlslsm), PARAMETER :: lai_min = 1.0e-3
                            ! Minimum LAI for the recovery terms.
REAL(KIND=real_jlslsm), PARAMETER :: c_per_mol_co2 = 12.0e-3
                            ! kg C per mol CO2.

! TRIFFID turnover rates (g_wood) are per 360-day year.
REAL(KIND=real_jlslsm), PARAMETER :: sec_per_trif_year = 360.0 * 86400.0
REAL(KIND=real_jlslsm), PARAMETER :: sec_per_year = 365.25 * 86400.0


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Damage-driving water potential.
SELECT CASE (ximpair_psi_driver(pft))
CASE (ximpair_driver_leaf)
  psi_x(:) = psi_leaf(:)
CASE (ximpair_driver_mean)
  psi_x(:) = 0.5 * (psi_leaf(:) + psi_root(:))
CASE (ximpair_driver_root)
  psi_x(:) = psi_root(:)
CASE (ximpair_driver_stem)
  psi_x(:) = psi_stem(:)
CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'ximpair_psi_driver should be leaf (1), mean (2), root (3) ' // &
               'or stem (4)')
END SELECT

kmax_pts(:)  = kmax_pft(pft)
kcrit_pts(:) = kcrit(pft)
b_pts(:)     = conductance_b_pft(pft)
c_pts(:)     = conductance_c_pft(pft)

!-----------------------------------------------------------------------------
! Recovery. With slow recovery the renewed fractions are summed in f_new
! and applied below; otherwise each is applied directly.
!-----------------------------------------------------------------------------
l_slow = ximpair_rec_years(pft) > 0.0
f_new(:) = 0.0

! MDK, 17/10/26
IF (l_ximpair_rec_nsc) THEN
  ! Xylem rebuilt from the TRIFFID store since the last update.
  CALL ximpair_memory_alloc( n_land_pts )
  IF (nsc_xylem_replace(pft) > 0) THEN
    ! Two pools: the stem's loss is the old pool's embolised share plus
    ! the rebuilt pool's share times its own worst loss. The new xylem
    ! (repair_pend, as a share of the stem) is built embolism free.
    ! MDK, 20/10/26
    DO l = 1, n_land_pts
      rep_l = ximpair_repair_pend(l,pft)
      IF (rep_l <= 0.0) CYCLE
      r_n   = ximpair_rebuilt(l,pft)
      m_n   = ximpair_lock_new(l,pft)
      e_new = r_n * m_n
      e_old = MAX(1.0 - kcap(l) / kmax_pts(l) - e_new, 0.0)
      IF (nsc_xylem_replace(pft) == 1) THEN
        ! Embolised xylem first: the old pool's, then the rebuilt pool's.
        d_old = MIN(rep_l, e_old)
        d_new = MIN(rep_l - d_old, e_new)
        e_old = e_old - d_old
        e_new = e_new - d_new
        r_n   = MIN(r_n + d_old, 1.0)
      ELSE
        ! Proportional: the old pool is replaced as it is, intact and
        ! embolised alike.
        d_old = 0.0
        d_new = 0.0
        IF (r_n < 1.0) d_old = MIN(rep_l, 1.0 - r_n) * e_old / (1.0 - r_n)
        e_old = e_old - d_old
        r_n   = MIN(r_n + rep_l, 1.0)
      END IF
      m_n = 0.0
      IF (r_n > 0.0) m_n = MIN(e_new / r_n, 1.0)
      ! Once the rebuilt xylem is the whole stem it becomes the old pool
      ! (with its damage), so later repair acts on it again.
      IF (r_n >= 1.0 - 1.0e-6) THEN
        e_old = e_old + e_new
        e_new = 0.0
        r_n   = 0.0
        m_n   = 0.0
      END IF
      ximpair_rebuilt(l,pft)  = r_n
      ximpair_lock_new(l,pft) = m_n
      kcap(l) = kmax_pts(l) * (1.0 - e_old - e_new)
    END DO
  ELSE
    kcap(:) = MIN(kcap(:) + kmax_pts(:) * ximpair_repair_pend(:,pft),         &
                  kmax_pts(:))
  END IF
  ximpair_repair_pend(:,pft) = 0.0
  ! Count the time the root zone (~predawn) is wetter than
  ! ximpair_psi_growth, for the growth gate on the next repair, as a
  ! fraction of the store period (the TRIFFID period, or a day without
  ! TRIFFID).
  store_days = 1.0
  IF (l_triffid) store_days = REAL(triffid_period)
  WHERE (psi_root(:) >= ximpair_psi_growth(pft))                               &
    ximpair_grow_acc(:,pft) = ximpair_grow_acc(:,pft)                          &
                              + REAL(timestep_len) / (86400.0 * store_days)
END IF

IF (l_ximpair_rec_growth) THEN
  ! New (fully conductive) wood carbon this timestep, c_new (kg C m-2).
  SELECT CASE (ximpair_growth_basis(pft))
  CASE (1)
    ! A fixed fraction of canopy net photosynthesis (mol CO2 m-2 s-1).
    c_new(:) = ximpair_wood_alloc(pft) * c_per_mol_co2                         &
               * MAX(anetc(:), 0.0) * REAL(timestep_len)
  CASE (2)
    ! A fixed fraction of NPP (previous timestep).
    IF (.NOT. ALLOCATED(ximpair_npp_prev)) THEN
      c_new(:) = 0.0
    ELSE
      c_new(:) = ximpair_wood_alloc(pft) * MAX(ximpair_npp_prev(:,pft), 0.0)          &
                 * REAL(timestep_len)
    END IF
  CASE (3)
    ! TRIFFID gross wood production: the change in allometric wood carbon,
    !   wood = a_ws eta_sl h lai_bal,
    !   lai_bal = (a_ws eta_sl h / a_wl)^(1/(b_wl-1))
    ! (as diagnosed from canopy height in sf_stom), plus the wood turnover
    ! it replaces, g_wood wood. Wood shrinkage does not remove conduits.
    IF (.NOT. l_triffid) THEN
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
                   'ximpair_growth_basis = 3 (TRIFFID wood) needs l_triffid')
    END IF
    wood(:) = a_ws(pft) * eta_sl(pft) * canht(:)                               &
              * lai_bal_from_ht(pft, MAX(canht(:), 0.0))   ! MDK, 20/10/26
    IF (.NOT. ALLOCATED(ximpair_wood_prev)) THEN
      ALLOCATE(ximpair_wood_prev(n_land_pts, npft))
      ximpair_wood_prev(:,:) = -1.0
    END IF
    WHERE (ximpair_wood_prev(:,pft) < 0.0) ximpair_wood_prev(:,pft) = wood(:)
    c_new(:) = MAX(wood(:) - ximpair_wood_prev(:,pft)                                  &
                   + g_wood(pft) * wood(:) * REAL(timestep_len)                &
                   / sec_per_trif_year, 0.0)
    ximpair_wood_prev(:,pft) = wood(:)
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
                 'ximpair_growth_basis should be anet (1), npp (2) or '     // &
                 'triffid wood (3)')
  END SELECT
  ! No new conducting xylem while the root zone (~predawn) is drier than
  ! ximpair_psi_growth: stem growth stops in drought (sink limitation),
  ! although photosynthesis continues.
  WHERE (psi_root(:) < ximpair_psi_growth(pft)) c_new(:) = 0.0
  ! Renew that fraction of the conducting (live stem) wood.
  f_renew(:) = c_new(:)                                                        &
               / MAX(live_stemwood(pft, canht(:), lai(:)),                     &
                     live_stemwood(pft, MAX(canht(:), 1.0), lai_min))   ! MDK, 20/10/26
  f_renew(:) = MIN(f_renew(:), 1.0)
  IF (l_slow .OR. l_ximpair_seg_memory) THEN
    f_new(:) = f_new(:) + f_renew(:)
  ELSE
    kcap(:) = kcap(:) + (kmax_pts(:) - kcap(:)) * f_renew(:)
  END IF
END IF

IF (l_slow .AND. l_ximpair_rec_growth) THEN
  CALL ximpair_memory_alloc( n_land_pts )
  ! Running mean of the renewed fraction (s-1), and its weight.
  a_mean = MIN(REAL(timestep_len) / (ximpair_rec_years(pft) * sec_per_year), &
               1.0)
  ximpair_renew_mean(:,pft) = ximpair_renew_mean(:,pft)                        &
                + a_mean * (f_new(:) / REAL(timestep_len)                      &
                            - ximpair_renew_mean(:,pft))
  ximpair_renew_wt(:,pft) = ximpair_renew_wt(:,pft)                            &
                + a_mean * (1.0 - ximpair_renew_wt(:,pft))
  ! Growth this timestep in units of a typical year's growth.
  WHERE (ximpair_renew_mean(:,pft) > 0.0 .AND. ximpair_renew_wt(:,pft) > 0.0)
    growth(:) = f_new(:) * ximpair_renew_wt(:,pft)                             &
                / (ximpair_renew_mean(:,pft) * sec_per_year)
  ELSEWHERE
    growth(:) = 0.0
  END WHERE
  ! Exponential recovery of the loss of conductivity with growth.
  IF (.NOT. l_ximpair_seg_memory) THEN
    plc(:) = 1.0 - kcap(:) / kmax_pts(:)
    plc(:) = plc(:) * EXP(-2.0 * growth(:) / ximpair_rec_years(pft))
    kcap(:) = kmax_pts(:) * (1.0 - plc(:))
  END IF
ELSE
  ! Per-segment memory without the growth clock: calendar years.
  growth(:) = REAL(timestep_len) / sec_per_year
END IF

IF (l_ximpair_seg_memory) THEN
  CALL ximpair_seg_state_alloc( n_land_pts )
  ximpair_psi_stem(:,pft) = psi_stem(:)
  ! Each segment's loss falls with its own e-folding time (the stem's
  ! instead by the store's repair with l_ximpair_rec_nsc).
  ! MDK, 17/10/26
  IF (ximpair_tau_stem(pft) > 0.0 .AND. .NOT. l_ximpair_rec_nsc) THEN
    plc(:) = 1.0 - kcap(:) / kmax_pts(:)
    plc(:) = plc(:) * EXP(-growth(:) / ximpair_tau_stem(pft))
    kcap(:) = kmax_pts(:) * (1.0 - plc(:))
  END IF
  ! With leaf renewal (l_ximpair_rec_leaf_turn) the leaf segment recovers
  ! by new leaves instead (triffid_nsc_turgor, phenol).
  ! MDK, 18/10/26
  IF (ximpair_tau_leaf(pft) > 0.0 .AND. .NOT. l_ximpair_rec_leaf_turn) THEN
    ximpair_kcap_leaf(:,pft) = 1.0 - (1.0 - ximpair_kcap_leaf(:,pft))         &
                               * EXP(-growth(:) / ximpair_tau_leaf(pft))
  END IF
END IF

IF (ximpair_tau_rec(pft) > 0.0) THEN
  f_rec = 1.0 - EXP(-REAL(timestep_len) / (ximpair_tau_rec(pft) * 86400.0))
  IF (l_ximpair_seg_memory) THEN
    WHERE (psi_stem(:) > ximpair_psi_refill(pft))
      kcap(:) = kcap(:) + (kmax_pts(:) - kcap(:)) * f_rec
    END WHERE
    WHERE (psi_leaf(:) > ximpair_psi_refill(pft))
      ximpair_kcap_leaf(:,pft) = ximpair_kcap_leaf(:,pft)                      &
                                 + (1.0 - ximpair_kcap_leaf(:,pft)) * f_rec
    END WHERE
  ELSE
    WHERE (psi_x(:) > ximpair_psi_refill(pft))
      kcap(:) = kcap(:) + (kmax_pts(:) - kcap(:)) * f_rec
    END WHERE
  END IF
END IF

IF (ximpair_reset_mmdd(pft) > 0) THEN
  IF (current_time%month * 100 + current_time%day == ximpair_reset_mmdd(pft)  &
      .AND. current_time%time < timestep_len) THEN
    kcap(:) = kmax_pts(:)
    IF (l_ximpair_seg_memory) ximpair_kcap_leaf(:,pft) = 1.0
    IF (ALLOCATED(ximpair_rebuilt)) THEN
      ximpair_rebuilt(:,pft)  = 0.0
      ximpair_lock_new(:,pft) = 0.0
    END IF
  END IF
END IF

!-----------------------------------------------------------------------------
! Damage: the cap can not exceed the intact conductance at the damage
! driver. Only points with open stomata (i.e. under tension) are damaged.
!-----------------------------------------------------------------------------
l_two = l_ximpair_rec_nsc .AND. nsc_xylem_replace(pft) > 0
IF (l_ximpair_seg_memory) THEN
  ! Per segment: the stem by psi_stem on the stem curve, the leaf by
  ! psi_leaf on the leaf curve (the segments' intact Weibull curves).
  DO j = 1, n_open_pts
    l = open_index(j)
    IF (l_two) THEN
      CALL damage_two_pool( l, 1.0 - seg_f(psi_stem(l), 2) )
    ELSE
      kcap(l) = MIN(kcap(l), kmax_pts(l) * seg_f(psi_stem(l), 2))
    END IF
    ! MDK, 17/10/26
    IF (l_ximpair_rec_nsc)                                                     &
      ximpair_fmin_acc(l,pft) = MIN(ximpair_fmin_acc(l,pft),                   &
                                    seg_f(psi_stem(l), 2))
    ximpair_kcap_leaf(l,pft) = MIN(ximpair_kcap_leaf(l,pft),                   &
                                   seg_f(psi_leaf(l), 3))
  END DO
  frac_crit = kcrit(pft) / MAX(kmax_pft(pft), TINY(1.0_real_jlslsm))
  ximpair_kcap_leaf(:,pft) = MAX(MIN(ximpair_kcap_leaf(:,pft), 1.0), frac_crit)
ELSE
  CALL leaf_conductance_jls( pft, n_land_pts, psi_x, kmax_pts, kcrit_pts,      &
                             b_pts, c_pts, k_x )
  DO j = 1, n_open_pts
    l = open_index(j)
    IF (l_two) THEN
      CALL damage_two_pool( l, 1.0 - k_x(l) / kmax_pts(l) )
    ELSE
      kcap(l) = MIN(kcap(l), k_x(l))
    END IF
    ! MDK, 17/10/26
    IF (l_ximpair_rec_nsc)                                                     &
      ximpair_fmin_acc(l,pft) = MIN(ximpair_fmin_acc(l,pft),                   &
                                    k_x(l) / kmax_pts(l))
  END DO
END IF

kcap(:) = MAX(MIN(kcap(:), kmax_pts(:)), kcrit_pts(:))

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
ximpair_lock(:,pft) = 1.0 - kcap(:) / kmax_pts(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

  ! Two-pool damage at point l from the intact loss plc_x at the current
  ! water potential: the old pool's embolised share can rise to its share
  ! of plc_x, the rebuilt pool's worst loss to plc_x.
  ! MDK, 20/10/26
  SUBROUTINE damage_two_pool( l, plc_x )
  INTEGER, INTENT(IN) :: l
  REAL(KIND=real_jlslsm), INTENT(IN) :: plc_x
  r_n   = ximpair_rebuilt(l,pft)
  m_n   = MAX(ximpair_lock_new(l,pft), plc_x)
  e_old = MAX(1.0 - kcap(l) / kmax_pts(l) - r_n * ximpair_lock_new(l,pft),   &
              0.0)
  e_old = MAX(e_old, (1.0 - r_n) * plc_x)
  ximpair_lock_new(l,pft) = m_n
  kcap(l) = kmax_pts(l) * (1.0 - e_old - r_n * m_n)
  END SUBROUTINE damage_two_pool

  ! Intact relative conductance of segment iseg at water potential psi.
  FUNCTION seg_f( psi, iseg ) RESULT( f )
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  INTEGER, INTENT(IN) :: iseg
  REAL(KIND=real_jlslsm) :: f
  f = EXP(-( ABS(MIN(psi, 0.0) / conductance_b_seg(pft,iseg)) )               &
          **conductance_c_seg(pft,iseg))
  END FUNCTION seg_f

END SUBROUTINE update_xylem_impairment_memory

! *********************************************************************
! Keep this timestep's NPP of a PFT for the next timestep's growth
! recovery (ximpair_growth_basis = 2). Called from sf_stom.
! *********************************************************************
SUBROUTINE ximpair_store_npp ( n_land_pts, pft, npp )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts, pft
REAL(KIND=real_jlslsm), INTENT(IN) :: npp(n_land_pts)

IF (.NOT. ALLOCATED(ximpair_npp_prev)) THEN
  ALLOCATE(ximpair_npp_prev(n_land_pts, npft))
  ximpair_npp_prev(:,:) = 0.0
END IF
ximpair_npp_prev(:,pft) = npp(:)

END SUBROUTINE ximpair_store_npp

END MODULE xylem_impairment_memory_mod
