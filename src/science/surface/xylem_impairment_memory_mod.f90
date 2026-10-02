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
          xylem_conductance_impaired_memory_stom_opt_jls,                      &
          leaf_psi_impaired_memory,                                            &
          update_xylem_impairment_memory,                                      &
          xylem_refit_weibull,                                                 &
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

! LAI and wood carbon of each PFT at the last update, for the leaf-area and
! growth recovery terms; -1 = not yet seen.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_lai_prev(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_wood_prev(:,:)

! Slow recovery (ximpair_rec_years > 0): loss of conductivity at the last
! damage, which sets the recovery rate, and the running mean (s-1) of the
! renewed fraction with its weight (for the start-up bias correction).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_plc_dam(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_renew_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_renew_wt(:,:)
! ximpair_npp_prev, ximpair_lock, ximpair_lai_prev, ximpair_wood_prev,
! ximpair_plc_dam, ximpair_renew_mean and ximpair_renew_wt are written to and
! read from dumps (ximpair_memory_alloc allocates them).

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
IF (.NOT. ALLOCATED(ximpair_lai_prev)) THEN
  ALLOCATE(ximpair_lai_prev(n_land_pts, npft))
  ximpair_lai_prev(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(ximpair_wood_prev)) THEN
  ALLOCATE(ximpair_wood_prev(n_land_pts, npft))
  ximpair_wood_prev(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(ximpair_plc_dam)) THEN
  ALLOCATE(ximpair_plc_dam(n_land_pts, npft))
  ximpair_plc_dam(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_renew_mean)) THEN
  ALLOCATE(ximpair_renew_mean(n_land_pts, npft))
  ximpair_renew_mean(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_renew_wt)) THEN
  ALLOCATE(ximpair_renew_wt(n_land_pts, npft))
  ximpair_renew_wt(:,:) = 0.0
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
!   - recovery: optional, from new leaf area (l_ximpair_rec_lai) and/or new
!               xylem grown from carbon gain (l_ximpair_rec_growth), plus
!               refilling (ximpair_tau_rec) and an annual reset
!               (ximpair_reset_mmdd).
!
! This is the target curve of Mackay et al. (2015, WRR, eqs 11-12; TREES):
! Kcav = Ksat f(psi_min) with the intact curve followed below psi_min. TREES
! then fits a new Weibull (b, c) to it; the kmax refit impairment model
! (pft_xylem_impairment_model = 4) does the same analytically (see
! xylem_refit_weibull), whereas this model uses the capped curve exactly.
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
! Xylem conductance on the capped curve for open points (stomatal
! optimisation layout).
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_impaired_memory_stom_opt_jls( pft,                &
                                                 n_water_potentials,           &
                                                 open_pnts,                    &
                                                 open_index,                   &
                                                 veg_index,                    &
                                                 land_pnts,                    &
                                                 water_potential,              &
                                                 kmax_ref,                     &
                                                 kcap,                         &
                                                 kcrit,                        &
                                                 conductance_b,                &
                                                 conductance_c,                &
                                               ! INTENT OUT
                                                 xylem_conductance             &
  )

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft, n_water_potentials, open_pnts, open_index(open_pnts), land_pnts,       &
  veg_index(land_pnts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(n_water_potentials, open_pnts)                               &
, kmax_ref(land_pnts), kcap(land_pnts), kcrit(land_pnts)                       &
, conductance_b(land_pnts), conductance_c(land_pnts)
                            ! See leaf_conductance_impaired_memory_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)

INTEGER :: i, l
REAL(KIND=real_jlslsm) ::                                                      &
  kmax_open(open_pnts), kcrit_open(open_pnts), b_open(open_pnts),             &
  c_open(open_pnts), kcap_open(open_pnts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_MEMORY_STOM_OPT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1, open_pnts
  l = veg_index(open_index(i))
  kmax_open(i)  = kmax_ref(l)
  kcrit_open(i) = kcrit(l)
  b_open(i)     = conductance_b(l)
  c_open(i)     = conductance_c(l)
  kcap_open(i)  = MAX(kcap(l), kcrit(l))
END DO

CALL xylem_conductance_jls( pft, n_water_potentials, open_pnts,                &
                            water_potential, kmax_open, kcrit_open, b_open,    &
                            c_open, xylem_conductance )
xylem_conductance = MIN(xylem_conductance,                                     &
                   SPREAD(kcap_open, DIM = 1, NCOPIES = n_water_potentials))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_impaired_memory_stom_opt_jls

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

INTEGER, PARAMETER :: max_nr_iter = 20
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Newton-Raphson convergence tolerance (Pa).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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
! Recovery (switches in jules_vegetation):
!   l_ximpair_rec_lai    - new leaf area (a rise in LAI) comes with undamaged
!                          xylem: k_cap is the leaf-area weighted mean of the
!                          existing cap and kmax for the new leaf area,
!                            k_cap' = (LAI_old k_cap + dLAI kmax) / LAI_new
!                          (no parameters; with prescribed LAI this captures
!                          e.g. new earlywood at leaf flush).
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
!
! Slow recovery (ximpair_rec_years > 0). Applied directly, the renewed
! fractions above recover the damage within ~1 year, as new leaf area and
! wood are compared only with the current (small) sapwood. Instead, new
! xylem replaces the damaged conduits over ximpair_rec_years years of growth
! (sapwood turnover; growth recovers over 3-5 years after drought, e.g.
! Anderegg et al. 2015, Kannenberg et al. 2019). The fraction renewed this
! timestep, f = f_lai + f_growth, is scaled by its running mean, <f>
! (e-folding time ximpair_rec_years, bias-corrected at the start), to the
! growth g = f / (<f> * 1 year) in units of a typical year's growth, and the
! loss of conductivity falls linearly with growth,
!   PLC' = MAX(PLC - PLC_dam * g / ximpair_rec_years, 0),
! with PLC_dam the loss at the last damage (reset each time damage raises
! PLC). So the loss recovers in ximpair_rec_years years of typical growth
! after the last damage, with the seasonal timing of leaf flush / growth,
! and more slowly after years of low growth.
!
! Also (per PFT, off by default): refilling with timescale ximpair_tau_rec
! while psi_x > ximpair_psi_refill, and an annual reset on
! ximpair_reset_mmdd.
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_memory ( n_land_pts                         &
,                                           n_open_pts                         &
,                                           open_index                         &
,                                           pft                                &
,                                           psi_leaf                           &
,                                           psi_root                           &
,                                           lai                                &
,                                           canht                              &
,                                           anetc                              &
,                                           kcap                               &
                                          )

USE pftparm, ONLY: kmax_pft, kcrit, conductance_b_pft, conductance_c_pft,      &
                   eta_sl, ximpair_psi_driver, ximpair_reset_mmdd,             &
                   ximpair_tau_rec, ximpair_psi_refill, ximpair_wood_alloc,    &
                   ximpair_growth_basis, ximpair_rec_years, a_wl, a_ws, b_wl
USE jules_vegetation_mod, ONLY: ximpair_driver_leaf, ximpair_driver_mean,      &
                                ximpair_driver_root, l_ximpair_rec_lai,        &
                                l_ximpair_rec_growth, l_triffid
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

REAL(KIND=real_jlslsm) ::                                                      &
  psi_x(n_land_pts), k_x(n_land_pts), kmax_pts(n_land_pts),                   &
  kcrit_pts(n_land_pts), b_pts(n_land_pts), c_pts(n_land_pts),                &
  f_renew(n_land_pts), c_new(n_land_pts), wood(n_land_pts),                  &
  f_new(n_land_pts), growth(n_land_pts), kcap_old(n_land_pts),                &
  plc(n_land_pts)
REAL(KIND=real_jlslsm) :: f_rec, a_mean
LOGICAL :: l_slow
                            ! Slow recovery (ximpair_rec_years > 0).

REAL(KIND=real_jlslsm), PARAMETER :: lai_min = 1.0e-3
                            ! Minimum LAI for the recovery terms.
REAL(KIND=real_jlslsm), PARAMETER :: c_per_mol_co2 = 12.0e-3
                            ! kg C per mol CO2.

! TRIFFID turnover rates (g_wood) are per 360-day year.
REAL(KIND=real_jlslsm), PARAMETER :: sec_per_trif_year = 360.0 * 86400.0
REAL(KIND=real_jlslsm), PARAMETER :: sec_per_year = 365.25 * 86400.0

                            ! LAI at the previous update, per PFT, for
                            ! l_ximpair_rec_lai. NOTE: not held in the dump,
                            ! so after a restart the first update sees no
                            ! change in LAI.

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
CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'ximpair_psi_driver should be leaf (1), mean (2) or root (3)')
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

IF (l_ximpair_rec_lai) THEN
  IF (.NOT. ALLOCATED(ximpair_lai_prev)) THEN
    ALLOCATE(ximpair_lai_prev(n_land_pts, npft))
    ximpair_lai_prev(:,:) = -1.0
  END IF
  DO l = 1, n_land_pts
    IF (ximpair_lai_prev(l,pft) >= 0.0 .AND. lai(l) > ximpair_lai_prev(l,pft) .AND.            &
        lai(l) > lai_min) THEN
      IF (l_slow) THEN
        f_new(l) = (lai(l) - ximpair_lai_prev(l,pft)) / lai(l)
      ELSE
        kcap(l) = ( MAX(ximpair_lai_prev(l,pft), 0.0) * kcap(l)                        &
                    + (lai(l) - ximpair_lai_prev(l,pft)) * kmax_pts(l) ) / lai(l)
      END IF
    END IF
    ximpair_lai_prev(l,pft) = lai(l)
  END DO
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
              * ( a_ws(pft) * eta_sl(pft) * MAX(canht(:), 0.0) / a_wl(pft) )   &
              ** (1.0 / (b_wl(pft) - 1.0))
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
  ! Renew that fraction of the conducting (live stem) wood.
  f_renew(:) = c_new(:)                                                        &
               / MAX(eta_sl(pft) * canht(:) * lai(:),                          &
                     eta_sl(pft) * MAX(canht(:), 1.0) * lai_min)
  f_renew(:) = MIN(f_renew(:), 1.0)
  IF (l_slow) THEN
    f_new(:) = f_new(:) + f_renew(:)
  ELSE
    kcap(:) = kcap(:) + (kmax_pts(:) - kcap(:)) * f_renew(:)
  END IF
END IF

IF (l_slow .AND. (l_ximpair_rec_lai .OR. l_ximpair_rec_growth)) THEN
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
  ! Linear recovery of the loss of conductivity.
  plc(:) = 1.0 - kcap(:) / kmax_pts(:)
  plc(:) = MAX(plc(:) - ximpair_plc_dam(:,pft) * growth(:)                     &
                        / ximpair_rec_years(pft), 0.0)
  kcap(:) = kmax_pts(:) * (1.0 - plc(:))
END IF

IF (ximpair_tau_rec(pft) > 0.0) THEN
  f_rec = 1.0 - EXP(-REAL(timestep_len) / (ximpair_tau_rec(pft) * 86400.0))
  WHERE (psi_x(:) > ximpair_psi_refill(pft))
    kcap(:) = kcap(:) + (kmax_pts(:) - kcap(:)) * f_rec
  END WHERE
END IF

IF (ximpair_reset_mmdd(pft) > 0) THEN
  IF (current_time%month * 100 + current_time%day == ximpair_reset_mmdd(pft)  &
      .AND. current_time%time < timestep_len) THEN
    kcap(:) = kmax_pts(:)
  END IF
END IF

!-----------------------------------------------------------------------------
! Damage: the cap can not exceed the intact conductance at the damage
! driver. Only points with open stomata (i.e. under tension) are damaged.
!-----------------------------------------------------------------------------
kcap_old(:) = kcap(:)
CALL leaf_conductance_jls( pft, n_land_pts, psi_x, kmax_pts, kcrit_pts,        &
                           b_pts, c_pts, k_x )
DO j = 1, n_open_pts
  l = open_index(j)
  kcap(l) = MIN(kcap(l), k_x(l))
END DO

kcap(:) = MAX(MIN(kcap(:), kmax_pts(:)), kcrit_pts(:))

! Slow recovery restarts from the loss at the latest damage.
IF (l_slow) THEN
  CALL ximpair_memory_alloc( n_land_pts )
  WHERE (kcap(:) < kcap_old(:))
    ximpair_plc_dam(:,pft) = 1.0 - kcap(:) / kmax_pts(:)
  END WHERE
END IF

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
ximpair_lock(:,pft) = 1.0 - kcap(:) / kmax_pts(:)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_memory

! ---------------------------------------------------------------------
! Refit the vulnerability curve to the capped curve, as in Mackay et al.
! (2015) / TREES: the refitted curve kcap * f(psi; b', c') matches
! MIN(kmax f(psi; b, c), kcap) relative to kcap at two points,
!   CW  (k = kmax exp(-(psi/b)^c)):   relative k = 1/e (defines b') and 0.12
!   SOX (k = kmax / (1+(psi/b)^c)):   relative k = 0.5 (defines b') and 0.12
! The CW b' is b (1 - ln r)^(1/c), r = kcap/kmax (as in the kmax model), but
! c' comes from the 0.12 point: the kmax model's c' = c (b/b')^c decreases
! with damage, whereas the capped curve, flat down to psi_min and then
! falling along the intact curve, needs a steeper (larger c) refit (Mackay
! et al. 2015, Fig. 2: c 4.08 -> 21.1).
! ---------------------------------------------------------------------
SUBROUTINE xylem_refit_weibull( pft, n, kcap, kmax, b, c, b_new, c_new )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance

INTEGER, INTENT(IN) :: pft, n
REAL(KIND=real_jlslsm), INTENT(IN) :: kcap(n), kmax, b, c
REAL(KIND=real_jlslsm), INTENT(OUT) :: b_new(n), c_new(n)

REAL(KIND=real_jlslsm) :: r(n), psi2(n)

r(:) = MIN(MAX(kcap(:) / kmax, 1.0e-6), 1.0)

IF (pft_conductance_model(pft) == CW_conductance) THEN
  b_new(:) = b * (1.0 - LOG(r(:)))**(1/c)
  psi2(:)  = b * (-LOG(0.12 * r(:)))**(1/c)
  c_new(:) = LOG(-LOG(0.12)) / LOG(psi2(:) / b_new(:))
ELSE
  b_new(:) = b * (2.0 / r(:) - 1.0)**(1/c)
  psi2(:)  = b * (1.0 / (0.12 * r(:)) - 1.0)**(1/c)
  c_new(:) = LOG(1.0 / 0.12 - 1.0) / LOG(psi2(:) / b_new(:))
END IF

END SUBROUTINE xylem_refit_weibull

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
