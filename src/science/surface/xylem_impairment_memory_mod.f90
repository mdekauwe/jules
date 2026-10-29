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

PRIVATE

PUBLIC :: ximpair_memory_alloc, ximpair_dilute_area,                           &
          leaf_conductance_impaired_memory_jls,                                &
          psi_leaf_impaired_memory,                                            &
          update_xylem_impairment_memory

PRIVATE :: k_intact, antideriv, psi_at_k

! Locked-in loss of conductivity, 1 - k_cap/kmax, of each PFT at each land
! point, kept for TRIFFID phenology (l_ximpair_leaf_loss). Not allocated
! until the first update (none during spin-up), and then read as zero
! where not yet set.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_lock(:,:)

! Stem recovery by ring renewal (l_ximpair_rec_nsc): the share of the
! stem's conductance renewed by the wood TRIFFID grew, queued at the
! TRIFFID step by triffid_nsc_ring_renewal and applied (then cleared) at
! the next update, where it dilutes the stem's loss.
!
! MDK, 26/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_repair_pend(:,:)

! ximpair_lock, ximpair_repair_pend and ximpair_area_w (with
! l_veg_compete) are written to and read from dumps (ximpair_memory_alloc
! allocates them). ximpair_npp_prev, ximpair_wood_prev, ximpair_renew_mean
! and ximpair_renew_wt went with the time-based recovery; ximpair_rebuilt,
! ximpair_lock_new (the two-pool memory) and ximpair_grow_acc (the growth
! gate) with the store-paid repair; ximpair_fmin_acc (the lowest intact
! stem fraction imposed) with the W_max hold of the TRIFFID store.
!
! MDK, 29/10/26

! Factor (0-1) by which TRIFFID's gain of PFT area has diluted the stem
! damage (ximpair_dilute_area) since the last update, applied to the
! impaired kmax (kmax_impaired_pft, held by the surface code) at the next
! update and then reset to 1. In the dump, so a restart between a TRIFFID
! step and the next update is exact (allocated by ximpair_memory_alloc).
!
! MDK, 25/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_area_w(:,:)

CONTAINS

! ---------------------------------------------------------------------
! The memory model's damage after TRIFFID's change of PFT area (lotka,
! with l_veg_compete). Area gained is new growth with intact xylem, so the
! damage of the PFT is spread over its new area: each loss of conductance
! is multiplied by w = frac_old / frac_new (0 for a PFT that had no area),
! as triffid_nsc_area spreads the store. Area lost keeps the damage of the
! area that is left. Diluted:
! - the stem's locked loss (ximpair_lock; the impaired kmax at the next
!   update, through ximpair_area_w);
! - the leaf segment's loss (1 - ximpair_kcap_leaf).
! The queued ring share (ximpair_repair_pend) is a multiplicative factor
! on the loss, L -> L (1 - rho), the same for old and new cover, so it is
! not diluted: scaling it too would give L w (1 - rho w) instead of the
! area-weighted L w (1 - rho).
! With hydraulic mortality (m_kill, the fraction of the old cover killed in
! the step) the trees killed carry the mean damage and the area regrown in
! their place is intact too, so w = frac_old (1 - m) / frac_new wherever
! the cover is above the survivors', even if the net change is 0 or less.
!
! (The survivors' cover frac_old (1 - m) is approximate where lotka does not
! remove exactly m frac_old: the old lotka without l_ht_compete, and the
! frac_seed floor. MDK, 25/10/26)
!
! MDK, 26/10/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_dilute_area( land_pts, trif_pts, trif_index, frac, dfrac,   &
                                m_kill )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype
USE jules_vegetation_mod, ONLY: xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: land_pts, trif_pts, trif_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  frac(land_pts,ntype),                                                        &
      ! Areal fractions after the TRIFFID step.
  dfrac(land_pts,npft)
      ! Change in the areal fractions over the step.
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: m_kill(land_pts,npft)
      ! Fraction of each PFT's old cover killed in the step (hydraulic
      ! mortality, l_trif_hydmort).
      !
      ! MDK, 25/10/26

INTEGER :: n, t, l
REAL(KIND=real_jlslsm) :: w, m, frac_surv

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XIMPAIR_DILUTE_AREA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( ALLOCATED(ximpair_lock) ) THEN
  CALL ximpair_memory_alloc( land_pts )
  ! Natural PFTs only: lotka changes their cover (crop cover is set
  ! elsewhere). MDK, 25/10/26
  DO n = 1,nnpft
    IF ( pft_xylem_impairment_model(n) /= xylem_impairment_memory ) CYCLE
    DO t = 1,trif_pts
      l = trif_index(t)
      IF ( frac(l,n) <= 0.0 ) CYCLE
      m = 0.0
      IF ( PRESENT(m_kill) ) m = m_kill(l,n)
      frac_surv = MAX(frac(l,n) - dfrac(l,n), 0.0) * (1.0 - m)
      IF ( m > 0.0 ) THEN
        IF ( frac(l,n) <= frac_surv ) CYCLE
      ELSE
        IF ( dfrac(l,n) <= 0.0 ) CYCLE
      END IF
      w = frac_surv / frac(l,n)
      ximpair_area_w(l,n)      = ximpair_area_w(l,n) * w
      ximpair_lock(l,n)        = ximpair_lock(l,n) * w
      IF ( ALLOCATED(ximpair_kcap_leaf) )                                      &
        ximpair_kcap_leaf(l,n) = 1.0 - (1.0 - ximpair_kcap_leaf(l,n)) * w
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ximpair_dilute_area

! ---------------------------------------------------------------------
! Allocate the impairment memory (n_land_pts, npft) with its start values,
! if not done yet (first update, or a dump / initial-condition read).
! ---------------------------------------------------------------------
SUBROUTINE ximpair_memory_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XIMPAIR_MEMORY_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
! MDK, 17/10/26
IF (.NOT. ALLOCATED(ximpair_repair_pend)) THEN
  ALLOCATE(ximpair_repair_pend(n_land_pts, npft))
  ximpair_repair_pend(:,:) = 0.0
END IF
! MDK, 25/10/26
IF (.NOT. ALLOCATED(ximpair_area_w)) THEN
  ALLOCATE(ximpair_area_w(n_land_pts, npft))
  ximpair_area_w(:,:) = 1.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

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
!   - recovery: optional, ring renewal of the stem by TRIFFID wood
!               growth (l_ximpair_rec_nsc), leaf turnover of the leaf
!               segment (ximpair_tau_leaf), and an annual reset
!               (ximpair_reset_mmdd). (The growth and refilling
!               recovery were retired 23/10/26, and the store-paid stem
!               repair replaced by ring renewal 26/10/26.)
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
! A(psi_1) - A(psi_2) (see the Newton-Raphson branches, psi_solver_newton,
! of psi_leaf_CW_jls / psi_leaf_SOX_jls).
!
! MDK, 24/10/26
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
                            ! Unimpaired maximum xylem conductance
                            ! (mol m-2 s-1 Pa-1)
, kcap(land_pnts)                                                              &
                            ! Conductance cap, kmax * (1 - PLC_mem)
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='LEAF_CONDUCTANCE_IMPAIRED_MEMORY_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL leaf_conductance_jls( pft, land_pnts, water_potential, kmax_ref, kcrit,   &
                           conductance_b, conductance_c, leaf_k )
leaf_k = MIN(leaf_k, MAX(kcap, kcrit))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_memory_jls


! ---------------------------------------------------------------------
! Leaf water potential from transpiration on the capped curve, for open
! points (stomatal optimisation layout). The transpiration integral is
! evaluated exactly (see the module header) and solved by Newton-Raphson
! (som_psi_solver = 2). With som_psi_solver = 3 the capped curve is
! inverted from the PFT's supply table instead.
!
! MDK, 25/10/26
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_impaired_memory( pft,                                      &
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
USE supply_table_mod, ONLY: psi_leaf_capped_lut_jls
! MDK, 18/10/26

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pnts, open_pnts, veg_index(land_pnts),                   &
  open_index(land_pnts)
                            ! See psi_leaf_impaired_jls.
                            !
                            ! MDK, 24/10/26

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 leaf s-1)
                            !
                            ! MDK, 24/10/26
, root_zone_psi(land_pnts)                                                     &
                            ! Water potential in the root zone (Pa)
, kmax_ref(land_pnts), kcap(land_pnts), kcrit(land_pnts)                       &
, conductance_b(land_pnts), conductance_c(land_pnts)
                            ! See leaf_conductance_impaired_memory_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

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
REAL(KIND=real_jlslsm) :: psi_vec(n_e_leaf), e_vec(n_e_leaf),                  &
                          k_vec(n_e_leaf), dpsi_vec(n_e_leaf)
LOGICAL :: l_capped(n_e_leaf)
                            ! Samples solved within the capped part.
LOGICAL :: l_nr_active(n_e_leaf)
                            ! Newton: samples still iterating (each stops
                            ! on its own convergence, so the answer for one
                            ! E does not depend on the rest of the call).
                            !
                            ! MDK, 24/10/26

INTEGER, PARAMETER :: max_nr_iter = 20
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Newton-Raphson convergence tolerance (Pa).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_IMPAIRED_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! som_psi_solver = 3: the same capped curve inverted exactly from the PFT's
! supply table (psi_leaf_capped_lut_jls). Every point is on the PFT curve
! (b, c), which the memory model does not change.
!
! MDK, 24/10/26
IF ( som_psi_solver == psi_solver_lut ) THEN
  CALL psi_leaf_capped_lut_jls( pft, n_e_leaf, land_pnts, open_pnts,           &
                                veg_index, open_index, e_leaf,                 &
                                root_zone_psi, kmax_ref, kcap, kcrit,          &
                                leaf_psi, leaf_k )
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
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

    l_nr_active(:) = .NOT. l_capped(:)
    DO i = 1, max_nr_iter
      e_vec(:) = e_capped + antideriv(pft, n_e_leaf, psi_vec, kmx, bb, cc)     &
                 - a_top(1)
      k_vec(:) = MAX(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc),            &
                     TINY(1.0_real_jlslsm))
      dpsi_vec(:) = (e_vec(:) - e_leaf(:,j)) / k_vec(:)
      WHERE (l_capped) dpsi_vec(:) = 0.0

      ! Only the samples on the intact part of the path are iterated (and
      ! bounded above by psi_top): the capped samples were solved exactly
      ! above and lie above psi_top.
      ! Each sample stops when its own step is below psi_tol.
      WHERE (l_nr_active)                                                      &
        psi_vec(:) = MIN(MAX(psi_vec(:) + dpsi_vec(:), psi_floor), psi_top)
      WHERE (ABS(dpsi_vec(:)) < psi_tol) l_nr_active(:) = .FALSE.
      IF ( .NOT. ANY(l_nr_active(:)) ) EXIT
    END DO
  END IF

  leaf_psi(:,j) = psi_vec(:)
  leaf_k(:,j) = MAX(MIN(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc), kc),    &
                    kcrit(l))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_impaired_memory

! ---------------------------------------------------------------------
! Update the conductance cap (leaf basis): recovery, then damage.
!
! Damage: k_cap = MIN(k_cap, k_intact(psi_x)) at points with open stomata,
! where psi_x is the damage-driving water potential (ximpair_psi_driver):
!   1: leaf, 2: mean of leaf and root zone, 3: root zone (~predawn, as used
!   for Kcav by Mackay et al. 2015), 4: stem (leaf-side end of the stem
!   segment).
!
! Recovery:
!   l_ximpair_rec_nsc - ring renewal: the stem's loss is diluted by the
!                       share of the stem renewed by TRIFFID's wood
!                       growth, queued at the TRIFFID step in
!                       ximpair_repair_pend (triffid_nsc_ring_renewal).
!   ximpair_tau_leaf  - with per-segment memory, the leaf segment recovers
!                       with leaf turnover (e-folding time in calendar
!                       years), or by new leaves with l_ximpair_rec_leaf_turn.
!   ximpair_reset_mmdd - an annual reset (off by default).
! The time-based recovery of the stem (l_ximpair_rec_growth, new xylem
! from carbon gain on a growth clock; ximpair_rec_years; ximpair_tau_stem;
! refilling with ximpair_tau_rec) and, earlier, a leaf-area term were
! retired in favour of the store repair, and the store repair (26/10/26)
! in favour of ring renewal. They are in the git history (e.g. branch
! ximpair_prune; 6133017 for the store repair).
!
! Per-segment memory (l_ximpair_seg_memory, with l_som_plant_segments):
! kcap is the stem segment's cap and ximpair_kcap_leaf the leaf segment's
! (as a fraction of the intact leaf segment). Each is damaged by its own
! leaf-side water potential on its own curve (conductance_b/c_seg): the stem
! by psi_stem, the leaf by psi_leaf (ximpair_psi_driver is not used). The
! leaf's loss falls with leaf turnover, g years after the damage,
!   PLC_leaf' = PLC_leaf EXP(-g / ximpair_tau_leaf),
! so a drought's damage to the distal, fast-renewed tissue fades within a
! few years while the stem keeps its damage until the rings renew it.
! The reset restores both.
!
! MDK, 26/10/26
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_memory ( n_land_pts                         &
,                                           n_open_pts                         &
,                                           open_index                         &
,                                           pft                                &
,                                           psi_leaf                           &
,                                           psi_root                           &
,                                           psi_stem                           &
,                                           kcap                               &
                                          )

USE pftparm, ONLY: kmax_pft, kcrit, conductance_b_pft, conductance_c_pft,      &
                   ximpair_psi_driver, ximpair_reset_mmdd,                     &
                   ximpair_tau_leaf, seg_rel_k
USE jules_vegetation_mod, ONLY: ximpair_driver_leaf, ximpair_driver_mean,      &
                                ximpair_driver_root, ximpair_driver_stem,      &
                                l_ximpair_seg_memory, l_ximpair_rec_nsc,       &
                                l_ximpair_rec_leaf_turn, l_trif_hydmort,       &
                                l_som_seg_whole
USE jules_vegetation_mod, ONLY: l_nsc_ring_additive   ! MDK, 29/10/26
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                   &
                                          ximpair_psi_stem,                    &
                                          ximpair_seg_state_alloc
USE triffid_hydmort_mod, ONLY: hydmort_accumulate
                            ! MDK, 21/10/26
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
, psi_stem(n_land_pts)
                            ! Stem water potential (Pa; leaf-side end of the
                            ! stem segment, = psi_root without segments)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  kcap(n_land_pts)
                            ! Conductance cap, kmax_pft * (1 - PLC_mem)
                            ! (leaf basis).

! Local variables
INTEGER :: j, l, errcode
REAL(KIND=real_jlslsm) :: rep_l
                            ! Queued ring share at one point.
                            !
                            ! MDK, 26/10/26

REAL(KIND=real_jlslsm) ::                                                      &
  psi_x(n_land_pts), k_x(n_land_pts), kmax_pts(n_land_pts),                    &
  kcrit_pts(n_land_pts), b_pts(n_land_pts), c_pts(n_land_pts)
                            ! Damage-driving water potential (Pa), the
                            ! intact conductance there, and the PFT's
                            ! intact curve (kmax, kcrit, b in Pa, c) on
                            ! each point.
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: dt_years
                            ! This timestep in years, for the leaf
                            ! segment's recovery.
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: frac_crit
                            ! kcrit / kmax (floor of the leaf segment's
                            ! cap).
                            !
                            ! MDK, 26/10/26

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

! Area gained by TRIFFID since the last update (ximpair_dilute_area): the
! stem's loss of conductance diluted by the new, intact area.
!
! MDK, 25/10/26
IF ( ALLOCATED(ximpair_area_w) ) THEN
  WHERE ( ximpair_area_w(:,pft) < 1.0 )                                        &
    kcap(:) = kmax_pts(:) - (kmax_pts(:) - kcap(:)) * ximpair_area_w(:,pft)
  ximpair_area_w(:,pft) = 1.0
END IF

!-----------------------------------------------------------------------------
! Recovery: the stem by ring renewal (l_ximpair_rec_nsc), the leaf segment
! by leaf turnover, and the annual reset. (The time-based recovery of the
! stem, l_ximpair_rec_growth, ximpair_rec_years, ximpair_tau_stem and
! ximpair_tau_rec, was retired. MDK, 23/10/26)
!
! MDK, 26/10/26
!-----------------------------------------------------------------------------

IF (l_ximpair_rec_nsc) THEN
  CALL ximpair_memory_alloc( n_land_pts )
  ! Ring renewal: repair_pend is the share of the stem's conductance
  ! renewed by the new rings since the last update, new and old xylem
  ! replaced alike (random), so the whole stem's loss is diluted by it,
  ! L -> L (1 - rho). Points with nothing queued are left untouched.
  !
  ! MDK, 26/10/26
  DO l = 1, n_land_pts
    rep_l = ximpair_repair_pend(l,pft)
    IF (rep_l <= 0.0) CYCLE
    IF (l_nsc_ring_additive) THEN
      ! The new ring replaces damaged xylem first (MEDFATE): the loss falls
      ! by the ring's share. MDK, 29/10/26
      kcap(l) = kmax_pts(l) * (1.0 - MAX(1.0 - kcap(l) / kmax_pts(l)          &
                                         - rep_l, 0.0))
    ELSE
      kcap(l) = kmax_pts(l) * (1.0 - MAX(1.0 - kcap(l) / kmax_pts(l), 0.0)    &
                                     * (1.0 - MIN(rep_l, 1.0)))
    END IF
  END DO
  ximpair_repair_pend(:,pft) = 0.0
END IF

! This timestep in years, for the leaf segment's recovery.
dt_years = REAL(timestep_len) / sec_per_year

! With the segments on the PFT curve (l_som_seg_whole) the stem psi is
! kept for the psi_stem_pft output only (the leaf state is not used).
!
! MDK, 26/10/26
IF (l_som_seg_whole) THEN
  IF (.NOT. ALLOCATED(ximpair_psi_stem)) THEN
    ALLOCATE(ximpair_psi_stem(n_land_pts, npft))
    ximpair_psi_stem(:,:) = 0.0
  END IF
  ximpair_psi_stem(:,pft) = psi_stem(:)
END IF

IF (l_ximpair_seg_memory) THEN
  CALL ximpair_seg_state_alloc( n_land_pts )
  ximpair_psi_stem(:,pft) = psi_stem(:)
  ! With leaf renewal (l_ximpair_rec_leaf_turn) the leaf segment recovers
  ! by new leaves instead (triffid_nsc_turgor, phenol).
  !
  ! MDK, 18/10/26
  IF (ximpair_tau_leaf(pft) > 0.0 .AND. .NOT. l_ximpair_rec_leaf_turn) THEN
    ximpair_kcap_leaf(:,pft) = 1.0 - (1.0 - ximpair_kcap_leaf(:,pft))          &
                               * EXP(-dt_years / ximpair_tau_leaf(pft))
  END IF
END IF

IF (ximpair_reset_mmdd(pft) > 0) THEN
  IF (current_time%month * 100 + current_time%day == ximpair_reset_mmdd(pft)   &
      .AND. current_time%time < timestep_len) THEN
    kcap(:) = kmax_pts(:)
    IF (l_ximpair_seg_memory) ximpair_kcap_leaf(:,pft) = 1.0
  END IF
END IF

!-----------------------------------------------------------------------------
! Damage: the cap can not exceed the intact conductance at the damage
! driver. Only points with open stomata (i.e. under tension) are damaged.
!-----------------------------------------------------------------------------
IF (l_ximpair_seg_memory) THEN
  ! Per segment: the stem by psi_stem on the stem curve, the leaf by
  ! psi_leaf on the leaf curve (the segments' intact Weibull curves).
  DO j = 1, n_open_pts
    l = open_index(j)
    kcap(l) = MIN(kcap(l), kmax_pts(l) * seg_rel_k(pft, 2, psi_stem(l)))
    ximpair_kcap_leaf(l,pft) = MIN(ximpair_kcap_leaf(l,pft),                   &
                                   seg_rel_k(pft, 3, psi_leaf(l)))
  END DO
  frac_crit = kcrit(pft) / MAX(kmax_pft(pft), TINY(1.0_real_jlslsm))
  ximpair_kcap_leaf(:,pft) = MAX(MIN(ximpair_kcap_leaf(:,pft), 1.0), frac_crit)
ELSE
  CALL leaf_conductance_jls( pft, n_land_pts, psi_x, kmax_pts, kcrit_pts,      &
                             b_pts, c_pts, k_x )
  DO j = 1, n_open_pts
    l = open_index(j)
    kcap(l) = MIN(kcap(l), k_x(l))
  END DO
END IF

kcap(:) = MAX(MIN(kcap(:), kmax_pts(:)), kcrit_pts(:))

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
ximpair_lock(:,pft) = 1.0 - kcap(:) / kmax_pts(:)
! Worst damage of the TRIFFID period, for hydraulic mortality.
!
! MDK, 21/10/26
IF (l_trif_hydmort) CALL hydmort_accumulate( n_land_pts, pft,                  &
                                             ximpair_lock(:,pft) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_memory

END MODULE xylem_impairment_memory_mod
