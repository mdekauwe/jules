! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE phenol_mod


IMPLICIT NONE


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHENOL_MOD'

CONTAINS

!----------------------------------------------------------------------------
! Subroutine PHENOL ---------------------------------------------------------
!
!
! Purpose :  Parametrizes leaf phenological changes and updates the
!            leaf area index and the leaf turnover rate.
!
! ---------------------------------------------------------------------------
SUBROUTINE phenol (land_pts, veg_pts, n, veg_index, dtime_phen, g_leaf, ht,    &
                   lai, g_leaf_phen,lai_bal_red)

USE um_types, ONLY: real_jlslsm
USE pftparm, ONLY: is_set, nsc_in_flush_window   ! MDK, 24/10/26
USE jules_vegetation_mod, ONLY: l_nitrogen, l_ximpair_leaf_loss,             &
                                l_ximpair_rec_leaf_turn, l_phen_psi,           &
                                l_triffid_nsc,                                 &
                                ! TEMP RED
                                l_red
USE triffid_nsc_mod, ONLY: phen_fbleaf, nsc_pft, nsc_regrow_acc
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf
USE calc_c_comps_triffid_mod, ONLY: lai_bal_from_ht   ! MDK, 20/10/26
USE pftparm, ONLY: a_wl, a_ws, b_wl, eta_sl, g_leaf_0, ximpair_leaf_sens,    &
                   sigl, nsc_flush, nsc_cmax, nsc_reserve_frac, phen_psi_tlp,   &
                   lma, nsc_flush_start_doy, nsc_flush_days
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE jules_vegetation_mod, ONLY: l_trait_phys
USE jules_surface_mod, ONLY: cmass
USE xylem_impairment_memory_mod, ONLY: ximpair_lock
USE trif, ONLY: g_grow
USE phen_budburst_mod, ONLY: phen_bb_on, phen_bb_blocked, phen_bb_day,         &
                             phen_dormant, phen_bb_in_leaf

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  veg_pts,                                                                     &
    ! Number of vegetated points.
  n
    ! Plant functional type.

INTEGER, INTENT(IN) ::                                                         &
  veg_index(land_pts)
    ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm) , INTENT(IN) ::                                         &
  dtime_phen
    ! Timestep (years).

REAL(KIND=real_jlslsm) , INTENT(IN) ::                                         &
  g_leaf(land_pts),                                                            &
    ! Rate of leaf turnover (/360days).
  ht(land_pts)
    ! Canopy height (m).
!-----------------------------------------------------------------------------
! Optional Arguments with INTENT(IN)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  lai_bal_red(land_pts)
    ! Balanced growth LAI from veg3.

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  lai(land_pts)
    ! Leaf area index.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  g_leaf_phen(land_pts)
    ! Rate of leaf turnover including leaf phenology (/360days).


!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  j,l       ! Loop counters

REAL(KIND=real_jlslsm) ::                                                      &
  dphen     ! Increment to phenological state.

REAL(KIND=real_jlslsm) ::                                                      &
  leafc_b, c_max_l, c_free
    ! Leaf carbon at the balanced LAI, store capacity and accessible store
    ! less queued regrowth (kg C m-2), for regrowth paid from the store.

LOGICAL :: l_pay
    ! Regrowth paid from the store this call (drought shedding, or inside
    ! the leaf-flush window).
REAL(KIND=real_jlslsm) :: doy
    ! Day of the year.

REAL(KIND=real_jlslsm) ::                                                      &
  lai_bal(land_pts),                                                           &
    ! Balanced growth LAI.
  phen(land_pts),                                                              &
    ! Phenological state.
  phen_max(land_pts)
    ! Largest phenological state the (damaged) xylem can supply.

! Dr Hook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHENOL'

!-----------------------------------------------------------------------------
!end of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Diagnose the phenological state
!-----------------------------------------------------------------------------

DO j = 1,veg_pts
  l = veg_index(j)
  ! If we're running l_red we use lai_bal as an input, otherwise estimate the
  ! classic TRIFFID way.
  IF (l_red .AND. PRESENT(lai_bal_red)) THEN
    lai_bal(l) = lai_bal_red(l)
  ELSE
    lai_bal(l) = lai_bal_from_ht(n, ht(l))   ! MDK, 20/10/26
  END IF
  phen(l)    = lai(l) / lai_bal(l)
END DO

!-----------------------------------------------------------------------------
! Lasting xylem damage (memory impairment models, l_ximpair_leaf_loss): the
! canopy cannot be fuller than the damaged xylem can supply,
!   phen <= 1 - ximpair_leaf_sens * (1 - k_cap/kmax),
! so it sheds leaves (at the rate of phenological leaf drop) when the
! locked-in loss rises and regrows as the xylem recovers. This is a cap on
! the canopy rather than an extra turnover rate: leaf drop below is a
! threshold on g_leaf, so amplifying the turnover would either do nothing
! or drop the whole canopy.
!-----------------------------------------------------------------------------
phen_max(:) = 1.0
IF (l_ximpair_leaf_loss .AND. ALLOCATED(ximpair_lock)) THEN
  DO j = 1,veg_pts
    l = veg_index(j)
    phen_max(l) = MAX(0.01, 1.0 - ximpair_leaf_sens(n)                         &
                                  * MAX(ximpair_lock(l,n), 0.0))
  END DO
END IF

! Drought leaf shedding (l_phen_psi): the canopy is also capped at the leaf
! biomass factor of the water status (Xu et al. 2016; triffid_nsc_mod).
! MDK, 18/10/26
IF (l_phen_psi .AND. is_set(phen_psi_tlp(n)) .AND.                             &
    ALLOCATED(phen_fbleaf)) THEN
  DO j = 1,veg_pts
    l = veg_index(j)
    phen_max(l) = MIN(phen_max(l), MAX(phen_fbleaf(l,n), 0.01))
  END DO
END IF

!-----------------------------------------------------------------------------
! Update the phenological state and output the leaf turnover rate in
! terms of the balanced growth LAI
!-----------------------------------------------------------------------------
DO j = 1,veg_pts
  l = veg_index(j)

  ! Budburst gate: no cold leaf drop in the first phen_bb_min_days after
  ! budburst.
  ! MDK, 20/10/26
  IF (g_leaf(l) >  2.0 * g_leaf_0(n) .AND. .NOT. phen_bb_in_leaf(l, n)) THEN
    dphen          = -dtime_phen * g_grow(n)
    dphen          = MAX(dphen,(0.01 - phen(l)))
    g_leaf_phen(l) = -dphen / dtime_phen

  ELSE IF (phen(l) > phen_max(l)) THEN
    ! Shed the leaf area the damaged xylem cannot supply; the shed leaves
    ! add to the leaf turnover.
    dphen          = -dtime_phen * g_grow(n)
    dphen          = MAX(dphen,(phen_max(l) - phen(l)))
    g_leaf_phen(l) = MAX(phen(l) * g_leaf(l), -dphen / dtime_phen)

  ELSE
    dphen          = dtime_phen * g_grow(n) * (phen_max(l) - phen(l))
    dphen          = MIN(dphen,(phen_max(l) - phen(l)))
    g_leaf_phen(l) = phen(l) * g_leaf(l)

    ! Budburst gate (phen_bb_doy or phen_bb_gdd): no regrowth while dormant
    ! or before this year's budburst.
    ! MDK, 20/10/26
    IF (phen_bb_blocked(l, n)) dphen = 0.0

    ! Drought leaf shedding with the store (l_phen_psi, l_triffid_nsc):
    ! regrowth is paid from the accessible store (queued for the next store
    ! step) and limited by it.
    ! MDK, 18/10/26
    ! With the leaf flush (nsc_flush_start_doy set) regrowth inside the
    ! flush window is paid from the store too, e.g. a deciduous canopy
    ! rebuilt each spring. Outside the window it stays free, as in TRIFFID:
    ! paying for the brief regrowth of mild winter spells would drain the
    ! store and collapse a deciduous canopy.
    ! MDK, 20/10/26
    l_pay = l_phen_psi .AND. is_set(phen_psi_tlp(n))
    IF (.NOT. l_pay .AND. is_set(nsc_flush_start_doy(n))) THEN
      doy = REAL(day_of_year(current_time%year, current_time%month,          &
                             current_time%day, l_360, l_leap))
      ! The flush window may cross the year end.
      l_pay = nsc_in_flush_window(n, NINT(doy),                                &
                          days_in_year(current_time%year, l_360, l_leap))
    END IF
    IF (l_triffid_nsc .AND. dphen > 0.0 .AND. ALLOCATED(nsc_pft) .AND. l_pay) THEN
     IF (nsc_pft(l,n) >= 0.0) THEN
      leafc_b = sigl(n) * lai_bal(l)
      IF (l_trait_phys) leafc_b = cmass * lma(n) * lai_bal(l)
      c_max_l = nsc_flush(n) * 2.0 * leafc_b
      IF (nsc_cmax(n) > 0.0) c_max_l = nsc_cmax(n)
      c_free  = MAX(nsc_pft(l,n) - nsc_reserve_frac(n) * c_max_l               &
                    - nsc_regrow_acc(l,n), 0.0)
      IF (leafc_b > 0.0) THEN
        dphen = MIN(dphen, c_free / leafc_b)
        nsc_regrow_acc(l,n) = nsc_regrow_acc(l,n) + dphen * leafc_b
      END IF
     END IF
    END IF

    ! In the Nitrogen scheme use g_leaf_phen to keep track of net change in
    ! phen.
    IF (l_nitrogen) g_leaf_phen(l) = - dphen / dtime_phen

  END IF

  !---------------------------------------------------------------------------
  ! Update the leaf area index
  !---------------------------------------------------------------------------
  phen(l) = phen(l) + dphen
  lai(l)  = phen(l) * lai_bal(l)

  ! Budburst gate: the canopy fallen to its floor after this year's
  ! budburst goes dormant until the next budburst.
  ! MDK, 20/10/26
  IF (phen_bb_on(n) .AND. ALLOCATED(phen_dormant)) THEN
    IF (dphen < 0.0 .AND. phen(l) <= 0.01 + 1.0e-6 .AND.                      &
        phen_bb_day(l,n) >= 0.0) phen_dormant(l,n) = 1.0
  END IF

  ! Leaf renewal (l_ximpair_rec_leaf_turn): leaves regrown after shedding
  ! arrive intact, so the leaf segment's loss falls by their share.
  ! MDK, 18/10/26
  IF (l_ximpair_rec_leaf_turn .AND. dphen > 0.0 .AND. phen(l) > 0.0) THEN
    IF (ALLOCATED(ximpair_kcap_leaf))                                          &
      ximpair_kcap_leaf(l,n) = 1.0 - (1.0 - ximpair_kcap_leaf(l,n))           &
                               * (1.0 - MIN(dphen / phen(l), 1.0))
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE phenol

END MODULE phenol_mod
