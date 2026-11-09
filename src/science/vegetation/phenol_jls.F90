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
USE jules_vegetation_mod, ONLY: l_nitrogen,                                    &
                                l_ximpair_rec_leaf_turn,                       &
                                l_triffid_nsc, l_trait_phys,                   &
                                ! TEMP RED
                                l_red
USE triffid_nsc_mod, ONLY: nsc_pft, nsc_regrow_acc
USE triffid_nsc_mod, ONLY: triffid_nsc_alloc, nsc_lai_sap_cap
USE jules_vegetation_mod, ONLY: nsc_damage_response, nsc_damage_trugman,       &
                                l_shed_regrow_flush
USE leaf_shed_mod, ONLY: leaf_shed_alloc, leaf_shed_caps, phen_shed_debt,     &
                         phen_regrow_c, leaf_shed_heat_clear
                           ! MDK, 09/11/26
USE model_time_mod, ONLY: is_spinup
                           ! MDK, 01/11/26
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf
USE calc_c_comps_triffid_mod, ONLY: lai_bal_from_ht   ! MDK, 20/10/26
USE pftparm, ONLY: a_wl, a_ws, b_wl, eta_sl, g_leaf_0,                         &
                   sigl, nsc_flush, nsc_cmax, nsc_reserve_frac,                &
                   lma, nsc_flush_start_doy, nsc_flush_days
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE jules_surface_mod, ONLY: cmass
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
    ! All regrowth is paid from the store this call (inside the leaf-flush
    ! window).
REAL(KIND=real_jlslsm) :: doy
    ! Day of the year.

LOGICAL :: l_cold
    ! The cold (or TRIFFID drought) leaf drop acted this call.

REAL(KIND=real_jlslsm) :: phen_old, dpay
    ! Phenological state before this call, and the part of the regrowth
    ! paid from the store.

REAL(KIND=real_jlslsm) ::                                                      &
  lai_bal(land_pts),                                                           &
    ! Balanced growth LAI.
  phen(land_pts),                                                              &
    ! Phenological state.
  phen_max(land_pts),                                                          &
    ! Largest phenological state the leaf shedding caps allow.
  phen_max_free(land_pts),                                                     &
    ! The same from the caps whose regrowth is free (damage).
  lai_cap(land_pts)
    ! Canopy the damage cap allows (m2 m-2; < 0 = no cap).
    !
    ! MDK, 09/11/26

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
! Leaf shedding (leaf_shed_mod): the caps on the phenological state set by
! drought, heat and xylem damage. Leaves above the cap are shed at the rate
! of phenological leaf drop and regrow as the cap rises. These are caps on
! the canopy rather than extra turnover rates: leaf drop below is a
! threshold on g_leaf, so amplifying the turnover would either do nothing
! or drop the whole canopy.
!
! MDK, 09/11/26
!-----------------------------------------------------------------------------
CALL leaf_shed_alloc( land_pts )
CALL leaf_shed_caps( land_pts, veg_pts, veg_index, n, lai_bal, phen_max,      &
                     phen_max_free, lai_cap )
IF (nsc_damage_response == nsc_damage_trugman) THEN
  CALL triffid_nsc_alloc( land_pts )
  DO j = 1,veg_pts
    l = veg_index(j)
    nsc_lai_sap_cap(l,n) = lai_cap(l)
  END DO
END IF

!-----------------------------------------------------------------------------
! Update the phenological state and output the leaf turnover rate in
! terms of the balanced growth LAI
!-----------------------------------------------------------------------------
DO j = 1,veg_pts
  l = veg_index(j)
  phen_old           = phen(l)
  l_cold             = .FALSE.
  phen_regrow_c(l,n) = 0.0

  ! Budburst gate: no cold leaf drop in the first phen_bb_min_days after
  ! budburst.
  !
  ! MDK, 20/10/26
  IF (g_leaf(l) >  2.0 * g_leaf_0(n) .AND. .NOT. phen_bb_in_leaf(l, n)) THEN
    dphen          = -dtime_phen * g_grow(n)
    dphen          = MAX(dphen,(0.01 - phen(l)))
    g_leaf_phen(l) = -dphen / dtime_phen
    l_cold         = .TRUE.

  ELSE IF (phen(l) > phen_max(l)) THEN
    ! Shed the leaf area above the leaf shedding caps; the shed leaves add
    ! to the leaf turnover.
    dphen          = -dtime_phen * g_grow(n)
    dphen          = MAX(dphen,(phen_max(l) - phen(l)))
    g_leaf_phen(l) = MAX(phen(l) * g_leaf(l), -dphen / dtime_phen)
    ! In the nitrogen scheme g_leaf_phen tracks the net change in phen, as
    ! in the other two branches.
    !
    ! MDK, 24/10/26
    IF (l_nitrogen) g_leaf_phen(l) = -dphen / dtime_phen

    ! The part shed below the free (damage) cap was shed by drought or
    ! heat: its regrowth is owed to the store.
    !
    ! MDK, 09/11/26
    phen_shed_debt(l,n) = phen_shed_debt(l,n)                                  &
      + MAX(MIN(phen_old, phen_max_free(l)) - (phen_old + dphen), 0.0)

    ! The leaves shed by the heat cap take their damage with them.
    !
    ! MDK, 09/11/26
    CALL leaf_shed_heat_clear( l, n, phen_max_free(l), phen_old,               &
                               phen_old + dphen )

  ELSE
    dphen          = dtime_phen * g_grow(n) * (phen_max(l) - phen(l))
    dphen          = MIN(dphen,(phen_max(l) - phen(l)))
    g_leaf_phen(l) = phen(l) * g_leaf(l)

    ! Budburst gate (phen_bb_doy or phen_bb_gdd): no regrowth while dormant
    ! or before this year's budburst.
    !
    ! MDK, 20/10/26
    IF (phen_bb_blocked(l, n)) dphen = 0.0

    ! Regrowth paid from the accessible store (queued for the next store
    ! step) and limited by it (l_triffid_nsc):
    ! - the leaves shed by drought or heat, up to the debt phen_shed_debt,
    !   in any season (with l_shed_regrow_flush only in the flush window);
    ! - with the leaf flush (nsc_flush_start_doy set), all regrowth inside
    !   the flush window, e.g. a deciduous canopy rebuilt each spring.
    ! Other regrowth stays free, as in TRIFFID: paying for the brief
    ! regrowth of mild winter spells, or for the spring canopy after the
    ! cold leaf drop, would drain the store and collapse a deciduous canopy.
    !
    ! MDK, 09/11/26
    l_pay = .FALSE.
    IF (is_set(nsc_flush_start_doy(n))) THEN
      doy = REAL(day_of_year(current_time%year, current_time%month,            &
                             current_time%day, l_360, l_leap))
      ! The flush window may cross the year end.
      l_pay = nsc_in_flush_window(n, NINT(doy),                                &
                          days_in_year(current_time%year, l_360, l_leap))
    END IF
    dpay = 0.0
    IF (dphen > 0.0) THEN
      IF (l_pay) THEN
        dpay = dphen
      ELSE
        dpay = MIN(dphen, phen_shed_debt(l,n))
        ! l_shed_regrow_flush: outside the flush window the shed leaves
        ! wait for the next flush.
        IF (l_shed_regrow_flush .AND. is_set(nsc_flush_start_doy(n))) THEN
          dphen = dphen - dpay
          dpay  = 0.0
        END IF
      END IF
    END IF
    IF (l_triffid_nsc .AND. dpay > 0.0 .AND. ALLOCATED(nsc_pft)) THEN
      IF (nsc_pft(l,n) >= 0.0) THEN
        leafc_b = sigl(n) * lai_bal(l)
        IF (l_trait_phys) leafc_b = cmass * lma(n) * lai_bal(l)
        c_max_l = nsc_flush(n) * 2.0 * leafc_b
        IF (nsc_cmax(n) > 0.0) c_max_l = nsc_cmax(n)
        c_free  = MAX(nsc_pft(l,n) - nsc_reserve_frac(n) * c_max_l             &
                      - nsc_regrow_acc(l,n), 0.0)
        IF (leafc_b > 0.0) THEN
          ! The free part of the regrowth goes ahead; the paid part only as
          ! far as the store allows.
          dphen = dphen - dpay
          dpay  = MIN(dpay, c_free / leafc_b)
          dphen = dphen + dpay
          nsc_regrow_acc(l,n) = nsc_regrow_acc(l,n) + dpay * leafc_b
          phen_regrow_c(l,n)  = dpay * leafc_b / dtime_phen
        END IF
      END IF
    END IF
    ! Only regrowth of the shed leaves clears the debt: dpay, the part of
    ! the regrowth counted against the debt (after the store's limit when
    ! the store is on; without the store it is regrown free). Free
    ! regrowth of other leaves (after a cold drop or under the damage cap)
    ! leaves the debt as it is.
    !
    ! MDK, 09/11/26
    IF (dpay > 0.0) phen_shed_debt(l,n) = MAX(phen_shed_debt(l,n) - dpay, 0.0)

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
  ! budburst goes dormant until the next budburst. Only the cold (or
  ! TRIFFID drought) leaf drop sends it dormant, not shedding by a cap.
  !
  ! MDK, 09/11/26
  IF (phen_bb_on(n) .AND. ALLOCATED(phen_dormant)) THEN
    IF (l_cold .AND. phen(l) <= 0.01 + 1.0e-6 .AND.                            &
        phen_bb_day(l,n) >= 0.0) phen_dormant(l,n) = 1.0
  END IF

  ! The debt cannot exceed the leaves missing from the full canopy.
  !
  ! MDK, 09/11/26
  phen_shed_debt(l,n) = MIN(phen_shed_debt(l,n), MAX(1.0 - phen(l), 0.0))

  ! Leaf renewal (l_ximpair_rec_leaf_turn): leaves regrown after shedding
  ! arrive intact, so the leaf segment's loss falls by their share. Not in
  ! spin-up, where the damage is held.
  !
  ! MDK, 04/11/26
  IF (l_ximpair_rec_leaf_turn .AND. .NOT. is_spinup .AND. dphen > 0.0 .AND.   &
      phen(l) > 0.0) THEN
    IF (ALLOCATED(ximpair_kcap_leaf))                                          &
      ximpair_kcap_leaf(l,n) = 1.0 - (1.0 - ximpair_kcap_leaf(l,n))            &
                               * (1.0 - MIN(dphen / phen(l), 1.0))
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE phenol

END MODULE phenol_mod
