! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_kmax_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_KMAX_MOD'

PUBLIC :: leaf_conductance_impaired_kmax_jls,                                  &
          xylem_conductance_impaired_kmax_stom_opt_jls,                        &
          leaf_psi_impaired_kmax,                                              &
          update_xylem_impairment_kmax,                                        &
          update_xylem_impairment_kmax_refit

CONTAINS

! *********************************************************************
! Contains routines used to apply the xylem impairment model based on
! reducing the maximum xylem conductance (kmax).
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calculate the leaf conductance from the leaf water
! potential.
! ---------------------------------------------------------------------

SUBROUTINE leaf_conductance_impaired_kmax_jls( pft,                               &
                                               land_pnts,                         &
                                               water_potential,                   &
                                               kmax,                              &
                                               kcrit,                             &
                                               conductance_b,                     &
                                               conductance_c,                     &
                                             ! INTENT OUT
                                               leaf_conductance                   &
)

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls

USE ereport_mod, ONLY: ereport
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
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_conductance(land_pnts)
                            ! Leaf conductance for each land point (m/s)

! Local variables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_KMAX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The impairment is handled by the kmax, conductance_b and conductance_c
! parameters so we just need to call the default leaf conductance code.
! NOTE: the per-point arrays are passed whole. This used to loop over the
!       points passing one point's scalar kmax(i)/b(i)/c(i) per call, but
!       each call overwrote the whole output array, so every point ended up
!       using the last point's curve.
CALL leaf_conductance_jls( pft,                                                &
                           land_pnts,                                          &
                           water_potential,                                    &
                           kmax,                                               &
                           kcrit,                                              &
                           conductance_b,                                      &
                           conductance_c,                                      &
                         ! INTENT OUT
                           leaf_conductance                                    &
                         )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_kmax_jls

! ---------------------------------------------------------------------
! Function to manage xylem conductance calculations for different
! impairment models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It asumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_impaired_kmax_stom_opt_jls( pft,                               &
                                            n_water_potentials,                &
                                            open_pnts,                         &
                                            open_index,                        &
                                            veg_index,                         &
                                            land_pnts,                         &
                                            water_potential,                   &
                                            kmax,                              &
                                            kcrit,                             &
                                            conductance_b,                     &
                                            conductance_c,                     &
                                          ! INTENT OUT
                                            xylem_conductance                  &
  )

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_water_potentials                                                           &
                            ! Number of water potentials per open point
, open_pnts                                                                    &
                            ! Number of open stomata
, open_index(open_pnts)                                                        &
                            ! Index of open stomata into veg_index
, land_pnts                                                                    &
                            ! Number of land points
, veg_index(land_pnts)
                            ! Index of vegetated points on the land grid

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(n_water_potentials, open_pnts)                               &
                            ! Water potentials for each open point (Pa)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: i, l

REAL(KIND=real_jlslsm) ::                                                      &
  kmax_open(open_pnts), kcrit_open(open_pnts), b_open(open_pnts),             &
  c_open(open_pnts)
                            ! Per-point curve parameters gathered onto the
                            ! open-point index.

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_KMAX_STOM_OPT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! The impairment is handled by the kmax, conductance_b and conductance_c
! parameters so we just need to call the default xylem conductance code,
! with the per-point parameters gathered onto the open-point index.
DO i = 1, open_pnts
  l = veg_index(open_index(i))
  kmax_open(i)  = kmax(l)
  kcrit_open(i) = kcrit(l)
  b_open(i)     = conductance_b(l)
  c_open(i)     = conductance_c(l)
END DO

CALL xylem_conductance_jls( pft,                                               &
                            n_water_potentials,                                &
                            open_pnts,                                         &
                            water_potential,                                   &
                            kmax_open,                                         &
                            kcrit_open,                                        &
                            b_open,                                            &
                            c_open,                                            &
                          ! INTENT OUT
                            xylem_conductance                                  &
                            )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_impaired_kmax_stom_opt_jls

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf water potential from the
! transpiration rate.
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_impaired_kmax( pft,                                        &
                                   n_e_leaf,                                   &
                                   land_pnts,                                  &
                                   open_pnts,                                  &
                                   veg_index,                                  &
                                   open_index,                                 &
                                   e_leaf,                                     &
                                   root_zone_psi,                              &
                                   kmax,                                       &
                                   kcrit,                                      &
                                   conductance_b,                              &
                                   conductance_c,                              &
                                ! INTENT OUT
                                   leaf_psi,                                   &
                                   leaf_k                                      &
  )

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_e_leaf                                                                     &
                            ! Number of transpirations per land point
, land_pnts                                                                    &
                            ! Number of land points
, open_pnts                                                                    &
                            ! Number of open stomata
, veg_index(land_pnts)                                                         &
                            ! Index of vegetation points on the land grid
, open_index(land_pnts)
                            ! Index of open stomata for each land point

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point (kg/m2/s)
, root_zone_psi(land_pnts)                                                     &
                             ! Water potential in the root zone (Pa)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: i, l

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_KMAX_STOM_OPT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The impairment is handled by the kmax, conductance_b and conductance_c
! parameters so we just need to call the default leaf water potential code.
! NOTE: the per-point arrays are passed whole - see the matching note in
!       leaf_conductance_impaired_kmax_jls above.
CALL leaf_psi_jls( pft,                                                        &
                   n_e_leaf,                                                   &
                   land_pnts,                                                  &
                   open_pnts,                                                  &
                   veg_index,                                                  &
                   open_index,                                                 &
                   e_leaf,                                                     &
                   root_zone_psi,                                              &
                   kmax,                                                       &
                   kcrit,                                                      &
                   conductance_b,                                              &
                   conductance_c,                                              &
                 ! INTENT OUT
                   leaf_psi,                                                   &
                   leaf_k                                                      &
                 )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_impaired_kmax

! ---------------------------------------------------------------------
! Function to update the xylem conductance model parameters.
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_kmax ( n_land_pts                           &
,                                         n_open_pts                           &
,                                         open_index                           &
,                                         pft                                  &
,                                         leaf_conductance                     &
,                                         root_conductance                     &
,                                         impaired_k_max                       &
,                                         impaired_conductance_b               &
,                                         impaired_conductance_c               &
                                        )

! Impairment model parameters
USE pftparm, ONLY: ximpair_leaf_weight, ximpair_new_kmax_weight,               &
                   ximpair_threshold

! Conductance model parameters
USE pftparm, ONLY: pft_conductance_model
USE pftparm, ONLY: kmax_pft, conductance_b_pft, conductance_c_pft, kcrit
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  n_land_pts                                                                   &
                            ! Number of land points
, n_open_pts                                                                   &
                            ! Number of open plant functional types
, open_index(n_open_pts)                                                       &
                            ! Index of open plant functional types
, pft
                            ! Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  leaf_conductance(n_land_pts)                                                 &
                            ! Leaf conductance (m/s)
, root_conductance(n_land_pts)
                            ! Root conductance (m/s)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  impaired_k_max(n_land_pts)
                            ! Maximum xylem conductance (m/s)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  impaired_conductance_b(n_land_pts)                                           &
                            ! Conductance model parameter b (MPa)
, impaired_conductance_c(n_land_pts)
                            ! Conductance model parameter c (unitless)

Real(KIND=real_jlslsm) :: impaired_k_max_new(n_land_pts)

INTEGER :: j

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_KMAX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the target value of kmax to a weighted average of the root zone
! and leaf conductance.
impaired_k_max_new = ximpair_leaf_weight(pft) * leaf_conductance            &
                    + (1-ximpair_leaf_weight(pft)) * root_conductance

! Calculate a value for the new kmax based on the leaf conductance.
! The value 0.12 determins the conductance loss threshold at which
! damage occurs.
!impaired_k_max_new = (1/(1-0.12)) * leaf_conductance
impaired_k_max_new = (impaired_k_max_new /                                   &
                      (impaired_k_max - ximpair_threshold(pft)))                 &
                      * impaired_k_max

! Take the weighted average of the new and old kmax values.
impaired_k_max_new = ximpair_new_kmax_weight(pft) * impaired_k_max_new       &
                     + (1 - ximpair_new_kmax_weight(pft))                    &
                       * impaired_k_max

! Limit the new maximum conductance to less than or equal to the previous
! maximum conductance.
impaired_k_max_new = MIN(impaired_k_max, impaired_k_max_new)

! Impose minimum k_max as the critical conductance.
impaired_k_max_new = MAX(impaired_k_max_new, kcrit(pft))

! Update the value of kmax only when stomata are open
! Loop over land points with open stomata
DO j = 1, n_open_pts
  impaired_k_max(open_index(j)) = impaired_k_max_new(open_index(j))
END DO

! Calculation of b and c depend on conductance model
SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )

  ! updat the value of conductance_b
  !
  !   b' = b * (1 - log(k_max' / k_max))^(1/c)
  !
  impaired_conductance_b = conductance_b_pft(pft)                                    &
                           * (1 - LOG(impaired_k_max / kmax_pft(pft)))           &
                              **(1/conductance_c_pft(pft))

  ! updat the value of conductance_c
  !
  !   c' = c * (b / b')^(c)
  !
  impaired_conductance_c = conductance_c_pft(pft)                                    &
                           * (conductance_b_pft(pft) / impaired_conductance_b)       &
                              **conductance_c_pft(pft)

CASE ( SOX_conductance )

  ! updat the value of conductance_b
  !
  !   b' = b * ( 2*(kmax/kmax') - 1 )^(1/c)
  !
  impaired_conductance_b = conductance_b_pft(pft)                                    &
                           * (2*(kmax_pft(pft) / impaired_k_max) - 1)            &
                              **(1/conductance_c_pft(pft))

  ! updat the value of conductance_c
  !
  !   c' = c * (kmax/kmax') * (b'/b)^c
  !
  impaired_conductance_c = conductance_c_pft(pft)                                    &
                           * (kmax_pft(pft) / impaired_k_max)                    &
                           * (impaired_conductance_b / conductance_b_pft(pft))       &
                               **conductance_c_pft(pft)

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_kmax

! ---------------------------------------------------------------------
! kmax refit impairment model (pft_xylem_impairment_model = 4): the kmax
! model made consistent with Mackay et al. (2015, WRR) / TREES, which it is
! based on:
!   - Kcav (the impaired kmax) = kmax_pft * f_intact(psi_min), with psi_min
!     the running minimum of the damage-driving water potential
!     (ximpair_psi_driver), i.e. evaluated on the intact curve with no
!     threshold or weights (the kmax model instead used a weighted leaf/root
!     conductance on the already-impaired curve, which feeds back on itself,
!     and a threshold term that goes negative once kmax < f * kmax_pft).
!     This, and the optional recovery, is the same state update as the
!     embolism memory model (update_xylem_impairment_memory).
!   - conductance_b/c are refitted to the capped curve as TREES does, with
!     xylem_refit_weibull (the kmax model's c' went the wrong way).
! The impaired curve is then kcap * f(psi; b', c'), using the kmax model's
! conductance/leaf water potential routines.
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_kmax_refit ( n_land_pts                     &
,                                               n_open_pts                     &
,                                               open_index                     &
,                                               pft                            &
,                                               psi_leaf                       &
,                                               psi_root                       &
,                                               lai                            &
,                                               canht                          &
,                                               anetc                          &
,                                               impaired_k_max                 &
,                                               impaired_conductance_b         &
,                                               impaired_conductance_c         &
                                              )

USE pftparm, ONLY: kmax_pft, conductance_b_pft, conductance_c_pft
USE xylem_impairment_memory_mod, ONLY: update_xylem_impairment_memory,         &
                                       xylem_refit_weibull

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  n_land_pts, n_open_pts, open_index(n_open_pts), pft
                            ! See update_xylem_impairment_memory.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  psi_leaf(n_land_pts), psi_root(n_land_pts), lai(n_land_pts),               &
  canht(n_land_pts), anetc(n_land_pts)
                            ! See update_xylem_impairment_memory.

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  impaired_k_max(n_land_pts)                                                   &
                            ! Impaired maximum xylem conductance, Kcav
                            ! (leaf basis).
, impaired_conductance_b(n_land_pts)                                           &
                            ! Refitted conductance parameter b (Pa)
, impaired_conductance_c(n_land_pts)
                            ! Refitted conductance parameter c

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_KMAX_REFIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL update_xylem_impairment_memory( n_land_pts, n_open_pts, open_index, pft,  &
                                     psi_leaf, psi_root, lai, canht, anetc,    &
                                     impaired_k_max )

CALL xylem_refit_weibull( pft, n_land_pts, impaired_k_max, kmax_pft(pft),      &
                          conductance_b_pft(pft), conductance_c_pft(pft),      &
                          impaired_conductance_b, impaired_conductance_c )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_kmax_refit


END MODULE xylem_impairment_kmax_mod