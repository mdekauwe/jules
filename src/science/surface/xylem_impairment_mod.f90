! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_MOD'

PRIVATE

PUBLIC :: leaf_conductance_impaired_jls,                                       &
          psi_leaf_impaired_jls,                                               &
          update_xylem_impairment

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
!
! NOTE on the curve parameter arguments passed to the routines below:
!   kmax_ref      : unimpaired maximum xylem conductance for each land
!                   point, i.e. kmax_pft(pft) with whatever canopy
!                   scaling the caller applies (e.g. * LAI for the
!                   big-leaf stomatal optimisation). Used, together with
!                   conductance_b_pft/conductance_c_pft, whenever no
!                   impairment is applied (spinup, xylem_impairment_none).
!   kmax          : the (possibly impaired) per-point kmax, with the same
!                   scaling as kmax_ref.
!   conductance_b, conductance_c
!                 : the PFT's vulnerability curve on each point. Neither
!                   impairment model changes it.
!   kcrit         : critical xylem conductance, with the same scaling as
!                   kmax_ref.
!
! MDK, 24/10/26
! *********************************************************************

! ---------------------------------------------------------------------
! Function to manage leaf conductance calculations for different
! impairment models.
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_impaired_jls( pft,                                 &
                                          land_pnts,                           &
                                          water_potential,                     &
                                          kmax_ref,                            &
                                          kmax,                                &
                                          kcrit,                               &
                                          conductance_b,                       &
                                          conductance_c,                       &
                                          psi_leaf_extreme,                    &
                                        ! INTENT OUT
                                          leaf_conductance                     &
  )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls
USE xylem_impairment_whole_trunk_mod,                                          &
    ONLY: leaf_conductance_impaired_whole_trunk_jls
USE xylem_impairment_memory_mod, ONLY: leaf_conductance_impaired_memory_jls

USE um_types, ONLY: real_jlslsm

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
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance for each
                            ! land point (mol m-2 s-1 Pa-1)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_conductance(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: errcode

REAL(KIND=real_jlslsm) :: b_ref(land_pnts), c_ref(land_pnts)
                            ! Unimpaired PFT curve parameters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_IMPAIRED_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and use the default
! conductance model when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  b_ref(:) = conductance_b_pft(pft)
  c_ref(:) = conductance_c_pft(pft)

  CALL leaf_conductance_jls( pft,                                              &
                             land_pnts,                                        &
                             water_potential,                                  &
                             kmax_ref,                                         &
                             kcrit,                                            &
                             b_ref,                                            &
                             c_ref,                                            &
                           ! INTENT OUT
                             leaf_conductance                                  &
                             )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE ( pft_xylem_impairment_model(pft) )

CASE ( xylem_impairment_whole_trunk )
  CALL leaf_conductance_impaired_whole_trunk_jls( pft,                         &
                                                 land_pnts,                    &
                                                 water_potential,              &
                                                 kmax,                         &
                                                 kcrit,                        &
                                                 conductance_b,                &
                                                 conductance_c,                &
                                                 psi_leaf_extreme,             &
                                               ! INTENT OUT
                                                 leaf_conductance              &
                                                 )

CASE ( xylem_impairment_memory )
  CALL leaf_conductance_impaired_memory_jls( pft,                              &
                                             land_pnts,                        &
                                             water_potential,                  &
                                             kmax_ref,                         &
                                             kmax,                             &
                                             kcrit,                            &
                                             conductance_b,                    &
                                             conductance_c,                    &
                                           ! INTENT OUT
                                             leaf_conductance                  &
                                             )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
    'pft_xylem_impairment_model should be none (0), whole_trunk (2) ' //       &
    'or memory (3)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_jls

! ---------------------------------------------------------------------
! Subroutine to manage leaf water potential calculations for different
! conductance models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_impaired_jls( pft,                                         &
                                  n_e_leaf,                                    &
                                  land_pnts,                                   &
                                  open_pnts,                                   &
                                  veg_index,                                   &
                                  open_index,                                  &
                                  e_leaf,                                      &
                                  root_zone_psi,                               &
                                  kmax_ref,                                    &
                                  kmax,                                        &
                                  kcrit,                                       &
                                  conductance_b,                               &
                                  conductance_c,                               &
                                  psi_leaf_extreme,                            &
                                  psi_root_extreme,                            &
                               ! INTENT OUT
                                  leaf_psi,                                    &
                                  leaf_k                                       &
  )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls, leaf_conductance_jls
USE xylem_impairment_whole_trunk_mod, ONLY: psi_leaf_impaired_whole_trunk
USE xylem_impairment_memory_mod, ONLY: psi_leaf_impaired_memory
USE jules_vegetation_mod, ONLY: l_som_rhizo_series
USE plant_hydraulic_path_mod, ONLY: ksr_path, som_psi_in_min
                            ! MDK, 21/10/26

USE um_types, ONLY: real_jlslsm

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
                            ! Transpiration rates for each open point
, root_zone_psi(land_pnts)                                                     &
                             ! Water potential in the root zone (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance for each
                            ! land point (mol m-2 s-1 Pa-1)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)                                                  &
                            ! Minimum leaf water potential for each land point
                            ! (Pa)
, psi_root_extreme(land_pnts)
                            ! Minimum root water potential for each land point
                            ! (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: errcode, i, j, l

REAL(KIND=real_jlslsm) :: b_ref(land_pnts), c_ref(land_pnts)
                            ! Unimpaired PFT curve parameters.
REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(land_pnts), k_inlet(land_pnts)                                        &
                            ! l_som_rhizo_series: root water potential (Pa)
                            ! and impaired conductance there, one sample.
, e1(1, open_pnts), psi1(1, open_pnts), k1(1, open_pnts)                       &
                            ! One sample per open point.
, k_s                                                                          &
                            ! Soil-to-root conductance of the path.
, psi_k(land_pnts)
                            ! Whole trunk: root-end potential the root-end
                            ! conductance is taken at (Pa).
                            !
                            ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_IMPAIRED_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and just call the base
! conductance models when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  b_ref(:) = conductance_b_pft(pft)
  c_ref(:) = conductance_c_pft(pft)

  CALL psi_leaf_jls( pft,                                                      &
                     n_e_leaf,                                                 &
                     land_pnts,                                                &
                     open_pnts,                                                &
                     veg_index,                                                &
                     open_index,                                               &
                     e_leaf,                                                   &
                     root_zone_psi,                                            &
                     kmax_ref,                                                 &
                     kcrit,                                                    &
                     b_ref,                                                    &
                     c_ref,                                                    &
                   ! INTENT OUT
                     leaf_psi,                                                 &
                     leaf_k                                                    &
            )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE ( pft_xylem_impairment_model(pft) )

  CASE (xylem_impairment_whole_trunk, xylem_impairment_memory)
    IF ( .NOT. l_som_rhizo_series ) THEN
      CALL impaired_path( n_e_leaf, e_leaf, root_zone_psi, leaf_psi, leaf_k )
    ELSE
      !-----------------------------------------------------------------------
      ! Soil-to-root conductance in series ahead of the impaired plant, as
      ! psi_leaf_jls: the root water potential of each sample is psi_in =
      ! psi_src - E/K_s, the impaired plant path runs from there, and the
      ! whole-path marginal conductance is k_p / (1 + k(psi_in) / K_s),
      ! k(psi_in) this model's impaired conductance at the root. The soil is not
      ! damaged.
      !-----------------------------------------------------------------------
      DO i = 1, n_e_leaf
        psi_in(:) = root_zone_psi(:)
        DO j = 1, open_pnts
          l = veg_index(open_index(j))
          k_s = ksr_path(l)
          IF ( ksr_path(l) < 0.0 ) THEN
            ! No soil link for this PFT (fsmc_mod /= 2).
          ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
            psi_in(l) = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s,              &
                            som_psi_in_min)
          ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
            psi_in(l) = som_psi_in_min
          END IF
          e1(1,j) = e_leaf(i,j)
        END DO
        CALL impaired_path( 1, e1, psi_in, psi1, k1 )
        IF ( pft_xylem_impairment_model(pft) ==                                &
             xylem_impairment_whole_trunk ) THEN
          ! Whole trunk: the root end of the path is where water enters, so its
          ! conductance is the unimpaired curve at the root's historic
          ! minimum (as psi_leaf_impaired_whole_trunk uses there), not the
          ! leaf's, which leaf_conductance_impaired_jls would use.
          !
          ! MDK, 21/10/26
          psi_k(:) = MIN(psi_in(:), psi_root_extreme(:))
          CALL leaf_conductance_jls( pft, land_pnts, psi_k, kmax, kcrit,       &
                                     conductance_b, conductance_c, k_inlet )
        ELSE
          CALL leaf_conductance_impaired_jls( pft, land_pnts, psi_in,          &
                                              kmax_ref, kmax, kcrit,           &
                                              conductance_b, conductance_c,    &
                                              psi_leaf_extreme, k_inlet )
        END IF
        DO j = 1, open_pnts
          l = veg_index(open_index(j))
          k_s = ksr_path(l)
                            ! MDK, 21/10/26
          leaf_psi(i,j) = psi1(1,j)
          IF ( ksr_path(l) < 0.0 ) THEN
            leaf_k(i,j) = k1(1,j)
          ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
            leaf_k(i,j) = k1(1,j) / (1.0 + k_inlet(l) / k_s)
          ELSE
            leaf_k(i,j) = 0.0
          END IF
        END DO
      END DO
    END IF
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), whole_trunk (2) ' //      &
     'or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! The impaired plant path (root to leaf), whole_trunk or memory model.
SUBROUTINE impaired_path( n_e, e_in, psi_root, psi_out, k_out )

INTEGER, INTENT(IN) :: n_e
REAL(KIND=real_jlslsm), INTENT(IN)  :: e_in(n_e, open_pnts),                   &
                                       psi_root(land_pnts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi_out(n_e, open_pnts),                &
                                       k_out(n_e, open_pnts)

SELECT CASE ( pft_xylem_impairment_model(pft) )
CASE (xylem_impairment_whole_trunk)
    CALL psi_leaf_impaired_whole_trunk( pft,                                   &
                                        n_e,                                   &
                                        land_pnts,                             &
                                        open_pnts,                             &
                                        veg_index,                             &
                                        open_index,                            &
                                        e_in,                                  &
                                        psi_root,                              &
                                        kmax,                                  &
                                        kcrit,                                 &
                                        conductance_b,                         &
                                        conductance_c,                         &
                                        psi_leaf_extreme,                      &
                                        psi_root_extreme,                      &
                                      ! INTENT OUT
                                        psi_out,                               &
                                        k_out                                  &
                                        )

CASE (xylem_impairment_memory)
    CALL psi_leaf_impaired_memory( pft,                                        &
                                   n_e,                                        &
                                   land_pnts,                                  &
                                   open_pnts,                                  &
                                   veg_index,                                  &
                                   open_index,                                 &
                                   e_in,                                       &
                                   psi_root,                                   &
                                   kmax_ref,                                   &
                                   kmax,                                       &
                                   kcrit,                                      &
                                   conductance_b,                              &
                                   conductance_c,                              &
                                 ! INTENT OUT
                                   psi_out,                                    &
                                   k_out                                       &
                                   )

END SELECT

END SUBROUTINE impaired_path

END SUBROUTINE psi_leaf_impaired_jls

! ---------------------------------------------------------------------
! Subroutine to manage which xylem impairment model is applied.
! Note the different types of conductance models are not handled here.
!
! MDK, 24/10/26
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment ( n_land_pts                                &
,                                    n_open_pts                                &
,                                    open_index                                &
,                                    pft                                       &
,                                    psi_leaf                                  &
,                                    psi_root                                  &
,                                    psi_stem                                  &
,                                    impaired_k_max                            &
,                                    psi_leaf_extreme                          &
,                                    psi_root_extreme                          &
,                                    lai                                       &
,                                    canht                                     &
,                                    anetc                                     &
                                   )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model

USE xylem_impairment_whole_trunk_mod, ONLY: update_xylem_impairment_whole_trunk
USE xylem_impairment_memory_mod, ONLY: update_xylem_impairment_memory

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  n_land_pts                                                                   &
                            ! Number of land points
, n_open_pts                                                                   &
                            ! Number of land points with open stomata
                            !
                            ! MDK, 24/10/26
, open_index(n_open_pts)                                                       &
                            ! Land point of each point with open stomata
                            !
                            ! MDK, 24/10/26
, pft
                            ! Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  psi_leaf(n_land_pts)                                                         &
                            ! Leaf water potential (Pa)
, psi_root(n_land_pts)                                                         &
                            ! Root water potential (Pa)
, psi_stem(n_land_pts)
                            ! Stem water potential (Pa; memory model driver 4)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
  impaired_k_max(n_land_pts)                                                   &
                            ! Maximum xylem conductance (mol m-2 s-1 Pa-1)
,  psi_leaf_extreme(n_land_pts)                                                &
                            ! Minimum historic leaf water potential (Pa)
, psi_root_extreme(n_land_pts)
                            ! Minimum historic root water potential (Pa)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  lai(n_land_pts)                                                              &
                            ! Leaf area index
, canht(n_land_pts)                                                            &
                            ! Canopy height (m)
, anetc(n_land_pts)
                            ! Canopy net photosynthesis (mol CO2 m-2 s-1)
                            ! (lai, canht and anetc are not used since the
                            ! time-based recovery was retired; kept so the
                            ! sf_stom calls are unchanged.)
                            !
                            ! MDK, 23/10/26


INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not update the xylem impairment model during spinup
IF (is_spinup) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE (pft_xylem_impairment_model(pft))
  CASE(xylem_impairment_none)
    ! No xylem impairment

  CASE(xylem_impairment_whole_trunk)
    CALL update_xylem_impairment_whole_trunk( n_land_pts,                      &
                                              n_open_pts,                      &
                                              open_index,                      &
                                              pft,                             &
                                              psi_leaf,                        &
                                              psi_root,                        &
                                              psi_leaf_extreme,                &
                                              psi_root_extreme                 &
                                            )

  CASE(xylem_impairment_memory)
    CALL update_xylem_impairment_memory( n_land_pts,                           &
                                         n_open_pts,                           &
                                         open_index,                           &
                                         pft,                                  &
                                         psi_leaf,                             &
                                         psi_root,                             &
                                         psi_stem,                             &
                                         impaired_k_max                        &
                                       )

  CASE DEFAULT
    errcode = 1
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), whole_trunk (2) ' //      &
     'or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment

! (canopy_impaired_psi_jls, a separate single-curve solve for the
! psi_canopy_pft and PLC_canopy_pft outputs, was removed: those outputs are
! now psi_leaf_pft and PLC_pft. MDK, 24/10/26)

END MODULE xylem_impairment_mod
