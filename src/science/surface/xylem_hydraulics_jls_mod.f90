! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_hydraulics_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_HYDRAULICS_JLS_MOD'

PUBLIC :: leaf_conductance_jls, xylem_conductance_jls, leaf_psi_jls

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
! *********************************************************************

! ---------------------------------------------------------------------
! Function to manage leaf conductance calculations for different
! conductance models, for all land points. JBaguley
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_jls( pft,                                          &
                                 land_pnts,                                    &
                                 water_potential,                              &
                                 kmax,                                         &
                                 kcrit,                                        &
                                 conductance_b,                                &
                                 conductance_c,                                &
                               ! INTENT OUT
                                 leaf_k                                        &
  )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_conductance_CW_jls
USE xylem_hydraulics_SOX_jls_mod, ONLY: leaf_conductance_SOX_jls

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
                            ! Water potentials for each land point (Pa)
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
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point (m/s)

! Local variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  CALL leaf_conductance_CW_jls( pft, land_pnts, water_potential,               &
                                kmax, kcrit, conductance_b, conductance_c,     &
                              ! INTENT OUT
                                leaf_k )

CASE ( SOX_conductance )
  CALL leaf_conductance_SOX_jls( pft, land_pnts, water_potential,              &
                                 kmax, kcrit, conductance_b, conductance_c,    &
                               ! INTENT OUT
                                 leaf_k )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE leaf_conductance_jls

! ---------------------------------------------------------------------
! Function to manage xylem conductance calculations for different
! conductance models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It assumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_jls( pft,                                         &
                                  n_water_potentials,                          &
                                  open_pnts,                                   &
                                  water_potential,                             &
                                  kmax,                                        &
                                  kcrit,                                       &
                                  conductance_b,                               &
                                  conductance_c,                               &
                                ! INTENT OUT
                                  xylem_conductance                            &
  )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance
USE xylem_hydraulics_CW_jls_mod, ONLY: xylem_conductance_CW_jls
USE xylem_hydraulics_SOX_jls_mod, ONLY: xylem_conductance_SOX_jls

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_water_potentials                                                           &
                            ! Number of water potentials per open point
, open_pnts
                            ! Number of open stomata

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(n_water_potentials, open_pnts)                               &
                            ! Water potentials for each open point (Pa)
, kmax(open_pnts)                                                              &
                            ! Maximum xylem conductance for each open point
                            ! (m/s).
, kcrit(open_pnts)                                                             &
                            ! Critical xylem conductance for each open point
                            ! (m/s).
, conductance_b(open_pnts)                                                     &
                            ! Conductance parameter b for each open point
                            ! (Pa). JBaguley
, conductance_c(open_pnts)
                            ! Conductance parameter c for each open point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  CALL xylem_conductance_CW_jls( pft,                                          &
                                 n_water_potentials,                           &
                                 open_pnts,                                    &
                                 water_potential,                              &
                                 kmax,                                         &
                                 kcrit,                                        &
                                 conductance_b,                                &
                                 conductance_c,                                &
                               ! INTENT OUT
                                 xylem_conductance                             &
                                 )
CASE ( SOX_conductance )
  CALL xylem_conductance_SOX_jls( pft,                                         &
                                  n_water_potentials,                          &
                                  open_pnts,                                   &
                                  water_potential,                             &
                                  kmax,                                        &
                                  kcrit,                                       &
                                  conductance_b,                               &
                                  conductance_c,                               &
                                ! INTENT OUT
                                  xylem_conductance                            &
                                  )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_jls

! ---------------------------------------------------------------------
! Subroutine to manage leaf water potential calculations for different
! conductance models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_jls( pft,                                                  &
                         n_e_leaf,                                             &
                         land_pts,                                             &
                         open_pnts,                                            &
                         veg_index,                                            &
                         open_index,                                           &
                         e_leaf,                                               &
                         root_zone_psi,                                        &
                         kmax,                                                 &
                         kcrit,                                                &
                         conductance_b,                                        &
                         conductance_c,                                        &
                      ! INTENT OUT
                         leaf_psi,                                             &
                         leaf_k                                                &
  )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance,             &
                                som_psi_aprox_method, psi_aprox_LUT,           &
                                l_som_plant_segments
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_CW_jls, leaf_psi_lut_jls
USE xylem_hydraulics_SOX_jls_mod, ONLY: leaf_psi_SOX_jls

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_e_leaf                                                                     &
                            ! Number of transpirations per land point
, land_pts                                                                     &
                            ! Number of land points
, open_pnts                                                                    &
                            ! Number of open stomata
, veg_index(land_pts)                                                          &
                            ! Index of vegetation points on the land grid
, open_index(land_pts)
                            ! Index of open stomata for each land point

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point (kg/m2/s)
, root_zone_psi(land_pts)                                                      &
                             ! Water potential in the root zone (Pa)
, kmax(land_pts)                                                               &
                            ! Maximum xylem conductance for each land point
                            ! (m/s).
, kcrit(land_pts)                                                              &
                            ! Critical xylem conductance for each land point
                            ! (m/s).
, conductance_b(land_pts)                                                      &
                            ! Conductance parameter b for each land point
                            ! (Pa). JBaguley
, conductance_c(land_pts)
                            ! Conductance parameter c for each land point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  IF ( som_psi_aprox_method == psi_aprox_LUT .AND.                             &
       .NOT. l_som_plant_segments ) THEN
    ! Direct call: skips leaf_psi_CW_jls's automatic work arrays.
    CALL leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,   &
                              open_index, e_leaf, root_zone_psi, kmax,         &
                              leaf_psi, leaf_k )
  ELSE
  CALL leaf_psi_CW_jls( pft,                                                   &
                        n_e_leaf,                                              &
                        land_pts,                                              &
                        open_pnts,                                             &
                        veg_index,                                             &
                        open_index,                                            &
                        e_leaf,                                                &
                        root_zone_psi,                                         &
                        kmax,                                                  &
                        kcrit,                                                 &
                        conductance_b,                                         &
                        conductance_c,                                         &
                     ! INTENT OUT
                        leaf_psi,                                              &
                        leaf_k                                                 &
                        )
  END IF

CASE ( SOX_conductance )
  IF ( som_psi_aprox_method == psi_aprox_LUT ) THEN
    CALL leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,      &
                           open_index, e_leaf, root_zone_psi, kmax,            &
                           leaf_psi, leaf_k )
  ELSE
  CALL leaf_psi_SOX_jls( pft,                                                  &
                         n_e_leaf,                                             &
                         land_pts,                                             &
                         open_pnts,                                            &
                         veg_index,                                            &
                         open_index,                                           &
                         e_leaf,                                               &
                         root_zone_psi,                                        &
                         kmax,                                                 &
                         kcrit,                                                &
                         conductance_b,                                        &
                         conductance_c,                                        &
                      ! INTENT OUT
                         leaf_psi,                                             &
                         leaf_k                                                &
                         )
  END IF

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_jls

END MODULE xylem_hydraulics_jls_mod