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

PUBLIC :: leaf_conductance_jls, xylem_conductance_jls, psi_leaf_jls

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
! *********************************************************************

! ---------------------------------------------------------------------
! Subroutine to manage leaf conductance calculations for different
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
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point
                            ! (dimensionless)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

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
! Subroutine to manage xylem conductance calculations for different
! conductance models.
! NOTE: This subroutine is designed for the stomatal optimisation model.
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
                            ! (mol m-2 s-1 Pa-1).
, kcrit(open_pnts)                                                             &
                            ! Critical xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1).
, conductance_b(open_pnts)                                                     &
                            ! Conductance parameter b for each open point
                            ! (Pa). JBaguley
, conductance_c(open_pnts)
                            ! Conductance parameter c for each open point
                            ! (dimensionless).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

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
!       It returns an output array that only contains values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_jls( pft,                                                  &
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
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance,               &
                                som_psi_solver, psi_solver_lut,                &
                                psi_solver_newton, l_som_plant_segments,       &
                                l_som_rhizo_series
USE xylem_hydraulics_CW_jls_mod, ONLY: psi_leaf_CW_jls, xylem_f_bc
USE supply_table_mod, ONLY: psi_leaf_lut_jls, lut_pft_curve
USE plant_hydraulic_path_mod, ONLY: ksr_path, som_psi_in_min,                  &
                                    psi_leaf_segments_jls,                     &
                                    psi_leaf_segments_lut_jls
USE xylem_hydraulics_SOX_jls_mod, ONLY: psi_leaf_SOX_jls

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
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 s-1)
, root_zone_psi(land_pts)                                                      &
                             ! Water potential in the root zone (Pa)
, kmax(land_pts)                                                               &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1).
, kcrit(land_pts)                                                              &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1).
, conductance_b(land_pts)                                                      &
                            ! Conductance parameter b for each land point
                            ! (Pa). JBaguley
, conductance_c(land_pts)
                            ! Conductance parameter c for each land point
                            ! (dimensionless).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Whole-path conductance -dE/dpsi_leaf for each
                            ! open point (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: errcode, i, j, l

REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(land_pts)                                                             &
                            ! Root water potential for one sample (Pa).
, e1(1, open_pnts), psi1(1, open_pnts), k1(1, open_pnts)                       &
                            ! One sample per open point.
, k_s, k_in
                            ! Soil-to-root conductance of the path and the
                            ! plant conductance at the root (kmax units).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. l_som_rhizo_series .OR. l_som_plant_segments ) THEN
  ! (With l_som_plant_segments the segment solvers add the soil link
  ! themselves, at the start of their chain.)
  CALL plant_path( n_e_leaf, e_leaf, root_zone_psi, leaf_psi, leaf_k )
ELSE
  !---------------------------------------------------------------------------
  ! Soil (rhizosphere) resistance in series with the plant (SPA, MAESPA).
  !
  ! The soil-to-root drop is linear in E within a step, so each sample has
  ! its own root water potential (where water enters the plant),
  !
  !   psi_in = psi_src - E / K_s
  !
  ! and the plant path runs from there to the leaf. With
  ! psi_l(E, psi_in(E)), the marginal conductance of the whole path is
  !
  !   dE / dpsi_l = k_p / (1 + k(psi_in) / K_s)
  !
  ! where
  !   k_p        plant marginal conductance at the leaf (leaf_k of the
  !              plant path)
  !   k(psi_in)  = kmax f(psi_in), the plant conductance at the root
  !   K_s        soil-to-root conductance
  !
  ! using dpsi_l / dpsi_in = k(psi_in) / k_p for a single conductance
  ! curve. At E = 0 this is the series conductance
  !
  !   1 / (1 / k(psi_src) + 1 / K_s)
  !
  ! MDK, 23/10/26
  !---------------------------------------------------------------------------
  ! (psi_in is set at the open points only: the solvers read no others, and
  ! setting the whole grid for every sample made a one-sample call cost
  ! O(land_pts). MDK, 24/10/26)
  DO i = 1, n_e_leaf
    DO j = 1, open_pnts
      l = veg_index(open_index(j))
      k_s = ksr_path(l)
      psi_in(l) = root_zone_psi(l)
      IF ( ksr_path(l) < 0.0 ) THEN
        ! No soil link for this PFT (fsmc_mod /= 2).
      ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        psi_in(l) = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s, som_psi_in_min)
      ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
        psi_in(l) = som_psi_in_min
      END IF
      e1(1,j) = e_leaf(i,j)
    END DO
    CALL plant_path( 1, e1, psi_in, psi1, k1 )
    DO j = 1, open_pnts
      l = veg_index(open_index(j))
      k_s  = ksr_path(l)
      ! The plant conductance at the root on this point's curve: the
      ! impaired one (conductance_b/c) under xylem impairment, else the PFT
      ! curve. (This used to take supply_lut_f when b and c equalled the
      ! PFT's, which gives the same value for psi <= 0 but built the supply
      ! table even with the Newton solver.)
      !
      ! MDK, 21/10/26
      k_in = kmax(l) * xylem_f_bc(pft, psi_in(l), conductance_b(l),            &
                                  conductance_c(l))
      leaf_psi(i,j) = psi1(1,j)
      IF ( ksr_path(l) < 0.0 ) THEN
        leaf_k(i,j) = k1(1,j)
      ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        leaf_k(i,j) = k1(1,j) / (1.0 + k_in / k_s)
      ELSE
        leaf_k(i,j) = 0.0
      END IF
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! The plant path (root to leaf) for n_e samples per open point.
SUBROUTINE plant_path( n_e, e_in, psi_root, psi_out, k_out )

INTEGER, INTENT(IN) :: n_e
REAL(KIND=real_jlslsm), INTENT(IN)  :: e_in(n_e, open_pnts),                   &
                                       psi_root(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi_out(n_e, open_pnts),                &
                                       k_out(n_e, open_pnts)

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  ! Root, stem and leaf segments in series (l_som_plant_segments) with the
  ! Newton or table solver: leaf psi from solving the segments downstream,
  ! and the whole-plant conductance -dE/dpsi_leaf in k_out. (This choice
  ! was made inside psi_leaf_CW_jls; it is here so that the Weibull module
  ! does not depend on the segment module. With the Taylor solver the
  ! single curve is used, as before.)
  !
  ! MDK, 23/10/26
  IF ( l_som_plant_segments .AND. som_psi_solver == psi_solver_newton ) THEN
    CALL psi_leaf_segments_jls( pft, n_e, land_pts, open_pnts, veg_index,      &
                                open_index, e_in, psi_root, kmax, kcrit,       &
                                psi_out, k_out )
  ELSE IF ( l_som_plant_segments .AND. som_psi_solver == psi_solver_lut ) THEN
    CALL psi_leaf_segments_lut_jls( pft, n_e, land_pts, open_pnts, veg_index,  &
                                    open_index, e_in, psi_root, kmax,          &
                                    psi_out, k_out )

  ! (Only on the PFT's curve: the table cannot follow a point's own b and
  ! c, which psi_leaf_CW_jls then solves by Newton. MDK, 21/10/26)
  ELSE IF ( som_psi_solver == psi_solver_lut .AND.                             &
       lut_pft_curve( pft, land_pts, open_pnts, veg_index, open_index,         &
                      conductance_b, conductance_c ) ) THEN
    ! Direct call: skips psi_leaf_CW_jls's automatic work arrays.
    CALL psi_leaf_lut_jls( pft, n_e, land_pts, open_pnts, veg_index,           &
                           open_index, e_in, psi_root, kmax,                   &
                           psi_out, k_out )
  ELSE
    CALL psi_leaf_CW_jls( pft, n_e, land_pts, open_pnts, veg_index,            &
                          open_index, e_in, psi_root, kmax, kcrit,             &
                          conductance_b, conductance_c,                        &
                       ! INTENT OUT
                          psi_out, k_out )
  END IF

CASE ( SOX_conductance )
  IF ( som_psi_solver == psi_solver_lut .AND.                                  &
       lut_pft_curve( pft, land_pts, open_pnts, veg_index, open_index,         &
                      conductance_b, conductance_c ) ) THEN   ! MDK, 21/10/26
    CALL psi_leaf_lut_jls( pft, n_e, land_pts, open_pnts, veg_index,           &
                           open_index, e_in, psi_root, kmax,                   &
                           psi_out, k_out )
  ELSE
    CALL psi_leaf_SOX_jls( pft, n_e, land_pts, open_pnts, veg_index,           &
                           open_index, e_in, psi_root, kmax, kcrit,            &
                           conductance_b, conductance_c,                       &
                        ! INTENT OUT
                           psi_out, k_out )
  END IF

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

END SUBROUTINE plant_path

END SUBROUTINE psi_leaf_jls

END MODULE xylem_hydraulics_jls_mod