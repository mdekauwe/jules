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

PUBLIC :: xylem_conductance_jls, leaf_psi_jls

! Soil-to-root conductance for l_som_rhizo_series, per land point, as a
! fraction of the whole-plant conductance kmax_pft * LAI: K_s(ground) /
! (kmax_pft * LAI). Set by physiol for the current PFT just before sf_stom;
! leaf_psi_jls gives each leaf path K_s = som_ksr_frac * kmax (its share of
! the soil conductance in proportion to its share of kmax).
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: som_ksr_frac(:)

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
! *********************************************************************

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
, kcrit(open_pnts)
                            ! Critical xylem conductance for each open point
                            ! (m/s).

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
                      ! INTENT OUT
                         leaf_psi,                                             &
                         leaf_k                                                &
  )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance,             &
                                som_psi_solver, psi_solver_lut,           &
                                l_som_plant_segments, l_som_rhizo_series
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_CW_jls, leaf_psi_lut_jls,     &
                                       supply_lut_f
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
, kcrit(land_pts)
                            ! Critical xylem conductance for each land point
                            ! (m/s).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: errcode, i, j, l

REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(land_pts)                                                             &
                            ! Root inlet potential for one sample (Pa).
, e1(1, open_pnts), psi1(1, open_pnts), k1(1, open_pnts)                       &
                            ! One sample per open point.
, k_s, k_in
                            ! Soil-to-root conductance of the path and the
                            ! plant conductance at the inlet (kmax units).

REAL(KIND=real_jlslsm), PARAMETER :: psi_in_min = -1.0e9
                            ! Floor on the inlet potential (Pa) when the
                            ! soil conductance is (near) zero.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. l_som_rhizo_series ) THEN
  CALL plant_path( n_e_leaf, e_leaf, root_zone_psi, leaf_psi, leaf_k )
ELSE
  !---------------------------------------------------------------------------
  ! Soil (rhizosphere) resistance in series with the plant (SPA, MAESPA).
  ! The soil link is linear in E within a step, so each sample has its own
  ! root inlet psi_in = psi_src - E / K_s, and the plant path runs from
  ! there. With psi_l(E, psi_in(E)), the marginal conductance of the whole
  ! path is
  !   dE/dpsi_l = k_p / (1 + k(psi_in) / K_s),
  ! with k_p the plant marginal conductance at the leaf (leaf_k of the plant
  ! path) and k(psi_in) = kmax f(psi_in), using dpsi_l/dpsi_in =
  ! k(psi_in) / k_p for a single conductance curve. At E = 0 this is the
  ! series conductance 1 / (1/k(psi_src) + 1/K_s).
  !---------------------------------------------------------------------------
  DO i = 1, n_e_leaf
    psi_in(:) = root_zone_psi(:)
    DO j = 1, open_pnts
      l = veg_index(open_index(j))
      k_s = som_ksr_frac(l) * kmax(l)
      IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        psi_in(l) = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s, psi_in_min)
      ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
        psi_in(l) = psi_in_min
      END IF
      e1(1,j) = e_leaf(i,j)
    END DO
    CALL plant_path( 1, e1, psi_in, psi1, k1 )
    DO j = 1, open_pnts
      l = veg_index(open_index(j))
      k_s  = som_ksr_frac(l) * kmax(l)
      k_in = kmax(l) * supply_lut_f(pft, psi_in(l))
      leaf_psi(i,j) = psi1(1,j)
      IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        leaf_k(i,j) = k1(1,j) / (1.0 + k_in / k_s)
      ELSE
        leaf_k(i,j) = 0.0
      END IF
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! The plant path (root inlet to leaf) for n_e samples per open point.
SUBROUTINE plant_path( n_e, e_in, psi_root, psi_out, k_out )

INTEGER, INTENT(IN) :: n_e
REAL(KIND=real_jlslsm), INTENT(IN)  :: e_in(n_e, open_pnts),                  &
                                       psi_root(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi_out(n_e, open_pnts),               &
                                       k_out(n_e, open_pnts)

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  IF ( som_psi_solver == psi_solver_lut .AND.                             &
       .NOT. l_som_plant_segments ) THEN
    ! Direct call: skips leaf_psi_CW_jls's automatic work arrays.
    CALL leaf_psi_lut_jls( pft, n_e, land_pts, open_pnts, veg_index,        &
                              open_index, e_in, psi_root, kmax,                &
                              psi_out, k_out )
  ELSE
  CALL leaf_psi_CW_jls( pft, n_e, land_pts, open_pnts, veg_index,             &
                        open_index, e_in, psi_root, kmax, kcrit,               &
                     ! INTENT OUT
                        psi_out, k_out )
  END IF

CASE ( SOX_conductance )
  IF ( som_psi_solver == psi_solver_lut ) THEN
    CALL leaf_psi_lut_jls( pft, n_e, land_pts, open_pnts, veg_index,           &
                           open_index, e_in, psi_root, kmax,                   &
                           psi_out, k_out )
  ELSE
  CALL leaf_psi_SOX_jls( pft, n_e, land_pts, open_pnts, veg_index,             &
                         open_index, e_in, psi_root, kmax, kcrit,              &
                      ! INTENT OUT
                         psi_out, k_out )
  END IF

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

END SUBROUTINE plant_path

END SUBROUTINE leaf_psi_jls

END MODULE xylem_hydraulics_jls_mod