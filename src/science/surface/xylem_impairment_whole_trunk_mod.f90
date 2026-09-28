! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_whole_trunk_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_MOD'

PUBLIC :: leaf_conductance_impaired_whole_trunk_jls,                           &
          xylem_conductance_impaired_whole_trunk_stom_opt_jls,                 &
          leaf_psi_impaired_whole_trunk,                                       &
          update_xylem_impairment_whole_trunk

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calcualte leaf conductance
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_impaired_whole_trunk_jls( pft,                    &
                                                      land_pnts,              &
                                                      leaf_psi,               &
                                                      kmax,                   &
                                                      kcrit,                  &
                                                      conductance_b,          &
                                                      conductance_c,          &
                                                      psi_leaf_extreme,       &
                                                   ! INTENT OUT
                                                      leaf_k                  &
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
  leaf_psi(land_pnts)                                                          &
                            ! Leaf water potential for each land point (Pa)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point (m/s)

! Local variables
REAL(KIND=real_jlslsm) :: new_psi_leaf_extreme(land_pnts)
                            ! Holds the historic minimum leaf water potential
                            ! for each land point.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_IMPAIRED_WHOLE_TRUNK_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the new historic minimum leaf water potentials, incase the leaf
! water potential is less than the historic minimum leaf water potential.
new_psi_leaf_extreme = MIN(leaf_psi, psi_leaf_extreme)

! The leaf conductance is just the conductance for the historic minimum leaf
! water potential.
! NOTE: kmax/conductance_b/conductance_c are the per-point arrays rather than
!       kmax_pft/conductance_b_pft/conductance_c_pft. The whole-trunk model
!       never updates them, so they hold the PFT values - but passed this way
!       they carry whatever canopy scaling (e.g. fpar) the caller applied.
CALL leaf_conductance_jls(pft,                                                &
                          land_pnts,                                          &
                          new_psi_leaf_extreme,                               &
                          kmax,                                               &
                          kcrit,                                              &
                          conductance_b,                                      &
                          conductance_c,                                      &
                        ! INTENT OUT
                          leaf_k                                              &
                         )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_whole_trunk_jls

! ---------------------------------------------------------------------
! Function to manage xylem conductance calculations for different
! impairment models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It asumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_impaired_whole_trunk_stom_opt_jls(                &
                                                   pft,                        &
                                                   n_water_potentials,         &
                                                   open_pnts,                  &
                                                   open_index,                 &
                                                   veg_index,                  &
                                                   land_pnts,                  &
                                                   water_potential,            &
                                                   kmax,                       &
                                                   kcrit,                      &
                                                   conductance_b,              &
                                                   conductance_c,              &
                                                   psi_extreme,                &
                                                 ! INTENT OUT
                                                   xylem_conductance           &
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
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_extreme(land_pnts)
                            ! Historic minimum water potential at a given point
                            ! in the tree for each land point (Pa).
                            ! To calculate the leaf (root) conductance this
                            ! value is the historic minm=imum leaf (root) water
                            ! potential.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
REAL(KIND=real_jlslsm) :: new_psi_extreme(n_water_potentials, open_pnts)
REAL(KIND=real_jlslsm) ::                                                      &
  kmax_open(open_pnts), kcrit_open(open_pnts), b_open(open_pnts),             &
  c_open(open_pnts)
                            ! Per-point curve parameters gathered onto the
                            ! open-point index.
INTEGER :: i, l

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_WHOLE_TRUNK_STOM_OPT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Fill the new_psi_extreme array with the historic minimum water potentials
! for each open point.
! NOTE: open_index indexes veg_index, not the land grid directly.
DO i = 1, open_pnts
  l = veg_index(open_index(i))
  new_psi_extreme(:,i) = psi_extreme(l)
  kmax_open(i)  = kmax(l)
  kcrit_open(i) = kcrit(l)
  b_open(i)     = conductance_b(l)
  c_open(i)     = conductance_c(l)
END DO

! Calculate the new historic minimum water potentials.
new_psi_extreme = MIN(water_potential, new_psi_extreme)

! Calculate the xylem conductance for each open point. See the note in
! leaf_conductance_impaired_whole_trunk_jls on why the per-point curve
! parameters are used rather than the PFT ones.
CALL xylem_conductance_jls(pft,                                                &
                           n_water_potentials,                                 &
                           open_pnts,                                          &
                           new_psi_extreme,                                    &
                           kmax_open,                                          &
                           kcrit_open,                                         &
                           b_open,                                             &
                           c_open,                                             &
                         ! INTENT OUT
                           xylem_conductance                                   &
                          )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_impaired_whole_trunk_stom_opt_jls

! ---------------------------------------------------------------------
! Subroutine to manage leaf water potential calculations for different
! conductance models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_impaired_whole_trunk( pft,                                 &
                                          n_e_leaf,                            &
                                          land_pnts,                           &
                                          open_pnts,                           &
                                          veg_index,                           &
                                          open_index,                          &
                                          e_leaf,                              &
                                          root_zone_psi,                       &
                                          kmax,                                &
                                          kcrit,                               &
                                          conductance_b,                       &
                                          conductance_c,                       &
                                          psi_leaf_extreme,                    &
                                          psi_root_extreme,                    &
                                       ! INTENT OUT
                                          leaf_psi,                            &
                                          leaf_k                               &
  )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance

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
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)                                                  &
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)
, psi_root_extreme(land_pnts)
                            ! Historic minimum root water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
REAL(KIND=real_jlslsm) :: new_psi_root_extreme(land_pnts)
                            ! Holds the historic minimum root water potential
                            ! for each land point.
REAL(KIND=real_jlslsm) :: new_psi_leaf_extreme(land_pnts)
                            ! Holds the historic minimum leaf water potential
                            ! for each land point.

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_WHOLE_TRUNK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Get the new historic minimum root water potentials.
! Note this is internal to the function and doesn't change the input.
! Changing the minimum root water potential is done in the subroutine
! update_xylem_impairment_whole_trunk.
new_psi_root_extreme = MIN(root_zone_psi, psi_root_extreme)

! The leaf water potential can never be greater than the root water potential,
! hence we need to impose that the historic minimum leaf water potential is
! less than the historic minimum root water potential.
new_psi_leaf_extreme = MIN(psi_leaf_extreme, new_psi_root_extreme)

! Switch between different conductance models
SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  CALL leaf_psi_impaired_CW_jls(pft,                                         &
                                n_e_leaf,                                    &
                                land_pnts,                                   &
                                open_pnts,                                   &
                                veg_index,                                   &
                                open_index,                                  &
                                e_leaf,                                      &
                                root_zone_psi,                               &
                                kmax,                                        &
                                kcrit,                                       &
                                conductance_b,                               &
                                conductance_c,                               &
                                new_psi_leaf_extreme,                        &
                                new_psi_root_extreme,                        &
                              ! INTENT OUT
                                leaf_psi,                                    &
                                leaf_k                                       &
                               )

CASE ( SOX_conductance )
  CALL leaf_psi_impaired_SOX_jls(pft,                                        &
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
                                 new_psi_leaf_extreme,                       &
                                 new_psi_root_extreme,                       &
                               ! INTENT OUT
                                 leaf_psi,                                   &
                                 leaf_k                                      &
                                )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_impaired_whole_trunk

! ------------------------------------------------------------------------
! Subroutine to update the historic minimum leaf and root water potentials
! ------------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_whole_trunk ( n_land_pts                    &
,                                                n_open_pts                    &
,                                                open_index                    &
,                                                pft                           &
,                                                leaf_psi                      &
,                                                root_zone_psi                 &
,                                                psi_leaf_extreme              &
,                                                psi_root_extreme              &
                                               )

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
  leaf_psi(n_land_pts)                                                         &
                            ! Leaf water potential (Pa)
, root_zone_psi(n_land_pts)
                            ! Root zone water potential (Pa)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  psi_leaf_extreme(n_land_pts)                                                 &
                            ! Historic minimum leaf water potential (Pa)
, psi_root_extreme(n_land_pts)
                            ! Historic minimum root water potential (Pa)

! Local variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_WHOLE_TRUNK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Update the historic minimum leaf and root water potentials
psi_leaf_extreme = MIN(leaf_psi, psi_leaf_extreme)
psi_root_extreme = MIN(root_zone_psi, psi_root_extreme)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_whole_trunk

! ---------------------------------------------------------------------
! Internal subroutines to calculate leaf water potential from
! transpiration and root zone water potential.
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Whole-trunk impairment: the conductance at each point along the xylem
! is that of the unimpaired curve at the lowest water potential that
! point has experienced. The historic minimum profile runs from
! psi_root_extreme (psi_re) at the roots to psi_leaf_extreme (psi_le) at
! the leaf, and is taken to have been a steady-state flow profile, so
! with G(psi) = integral of k(psi') dpsi' from psi to psi_re (on the
! unimpaired curve):
!
!  Case 1 - leaf stays above its historic minimum (psi_l >= psi_le): the
!   whole path keeps its historic conductance, so
!     e = keff * (psi_r - psi_l),   keff = G(psi_le) / (psi_re - psi_le)
!   (keff -> k(psi_re) as psi_le -> psi_re).
!
!  Case 2 - leaf goes below its historic minimum (psi_l < psi_le): the
!   historic profile is reset to run from psi_re to psi_l, so
!     E(psi_l) = G(psi_l) * (psi_r - psi_l) / (psi_re - psi_l)
!   solved for E(psi_l) = e by Newton-Raphson, with
!     dE/dpsi_l = - k(psi_l) (psi_r - psi_l) / (psi_re - psi_l)
!                 + G(psi_l) (psi_r - psi_re) / (psi_re - psi_l)^2
!
! Case 1 applies when its solution satisfies psi_l >= psi_le, otherwise
! case 2 (E is continuous at psi_l = psi_le, where both give the same e).
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf water potential for the Cumulative
! Weibull model.
! ---------------------------------------------------------------------

SUBROUTINE leaf_psi_impaired_CW_jls( pft,                                      &
                                     n_e_leaf,                                 &
                                     land_pnts,                                &
                                     open_pnts,                                &
                                     veg_index,                                &
                                     open_index,                               &
                                     e_leaf,                                   &
                                     root_zone_psi,                            &
                                     kmax_pft,                                 &
                                     kcrit,                                    &
                                     conductance_b,                            &
                                     conductance_c,                            &
                                     psi_leaf_extreme,                         &
                                     psi_root_extreme,                         &
                                  ! INTENT OUT
                                     leaf_psi,                                 &
                                     leaf_k                                    &
  )

USE xylem_hydraulics_CW_jls_mod, ONLY: incomplete_gamma

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
, kmax_pft(land_pnts)                                                          &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)                                                  &
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)
, psi_root_extreme(land_pnts)
                            ! Historic minimum root water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
REAL(KIND=real_jlslsm) :: gamma_re
                            ! Lower incomplete gamma function at psi_re.
REAL(KIND=real_jlslsm) :: psi_r, psi_re, psi_le
                            ! Root zone, historic minimum root and historic
                            ! minimum leaf water potentials for this point
                            ! (Pa).
REAL(KIND=real_jlslsm) :: g_le
                            ! G(psi_le): transpiration along the historic
                            ! minimum profile.
REAL(KIND=real_jlslsm) :: keff
                            ! Effective whole-path conductance for case 1.
REAL(KIND=real_jlslsm) :: psi_c1(n_e_leaf)
                            ! Case 1 leaf water potential (Pa).
REAL(KIND=real_jlslsm) :: psi_vec(n_e_leaf), g_vec(n_e_leaf),               &
                          k_vec(n_e_leaf), d_vec(n_e_leaf),                 &
                          e_vec(n_e_leaf), de_vec(n_e_leaf),                &
                          dpsi_vec(n_e_leaf)
                            ! Case 2 Newton-Raphson work arrays.
REAL(KIND=real_jlslsm) :: psi_floor
                            ! Lower bound on psi_l for case 2, to avoid
                            ! runaway steps for infeasible samples (as in
                            ! the base leaf_psi routines).

INTEGER, PARAMETER :: max_nr_iter = 10
                            ! Cap on case 2 Newton-Raphson iterations.
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Case 2 convergence tolerance on psi_l (Pa).
REAL(KIND=real_jlslsm), PARAMETER :: dpsi_min = 1.0
                            ! Minimum psi_re - psi_le (Pa) below which the
                            ! historic profile is treated as flat.

INTEGER :: i, j, l

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j = 1, open_pnts
  l = veg_index(open_index(j))

  psi_r  = root_zone_psi(l)
  psi_re = psi_root_extreme(l)
  psi_le = psi_leaf_extreme(l)

  gamma_re = incomplete_gamma_scalar(1/conductance_c(l),                   &
                   (psi_re/conductance_b(l))**conductance_c(l))

  !-------------------------------------------------------------------------
  ! Case 1: leaf above its historic minimum.
  !-------------------------------------------------------------------------
  IF (psi_re - psi_le > dpsi_min) THEN
    g_le = kmax_pft(l) * (-conductance_b(l)/conductance_c(l))              &
           * ( incomplete_gamma_scalar(1/conductance_c(l),                 &
                   (psi_le/conductance_b(l))**conductance_c(l)) - gamma_re )
    keff = g_le / (psi_re - psi_le)
  ELSE
    ! Flat historic profile: the whole path has the conductance at psi_re.
    keff = kmax_pft(l) * EXP( -(psi_re/conductance_b(l))**conductance_c(l) )
  END IF
  keff = MAX(keff, TINY(1.0_real_jlslsm))

  psi_c1(:) = psi_r - e_leaf(:,j) / keff

  !-------------------------------------------------------------------------
  ! Case 2: leaf below its historic minimum. Solved for every sample and
  ! only used where case 1 does not apply.
  !-------------------------------------------------------------------------
  psi_floor = psi_r - 5.0 * ABS(conductance_b(l))
  psi_vec(:) = MAX(MIN(psi_c1(:), psi_le), psi_floor)

  DO i = 1, max_nr_iter
    g_vec(:) = kmax_pft(l) * (-conductance_b(l)/conductance_c(l))            &
             * ( incomplete_gamma(n_e_leaf, 1/conductance_c(l),              &
                     (psi_vec(:)/conductance_b(l))**conductance_c(l))        &
                 - gamma_re )
    k_vec(:) = kmax_pft(l)                                                   &
             * EXP( -(psi_vec(:)/conductance_b(l))**conductance_c(l) )

    ! Distance below psi_re, floored to avoid dividing by zero where the
    ! leaf sits at psi_re (zero transpiration).
    d_vec(:) = MAX(psi_re - psi_vec(:), dpsi_min)

    e_vec(:)  = g_vec(:) * (psi_r - psi_vec(:)) / d_vec(:)
    de_vec(:) = - k_vec(:) * (psi_r - psi_vec(:)) / d_vec(:)                 &
                + g_vec(:) * (psi_r - psi_re) / d_vec(:)**2

    ! dE/dpsi_l < 0 for any physical solution; guard the step against a
    ! non-negative derivative.
    de_vec(:) = MIN(de_vec(:), -TINY(1.0_real_jlslsm))

    dpsi_vec(:) = - (e_vec(:) - e_leaf(:,j)) / de_vec(:)
    psi_vec(:)  = MIN(MAX(psi_vec(:) + dpsi_vec(:), psi_floor), psi_le)

    IF (MAXVAL(ABS(dpsi_vec(:))) < psi_tol) EXIT
  END DO

  ! Case 1 applies where its solution stays at or above psi_le.
  WHERE (psi_c1(:) >= psi_le)
    leaf_psi(:,j) = psi_c1(:)
  ELSEWHERE
    leaf_psi(:,j) = psi_vec(:)
  END WHERE
END DO

CALL xylem_conductance_impaired_whole_trunk_stom_opt_jls(pft,                  &
                                            n_e_leaf,                          &
                                            open_pnts,                         &
                                            open_index,                        &
                                            veg_index,                         &
                                            land_pnts,                         &
                                            leaf_psi,                          &
                                            kmax_pft,                          &
                                            kcrit,                             &
                                            conductance_b,                     &
                                            conductance_c,                     &
                                            psi_leaf_extreme,                  &
                                          ! INTENT OUT
                                            leaf_k                             &
                                           )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_impaired_CW_jls

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf water potential for the SOX model.
! ---------------------------------------------------------------------

SUBROUTINE leaf_psi_impaired_SOX_jls( pft,                                    &
                                      n_e_leaf,                               &
                                      land_pnts,                              &
                                      open_pnts,                              &
                                      veg_index,                              &
                                      open_index,                             &
                                      e_leaf,                                 &
                                      root_zone_psi,                          &
                                      kmax_pft,                               &
                                      kcrit,                                  &
                                      conductance_b,                          &
                                      conductance_c,                          &
                                      psi_leaf_extreme,                       &
                                      psi_root_extreme,                       &
                                   ! INTENT OUT
                                      leaf_psi,                               &
                                      leaf_k                                  &
  )

USE xylem_hydraulics_SOX_jls_mod, ONLY: SOX_2F1

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
, kmax_pft(land_pnts)                                                          &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)                                                  &
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)
, psi_root_extreme(land_pnts)
                            ! Historic minimum root water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
REAL(KIND=real_jlslsm) :: gamma_re
                            ! psi_re * 2F1(...) at psi_re (the psi_re end of
                            ! the transpiration integral).
REAL(KIND=real_jlslsm) :: psi_r, psi_re, psi_le
                            ! Root zone, historic minimum root and historic
                            ! minimum leaf water potentials for this point
                            ! (Pa).
REAL(KIND=real_jlslsm) :: g_le
                            ! G(psi_le): transpiration along the historic
                            ! minimum profile.
REAL(KIND=real_jlslsm) :: keff
                            ! Effective whole-path conductance for case 1.
REAL(KIND=real_jlslsm) :: psi_c1(n_e_leaf)
                            ! Case 1 leaf water potential (Pa).
REAL(KIND=real_jlslsm) :: psi_vec(n_e_leaf), g_vec(n_e_leaf),               &
                          k_vec(n_e_leaf), d_vec(n_e_leaf),                 &
                          e_vec(n_e_leaf), de_vec(n_e_leaf),                &
                          dpsi_vec(n_e_leaf)
                            ! Case 2 Newton-Raphson work arrays.
REAL(KIND=real_jlslsm) :: psi_floor
                            ! Lower bound on psi_l for case 2, to avoid
                            ! runaway steps for infeasible samples (as in
                            ! the base leaf_psi routines).

INTEGER, PARAMETER :: max_nr_iter = 10
                            ! Cap on case 2 Newton-Raphson iterations.
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Case 2 convergence tolerance on psi_l (Pa).
REAL(KIND=real_jlslsm), PARAMETER :: dpsi_min = 1.0
                            ! Minimum psi_re - psi_le (Pa) below which the
                            ! historic profile is treated as flat.

INTEGER :: i, j, l

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO j = 1, open_pnts
  l = veg_index(open_index(j))

  psi_r  = root_zone_psi(l)
  psi_re = psi_root_extreme(l)
  psi_le = psi_leaf_extreme(l)

  gamma_re = psi_re * sox_2f1_scalar(1 + 1/conductance_c(l),               &
                   -(psi_re/conductance_b(l))**conductance_c(l))

  !-------------------------------------------------------------------------
  ! Case 1: leaf above its historic minimum.
  !-------------------------------------------------------------------------
  IF (psi_re - psi_le > dpsi_min) THEN
    g_le = - kmax_pft(l)                                                   &
           * ( psi_le * sox_2f1_scalar(1 + 1/conductance_c(l),             &
                   -(psi_le/conductance_b(l))**conductance_c(l)) - gamma_re )
    keff = g_le / (psi_re - psi_le)
  ELSE
    ! Flat historic profile: the whole path has the conductance at psi_re.
    keff = kmax_pft(l) / (1 + (psi_re/conductance_b(l))**conductance_c(l))
  END IF
  keff = MAX(keff, TINY(1.0_real_jlslsm))

  psi_c1(:) = psi_r - e_leaf(:,j) / keff

  !-------------------------------------------------------------------------
  ! Case 2: leaf below its historic minimum. Solved for every sample and
  ! only used where case 1 does not apply.
  !-------------------------------------------------------------------------
  psi_floor = psi_r - 5.0 * ABS(conductance_b(l))
  psi_vec(:) = MAX(MIN(psi_c1(:), psi_le), psi_floor)

  DO i = 1, max_nr_iter
    g_vec(:) = - kmax_pft(l)                                                 &
             * ( psi_vec(:) * SOX_2F1(n_e_leaf, 1 + 1/conductance_c(l),      &
                     -(psi_vec(:)/conductance_b(l))**conductance_c(l))       &
                 - gamma_re )
    k_vec(:) = kmax_pft(l)                                                   &
             / (1 + (psi_vec(:)/conductance_b(l))**conductance_c(l))

    ! Distance below psi_re, floored to avoid dividing by zero where the
    ! leaf sits at psi_re (zero transpiration).
    d_vec(:) = MAX(psi_re - psi_vec(:), dpsi_min)

    e_vec(:)  = g_vec(:) * (psi_r - psi_vec(:)) / d_vec(:)
    de_vec(:) = - k_vec(:) * (psi_r - psi_vec(:)) / d_vec(:)                 &
                + g_vec(:) * (psi_r - psi_re) / d_vec(:)**2

    ! dE/dpsi_l < 0 for any physical solution; guard the step against a
    ! non-negative derivative.
    de_vec(:) = MIN(de_vec(:), -TINY(1.0_real_jlslsm))

    dpsi_vec(:) = - (e_vec(:) - e_leaf(:,j)) / de_vec(:)
    psi_vec(:)  = MIN(MAX(psi_vec(:) + dpsi_vec(:), psi_floor), psi_le)

    IF (MAXVAL(ABS(dpsi_vec(:))) < psi_tol) EXIT
  END DO

  ! Case 1 applies where its solution stays at or above psi_le.
  WHERE (psi_c1(:) >= psi_le)
    leaf_psi(:,j) = psi_c1(:)
  ELSEWHERE
    leaf_psi(:,j) = psi_vec(:)
  END WHERE
END DO

CALL xylem_conductance_impaired_whole_trunk_stom_opt_jls(                       &
                                            pft,                                &
                                            n_e_leaf,                           &
                                            open_pnts,                          &
                                            open_index,                         &
                                            veg_index,                          &
                                            land_pnts,                          &
                                            leaf_psi,                           &
                                            kmax_pft,                           &
                                            kcrit,                              &
                                            conductance_b,                      &
                                            conductance_c,                      &
                                            psi_leaf_extreme,                   &
                                          ! INTENT OUT
                                            leaf_k                              &
                                           )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_impaired_SOX_jls

! ---------------------------------------------------------------------
! Scalar wrappers around the array series functions used by the
! whole-trunk solvers above.
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma_scalar( a, x ) RESULT( gamma )

USE xylem_hydraulics_CW_jls_mod, ONLY: incomplete_gamma

REAL(KIND=real_jlslsm), INTENT(IN) :: a, x
REAL(KIND=real_jlslsm) :: gamma
REAL(KIND=real_jlslsm) :: gamma_arr(1)

gamma_arr = incomplete_gamma(1, a, (/ x /))
gamma = gamma_arr(1)

END FUNCTION incomplete_gamma_scalar

FUNCTION sox_2f1_scalar( cf, xf ) RESULT( sum )

USE xylem_hydraulics_SOX_jls_mod, ONLY: SOX_2F1

REAL(KIND=real_jlslsm), INTENT(IN) :: cf, xf
REAL(KIND=real_jlslsm) :: sum
REAL(KIND=real_jlslsm) :: sum_arr(1)

sum_arr = SOX_2F1(1, cf, (/ xf /))
sum = sum_arr(1)

END FUNCTION sox_2f1_scalar

END MODULE xylem_impairment_whole_trunk_mod