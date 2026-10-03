! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_hydraulics_CW_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='XYLEM_HYDRAULICS_CW_JLS_MOD'

PUBLIC :: xylem_conductance_CW_jls, leaf_psi_CW_jls, leaf_psi_lut_jls,        &
          supply_lut_psi, supply_lut_e_crit, supply_lut_f

! ---------------------------------------------------------------------
! Supply-function lookup table (som_psi_solver = psi_solver_lut),
! for either conductance model of the PFT (pft_conductance_model):
!   cumulative Weibull: f(psi) = exp(-(psi/b)^c)
!   SOX:                f(psi) = 1 / (1 + (psi/b)^c)
! With k(psi) = kmax * f(psi), the transpiration between root zone and
! leaf is E = kmax * (S(psi_l) - S(psi_r)), S(psi) = integral(f, psi, 0).
! S depends only on the PFT's b and c (not on kmax, which can vary in
! space and time), so it is tabulated once per PFT on a uniform psi grid
! from 0 down to where f has fallen to well below kcrit/kmax. psi_leaf is
! then found by inverting S directly (binary search plus linear
! interpolation) instead of Newton-Raphson on the incomplete gamma or
! hypergeometric series.
! ---------------------------------------------------------------------
INTEGER, PARAMETER, PRIVATE :: n_lut = 4001
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor = 1.0e-6
                            ! Bottom of the Weibull table (f falls fast).
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor_kcrit = 0.1
                            ! Bottom of the SOX table, as a fraction of
                            ! kcrit/kmax (f has a power-law tail, so a fixed
                            ! tiny floor would spread the table too thin).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE :: lut_s(:,:)
                            ! S(psi_i) (Pa), psi_i = -(i-1)*lut_dpsi.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE :: lut_dpsi(:)
                            ! Grid spacing for each PFT (Pa).
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_ready(:)
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_sox(:)
                            ! Table is for the SOX curve (else Weibull).
! Root / stem / leaf segment tables (l_som_plant_segments with
! psi_solver_lut): the same S(psi) for each segment's own Weibull
! (conductance_b_seg, conductance_c_seg), so the segments in series are
! solved by three table inversions (leaf_psi_segments_lut_jls).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE :: lut_seg_s(:,:,:)
                            ! S(psi_i) (Pa) for (i, pft, segment).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE :: lut_seg_dpsi(:,:)
                            ! Grid spacing for (pft, segment) (Pa).
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_seg_ready(:)

CONTAINS

! *********************************************************************
! Contains routines used to calculate conductance and leaf water
! potential for a cumulative Weibull conductance model.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calculate xylem conductance using a cumulative weibull
! distribution.
! NOTE: This function is designed for the stomatal optimisation model.
!       It assumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_CW_jls( pft,                                      &
                                     n_water_potentials,                       &
                                     open_pnts,                                &
                                     water_potential,                          &
                                     kmax,                                     &
                                     kcrit,                                    &
                                  ! INTENT OUT
                                     xylem_conductance                         &
  )

USE pftparm, ONLY:                                                             &
        conductance_b, conductance_c

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
                            ! Water potentials for each open point)
, kmax(open_pnts)                                                              &
                            ! Maximum xylem conductance for each open point
                            ! (m/s). Lets a layer- or PFT-specific maximum
                            ! be used rather than always reading kmax_pft.
, kcrit(open_pnts)
                            ! Critical xylem conductance for each open point
                            ! (m/s), scaled to match whatever kmax is being
                            ! used (e.g. kcrit_per_lyr alongside kmax_per_lyr).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

xylem_conductance = SPREAD(kmax, DIM = 1, NCOPIES = n_water_potentials)        &
     * EXP(-ABS(water_potential / conductance_b(pft))                          &
                ** conductance_c(pft))

! Apply minimum conductance limit
xylem_conductance = MAX(xylem_conductance,                                     &
                        SPREAD(kcrit, DIM = 1, NCOPIES = n_water_potentials))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_CW_jls



! ---------------------------------------------------------------------
! Function to calculate the leaf water potential from the transpiration
! rate and rootzone water potential. Uses a cumulative Weibull
! distribution to model conductance.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! The code bellow uses one of two approximations to estimate the leaf
! water potential from the transpiration rate. This is done to avoid
! the need to integrate the conductance equation each iteration.
!   1: Applies a zeroth order taylor expansion to the conductance
!      equation (psi_solver_taylor).
!   2: Uses a Newton Raphson approximation to find the leaf water
!      potential (psi_solver_newton).
!
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_CW_jls( pft,                                               &
                            n_e_leaf,                                          &
                            land_pts,                                          &
                            open_pnts,                                         &
                            veg_index,                                         &
                            open_index,                                        &
                            e_leaf,                                            &
                            root_zone_psi,                                     &
                            kmax,                                              &
                            kcrit,                                             &
                         ! INTENT OUT
                            leaf_psi,                                          &
                            leaf_k                                             &
  )

USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_taylor,            &
                                psi_solver_newton, psi_solver_lut,                   &
                                l_som_plant_segments

USE pftparm, ONLY: conductance_b, conductance_c

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_e_leaf                                                                     &
                            ! Number of transpiration values per land point
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
, kmax(land_pts)                                                                &
                            ! Maximum xylem conductance for each land point
                            ! (m/s), e.g. kmax_pft(pft) or a layer-scaled
                            ! value such as kmax_per_lyr.
, kcrit(land_pts)
                            ! Critical xylem conductance for each land point
                            ! (m/s), scaled to match whatever kmax is being
                            ! used (e.g. kcrit_per_lyr alongside kmax_per_lyr).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: i, j, k, l

! Taylor expansion variables
REAL(KIND=real_jlslsm) :: leaf_k_new
REAL(KIND=real_jlslsm) :: k_conversion_limit
REAL(KIND=real_jlslsm) :: reference_psi
REAL(KIND=real_jlslsm) :: kmax_open(open_pnts)
                            ! kmax gathered onto the open-point index, for
                            ! passing into xylem_conductance_CW_jls.
REAL(KIND=real_jlslsm) :: kcrit_open(open_pnts)
                            ! kcrit gathered onto the open-point index, for
                            ! passing into xylem_conductance_CW_jls.

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (kg/m2/s)
REAL(KIND=real_jlslsm) :: gamma_root_psi(land_pts)
                            ! Lower incomplete gamma function for the root
                            ! zone water potential. Used in Newton-Raphson
                            ! method.
REAL(KIND=real_jlslsm) :: leaf_k_prev(n_e_leaf, open_pnts)
                            ! leaf_k from the previous Newton-Raphson
                            ! iteration, used for the convergence check.
INTEGER, PARAMETER :: max_nr_iter = 4
                            ! Cap on Newton-Raphson iterations. This used to
                            ! be a fixed 2 iterations with no convergence
                            ! check: accurate near well-watered conditions,
                            ! but leaves a large, one-sided error under water
                            ! stress (leaf_psi ends up not negative enough
                            ! for the demanded e_leaf, which overstates the
                            ! leaf_k reported downstream) - most severely
                            ! once e_leaf approaches the finite maximum this
                            ! vulnerability curve can sustain from
                            ! root_zone_psi, where the true solution runs to
                            ! very negative psi and 2 fixed iterations barely
                            ! move off the initial guess.
                            ! Checked numerically against the exact (gamma-
                            ! function) solution across root/leaf psi from
                            ! -1 to -9.5 MPa: 4 iterations is within ~0.003
                            ! MPa of exact for every case that is actually
                            ! achievable, and the only case that still has
                            ! residual error beyond that (right at this
                            ! curve's finite supply ceiling) is already
                            ! correctly excluded by the kl_sample > kcrit
                            ! mask in stom_opt_mod by iteration 2, so more
                            ! iterations there buy nothing. Originally tried
                            ! 15, which was correct but ~2.3x slower overall
                            ! because the exit check below is array-wide
                            ! (blocks on the whole 800-sample vector, not
                            ! per-sample) - see the exit check's own comment.

INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Select the method to approximate the leaf water potential
SELECT CASE (som_psi_solver)

! ---------------------------------------------------------------------
! Apply a zeroth order taylor expansion to simplify the integration
!
! Note: The taylor series expansion is approximated at the reference
!       water potential (reference_psi). The reference water potential
!       is the mean (midpoint) of the root zone and leaf water
!       potentials.
! ---------------------------------------------------------------------
CASE(psi_solver_taylor)

  DO j = 1,open_pnts ! Iterate over the open points.
      l = veg_index(open_index(j))
      kmax_open(j) = kmax(l)
      kcrit_open(j) = kcrit(l)

      ! Calculate the conductance conversion limit
      k_conversion_limit = 0.1 * kcrit(l)

      ! Set the initial leaf water potential to that of the root zone
      leaf_psi(:,j) = root_zone_psi(l)

      ! Calculate the initial xylem conductance for the initial
      !  reference water potential (this is equal to the root
      !  zone water potential).
      !  k(psi) = kmax * exp(-(psi/b) ** c)
      leaf_k(:,j) = kmax(l)                                               &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

      !  Iterate over the transpiration for each open point.
      DO i = 1,n_e_leaf

        DO k = 1,100 ! Iterate untill the change in xylem conductance is
                     ! less than k_conversion_limit. or 100 iterations.

          ! Calculate the new leaf water potential.
          leaf_psi(i,j) = root_zone_psi(l) - e_leaf(i,j) / leaf_k(i,j)

          ! Calculate the reference psi at which to apply the zeroth order
          ! Taylor series expansion.
          !  psi_ref = mean(leaf_psi, root_zone_psi)
          reference_psi = 0.5 * (leaf_psi(i,j) + root_zone_psi(l))

          ! Calculate the new xylem conductance
          !  k(psi) = kmax * exp(-(psi/b) ** c)
          leaf_k_new = kmax(l)                                              &
             * EXP( -(reference_psi/conductance_b(pft))**conductance_c(pft) )

          ! If the change in xylem conductance is less than
          !  k_conversion_limit then break the loop.
          IF (ABS(leaf_k_new - leaf_k(i,j)) < k_conversion_limit) THEN
            !leaf_k(i,j) = leaf_k_new
            EXIT
          END IF

          ! Set the xylem conductance to the new value for the next
          ! iteration.
          leaf_k(i,j) = leaf_k_new

        END DO ! Iterate untill conductance converges.
      END DO ! Iterate over the sample ci values.
    END DO ! Iterate over the open points.

    ! Calculate conductance at th leaf water potential. The value calculated
    !  above is the conductance at the reference psi.
    CALL xylem_conductance_CW_jls(pft, n_e_leaf, open_pnts, leaf_psi,         &
                                  kmax_open, kcrit_open, leaf_k)

    leaf_k(:,:) = MAX(leaf_k(:,:),                                            &
                      SPREAD(kcrit_open, DIM = 1, NCOPIES = n_e_leaf))

! ---------------------------------------------------------------------
! Apply a Newton-Raphson method to approximate the leaf water potential
! from the transpiration rate. The transpiration rate is given by the
! integral of the conductance model, k(psi), over the range
! [psi_r, psi_l]. Here the conductance model is given by the cumulative
! Weibull distribution. The integral of k(psi) with respect to psi is;
!
!  E(psi_l, psi_r) = kmax (b/c)
!                     * [ gamma(1/c, (psi_r/b)^c)
!                         - gamma(1/c, (psi_l/b)^c) ]
!
! Where gamma is the lower incomplete gamma function.
!
! The leaf water potential is then calculated by solving the equation;
!
!  0 = E(psi_l, psi_r) - e_leaf
!
! This means that the updated prediction of the leaf water potential
! (psi_l') is given by;
!
!  psi_l' = psi_r - (e_leaf - E(psi_l, psi_r))/k(psi_l)
!
! ---------------------------------------------------------------------
CASE(psi_solver_newton)
  ! Root / stem / leaf segments in series (l_som_plant_segments): leaf psi
  ! from solving the segments downstream, and the whole-plant conductance
  ! -dE/dpsi_leaf returned in leaf_k. The single-segment code below is left
  ! unchanged for the default case.
  IF ( l_som_plant_segments ) THEN
    CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,  &
                                open_index, e_leaf, root_zone_psi, kmax, kcrit, &
                                leaf_psi, leaf_k )
  ELSE

  ! Calculate the lower incomplete gamma function for the root zone water
  ! potential, for all land points.
  gamma_root_psi = incomplete_gamma(land_pts, 1/conductance_c(pft),            &
                                    (root_zone_psi/conductance_b(pft))         &
                                    **conductance_c(pft))

  DO j = 1, open_pnts ! Iterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        * EXP( -(root_zone_psi(l)/conductance_b(pft))**conductance_c(pft) )

    ! The initial guess for the leaf water potential is that of the root zone
    ! meaning that the transpiration rate is zero. This simplifies the first
    ! prediction of the leaf water potential to;
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j) / leaf_k(:,j)

    ! Convergence tolerance for this point, matching the psi_solver_taylor
    ! tolerance above (10% of kcrit).
    k_conversion_limit = 0.1 * kcrit(l)

    ! Iterate the water potential until converged, or until max_nr_iter is
    ! reached. A sample whose e_leaf demand exceeds what this vulnerability
    ! curve can ever supply from root_zone_psi (i.e. no finite leaf_psi
    ! satisfies it) will not converge: leaf_k keeps shrinking each pass and
    ! leaf_psi keeps stepping further negative, which is the correct
    ! behaviour here - it drives leaf_k below kcrit so stom_opt_mod's
    ! feasibility mask (kl_sample > kcrit) correctly rejects the sample,
    ! rather than silently reporting the stalled, too-high leaf_k that the
    ! old fixed-2-iteration version left behind.
    DO i = 1, max_nr_iter

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k_prev(:,j) = leaf_k(:,j)
      leaf_k(:,j) = kmax(l)                                                    &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

      ! Calculate the transpiration rate for the current leaf water potential.
      ! First get the incomplete gamma function for the leaf water potential.
      e_leaf_current(:,j) = incomplete_gamma(n_e_leaf,                         &
                                             1/conductance_c(pft),             &
                                             (leaf_psi(:,j)/conductance_b(pft))&
                                              **conductance_c(pft))

      ! Calculate the difference between the incomplete gamma function at the
      ! leaf and root zone water potentials.
      e_leaf_current(:,j) = e_leaf_current(:,j) - gamma_root_psi(l)

      ! Multiply by kmax * (-b/c) to get the transpiration rate
      ! Note: The negative sign is present because conductance_b is negative.
      e_leaf_current(:,j) = e_leaf_current(:,j) * kmax(l)                      &
          * (-conductance_b(pft)/conductance_c(pft))

      ! Calculate the new leaf water potential. leaf_k is floored in the
      ! denominator to avoid a literal divide-by-zero for a sample that has
      ! already collapsed leaf_k to (numerically) 0; MAX(...) here does not
      ! change the update for any sample that is still within reach of a
      ! finite solution.
      leaf_psi(:,j) = leaf_psi(:,j)                                            &
             - (e_leaf(:,j) - e_leaf_current(:,j))                             &
               / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

      ! Guard against runaway steps for an infeasible sample: once leaf_k
      ! has collapsed there is no need to keep tracking the true (very
      ! negative) asymptotic solution, only to land clearly below kcrit
      ! without risking overflow in (leaf_psi/b)**c on the next pass.
      leaf_psi(:,j) = MAX(leaf_psi(:,j),                                       &
                          root_zone_psi(l) - 5.0 * ABS(conductance_b(pft)))

      IF (MAXVAL(ABS(leaf_k(:,j) - leaf_k_prev(:,j))) < k_conversion_limit) EXIT

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax(l)                                                    &
          * EXP( -(leaf_psi(:,j)/conductance_b(pft))**conductance_c(pft) )

  END DO

  END IF  ! l_som_plant_segments

! ---------------------------------------------------------------------
! Invert the tabulated supply function (see lut_s at the top of this
! module): S(psi_l) = S(psi_r) + e_leaf/kmax. Exact up to the table's
! linear interpolation, with no iteration. A demand beyond what the curve
! can supply from psi_r puts the leaf at the bottom of the table, where
! k = lut_f_floor * kmax < kcrit, so stom_opt_mod rejects it as before.
! ---------------------------------------------------------------------
CASE(psi_solver_lut)
  IF ( l_som_plant_segments ) THEN
    CALL leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pnts,         &
                                    veg_index, open_index, e_leaf,              &
                                    root_zone_psi, kmax, leaf_psi, leaf_k )
  ELSE
    CALL leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,    &
                              open_index, e_leaf, root_zone_psi, kmax,          &
                              leaf_psi, leaf_k )
  END IF  ! l_som_plant_segments

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_solver should be psi_solver_taylor, psi_solver_newton ' //       &
        'or psi_solver_lut')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_CW_jls


! *****************************************************************************
! Leaf water potential and conductance by inverting the tabulated supply
! function (see lut_s at the top of this module):
!   S(psi_l) = S(psi_r) + e_leaf/kmax.
! Exact up to the table's linear interpolation, with no iteration. A demand
! beyond what the curve can supply from psi_r puts the leaf at the bottom of
! the table, where k = lut_f_floor * kmax < kcrit, so stom_opt_mod rejects it
! as for the Newton-Raphson solve. Has no automatic arrays, so it is cheap to
! call for a single sample per point.
! *****************************************************************************
SUBROUTINE leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index, &
                                open_index, e_leaf, root_zone_psi, kmax,       &
                                leaf_psi, leaf_k )

USE pftparm, ONLY: conductance_b, conductance_c

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pts, open_pnts, veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point.
, root_zone_psi(land_pts)                                                      &
                            ! Water potential in the root zone (Pa).
, kmax(land_pts)
                            ! Maximum xylem conductance (same units as
                            ! e_leaf per Pa).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential (Pa).
, leaf_k(n_e_leaf, open_pnts)
                            ! Xylem conductance at the leaf water potential.

INTEGER :: i, j, l

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  DO i = 1, n_e_leaf
    leaf_psi(i,j) = supply_lut_psi(pft, root_zone_psi(l),                    &
                    e_leaf(i,j) / MAX(kmax(l), TINY(1.0_real_jlslsm)))
    leaf_k(i,j) = kmax(l) * supply_lut_f(pft, leaf_psi(i,j))
  END DO
END DO

END SUBROUTINE leaf_psi_lut_jls


! *****************************************************************************
! Leaf water potential (Pa) from the supply-function table, for root zone
! potential psi_root (Pa) and normalised demand e_over_kmax = e_leaf/kmax (Pa):
! solves S(psi_l) = S(psi_r) + e_over_kmax. See leaf_psi_lut_jls.
! *****************************************************************************
FUNCTION supply_lut_psi( pft, psi_root, e_over_kmax ) RESULT( psi_l )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, e_over_kmax
REAL(KIND=real_jlslsm) :: psi_l

INTEGER :: lo, hi, mid, i_root
REAL(KIND=real_jlslsm) :: dpsi, r_root, s_target

CALL build_supply_lut(pft)
dpsi = lut_dpsi(pft)

! S at the root zone water potential plus the demand, and the table cell
! psi_root is in (the search starts there).
r_root = MIN(MAX(-psi_root / dpsi, 0.0), REAL(n_lut - 1))
i_root = MIN(INT(r_root) + 1, n_lut - 1)
s_target = lut_s_at(pft, psi_root) + e_over_kmax

IF ( s_target >= lut_s(n_lut,pft) ) THEN
  psi_l = MIN(-(n_lut - 1) * dpsi, psi_root)
ELSE
  ! lut_s(lo) <= s_target < lut_s(hi); S is strictly increasing.
  lo = i_root
  hi = n_lut
  DO WHILE ( hi - lo > 1 )
    mid = (lo + hi) / 2
    IF ( lut_s(mid,pft) <= s_target ) THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO
  psi_l = -(lo - 1) * dpsi - dpsi * (s_target - lut_s(lo,pft))                 &
          / (lut_s(hi,pft) - lut_s(lo,pft))
  psi_l = MIN(psi_l, psi_root)
END IF

END FUNCTION supply_lut_psi


! *****************************************************************************
! Largest normalised transpiration, e/kmax (Pa), the curve can supply from
! psi_root while keeping k >= kcrit_frac * kmax: S(psi_crit) - S(psi_root),
! with k(psi_crit) = kcrit_frac * kmax. Zero if psi_root is already beyond
! psi_crit.
! *****************************************************************************
FUNCTION supply_lut_e_crit( pft, psi_root, kcrit_frac ) RESULT( e_crit )

USE pftparm, ONLY: conductance_b, conductance_c

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, kcrit_frac
REAL(KIND=real_jlslsm) :: e_crit

REAL(KIND=real_jlslsm) :: psi_crit, r

CALL build_supply_lut(pft)
! f(psi_crit) = kcrit_frac, kept inside (0, 1) (e.g. kcrit >= kmax would
! otherwise give a NaN).
r = MIN(MAX(kcrit_frac, 1.0e-6_real_jlslsm), 1.0_real_jlslsm - 1.0e-6_real_jlslsm)
IF ( lut_sox(pft) ) THEN
  psi_crit = conductance_b(pft)                                                &
             * (1.0 / r - 1.0)**(1.0 / conductance_c(pft))
ELSE
  psi_crit = conductance_b(pft)                                                &
             * (-LOG(r))**(1.0 / conductance_c(pft))
END IF
e_crit = MAX(lut_s_at(pft, psi_crit) - lut_s_at(pft, psi_root), 0.0)

END FUNCTION supply_lut_e_crit

! Relative conductance f(psi) = k/kmax of the PFT's conductance model.
FUNCTION supply_lut_f( pft, psi ) RESULT( f )

USE pftparm, ONLY: conductance_b, conductance_c

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: f

CALL build_supply_lut(pft)
IF ( lut_sox(pft) ) THEN
  f = 1.0 / (1.0 + (ABS(psi / conductance_b(pft)))**conductance_c(pft))
ELSE
  f = EXP(-(ABS(psi / conductance_b(pft)))**conductance_c(pft))
END IF

END FUNCTION supply_lut_f

! S(psi) by linear interpolation in the table (the table must be built).
FUNCTION lut_s_at( pft, psi ) RESULT( s_psi )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: s_psi

INTEGER :: i0
REAL(KIND=real_jlslsm) :: r

r = MIN(MAX(-psi / lut_dpsi(pft), 0.0), REAL(n_lut - 1))
i0 = MIN(INT(r) + 1, n_lut - 1)
s_psi = lut_s(i0,pft) + (r - (i0 - 1)) * (lut_s(i0+1,pft) - lut_s(i0,pft))

END FUNCTION lut_s_at

! *****************************************************************************
! Build (once per PFT) the supply-function table used by psi_solver_lut:
! lut_s(i,pft) = integral(f, psi_i, 0) by Simpson's rule on each grid cell,
! for the PFT's conductance model (see lut_s above). Built by numerical
! integration rather than the closed forms (incomplete gamma /
! hypergeometric series), whose series are poorly conditioned deep into the
! curve.
! *****************************************************************************
SUBROUTINE build_supply_lut( pft )

USE max_dimensions, ONLY: npft_max
USE pftparm, ONLY: conductance_b, conductance_c, pft_conductance_model,       &
                   kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: SOX_conductance

INTEGER, INTENT(IN) :: pft

INTEGER :: i
REAL(KIND=real_jlslsm) :: b, c, dpsi, p0, pm, p1, f_floor
LOGICAL :: l_sox

IF ( ALLOCATED(lut_ready) ) THEN
  IF ( lut_ready(pft) ) RETURN
END IF

!$OMP CRITICAL (som_supply_lut)
IF ( .NOT. ALLOCATED(lut_ready) ) THEN
  ALLOCATE( lut_s(n_lut, npft_max), lut_dpsi(npft_max), lut_ready(npft_max), &
            lut_sox(npft_max) )
  lut_ready(:) = .FALSE.
END IF

IF ( .NOT. lut_ready(pft) ) THEN
  b = conductance_b(pft)
  c = conductance_c(pft)
  l_sox = pft_conductance_model(pft) == SOX_conductance
  ! Bottom of the table: where f = f_floor.
  IF ( l_sox ) THEN
    f_floor = lut_f_floor_kcrit * (1.0 - kcrit_fractional_loss(pft))
    dpsi = ABS(b) * (1.0 / f_floor - 1.0)**(1.0 / c) / (n_lut - 1)
  ELSE
    dpsi = ABS(b) * (-LOG(lut_f_floor))**(1.0 / c) / (n_lut - 1)
  END IF
  lut_dpsi(pft) = dpsi
  lut_s(1,pft) = 0.0
  DO i = 2, n_lut
    p0 = -(i - 2) * dpsi
    p1 = -(i - 1) * dpsi
    pm = 0.5 * (p0 + p1)
    lut_s(i,pft) = lut_s(i-1,pft) + dpsi / 6.0                                 &
                   * ( f_of(p0) + 4.0 * f_of(pm) + f_of(p1) )
  END DO
  lut_sox(pft) = l_sox
  lut_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut)

CONTAINS

  REAL(KIND=real_jlslsm) FUNCTION f_of(p)
  REAL(KIND=real_jlslsm), INTENT(IN) :: p
  IF ( l_sox ) THEN
    f_of = 1.0 / (1.0 + (ABS(p / b))**c)
  ELSE
    f_of = EXP(-(ABS(p / b))**c)
  END IF
  END FUNCTION f_of

END SUBROUTINE build_supply_lut

! *****************************************************************************
! Leaf water potential and whole-plant conductance for the plant as root,
! stem and leaf segments in series (l_som_plant_segments), as in GDAY
! gs_opt.c (setup_plant / plant_supply) and Wang et al. (2019).
!
! Segment s has maximum conductance kmax * seg_kfac(pft,s) (the segments
! share the whole-plant resistance, so in series they give kmax when well
! watered) and its own cumulative Weibull curve (conductance_b_seg,
! conductance_c_seg). For each sampled transpiration the outlet potential of
! each segment is solved downstream from psi_root_zone with the same
! Newton-Raphson as the single segment (psi_solver_newton), and becomes the next
! segment's inlet. The last outlet is the leaf water potential.
!
! The conductance the hydraulic cost uses is the whole-plant k = -dE/dpsi_leaf
! (Sperry et al. 2017), from E = int_{psi_out}^{psi_in} k_s dpsi:
!   dpsi_out/dE = (k_s(psi_in) * dpsi_in/dE - 1) / k_s(psi_out),
! starting from dpsi_root_zone/dE = 0. With no flow 1/k = sum_s 1/k_s(psi);
! a single segment gives k = k(psi_leaf) (the default model).
! *****************************************************************************
SUBROUTINE leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pnts,          &
                                  veg_index, open_index, e_leaf,               &
                                  root_zone_psi, kmax, kcrit,                  &
                                  leaf_psi, leaf_k )

USE pftparm, ONLY: seg_kfac, conductance_b_seg, conductance_c_seg

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts),                                                 &
      ! Sampled transpiration (same units as kmax * psi).
  root_zone_psi(land_pts),                                                     &
      ! Root-zone water potential (Pa).
  kmax(land_pts),                                                              &
      ! Whole-plant maximum conductance.
  kcrit(land_pts)
      ! Whole-plant critical conductance.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts),                                               &
      ! Leaf water potential (Pa).
  leaf_k(n_e_leaf, open_pnts)
      ! Whole-plant conductance -dE/dpsi_leaf.

INTEGER, PARAMETER :: n_seg = 3, max_nr_iter = 4

REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12
      ! Floor on a segment's conductance, as a fraction of its kmax.

INTEGER :: j, l, iseg, it

REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(n_e_leaf), psi_out(n_e_leaf),                                         &
      ! Segment inlet / outlet water potential (Pa).
  k_in(n_e_leaf), k_cur(n_e_leaf), k_prev(n_e_leaf),                           &
      ! Segment conductance at the inlet / current outlet / previous pass.
  g_in(n_e_leaf), e_cur(n_e_leaf),                                             &
      ! Incomplete gamma at the inlet; flow for the current outlet.
  dpsi_de(n_e_leaf),                                                           &
      ! d(psi_out)/dE carried down the segments.
  g_one(1),                                                                    &
      ! Incomplete gamma at psi_root_zone (root inlet, same for all samples).
  kmx, kcr, bs, cs
      ! Segment kmax, kcrit, Weibull b (Pa) and c.

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  psi_in(:)  = root_zone_psi(l)
  dpsi_de(:) = 0.0

  DO iseg = 1, n_seg
    kmx = kmax(l) * seg_kfac(pft,iseg)
    kcr = kcrit(l) * seg_kfac(pft,iseg)
    bs  = conductance_b_seg(pft,iseg)
    cs  = conductance_c_seg(pft,iseg)

    k_in(:) = kmx * EXP( -(psi_in(:)/bs)**cs )
    IF ( iseg == 1 ) THEN
      ! The root inlet is psi_root_zone for every sample: one evaluation.
      g_one(:) = incomplete_gamma(1, 1.0/cs, [(psi_in(1)/bs)**cs])
      g_in(:)  = g_one(1)
    ELSE
      g_in(:) = incomplete_gamma(n_e_leaf, 1.0/cs, (psi_in(:)/bs)**cs)
    END IF

    ! First guess: the inlet conductance over the whole drop.
    psi_out(:) = psi_in(:) - e_leaf(:,j) / MAX(k_in(:), TINY(1.0_real_jlslsm))
    psi_out(:) = MAX(psi_out(:), psi_in(:) - 5.0 * ABS(bs))
    k_cur(:)   = k_in(:)

    DO it = 1, max_nr_iter
      k_prev(:) = k_cur(:)
      k_cur(:)  = kmx * EXP( -(psi_out(:)/bs)**cs )
      e_cur(:)  = ( incomplete_gamma(n_e_leaf, 1.0/cs, (psi_out(:)/bs)**cs)    &
                    - g_in(:) ) * kmx * (-bs/cs)
      psi_out(:) = psi_out(:) - (e_leaf(:,j) - e_cur(:))                       &
                                / MAX(k_cur(:), TINY(1.0_real_jlslsm))
      ! Same runaway guard as the single segment (infeasible demand).
      psi_out(:) = MAX(psi_out(:), psi_in(:) - 5.0 * ABS(bs))
      IF ( MAXVAL(ABS(k_cur(:) - k_prev(:))) < 0.1 * kcr ) EXIT
    END DO

    ! Segment conductances are floored at a tiny fraction of the segment
    ! kmax (not TINY): a collapsed segment has k -> 0 by underflow, and with
    ! a TINY floor the recursion overflows to -Inf and then 0 * -Inf = NaN
    ! in the next segment, which previously came out as a huge whole-plant k
    ! and let an infeasible state (psi_leaf ~ -70 MPa) pass the k > kcrit
    ! mask. With this floor dpsi_de stays finite and k correctly ~ 0.
    k_cur(:)   = MAX(kmx * EXP( -(psi_out(:)/bs)**cs ), k_floor * kmx)
    dpsi_de(:) = ( MAX(k_in(:), k_floor * kmx) * dpsi_de(:) - 1.0 ) / k_cur(:)
    psi_in(:)  = psi_out(:)
  END DO

  leaf_psi(:,j) = psi_in(:)
  ! dpsi_de < 0 and finite here; anything else (defensive) is infeasible.
  WHERE ( dpsi_de(:) < 0.0 .AND. ABS(dpsi_de(:)) < HUGE(1.0_real_jlslsm) )
    leaf_k(:,j) = -1.0 / dpsi_de(:)
  ELSEWHERE
    leaf_k(:,j) = 0.0
  END WHERE
END DO

END SUBROUTINE leaf_psi_segments_jls

! *****************************************************************************
! As leaf_psi_segments_jls, with each segment solved from its supply-function
! table instead of Newton-Raphson (psi_solver_lut): segment s carries
!   E = kmax_s * ( S_s(psi_out) - S_s(psi_in) ),  kmax_s = kmax * seg_kfac,
! so psi_out = S_s^-1( S_s(psi_in) + E / kmax_s ), downstream from
! psi_root_zone. The whole-plant conductance -dE/dpsi_leaf comes from the
! same recursion, dpsi_out/dE = (k_s(psi_in) dpsi_in/dE - 1) / k_s(psi_out),
! with the same floor on segment conductances. A demand beyond what a
! segment can supply puts its outlet at the bottom of its table, where k is
! ~0, so the state is rejected as infeasible (k <= kcrit).
! *****************************************************************************
SUBROUTINE leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pnts,      &
                                      veg_index, open_index, e_leaf,           &
                                      root_zone_psi, kmax, leaf_psi, leaf_k )

USE pftparm, ONLY: seg_kfac, conductance_b_seg, conductance_c_seg

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts), root_zone_psi(land_pts), kmax(land_pts)
      ! As leaf_psi_segments_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts), leaf_k(n_e_leaf, open_pnts)
      ! Leaf water potential (Pa) and whole-plant -dE/dpsi_leaf.

INTEGER, PARAMETER :: n_seg = 3

REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12
      ! Floor on a segment's conductance, as a fraction of its kmax (as
      ! leaf_psi_segments_jls).

INTEGER :: i, j, l, iseg
REAL(KIND=real_jlslsm) :: kmx, bs, cs, psi_in, psi_out, k_in, k_out, dpsi_de

CALL build_supply_lut_seg(pft)

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  DO i = 1, n_e_leaf
    psi_in  = root_zone_psi(l)
    dpsi_de = 0.0
    DO iseg = 1, n_seg
      kmx = kmax(l) * seg_kfac(pft,iseg)
      bs  = conductance_b_seg(pft,iseg)
      cs  = conductance_c_seg(pft,iseg)
      psi_out = seg_lut_psi(pft, iseg, psi_in,                                 &
                            e_leaf(i,j) / MAX(kmx, TINY(1.0_real_jlslsm)))
      k_in  = MAX(kmx * EXP( -(ABS(psi_in /bs))**cs ), k_floor * kmx)
      k_out = MAX(kmx * EXP( -(ABS(psi_out/bs))**cs ), k_floor * kmx)
      dpsi_de = ( k_in * dpsi_de - 1.0 ) / k_out
      psi_in  = psi_out
    END DO
    leaf_psi(i,j) = psi_in
    IF ( dpsi_de < 0.0 .AND. ABS(dpsi_de) < HUGE(1.0_real_jlslsm) ) THEN
      leaf_k(i,j) = -1.0 / dpsi_de
    ELSE
      leaf_k(i,j) = 0.0
    END IF
  END DO
END DO

END SUBROUTINE leaf_psi_segments_lut_jls

! *****************************************************************************
! Outlet water potential (Pa) of segment iseg for inlet psi_in (Pa) and
! normalised flow e_over_k = E / kmax_s (Pa): S_s(psi_out) = S_s(psi_in) +
! e_over_k, by binary search and linear interpolation (as supply_lut_psi).
! *****************************************************************************
FUNCTION seg_lut_psi( pft, iseg, psi_in, e_over_k ) RESULT( psi_out )

INTEGER, INTENT(IN) :: pft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

INTEGER :: lo, hi, mid, i0
REAL(KIND=real_jlslsm) :: dpsi, r, s_in, s_target

dpsi = lut_seg_dpsi(pft,iseg)
r  = MIN(MAX(-psi_in / dpsi, 0.0), REAL(n_lut - 1))
i0 = MIN(INT(r) + 1, n_lut - 1)
s_in = lut_seg_s(i0,pft,iseg)                                                  &
       + (r - (i0 - 1)) * (lut_seg_s(i0+1,pft,iseg) - lut_seg_s(i0,pft,iseg))
s_target = s_in + MAX(e_over_k, 0.0)

IF ( s_target >= lut_seg_s(n_lut,pft,iseg) ) THEN
  psi_out = MIN(-(n_lut - 1) * dpsi, psi_in)
ELSE
  lo = i0
  hi = n_lut
  DO WHILE ( hi - lo > 1 )
    mid = (lo + hi) / 2
    IF ( lut_seg_s(mid,pft,iseg) <= s_target ) THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO
  psi_out = -(lo - 1) * dpsi - dpsi * (s_target - lut_seg_s(lo,pft,iseg))      &
            / (lut_seg_s(hi,pft,iseg) - lut_seg_s(lo,pft,iseg))
  psi_out = MIN(psi_out, psi_in)
END IF

END FUNCTION seg_lut_psi

! *****************************************************************************
! Build (once per PFT) the segment supply-function tables: as
! build_supply_lut for each segment's Weibull, from psi = 0 down to where the
! normalised conductance falls to lut_f_floor.
! *****************************************************************************
SUBROUTINE build_supply_lut_seg( pft )

USE max_dimensions, ONLY: npft_max
USE pftparm, ONLY: conductance_b_seg, conductance_c_seg

INTEGER, INTENT(IN) :: pft

INTEGER :: i, iseg
REAL(KIND=real_jlslsm) :: b, c, dpsi, p0, pm, p1

IF ( ALLOCATED(lut_seg_ready) ) THEN
  IF ( lut_seg_ready(pft) ) RETURN
END IF

!$OMP CRITICAL (som_supply_lut_seg)
IF ( .NOT. ALLOCATED(lut_seg_ready) ) THEN
  ALLOCATE( lut_seg_s(n_lut, npft_max, 3), lut_seg_dpsi(npft_max, 3),         &
            lut_seg_ready(npft_max) )
  lut_seg_ready(:) = .FALSE.
END IF

IF ( .NOT. lut_seg_ready(pft) ) THEN
  DO iseg = 1, 3
    b = conductance_b_seg(pft,iseg)
    c = conductance_c_seg(pft,iseg)
    dpsi = ABS(b) * (-LOG(lut_f_floor))**(1.0 / c) / (n_lut - 1)
    lut_seg_dpsi(pft,iseg) = dpsi
    lut_seg_s(1,pft,iseg) = 0.0
    DO i = 2, n_lut
      p0 = -(i - 2) * dpsi
      p1 = -(i - 1) * dpsi
      pm = 0.5 * (p0 + p1)
      lut_seg_s(i,pft,iseg) = lut_seg_s(i-1,pft,iseg) + dpsi / 6.0             &
                              * ( EXP(-(ABS(p0 / b))**c)                       &
                                  + 4.0 * EXP(-(ABS(pm / b))**c)               &
                                  + EXP(-(ABS(p1 / b))**c) )
    END DO
  END DO
  lut_seg_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut_seg)

END SUBROUTINE build_supply_lut_seg

! ---------------------------------------------------------------------
! Function to calculate the lower incomplete gamma function. The
! lower incomplete gamma function is defined as;
!
! gamma(a, x) = integral( t^(a-1) * e^(-t) dt, 0, x)
!
! The code here is based on the implementation in Numerical Recipes
! in FORTRAN 90, page 211, and makes use of the series representation
! of the lower incomplete gamma function.
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma( n_x, a, x ) RESULT( gamma )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_x ! Size of x array

REAL(KIND=real_jlslsm), INTENT(IN) :: a
REAL(KIND=real_jlslsm), INTENT(IN) :: x(n_x)

REAL(KIND=real_jlslsm) :: gamma(n_x)

! Local variables
REAL(KIND=real_jlslsm) :: step(n_x)
INTEGER :: i

! Series convergence controls. NOTE: this series solution is only well
! conditioned for x < a+1 or so (see Numerical Recipes in Fortran 90,
! section 6.2): for larger x - i.e. deep into the vulnerability curve,
! close to hydraulic failure - it converges slowly, and the intermediate
! terms can grow very large before being cancelled by the exp(-x) factor
! below, risking overflow. If this routine is ever exercised with
! som_psi_solver = psi_solver_newton under severe water stress and that
! becomes a problem in practice, this branch should be paired with a
! continued-fraction evaluation of the upper incomplete gamma function
! for x >= a+1 (NR's gcf), rather than raising max_terms further.
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 40
                            ! Checked numerically against this branch's own
                            ! kcrit_fractional_loss=0.95 cutoff: with the
                            ! corrected P50/P88 constants, k/kmax crosses
                            ! below kcrit at x corresponding to psi ~ -8.2
                            ! MPa, which needs ~17-21 series terms here. 40
                            ! terms covers everything up to psi ~ -10.5 MPa
                            ! with real margin past that cutoff - anything
                            ! needing more terms than that is already,
                            ! separately, excluded by the kl_sample > kcrit
                            ! mask in stom_opt_mod regardless of how
                            ! precisely this series converges for it, so
                            ! spending up to 100 terms chasing it (as the
                            ! psi_solver_newton loop, capped at 4 iterations,
                            ! calls this every pass) was wasted work. This
                            ! does not touch the psi_solver_taylor branch or its
                            ! own convergence tolerance.

! zhook variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INCOMPLETE_GAMMA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial values (the n=0 term of the series)
step(:) = 1/a
gamma(:) = step

! Calculate series, iterating until every point's remaining term is
! negligible relative to its running sum, rather than a fixed number of
! terms (previously hard-coded to 2 extra terms, which is only accurate
! for x close to 0 - see max_terms/series_tol comment above for the
! remaining caveat at large x).
DO i = 1, max_terms
  step = step * x / (a + i)
  gamma = gamma + step

  IF (ALL(ABS(step) < series_tol * ABS(gamma))) EXIT
END DO

! Multiply by terms outside the sum
gamma = gamma * EXP(-x) * x**a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma

END MODULE xylem_hydraulics_CW_jls_mod