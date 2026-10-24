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

! The cumulative Weibull vulnerability curve: conductance, the leaf water
! potential solvers on the single curve (Taylor and Newton; the table
! solver is called from here), the curve fractions xylem_f and xylem_f_bc,
! and the incomplete gamma function. The supply tables are in
! supply_table_mod and the segments, the soil-to-root link and the
! zero-flow conductance of the path in plant_hydraulic_path_mod (moved
! there 23/10/26).
!
! MDK, 23/10/26
PRIVATE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='XYLEM_HYDRAULICS_CW_JLS_MOD'

PUBLIC :: leaf_conductance_CW_jls                                              &
,         xylem_conductance_CW_jls                                             &
,         psi_leaf_CW_jls                                                      &
,         incomplete_gamma                                                     &
,         incomplete_gamma_different_a                                         &
,         xylem_f                                                              &
,         xylem_f_bc

CONTAINS

! *****************************************************************************
! Fraction of kmax left on the PFT's whole-plant vulnerability curve at psi
! (Pa), as pft_conductance_model. psi >= 0 gives 1. Xylem only (no soil).
! *****************************************************************************
FUNCTION xylem_f( pft, psi ) RESULT( f )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: f

f = xylem_f_bc( pft, psi, conductance_b_pft(pft), conductance_c_pft(pft) )

END FUNCTION xylem_f

! *****************************************************************************
! As xylem_f, on the curve with parameters b (Pa, negative) and c
! (dimensionless), e.g. a point's impaired curve: cumulative Weibull
! exp(-(psi/b)^c), or SOX 1 / (1 + (psi/b)^c), as pft_conductance_model.
! psi >= 0 gives 1.
!
! MDK, 21/10/26
! *****************************************************************************
FUNCTION xylem_f_bc( pft, psi, b, c ) RESULT( f )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: SOX_conductance

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi, b, c
REAL(KIND=real_jlslsm) :: f, x

x = ( ABS(MIN(psi, 0.0) / b) )**c
IF ( pft_conductance_model(pft) == SOX_conductance ) THEN
  f = 1.0 / (1.0 + x)
ELSE
  f = EXP(-x)
END IF

END FUNCTION xylem_f_bc

! *********************************************************************
! Contains routines used to calculate conductance and leaf water
! potential for a cumulative Weibull conductance model.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! *********************************************************************

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf conductance from the leaf water
! potential for all land points, floored at kcrit: it returns
! MAX(k(psi), kcrit), not k(psi). JBaguley
!
! MDK, 21/10/26
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_CW_jls( pft,                                       &
                                    land_pnts,                                 &
                                    water_potential,                           &
                                    kmax,                                      &
                                    kcrit,                                     &
                                    conductance_b,                             &
                                    conductance_c,                             &
                                 ! INTENT OUT
                                    leaf_conductance                           &
  )

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
  leaf_conductance(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the leaf conductance
leaf_conductance = kmax                                                        &
     * EXP(-ABS(water_potential / conductance_b)                               &
              ** conductance_c)

! Apply minimum conductance limit
leaf_conductance = MAX(leaf_conductance, kcrit)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_CW_jls

! ---------------------------------------------------------------------
! Subroutine to calculate xylem conductance using a cumulative Weibull
! distribution, floored at kcrit: it returns MAX(k(psi), kcrit), not
! k(psi).
! NOTE: This subroutine is designed for the stomatal optimisation model.
!       It assumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
!
! MDK, 21/10/26
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_CW_jls( pft,                                      &
                                     n_water_potentials,                       &
                                     open_pnts,                                &
                                     water_potential,                          &
                                     kmax,                                     &
                                     kcrit,                                    &
                                     conductance_b,                            &
                                     conductance_c,                            &
                                  ! INTENT OUT
                                     xylem_conductance                         &
  )

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
                            ! (mol m-2 s-1 Pa-1): the kmax of this call (whole
                            ! canopy, the sunlit or shaded leaves, or a
                            ! layer's leaves), not always kmax_pft.
, kcrit(open_pnts)                                                             &
                            ! Critical xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1), on the same basis as kmax.
, conductance_b(open_pnts)                                                     &
                            ! Conductance parameter b for each open point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(open_pnts)
                            ! Conductance parameter c for each open point
                            ! (dimensionless).
                            !
                            ! MDK, 21/10/26

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

INTEGER :: i, j
                            ! Water potential and open point.
                            !
                            ! MDK, 23/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! k(psi) = kmax exp(-|psi/b|^c), floored at kcrit.
!
! MDK, 23/10/26
DO j = 1, open_pnts
  DO i = 1, n_water_potentials
    xylem_conductance(i,j) = kmax(j)                                           &
         * EXP(-ABS(water_potential(i,j) / conductance_b(j))                   &
               ** conductance_c(j))
    xylem_conductance(i,j) = MAX(xylem_conductance(i,j), kcrit(j))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_CW_jls

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf water potential from the transpiration
! rate and root zone water potential. Uses a cumulative Weibull
! distribution to model conductance.
!
!                  k(psi) = kmax * e^(-(psi / b) ^ c)
!
! The code below uses one of three methods to estimate the leaf water
! potential from the transpiration rate (som_psi_solver):
!   1: Applies a zeroth order Taylor expansion to the conductance
!      equation (psi_solver_taylor).
!   2: Uses a Newton-Raphson solve of the integrated conductance
!      equation (psi_solver_newton).
!   3: Inverts the tabulated supply function (psi_solver_lut; see
!      psi_leaf_lut_jls).
! With l_som_plant_segments, methods 2 and 3 are not done here: plant_path
! (xylem_hydraulics_jls_mod) calls the segment solvers
! (psi_leaf_segments_jls, psi_leaf_segments_lut_jls) instead; method 1
! with segments is rejected by the namelist check.
!
! The conductance returned, leaf_k, is the curve at the leaf, k(psi_leaf).
! That is also the plant's marginal conductance, the one the hydraulic
! cost needs: E(psi_leaf) = integral( k(psi) dpsi, psi_leaf, psi_root ), so
!
!   -dE / dpsi_leaf = k(psi_leaf)
!
! for one conductance curve. With segments or a soil link in series it is
! not, and those solvers carry dpsi_leaf/dE down the path instead.
!
! NOTE: This subroutine is designed for the stomatal optimisation model.
!       It returns an output array that only contains values for points
!       with open stomata.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_CW_jls( pft,                                               &
                            n_e_leaf,                                          &
                            land_pts,                                          &
                            open_pnts,                                         &
                            veg_index,                                         &
                            open_index,                                        &
                            e_leaf,                                            &
                            root_zone_psi,                                     &
                            kmax,                                              &
                            kcrit,                                             &
                            conductance_b,                                     &
                            conductance_c,                                     &
                         ! INTENT OUT
                            leaf_psi,                                          &
                            leaf_k                                             &
  )

USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_taylor,             &
                                psi_solver_newton, psi_solver_lut
USE supply_table_mod, ONLY: psi_leaf_lut_jls, lut_pft_curve

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
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 s-1)
, root_zone_psi(land_pts)                                                      &
                            ! Water potential in the root zone (Pa)
, kmax(land_pts)                                                               &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1): the kmax of this call
                            ! (whole canopy, the sunlit or shaded leaves, or
                            ! a layer's leaves).
, kcrit(land_pts)                                                              &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1), on the same basis as kmax.
, conductance_b(land_pts)                                                      &
                            ! Conductance parameter b for each land point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(land_pts)
                            ! Conductance parameter c for each land point
                            ! (dimensionless).
                            !
                            ! MDK, 21/10/26

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)
                            !
                            ! MDK, 21/10/26

! Local variables
INTEGER :: i, j, k, l

! Taylor expansion variables
REAL(KIND=real_jlslsm) :: leaf_k_new
REAL(KIND=real_jlslsm) :: k_conversion_limit
LOGICAL :: l_nr_active(n_e_leaf)
                   ! Newton: samples still iterating (each sample stops on
                   ! its own convergence, so the answer for one E does not
                   ! depend on what else is in the call)
                   !
                   ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: reference_psi
REAL(KIND=real_jlslsm) :: kmax_open(open_pnts)
                            ! kmax gathered onto the open-point index, for
                            ! passing into xylem_conductance_CW_jls.
REAL(KIND=real_jlslsm) :: kcrit_open(open_pnts)
                            ! kcrit gathered onto the open-point index, for
                            ! passing into xylem_conductance_CW_jls.
REAL(KIND=real_jlslsm) :: b_open(open_pnts), c_open(open_pnts)
                            ! conductance_b/c gathered onto the open-point
                            ! index, for passing into xylem_conductance_CW_jls.

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (mol H2O m-2 s-1)
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm) :: gamma_root_psi(land_pts)
                            ! Lower incomplete gamma function for the root
                            ! zone water potential. Used in Newton-Raphson
                            ! method.
REAL(KIND=real_jlslsm) :: leaf_k_prev(n_e_leaf, open_pnts)
                            ! leaf_k from the previous Newton-Raphson
                            ! iteration, used for the convergence check.
INTEGER, PARAMETER :: max_nr_iter = 4
                            ! Cap on Newton-Raphson iterations. Checked
                            ! against the exact (gamma-function) solution for
                            ! root and leaf psi from -1 to -9.5 MPa: 4 is
                            ! within ~0.003 MPa for every reachable case.
                            ! Samples that do not converge (demand beyond what
                            ! the curve can supply) end with k below kcrit
                            ! and are rejected by the feasibility mask.
                            !
                            ! MDK, 21/10/26

INTEGER :: errcode

INTEGER :: psi_solver
                            ! som_psi_solver, or Newton where the supply
                            ! table cannot follow a point's own b and c.
                            !
                            ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The table is keyed by PFT, so a point with its own curve (b, c) is solved
! by Newton. (Segments with the Newton or table solver do not come here:
! plant_path in xylem_hydraulics_jls_mod calls the segment solvers.)
!
! MDK, 23/10/26
psi_solver = som_psi_solver
IF ( psi_solver == psi_solver_lut ) THEN
  IF ( .NOT. lut_pft_curve( pft, land_pts, open_pnts, veg_index, open_index,   &
                            conductance_b, conductance_c ) )                   &
    psi_solver = psi_solver_newton
END IF

! Select the method to approximate the leaf water potential
SELECT CASE (psi_solver)

! ---------------------------------------------------------------------
! Apply a zeroth order Taylor expansion to simplify the integration
!
! Note: The Taylor series expansion is approximated at the reference
!       water potential (reference_psi). The reference water potential
!       is the mean (midpoint) of the root zone and leaf water
!       potentials.
! psi is held at <= 0 in (psi/b)**c: a positive psi gave a negative base
! and NaN. Unchanged for psi <= 0.
!
! In plain terms: "zeroth order" means the conductance is treated as
! constant over the drop from the root zone to the leaf, taken at the
! midpoint (the midpoint rule for the supply integral). The leaf water
! potential then follows Ohm's law,
!
!   psi_leaf = psi_root - E / k(psi_mid),   psi_mid = (psi_leaf + psi_root) / 2
!
! which is repeated (psi_leaf -> psi_mid -> k) until k changes by less than
! 0.1 kcrit. It is cheap, but on a steep curve it overestimates k (and so
! underestimates the drop). The k returned is then re-evaluated at
! psi_leaf and floored at kcrit.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
CASE(psi_solver_taylor)

  DO j = 1,open_pnts ! Iterate over the open points.
    l = veg_index(open_index(j))
    kmax_open(j) = kmax(l)
    kcrit_open(j) = kcrit(l)
    b_open(j) = conductance_b(l)
    c_open(j) = conductance_c(l)

    ! Calculate the conductance conversion limit
    k_conversion_limit = 0.1 * kcrit(l)

    ! Set the initial leaf water potential to that of the root zone
    leaf_psi(:,j) = root_zone_psi(l)

    ! Calculate the initial xylem conductance for the initial
    !  reference water potential (this is equal to the root
    !  zone water potential).
    !  k(psi) = kmax * exp(-(psi/b) ** c)
    leaf_k(:,j) = kmax(l)                                                      &
        * EXP( -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))          &
               **conductance_c(l) )

    !  Iterate over the transpiration for each open point.
    DO i = 1,n_e_leaf

      DO k = 1,100 ! Iterate until the change in xylem conductance is
                   ! less than k_conversion_limit, or 100 iterations.

        ! Calculate the new leaf water potential. leaf_k is floored at TINY,
        ! as in the Newton solver: deep in the curve kmax exp(-(psi/b)^c)
        ! underflows to 0, which gave -Inf (E > 0) or NaN (E = 0).
        !
        ! MDK, 23/10/26
        leaf_psi(i,j) = root_zone_psi(l)                                       &
                        - e_leaf(i,j) / MAX(leaf_k(i,j), TINY(1.0_real_jlslsm))

        ! Calculate the reference psi at which to apply the zeroth order
        ! Taylor series expansion.
        !  psi_ref = mean(leaf_psi, root_zone_psi)
        reference_psi = 0.5 * (leaf_psi(i,j) + root_zone_psi(l))

        ! Calculate the new xylem conductance
        !  k(psi) = kmax * exp(-(psi/b) ** c)
        leaf_k_new = kmax(l)                                                   &
           * EXP( -(MIN(reference_psi,0.0_real_jlslsm)/conductance_b(l))       &
                  **conductance_c(l) )

        ! If the change in xylem conductance is less than
        !  k_conversion_limit then break the loop.
        IF (ABS(leaf_k_new - leaf_k(i,j)) < k_conversion_limit) THEN
          EXIT
        END IF

        ! Set the xylem conductance to the new value for the next
        ! iteration.
        leaf_k(i,j) = leaf_k_new

      END DO ! Iterate until conductance converges.
    END DO ! Iterate over the sample transpiration values.
  END DO ! Iterate over the open points.

  ! Calculate conductance at the leaf water potential. The value calculated
  !  above is the conductance at the reference psi.
  CALL xylem_conductance_CW_jls(pft, n_e_leaf, open_pnts, leaf_psi,            &
                                kmax_open, kcrit_open, b_open, c_open,         &
                                leaf_k)

  leaf_k(:,:) = MAX(leaf_k(:,:),                                               &
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
! Newton-Raphson on F(psi_l) = E(psi_l, psi_r) - e_leaf. Since E is the
! integral of k down to psi_l, dF/dpsi_l = -k(psi_l), so each step moves
! the current estimate psi_l by the missing flow over the conductance there:
!
!  psi_l' = psi_l - (e_leaf - E(psi_l, psi_r)) / k(psi_l)
!
! (corrected 23/10/26: this read psi_r - ..., but the code has always
! stepped from psi_l). The first guess is psi_l = psi_r - e_leaf / k(psi_r);
! each step is kept no more than 5 |b| below psi_r, and the loop stops when
! k changes by less than k_conversion_limit (or after max_nr_iter steps).
!
! psi is held at <= 0 in (psi/b)**c, as for the Taylor solver.
!
! MDK, 21/10/26
! ---------------------------------------------------------------------
CASE(psi_solver_newton)
  ! Calculate the lower incomplete gamma function for the root zone water
  ! potential, for all land points.
  gamma_root_psi = incomplete_gamma_different_a(land_pts, 1/conductance_c,     &
                     (MIN(root_zone_psi,0.0_real_jlslsm)/conductance_b)        &
                     **conductance_c)

  DO j = 1, open_pnts ! Iterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        * EXP( -(MIN(root_zone_psi(l),0.0_real_jlslsm)/conductance_b(l))       &
               **conductance_c(l) )

    ! The initial guess for the leaf water potential is that of the root
    ! zone meaning that the transpiration rate is zero. This simplifies the
    ! first prediction of the leaf water potential to;
    ! (leaf_k floored as in the update below: in very dry soil it
    ! underflows to 0, and the first guess was then -Inf and the sample
    ! NaN. Only that case changes.)
    !
    ! MDK, 21/10/26
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j)                             &
                    / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

    ! Convergence tolerance for this point, matching the psi_solver_taylor
    ! tolerance above (10% of kcrit).
    k_conversion_limit = 0.1 * kcrit(l)

    ! Iterate the water potential until converged, or until max_nr_iter
    ! is reached. A sample whose e_leaf demand exceeds what this
    ! vulnerability curve can ever supply from root_zone_psi (no finite
    ! leaf_psi satisfies it) does not converge: leaf_k keeps shrinking
    ! and leaf_psi keeps stepping further negative, which drives leaf_k
    ! below kcrit so the feasibility mask (k > kcrit) rejects the sample.
    ! Each sample stops when its own k has settled.
    !
    ! MDK, 24/10/26
    l_nr_active(:) = .TRUE.
    DO i = 1, max_nr_iter

      ! Calculate the xylem conductance for the current leaf water
      ! potential.
      leaf_k_prev(:,j) = leaf_k(:,j)
      leaf_k(:,j) = kmax(l)                                                    &
          * EXP( -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))        &
                 **conductance_c(l) )

      ! Calculate the transpiration rate for the current leaf water
      ! potential. First get the incomplete gamma function for the leaf
      ! water potential.
      e_leaf_current(:,j) = incomplete_gamma(n_e_leaf,                         &
                              1/conductance_c(l),                              &
                              (MIN(leaf_psi(:,j),0.0_real_jlslsm)              &
                               /conductance_b(l))**conductance_c(l))

      ! Calculate the difference between the incomplete gamma function at
      ! the leaf and root zone water potentials.
      e_leaf_current(:,j) = e_leaf_current(:,j) - gamma_root_psi(l)

      ! Multiply by kmax * (-b/c) to get the transpiration rate
      ! Note: The negative sign is present because conductance_b is
      ! negative.
      e_leaf_current(:,j) = e_leaf_current(:,j) * kmax(l)                      &
          * (-conductance_b(l)/conductance_c(l))

      ! Calculate the new leaf water potential. leaf_k is floored in the
      ! denominator to avoid a literal divide-by-zero for a sample that has
      ! already collapsed leaf_k to (numerically) 0; MAX(...) here does not
      ! change the update for any sample that is still within reach of a
      ! finite solution.
      WHERE (l_nr_active(:))                                                   &
        leaf_psi(:,j) = leaf_psi(:,j)                                          &
               - (e_leaf(:,j) - e_leaf_current(:,j))                           &
                 / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

      ! Guard against runaway steps for an infeasible sample: once leaf_k
      ! has collapsed there is no need to keep tracking the true (very
      ! negative) asymptotic solution, only to land clearly below kcrit
      ! without risking overflow in (leaf_psi/b)**c on the next pass.
      leaf_psi(:,j) = MAX(leaf_psi(:,j),                                       &
                          root_zone_psi(l) - 5.0 * ABS(conductance_b(l)))

      WHERE (ABS(leaf_k(:,j) - leaf_k_prev(:,j)) < k_conversion_limit)         &
        l_nr_active(:) = .FALSE.
      IF ( .NOT. ANY(l_nr_active(:)) ) EXIT

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax(l)                                                      &
          * EXP( -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))        &
                 **conductance_c(l) )

  END DO


! ---------------------------------------------------------------------
! Invert the tabulated supply function; see psi_leaf_lut_jls.
! ---------------------------------------------------------------------
CASE(psi_solver_lut)
  CALL psi_leaf_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,        &
                         open_index, e_leaf, root_zone_psi, kmax,              &
                         leaf_psi, leaf_k )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_solver should be psi_solver_taylor, psi_solver_newton ' //    &
        'or psi_solver_lut')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_CW_jls

! ---------------------------------------------------------------------
! Function to calculate the lower incomplete gamma function, with its own
! shape parameter a for each x (incomplete_gamma below is the same with
! one a for all x). The lower incomplete gamma function is defined as;
!
! gamma(a, x) = integral( t^(a-1) * e^(-t) dt, 0, x)
!
! The code here is based on the implementation in Numerical Recipes
! in FORTRAN 90, page 211, and makes use of the series representation
! of the lower incomplete gamma function. JBaguley
!
! The series (Numerical Recipes, section 6.2) is
!
!   gamma(a, x) = x^a e^-x  SUM_n=0..inf  x^n / ( a (a+1) ... (a+n) )
!
! Each term is the previous one times x / (a + n), which is what the loop
! below accumulates (step), before multiplying by x^a e^-x at the end.
!
! What it is for, in plain terms. The water a plant path can move is the
! area under its vulnerability curve: conductance k(psi) summed over the
! drop in water potential from the soil to the leaf (Sperry et al. 2017),
!
!   E = integral( k(psi) dpsi, psi_leaf, psi_soil )
!
! With the cumulative Weibull curve k(psi) = kmax exp(-(psi/b)^c) this area
! has no simple formula, but the change of variable t = (psi/b)^c turns it
! into the standard integral above, with a = 1/c and x = (psi/b)^c:
!
!   integral( exp(-(psi'/b)^c) dpsi', psi, 0 ) = (-b/c) gamma(1/c, x)
!
! So kmax (-b/c) gamma(1/c, (psi/b)^c) is the water the path can carry from
! psi = 0 down to psi (the "supply integral"), and E between two water
! potentials is the difference of two such values. The leaf water
! potential is then found by asking which psi_leaf gives the required E.
! gamma rises from 0 (psi = 0) towards Gamma(1/c) as psi falls, so the
! supply integral levels off: past a point the path cannot carry more
! water however far psi_leaf falls.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma_different_a( n_x, a, x ) RESULT( gamma )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_x ! Size of x array

REAL(KIND=real_jlslsm), INTENT(IN) :: a(n_x)
                            ! Shape parameter of each value, 1/c
                            ! (dimensionless).
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm), INTENT(IN) :: x(n_x)
                            ! Upper limit of the integral, (psi/b)^c
                            ! (dimensionless).

REAL(KIND=real_jlslsm) :: gamma(n_x)
                            ! Lower incomplete gamma of each x
                            ! (dimensionless).
                            !
                            ! MDK, 21/10/26

! Local variables
REAL(KIND=real_jlslsm) :: step(n_x)
LOGICAL :: l_active(n_x)
                            ! .TRUE. while x's own series has not converged.
                            !
                            ! MDK, 23/10/26
INTEGER :: i

! Series convergence controls. NOTE: the series is only well conditioned
! for x < a+1 or so (Numerical Recipes in Fortran 90, section 6.2). Deep
! into the vulnerability curve it converges slowly and its terms grow
! large. If that becomes a problem, pair it with a continued-fraction
! evaluation for x >= a+1 (NR's gcf) rather than raising max_terms.
!
! MDK, 21/10/26
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 40
                            ! Enough for the feasible range of any PFT. The
                            ! terms shrink once n > x, and at the critical
                            ! conductance x = (psi/b)^c = -ln(k_crit/k_max)
                            ! whatever b and c are (about 3 for a 95 % loss),
                            ! so a few tens of terms converge there. Only
                            ! samples far below kcrit need more, and those
                            ! are rejected by the k > kcrit test anyway.
                            ! (For Q. ilex: down to about -10.5 MPa against
                            ! a kcrit at about -8.2 MPa, 17-21 terms.)
                            !
                            ! MDK, 23/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INCOMPLETE_GAMMA_DIFFERENT_A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial values (the n=0 term of the series)
step(:) = 1/a
gamma(:) = step

! Calculate the series. Each x stops adding terms once its own remaining
! term is negligible relative to its running sum (see max_terms and
! series_tol above), so its value does not depend on which other values
! share the call. (It used to stop all of them together, when the last
! had converged, so the same x could differ in the last bits between
! calls.)
!
! MDK, 23/10/26
l_active(:) = .TRUE.
DO i = 1, max_terms
  WHERE ( l_active )
    step = step * x / (a + i)
    gamma = gamma + step
    l_active = ABS(step) >= series_tol * ABS(gamma)
  END WHERE

  IF ( .NOT. ANY(l_active) ) EXIT
END DO

! Multiply by terms outside the sum
gamma = gamma * EXP(-x) * x**a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma_different_a

! ---------------------------------------------------------------------
! The lower incomplete gamma function with one shape parameter a for all
! x: incomplete_gamma_different_a with a repeated. (This was a second copy
! of the series.)
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma( n_x, a, x ) RESULT( gamma )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_x ! Size of x array

REAL(KIND=real_jlslsm), INTENT(IN) :: a
                            ! Shape parameter, 1/c (dimensionless).
REAL(KIND=real_jlslsm), INTENT(IN) :: x(n_x)
                            ! Upper limit of the integral, (psi/b)^c
                            ! (dimensionless).

REAL(KIND=real_jlslsm) :: gamma(n_x)
                            ! Lower incomplete gamma of each x
                            ! (dimensionless).
                            !
                            ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INCOMPLETE_GAMMA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

gamma = incomplete_gamma_different_a( n_x, SPREAD(a, 1, n_x), x )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma

END MODULE xylem_hydraulics_CW_jls_mod
