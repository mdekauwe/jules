! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_hydraulics_SOX_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_HYDRAULICS_SOX_JLS_MOD'

PUBLIC :: leaf_conductance_SOX_jls                                             &
,         xylem_conductance_SOX_jls                                            &
,         leaf_psi_SOX_jls                                                     &
,         SOX_2F1                                                              &
,         SOX_2F1_multiple_cf

CONTAINS

! *********************************************************************
! Contains routines used to calculate conductance and leaf water
! potential using the conductance model from SOX, Eller CB et al. 2018.
!
!                                kmax
!                  k(psi) = ---------------
!                            1 + (psi/b)^c
! *********************************************************************

! ---------------------------------------------------------------------
! Function to calculate the leaf conductance using the SOX conductance
! model, Eller CB et al. 2018 for all land points. JBaguley
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_SOX_jls( pft,                                      &
                                     land_pnts,                                &
                                     water_potential,                          &
                                     kmax,                                     &
                                     kcrit,                                    &
                                     conductance_b,                            &
                                     conductance_c,                            &
                                    ! INTENT OUT
                                     leaf_k                                    &
  )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                        &
  pft                                                                         &
                            ! Plant functional type index
, land_pnts
                            ! Number of land points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  water_potential(land_pnts)                                                  &
                            ! Leaf water potential for each land point (Pa)
, kmax(land_pnts)                                                             &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                            &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                    &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point (m/s)

! Local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the leaf conductance for each land point.
leaf_k = kmax / (1 + (water_potential/conductance_b)**conductance_c)

! Apply minimum conductance limit
leaf_k = MAX(leaf_k, kcrit)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_SOX_jls

! ---------------------------------------------------------------------
! Function to calculate xylem conductance using the SOX conductance
! model, Eller CB et al. 2018.
! NOTE: This function is designed for the stomatal optimisation model.
!       It assumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_SOX_jls( pft,                                    &
                                      n_water_potentials,                     &
                                      open_pnts,                              &
                                      water_potential,                        &
                                      kmax,                                   &
                                      kcrit,                                  &
                                      conductance_b,                          &
                                      conductance_c,                          &
                                   ! INTENT OUT
                                      xylem_conductance                       &
  )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                        &
  pft                                                                         &
                            ! Plant functional type index
, n_water_potentials                                                          &
                            ! Number of water potentials per open point
, open_pnts
                            ! Number of open stomata

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  water_potential(n_water_potentials, open_pnts)                              &
                            ! Water potentials for each open point)
, kmax(open_pnts)                                                             &
                            ! Maximum xylem conductance for each open point
                            ! (m/s). Lets a layer- or PFT-specific maximum
                            ! be used rather than always reading kmax_pft.
, kcrit(open_pnts)                                                            &
                            ! Critical xylem conductance for each open point
                            ! (m/s), scaled to match whatever kmax is being
                            ! used (e.g. kcrit_per_lyr alongside kmax_per_lyr).
, conductance_b(open_pnts)                                                    &
                            ! Conductance parameter b for each open point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(open_pnts)
                            ! Conductance parameter c for each open point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


xylem_conductance = SPREAD(kmax, DIM = 1, NCOPIES = n_water_potentials)        &
      / (1 + (water_potential                                                  &
              / SPREAD(conductance_b, DIM = 1, NCOPIES = n_water_potentials))  &
             **SPREAD(conductance_c, DIM = 1, NCOPIES = n_water_potentials))


! Apply minimum conductance limit
xylem_conductance = MAX(xylem_conductance,                                     &
                        SPREAD(kcrit, DIM = 1, NCOPIES = n_water_potentials))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_SOX_jls



! ---------------------------------------------------------------------
! Function to calculate the leaf water potential from the transpiration
! rate and rootzone water potential. Uses a cumulative Weibull
! distribution to model conductance.
!
!                                kmax
!                  k(psi) = ---------------
!                            1 + (psi/b)^c
!
! The code bellow uses one of two approximations to estimate the leaf
! water potential from the transpiration rate. This is done to avoid
! the need to integrate the conductance equation each iteration.
!   1: Applies a zeroth order taylor expansion to the conductance
!      equation (psi_solver_taylor).
!   2: Uses a precalculated look up table (psi_solver_newton).
!
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_SOX_jls( pft,                                              &
                             n_e_leaf,                                         &
                             land_pts,                                         &
                             open_pnts,                                        &
                             veg_index,                                        &
                             open_index,                                       &
                             e_leaf,                                           &
                             root_zone_psi,                                    &
                             kmax,                                             &
                             kcrit,                                            &
                             conductance_b,                                    &
                             conductance_c,                                    &
                          ! INTENT OUT
                             leaf_psi,                                         &
                             leaf_k                                            &
  )

USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_taylor,            &
                                psi_solver_newton, psi_solver_lut
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_lut_jls

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
, kcrit(land_pts)                                                              &
                            ! Critical xylem conductance for each land point
                            ! (m/s), scaled to match whatever kmax is being
                            ! used (e.g. kcrit_per_lyr alongside kmax_per_lyr).
, conductance_b(land_pts)                                                      &
                            ! Conductance parameter b for each land point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(land_pts)
                            ! Conductance parameter c for each land point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: i, j, k, l

! Taylor series variables
REAL(KIND=real_jlslsm) :: leaf_k_new
REAL(KIND=real_jlslsm) :: k_conversion_limit
REAL(KIND=real_jlslsm) :: reference_psi
REAL(KIND=real_jlslsm) :: kmax_open(open_pnts)
                            ! kmax gathered onto the open-point index, for
                            ! passing into xylem_conductance_SOX_jls.
REAL(KIND=real_jlslsm) :: kcrit_open(open_pnts)
                            ! kcrit gathered onto the open-point index, for
                            ! passing into xylem_conductance_SOX_jls.
REAL(KIND=real_jlslsm) :: b_open(open_pnts), c_open(open_pnts)
                            ! conductance_b/c gathered onto the open-point
                            ! index, for passing into xylem_conductance_SOX_jls.

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (kg/m2/s)
REAL(KIND=real_jlslsm) :: root_psi_component(land_pts)
                            ! Lower incomplete gamma function for the root
                            ! zone water potential. Used in Newton-Raphson
                            ! method.
REAL(KIND=real_jlslsm) :: leaf_k_prev(n_e_leaf, open_pnts)
                            ! leaf_k from the previous Newton-Raphson
                            ! iteration, used for the convergence check.
INTEGER, PARAMETER :: max_nr_iter = 4
                            ! Cap on Newton-Raphson iterations. Mirrors the
                            ! equivalent fix and its reasoning in
                            ! leaf_psi_CW_jls (xylem_hydraulics_cumulative_
                            ! weibull_jls_mod.f90): this used to be a fixed
                            ! 2 iterations with no convergence check, which
                            ! is inaccurate under water stress and does not
                            ! correctly drive an infeasible sample's leaf_k
                            ! below kcrit. Not independently numerically
                            ! re-verified for this conductance model's
                            ! k(psi) = kmax/(1+(psi/b)^c) form specifically -
                            ! same value used on the assumption the same
                            ! convergence behaviour holds, since the
                            ! functional form is qualitatively the same
                            ! shape (monotonic, collapsing to 0 as psi
                            ! becomes very negative). Re-check numerically
                            ! before relying on this if SOX is put into use.

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
      b_open(j) = conductance_b(l)
      c_open(j) = conductance_c(l)

      ! Calculate the conductance conversion limit
      k_conversion_limit = 0.1 * kcrit(l)

      ! Set the initial leaf water potential to that of the root zone
      leaf_psi(:,j) = root_zone_psi(l)

      ! Calculate the initial xylem conductance for the initial
      !  reference water potential (this is equal to the root
      !  zone water potential).
      leaf_k(:,j) = kmax(l)                                                    &
                 / (1 + (leaf_psi(:,j)/conductance_b(l))**conductance_c(l))

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
          leaf_k_new = kmax(l)                                                 &
                 / (1 + (reference_psi/conductance_b(l))**conductance_c(l))

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
    CALL xylem_conductance_SOX_jls(pft, n_e_leaf, open_pnts, leaf_psi,        &
                                   kmax_open, kcrit_open, b_open, c_open,     &
                                   leaf_k)

    leaf_k(:,:) = MAX(leaf_k(:,:),                                            &
                      SPREAD(kcrit_open, DIM = 1, NCOPIES = n_e_leaf))

! ---------------------------------------------------------------------
! Apply a Newton-Raphson method to approximate the leaf water potential
! from the transpiration rate. The transpiration rate is given by the
! integral of the conductance model, k(psi), over the range
! [psi_r, psi_l]. Here the conductance model is that used in the SOX
! 2018 paper. The integral of k(psi) with respect to psi is;
!
!  E(psi_l, psi_r) = kmax (b/c)
!                     * [ psi_l * 2F1(1, 1/c; 1 + 1/c; -(psi_l/b)^c)
!                        - psi_r * 2F1(1, 1/c; 1 + 1/c; -(psi_r/b)^c)]
!
! Where 2F1() is the Gaussian Hypergeometric series. (This integral
! was calculated using Wolfram Alpha's online integrator)
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

  ! Calculate the lower 2F1 function for the root zone water
  ! potential, for all land points.
  root_psi_component = SOX_2F1_multiple_cf(land_pts,                           &
                               1 + 1/conductance_c,                            &
                -(root_zone_psi/conductance_b)**conductance_c)

  root_psi_component = root_zone_psi * root_psi_component

  DO j = 1, open_pnts ! Iterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop. By setting the initial guess
    ! to the root zone water potential, the transpiration rate is zero.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        / (1 + (root_zone_psi(l)/conductance_b(l))**conductance_c(l))

    ! The initial guess for the leaf water potential is that of the root zone
    ! meaning that the transpiration rate is zero. This simplifies the first
    ! prediction of the leaf water potential to;
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j) / leaf_k(:,j)

    !WRITE(*,*) ' '
    !WRITE(*,*) 'root_psi', root_zone_psi(l)
    !WRITE(*,*) 'leaf_psi', leaf_psi(400,j)

    ! Convergence tolerance for this point, matching the psi_solver_taylor
    ! tolerance above (10% of kcrit).
    k_conversion_limit = 0.1 * kcrit(l)

    ! Iterate the water potential until converged, or until max_nr_iter is
    ! reached. See the max_nr_iter comment above for why a fixed 2
    ! iterations with no convergence check was wrong.
    DO i = 1, max_nr_iter

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k_prev(:,j) = leaf_k(:,j)
      leaf_k(:,j) = kmax(l)                                                    &
        / (1 + (leaf_psi(:,j)/conductance_b(l))**conductance_c(l))

      ! Calculate the transpiration rate for the current leaf water potential.
      ! First get the incomplete gamma function for the leaf water potential.
      e_leaf_current(:,j) = SOX_2F1(n_e_leaf,                                  &
                                    1 + 1/conductance_c(l),                  &
                   -(leaf_psi(:,j)/conductance_b(l))**conductance_c(l))

      e_leaf_current(:,j) = leaf_psi(:,j) * e_leaf_current(:,j)

      ! Calculate the difference between the incomplete gamma function at the
      ! leaf and root zone water potentials.
      e_leaf_current(:,j) = e_leaf_current(:,j) - root_psi_component(l)

      ! Multiply by kmax to get the transpiration rate.
      e_leaf_current(:,j) = - e_leaf_current(:,j) * kmax(l)

      ! Calculate the new leaf water potential. leaf_k is floored in the
      ! denominator to avoid a literal divide-by-zero for a sample that has
      ! already collapsed leaf_k to (numerically) 0; see leaf_psi_CW_jls for
      ! the equivalent reasoning.
      ! (Was "+": the update must lower psi_leaf when the demand exceeds the
      ! supply, as in leaf_psi_CW_jls; with "+" it stepped away from the
      ! root, e.g. -5.0 instead of -5.5 MPa at high demand.)
      leaf_psi(:,j) = leaf_psi(:,j)                                            &
             - (e_leaf(:,j) - e_leaf_current(:,j))                             &
               / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

      ! Guard against runaway steps for an infeasible sample, as in
      ! leaf_psi_CW_jls.
      leaf_psi(:,j) = MAX(leaf_psi(:,j),                                       &
                          root_zone_psi(l) - 5.0 * ABS(conductance_b(l)))

      !WRITE(*,*) 'leaf_psi', leaf_psi(400,j)

      IF (MAXVAL(ABS(leaf_k(:,j) - leaf_k_prev(:,j))) < k_conversion_limit) EXIT

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        / (1 + (leaf_psi(:,j)/conductance_b(l))**conductance_c(l))

  END DO

! Supply-function lookup table (see xylem_hydraulics_CW_jls_mod).
CASE(psi_solver_lut)
  CALL leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,        &
                         open_index, e_leaf, root_zone_psi, kmax,              &
                         leaf_psi, leaf_k )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_solver should be psi_solver_taylor, psi_solver_newton or ' //     &
        'psi_solver_lut')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_SOX_jls

! ---------------------------------------------------------------------
! Function to calculate the Gaussian Hypergeometric series
! 2F1(af,bf;cf;xf). THis function is used to calculate the
! transpiration rate from the leaf water potential for the SOX
! conductance model.
!
! For the SOX conductance model the inputs to the 2F1 function are
! always;
!   af = 1
!   bf = 1/c
!   cf = 1 + 1/c
!   xf = -(psi/b)^c <= 0
!
! The series solution for 2F1(af,bf;cf;xf) is;
!
!                            (af)_n * (bf)_n     xf^n
!   2F1(af,bf;cf;xf) = SUM[ ----------------- * ------ , n=0,inf]
!                                (cf)_n           n!
!
! over the range |xf| < 1. The notation (q)_n is the Pochhammer symbol,
! usd to denote the rising factorial;
!
!    (q)_n = q * (q+1) * (q+2) * ... * (q+n-1)
!
!    (q)_0 = 1
!
! To allow for the calculation of 2F1 for all xf <=0 the code below uses
! the analytic continuation of 2F1 described in equation B.5 of
! W. Becken, P. Schmelcher 1997 (The analytic continuation of the
! Gaussian hypergeometric function 2F1(a,b;c;z) for arbitrary
! parameters). Equation B.5 takes the form;
!
! 2F1(af,bf;cf;xf) / Gamma(cf) =
!      (1 - xf)^(-af) * 2F1(af,cf-bf;cf;xf/(xf-1))/Gamma(cf)
!
! Note that this equation can be used because cf - af - bf = 0.
! Canceling out the two Gamma functions and substituting for the values
! of, bf and cf.
!
! 2F1(1,1/c;1+1/c;xf) =
!     (1 - xf)^(-1) * 2F1(1,1;1+1/c;xf/(xf - 1))
!
! (NOTE: an earlier version of this comment had (1 + xf)^(-1) here, which
! does not match the general formula above with af=1, (1-xf)^(-af). The
! sign was corrected; see also the note below the series recursion.)
!
! Given that xf <= 0, xf/(xf-1) exists over the range [0,1). Meaning
! that the series solution for 2F1(1,1;1+1/c;xf/(xf-1)) can be
! calculated using the series solution for 2F1.
!
!  2F1(1,1;1+1/c;xf/(xf-1)) =
!
!          (1)_n * (1)_n       xf         1
!    SUM[ --------------- * (------)^n * ---- , n=0,inf]
!            (1+1/c)_n        xf-1        n!
!
! Since (1)_n = n! this simplifies to;
!
!  2F1(1,1;1+1/c;xf/(xf-1)) =
!
!              n!          xf
!    SUM[ ----------- * (------)^n, n=0,inf]
!          (1+1/c)_n      xf-1
!
! The value within this sum would be expensive to calculate for large
! n. However, it is possible to re-write the value within the sum such
! that it is a multiple of the previous iteration (n-1). Writing
! z = xf/(xf-1) and Tn for the nth term (T0 = 1), the term-to-term ratio
! is;
!
!            n! / (1+1/c)_n * z^n              n
!  Tn/Tn-1 = ----------------------- = ------------------ * z
!            (n-1)! / (1+1/c)_(n-1) * z^(n-1)   n - 1 + cf
!
! since (1+1/c)_n = (1+1/c)_(n-1) * (n - 1 + cf), where cf = 1 + 1/c (the
! argument passed into this function). So;
!
!                  n
!   Tn = ------------------ * z * Tn-1
!         n - 1 + cf
!
! Where Tn-1 is the previous term in the sum.
!
! (NOTE: an earlier version of this routine used Tn = n * cf * z * Tn-1,
! i.e. multiplying by cf instead of dividing by (n-1+cf). That recursion
! does not converge - for any fixed xf < 0 the terms grow without bound
! as n increases, rather than shrinking, so taking more terms made the
! result worse rather than more accurate. It was also missing the
! (1-xf)^(-1) prefactor derived above. Both are fixed in the code below.)
! ---------------------------------------------------------------------
FUNCTION SOX_2F1(n_xf, cf, xf) RESULT( sum )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_xf ! Size of xf array

REAL(KIND=real_jlslsm), INTENT(IN) :: cf
REAL(KIND=real_jlslsm), INTENT(IN) :: xf(n_xf)

REAL(KIND=real_jlslsm) :: sum(n_xf)

! Local variables
REAL(KIND=real_jlslsm) :: current_term(n_xf)
REAL(KIND=real_jlslsm) :: z(n_xf)
INTEGER :: i

! Series convergence controls, mirroring the equivalent fix to
! incomplete_gamma in xylem_hydraulics_cumulative_weibull_jls_mod.f90.
! Same caveat applies: this series is only well conditioned while z is
! not too close to 1 (i.e. xf not too large in magnitude - deep into the
! vulnerability curve, close to hydraulic failure); very extreme xf would
! need more terms than max_terms to converge fully.
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 100

! zhook variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SOX_2F1'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

z = xf / (xf - 1)

! Calculate the 0th term in the sum.
current_term = 1.0
sum = 1.0

! Iterate over the sum: Tn = [n / (n-1+cf)] * z * Tn-1 (see derivation
! above), until each point's remaining term is negligible relative to its
! running sum.
DO i = 1, max_terms
  current_term = current_term * i * z / (i - 1 + cf)
  sum = sum + current_term

  IF (ALL(ABS(current_term) < series_tol * ABS(sum))) EXIT
END DO

! Apply the (1-xf)^(-1) prefactor from the analytic continuation identity
! above to convert the z-series (2F1(1,1;1+1/c;xf/(xf-1))) into the
! function this routine is documented to return, 2F1(1,1/c;1+1/c;xf).
sum = sum / (1 - xf)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION SOX_2F1

! ---------------------------------------------------------------------
! Same as the function above but with cf as an array input, for use with
! per-point (e.g. impaired) conductance curves. JBaguley
! ---------------------------------------------------------------------
FUNCTION SOX_2F1_multiple_cf(n_xf, cf, xf) RESULT( sum )

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_xf ! Size of xf array

REAL(KIND=real_jlslsm), INTENT(IN) :: cf(n_xf)
REAL(KIND=real_jlslsm), INTENT(IN) :: xf(n_xf)

REAL(KIND=real_jlslsm) :: sum(n_xf)

! Local variables
REAL(KIND=real_jlslsm) :: current_term(n_xf)
REAL(KIND=real_jlslsm) :: z(n_xf)
INTEGER :: i

! Series convergence controls - see SOX_2F1 above.
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 100

! zhook variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SOX_2F1_MULTIPLE_CF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

z = xf / (xf - 1)

! Calculate the 0th term in the sum.
current_term = 1.0
sum = 1.0

! Tn = [n / (n-1+cf)] * z * Tn-1 - see the derivation above SOX_2F1.
DO i = 1, max_terms
  current_term = current_term * i * z / (i - 1 + cf)
  sum = sum + current_term

  IF (ALL(ABS(current_term) < series_tol * ABS(sum))) EXIT
END DO

! Apply the (1-xf)^(-1) prefactor - see SOX_2F1 above.
sum = sum / (1 - xf)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION SOX_2F1_multiple_cf

END MODULE xylem_hydraulics_SOX_jls_mod