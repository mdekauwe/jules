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

CHARACTER(LEN=*), PARAMETER, PRIVATE ::                                        &
  ModuleName='XYLEM_HYDRAULICS_SOX_JLS_MOD'

PUBLIC :: leaf_conductance_SOX_jls                                             &
,         xylem_conductance_SOX_jls                                            &
,         psi_leaf_SOX_jls                                                     &
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

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, land_pnts
                            ! Number of land points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(land_pnts)                                                   &
                            ! Leaf water potential for each land point (Pa)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c for each land point

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

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
SUBROUTINE xylem_conductance_SOX_jls( pft,                                     &
                                      n_water_potentials,                      &
                                      open_pnts,                               &
                                      water_potential,                         &
                                      kmax,                                    &
                                      kcrit,                                   &
                                      conductance_b,                           &
                                      conductance_c,                           &
                                   ! INTENT OUT
                                      xylem_conductance                        &
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
                            ! (mol m-2 s-1 Pa-1). Lets a layer- or PFT-specific
                            ! maximum be used rather than always reading
                            ! kmax_pft.
, kcrit(open_pnts)                                                             &
                            ! Critical xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1), scaled to match whatever kmax
                            ! is being used (e.g. kcrit_per_lyr alongside
                            ! kmax_per_lyr).
, conductance_b(open_pnts)                                                     &
                            ! Conductance parameter b for each open point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(open_pnts)
                            ! Conductance parameter c for each open point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: i, j
                            ! Water potential and open point.
                            !
                            ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! k(psi) = kmax / (1 + (psi/b)^c), floored at kcrit (a double loop rather
! than SPREAD copies, as in xylem_conductance_CW_jls).
!
! MDK, 24/10/26
DO j = 1, open_pnts
  DO i = 1, n_water_potentials
    xylem_conductance(i,j) = kmax(j)                                           &
          / (1 + (water_potential(i,j) / conductance_b(j))**conductance_c(j))
    xylem_conductance(i,j) = MAX(xylem_conductance(i,j), kcrit(j))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_SOX_jls

! ---------------------------------------------------------------------
! Subroutine to calculate the leaf water potential from the transpiration
! rate and root zone water potential, with the SOX conductance curve
! (Eller et al. 2018).
!
!                                kmax
!                  k(psi) = ---------------
!                            1 + (psi/b)^c
!
! The code below uses one of two methods to estimate the leaf
! water potential from the transpiration rate. (Method 1, a zeroth-order
! Taylor expansion about the mid-path psi, was retired: it overstated k
! near E_crit. MDK, 25/10/26)
!   2: Uses a Newton-Raphson solve of the integrated conductance
!      equation (psi_solver_newton).
!   3: Inverts the tabulated supply function (psi_solver_lut; see
!      psi_leaf_lut_jls), on the PFT's own curve only (else Newton).
!
! NOTE: This subroutine is designed for the stomatal optimisation model.
!       It returns an output array that only contains values for points
!       with open stomata.
!
! MDK, 24/10/26
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_SOX_jls( pft,                                              &
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

USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_newton,             &
                                psi_solver_lut
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
                            !
                            ! MDK, 24/10/26
, root_zone_psi(land_pts)                                                      &
                            ! Water potential in the root zone (Pa)
, kmax(land_pts)                                                               &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1), e.g. kmax_pft(pft) or a
                            ! layer-scaled value such as kmax_per_lyr.
, kcrit(land_pts)                                                              &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1), scaled to match whatever kmax
                            ! is being used (e.g. kcrit_per_lyr alongside
                            ! kmax_per_lyr).
, conductance_b(land_pts)                                                      &
                            ! Conductance parameter b for each land point
                            ! (Pa). Per point rather than per PFT so an
                            ! impaired vulnerability curve can be used.
                            ! JBaguley
, conductance_c(land_pts)
                            ! Conductance parameter c for each land point.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: i, j, l

REAL(KIND=real_jlslsm) :: k_conv_tol
                            ! Convergence tolerance on k, 0.1 kcrit
                            ! (mol m-2 s-1 Pa-1).
                            !
                            ! MDK, 24/10/26
LOGICAL :: l_nr_active(n_e_leaf)
                   ! Newton: samples still iterating (each sample stops on
                   ! its own convergence)
                   !
                   ! MDK, 24/10/26

! Newton-Raphson variables
REAL(KIND=real_jlslsm) :: e_leaf_current(n_e_leaf, open_pnts)
                            ! Current transpiration rates for each open point
                            ! calculated from current leaf_psi. Used in
                            ! Newton-Raphson method. (mol H2O m-2 s-1)
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: root_psi_component(open_pnts)
                            ! psi_r 2F1 at the root zone water potential of
                            ! each open point. Used in Newton-Raphson method.
REAL(KIND=real_jlslsm) :: cf_open(open_pnts), xf_open(open_pnts)
                            ! Its arguments at the open points: only these
                            ! are solved, so a one-sample call no longer
                            ! sums the series over the whole grid.
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: leaf_k_prev(n_e_leaf, open_pnts)
                            ! leaf_k from the previous Newton-Raphson
                            ! iteration, used for the convergence check.
INTEGER, PARAMETER :: max_nr_iter = 4
                            ! Cap on Newton-Raphson iterations. Mirrors the
                            ! equivalent fix and its reasoning in
                            ! psi_leaf_CW_jls (xylem_hydraulics_cumulative_
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

INTEGER :: psi_solver
                            ! som_psi_solver, or Newton where the supply
                            ! table cannot follow a point's own b and c.
                            !
                            ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_SOX_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The table is keyed by PFT, so a point with its own curve (b, c) is solved
! by Newton.
!
! MDK, 21/10/26
psi_solver = som_psi_solver
IF ( psi_solver == psi_solver_lut ) THEN
  IF ( .NOT. lut_pft_curve( pft, land_pts, open_pnts, veg_index, open_index,   &
                            conductance_b, conductance_c ) )                   &
    psi_solver = psi_solver_newton
END IF

! Select the method to approximate the leaf water potential.
! psi is held at <= 0 in (psi/b)**c: a positive psi gave a negative base and
! NaN. Unchanged for psi <= 0.
!
! MDK, 24/10/26
SELECT CASE (psi_solver)

! ---------------------------------------------------------------------
! Apply a Newton-Raphson method to approximate the leaf water potential
! from the transpiration rate. The transpiration rate is given by the
! integral of the conductance model, k(psi), over the range
! [psi_l, psi_r]. Here the conductance model is that used in the SOX
! 2018 paper. The integral of k(psi) with respect to psi is;
!
!  E(psi_l, psi_r) = kmax * [ psi_r * 2F1(1, 1/c; 1 + 1/c; -(psi_r/b)^c)
!                           - psi_l * 2F1(1, 1/c; 1 + 1/c; -(psi_l/b)^c) ]
!
! Where 2F1() is the Gaussian Hypergeometric series. There is no b/c
! factor: psi * 2F1(1, 1/c; 1 + 1/c; -(psi/b)^c) is itself the integral
! of 1 / (1 + (psi'/b)^c) from 0 to psi. This is what the code computes.
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
! k changes by less than k_conv_tol (or after max_nr_iter steps). Each
! sample stops on its own convergence (l_nr_active).
!
! MDK, 24/10/26
! ---------------------------------------------------------------------
CASE(psi_solver_newton)

  ! psi_r * 2F1 at the root zone water potential, for the open points
  ! (each value converges on its own).
  !
  ! MDK, 24/10/26
  DO j = 1, open_pnts
    l = veg_index(open_index(j))
    cf_open(j) = 1 + 1/conductance_c(l)
    xf_open(j) = -(MIN(root_zone_psi(l),0.0_real_jlslsm)/conductance_b(l))     &
                 **conductance_c(l)
  END DO
  root_psi_component = SOX_2F1_multiple_cf(open_pnts, cf_open, xf_open)
  DO j = 1, open_pnts
    l = veg_index(open_index(j))
    root_psi_component(j) = root_zone_psi(l) * root_psi_component(j)
  END DO

  DO j = 1, open_pnts ! Iterate over the open points.
    l = veg_index(open_index(j))

    ! Do the first iteration outside the loop. By setting the initial guess
    ! to the root zone water potential, the transpiration rate is zero.

    ! Calculate the xylem conductance at the root zone water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        / (1 + (MIN(root_zone_psi(l),0.0_real_jlslsm)/conductance_b(l))        &
           **conductance_c(l))

    ! The initial guess for the leaf water potential is that of the root zone
    ! meaning that the transpiration rate is zero. This simplifies the first
    ! prediction of the leaf water potential to;
    ! (leaf_k floored: it underflows to 0 in very dry soil, which gave a
    ! -Inf first guess and then NaN. Only that case changes. MDK, 22/10/26)
    leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j)                             &
                    / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

    ! Convergence tolerance for this point (10% of kcrit).
    k_conv_tol = 0.1 * kcrit(l)

    ! Iterate the water potential until converged, or until max_nr_iter is
    ! reached. See the max_nr_iter comment above for why a fixed 2
    ! iterations with no convergence check was wrong. Each sample stops
    ! when its own k has settled.
    !
    ! MDK, 24/10/26
    l_nr_active(:) = .TRUE.
    DO i = 1, max_nr_iter

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k_prev(:,j) = leaf_k(:,j)
      leaf_k(:,j) = kmax(l)                                                    &
        / (1 + (MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))           &
           **conductance_c(l))

      ! Calculate the transpiration rate for the current leaf water potential.
      ! First get the 2F1 term for the leaf water potential.
      !
      ! MDK, 24/10/26
      e_leaf_current(:,j) = SOX_2F1(n_e_leaf,                                  &
                                    1 + 1/conductance_c(l),                    &
                   -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))      &
                   **conductance_c(l))

      e_leaf_current(:,j) = leaf_psi(:,j) * e_leaf_current(:,j)

      ! Calculate the difference between psi * 2F1 at the leaf and root
      ! zone water potentials.
      !
      ! MDK, 24/10/26
      e_leaf_current(:,j) = e_leaf_current(:,j) - root_psi_component(j)

      ! Multiply by -kmax to get the transpiration rate,
      !  E = kmax * (psi_r 2F1_r - psi_l 2F1_l).
      !
      ! MDK, 24/10/26
      e_leaf_current(:,j) = - e_leaf_current(:,j) * kmax(l)

      ! Calculate the new leaf water potential. leaf_k is floored in the
      ! denominator to avoid a literal divide-by-zero for a sample that has
      ! already collapsed leaf_k to (numerically) 0; see psi_leaf_CW_jls for
      ! the equivalent reasoning.
      ! (Was "+": the update must lower psi_leaf when the demand exceeds the
      ! supply, as in psi_leaf_CW_jls; with "+" it stepped away from the
      ! root, e.g. -5.0 instead of -5.5 MPa at high demand.)
      WHERE (l_nr_active(:))                                                   &
        leaf_psi(:,j) = leaf_psi(:,j)                                          &
               - (e_leaf(:,j) - e_leaf_current(:,j))                           &
                 / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

      ! Guard against runaway steps for an infeasible sample, as in
      ! psi_leaf_CW_jls.
      leaf_psi(:,j) = MAX(leaf_psi(:,j),                                       &
                          root_zone_psi(l) - 5.0 * ABS(conductance_b(l)))

      WHERE (ABS(leaf_k(:,j) - leaf_k_prev(:,j)) < k_conv_tol)                 &
        l_nr_active(:) = .FALSE.
      IF ( .NOT. ANY(l_nr_active(:)) ) EXIT

    END DO

    ! Calculate the xylem conductance for the current leaf water potential.
    leaf_k(:,j) = kmax(l)                                                      &
        / (1 + (MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))           &
           **conductance_c(l))

  END DO

! Supply-function lookup table (see psi_leaf_lut_jls in supply_table_mod).
!
! MDK, 24/10/26
CASE(psi_solver_lut)
  CALL psi_leaf_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,        &
                         open_index, e_leaf, root_zone_psi, kmax,              &
                         leaf_psi, leaf_k )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_solver should be psi_solver_newton (2) or ' //               &
        'psi_solver_lut (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_SOX_jls

! ---------------------------------------------------------------------
! Function to calculate the Gaussian Hypergeometric series
! 2F1(af,bf;cf;xf). This function is used to calculate the
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
! used to denote the rising factorial;
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
FUNCTION SOX_2F1(n_xf, cf, xf) RESULT( f21 )

! One cf for all values: the per-value series below (SOX_2F1_multiple_cf),
! as incomplete_gamma does.
!
! MDK, 24/10/26

INTEGER, INTENT(IN) :: n_xf ! Size of xf array

REAL(KIND=real_jlslsm), INTENT(IN) :: cf
REAL(KIND=real_jlslsm), INTENT(IN) :: xf(n_xf)

REAL(KIND=real_jlslsm) :: f21(n_xf)

f21 = SOX_2F1_multiple_cf(n_xf, SPREAD(cf, 1, n_xf), xf)

END FUNCTION SOX_2F1

! ---------------------------------------------------------------------
! As SOX_2F1, with cf per value (per-point, e.g. impaired, curves); the
! series itself. JBaguley
! ---------------------------------------------------------------------
FUNCTION SOX_2F1_multiple_cf(n_xf, cf, xf) RESULT( f21 )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_xf ! Size of xf array

REAL(KIND=real_jlslsm), INTENT(IN) :: cf(n_xf)
REAL(KIND=real_jlslsm), INTENT(IN) :: xf(n_xf)

REAL(KIND=real_jlslsm) :: f21(n_xf)
                            ! (Was named sum, which hid the intrinsic.
                            ! MDK, 24/10/26)

! Local variables
REAL(KIND=real_jlslsm) :: current_term(n_xf)
LOGICAL :: l_active(n_xf)
                            ! .TRUE. while this value's own series has not
                            ! converged (as incomplete_gamma).
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm) :: z(n_xf)
INTEGER :: i

! Series convergence controls, mirroring the equivalent fix to
! incomplete_gamma in xylem_hydraulics_cumulative_weibull_jls_mod.f90.
! Same caveat applies: this series is only well conditioned while z is
! not too close to 1 (i.e. xf not too large in magnitude - deep into the
! vulnerability curve, close to hydraulic failure); very extreme xf would
! need more terms than max_terms to converge fully.
!
! How far off: z = xf / (xf - 1) tends to 1 as |psi/b|^c grows, and the
! terms shrink about as z^n. At xf = -25, z = 0.96 and 0.96^100 is about
! 0.017, so E is a few per cent out near and beyond kcrit for curves with
! c <= 2. This matters only with the series-based solvers; for SOX
! curves use the lookup table (som_psi_solver = psi_solver_lut), which
! integrates k(psi) directly.
!
! MDK, 23/10/26
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 100

! zhook variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SOX_2F1_MULTIPLE_CF'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

z = xf / (xf - 1)

! Calculate the 0th term in the sum.
current_term = 1.0
f21 = 1.0

! Tn = [n / (n-1+cf)] * z * Tn-1 (see the derivation above SOX_2F1),
! each value until its own remaining term is negligible relative to its
! running sum (as incomplete_gamma; MDK, 23/10/26).
l_active(:) = .TRUE.
DO i = 1, max_terms
  WHERE ( l_active )
    current_term = current_term * i * z / (i - 1 + cf)
    f21 = f21 + current_term
    l_active = ABS(current_term) >= series_tol * ABS(f21)
  END WHERE

  IF ( .NOT. ANY(l_active) ) EXIT
END DO

! Apply the (1-xf)^(-1) prefactor - see SOX_2F1 above.
f21 = f21 / (1 - xf)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END FUNCTION SOX_2F1_multiple_cf

END MODULE xylem_hydraulics_SOX_jls_mod