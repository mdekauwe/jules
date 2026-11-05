! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE
MODULE supply_table_mod

USE um_types, ONLY: real_jlslsm
USE jules_vegetation_mod, ONLY: k_ground_ratio_max
USE jules_vegetation_mod, ONLY: l_nsc_k_sapwood

IMPLICIT NONE

! Lookup tables of the plant's water supply, used when
! som_psi_solver = psi_solver_lut, and the routines that read them.
!
! A table holds, for one vulnerability curve, how much water the path can
! carry per unit of maximum conductance from psi = 0 down to each water
! potential on a fine grid (S(psi), below). There is one table for the
! whole plant of each PFT and, for plants made of root, stem and leaf
! segments, one for each segment. Reading a table gives the leaf water
! potential for a transpiration rate without iterating.
!
! MDK, 31/10/26
PRIVATE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='SUPPLY_TABLE_MOD'

PUBLIC :: psi_leaf_lut_jls                                                     &
,         psi_leaf_capped_lut_jls                                              &
,         supply_lut_psi                                                       &
,         supply_lut_e_crit                                                    &
,         supply_lut_f                                                         &
,         lut_pft_curve                                                        &
,         init_supply_luts                                                     &
,         seg_lut_psi                                                          &
,         build_supply_lut_seg                                                 &
,         newton_psi_drop                                                      &
,         seg_lut_psi_back                                                     &
,         seg_lut_s

! ---------------------------------------------------------------------
! The supply function. The PFT's vulnerability curve gives the fraction
! of the maximum conductance kmax left at water potential psi (Pa,
! negative), f(psi) = k(psi) / kmax, in one of two forms
! (pft_conductance_model):
!
!   cumulative Weibull: f(psi) = exp(-(psi/b)^c)
!   SOX:                f(psi) = 1 / (1 + (psi/b)^c)
!
! where b (Pa, negative) and c (dimensionless) set the curve. The water
! the path carries from the root zone at psi_r to the leaf at psi_l is the
! area under the conductance curve between them,
!
!   E = kmax (S(psi_l) - S(psi_r))         (mol H2O m-2 s-1)
!
! where S(psi) (Pa) is the area under f from psi up to 0. S depends only
! on b and c, not on kmax, which changes from point to point and in time,
! so it is worked out once per PFT on an evenly spaced psi grid from 0
! down to where f is far below the critical conductance (lut_grid). The
! leaf water potential for a given E is the psi_l where
! S(psi_l) = S(psi_r) + E / kmax, found by searching the table and
! interpolating along a straight line between grid points.
! ---------------------------------------------------------------------
! Nodes in a table that ends at its usual floor (below). A table extended down
! to a lower critical conductance keeps the same spacing and has more nodes, up
! to n_lut_max (lut_grid).
INTEGER, PARAMETER, PRIVATE :: n_lut = 4001
! Most nodes in an extended table; beyond it the spacing is widened instead.
INTEGER, PARAMETER, PRIVATE :: n_lut_max = 16 * (n_lut - 1) + 1
! Bottom of the Weibull table (f falls fast).
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor = 1.0e-6
! Smallest sapwood-area factor on kmax (nsc_k_area in triffid_nsc_mod),
! for the depth of the segment tables with l_nsc_k_sapwood.
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: f_a_min_lut = 0.01
! Bottom of the SOX table, as a fraction of kcrit/kmax (f has a power-law tail,
! so a fixed tiny floor would spread the table too thin).
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor_kcrit = 0.1
! One table. The node count is per table, so a table extended to a low
! sl_kfail_frac does not enlarge the others.
TYPE, PRIVATE :: supply_table
  ! Number of nodes.
  INTEGER :: n = 0
  ! Grid spacing (Pa).
  REAL(KIND=real_jlslsm) :: dpsi = 0.0
  ! S(psi_i) (Pa), psi_i = -(i-1)*dpsi, i = 1..n.
  REAL(KIND=real_jlslsm), ALLOCATABLE :: s(:)
END TYPE supply_table
! Whole-plant table for each PFT.
TYPE(supply_table), ALLOCATABLE, SAVE, PRIVATE :: lut_tab(:)
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_ready(:)
! Table is for the SOX curve (else Weibull).
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_sox(:)
! Root / stem / leaf segment tables (l_som_plant_segments with
! psi_solver_lut): the same S(psi) for each segment's own Weibull
! (conductance_b_seg, conductance_c_seg), so the segments in series are
! solved by three table inversions (psi_leaf_segments_lut_jls).
! Table for (pft, segment).
TYPE(supply_table), ALLOCATABLE, SAVE, PRIVATE :: lut_seg_tab(:,:)
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_seg_ready(:)

CONTAINS

! *****************************************************************************
! Check whether the PFT's supply table can be used for these points: it
! returns .TRUE. if every open point uses the PFT's own vulnerability
! curve (the PFT's b and c). The supply table is built from the PFT's
! curve, so it holds only when this is true. A point with its own curve,
! for example one changed by xylem impairment, is solved by the Newton
! solver instead, which finds psi by repeatedly improving a guess using
! the slope of the flow.
!
! MDK, 21/10/26
! *****************************************************************************
LOGICAL FUNCTION lut_pft_curve( pft, land_pts, open_pnts, veg_index,           &
                                open_index, b, c )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, land_pts, open_pnts, veg_index(land_pts),          &
                       open_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: b(land_pts), c(land_pts)
INTEGER :: j, l

lut_pft_curve = .TRUE.
DO j = 1, open_pnts
  l = veg_index(open_index(j))
  IF ( b(l) /= conductance_b_pft(pft) .OR.                                     &
       c(l) /= conductance_c_pft(pft) ) THEN
    lut_pft_curve = .FALSE.
    RETURN
  END IF
END DO

END FUNCTION lut_pft_curve

! *****************************************************************************
! For each candidate transpiration rate of the stomatal optimiser, read the
! leaf water potential and conductance from the PFT's supply table, which
! is quicker than solving the vulnerability curve each time. For root
! zone water potential
! psi_r (Pa) and transpiration E, the leaf water potential psi_l (Pa) is
! where
!
!   S(psi_l) = S(psi_r) + E / kmax
!
! with S the supply function (top of this module) and kmax the point's
! maximum conductance (mol m-2 s-1 Pa-1). The answer is exact apart from
! the straight-line interpolation between grid points, and needs no
! iteration. The conductance returned is the curve at the leaf,
! k = kmax f(psi_l), which for a single curve is also -dE/dpsi_l.
!
! If E is more than the curve can carry from psi_r, the leaf is put at
! the bottom of the table, where k is at most half the critical
! conductance, so stom_opt_mod rejects the sample, as it does with the
! Newton solver. The routine uses no work arrays, so a call for one
! sample per point is cheap.
!
! MDK, 18/10/26
! *****************************************************************************
SUBROUTINE psi_leaf_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,    &
                                open_index, e_leaf, root_zone_psi, kmax,       &
                                leaf_psi, leaf_k )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pts, open_pnts, veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Transpiration rates for each open point (mol H2O m-2 s-1).
  e_leaf(n_e_leaf, open_pnts)                                                  &
! Water potential in the root zone (Pa).
, root_zone_psi(land_pts)                                                      &
! Maximum xylem conductance (same units as e_leaf per Pa).
, kmax(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Leaf water potential (Pa).
  leaf_psi(n_e_leaf, open_pnts)                                                &
! Xylem conductance at the leaf water potential.
, leaf_k(n_e_leaf, open_pnts)

INTEGER :: i, j, l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_LUT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! For each open point and each candidate transpiration rate: divide E by the
! point's kmax (the table is per unit kmax, so one table serves every point),
! invert the table for psi_leaf, then take k at that psi_leaf from the curve.
! The MAX guards a zero kmax.
DO j = 1, open_pnts
  l = veg_index(open_index(j))
  DO i = 1, n_e_leaf
    leaf_psi(i,j) = supply_lut_psi(pft, root_zone_psi(l),                      &
                    e_leaf(i,j) / MAX(kmax(l), TINY(1.0_real_jlslsm)))
    leaf_k(i,j) = kmax(l) * supply_lut_f(pft, leaf_psi(i,j))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_lut_jls

! *****************************************************************************
! For each candidate transpiration rate of the stomatal optimiser, read the
! leaf water potential and conductance from the supply table for a plant
! whose xylem has been damaged by drought, under the memory model of
! xylem impairment (pft_xylem_impairment_model = 3), for a plant without
! segments. The memory model caps the conductance at kc,
! the most the damaged xylem can carry:
!
!   k(psi) = MIN(kmax f(psi), kc)
!
! The memory model does not change b or c, so the curve f is the PFT's
! and the intact table can be used. psi_r is the root zone water
! potential and psi_l the leaf's (Pa). Let psi_s be where the intact curve
! falls to the cap, f(psi_s) = kc / kmax. Between psi_r and psi_s the
! conductance is the constant kc, so that part of the path carries at
! most
!
!   e_cap = kc (psi_r - psi_s)   if psi_r > psi_s, else 0.
!
! A transpiration E up to e_cap gives psi_l = psi_r - E / kc. A larger E
! continues down the intact curve from psi_top = MIN(psi_r, psi_s):
!
!   S(psi_l) = S(psi_top) + (E - e_cap) / kmax
!
! which is one table look-up. This is the equation the Newton solver in
! psi_leaf_impaired_memory solves. The conductance returned is the capped
! curve at psi_l, floored at the critical conductance kcrit, as there.
!
! MDK, 18/10/26
! *****************************************************************************
SUBROUTINE psi_leaf_capped_lut_jls( pft, n_e_leaf, land_pts, open_pnts,        &
                                    veg_index, open_index, e_leaf,             &
                                    root_zone_psi, kmax_ref, kcap, kcrit,      &
                                    leaf_psi, leaf_k )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pts, open_pnts, veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Transpiration rates for each open point (mol H2O m-2 s-1).
  e_leaf(n_e_leaf, open_pnts)                                                  &
! Water potential in the root zone (Pa).
, root_zone_psi(land_pts)                                                      &
! Intact kmax, the cap and kcrit (same units as e_leaf per Pa).
, kmax_ref(land_pts), kcap(land_pts), kcrit(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Leaf water potential (Pa).
  leaf_psi(n_e_leaf, open_pnts)                                                &
! Capped conductance at leaf_psi (>= kcrit).
, leaf_k(n_e_leaf, open_pnts)

INTEGER :: i, j, l
REAL(KIND=real_jlslsm) :: kmx, kc, psi_r, psi_s, psi_top, e_cap, r

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_CAPPED_LUT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL build_supply_lut(pft)

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  kmx   = kmax_ref(l)
  kc    = MAX(MIN(kcap(l), kmx), kcrit(l))
  psi_r = root_zone_psi(l)

  ! psi_s exactly from the curve (as psi_at_k; +HUGE: no cap).
  IF ( kc < kmx * (1.0 - 1.0e-6) ) THEN
    r = kc / kmx
    IF ( lut_sox(pft) ) THEN
      psi_s = conductance_b_pft(pft)                                           &
              * (1.0 / r - 1.0)**(1.0 / conductance_c_pft(pft))
    ELSE
      psi_s = conductance_b_pft(pft) * (-LOG(r))**(1.0 / conductance_c_pft(pft))
    END IF
  ELSE
    psi_s = HUGE(1.0_real_jlslsm)
  END IF

  IF ( psi_r > psi_s ) THEN
    e_cap = kc * (psi_r - psi_s)
  ELSE
    e_cap = 0.0
  END IF
  psi_top = MIN(psi_r, psi_s)

  DO i = 1, n_e_leaf
    IF ( psi_r > psi_s .AND. e_leaf(i,j) <= e_cap ) THEN
      leaf_psi(i,j) = psi_r - e_leaf(i,j) / kc
    ELSE
      leaf_psi(i,j) = supply_lut_psi(pft, psi_top,                             &
                      (e_leaf(i,j) - e_cap) / MAX(kmx, TINY(1.0_real_jlslsm)))
    END IF
    leaf_k(i,j) = MAX(MIN(kmx * supply_lut_f(pft, leaf_psi(i,j)), kc), kcrit(l))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_capped_lut_jls

! *****************************************************************************
! Return the leaf water potential psi_l (Pa) for one transpiration rate,
! read from the PFT's supply table, for root zone water potential psi_root
! (Pa) and transpiration per unit maximum conductance,
! e_over_kmax = E / kmax (Pa). It solves
!
!   S(psi_l) = S(psi_root) + e_over_kmax
!
! (see psi_leaf_lut_jls). The table is built first if it is not ready.
! *****************************************************************************
FUNCTION supply_lut_psi( pft, psi_root, e_over_kmax ) RESULT( psi_l )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, e_over_kmax
REAL(KIND=real_jlslsm) :: psi_l

CALL build_supply_lut(pft)
psi_l = table_psi(lut_tab(pft), psi_root, e_over_kmax)

END FUNCTION supply_lut_psi

! *****************************************************************************
! Return the most transpiration per unit maximum conductance (E / kmax,
! Pa) the PFT's curve can carry from the root zone at psi_root (Pa) before
! the conductance at the leaf falls to the fraction kcrit_frac of kmax,
! which sets the top of the stomatal optimiser's search. With
! psi_crit the water potential where f(psi_crit) = kcrit_frac, found
! from the curve's formula,
!
!   e_crit = S(psi_crit) - S(psi_root)
!
! It is 0 when the root zone is already drier than psi_crit.
!
! MDK, 18/10/26
! *****************************************************************************
FUNCTION supply_lut_e_crit( pft, psi_root, kcrit_frac ) RESULT( e_crit )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, kcrit_frac
REAL(KIND=real_jlslsm) :: e_crit

REAL(KIND=real_jlslsm) :: psi_crit, r

CALL build_supply_lut(pft)

! f(psi_crit) = kcrit_frac, kept inside (0, 1) (e.g. kcrit >= kmax would
! otherwise give a NaN).
r = MIN(MAX(kcrit_frac, 1.0e-6_real_jlslsm),                                   &
        1.0_real_jlslsm - 1.0e-6_real_jlslsm)

! psi_crit by inverting the curve exactly (not from the table):
!
!   SOX:      f = 1/(1 + (psi/b)^c)  gives  psi = b (1/f - 1)^(1/c),
!   Weibull:  f = exp(-(psi/b)^c)    gives  psi = b (-ln f)^(1/c).
!
IF ( lut_sox(pft) ) THEN
  psi_crit = conductance_b_pft(pft)                                            &
             * (1.0 / r - 1.0)**(1.0 / conductance_c_pft(pft))
ELSE
  psi_crit = conductance_b_pft(pft)                                            &
             * (-LOG(r))**(1.0 / conductance_c_pft(pft))
END IF

! The flow (per unit kmax) that takes the leaf from psi_root down to
! psi_crit. It uses the same interpolation as supply_lut_psi, so inverting
! e_crit lands back on psi_crit. Zero if the root zone is already drier
! than psi_crit: then no flow keeps k above kcrit.
e_crit = MAX(table_s_at(lut_tab(pft), psi_crit)                                &
             - table_s_at(lut_tab(pft), psi_root), 0.0)

END FUNCTION supply_lut_e_crit

! Return the fraction of the maximum conductance left at water potential
! psi (Pa), f = k / kmax, on the PFT's curve (Weibull or SOX), which the
! stomatal optimiser uses for the conductance at one water potential. It
! is worked out from the curve's formula; the table is built only to
! record which form the PFT uses.
!
! MDK, 04/11/26
FUNCTION supply_lut_f( pft, psi ) RESULT( f )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: f

! Evaluated exactly from the curve, not interpolated from the table; the
! build call only makes sure lut_sox(pft) is set. Callers use it for k at
! the psi_leaf from supply_lut_psi (k = kmax * f). b is negative, and psi
! is held at <= 0 as in the Newton solver (psi_leaf_segments_jls), so
! psi/b >= 0.
CALL build_supply_lut(pft)
IF ( lut_sox(pft) ) THEN
  f = 1.0 / (1.0 + (MIN(psi, 0.0_real_jlslsm) / conductance_b_pft(pft))       &
                   **conductance_c_pft(pft))
ELSE
  f = EXP(-(MIN(psi, 0.0_real_jlslsm) / conductance_b_pft(pft))               &
          **conductance_c_pft(pft))
END IF

END FUNCTION supply_lut_f

! *****************************************************************************
! Routines shared by the whole-plant tables (lut_tab) and the segment
! tables (lut_seg_tab), so every table is read, inverted and built the
! same way:
!
!   table_s_at:  S(psi) at any psi, interpolated between grid points;
!   table_psi:   the psi_out where S(psi_out) = S(psi_in) + E / kmax;
!   table_build: S at every grid point, by Simpson's rule.
! *****************************************************************************

! Return the supply function S(psi) (Pa) at water potential psi (Pa), read
! from table tab by straight-line interpolation between its grid points.
! The table must be built. psi above 0
! gives S = 0, and psi below the bottom of the table gives S at the
! bottom.
!
! MDK, 20/10/26
FUNCTION table_s_at( tab, psi ) RESULT( s_psi )

TYPE(supply_table), INTENT(IN) :: tab
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: s_psi

INTEGER :: i0
REAL(KIND=real_jlslsm) :: r

! r is psi's fractional position in the table: 0 at psi = 0, n - 1 at
! the bottom. It is clamped, so psi > 0 gives S = 0 and psi below the table
! gives S at the bottom (no extrapolation). i0 is the upper (wetter) node of
! psi's cell, and S is interpolated linearly between nodes i0 and i0+1.
r = MIN(MAX(-psi / tab%dpsi, 0.0), REAL(tab%n - 1))
i0 = MIN(INT(r) + 1, tab%n - 1)
s_psi = tab%s(i0) + (r - (i0 - 1)) * (tab%s(i0+1) - tab%s(i0))

END FUNCTION table_s_at

! *****************************************************************************
! Return the water potential psi_out (Pa) at the leaf-side end of a path
! that carries a given flow, read from table tab. psi_in (Pa) is the water
! potential at the root-side end and e_over_k = E / kmax (Pa) the flow per
! unit maximum conductance:
!
!   S(psi_out) = S(psi_in) + e_over_k
!
! The grid cell that holds the answer is found by stepping down from
! psi_in's cell in steps that double (1, 2, 4, ... grid points), then
! halving the range until one cell is left. psi_out is interpolated along
! a straight line inside that cell. Most flows move psi_out only a few
! cells, so this takes few steps. psi_out is never wetter than psi_in. A
! flow more than the table can carry puts psi_out at the bottom of the
! table.
!
! MDK, 04/11/26
! *****************************************************************************
FUNCTION table_psi( tab, psi_in, e_over_k ) RESULT( psi_out )

TYPE(supply_table), INTENT(IN) :: tab
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

INTEGER :: lo, hi, mid, i_in, n, step
REAL(KIND=real_jlslsm) :: dpsi, r_in, s_target

dpsi = tab%dpsi
n = tab%n

! s_target is the S the far end must reach: E = kmax * (S(psi_out) -
! S(psi_in)), rearranged. r_in is psi_in's fractional position in the table
! and i_in the upper (wetter) node of its cell, clamped as in table_s_at.
r_in = MIN(MAX(-psi_in / dpsi, 0.0), REAL(n - 1))
i_in = MIN(INT(r_in) + 1, n - 1)
s_target = table_s_at(tab, psi_in) + e_over_k

! Flow at or beyond all the table can carry: put the far end at the bottom
! of the table (or at psi_in, if that is lower). k there is at most half
! of kcrit (and of sl_kfail_frac * kmax, when set), so the sample is
! rejected (stom_opt_mod) or counted as failed (supply-loss).
IF ( s_target >= tab%s(n) ) THEN
  psi_out = MIN(-(n - 1) * dpsi, psi_in)
ELSE
  ! S(lo) <= s_target < S(hi). S does not decrease with depth, but in
  ! single precision it is flat over some cells deep in the tail (below
  ! kcrit), so it is not strictly increasing. The cell is still unique (lo
  ! the last node with S <= s_target), and S(hi) > S(lo) there.
  ! Search for that cell starting from psi_in's cell (the far end is never
  ! wetter). The bracket is first
  ! widened from psi_in's cell in steps of 1, 2, 4, ... nodes, so a flow
  ! that takes the far end only a few cells down (most of them) needs a few
  ! steps on nearby nodes instead of about log2(n) = 12 across the table;
  ! then a binary search within the bracket. The cell, and so psi_out, is
  ! the same as from a binary search over the whole table.
  lo = i_in
  step = 1
  DO WHILE ( lo + step < n )
    IF ( tab%s(lo + step) > s_target ) EXIT
    lo = lo + step
    step = 2 * step
  END DO
  hi = MIN(lo + step, n)
  DO WHILE ( hi - lo > 1 )
    mid = (lo + hi) / 2
    IF ( tab%s(mid) <= s_target ) THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO

  ! Linear interpolation of psi within the cell (the inverse of table_s_at):
  ! node lo is at psi = -(lo-1)*dpsi, and the fraction of the cell's S rise
  ! that s_target covers sets how far below it psi_out lies. The cap at
  ! psi_in stops rounding putting the far end above the near end.
  psi_out = -(lo - 1) * dpsi - dpsi * (s_target - tab%s(lo))                   &
            / (tab%s(hi) - tab%s(lo))
  psi_out = MIN(psi_out, psi_in)
END IF

END FUNCTION table_psi

! *****************************************************************************
! Fill table tab for the curve with parameters b (Pa) and c (SOX if l_sox,
! else Weibull). lut_grid sets the grid, which ends where f = f_floor, or
! further down where f = f_need if that is smaller. S starts at S(0) = 0,
! and each grid cell from p0 down to p1 adds the area under f across it by
! Simpson's rule (f at the two ends and the middle, weighted 1:4:1):
!
!   S(p1) = S(p0) + dpsi / 6 (f(p0) + 4 f(pm) + f(p1))
!
! with dpsi the grid spacing (Pa) and pm the cell's midpoint. The sum and
! f are worked out in double precision and each S rounded once to the
! table's precision, so the error does not build up along the table. f > 0,
! so S never falls as psi falls, which the search in table_psi needs (deep
! in the tail, where a cell adds less than S's last digit, neighbouring
! entries can be equal).
!
! MDK, 04/11/26
! *****************************************************************************
SUBROUTINE table_build( tab, b, c, l_sox, f_floor, f_need )

TYPE(supply_table), INTENT(IN OUT) :: tab
REAL(KIND=real_jlslsm), INTENT(IN) :: b, c, f_floor, f_need
LOGICAL, INTENT(IN) :: l_sox

! Kind of the double-precision sum.
INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15, 307)

INTEGER :: i, n
REAL(KIND=real_jlslsm) :: dpsi
REAL(KIND=r8) :: dpsi8, p0, pm, p1, s_sum

CALL lut_grid( b, c, l_sox, f_floor, f_need, dpsi, n )
tab%dpsi = dpsi
tab%n = n
ALLOCATE( tab%s(n) )
tab%s(1) = 0.0
dpsi8 = REAL(dpsi, r8)
s_sum = 0.0_r8
DO i = 2, n
  p0 = -REAL(i - 2, r8) * dpsi8
  p1 = -REAL(i - 1, r8) * dpsi8
  pm = 0.5_r8 * (p0 + p1)
  s_sum = s_sum + dpsi8 / 6.0_r8 * ( f_of(p0) + 4.0_r8 * f_of(pm) + f_of(p1) )
  tab%s(i) = REAL(s_sum, real_jlslsm)
END DO

CONTAINS

  ! Return the fraction of the maximum conductance left at water potential
  ! p (Pa), f = k / kmax, on the curve being tabulated (as supply_lut_f),
  ! for the integration that fills the table, in double precision.
  REAL(KIND=r8) FUNCTION f_of(p)
  ! Water potential (Pa, <= 0).
  REAL(KIND=r8), INTENT(IN) :: p
  IF ( l_sox ) THEN
    f_of = 1.0_r8 / (1.0_r8 + (MIN(p, 0.0_r8) / REAL(b, r8))**REAL(c, r8))
  ELSE
    f_of = EXP(-(MIN(p, 0.0_r8) / REAL(b, r8))**REAL(c, r8))
  END IF
  END FUNCTION f_of

END SUBROUTINE table_build

! *****************************************************************************
! Build the supply tables of every PFT, and the segment tables when the
! plant is made of root, stem and leaf segments (l_som_plant_segments).
! It is called once at initialisation, after the PFT parameters are read
! and checked (init_pftparm_jules). Building the tables here means none is
! built inside a parallel (OpenMP) region, and every "ready" flag is set
! before any thread reads it. Tables from an earlier call are discarded
! first, so a second initialisation with other parameters rebuilds them
! all. The table routines still build a missing table on first use, for
! callers that do not come through here.
!
! MDK, 21/10/26
! *****************************************************************************
SUBROUTINE init_supply_luts( npft )

USE jules_vegetation_mod, ONLY: l_som_plant_segments

INTEGER, INTENT(IN) :: npft
INTEGER :: pft

IF ( ALLOCATED(lut_tab) ) DEALLOCATE( lut_tab )
IF ( ALLOCATED(lut_ready) ) DEALLOCATE( lut_ready )
IF ( ALLOCATED(lut_sox) ) DEALLOCATE( lut_sox )
IF ( ALLOCATED(lut_seg_tab) ) DEALLOCATE( lut_seg_tab )
IF ( ALLOCATED(lut_seg_ready) ) DEALLOCATE( lut_seg_ready )
DO pft = 1, npft
  CALL build_supply_lut(pft)
  IF ( l_som_plant_segments ) CALL build_supply_lut_seg(pft)
END DO

END SUBROUTINE init_supply_luts

! *****************************************************************************
! Build the supply table of one PFT, once, for the PFT's curve (Weibull or
! SOX, pft_conductance_model):
!
!   lut_tab(pft)%s(i) = integral of f(psi) from psi_i up to 0     (Pa)
!
! with psi_i the i-th grid point (Pa, negative). Each grid cell's part is
! found by Simpson's rule (table_build). The closed formulas for this
! integral (the incomplete gamma function for Weibull, a hypergeometric
! series for SOX) lose accuracy deep into the curve, so they are not used.
!
! In plain terms: the water a plant path can carry depends on how its
! conductance falls as the water potential drops (the vulnerability
! curve, f = k / kmax). Adding up f from psi = 0 down to some psi gives
! S(psi), the water the path can move across that drop per unit kmax. S
! depends only on the shape of the PFT's curve, not on the weather, the
! soil or kmax, so it is worked out once on a fine grid and stored. Then:
!
!   - transpiration between two water potentials is
!       E = kmax (S(psi_leaf) - S(psi_soil));
!   - for a given E the leaf water potential is the psi where
!       S(psi_leaf) = S(psi_soil) + E / kmax,
!     found by searching the table (supply_lut_psi);
!   - the most the path can carry before the conductance at the leaf
!     falls to a given fraction of kmax comes from the same table
!     (supply_lut_e_crit).
!
! S is smooth and rises steadily as psi falls, so straight-line
! interpolation between close grid points is accurate, and much cheaper
! than working out the integral again for every sample at every step.
!
! The grid ends where f falls to a floor: lut_f_floor for Weibull, and
! lut_f_floor_kcrit times the critical fraction kcrit / kmax for SOX,
! whose curve has a long tail. It goes further down when the failure
! conductance is lower (lut_grid). The first call allocates the tables of
! all PFTs. A table is built inside a critical section, so two threads
! cannot build it at once, and is marked ready only when it is complete.
!
! MDK, 23/10/26
! *****************************************************************************
SUBROUTINE build_supply_lut( pft )

USE max_dimensions, ONLY: npft_max
USE pftparm, ONLY: conductance_b_pft, conductance_c_pft,                       &
                   pft_conductance_model, kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: SOX_conductance
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Subroutine arguments
!-----------------------------------------------------------------------------
! Plant functional type whose table is built.
INTEGER, INTENT(IN) :: pft

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
! Vulnerability curve parameter b (Pa).
REAL(KIND=real_jlslsm) :: b
! Vulnerability curve shape parameter c (dimensionless).
REAL(KIND=real_jlslsm) :: c
! Conductance fraction (k / kmax) at the bottom of the table.
REAL(KIND=real_jlslsm) :: f_floor

! .TRUE.: SOX conductance curve; .FALSE.: cumulative Weibull.
LOGICAL :: l_sox

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='BUILD_SUPPLY_LUT'

! Fast exit when this PFT's table is already built. Every supply_lut_*
! call passes through here, so after the first call this is the only work
! done. The check is outside the critical section so that it costs no
! locking.
IF ( ALLOCATED(lut_ready) ) THEN
  IF ( lut_ready(pft) ) RETURN
END IF

! (The DrHook calls come after the fast exit, which is taken on almost
! every call, so that every IN has a matching OUT.)
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First call for any PFT: allocate the tables for all PFTs (npft_max) at
! once, none of them built yet. Inside the critical section, so two threads
! cannot both allocate or both build the same PFT's table.
!$OMP CRITICAL (som_supply_lut)
IF ( .NOT. ALLOCATED(lut_ready) ) THEN
  ! (Each PFT's S array is allocated when its table is built.)
  ALLOCATE( lut_tab(npft_max), lut_ready(npft_max), lut_sox(npft_max) )
  lut_ready(:) = .FALSE.
END IF

IF ( .NOT. lut_ready(pft) ) THEN
  ! Build from the PFT's own vulnerability curve (b, c) and conductance
  ! model. kmax is not needed: the table holds S for f = k/kmax, and the
  ! callers divide E by each point's kmax.
  b = conductance_b_pft(pft)
  c = conductance_c_pft(pft)
  l_sox = pft_conductance_model(pft) == SOX_conductance

  ! Bottom of the table: where f = f_floor.
  ! The grid runs from psi = 0 (node 1) down to that psi (node n_lut) in
  ! n_lut - 1 equal cells of dpsi. The psi comes from inverting the curve
  ! at f_floor (as in supply_lut_e_crit). The SOX floor scales with
  ! kcrit/kmax (lut_f_floor_kcrit); Weibull uses the fixed lut_f_floor.
  ! lut_grid extends the table, at the same spacing, if the critical
  ! conductance (kcrit or sl_kfail_frac) lies below that floor.
  IF ( l_sox ) THEN
    f_floor = lut_f_floor_kcrit * (1.0 - kcrit_fractional_loss(pft))
  ELSE
    f_floor = lut_f_floor
  END IF
  CALL table_build( lut_tab(pft), b, c, l_sox, f_floor,                        &
                    lut_f_crit(1.0 - kcrit_fractional_loss(pft),               &
                               1.0_real_jlslsm) )

  ! Mark the table ready only once it is complete, so the fast exit above
  ! never returns a half-built table.
  lut_sox(pft) = l_sox
  lut_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE build_supply_lut

! *****************************************************************************
! Return the leaf-side water potential psi_out (Pa) of segment iseg
! (1 root, 2 stem, 3 leaf) carrying a given flow, read from the segment's
! supply table, for root-side water potential psi_in (Pa) and flow per
! unit of the segment's maximum conductance, e_over_k = E / kmax_s (Pa):
!
!   S_s(psi_out) = S_s(psi_in) + e_over_k
!
! where S_s is the segment's supply table. A negative flow is treated as
! zero.
!
! MDK, 23/10/26
! *****************************************************************************
FUNCTION seg_lut_psi( pft, iseg, psi_in, e_over_k ) RESULT( psi_out )

INTEGER, INTENT(IN) :: pft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

! A negative flow is treated as zero, so the leaf-side end is never wetter
! than the root-side end.
psi_out = table_psi(lut_seg_tab(pft,iseg), psi_in, MAX(e_over_k, 0.0))

END FUNCTION seg_lut_psi

! *****************************************************************************
! Return the leaf-side water potential psi_out (Pa) of segment iseg when
! water flows back towards the soil, e_over_k = E / kmax_s < 0 (Pa), read
! from the segment's supply table. This happens on the
! root and stem when the plant capacitance store gives more water than
! the leaves transpire (l_store_q_path): they carry E - q_store < 0, where
! q_store is the store's flow out. It solves
!
!   S_s(psi_out) = S_s(psi_in) + e_over_k
!
! upwards from psi_in, so psi_out is wetter than psi_in. If the flow is
! more than the segment can carry from psi_in up to psi = 0, psi_out is 0,
! as in the Newton solver (psi_leaf_segments_jls).
!
! MDK, 31/10/26
! *****************************************************************************
FUNCTION seg_lut_psi_back( pft, iseg, psi_in, e_over_k ) RESULT( psi_out )

INTEGER, INTENT(IN) :: pft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

INTEGER :: lo, hi, mid, n
REAL(KIND=real_jlslsm) :: dpsi, r_in, s_target

dpsi = lut_seg_tab(pft,iseg)%dpsi
n = lut_seg_tab(pft,iseg)%n
s_target = table_s_at(lut_seg_tab(pft,iseg), psi_in) + MIN(e_over_k, 0.0)

IF ( s_target <= 0.0 ) THEN
  psi_out = 0.0
ELSE
  ! S(lo) <= s_target < S(hi), searched between psi = 0 and psi_in's cell
  ! (the far end is never drier than psi_in).
  r_in = MIN(MAX(-psi_in / dpsi, 0.0), REAL(n - 1))
  lo = 1
  hi = MIN(INT(r_in) + 2, n)
  DO WHILE ( hi - lo > 1 )
    mid = (lo + hi) / 2
    IF ( lut_seg_tab(pft,iseg)%s(mid) <= s_target ) THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO
  psi_out = -(lo - 1) * dpsi - dpsi * (s_target - lut_seg_tab(pft,iseg)%s(lo)) &
            / (lut_seg_tab(pft,iseg)%s(hi) - lut_seg_tab(pft,iseg)%s(lo))
  psi_out = MIN(MAX(psi_out, psi_in), 0.0_real_jlslsm)
END IF

END FUNCTION seg_lut_psi_back

! *****************************************************************************
! Return the supply function S_s(psi) (Pa) of segment iseg at water
! potential psi (Pa), by straight-line interpolation in the segment's
! table.
!
! MDK, 31/10/26
! *****************************************************************************
FUNCTION seg_lut_s( pft, iseg, psi ) RESULT( s_psi )

INTEGER, INTENT(IN) :: pft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: s_psi

s_psi = table_s_at(lut_seg_tab(pft,iseg), psi)

END FUNCTION seg_lut_s

! *****************************************************************************
! Build the three segment supply tables (root, stem, leaf) of one PFT,
! once, as build_supply_lut does for the whole plant. Each segment has its
! own Weibull curve (conductance_b_seg, conductance_c_seg). A table runs
! from psi = 0 down to where f falls to lut_f_floor, or further if the
! segment's failure conductance lies below that (lut_f_crit, lut_grid).
! A segment's share of the plant's conductance changes with LAI
! (kplant_model = 1) and with the sapwood area (l_nsc_k_sapwood), so each
! table is built deep enough for the largest share these allow.
!
! MDK, 31/10/26
! *****************************************************************************
SUBROUTINE build_supply_lut_seg( pft )

USE max_dimensions, ONLY: npft_max
USE pftparm, ONLY: conductance_b_seg, conductance_c_seg, seg_kfac,             &
                   kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: kplant_model, kplant_model_ground
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: pft

INTEGER :: iseg
REAL(KIND=real_jlslsm) :: b, c
REAL(KIND=real_jlslsm) :: kfac_lut   ! segment share the table is built for

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='BUILD_SUPPLY_LUT_SEG'

! Same structure as build_supply_lut: fast exit once built, allocate for all
! PFTs inside the critical section on the first call, then build this PFT's
! three tables (root, stem, leaf). Segments are Weibull only, so the floor
! is always lut_f_floor.
IF ( ALLOCATED(lut_seg_ready) ) THEN
  IF ( lut_seg_ready(pft) ) RETURN
END IF

! DrHook after the fast exit, as in build_supply_lut.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP CRITICAL (som_supply_lut_seg)
IF ( .NOT. ALLOCATED(lut_seg_ready) ) THEN
  ! (Each table's S array is allocated when it is built.)
  ALLOCATE( lut_seg_tab(npft_max, 3), lut_seg_ready(npft_max) )
  lut_seg_ready(:) = .FALSE.
END IF

IF ( .NOT. lut_seg_ready(pft) ) THEN
  DO iseg = 1, 3
    b = conductance_b_seg(pft,iseg)
    c = conductance_c_seg(pft,iseg)

    ! A segment's k is f kmax seg_kfac, and the whole-path k is at most
    ! the leaf segment's, so the critical fraction of the path (kcrit or
    ! sl_kfail_frac) is a fraction 1/seg_kfac of it on the segment's own
    ! curve. seg_kfac >= 1 for segments in series; the MAX guards others.
    ! With kplant_model = 1 a segment's share also moves with LAI
    ! (seg_kscale), so the table is built for the largest share the
    ! LAI_ref / LAI cap allows, k_ground_ratio_max x seg_kfac.
    ! With the sapwood-area factor (l_nsc_k_sapwood) the stem's share can
    ! fall to f_a_min and the other segments' shares rise by up to its
    ! inverse, so the table is built that far too.
    kfac_lut = MAX(seg_kfac(pft,iseg), 1.0_real_jlslsm)
    IF ( kplant_model == kplant_model_ground )                                 &
      kfac_lut = kfac_lut * k_ground_ratio_max
    IF ( l_nsc_k_sapwood ) kfac_lut = kfac_lut / f_a_min_lut
    CALL table_build( lut_seg_tab(pft,iseg), b, c, .FALSE., lut_f_floor,       &
                      lut_f_crit(1.0 - kcrit_fractional_loss(pft), kfac_lut) )
  END DO
  lut_seg_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut_seg)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE build_supply_lut_seg

! *****************************************************************************
! Return how far down the vulnerability curve a supply table must reach,
! as the lowest conductance fraction f = k / kmax it must hold, so the
! table covers every flow the stomatal optimiser may try. It is half the
! path's failure fraction, divided by kfac, the segment's maximum
! conductance over the whole path's (1 for the whole-plant table). The
! failure fraction is kcrit / kmax = 1 - kcrit_fractional_loss, or
! sl_kfail_frac when that is set and lower. sl_kfail_frac
! (stomata_model = 9) can be far below kcrit / kmax (Sperry et al. 2016
! use 0.0005), so a table that stopped at the usual floor (lut_f_floor,
! lut_f_floor_kcrit) could end before the failure point.
!
! MDK, 20/10/26
! *****************************************************************************
FUNCTION lut_f_crit( f_kcrit, kfac ) RESULT( f_need )

USE jules_vegetation_mod, ONLY: sl_kfail_frac

REAL(KIND=real_jlslsm), INTENT(IN) :: f_kcrit, kfac
REAL(KIND=real_jlslsm) :: f_need

f_need = f_kcrit
IF ( sl_kfail_frac > 0.0 ) f_need = MIN(f_need, sl_kfail_frac)
f_need = 0.5 * f_need / kfac

END FUNCTION lut_f_crit

! *****************************************************************************
! Return the furthest the Newton solvers may move the leaf-side water
! potential below the root-side one (Pa, positive), on the curve with b
! (Pa) and c
! (SOX if l_sox, else Weibull). The Newton solvers find psi by repeatedly
! improving a guess using the slope of the flow; for a demand the curve
! cannot supply, the guesses would run away to ever more negative psi,
! and this floor stops them.
!
! A sample stopped at the floor is rejected only if k there is below the
! failure conductance, so the floor must lie below the point where f
! falls to the failure fraction (f_kcrit = kcrit / kmax, or sl_kfail_frac
! when that is set and lower). That point is at
!
!   Weibull: psi_f = b (-ln f)^(1/c),   SOX: psi_f = b (1/f - 1)^(1/c)
!
! and the drop is 1.2 |psi_f| or 5 |b|, whichever is larger. 5 |b| alone
! lies above the failure point for wide curves (Weibull c below about
! 0.7, SOX c below about 1.8, or a small sl_kfail_frac).
!
! MDK, 25/10/26
! *****************************************************************************
ELEMENTAL FUNCTION newton_psi_drop( b, c, f_kcrit, l_sox ) RESULT( drop )

USE jules_vegetation_mod, ONLY: sl_kfail_frac

REAL(KIND=real_jlslsm), INTENT(IN) :: b, c, f_kcrit
LOGICAL, INTENT(IN) :: l_sox
REAL(KIND=real_jlslsm) :: drop

REAL(KIND=real_jlslsm) :: f, x_f

f = f_kcrit
IF ( sl_kfail_frac > 0.0 ) f = MIN(f, sl_kfail_frac)
f = MIN(MAX(f, 1.0e-12_real_jlslsm), 0.999_real_jlslsm)
IF ( l_sox ) THEN
  x_f = 1.0 / f - 1.0
ELSE
  x_f = -LOG(f)
END IF
drop = ABS(b) * MAX(5.0_real_jlslsm, 1.2 * x_f**(1.0 / c))

END FUNCTION newton_psi_drop

! *****************************************************************************
! Choose the grid of a supply table on the curve with b (Pa) and c (SOX if
! l_sox, else Weibull): the spacing dpsi (Pa) and the number of grid
! points n.
! The usual table has n_lut points and ends where f = f_floor. If f_need
! (lut_f_crit) is lower, the table is made longer at the same spacing,
! down to where f = f_need, so the accuracy near psi = 0, where most of
! the flow is, does not change. Only beyond n_lut_max points (16 times the
! usual length, for example a SOX curve with c near 1 and a very small
! sl_kfail_frac) is the spacing widened instead.
!
! MDK, 20/10/26
! *****************************************************************************
SUBROUTINE lut_grid( b, c, l_sox, f_floor, f_need, dpsi, n )

REAL(KIND=real_jlslsm), INTENT(IN) :: b, c, f_floor, f_need
LOGICAL, INTENT(IN) :: l_sox
REAL(KIND=real_jlslsm), INTENT(OUT) :: dpsi
INTEGER, INTENT(OUT) :: n

REAL(KIND=real_jlslsm) :: psi_need, r

! Spacing of the table that ends at f = f_floor.
IF ( l_sox ) THEN
  dpsi = ABS(b) * (1.0 / f_floor - 1.0)**(1.0 / c) / (n_lut - 1)
ELSE
  dpsi = ABS(b) * (-LOG(f_floor))**(1.0 / c) / (n_lut - 1)
END IF
n = n_lut

IF ( f_need > 0.0 .AND. f_need < f_floor ) THEN
  ! |psi| where f = f_need, in cells of dpsi.
  IF ( l_sox ) THEN
    psi_need = ABS(b) * (1.0 / f_need - 1.0)**(1.0 / c)
  ELSE
    psi_need = ABS(b) * (-LOG(f_need))**(1.0 / c)
  END IF
  r = psi_need / dpsi
  IF ( r < REAL(n_lut_max - 1, real_jlslsm) ) THEN
    n = MAX(n_lut, CEILING(r) + 1)
  ELSE
    n = n_lut_max
    dpsi = psi_need / (n_lut_max - 1)
  END IF
END IF

END SUBROUTINE lut_grid

END MODULE supply_table_mod
