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

IMPLICIT NONE

! Supply-function lookup tables for the plant hydraulics
! (som_psi_solver = psi_solver_lut): the whole-plant table of each PFT and
! the root, stem and leaf segment tables, and the solvers that read them.
! (Moved here from xylem_hydraulics_CW_jls_mod, which now keeps the
! cumulative Weibull routines.)
!
! MDK, 23/10/26
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
,         build_supply_lut_seg

! ---------------------------------------------------------------------
! Supply-function lookup table (som_psi_solver = psi_solver_lut),
! for either conductance model of the PFT (pft_conductance_model):
!   cumulative Weibull: f(psi) = exp(-(psi/b)^c)
!   SOX:                f(psi) = 1 / (1 + (psi/b)^c)
! With k(psi) = kmax * f(psi), the transpiration between root zone and
! leaf is E = kmax * (S(psi_l) - S(psi_r)), S(psi) = integral(f, psi, 0).
! S depends only on the PFT's b and c (not on kmax, which can vary in
! space and time), so it is tabulated once per PFT on a uniform psi grid
! from 0 down to where f has fallen to well below kcrit/kmax (and below
! sl_kfail_frac; see lut_grid). psi_leaf is
! then found by inverting S directly (binary search plus linear
! interpolation) instead of Newton-Raphson on the incomplete gamma or
! hypergeometric series.
! ---------------------------------------------------------------------
INTEGER, PARAMETER, PRIVATE :: n_lut = 4001
                            ! Nodes in a table that ends at its usual floor
                            ! (below). A table extended down to a lower
                            ! critical conductance keeps the same spacing and
                            ! has more nodes, up to n_lut_max (lut_grid).
INTEGER, PARAMETER, PRIVATE :: n_lut_max = 16 * (n_lut - 1) + 1
                            ! Most nodes in an extended table; beyond it the
                            ! spacing is widened instead.
                            !
                            ! MDK, 20/10/26
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor = 1.0e-6
                            ! Bottom of the Weibull table (f falls fast).
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: lut_f_floor_kcrit = 0.1
                            ! Bottom of the SOX table, as a fraction of
                            ! kcrit/kmax (f has a power-law tail, so a fixed
                            ! tiny floor would spread the table too thin).
TYPE, PRIVATE :: supply_table
  INTEGER :: n = 0
                            ! Number of nodes.
  REAL(KIND=real_jlslsm) :: dpsi = 0.0
                            ! Grid spacing (Pa).
  REAL(KIND=real_jlslsm), ALLOCATABLE :: s(:)
                            ! S(psi_i) (Pa), psi_i = -(i-1)*dpsi, i = 1..n.
END TYPE supply_table
                            ! One table. The node count is per table, so a
                            ! table extended to a low sl_kfail_frac does not
                            ! enlarge the others.
                            !
                            ! MDK, 20/10/26
TYPE(supply_table), ALLOCATABLE, SAVE, PRIVATE :: lut_tab(:)
                            ! Whole-plant table for each PFT.
                            !
                            ! MDK, 20/10/26
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_ready(:)
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_sox(:)
                            ! Table is for the SOX curve (else Weibull).
! Root / stem / leaf segment tables (l_som_plant_segments with
! psi_solver_lut): the same S(psi) for each segment's own Weibull
! (conductance_b_seg, conductance_c_seg), so the segments in series are
! solved by three table inversions (psi_leaf_segments_lut_jls).
TYPE(supply_table), ALLOCATABLE, SAVE, PRIVATE :: lut_seg_tab(:,:)
                            ! Table for (pft, segment).
                            !
                            ! MDK, 20/10/26
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_seg_ready(:)

CONTAINS

! *****************************************************************************
! .TRUE. if every open point's curve (b, c) is the PFT's, so the supply
! table (keyed by PFT) holds for all of them. A point with its own curve
! (e.g. impaired b and c) must be solved by Newton instead; the table
! solvers used to ignore b and c.
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
! Leaf water potential and conductance by inverting the tabulated supply
! function (the table lut_tab, see the top of this module):
!   S(psi_l) = S(psi_r) + e_leaf/kmax.
! Exact up to the table's linear interpolation, with no iteration. A demand
! beyond what the curve can supply from psi_r puts the leaf at the bottom of
! the table, where k = lut_f_floor * kmax < kcrit, so stom_opt_mod rejects it
! as for the Newton-Raphson solve. Has no automatic arrays, so it is cheap to
! call for a single sample per point.
! *****************************************************************************
SUBROUTINE psi_leaf_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,    &
                                open_index, e_leaf, root_zone_psi, kmax,       &
                                leaf_psi, leaf_k )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pts, open_pnts, veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 s-1).
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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_LUT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! For each open point and each candidate transpiration rate: divide E by the
! point's kmax (the table is per unit kmax, so one table serves every point),
! invert the table for psi_leaf, then take k at that psi_leaf from the curve.
! The MAX guards a zero kmax.
!
! MDK, 18/10/26
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
! Memory xylem impairment model (pft_xylem_impairment_model = 3) without
! segments, from the supply-function table (psi_solver_lut): the capped
! curve k(psi) = MIN(kmax f(psi), kc) of psi_leaf_impaired_memory, inverted
! exactly from the intact PFT table. With psi_s where the intact curve
! meets the cap (f(psi_s) = kc/kmax), the flow from psi_r is carried at the
! constant conductance kc down to psi_s, e_cap = kc (psi_r - psi_s) if
! psi_r > psi_s (else 0), so
!   E <= e_cap:  psi_l = psi_r - E / kc                     (exact)
!   E >  e_cap:  S(psi_l) = S(psi_top) + (E - e_cap) / kmax, psi_top =
!                MIN(psi_r, psi_s)                          (one inversion)
! as the Newton-Raphson of psi_leaf_impaired_memory solves. leaf_k is the
! capped conductance at psi_l floored at kcrit, as there. The curve must be
! the PFT's (the table's): b and c are not changed by the memory model.
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
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 s-1).
, root_zone_psi(land_pts)                                                      &
                            ! Water potential in the root zone (Pa).
, kmax_ref(land_pts), kcap(land_pts), kcrit(land_pts)
                            ! Intact kmax, the cap and kcrit (same units as
                            ! e_leaf per Pa).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential (Pa).
, leaf_k(n_e_leaf, open_pnts)
                            ! Capped conductance at leaf_psi (>= kcrit).

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
! Leaf water potential (Pa) from the supply-function table, for root zone
! potential psi_root (Pa) and normalised demand e_over_kmax = e_leaf/kmax (Pa):
! solves S(psi_l) = S(psi_r) + e_over_kmax. See psi_leaf_lut_jls.
! *****************************************************************************
FUNCTION supply_lut_psi( pft, psi_root, e_over_kmax ) RESULT( psi_l )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, e_over_kmax
REAL(KIND=real_jlslsm) :: psi_l

CALL build_supply_lut(pft)
psi_l = table_psi(lut_tab(pft), psi_root, e_over_kmax)

END FUNCTION supply_lut_psi

! *****************************************************************************
! Largest normalised transpiration, e/kmax (Pa), the curve can supply from
! psi_root while keeping k >= kcrit_frac * kmax: S(psi_crit) - S(psi_root),
! with k(psi_crit) = kcrit_frac * kmax. Zero if psi_root is already beyond
! psi_crit.
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
! SOX f = 1/(1 + (psi/b)^c) gives psi = b (1/f - 1)^(1/c);
! Weibull f = exp(-(psi/b)^c) gives psi = b (-ln f)^(1/c).
!
! MDK, 18/10/26
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
!
! MDK, 18/10/26
e_crit = MAX(table_s_at(lut_tab(pft), psi_crit)                                &
             - table_s_at(lut_tab(pft), psi_root), 0.0)

END FUNCTION supply_lut_e_crit

! Relative conductance f(psi) = k/kmax of the PFT's conductance model.
FUNCTION supply_lut_f( pft, psi ) RESULT( f )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: f

! Evaluated exactly from the curve, not interpolated from the table; the
! build call only makes sure lut_sox(pft) is set. Callers use it for k at
! the psi_leaf from supply_lut_psi (k = kmax * f). b and psi are both
! negative, so psi/b >= 0; ABS only guards a psi > 0.
!
! MDK, 18/10/26
CALL build_supply_lut(pft)
IF ( lut_sox(pft) ) THEN
  f = 1.0 / (1.0 + (ABS(psi / conductance_b_pft(pft)))**conductance_c_pft(pft))
ELSE
  f = EXP(-(ABS(psi / conductance_b_pft(pft)))**conductance_c_pft(pft))
END IF

END FUNCTION supply_lut_f

! *****************************************************************************
! Helpers shared by the whole-plant tables (lut_tab) and the segment tables
! (lut_seg_tab), so that a table is read, inverted and built in one place:
!
!   table_s_at:  S(psi) by linear interpolation;
!   table_psi:   the psi_out where S(psi_out) = S(psi_in) + e_over_k;
!   table_build: S at every node, by Simpson's rule.
!
! MDK, 23/10/26
! *****************************************************************************

! S(psi) by linear interpolation in table tab (which must be built).
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
!
! MDK, 20/10/26
r = MIN(MAX(-psi / tab%dpsi, 0.0), REAL(tab%n - 1))
i0 = MIN(INT(r) + 1, tab%n - 1)
s_psi = tab%s(i0) + (r - (i0 - 1)) * (tab%s(i0+1) - tab%s(i0))

END FUNCTION table_s_at

! *****************************************************************************
! Water potential psi_out (Pa) at the far end of a path on table tab, from
! psi_in (Pa) at the near end and normalised flow e_over_k = E / kmax (Pa):
! solves S(psi_out) = S(psi_in) + e_over_k. psi_out is never wetter than
! psi_in.
! *****************************************************************************
FUNCTION table_psi( tab, psi_in, e_over_k ) RESULT( psi_out )

TYPE(supply_table), INTENT(IN) :: tab
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

INTEGER :: lo, hi, mid, i_in, n
REAL(KIND=real_jlslsm) :: dpsi, r_in, s_target

dpsi = tab%dpsi
n = tab%n

! s_target is the S the far end must reach: E = kmax * (S(psi_out) -
! S(psi_in)), rearranged. r_in is psi_in's fractional position in the table
! and i_in the upper (wetter) node of its cell, clamped as in table_s_at.
!
! MDK, 18/10/26
r_in = MIN(MAX(-psi_in / dpsi, 0.0), REAL(n - 1))
i_in = MIN(INT(r_in) + 1, n - 1)
s_target = table_s_at(tab, psi_in) + e_over_k

! Flow at or beyond all the table can carry: put the far end at the bottom
! of the table (or at psi_in, if that is lower). k there is at most half
! of kcrit (and of sl_kfail_frac * kmax, when set), so the sample is
! rejected (stom_opt_mod) or counted as failed (supply-loss).
!
! MDK, 20/10/26
IF ( s_target >= tab%s(n) ) THEN
  psi_out = MIN(-(n - 1) * dpsi, psi_in)
ELSE
  ! S(lo) <= s_target < S(hi); S is strictly increasing.
  ! Binary search for that cell, starting from psi_in's cell (the far end
  ! is never wetter). S rises with depth, so the cell is unique; the search
  ! takes about log2(n) = 12 steps (16 for n_lut_max).
  !
  ! MDK, 18/10/26
  lo = i_in
  hi = n
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
  !
  ! MDK, 18/10/26
  psi_out = -(lo - 1) * dpsi - dpsi * (s_target - tab%s(lo))                   &
            / (tab%s(hi) - tab%s(lo))
  psi_out = MIN(psi_out, psi_in)
END IF

END FUNCTION table_psi

! *****************************************************************************
! Fill table tab on the curve (b, c; SOX if l_sox, else Weibull): the grid
! from lut_grid (bottom at f_floor, extended down to f_need), then S
! accumulated downwards from S(0) = 0. Each cell [p1, p0] adds the integral
! of f over it by Simpson's rule (ends and midpoint weighted 1:4:1, times
! dpsi/6). f > 0, so S rises strictly with depth, which the binary search in
! table_psi relies on.
!
! MDK, 23/10/26
! *****************************************************************************
SUBROUTINE table_build( tab, b, c, l_sox, f_floor, f_need )

TYPE(supply_table), INTENT(IN OUT) :: tab
REAL(KIND=real_jlslsm), INTENT(IN) :: b, c, f_floor, f_need
LOGICAL, INTENT(IN) :: l_sox

INTEGER :: i, n
REAL(KIND=real_jlslsm) :: dpsi, p0, pm, p1

CALL lut_grid( b, c, l_sox, f_floor, f_need, dpsi, n )
tab%dpsi = dpsi
tab%n = n
ALLOCATE( tab%s(n) )
tab%s(1) = 0.0
DO i = 2, n
  p0 = -(i - 2) * dpsi
  p1 = -(i - 1) * dpsi
  pm = 0.5 * (p0 + p1)
  tab%s(i) = tab%s(i-1) + dpsi / 6.0 * ( f_of(p0) + 4.0 * f_of(pm) + f_of(p1) )
END DO

CONTAINS

  ! f = k/kmax at psi = p on the curve (as supply_lut_f).
  REAL(KIND=real_jlslsm) FUNCTION f_of(p)
  REAL(KIND=real_jlslsm), INTENT(IN) :: p
                            ! Water potential (Pa, <= 0).
  IF ( l_sox ) THEN
    f_of = 1.0 / (1.0 + (ABS(p / b))**c)
  ELSE
    f_of = EXP(-(ABS(p / b))**c)
  END IF
  END FUNCTION f_of

END SUBROUTINE table_build

! *****************************************************************************
! Build the supply tables of every PFT (and the segment tables with
! l_som_plant_segments) at initialisation, once the PFT parameters are read
! and checked (init_pftparm_jules), instead of on first use inside the
! timestep. So no table is built inside an OpenMP region, and the ready
! flags the table functions test outside the critical sections are all set
! before any thread reads them. Tables from an earlier call are discarded
! first, so a re-initialisation with other parameters cannot keep stale
! ones. (The table functions still build a missing table on first use, for
! callers that do not come through here.)
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
! Build (once per PFT) the supply-function table used by psi_solver_lut:
!
!   lut_tab(pft)%s(i) = integral( f(psi) dpsi, psi_i, 0 )
!
! by Simpson's rule on each grid cell, for the PFT's conductance model (see
! lut_tab at the top of this module). Built by numerical integration rather
! than the closed forms (incomplete gamma / hypergeometric series), whose
! series are poorly conditioned deep into the curve.
!
! In plain terms. The water a plant path can carry depends on how its
! conductance falls as the water potential drops (the vulnerability curve,
! f = k / kmax). Adding up f from psi = 0 down to some psi gives S(psi): how
! much water (per unit kmax) the path can move across that drop. This only
! depends on the shape of the PFT's curve, not on the weather, the soil or
! kmax, so it is worked out once, on a fine grid of psi values from 0 down
! to where the conductance has almost gone, and stored.
!
! Every later call then just reads the table:
!   - transpiration between two water potentials is
!       E = kmax (S(psi_leaf) - S(psi_soil))
!   - so for a given E the leaf water potential is the psi where
!       S(psi_leaf) = S(psi_soil) + E / kmax
!     found by searching the table and interpolating between grid points
!     (supply_lut_psi);
!   - the conductance at any psi is read from the same grid
!     (supply_lut_f), and so is E_crit (supply_lut_e_crit).
! S is smooth and only ever increases with depth, so straight-line
! interpolation between closely spaced grid points is very accurate, and it
! is far cheaper than solving the integral again for every Ci sample at
! every timestep.
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
INTEGER, INTENT(IN) :: pft
                            ! Plant functional type whose table is built.

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) :: b
                            ! Vulnerability curve parameter b (Pa).
REAL(KIND=real_jlslsm) :: c
                            ! Vulnerability curve shape parameter c
                            ! (dimensionless).
REAL(KIND=real_jlslsm) :: f_floor
                            ! Conductance fraction (k / kmax) at the bottom
                            ! of the table.

LOGICAL :: l_sox
                            ! .TRUE.: SOX conductance curve; .FALSE.:
                            ! cumulative Weibull.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='BUILD_SUPPLY_LUT'

! Fast exit when this PFT's table is already built. Every supply_lut_*
! call passes through here, so after the first call this is the only work
! done. The check is outside the critical section so that it costs no
! locking.
!
! MDK, 18/10/26
IF ( ALLOCATED(lut_ready) ) THEN
  IF ( lut_ready(pft) ) RETURN
END IF

! (The DrHook calls come after the fast exit, which is taken on almost
! every call, so that every IN has a matching OUT.)
!
! MDK, 23/10/26
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! First call for any PFT: allocate the tables for all PFTs (npft_max) at
! once, none of them built yet. Inside the critical section, so two threads
! cannot both allocate or both build the same PFT's table.
!
! MDK, 18/10/26
!$OMP CRITICAL (som_supply_lut)
IF ( .NOT. ALLOCATED(lut_ready) ) THEN
  ! (Each PFT's S array is allocated when its table is built; MDK, 20/10/26.)
  ALLOCATE( lut_tab(npft_max), lut_ready(npft_max), lut_sox(npft_max) )
  lut_ready(:) = .FALSE.
END IF

IF ( .NOT. lut_ready(pft) ) THEN
  ! Build from the PFT's own vulnerability curve (b, c) and conductance
  ! model. kmax is not needed: the table holds S for f = k/kmax, and the
  ! callers divide E by each point's kmax.
  !
  ! MDK, 18/10/26
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
  !
  ! MDK, 20/10/26
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
  !
  ! MDK, 18/10/26
  lut_sox(pft) = l_sox
  lut_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE build_supply_lut

! *****************************************************************************
! Leaf-side water potential (Pa) of segment iseg for root-side psi_in (Pa) and
! normalised flow e_over_k = E / kmax_s (Pa): S_s(psi_out) = S_s(psi_in) +
! e_over_k (table_psi on the segment's table).
! *****************************************************************************
FUNCTION seg_lut_psi( pft, iseg, psi_in, e_over_k ) RESULT( psi_out )

INTEGER, INTENT(IN) :: pft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_in, e_over_k
REAL(KIND=real_jlslsm) :: psi_out

! A negative flow is treated as zero, so the leaf-side end is never wetter
! than the root-side end.
!
! MDK, 23/10/26
psi_out = table_psi(lut_seg_tab(pft,iseg), psi_in, MAX(e_over_k, 0.0))

END FUNCTION seg_lut_psi

! *****************************************************************************
! Build (once per PFT) the segment supply-function tables: as
! build_supply_lut for each segment's Weibull, from psi = 0 down to where the
! normalised conductance falls to lut_f_floor (or lower, see lut_grid).
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
                                     ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='BUILD_SUPPLY_LUT_SEG'

! Same structure as build_supply_lut: fast exit once built, allocate for all
! PFTs inside the critical section on the first call, then build this PFT's
! three tables (root, stem, leaf). Segments are Weibull only, so the floor
! is always lut_f_floor.
!
! MDK, 18/10/26
IF ( ALLOCATED(lut_seg_ready) ) THEN
  IF ( lut_seg_ready(pft) ) RETURN
END IF

! DrHook after the fast exit, as in build_supply_lut.
!
! MDK, 24/10/26
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP CRITICAL (som_supply_lut_seg)
IF ( .NOT. ALLOCATED(lut_seg_ready) ) THEN
  ! (Each table's S array is allocated when it is built; MDK, 20/10/26.)
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
    !
    ! MDK, 21/10/26
    kfac_lut = MAX(seg_kfac(pft,iseg), 1.0_real_jlslsm)
    IF ( kplant_model == kplant_model_ground )                                 &
      kfac_lut = kfac_lut * k_ground_ratio_max
    CALL table_build( lut_seg_tab(pft,iseg), b, c, .FALSE., lut_f_floor,       &
                      lut_f_crit(1.0 - kcrit_fractional_loss(pft), kfac_lut) )
  END DO
  lut_seg_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut_seg)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE build_supply_lut_seg

! *****************************************************************************
! Lowest normalised conductance f = k/kmax the supply tables must reach:
! half the critical fraction of the path (kcrit/kmax = 1 -
! kcrit_fractional_loss, or sl_kfail_frac when set and lower), divided by
! kfac (a segment's kmax over the path's; 1 for the whole-plant table).
! sl_kfail_frac (stomata_model = 9) can be far below kcrit (Sperry et al.
! 2016 use 0.0005), and the tables used to end at a fixed fraction (1e-6
! Weibull, 0.1 kcrit/kmax SOX), so the failure point could lie below the
! bottom of the table.
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
! Grid of a supply table on the curve (b, c; SOX or Weibull): spacing dpsi
! (Pa) and node count n. As before, the table ends where f = f_floor, in
! n_lut - 1 cells. If f_need (lut_f_crit) is lower, the table is extended
! down to where f = f_need with the same spacing and more nodes, so the
! accuracy near psi = 0, where most of the flow is, is unchanged and a
! table that needs no extension is exactly as before (bit-for-bit). Only
! past n_lut_max nodes (a 16-fold longer table, e.g. a SOX curve with c ~ 1
! and a tiny sl_kfail_frac) is the spacing widened instead.
!
! MDK, 20/10/26
! *****************************************************************************
SUBROUTINE lut_grid( b, c, l_sox, f_floor, f_need, dpsi, n )

REAL(KIND=real_jlslsm), INTENT(IN) :: b, c, f_floor, f_need
LOGICAL, INTENT(IN) :: l_sox
REAL(KIND=real_jlslsm), INTENT(OUT) :: dpsi
INTEGER, INTENT(OUT) :: n

REAL(KIND=real_jlslsm) :: psi_need, r

! Spacing as before (same expressions, so the same bits).
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
