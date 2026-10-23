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

! Everything is private unless listed here or marked PUBLIC where it is
! declared (the module variables ksr, ksr_path, kplant_scale, seg_kscale,
! k_ground_ratio_max and som_psi_in_min). The supply-table helpers
! (lut_s_at, build_supply_lut, build_supply_lut_seg, seg_lut_psi,
! lut_f_crit, lut_grid) are internal.
!
! MDK, 21/10/26
PRIVATE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='XYLEM_HYDRAULICS_CW_JLS_MOD'

PUBLIC :: leaf_conductance_CW_jls                                              &
,         xylem_conductance_CW_jls                                             &
,         leaf_psi_CW_jls                                                      &
,         incomplete_gamma                                                     &
,         incomplete_gamma_different_a                                         &
,         leaf_psi_lut_jls                                                     &
,         leaf_psi_capped_lut_jls                                              &
,         supply_lut_psi                                                       &
,         supply_lut_e_crit                                                    &
,         supply_lut_f                                                         &
,         leaf_psi_segments_jls                                                &
,         leaf_psi_segments_lut_jls                                            &
,         set_ksr_path                                                         &
,         xylem_f                                                              &
,         xylem_f_bc                                                           &
,         k_path_zero_flow                                                     &
,         set_k_ground_basis                                                   &
,         seg_kmax_fac                                                         &
,         lut_pft_curve                                                        &
,         init_supply_luts   ! MDK, 21/10/26

! Soil-to-root (rhizosphere) conductance Ksr for l_som_rhizo_series, a
! conductance (not the soil conductivity soil_k), kept separate from the
! plant xylem kmax so that anything that changes kmax (e.g. xylem
! impairment) leaves the soil term alone; the two meet only in series in
! the solvers.
!   ksr      : the PFT's Ksr summed over the layers, per m2 ground
!              (mol m-2 s-1 Pa-1, the units of kmax_pft * LAI). Set by
!              physiol for the current PFT just before sf_stom. < 0: no
!              soil link (fsmc_mod /= 2).
!   ksr_path : this leaf path's share of ksr, in the units of the kmax the
!              solver is given. Set by sf_stom (set_ksr_path) before each
!              solver call, from the path's share of the undamaged plant
!              conductance (share_sup). < 0: no soil link.
! Used by leaf_psi_jls (single conductance curve) and the segment solvers
! below.
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: ksr(:), ksr_path(:)

! Root and stem conductance per unit ground area (kplant_model = 1, with
! l_som_plant_segments). On the leaf-area basis of kmax_pft the root and
! stem segments have kmax_pft * seg_kfac_s * r, r = LAI_ref / LAI (capped
! at k_ground_ratio_max), and the leaf segment kmax_pft * seg_kfac_leaf, so
! the well-watered whole plant (segments in series) is
!   kmax_pft * g,  1 / g = sum_s 1 / (seg_kfac_s r_s)  (r_leaf = 1).
! sf_stom gives the optimiser kmax, its intact reference and kcrit times g
! (kplant_scale), so every use of the whole-plant kmax (cost, kcrit test,
! closed leaf_k) sees the plant actually present, and the segment solvers
! use seg_kfac_s * r_s / g (seg_kscale) of that kmax, which keeps the
! series sum of the segments equal to the kmax they are given (as seg_kfac
! does). Xylem impairment caps each segment at a fraction of its intact
! (scaled) conductance, so the stored damage does not change with LAI.
! Both are set per PFT by set_k_ground_basis before the leaf paths.
!   kplant_scale : g (land_pts).
!   seg_kscale   : r_s / g (land_pts, 3).
!
! MDK, 21/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: kplant_scale(:),               &
                                               seg_kscale(:,:)

REAL(KIND=real_jlslsm), PARAMETER, PUBLIC :: k_ground_ratio_max = 10.0
      ! Cap on LAI_ref / LAI (kplant_model = 1), also used as LAI -> 0.
      !
      ! MDK, 21/10/26

REAL(KIND=real_jlslsm), PARAMETER, PUBLIC :: som_psi_in_min = -1.0e9
      ! Floor on the root inlet potential (Pa) with l_som_rhizo_series when
      ! the soil conductance is (near) zero.

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
! solved by three table inversions (leaf_psi_segments_lut_jls).
TYPE(supply_table), ALLOCATABLE, SAVE, PRIVATE :: lut_seg_tab(:,:)
                            ! Table for (pft, segment).
                            !
                            ! MDK, 20/10/26
LOGICAL, ALLOCATABLE, SAVE, PRIVATE :: lut_seg_ready(:)

CONTAINS

! *****************************************************************************
! l_som_rhizo_series: set ksr_path, the soil-to-root conductance of the next
! leaf path, as share * ksr. share is the path's share of the plant (its
! kmax over kmax_pft * LAI, in the units of that kmax). It must come from
! the undamaged plant conductance, so that xylem impairment, which lowers
! the kmax given to the solver, does not lower the soil term.
! *****************************************************************************
SUBROUTINE set_ksr_path( land_pts, veg_pts, veg_index, share )

USE jules_vegetation_mod, ONLY: l_som_rhizo_series

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: share(land_pts)

INTEGER :: l, m

IF ( .NOT. l_som_rhizo_series ) RETURN
IF ( .NOT. ALLOCATED(ksr_path) ) ALLOCATE(ksr_path(land_pts))
ksr_path(:) = -1.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( ksr(l) >= 0.0 ) ksr_path(l) = ksr(l) * share(l)
END DO

END SUBROUTINE set_ksr_path

! *****************************************************************************
! kplant_model = 1: set kplant_scale and seg_kscale (see their declaration)
! for PFT pft from the current LAI and the reference LAI lai_ref (both m2
! leaf m-2 ground). Root and stem keep the conductance they have at
! lai_ref, per unit ground area; the leaf segment scales with LAI. With
! kplant_model = 0 both are 1 and are not used.
!
! MDK, 21/10/26
! *****************************************************************************
SUBROUTINE set_k_ground_basis( pft, land_pts, veg_pts, veg_index, lai,        &
                               lai_ref )

USE pftparm, ONLY: seg_kfac
USE jules_vegetation_mod, ONLY: kplant_model, kplant_model_ground

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), lai_ref(land_pts)

INTEGER :: l, m, iseg
REAL(KIND=real_jlslsm) :: r(3), rinv

IF ( .NOT. ALLOCATED(kplant_scale) ) THEN
  ALLOCATE(kplant_scale(land_pts))
  ALLOCATE(seg_kscale(land_pts,3))
END IF
kplant_scale(:) = 1.0
seg_kscale(:,:) = 1.0
IF ( kplant_model /= kplant_model_ground ) RETURN

DO m = 1,veg_pts
  l = veg_index(m)
  ! Root and stem: LAI_ref / LAI, capped (and the cap as LAI -> 0).
  IF ( lai(l) * k_ground_ratio_max > lai_ref(l) ) THEN
    r(1:2) = lai_ref(l) / lai(l)
  ELSE
    r(1:2) = k_ground_ratio_max
  END IF
  ! and floored at 1 / k_ground_ratio_max, so LAI_ref = 0 (no height) or
  ! LAI >> LAI_ref cannot give a zero segment (NaN shares).
  !
  ! MDK, 22/10/26
  r(1:2) = MAX(r(1:2), 1.0 / k_ground_ratio_max)
  r(3) = 1.0
  ! Whole-plant (series) conductance relative to kmax: 1 / sum_s 1/k_s.
  rinv = 0.0
  DO iseg = 1, 3
    rinv = rinv + 1.0 / (seg_kfac(pft,iseg) * r(iseg))
  END DO
  kplant_scale(l) = 1.0 / rinv
  DO iseg = 1, 3
    seg_kscale(l,iseg) = r(iseg) * rinv
  END DO
END DO

END SUBROUTINE set_k_ground_basis

! *****************************************************************************
! Maximum conductance of segment iseg as a multiple of the whole-plant kmax
! given to the solvers at point l: seg_kfac, times seg_kscale with
! kplant_model = 1 (see set_k_ground_basis).
!
! MDK, 21/10/26
! *****************************************************************************
FUNCTION seg_kmax_fac( pft, l, iseg ) RESULT( fac )

USE pftparm, ONLY: seg_kfac
USE jules_vegetation_mod, ONLY: kplant_model, kplant_model_ground

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, l, iseg
REAL(KIND=real_jlslsm) :: fac

fac = seg_kfac(pft,iseg)
IF ( kplant_model == kplant_model_ground ) fac = fac * seg_kscale(l,iseg)

END FUNCTION seg_kmax_fac

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
! .TRUE. if every open point's curve (b, c) is the PFT's, so the supply
! table (keyed by PFT) holds for all of them. A point with its own curve
! (e.g. impaired b and c) must be solved by Newton instead; the table
! solvers used to ignore b and c.
!
! MDK, 21/10/26
! *****************************************************************************
LOGICAL FUNCTION lut_pft_curve( pft, land_pts, open_pnts, veg_index,          &
                                open_index, b, c )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, land_pts, open_pnts, veg_index(land_pts),         &
                       open_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: b(land_pts), c(land_pts)
INTEGER :: j, l

lut_pft_curve = .TRUE.
DO j = 1, open_pnts
  l = veg_index(open_index(j))
  IF ( b(l) /= conductance_b_pft(pft) .OR.                                    &
       c(l) /= conductance_c_pft(pft) ) THEN
    lut_pft_curve = .FALSE.
    RETURN
  END IF
END DO

END FUNCTION lut_pft_curve

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

! *****************************************************************************
! Whole-path conductance -dE/dpsi_leaf (in the units of kmax) when the
! plant is not transpiring. With no flow there is no drop in water
! potential along the path, so every element sits at psi_root and the
! elements add in series:
!   1 / K = 1 / ksr_path + sum over segments of 1 / k_s(psi_root)
! This is what the solvers give at E = 0, and is used for points the
! optimiser does not visit (closed stomata).
!
! The elements:
! - Segments (l_som_plant_segments): root, stem and leaf, each with maximum
!   conductance kmax * seg_kmax_fac (seg_kfac, times the LAI scaling
!   seg_kscale with kplant_model = 1) and its own curve. Without segments,
!   one element on the PFT curve, or on the curve b, c if given (a point's
!   impaired curve).
! - kcap_frac (xylem impairment, memory model with segments): caps the stem
!   and leaf segments, as leaf_psi_segments_jls; kmax is then the intact
!   value. With l_ximpair_seg_memory the leaf segment has its own cap.
! - kx: the plant (xylem) conductance at psi_root, if already known (e.g.
!   from an impairment model); replaces the curve.
! - The soil-to-root link ksr_path (l_som_rhizo_series, ksr_path >= 0),
!   which must be set by set_ksr_path for this path. Left out when
!   l_plant_only is .TRUE. (e.g. CGain's k_max).
!
! MDK, 21/10/26
! *****************************************************************************
FUNCTION k_path_zero_flow( pft, l, kmax, psi_root, b, c, kcap_frac, kx,       &
                           l_plant_only ) RESULT( k0 )

USE pftparm, ONLY: conductance_b_seg, conductance_c_seg
USE jules_vegetation_mod, ONLY: l_som_plant_segments, l_som_rhizo_series
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, l
                            ! Plant functional type and land point.
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax, psi_root
                            ! Whole-plant maximum conductance
                            ! (mol m-2 s-1 Pa-1) and the root zone water
                            ! potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: b, c, kcap_frac, kx
                            ! Curve parameters b (Pa) and c
                            ! (dimensionless); the segment cap fraction
                            ! (dimensionless); the known plant conductance
                            ! at psi_root (mol m-2 s-1 Pa-1).
LOGICAL, INTENT(IN), OPTIONAL :: l_plant_only
                            ! .TRUE.: leave out the soil-to-root link.
REAL(KIND=real_jlslsm) :: k0, r, k_el, k_s
                            ! The result (mol m-2 s-1 Pa-1); the series
                            ! resistance sum; an element's and the soil's
                            ! conductance.
                            !
                            ! MDK, 21/10/26
INTEGER :: iseg
LOGICAL :: l_soil

l_soil = l_som_rhizo_series
IF ( PRESENT(l_plant_only) ) l_soil = l_soil .AND. .NOT. l_plant_only

k0 = 0.0
IF ( kmax <= 0.0 ) RETURN

IF ( PRESENT(kx) ) THEN
  IF ( kx <= TINY(1.0_real_jlslsm) ) RETURN
  r = 1.0 / kx
ELSE IF ( l_som_plant_segments ) THEN
  r = 0.0
  DO iseg = 1, 3
    k_el = kmax * seg_kmax_fac(pft,l,iseg)                                     &
           * EXP(-( ABS(MIN(psi_root, 0.0) / conductance_b_seg(pft,iseg)) )    &
                 **conductance_c_seg(pft,iseg))
    IF ( PRESENT(kcap_frac) .AND. iseg >= 2 )                                  &
      k_el = MIN(k_el, MIN(MAX(ximpair_seg_frac(pft, l, iseg, kcap_frac),      &
                               0.0), 1.0) * kmax * seg_kmax_fac(pft,l,iseg))
    IF ( k_el <= TINY(1.0_real_jlslsm) ) RETURN
    r = r + 1.0 / k_el
  END DO
ELSE
  IF ( PRESENT(b) .AND. PRESENT(c) ) THEN
    k_el = kmax * xylem_f_bc(pft, psi_root, b, c)
  ELSE
    k_el = kmax * xylem_f(pft, psi_root)
  END IF
  IF ( k_el <= TINY(1.0_real_jlslsm) ) RETURN
  r = 1.0 / k_el
END IF

IF ( l_soil ) THEN
  IF ( ksr_path(l) >= 0.0 ) THEN
    k_s = ksr_path(l)
    IF ( k_s <= TINY(1.0_real_jlslsm) ) RETURN
    r = r + 1.0 / k_s
  END IF
END IF

k0 = 1.0 / r

END FUNCTION k_path_zero_flow

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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

xylem_conductance = SPREAD(kmax, DIM = 1, NCOPIES = n_water_potentials)        &
     * EXP(-ABS(water_potential                                                &
                / SPREAD(conductance_b, DIM = 1, NCOPIES = n_water_potentials))&
                ** SPREAD(conductance_c, DIM = 1, NCOPIES = n_water_potentials))

! Apply minimum conductance limit
xylem_conductance = MAX(xylem_conductance,                                     &
                        SPREAD(kcrit, DIM = 1, NCOPIES = n_water_potentials))

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
!      leaf_psi_lut_jls).
! With l_som_plant_segments, methods 2 and 3 hand over to the segment
! solvers (leaf_psi_segments_jls, leaf_psi_segments_lut_jls); method 1 with
! segments is rejected by the namelist check.
!
! NOTE: This subroutine is designed for the stomatal optimisation model.
!       It returns an output array that only contains values for points
!       with open stomata.
!
! MDK, 21/10/26
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
                            conductance_b,                                     &
                            conductance_c,                                     &
                         ! INTENT OUT
                            leaf_psi,                                          &
                            leaf_k                                             &
  )

USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_taylor,            &
                                psi_solver_newton, psi_solver_lut,             &
                                l_som_plant_segments

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
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_CW_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! The table is keyed by PFT, so a point with its own curve (b, c) is solved
! by Newton. (The segment solvers take their curves from the PFT's segment
! parameters, not b and c, so they keep the table.)
!
! MDK, 21/10/26
psi_solver = som_psi_solver
IF ( psi_solver == psi_solver_lut .AND. .NOT. l_som_plant_segments ) THEN
  IF ( .NOT. lut_pft_curve( pft, land_pts, open_pnts, veg_index, open_index,  &
                            conductance_b, conductance_c ) )                  &
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
! MDK, 21/10/26
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

        ! Calculate the new leaf water potential.
        leaf_psi(i,j) = root_zone_psi(l) - e_leaf(i,j) / leaf_k(i,j)

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
  CALL xylem_conductance_CW_jls(pft, n_e_leaf, open_pnts, leaf_psi,           &
                                kmax_open, kcrit_open, b_open, c_open,        &
                                leaf_k)

  leaf_k(:,:) = MAX(leaf_k(:,:),                                              &
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
! psi is held at <= 0 in (psi/b)**c, as for the Taylor solver.
!
! MDK, 21/10/26
! ---------------------------------------------------------------------
CASE(psi_solver_newton)
  ! Root / stem / leaf segments in series (l_som_plant_segments): leaf psi
  ! from solving the segments downstream, and the whole-plant conductance
  ! -dE/dpsi_leaf returned in leaf_k. The single-segment code below is left
  ! unchanged for the default case.
  IF ( l_som_plant_segments ) THEN
    CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pnts,           &
                                veg_index, open_index, e_leaf, root_zone_psi,  &
                                kmax, kcrit, leaf_psi, leaf_k )
  ELSE

    ! Calculate the lower incomplete gamma function for the root zone water
    ! potential, for all land points.
    gamma_root_psi = incomplete_gamma_different_a(land_pts, 1/conductance_c,  &
                       (MIN(root_zone_psi,0.0_real_jlslsm)/conductance_b)      &
                       **conductance_c)

    DO j = 1, open_pnts ! Iterate over the open points.
      l = veg_index(open_index(j))

      ! Do the first iteration outside the loop.

      ! Calculate the xylem conductance at the root zone water potential.
      leaf_k(:,j) = kmax(l)                                                    &
          * EXP( -(MIN(root_zone_psi(l),0.0_real_jlslsm)/conductance_b(l))     &
                 **conductance_c(l) )

      ! The initial guess for the leaf water potential is that of the root
      ! zone meaning that the transpiration rate is zero. This simplifies the
      ! first prediction of the leaf water potential to;
      ! (leaf_k floored as in the update below: in very dry soil it
      ! underflows to 0, and the first guess was then -Inf and the sample
      ! NaN. Only that case changes.)
      !
      ! MDK, 21/10/26
      leaf_psi(:,j) = root_zone_psi(l) - e_leaf(:,j)                           &
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
      DO i = 1, max_nr_iter

        ! Calculate the xylem conductance for the current leaf water
        ! potential.
        leaf_k_prev(:,j) = leaf_k(:,j)
        leaf_k(:,j) = kmax(l)                                                  &
            * EXP( -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))      &
                   **conductance_c(l) )

        ! Calculate the transpiration rate for the current leaf water
        ! potential. First get the incomplete gamma function for the leaf
        ! water potential.
        e_leaf_current(:,j) = incomplete_gamma(n_e_leaf,                       &
                                1/conductance_c(l),                            &
                                (MIN(leaf_psi(:,j),0.0_real_jlslsm)            &
                                 /conductance_b(l))**conductance_c(l))

        ! Calculate the difference between the incomplete gamma function at
        ! the leaf and root zone water potentials.
        e_leaf_current(:,j) = e_leaf_current(:,j) - gamma_root_psi(l)

        ! Multiply by kmax * (-b/c) to get the transpiration rate
        ! Note: The negative sign is present because conductance_b is
        ! negative.
        e_leaf_current(:,j) = e_leaf_current(:,j) * kmax(l)                    &
            * (-conductance_b(l)/conductance_c(l))

        ! Calculate the new leaf water potential. leaf_k is floored in the
        ! denominator to avoid a literal divide-by-zero for a sample that has
        ! already collapsed leaf_k to (numerically) 0; MAX(...) here does not
        ! change the update for any sample that is still within reach of a
        ! finite solution.
        leaf_psi(:,j) = leaf_psi(:,j)                                          &
               - (e_leaf(:,j) - e_leaf_current(:,j))                           &
                 / MAX(leaf_k(:,j), TINY(1.0_real_jlslsm))

        ! Guard against runaway steps for an infeasible sample: once leaf_k
        ! has collapsed there is no need to keep tracking the true (very
        ! negative) asymptotic solution, only to land clearly below kcrit
        ! without risking overflow in (leaf_psi/b)**c on the next pass.
        leaf_psi(:,j) = MAX(leaf_psi(:,j),                                     &
                            root_zone_psi(l) - 5.0 * ABS(conductance_b(l)))

        IF (MAXVAL(ABS(leaf_k(:,j) - leaf_k_prev(:,j)))                        &
            < k_conversion_limit) EXIT

      END DO

      ! Calculate the xylem conductance for the current leaf water potential.
      leaf_k(:,j) = kmax(l)                                                    &
            * EXP( -(MIN(leaf_psi(:,j),0.0_real_jlslsm)/conductance_b(l))      &
                   **conductance_c(l) )

    END DO

  END IF  ! l_som_plant_segments

! ---------------------------------------------------------------------
! Invert the tabulated supply function; see leaf_psi_lut_jls.
! ---------------------------------------------------------------------
CASE(psi_solver_lut)
  IF ( l_som_plant_segments ) THEN
    CALL leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pnts,       &
                                    veg_index, open_index, e_leaf,            &
                                    root_zone_psi, kmax, leaf_psi, leaf_k )
  ELSE
    CALL leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index,     &
                           open_index, e_leaf, root_zone_psi, kmax,           &
                           leaf_psi, leaf_k )
  END IF  ! l_som_plant_segments

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
        'som_psi_solver should be psi_solver_taylor, psi_solver_newton ' //    &
        'or psi_solver_lut')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_psi_CW_jls


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
SUBROUTINE leaf_psi_lut_jls( pft, n_e_leaf, land_pts, open_pnts, veg_index, &
                                open_index, e_leaf, root_zone_psi, kmax,       &
                                leaf_psi, leaf_k )

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

! For each open point and each candidate transpiration rate: divide E by the
! point's kmax (the table is per unit kmax, so one table serves every point),
! invert the table for psi_leaf, then take k at that psi_leaf from the curve.
! The MAX guards a zero kmax.
!
! MDK, 18/10/26
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
! Memory xylem impairment model (pft_xylem_impairment_model = 3) without
! segments, from the supply-function table (psi_solver_lut): the capped
! curve k(psi) = MIN(kmax f(psi), kc) of leaf_psi_impaired_memory, inverted
! exactly from the intact PFT table. With psi_s where the intact curve
! meets the cap (f(psi_s) = kc/kmax), the flow from psi_r is carried at the
! constant conductance kc down to psi_s, e_cap = kc (psi_r - psi_s) if
! psi_r > psi_s (else 0), so
!   E <= e_cap:  psi_l = psi_r - E / kc                     (exact)
!   E >  e_cap:  S(psi_l) = S(psi_top) + (E - e_cap) / kmax, psi_top =
!                MIN(psi_r, psi_s)                          (one inversion)
! as the Newton-Raphson of leaf_psi_impaired_memory solves. leaf_k is the
! capped conductance at psi_l floored at kcrit, as there. The curve must be
! the PFT's (the table's): b and c are not changed by the memory model.
!
! MDK, 18/10/26
! *****************************************************************************
SUBROUTINE leaf_psi_capped_lut_jls( pft, n_e_leaf, land_pts, open_pnts,        &
                                    veg_index, open_index, e_leaf,             &
                                    root_zone_psi, kmax_ref, kcap, kcrit,      &
                                    leaf_psi, leaf_k )

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

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

END SUBROUTINE leaf_psi_capped_lut_jls

! *****************************************************************************
! Leaf water potential (Pa) from the supply-function table, for root zone
! potential psi_root (Pa) and normalised demand e_over_kmax = e_leaf/kmax (Pa):
! solves S(psi_l) = S(psi_r) + e_over_kmax. See leaf_psi_lut_jls.
! *****************************************************************************
FUNCTION supply_lut_psi( pft, psi_root, e_over_kmax ) RESULT( psi_l )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, e_over_kmax
REAL(KIND=real_jlslsm) :: psi_l

INTEGER :: lo, hi, mid, i_root, n
REAL(KIND=real_jlslsm) :: dpsi, r_root, s_target

CALL build_supply_lut(pft)
! The table's spacing and node count (n is n_lut unless the table was
! extended; see lut_grid).
!
! MDK, 20/10/26
dpsi = lut_tab(pft)%dpsi
n = lut_tab(pft)%n

! S at the root zone water potential plus the demand, and the table cell
! psi_root is in (the search starts there).
! s_target is the S the leaf must reach: E = kmax * (S(psi_l) - S(psi_r)),
! rearranged. r_root is psi_root's fractional position in the table and
! i_root the upper (wetter) node of its cell, clamped as in lut_s_at.
!
! MDK, 18/10/26
r_root = MIN(MAX(-psi_root / dpsi, 0.0), REAL(n - 1))
i_root = MIN(INT(r_root) + 1, n - 1)
s_target = lut_s_at(pft, psi_root) + e_over_kmax

! Demand at or beyond all the table can supply: put the leaf at the bottom
! of the table (or at psi_root, if that is lower). k there is at most half
! of kcrit (and of sl_kfail_frac * kmax, when set), so the sample is
! rejected (stom_opt_mod) or counted as failed (supply-loss).
!
! MDK, 20/10/26
IF ( s_target >= lut_tab(pft)%s(n) ) THEN
  psi_l = MIN(-(n - 1) * dpsi, psi_root)
ELSE
  ! S(lo) <= s_target < S(hi); S is strictly increasing.
  ! Binary search for that cell, starting from psi_root's cell (the leaf is
  ! never wetter than the root zone). S rises with depth, so the cell is
  ! unique; the search takes about log2(n) = 12 steps (16 for n_lut_max).
  !
  ! MDK, 18/10/26
  lo = i_root
  hi = n
  DO WHILE ( hi - lo > 1 )
    mid = (lo + hi) / 2
    IF ( lut_tab(pft)%s(mid) <= s_target ) THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO
  ! Linear interpolation of psi within the cell (the inverse of lut_s_at):
  ! node lo is at psi = -(lo-1)*dpsi, and the fraction of the cell's S rise
  ! that s_target covers sets how far below it psi_l lies. The cap at
  ! psi_root stops rounding putting the leaf above the root zone.
  !
  ! MDK, 18/10/26
  psi_l = -(lo - 1) * dpsi - dpsi * (s_target - lut_tab(pft)%s(lo))          &
          / (lut_tab(pft)%s(hi) - lut_tab(pft)%s(lo))
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

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root, kcrit_frac
REAL(KIND=real_jlslsm) :: e_crit

REAL(KIND=real_jlslsm) :: psi_crit, r

CALL build_supply_lut(pft)
! f(psi_crit) = kcrit_frac, kept inside (0, 1) (e.g. kcrit >= kmax would
! otherwise give a NaN).
r = MIN(MAX(kcrit_frac, 1.0e-6_real_jlslsm),                                  &
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
e_crit = MAX(lut_s_at(pft, psi_crit) - lut_s_at(pft, psi_root), 0.0)

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

! S(psi) by linear interpolation in the table (the table must be built).
FUNCTION lut_s_at( pft, psi ) RESULT( s_psi )

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: s_psi

INTEGER :: i0, n
REAL(KIND=real_jlslsm) :: r

! r is psi's fractional position in the table: 0 at psi = 0, n - 1 at
! the bottom. It is clamped, so psi > 0 gives S = 0 and psi below the table
! gives S at the bottom (no extrapolation). i0 is the upper (wetter) node of
! psi's cell, and S is interpolated linearly between nodes i0 and i0+1.
! n is the table's node count.
!
! MDK, 20/10/26
n = lut_tab(pft)%n
r = MIN(MAX(-psi / lut_tab(pft)%dpsi, 0.0), REAL(n - 1))
i0 = MIN(INT(r) + 1, n - 1)
s_psi = lut_tab(pft)%s(i0)                                                     &
        + (r - (i0 - 1)) * (lut_tab(pft)%s(i0+1) - lut_tab(pft)%s(i0))

END FUNCTION lut_s_at

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
! lut_tab(pft)%s(i) = integral(f, psi_i, 0) by Simpson's rule on each grid
! cell, for the PFT's conductance model (see lut_tab at the top of this
! module). Built by numerical
! integration rather than the closed forms (incomplete gamma /
! hypergeometric series), whose series are poorly conditioned deep into the
! curve.
! *****************************************************************************
SUBROUTINE build_supply_lut( pft )

USE max_dimensions, ONLY: npft_max
USE pftparm, ONLY: conductance_b_pft, conductance_c_pft,                       &
                   pft_conductance_model, kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: SOX_conductance

INTEGER, INTENT(IN) :: pft

INTEGER :: i, n
REAL(KIND=real_jlslsm) :: b, c, dpsi, p0, pm, p1, f_floor
LOGICAL :: l_sox

! Fast exit when this PFT's table is already built. Every supply_lut_*
! call passes through here, so after the first call this is the only work
! done. The check is outside the critical section so that it costs no
! locking.
!
! MDK, 18/10/26
IF ( ALLOCATED(lut_ready) ) THEN
  IF ( lut_ready(pft) ) RETURN
END IF

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
  CALL lut_grid( b, c, l_sox, f_floor,                                         &
                 lut_f_crit(1.0 - kcrit_fractional_loss(pft), 1.0_real_jlslsm),&
                 dpsi, n )
  lut_tab(pft)%dpsi = dpsi
  lut_tab(pft)%n = n
  ALLOCATE( lut_tab(pft)%s(n) )
  ! Accumulate S downwards from S(0) = 0: each cell [p1, p0] adds the
  ! integral of f over it by Simpson's rule (ends and midpoint weighted
  ! 1:4:1, times dpsi/6). f > 0, so S rises strictly with depth, which the
  ! binary search in supply_lut_psi relies on.
  !
  ! MDK, 18/10/26
  lut_tab(pft)%s(1) = 0.0
  DO i = 2, n
    p0 = -(i - 2) * dpsi
    p1 = -(i - 1) * dpsi
    pm = 0.5 * (p0 + p1)
    lut_tab(pft)%s(i) = lut_tab(pft)%s(i-1) + dpsi / 6.0                       &
                        * ( f_of(p0) + 4.0 * f_of(pm) + f_of(p1) )
  END DO
  ! Mark the table ready only once it is complete, so the fast exit above
  ! never returns a half-built table.
  !
  ! MDK, 18/10/26
  lut_sox(pft) = l_sox
  lut_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut)

CONTAINS

  ! f = k/kmax at psi = p on this PFT's curve (as supply_lut_f), using the
  ! b, c and l_sox set above.
  !
  ! MDK, 18/10/26
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
! gs_opt.c (setup_plant / plant_supply).
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
                                  leaf_psi, leaf_k,                            &
                                  kcap_frac, leaf_k_intact, stem_psi,          &
                                  l_exact )

! Xylem impairment (memory model, pft_xylem_impairment_model = 3) with
! segments: with kcap_frac present, the stem and leaf segments (which use
! the PFT curve the impairment model damages) are capped,
!   k_s(psi) = MIN( kmax_s f_s(psi), kcap_frac * kmax_s ),
! while the root segment stays intact (no root damage). With
! l_ximpair_seg_memory kcap_frac caps the stem only, and the leaf segment
! takes its own cap (ximpair_kcap_leaf). Each segment's flow
! is then E = A_s(psi_out) - A_s(psi_in), A_s the integral of the capped
! curve from psi to 0: kcap_s * (-psi) above psi_cap (where the intact curve
! meets the cap) and the intact integral below it. leaf_k_intact, if
! present, is the whole-plant conductance of the intact segments at the
! same (impaired) water potentials, for the hydraulic cost.

USE pftparm, ONLY: conductance_b_seg, conductance_c_seg
USE jules_vegetation_mod, ONLY: l_som_rhizo_series
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts),                                                 &
      ! Sampled transpiration (same units as kmax * psi).
  root_zone_psi(land_pts),                                                     &
      ! Root-zone water potential (Pa).
  kmax(land_pts),                                                              &
      ! Whole-plant maximum conductance (intact).
  kcrit(land_pts)
      ! Whole-plant critical conductance.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts),                                               &
      ! Leaf water potential (Pa).
  leaf_k(n_e_leaf, open_pnts)
      ! Whole-plant conductance -dE/dpsi_leaf.

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  kcap_frac(land_pts)
      ! Impaired / intact maximum conductance of the stem and leaf segments.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  leaf_k_intact(n_e_leaf, open_pnts),                                          &
      ! Whole-plant conductance of the intact segments at leaf_psi.
  stem_psi(n_e_leaf, open_pnts)
      ! Water potential at the outlet of the stem segment (Pa).
LOGICAL, INTENT(IN), OPTIONAL :: l_exact
      ! .TRUE.: each segment's outlet psi is refined by bisection on the
      ! supply integral after the Newton steps, so it is exact (to
      ! rounding) also near the most the segment can carry, where the few
      ! Newton steps may not converge; a flow above that most is
      ! infeasible (leaf_psi = som_psi_in_min, leaf_k = 0). Used for the
      ! cuticular loss, whose psi is used whatever its k.
      !
      ! MDK, 22/10/26

INTEGER, PARAMETER :: n_seg = 3, max_nr_iter = 4

REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12
      ! Floor on a segment's conductance, as a fraction of its kmax.

INTEGER :: j, l, iseg, it, i
LOGICAL :: l_closed(n_e_leaf)
      ! l_som_rhizo_series: no soil conductance and E > 0 (infeasible).
REAL(KIND=real_jlslsm) :: k_s
      ! l_som_rhizo_series: this path's soil-to-root conductance.

REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(n_e_leaf), psi_out(n_e_leaf),                                         &
      ! Segment inlet / outlet water potential (Pa).
  k_in(n_e_leaf), k_cur(n_e_leaf), k_prev(n_e_leaf),                           &
      ! Segment conductance at the inlet / current outlet / previous pass.
  a_in(n_e_leaf), e_cur(n_e_leaf),                                             &
      ! Supply integral at the inlet; flow for the current outlet.
  dpsi_de(n_e_leaf), dpsi_de_int(n_e_leaf),                                    &
      ! d(psi_out)/dE carried down the segments (as solved / intact curves).
  a_one(1),                                                                    &
      ! Supply integral at psi_root_zone (root inlet, same for all samples).
  kmx, kcr, bs, cs, fr, psi_cap, a_cap
      ! Segment kmax, kcrit, Weibull b (Pa) and c; cap fraction, the psi
      ! where the intact curve meets the cap and the intact integral there.
LOGICAL :: l_cap, l_impair, l_exact_on

l_impair = PRESENT(kcap_frac)
l_exact_on = .FALSE.
IF ( PRESENT(l_exact) ) l_exact_on = l_exact

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  psi_in(:)  = root_zone_psi(l)
  dpsi_de(:) = 0.0
  dpsi_de_int(:) = 0.0
  l_closed(:) = .FALSE.
  ! l_som_rhizo_series: the soil (rhizosphere) link in series ahead of the
  ! root segment. Linear in E: psi_root = psi_src - E / K_s, and the chain
  ! starts from dpsi_root/dE = -1/K_s, so leaf_k is the conductance of the
  ! whole soil-to-leaf path. The soil is not damaged, so the intact chain
  ! (leaf_k_intact, for the hydraulic cost under impairment: the intact
  ! segments at the same water potentials) starts the same way.
  IF ( l_som_rhizo_series .AND. ksr_path(l) >= 0.0 ) THEN
    k_s = ksr_path(l)
    IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
      psi_in(:)  = MAX(root_zone_psi(l) - e_leaf(:,j) / k_s,                   &
                       som_psi_in_min)
      dpsi_de(:) = -1.0 / k_s
      dpsi_de_int(:) = -1.0 / k_s
    ELSE
      l_closed(:) = e_leaf(:,j) > 0.0
    END IF
  END IF

  DO iseg = 1, n_seg
    kmx = kmax(l) * seg_kmax_fac(pft,l,iseg)   ! MDK, 21/10/26
    kcr = kcrit(l) * seg_kmax_fac(pft,l,iseg)
    bs  = conductance_b_seg(pft,iseg)
    cs  = conductance_c_seg(pft,iseg)

    ! Cap on the stem and leaf segments (not the root).
    fr = 1.0
    IF ( l_impair .AND. iseg >= 2 )                                            &
      fr = MIN(MAX(ximpair_seg_frac(pft, l, iseg, kcap_frac(l)), k_floor), 1.0)
    l_cap = fr < 1.0 - 1.0e-6
    IF ( l_cap ) THEN
      psi_cap = bs * (-LOG(fr))**(1.0 / cs)
      a_one(:) = incomplete_gamma(1, 1.0/cs, [(psi_cap/bs)**cs])
      a_cap = kmx * (-bs/cs) * a_one(1)
    ELSE
      psi_cap = 0.0
      a_cap = 0.0
    END IF

    k_in(:) = seg_k(psi_in)
    IF ( iseg == 1 .AND. .NOT. ( l_som_rhizo_series .AND.                     &
                                 ksr_path(l) >= 0.0 ) ) THEN
      ! The root inlet is psi_root_zone for every sample: one evaluation.
      ! (With l_som_rhizo_series it falls with E, so not here.)
      a_one(:) = seg_a(1, psi_in(1:1))
      a_in(:)  = a_one(1)
    ELSE
      a_in(:) = seg_a(n_e_leaf, psi_in)
    END IF

    ! First guess: the inlet conductance over the whole drop.
    psi_out(:) = psi_in(:) - e_leaf(:,j) / MAX(k_in(:), TINY(1.0_real_jlslsm))
    psi_out(:) = MAX(psi_out(:), psi_in(:) - 5.0 * ABS(bs))
    k_cur(:)   = k_in(:)

    DO it = 1, max_nr_iter
      k_prev(:) = k_cur(:)
      k_cur(:)  = seg_k(psi_out)
      e_cur(:)  = seg_a(n_e_leaf, psi_out) - a_in(:)
      psi_out(:) = psi_out(:) - (e_leaf(:,j) - e_cur(:))                       &
                                / MAX(k_cur(:), TINY(1.0_real_jlslsm))
      ! Same runaway guard as the single segment (infeasible demand).
      psi_out(:) = MAX(psi_out(:), psi_in(:) - 5.0 * ABS(bs))
      IF ( MAXVAL(ABS(k_cur(:) - k_prev(:))) < 0.1 * kcr ) EXIT
    END DO
    IF ( l_exact_on ) THEN
      DO i = 1, n_e_leaf
        CALL seg_bisect( i )
      END DO
    END IF

    ! Segment conductances are floored at a tiny fraction of the segment
    ! kmax (not TINY): a collapsed segment has k -> 0 by underflow, and with
    ! a TINY floor the recursion overflows to -Inf and then 0 * -Inf = NaN
    ! in the next segment, which previously came out as a huge whole-plant k
    ! and let an infeasible state (psi_leaf ~ -70 MPa) pass the k > kcrit
    ! mask. With this floor dpsi_de stays finite and k correctly ~ 0.
    k_cur(:)   = MAX(seg_k(psi_out), k_floor * kmx)
    dpsi_de(:) = ( MAX(k_in(:), k_floor * kmx) * dpsi_de(:) - 1.0 ) / k_cur(:)
    IF ( PRESENT(leaf_k_intact) ) THEN
      dpsi_de_int(:) = ( MAX(kmx * EXP( -(psi_in(:)/bs)**cs ), k_floor * kmx) &
                         * dpsi_de_int(:) - 1.0 )                              &
                       / MAX(kmx * EXP( -(psi_out(:)/bs)**cs ), k_floor * kmx)
    END IF
    psi_in(:)  = psi_out(:)
    IF ( iseg == 2 .AND. PRESENT(stem_psi) ) stem_psi(:,j) = psi_out(:)
  END DO

  leaf_psi(:,j) = psi_in(:)
  ! dpsi_de < 0 and finite here; anything else (defensive) is infeasible.
  WHERE ( dpsi_de(:) < 0.0 .AND. ABS(dpsi_de(:)) < HUGE(1.0_real_jlslsm) )
    leaf_k(:,j) = -1.0 / dpsi_de(:)
  ELSEWHERE
    leaf_k(:,j) = 0.0
  END WHERE
  IF ( PRESENT(leaf_k_intact) ) THEN
    WHERE ( dpsi_de_int(:) < 0.0 .AND.                                         &
            ABS(dpsi_de_int(:)) < HUGE(1.0_real_jlslsm) )
      leaf_k_intact(:,j) = -1.0 / dpsi_de_int(:)
    ELSEWHERE
      leaf_k_intact(:,j) = 0.0
    END WHERE
  END IF
  WHERE ( l_closed(:) )
    leaf_psi(:,j) = som_psi_in_min
    leaf_k(:,j)   = 0.0
  END WHERE
  IF ( PRESENT(leaf_k_intact) ) THEN
    WHERE ( l_closed(:) ) leaf_k_intact(:,j) = 0.0
  END IF
  IF ( PRESENT(stem_psi) ) THEN
    WHERE ( l_closed(:) ) stem_psi(:,j) = som_psi_in_min
  END IF
END DO

CONTAINS

  ! Segment conductance (capped where l_cap).
  FUNCTION seg_k( psi ) RESULT( k )
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n_e_leaf)
  REAL(KIND=real_jlslsm) :: k(n_e_leaf)
  k(:) = kmx * EXP( -(psi(:)/bs)**cs )
  IF ( l_cap ) k(:) = MIN(k(:), fr * kmx)
  END FUNCTION seg_k

  ! l_exact: the outlet psi of sample i by bisection, seg_a(psi_out) =
  ! a_in + E (seg_a rises as psi falls), between the inlet and the psi
  ! where (psi/b)**c = 12 (f = 6e-6: the curve has carried all it can, and
  ! the series in incomplete_gamma is still converged there; deeper, its
  ! 40 terms are not enough and seg_a is no longer monotone). A flow above
  ! that is infeasible.
  !
  ! MDK, 22/10/26
  SUBROUTINE seg_bisect( i )
  INTEGER, INTENT(IN) :: i
  INTEGER :: ib
  REAL(KIND=real_jlslsm) :: lo(1), hi(1), mid(1), a_t, a1(1)
  a_t   = a_in(i) + e_leaf(i,j)
  hi(1) = psi_in(i)
  lo(1) = -ABS(bs) * 12.0**(1.0 / cs)
  IF ( l_closed(i) ) RETURN
  IF ( lo(1) >= hi(1) ) THEN
    l_closed(i) = .TRUE.
    psi_out(i)  = lo(1)
    RETURN
  END IF
  a1 = seg_a(1, lo)
  IF ( a1(1) < a_t ) THEN
    l_closed(i) = .TRUE.
    psi_out(i)  = lo(1)
    RETURN
  END IF
  DO ib = 1, 60
    mid(1) = 0.5 * (lo(1) + hi(1))
    a1 = seg_a(1, mid)
    IF ( a1(1) < a_t ) THEN
      hi(1) = mid(1)
    ELSE
      lo(1) = mid(1)
    END IF
  END DO
  psi_out(i) = 0.5 * (lo(1) + hi(1))
  END SUBROUTINE seg_bisect

  ! Integral of the segment conductance from psi to 0 (so that the flow
  ! from psi_in to psi_out is seg_a(psi_out) - seg_a(psi_in)).
  FUNCTION seg_a( n, psi ) RESULT( a )
  INTEGER, INTENT(IN) :: n
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n)
  REAL(KIND=real_jlslsm) :: a(n)
  ! Note: the negative sign is present because bs is negative.
  a(:) = kmx * (-bs/cs) * incomplete_gamma(n, 1.0/cs, (psi(:)/bs)**cs)
  IF ( l_cap ) THEN
    WHERE ( psi(:) >= psi_cap )
      a(:) = fr * kmx * (-psi(:))
    ELSEWHERE
      a(:) = fr * kmx * (-psi_cap) + a(:) - a_cap
    END WHERE
  END IF
  END FUNCTION seg_a

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
!
! Xylem impairment (memory model with segments; kcap_frac present): as
! leaf_psi_segments_jls, the stem and leaf segments are capped at a fraction
! fr of their kmax_s (ximpair_seg_frac), the root is intact. The capped
! integral is
!   A_s(psi) = fr kmax_s (-psi)                                psi >= psi_cap
!            = fr kmax_s (-psi_cap) + kmax_s (S_s(psi) - S_s(psi_cap))  below,
! psi_cap where the intact curve meets the cap, f_s(psi_cap) = fr. So for an
! inlet above psi_cap the flow down to psi_cap, e_cap = fr kmax_s (psi_in -
! psi_cap), is carried at the constant conductance fr kmax_s: an E <= e_cap
! gives psi_out = psi_in - E / (fr kmax_s) exactly, and a larger E continues
! on the intact curve from psi_cap, S_s(psi_out) = S_s(psi_cap) + (E - e_cap)
! / kmax_s, one inversion of the unchanged table. An inlet below psi_cap is
! on the intact curve throughout. The caps are per-point numbers, so the
! tables need no rebuilding. leaf_k uses the capped segment conductances in
! the recursion; leaf_k_intact (if present) the intact ones at the same water
! potentials (the hydraulic cost), and stem_psi the stem outlet.
!
! MDK, 18/10/26
! *****************************************************************************
SUBROUTINE leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pnts,      &
                                      veg_index, open_index, e_leaf,           &
                                      root_zone_psi, kmax, leaf_psi, leaf_k,   &
                                      kcap_frac, leaf_k_intact, stem_psi )

USE pftparm, ONLY: conductance_b_seg, conductance_c_seg
USE jules_vegetation_mod, ONLY: l_som_rhizo_series
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts), root_zone_psi(land_pts), kmax(land_pts)
      ! As leaf_psi_segments_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts), leaf_k(n_e_leaf, open_pnts)
      ! Leaf water potential (Pa) and whole-plant -dE/dpsi_leaf.

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  kcap_frac(land_pts)
      ! Impaired / intact maximum conductance of the stem (and leaf) segment,
      ! as leaf_psi_segments_jls (kmax is then the intact value).
      !
      ! MDK, 18/10/26
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  leaf_k_intact(n_e_leaf, open_pnts),                                          &
      ! Whole-plant conductance of the intact segments at leaf_psi.
  stem_psi(n_e_leaf, open_pnts)
      ! Water potential at the outlet of the stem segment (Pa).
      !
      ! MDK, 18/10/26

INTEGER, PARAMETER :: n_seg = 3

REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12
      ! Floor on a segment's conductance, as a fraction of its kmax (as
      ! leaf_psi_segments_jls).

INTEGER :: i, j, l, iseg
LOGICAL :: l_cap, l_impair
REAL(KIND=real_jlslsm) :: kmx, bs, cs, psi_in, psi_out, k_in, k_out, dpsi_de, &
                          k_s, fr, psi_cap, e_cap, e_pos, dpsi_de_int
      ! fr, psi_cap: the segment's cap fraction and where the intact curve
      ! meets it; e_cap the flow the capped part above psi_cap carries;
      ! e_pos the flow (>= 0); dpsi_de_int d(psi)/dE on the intact segments.
      !
      ! MDK, 18/10/26

CALL build_supply_lut_seg(pft)
l_impair = PRESENT(kcap_frac)

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  DO i = 1, n_e_leaf
    psi_in  = root_zone_psi(l)
    dpsi_de = 0.0
    dpsi_de_int = 0.0
    ! l_som_rhizo_series: soil link ahead of the root segment (as
    ! leaf_psi_segments_jls).
    IF ( l_som_rhizo_series .AND. ksr_path(l) >= 0.0 ) THEN
      k_s = ksr_path(l)
      IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        psi_in  = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s, som_psi_in_min)
        dpsi_de = -1.0 / k_s
        dpsi_de_int = -1.0 / k_s
      ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
        leaf_psi(i,j) = som_psi_in_min
        leaf_k(i,j)   = 0.0
        ! MDK, 18/10/26
        IF ( PRESENT(leaf_k_intact) ) leaf_k_intact(i,j) = 0.0
        IF ( PRESENT(stem_psi) ) stem_psi(i,j) = som_psi_in_min
        CYCLE
      END IF
    END IF
    ! Root, stem then leaf in series. Each segment carries the same E, so its
    ! outlet psi comes from its own table with E divided by its own kmax, and
    ! becomes the next segment's inlet. dpsi_de carries dpsi/dE down the chain:
    ! from E = integral(k), k_out dpsi_out = k_in dpsi_in - dE.
    !
    ! MDK, 18/10/26
    DO iseg = 1, n_seg
      kmx = kmax(l) * seg_kmax_fac(pft,l,iseg)   ! MDK, 21/10/26
      bs  = conductance_b_seg(pft,iseg)
      cs  = conductance_c_seg(pft,iseg)
      ! Cap on the stem and leaf segments (not the root), as
      ! leaf_psi_segments_jls.
      !
      ! MDK, 18/10/26
      l_cap = .FALSE.
      fr = 1.0
      IF ( l_impair .AND. iseg >= 2 ) THEN
        fr = MIN(MAX(ximpair_seg_frac(pft, l, iseg, kcap_frac(l)), k_floor),  &
                 1.0)
        l_cap = fr < 1.0 - 1.0e-6
      END IF
      IF ( l_cap ) THEN
        ! Exact inversion of the capped integral (see above): the straight
        ! capped part down to psi_cap, then the intact table from psi_cap.
        !
        ! MDK, 18/10/26
        psi_cap = bs * (-LOG(fr))**(1.0 / cs)
        e_pos   = MAX(e_leaf(i,j), 0.0)
        IF ( psi_in > psi_cap ) THEN
          e_cap = fr * kmx * (psi_in - psi_cap)
          IF ( e_pos <= e_cap ) THEN
            psi_out = psi_in - e_pos / MAX(fr * kmx, TINY(1.0_real_jlslsm))
          ELSE
            psi_out = seg_lut_psi(pft, iseg, psi_cap,                          &
                                  (e_pos - e_cap)                              &
                                  / MAX(kmx, TINY(1.0_real_jlslsm)))
          END IF
        ELSE
          psi_out = seg_lut_psi(pft, iseg, psi_in,                             &
                                e_pos / MAX(kmx, TINY(1.0_real_jlslsm)))
        END IF
        k_in  = MAX(MIN(kmx * EXP( -(ABS(psi_in /bs))**cs ), fr * kmx),        &
                    k_floor * kmx)
        k_out = MAX(MIN(kmx * EXP( -(ABS(psi_out/bs))**cs ), fr * kmx),        &
                    k_floor * kmx)
      ELSE
        psi_out = seg_lut_psi(pft, iseg, psi_in,                               &
                              e_leaf(i,j) / MAX(kmx, TINY(1.0_real_jlslsm)))
        k_in  = MAX(kmx * EXP( -(ABS(psi_in /bs))**cs ), k_floor * kmx)
        k_out = MAX(kmx * EXP( -(ABS(psi_out/bs))**cs ), k_floor * kmx)
      END IF
      dpsi_de = ( k_in * dpsi_de - 1.0 ) / k_out
      ! The intact chain at the same water potentials (hydraulic cost).
      !
      ! MDK, 18/10/26
      IF ( PRESENT(leaf_k_intact) )                                            &
        dpsi_de_int = ( MAX(kmx * EXP( -(ABS(psi_in /bs))**cs ),              &
                            k_floor * kmx)                                     &
                        * dpsi_de_int - 1.0 )                                  &
                      / MAX(kmx * EXP( -(ABS(psi_out/bs))**cs ), k_floor * kmx)
      psi_in  = psi_out
      IF ( iseg == 2 .AND. PRESENT(stem_psi) ) stem_psi(i,j) = psi_out
    END DO
    ! The leaf outlet is the leaf psi. The whole-plant conductance is
    ! -dE/dpsi_leaf = -1/dpsi_de; it is zero if the chain gave no usable
    ! slope (e.g. a segment at its conductance floor).
    !
    ! MDK, 18/10/26
    leaf_psi(i,j) = psi_in
    IF ( dpsi_de < 0.0 .AND. ABS(dpsi_de) < HUGE(1.0_real_jlslsm) ) THEN
      leaf_k(i,j) = -1.0 / dpsi_de
    ELSE
      leaf_k(i,j) = 0.0
    END IF
    IF ( PRESENT(leaf_k_intact) ) THEN
      IF ( dpsi_de_int < 0.0 .AND.                                             &
           ABS(dpsi_de_int) < HUGE(1.0_real_jlslsm) ) THEN
        leaf_k_intact(i,j) = -1.0 / dpsi_de_int
      ELSE
        leaf_k_intact(i,j) = 0.0
      END IF
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

INTEGER :: lo, hi, mid, i0, n
REAL(KIND=real_jlslsm) :: dpsi, r, s_in, s_target

! S at the inlet by linear interpolation (as lut_s_at, for this segment's
! table), plus the normalised flow. A negative flow is treated as zero, so
! the outlet is never wetter than the inlet. n is the table's node count.
!
! MDK, 20/10/26
dpsi = lut_seg_tab(pft,iseg)%dpsi
n = lut_seg_tab(pft,iseg)%n
r  = MIN(MAX(-psi_in / dpsi, 0.0), REAL(n - 1))
i0 = MIN(INT(r) + 1, n - 1)
s_in = lut_seg_tab(pft,iseg)%s(i0)                                             &
       + (r - (i0 - 1)) * (lut_seg_tab(pft,iseg)%s(i0+1)                       &
                           - lut_seg_tab(pft,iseg)%s(i0))
s_target = s_in + MAX(e_over_k, 0.0)

! Then as supply_lut_psi: the bottom of the table if the flow is beyond all
! the segment can carry, otherwise a binary search from the inlet's cell and
! linear interpolation within the cell, capped at the inlet psi.
!
! MDK, 18/10/26
IF ( s_target >= lut_seg_tab(pft,iseg)%s(n) ) THEN
  psi_out = MIN(-(n - 1) * dpsi, psi_in)
ELSE
  lo = i0
  hi = n
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
  psi_out = MIN(psi_out, psi_in)
END IF

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
                                ! MDK, 21/10/26

INTEGER, INTENT(IN) :: pft

INTEGER :: i, iseg, n
REAL(KIND=real_jlslsm) :: b, c, dpsi, p0, pm, p1
REAL(KIND=real_jlslsm) :: kfac_lut   ! segment share the table is built for
                                     ! MDK, 21/10/26

! Same structure as build_supply_lut: fast exit once built, allocate for all
! PFTs inside the critical section on the first call, then build this PFT's
! three tables (root, stem, leaf). Segments are Weibull only, so the floor
! is always lut_f_floor and f is written out in the Simpson sum.
!
! MDK, 18/10/26
IF ( ALLOCATED(lut_seg_ready) ) THEN
  IF ( lut_seg_ready(pft) ) RETURN
END IF

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
    CALL lut_grid( b, c, .FALSE., lut_f_floor,                                 &
                   lut_f_crit(1.0 - kcrit_fractional_loss(pft), kfac_lut),     &
                   dpsi, n )
    lut_seg_tab(pft,iseg)%dpsi = dpsi
    lut_seg_tab(pft,iseg)%n = n
    ALLOCATE( lut_seg_tab(pft,iseg)%s(n) )
    lut_seg_tab(pft,iseg)%s(1) = 0.0
    DO i = 2, n
      p0 = -(i - 2) * dpsi
      p1 = -(i - 1) * dpsi
      pm = 0.5 * (p0 + p1)
      lut_seg_tab(pft,iseg)%s(i) = lut_seg_tab(pft,iseg)%s(i-1) + dpsi / 6.0   &
                              * ( EXP(-(ABS(p0 / b))**c)                       &
                                  + 4.0 * EXP(-(ABS(pm / b))**c)               &
                                  + EXP(-(ABS(p1 / b))**c) )
    END DO
  END DO
  lut_seg_ready(pft) = .TRUE.
END IF
!$OMP END CRITICAL (som_supply_lut_seg)

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

! ---------------------------------------------------------------------
! Function to calculate the lower incomplete gamma function. The
! lower incomplete gamma function is defined as;
!
! gamma(a, x) = integral( t^(a-1) * e^(-t) dt, 0, x)
!
! The code here is based on the implementation in Numerical Recipes
! in FORTRAN 90, page 211, and makes use of the series representation
! of the lower incomplete gamma function. JBaguley
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

! Local variables
REAL(KIND=real_jlslsm) :: step(n_x)
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
                            ! Covers psi down to about -10.5 MPa, past the
                            ! kcrit cut-off (about -8.2 MPa, 17-21 terms).
                            ! Samples that need more are rejected by the
                            ! k > kcrit mask anyway.
                            !
                            ! MDK, 21/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INCOMPLETE_GAMMA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial values (the n=0 term of the series)
step(:) = 1/a
gamma(:) = step

! Calculate series, iterating until every point's remaining term is
! negligible relative to its running sum (see max_terms and series_tol
! above).
!
! MDK, 21/10/26
DO i = 1, max_terms
  step = step * x / (a + i)
  gamma = gamma + step

  IF (ALL(ABS(step) < series_tol * ABS(gamma))) EXIT
END DO

! Multiply by terms outside the sum
gamma = gamma * EXP(-x) * x**a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma

! ---------------------------------------------------------------------
! The same as the above function but with an array input for a, for use
! with per-point (e.g. impaired) conductance curves. JBaguley
! ---------------------------------------------------------------------
FUNCTION incomplete_gamma_different_a( n_x, a, x ) RESULT( gamma )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_x ! Size of x array

REAL(KIND=real_jlslsm), INTENT(IN) :: a(n_x)
                            ! Shape parameter of each value, 1/c
                            ! (dimensionless).
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
INTEGER :: i

! Series convergence controls - see incomplete_gamma above.
REAL(KIND=real_jlslsm), PARAMETER :: series_tol = 1.0e-8
INTEGER, PARAMETER :: max_terms = 40

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='INCOMPLETE_GAMMA_DIFFERENT_A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set initial values (the n=0 term of the series)
step(:) = 1/a
gamma(:) = step

! Calculate series until converged - see incomplete_gamma above.
DO i = 1, max_terms
  step = step * x / (a + i)
  gamma = gamma + step

  IF (ALL(ABS(step) < series_tol * ABS(gamma))) EXIT
END DO

! Multiply by terms outside the sum
gamma = gamma * EXP(-x) * x**a

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END FUNCTION incomplete_gamma_different_a

END MODULE xylem_hydraulics_CW_jls_mod
