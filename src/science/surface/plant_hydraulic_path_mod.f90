! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE
MODULE plant_hydraulic_path_mod

USE um_types, ONLY: real_jlslsm
USE xylem_hydraulics_CW_jls_mod, ONLY: xylem_f, xylem_f_bc, incomplete_gamma_1
USE supply_table_mod, ONLY: seg_lut_psi, build_supply_lut_seg,                &
                            newton_psi_drop, seg_lut_psi_back, seg_lut_s
USE jules_vegetation_mod, ONLY: k_ground_ratio_max, l_som_seg_whole
USE jules_vegetation_mod, ONLY: l_som_plant_segments

IMPLICIT NONE

! The plant's water path from the soil to the leaf, around the
! vulnerability curves: the soil-to-root link (ksr, ksr_path), the root,
! stem and leaf segments in series and their solvers, the zero-flow
! conductance of the path, and the kplant_model = 1 scaling.
!
! MDK, 30/10/26
PRIVATE

CHARACTER(LEN=*),PARAMETER,PRIVATE :: ModuleName='PLANT_HYDRAULIC_PATH_MOD'

PUBLIC :: set_ksr_path                                                         &
,         set_store_q_path                                                     &
,         clear_store_q_path                                                   &
,         set_shared_stem                                                      &
,         clear_shared_stem                                                    &
,         set_k_ground_basis                                                   &
,         k_path_zero_flow                                                     &
,         psi_leaf_segments_jls                                                &
,         psi_leaf_segments_lut_jls                                            &
,         kcrit_of_path

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
! Used by psi_leaf_jls (single conductance curve), the segment solvers
! below, k_path_zero_flow, and som_path_psi and the zero-flow cost
! reference in stom_opt_mod.
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: ksr(:), ksr_path(:)

! Plant capacitance coupled to the optimiser (the capacitance store): the
! store's flow out at the stem node, held fixed for this leaf path while
! the optimiser searches (q_f, > 0: the store gives; < 0: it refills).
! With it set, the soil link, the root and the stem segments carry
!
!   E_r = E - q_f
!
! (not floored at 0, see psi_leaf_segments_jls) and the leaf segment
! carries E, so a full store leaves the stem and leaf wetter at the same
! E, and the hydraulic cost (k = -dE/dpsi_leaf) is that of the xylem
! alone: q_f does not move with E, so the store adds no conductance of its
! own. The cost is measured against the unshifted path at zero flow
! (stom_opt_mod, kl_hc_max). With q_f = 0 the chain is the steady-state
! one.
!   store_q_path   : this leaf path's share of q_f, in the units of the
!                    transpiration the solver is given. Set by sf_stom
!                    (set_store_q_path) before the optimiser's call.
!   l_store_q_path : .TRUE. while store_q_path applies (cleared after the
!                    call, so no other caller of the segment solver sees
!                    it).
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: store_q_path(:)
LOGICAL, PUBLIC :: l_store_q_path = .FALSE.

! Shared root and stem (l_som_shared_stem, two-leaf canopy): the sunlit and
! shaded leaf paths draw on one root and one stem. While a class is
! searched, its root and stem segments (and the soil-to-root link) are the
! whole plant's, and carry the whole canopy's flow,
!
!   E_r = E + E_o - q_f
!
! (E the class's transpiration, E_o the other class's, held fixed, and q_f
! the whole store's flow out with the store coupled), while its leaf
! segment keeps the class's own conductance and carries E. So the stem
! water potential the class sees is that of the shared stem.
!   shared_k_fac   : the whole plant's root and stem conductance over the
!                    class's (1 / its share of the plant), applied to the
!                    root and stem segments in seg_kmax_fac.
!   shared_e_other : E_o, in the units of the transpiration the solver is
!                    given.
!   l_shared_stem_path : .TRUE. while they apply (set by sf_stom with
!                    set_shared_stem before a class's search, cleared after
!                    it, as store_q_path).
! The hydraulic cost of either class is that of the whole plant at the
! whole canopy's flow: the shared root and stem, then one leaf segment
! with the sunlit and shaded leaves' conductance together, all carrying
! E + E_o (less q_f on the root and stem). The search solves that path
! with l_shared_whole_path set, which makes the leaf segment the whole
! plant's too (shared_k_fac on all three segments) and puts E_o on it.
!   l_shared_whole_path : .TRUE. while the whole-plant path is solved.
!   shared_kcrit_plant  : kcrit of the whole plant (kcrit_of_path of its
!                    intact conductance), in the units of the class's kmax,
!                    for the cost.
! The gain of a class may also be on the whole plant's basis
! (l_som_shared_gain): A over the sum of both classes' A at their
! hydraulic edges.
!   l_shared_gain_path : .TRUE. for that gain.
!   shared_amax_other  : the other class's A at its hydraulic edge (as
!                    max_al in stom_opt_bounded_search), held fixed.
!   shared_amax_self   : this class's, set by the search (0 where it does
!                    not search), for the other class's next search.
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: shared_k_fac(:),                &
                                               shared_e_other(:),              &
                                               shared_kcrit_plant(:),          &
                                               shared_amax_other(:),           &
                                               shared_amax_self(:)
LOGICAL, PUBLIC :: l_shared_stem_path = .FALSE.
LOGICAL, PUBLIC :: l_shared_whole_path = .FALSE.
LOGICAL, PUBLIC :: l_shared_gain_path = .FALSE.

! The whole-plant solve differs from the class's solve of the same E only
! in the leaf segment, so with one sample at one point the class's solve
! keeps the state at the stem's leaf-side end (water potential, dpsi/dE on
! the solved and on the intact segments, infeasible flag), with the inputs
! it came from, and the whole-plant solve of the same inputs starts from
! it at the leaf segment. The root and stem are then not solved twice;
! the result is the same to the bit. Cleared by set_shared_stem and
! clear_shared_stem.
LOGICAL :: l_rs_cache = .FALSE.
LOGICAL :: rs_intact, rs_exact, rs_closed
INTEGER :: rs_l
REAL(KIND=real_jlslsm) :: rs_e, rs_kmax, rs_kcrit, rs_psi0, rs_kcap,           &
                          rs_psi, rs_dpsi, rs_dpsi_int

! Root and stem conductance per unit ground area (kplant_model = 1, with
! l_som_plant_segments). On the leaf-area basis of kmax_pft the root and
! stem segments have kmax_pft * seg_kfac_s * r, r = LAI_ref / LAI (capped
! at k_ground_ratio_max), and the leaf segment kmax_pft * seg_kfac_leaf, so
! the well-watered whole plant (segments in series) is
!
!   kmax_pft g,   1 / g = sum_s 1 / (seg_kfac_s r_s)   (r_leaf = 1).
!
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
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: kplant_scale(:),                &
                                               seg_kscale(:,:)

! .TRUE. while kplant_scale and seg_kscale act for the current PFT
! (kplant_model = 1, or the sapwood-area factor of l_nsc_k_sapwood), set by
! set_k_ground_basis.
LOGICAL, PUBLIC :: l_kscale = .FALSE.

! Floor on the root water potential (Pa) with l_som_rhizo_series when
! the soil conductance is (near) zero.
REAL(KIND=real_jlslsm), PARAMETER, PUBLIC :: som_psi_in_min = -1.0e9

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
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: share(land_pts)

INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_KSR_PATH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. l_som_rhizo_series ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
IF ( .NOT. ALLOCATED(ksr_path) ) ALLOCATE(ksr_path(land_pts))
ksr_path(:) = -1.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( ksr(l) >= 0.0 ) ksr_path(l) = ksr(l) * share(l)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_ksr_path

! *****************************************************************************
! the capacitance store: set store_q_path, the store's flow out held fixed
! for the next leaf path, as share * q (q per m2 ground, share as in
! set_ksr_path), and switch it on for the segment solver.
!
! MDK, 24/10/26
! *****************************************************************************
SUBROUTINE set_store_q_path( land_pts, veg_pts, veg_index, share, q )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: share(land_pts), q(land_pts)

INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_STORE_Q_PATH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ALLOCATED(store_q_path) ) ALLOCATE(store_q_path(land_pts))
store_q_path(:) = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  store_q_path(l) = q(l) * share(l)
END DO
l_store_q_path = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_store_q_path

! *****************************************************************************
! the capacitance store: switch store_q_path off again after the call it
! was set for.
!
! MDK, 24/10/26
! *****************************************************************************
SUBROUTINE clear_store_q_path()

IMPLICIT NONE

l_store_q_path = .FALSE.

END SUBROUTINE clear_store_q_path

! *****************************************************************************
! l_som_shared_stem: set the shared root and stem for the next leaf path
! (one class of the two-leaf canopy). share is the class's share of the
! undamaged plant conductance (as in set_ksr_path), so its root and stem
! segments are 1 / share times its own; e_other is the other class's
! transpiration (same units as the class's). A class with no share keeps
! its own root and stem. kcrit_plant is the whole plant's kcrit (units of
! the class's kmax), amax_other the other class's A at its hydraulic edge
! and l_gain the whole-plant gain (l_som_shared_gain).
!
! MDK, 30/10/26
! *****************************************************************************
SUBROUTINE set_shared_stem( land_pts, veg_pts, veg_index, share, e_other,      &
                            kcrit_plant, amax_other, l_gain )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: share(land_pts), e_other(land_pts),     &
                                      kcrit_plant(land_pts),                   &
                                      amax_other(land_pts)
LOGICAL, INTENT(IN) :: l_gain

INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_SHARED_STEM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ALLOCATED(shared_k_fac) ) THEN
  ALLOCATE(shared_k_fac(land_pts))
  ALLOCATE(shared_e_other(land_pts))
  ALLOCATE(shared_kcrit_plant(land_pts))
  ALLOCATE(shared_amax_other(land_pts))
  ALLOCATE(shared_amax_self(land_pts))
END IF
shared_k_fac(:)       = 1.0
shared_e_other(:)     = 0.0
shared_kcrit_plant(:) = 0.0
shared_amax_other(:)  = 0.0
shared_amax_self(:)   = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  shared_kcrit_plant(l) = kcrit_plant(l)
  shared_amax_other(l)  = MAX(amax_other(l), 0.0_real_jlslsm)
  IF ( share(l) > TINY(1.0_real_jlslsm) ) THEN
    shared_k_fac(l)   = 1.0 / share(l)
    shared_e_other(l) = MAX(e_other(l), 0.0_real_jlslsm)
  END IF
END DO
l_shared_stem_path  = .TRUE.
l_shared_whole_path = .FALSE.
l_rs_cache          = .FALSE.
l_shared_gain_path  = l_gain

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_shared_stem

! *****************************************************************************
! l_som_shared_stem: switch the shared root and stem off again after the
! call they were set for.
!
! MDK, 30/10/26
! *****************************************************************************
SUBROUTINE clear_shared_stem()

IMPLICIT NONE

l_shared_stem_path  = .FALSE.
l_shared_whole_path = .FALSE.
l_shared_gain_path  = .FALSE.
l_rs_cache          = .FALSE.

END SUBROUTINE clear_shared_stem

! *****************************************************************************
! kplant_model = 1: set kplant_scale and seg_kscale (see their declaration)
! for PFT pft from the current LAI and the reference LAI lai_ref (both m2
! leaf m-2 ground). Root and stem keep the conductance they have at
! lai_ref, per unit ground area; the leaf segment scales with LAI. With
! the sapwood-area factor k_area (l_nsc_k_sapwood) the stem segment's
! conductance is also scaled by it, or without segments the whole plant's.
! With kplant_model = 0 and no k_area both are 1 and are not used
! (l_kscale = .FALSE.).
!
! MDK, 31/10/26
! *****************************************************************************
SUBROUTINE set_k_ground_basis( pft, land_pts, veg_pts, veg_index, lai,         &
                               lai_ref, k_area )

USE pftparm, ONLY: seg_kfac
USE jules_vegetation_mod, ONLY: kplant_model, kplant_model_ground
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), lai_ref(land_pts)
! Sapwood-area factor on kmax (0-1; l_nsc_k_sapwood).
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: k_area(land_pts)

INTEGER :: l, m, iseg
REAL(KIND=real_jlslsm) :: r(3), rinv
! The plant has root, stem and leaf segments.
LOGICAL :: l_seg

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_K_GROUND_BASIS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ALLOCATED(kplant_scale) ) THEN
  ALLOCATE(kplant_scale(land_pts))
  ALLOCATE(seg_kscale(land_pts,3))
END IF
kplant_scale(:) = 1.0
seg_kscale(:,:) = 1.0
l_kscale = kplant_model == kplant_model_ground .OR. PRESENT(k_area)
IF ( .NOT. l_kscale ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
l_seg = l_som_plant_segments .OR. l_som_seg_whole

DO m = 1,veg_pts
  l = veg_index(m)
  r(:) = 1.0
  IF ( kplant_model == kplant_model_ground ) THEN
    ! Root and stem: LAI_ref / LAI, capped (and the cap as LAI -> 0).
    IF ( lai(l) * k_ground_ratio_max > lai_ref(l) ) THEN
      r(1:2) = lai_ref(l) / lai(l)
    ELSE
      r(1:2) = k_ground_ratio_max
    END IF
    ! and floored at 1 / k_ground_ratio_max, so LAI_ref = 0 (no height) or
    ! LAI >> LAI_ref cannot give a zero segment (NaN shares).
    r(1:2) = MAX(r(1:2), 1.0 / k_ground_ratio_max)
    r(3) = 1.0
  END IF

  ! Sapwood-area factor: on the stem segment, or without segments on the
  ! whole plant. A full area (1) with kplant_model = 0 leaves the point
  ! untouched (factors 1).
  IF ( PRESENT(k_area) ) THEN
    IF ( kplant_model /= kplant_model_ground .AND. k_area(l) >= 1.0 ) CYCLE
    IF ( .NOT. l_seg ) THEN
      kplant_scale(l) = k_area(l)
      CYCLE
    END IF
    r(2) = r(2) * k_area(l)
  END IF

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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_k_ground_basis

! *****************************************************************************
! Maximum conductance of segment iseg as a multiple of the whole-plant kmax
! given to the solvers at point l: seg_kfac, times seg_kscale with
! kplant_model = 1 or the sapwood-area factor (l_kscale; see
! set_k_ground_basis), and for the root and stem
! times shared_k_fac while the shared root and stem apply (and the leaf
! too while the whole-plant path is solved).
!
! MDK, 30/10/26
! *****************************************************************************
FUNCTION seg_kmax_fac( pft, l, iseg ) RESULT( fac )

USE pftparm, ONLY: seg_kfac

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, l, iseg
REAL(KIND=real_jlslsm) :: fac

fac = seg_kfac(pft,iseg)
IF ( l_kscale ) fac = fac * seg_kscale(l,iseg)
! The shared root and stem (l_som_shared_stem): the whole plant's.
IF ( l_shared_stem_path .AND. ( iseg <= 2 .OR. l_shared_whole_path ) )        &
  fac = fac * shared_k_fac(l)

END FUNCTION seg_kmax_fac

! *****************************************************************************
! Whole-path conductance -dE/dpsi_leaf (in the units of kmax) when the
! plant is not transpiring. With no flow there is no drop in water
! potential along the path, so every element sits at psi_root and the
! elements add in series:
!
!   1 / K = 1 / ksr_path + sum over segments of 1 / k_s(psi_root).
!
! This is what the solvers give at E = 0. It is used for points the
! optimiser does not visit (closed stomata), for the zero-flow reference
! of the hydraulic cost, for SOX's k(psi_bar) and for CGain's k_max
! (all in stom_opt_mod).
!
! The elements:
! - Segments (l_som_plant_segments): root, stem and leaf, each with maximum
!   conductance kmax * seg_kmax_fac (seg_kfac, times the LAI scaling
!   seg_kscale with kplant_model = 1) and its own curve. Without segments,
!   one element on the PFT curve, or on the curve b, c if given (a point's
!   impaired curve).
! - kcap_frac (xylem impairment, memory model with segments): caps the stem
!   and leaf segments, as psi_leaf_segments_jls; kmax is then the intact
!   value. With l_ximpair_seg_memory the leaf segment has its own cap.
! - kx: the plant (xylem) conductance at psi_root, if already known (e.g.
!   from an impairment model); replaces the curve.
! - The soil-to-root link ksr_path (l_som_rhizo_series, ksr_path >= 0),
!   which must be set by set_ksr_path for this path. Left out when
!   l_plant_only is .TRUE. (e.g. CGain's k_max).
!
! MDK, 26/10/26
! *****************************************************************************
FUNCTION k_path_zero_flow( pft, l, kmax, psi_root, b, c, kcap_frac, kx,        &
                           l_plant_only ) RESULT( k0 )

USE pftparm, ONLY: seg_rel_k
USE jules_vegetation_mod, ONLY: l_som_plant_segments, l_som_rhizo_series
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac

IMPLICIT NONE

! Plant functional type and land point.
INTEGER, INTENT(IN) :: pft, l
! Whole-plant maximum conductance (mol m-2 s-1 Pa-1) and the root zone water
! potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax, psi_root
! Curve parameters b (Pa) and c (dimensionless); the segment cap fraction
! (dimensionless); the known plant conductance at psi_root (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: b, c, kcap_frac, kx
! .TRUE.: leave out the soil-to-root link.
LOGICAL, INTENT(IN), OPTIONAL :: l_plant_only
! The result (mol m-2 s-1 Pa-1); the series resistance sum; an element's and the
! soil's conductance.
REAL(KIND=real_jlslsm) :: k0, r, k_el, k_s
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
    k_el = kmax * seg_kmax_fac(pft,l,iseg) * seg_rel_k(pft, iseg, psi_root)
    ! The cap on the stem and leaf (and the root with l_som_seg_whole).
    IF ( PRESENT(kcap_frac) .AND. ( iseg >= 2 .OR. l_som_seg_whole ) )         &
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

! *****************************************************************************
! Leaf water potential and whole-plant conductance for the plant as root,
! stem and leaf segments in series (l_som_plant_segments), as in GDAY
! gs_opt.c (setup_plant / plant_supply).
!
! Segment s has maximum conductance
!
!     kmax_s = kmax * seg_kfac(pft,s)
!
! (the segments share the whole-plant resistance, so in series they give
! kmax when well watered) and its own cumulative Weibull curve
! (conductance_b_seg, conductance_c_seg).
!
! For each sampled transpiration the leaf-side water potential of each
! segment is solved downstream from psi_root_zone with the same
! Newton-Raphson as the single segment (psi_solver_newton), and becomes the
! next segment's root-side end. The last leaf-side end is the leaf water
! potential.
!
! The conductance the hydraulic cost uses is the whole-plant conductance
! (Sperry et al. 2017)
!
!     k = -dE / dpsi_leaf
!
! Each segment carries the same flow
!
!     E = integral from psi_out to psi_in of k_s(psi) dpsi
!
! so differentiating with respect to E gives, segment by segment,
!
!     dpsi_out/dE = ( k_s(psi_in) * dpsi_in/dE - 1 ) / k_s(psi_out)
!
! starting from dpsi_root_zone/dE = 0 at the root end. With no flow the
! segments add as resistances in series:
!
!     1/k = sum over s of 1/k_s(psi)
!
! A single segment gives k = k(psi_leaf) (the default model).
!
! MDK, 30/10/26
! *****************************************************************************
SUBROUTINE psi_leaf_segments_jls( pft, n_e_leaf, land_pts, open_pnts,          &
                                  veg_index, open_index, e_leaf,               &
                                  root_zone_psi, kmax, kcrit,                  &
                                  leaf_psi, leaf_k,                            &
                                  kcap_frac, leaf_k_intact, stem_psi,          &
                                  l_exact )

! Xylem impairment (memory model, pft_xylem_impairment_model = 3) with
! segments: with kcap_frac present, the stem and leaf segments (which use
! the PFT curve the impairment model damages) are capped,
!
!   k_s(psi) = MIN( kmax_s f_s(psi), kcap_frac kmax_s ),
!
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
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Sampled transpiration (same units as kmax * psi).
  e_leaf(n_e_leaf, open_pnts),                                                 &
  ! Root-zone water potential (Pa).
  root_zone_psi(land_pts),                                                     &
  ! Whole-plant maximum conductance (intact).
  kmax(land_pts),                                                              &
  ! Whole-plant critical conductance.
  kcrit(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Leaf water potential (Pa).
  leaf_psi(n_e_leaf, open_pnts),                                               &
  ! Whole-plant conductance -dE/dpsi_leaf.
  leaf_k(n_e_leaf, open_pnts)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  ! Impaired / intact maximum conductance of the stem and leaf segments.
  kcap_frac(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  ! Whole-plant conductance of the intact segments at leaf_psi.
  leaf_k_intact(n_e_leaf, open_pnts),                                          &
  ! Water potential at the leaf-side end of the stem segment (Pa).
  stem_psi(n_e_leaf, open_pnts)
! .TRUE.: each segment's leaf-side psi is refined by bisection on the
! supply integral after the Newton steps, so it is exact (to
! rounding) also near the most the segment can carry, where the few
! Newton steps may not converge; a flow above that most is
! infeasible (leaf_psi = som_psi_in_min, leaf_k = 0). Used for the
! cuticular loss, whose psi is used whatever its k.
LOGICAL, INTENT(IN), OPTIONAL :: l_exact

INTEGER, PARAMETER :: n_seg = 3, max_nr_iter = 4

! Floor on a segment's conductance, as a fraction of its kmax.
REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12

INTEGER :: j, l, iseg, it, i
! l_som_rhizo_series: no soil conductance and E > 0 (infeasible).
LOGICAL :: l_closed(n_e_leaf)
! l_som_rhizo_series: this path's soil-to-root conductance.
REAL(KIND=real_jlslsm) :: k_s

REAL(KIND=real_jlslsm) ::                                                      &
  ! Segment water potential at the root-side / leaf-side end (Pa).
  psi_in(n_e_leaf), psi_out(n_e_leaf),                                         &
  ! Segment conductance at the root-side end / current leaf-side end /
  ! previous pass.
  k_in(n_e_leaf), k_cur(n_e_leaf), k_prev(n_e_leaf),                           &
  ! Supply integral at the root-side end; flow for the current leaf-side
  ! end.
  a_in(n_e_leaf), e_cur(n_e_leaf),                                             &
  ! d(psi_out)/dE carried down the segments (as solved / intact curves).
  dpsi_de(n_e_leaf), dpsi_de_int(n_e_leaf),                                    &
  ! Segment kmax, kcrit, Weibull b (Pa) and c; cap fraction, the psi
  ! where the intact curve meets the cap and the intact integral there.
  kmx, kcr, bs, cs, fr, psi_cap, a_cap
! Deepest Newton step below psi_in on the segment's curve (Pa,
! newton_psi_drop).
REAL(KIND=real_jlslsm) :: drop_nr
! The critical conductance of the path the segment belongs to, for
! drop_nr: kcrit, or the whole plant's for the shared root and stem.
REAL(KIND=real_jlslsm) :: kcr_drop
! The flow the current segment carries: the leaf's E, or with the
! store's flow held fixed (l_store_q_path) E_r = E - q_f for the soil
! link, the root and the stem (dE_r/dE = 1 either way).
REAL(KIND=real_jlslsm) :: e_seg(n_e_leaf)
LOGICAL :: l_cap, l_impair, l_exact_on
! The shared root and stem, one sample at one point: the stem state
! is kept by a class solve and reused by the whole-plant solve of the
! same inputs (see l_rs_cache).
LOGICAL :: l_one
INTEGER :: iseg0
REAL(KIND=real_jlslsm) :: kc_key
! Newton: samples still iterating (each stops on its own convergence).
LOGICAL :: l_nr_active(n_e_leaf)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_SEGMENTS_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_impair = PRESENT(kcap_frac)
l_exact_on = .FALSE.
IF ( PRESENT(l_exact) ) l_exact_on = l_exact
l_one = l_shared_stem_path .AND. n_e_leaf == 1 .AND. open_pnts == 1
IF ( l_one .AND. .NOT. l_shared_whole_path ) l_rs_cache = .FALSE.

! For each open point, walk every sampled transpiration down the chain of
! segments (root, then stem, then leaf). Each segment takes the water
! potential at its root-side end (psi_in), solves for the water potential
! at its leaf-side end (psi_out) that lets it carry E, and hands psi_out on
! as the next segment's psi_in. Alongside, dpsi_de carries how fast that
! water potential falls with E, which gives the whole-plant conductance at
! the end.
DO j = 1, open_pnts
  l = veg_index(open_index(j))

  ! Start of the chain: the root zone, which does not move with E.
  psi_in(:)  = root_zone_psi(l)
  dpsi_de(:) = 0.0
  dpsi_de_int(:) = 0.0
  l_closed(:) = .FALSE.

  ! The flow up to the stem node. With the store's flow held fixed
  ! (l_store_q_path) the soil, root and stem carry E_r = E - q_f. It is not
  ! floored at 0: below q_f it runs back towards the soil, so the supply
  ! curve is the steady one shifted by q_f, smooth in E. (A floor at 0
  ! would leave only the leaf segment below E = q_f, a conductance far
  ! above the whole plant's, so the hydraulic cost would be negative there
  ! and jump at q_f, pinning the optimum near q_f.) The store's balance
  ! (sf_stom, store_cpl_q) keeps E_r >= 0 in what is committed.
  e_seg(:) = e_leaf(:,j)
  IF ( l_store_q_path ) e_seg(:) = e_leaf(:,j) - store_q_path(l)

  ! The shared root and stem (l_shared_stem_path): the soil link, the root
  ! and the stem also carry the other class's E, held fixed, so
  !
  !   E_r = E + E_o - q_f
  !
  ! and dE_r/dE is still 1.
  IF ( l_shared_stem_path ) e_seg(:) = e_seg(:) + shared_e_other(l)

  ! l_som_rhizo_series: the soil (rhizosphere) link in series ahead of the
  ! root segment. Linear in E: psi_root = psi_src - E / K_s, and the chain
  ! starts from dpsi_root/dE = -1/K_s, so leaf_k is the conductance of the
  ! whole soil-to-leaf path. The soil is not damaged, so the intact chain
  ! (leaf_k_intact, for the hydraulic cost under impairment: the intact
  ! segments at the same water potentials) starts the same way.
  IF ( l_som_rhizo_series .AND. ksr_path(l) >= 0.0 ) THEN
    k_s = ksr_path(l)
    IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
      psi_in(:)  = MAX(root_zone_psi(l) - e_seg(:) / k_s,                      &
                       som_psi_in_min)
      dpsi_de(:) = -1.0 / k_s
      dpsi_de_int(:) = -1.0 / k_s
    ELSE
      l_closed(:) = e_seg(:) > 0.0
    END IF
  END IF

  ! The whole-plant solve of the same inputs as the last class solve:
  ! start at the leaf segment from the kept stem state (l_rs_cache).
  iseg0 = 1
  kc_key = -1.0
  IF ( l_impair ) kc_key = kcap_frac(l)
  IF ( l_one .AND. l_shared_whole_path .AND. l_rs_cache ) THEN
    IF ( rs_l == l .AND. rs_e == e_leaf(1,j) .AND. rs_kmax == kmax(l) .AND.   &
         rs_kcrit == kcrit(l) .AND. rs_psi0 == root_zone_psi(l) .AND.         &
         rs_kcap == kc_key .AND. ( rs_intact .EQV. PRESENT(leaf_k_intact) )   &
         .AND. ( rs_exact .EQV. l_exact_on ) ) THEN
      psi_in(1)      = rs_psi
      dpsi_de(1)     = rs_dpsi
      dpsi_de_int(1) = rs_dpsi_int
      l_closed(1)    = rs_closed
      IF ( PRESENT(stem_psi) ) stem_psi(1,j) = rs_psi
      iseg0 = n_seg
    END IF
  END IF

  DO iseg = iseg0, n_seg
    ! The leaf segment carries the leaf's E, whatever the store gives, and
    ! not the other class's unless the whole plant's leaf is solved.
    IF ( iseg == n_seg .AND. ( l_store_q_path .OR. l_shared_stem_path ) )     &
      e_seg(:) = e_leaf(:,j)
    IF ( iseg == n_seg .AND. l_shared_stem_path .AND. l_shared_whole_path )    &
      e_seg(:) = e_seg(:) + shared_e_other(l)

    ! This segment's maximum and critical conductance (its share of the
    ! whole-plant values) and its own vulnerability curve (b, c).
    kmx = kmax(l) * seg_kmax_fac(pft,l,iseg)
    kcr = kcrit(l) * seg_kmax_fac(pft,l,iseg)
    bs  = conductance_b_seg(pft,iseg)
    cs  = conductance_c_seg(pft,iseg)
    ! Deepest Newton step below psi_in: k there is below the path's kcrit
    ! on this segment's curve, a fraction kcrit / kmx of the segment's own
    ! kmax (newton_psi_drop; 5 |bs| unless the curve is wide). For the
    ! shared root and stem the path is the whole plant's, with its kcrit in
    ! proportion to its kmax.
    kcr_drop = kcrit(l)
    IF ( l_shared_stem_path .AND. ( iseg <= 2 .OR. l_shared_whole_path ) )     &
      kcr_drop = kcr_drop * shared_k_fac(l)
    drop_nr = newton_psi_drop(bs, cs,                                          &
                              kcr_drop / MAX(kmx, TINY(1.0_real_jlslsm)),      &
                              .FALSE.)

    ! Cap on the stem and leaf segments (not the root) under xylem
    ! impairment: their conductance cannot exceed fr * kmx. psi_cap is
    ! where the intact curve falls to the cap,
    !
    !   kmx exp(-(psi_cap / b)^c) = fr kmx,   psi_cap = b (-ln fr)^(1/c)
    !
    ! so above psi_cap (wetter) the cap binds and below it the intact curve
    ! does. a_cap is the intact supply integral at psi_cap, used to join the
    ! two parts of the capped integral (seg_a).
    ! With l_som_seg_whole the root is capped too (whole-plant damage).
    fr = 1.0
    IF ( l_impair .AND. ( iseg >= 2 .OR. l_som_seg_whole ) )                   &
      fr = MIN(MAX(ximpair_seg_frac(pft, l, iseg, kcap_frac(l)), k_floor), 1.0)
    l_cap = fr < 1.0 - 1.0e-6
    IF ( l_cap ) THEN
      psi_cap = bs * (-LOG(fr))**(1.0 / cs)
      a_cap = kmx * (-bs/cs) * incomplete_gamma_1(1.0/cs, (psi_cap/bs)**cs)
    ELSE
      psi_cap = 0.0
      a_cap = 0.0
    END IF

    ! Conductance and supply integral at the root-side end of the segment.
    k_in(:) = seg_k(psi_in)
    IF ( iseg == 1 .AND. .NOT. ( l_som_rhizo_series .AND.                      &
                                 ksr_path(l) >= 0.0 ) ) THEN
      ! The root water potential is psi_root_zone for every sample: one
      ! evaluation.
      ! (With l_som_rhizo_series it falls with E, so not here.)
      a_in(:) = seg_a(psi_in(1))
    ELSE
      a_in(:) = seg_a(psi_in)
    END IF

    ! First guess: the root-side conductance over the whole drop (Ohm's
    ! law, psi_out = psi_in - E / k_in), no more than drop_nr below psi_in.
    psi_out(:) = psi_in(:) - e_seg(:) / MAX(k_in(:), TINY(1.0_real_jlslsm))
    psi_out(:) = MAX(psi_out(:), psi_in(:) - drop_nr)
    k_cur(:)   = k_in(:)

    ! Newton-Raphson for the leaf-side water potential: the segment carries
    !
    !   E(psi_out) = seg_a(psi_out) - a_in
    !
    ! and dE/dpsi_out = -k(psi_out), so each step moves psi_out by the
    ! missing flow divided by the conductance there. Each sample stops
    ! when its own k has settled (changes by less than 0.1 kcrit), or
    ! after max_nr_iter steps.
    l_nr_active(:) = .TRUE.
    DO it = 1, max_nr_iter
      k_prev(:) = k_cur(:)
      k_cur(:)  = seg_k(psi_out)
      e_cur(:)  = seg_a(psi_out) - a_in(:)
      WHERE (l_nr_active(:))                                                   &
        psi_out(:) = psi_out(:) - (e_seg(:) - e_cur(:))                        &
                                  / MAX(k_cur(:), TINY(1.0_real_jlslsm))
      ! Same runaway guard as the single segment (infeasible demand).
      psi_out(:) = MAX(psi_out(:), psi_in(:) - drop_nr)
      ! A flow back from the store (E_r < 0) larger than the segment can
      ! carry up to psi = 0: the leaf-side end at 0.
      IF ( l_store_q_path ) psi_out(:) = MIN(psi_out(:), 0.0_real_jlslsm)
      WHERE (ABS(k_cur(:) - k_prev(:)) < 0.1 * kcr) l_nr_active(:) = .FALSE.
      IF ( .NOT. ANY(l_nr_active(:)) ) EXIT
    END DO

    ! Optionally make psi_out exact by bisection (see seg_bisect).
    IF ( l_exact_on ) THEN
      DO i = 1, n_e_leaf
        CALL seg_bisect( i )
      END DO
    END IF

    ! Segment conductances are floored at a tiny fraction of the segment
    ! kmax (not TINY): a collapsed segment has k -> 0 by underflow, and with
    ! a TINY floor the recursion overflows to -Inf and then 0 * -Inf = NaN
    ! in the next segment, a huge whole-plant k that would let an infeasible
    ! state (psi_leaf ~ -70 MPa) pass the k > kcrit mask. With this floor
    ! dpsi_de stays finite and k correctly ~ 0.
    !
    ! Carry dpsi/dE down to the leaf-side end of this segment (header):
    !
    !   dpsi_out/dE = (k(psi_in) dpsi_in/dE - 1) / k(psi_out)
    !
    k_cur(:)   = MAX(seg_k(psi_out), k_floor * kmx)
    dpsi_de(:) = ( MAX(k_in(:), k_floor * kmx) * dpsi_de(:) - 1.0 ) / k_cur(:)

    ! The same for the intact segments at the same water potentials (the
    ! conductance the hydraulic cost uses under impairment).
    IF ( PRESENT(leaf_k_intact) ) THEN
      dpsi_de_int(:) = ( MAX(kmx * EXP( -(MIN(psi_in(:), 0.0_real_jlslsm)      &
                                     / bs)**cs ), k_floor * kmx)               &
                         * dpsi_de_int(:) - 1.0 )                              &
                       / MAX(kmx * EXP( -(MIN(psi_out(:), 0.0_real_jlslsm)     &
                                        / bs)**cs ), k_floor * kmx)
    END IF

    ! Hand on to the next segment; keep the stem segment's leaf-side psi.
    psi_in(:)  = psi_out(:)
    IF ( iseg == 2 .AND. PRESENT(stem_psi) ) stem_psi(:,j) = psi_out(:)
    ! Keep the stem state of a class solve (l_rs_cache).
    IF ( iseg == 2 .AND. l_one .AND. .NOT. l_shared_whole_path ) THEN
      l_rs_cache  = .TRUE.
      rs_l        = l
      rs_e        = e_leaf(1,j)
      rs_kmax     = kmax(l)
      rs_kcrit    = kcrit(l)
      rs_psi0     = root_zone_psi(l)
      rs_kcap     = kc_key
      rs_intact   = PRESENT(leaf_k_intact)
      rs_exact    = l_exact_on
      rs_psi      = psi_out(1)
      rs_dpsi     = dpsi_de(1)
      rs_dpsi_int = dpsi_de_int(1)
      rs_closed   = l_closed(1)
    END IF
  END DO

  ! The leaf-side end of the leaf segment is the leaf water potential, and
  ! the whole-plant conductance is k = -dE/dpsi_leaf = -1 / (dpsi_leaf/dE).
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

  ! Samples the path cannot carry (no soil conductance, or a flow above the
  ! most a segment can carry with l_exact): infeasible, k = 0.
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

  ! Segment conductance (capped where l_cap). psi is held at <= 0 here and in
  ! seg_a: a root-zone psi a hair above 0 (rounding) would give a negative
  ! base to a non-integer power (NaN) and close a wet point. ELEMENTAL (one
  ! psi at a time), so the calls make no array temporaries.
  ELEMENTAL FUNCTION seg_k( psi ) RESULT( k )
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  REAL(KIND=real_jlslsm) :: k
  k = kmx * EXP( -(MIN(psi, 0.0_real_jlslsm)/bs)**cs )
  IF ( l_cap ) k = MIN(k, fr * kmx)
  END FUNCTION seg_k

  ! l_exact: the leaf-side psi of sample i by bisection, seg_a(psi_out) =
  ! a_in + E (seg_a rises as psi falls), between the root-side psi and the psi
  ! where (psi/b)**c = 12 (f = 6e-6: the curve has carried all it can, and
  ! the series in incomplete_gamma is still converged there; deeper, its
  ! 40 terms are not enough and seg_a is no longer monotone). A flow above
  ! that is infeasible.
  SUBROUTINE seg_bisect( i )
  INTEGER, INTENT(IN) :: i
  INTEGER :: ib
  REAL(KIND=real_jlslsm) :: lo, hi, mid, a_t
  a_t = a_in(i) + e_seg(i)
  hi  = psi_in(i)
  lo  = -ABS(bs) * 12.0**(1.0 / cs)
  IF ( l_closed(i) ) RETURN
  IF ( lo >= hi ) THEN
    l_closed(i) = .TRUE.
    psi_out(i)  = lo
    RETURN
  END IF
  IF ( seg_a(lo) < a_t ) THEN
    l_closed(i) = .TRUE.
    psi_out(i)  = lo
    RETURN
  END IF
  DO ib = 1, 60
    mid = 0.5 * (lo + hi)
    IF ( seg_a(mid) < a_t ) THEN
      hi = mid
    ELSE
      lo = mid
    END IF
  END DO
  psi_out(i) = 0.5 * (lo + hi)
  END SUBROUTINE seg_bisect

  ! Integral of the segment conductance from psi to 0 (so that the flow
  ! from psi_in to psi_out is seg_a(psi_out) - seg_a(psi_in)).
  ! ELEMENTAL, as seg_k.
  ELEMENTAL FUNCTION seg_a( psi ) RESULT( a )
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  REAL(KIND=real_jlslsm) :: a
  ! Note: the negative sign is present because bs is negative.
  a = kmx * (-bs/cs) * incomplete_gamma_1(1.0/cs,                              &
                                          (MIN(psi, 0.0_real_jlslsm)/bs)**cs)
  IF ( l_cap ) THEN
    IF ( psi >= psi_cap ) THEN
      a = fr * kmx * (-psi)
    ELSE
      a = fr * kmx * (-psi_cap) + a - a_cap
    END IF
  END IF
  END FUNCTION seg_a

END SUBROUTINE psi_leaf_segments_jls

! *****************************************************************************
! As psi_leaf_segments_jls, with each segment solved from its supply-function
! table instead of Newton-Raphson (psi_solver_lut): segment s carries
!
!   E = kmax_s ( S_s(psi_out) - S_s(psi_in) ),   kmax_s = kmax seg_kfac,
!
! so psi_out = S_s^-1( S_s(psi_in) + E / kmax_s ), downstream from
! psi_root_zone. The whole-plant conductance -dE/dpsi_leaf comes from the
! same recursion,
!
!   dpsi_out/dE = (k_s(psi_in) dpsi_in/dE - 1) / k_s(psi_out),
!
! with the same floor on segment conductances. A demand beyond what a
! segment can supply puts its leaf-side end at the bottom of its table, where k
! is ~0, so the state is rejected as infeasible (k <= kcrit).
!
! Xylem impairment (memory model with segments; kcap_frac present): as
! psi_leaf_segments_jls, the stem and leaf segments are capped at a fraction
! fr of their kmax_s (ximpair_seg_frac), the root is intact. The capped
! integral is
!
!   A_s(psi) = fr kmax_s (-psi)                                psi >= psi_cap
!            = fr kmax_s (-psi_cap) + kmax_s (S_s(psi) - S_s(psi_cap))  below,
!
! psi_cap where the intact curve meets the cap, f_s(psi_cap) = fr. So for an
! root-side psi above psi_cap the flow down to psi_cap, e_cap = fr kmax_s
! (psi_in - psi_cap), is carried at the constant conductance fr kmax_s: an E <=
! e_cap gives psi_out = psi_in - E / (fr kmax_s) exactly, and a larger E
! continues on the intact curve from psi_cap, S_s(psi_out) = S_s(psi_cap) + (E -
! e_cap) / kmax_s, one inversion of the unchanged table. A root-side psi below
! psi_cap is on the intact curve throughout. The caps are per-point numbers, so
! the tables need no rebuilding. leaf_k uses the capped segment conductances in
! the recursion; leaf_k_intact (if present) the intact ones at the same water
! potentials (the hydraulic cost), and stem_psi the stem's leaf-side psi.
!
! The plant capacitance store (l_store_q_path) and the shared root and stem
! (l_shared_stem_path) as psi_leaf_segments_jls: the soil link, the root
! and the stem carry E_r = E - q_f + E_o, the leaf E (plus E_o on the
! whole-plant path). E_r < 0 (the store giving more than E) is a flow back
! towards the soil, found upwards in the table (seg_lut_psi_back) and capped
! at psi = 0, as the Newton solver's cap.
!
! MDK, 31/10/26
! *****************************************************************************
SUBROUTINE psi_leaf_segments_lut_jls( pft, n_e_leaf, land_pts, open_pnts,      &
                                      veg_index, open_index, e_leaf,           &
                                      root_zone_psi, kmax, leaf_psi, leaf_k,   &
                                      kcap_frac, leaf_k_intact, stem_psi )

USE pftparm, ONLY: conductance_b_seg, conductance_c_seg
USE jules_vegetation_mod, ONLY: l_som_rhizo_series
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! As psi_leaf_segments_jls.
  e_leaf(n_e_leaf, open_pnts), root_zone_psi(land_pts), kmax(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Leaf water potential (Pa) and whole-plant -dE/dpsi_leaf.
  leaf_psi(n_e_leaf, open_pnts), leaf_k(n_e_leaf, open_pnts)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  ! Impaired / intact maximum conductance of the stem (and leaf) segment,
  ! as psi_leaf_segments_jls (kmax is then the intact value).
  kcap_frac(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  ! Whole-plant conductance of the intact segments at leaf_psi.
  leaf_k_intact(n_e_leaf, open_pnts),                                          &
  ! Water potential at the leaf-side end of the stem segment (Pa).
  stem_psi(n_e_leaf, open_pnts)

INTEGER, PARAMETER :: n_seg = 3

! Floor on a segment's conductance, as a fraction of its kmax (as
! psi_leaf_segments_jls).
REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12

INTEGER :: i, j, l, iseg
LOGICAL :: l_cap, l_impair
REAL(KIND=real_jlslsm) :: kmx, bs, cs, psi_in, psi_out, k_in, k_out, dpsi_de,  &
                          k_s, fr, psi_cap, e_cap, e_pos, dpsi_de_int
      ! fr, psi_cap: the segment's cap fraction and where the intact curve
      ! meets it; e_cap the flow the capped part above psi_cap carries;
      ! e_pos the flow (>= 0); dpsi_de_int d(psi)/dE on the intact segments.
! The flow of the soil link, root and stem (E_r) and of the current
! segment; with a flow back, the flow per unit kmax from psi_in up to psi_cap
! (< 0) and the flow per unit kmax.
REAL(KIND=real_jlslsm) :: e_r, e_s, x_cap, x_s

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_SEGMENTS_LUT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL build_supply_lut_seg(pft)
l_impair = PRESENT(kcap_frac)

DO j = 1, open_pnts
  l = veg_index(open_index(j))
  DO i = 1, n_e_leaf
    psi_in  = root_zone_psi(l)
    dpsi_de = 0.0
    dpsi_de_int = 0.0
    ! The flow up to the stem node, as psi_leaf_segments_jls: less the
    ! store's flow held fixed (not floored at 0), plus the other class's
    ! flow on the shared root and stem.
    e_r = e_leaf(i,j)
    IF ( l_store_q_path ) e_r = e_r - store_q_path(l)
    IF ( l_shared_stem_path ) e_r = e_r + shared_e_other(l)
    ! l_som_rhizo_series: soil link ahead of the root segment (as
    ! psi_leaf_segments_jls).
    IF ( l_som_rhizo_series .AND. ksr_path(l) >= 0.0 ) THEN
      k_s = ksr_path(l)
      IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
        psi_in  = MAX(root_zone_psi(l) - e_r / k_s, som_psi_in_min)
        dpsi_de = -1.0 / k_s
        dpsi_de_int = -1.0 / k_s
      ELSE IF ( e_r > 0.0 ) THEN
        leaf_psi(i,j) = som_psi_in_min
        leaf_k(i,j)   = 0.0
        IF ( PRESENT(leaf_k_intact) ) leaf_k_intact(i,j) = 0.0
        IF ( PRESENT(stem_psi) ) stem_psi(i,j) = som_psi_in_min
        CYCLE
      END IF
    END IF

    ! Root, stem then leaf in series. Each segment carries the same E, so its
    ! leaf-side psi comes from its own table with E divided by its own kmax, and
    ! becomes the next segment's root-side psi. dpsi_de carries dpsi/dE down the
    ! chain: from E = integral(k), k_out dpsi_out = k_in dpsi_in - dE.
    DO iseg = 1, n_seg
      kmx = kmax(l) * seg_kmax_fac(pft,l,iseg)
      bs  = conductance_b_seg(pft,iseg)
      cs  = conductance_c_seg(pft,iseg)

      ! The leaf carries the leaf's E (and the other class's on the
      ! whole-plant path); the soil link, root and stem E_r.
      e_s = e_r
      IF ( iseg == n_seg ) e_s = e_leaf(i,j)
      IF ( iseg == n_seg .AND. l_shared_stem_path .AND. l_shared_whole_path )  &
        e_s = e_s + shared_e_other(l)

      ! Cap on the stem and leaf segments (not the root; all three with
      ! l_som_seg_whole), as psi_leaf_segments_jls.
      l_cap = .FALSE.
      fr = 1.0
      IF ( l_impair .AND. ( iseg >= 2 .OR. l_som_seg_whole ) ) THEN
        fr = MIN(MAX(ximpair_seg_frac(pft, l, iseg, kcap_frac(l)), k_floor),   &
                 1.0)
        l_cap = fr < 1.0 - 1.0e-6
      END IF

      IF ( l_cap .AND. l_store_q_path .AND. e_s < 0.0 ) THEN
        ! A flow back towards the soil on a capped segment: up the intact
        ! curve to psi_cap, then up the straight capped part, to at most
        ! psi = 0.
        psi_cap = bs * (-LOG(fr))**(1.0 / cs)
        x_s = e_s / MAX(kmx, TINY(1.0_real_jlslsm))
        IF ( psi_in >= psi_cap ) THEN
          psi_out = psi_in - e_s / MAX(fr * kmx, TINY(1.0_real_jlslsm))
        ELSE
          x_cap = seg_lut_s(pft, iseg, psi_cap) - seg_lut_s(pft, iseg, psi_in)
          IF ( x_s >= x_cap ) THEN
            psi_out = seg_lut_psi_back(pft, iseg, psi_in, x_s)
          ELSE
            psi_out = psi_cap - (x_s - x_cap) / fr
          END IF
        END IF
        psi_out = MIN(psi_out, 0.0_real_jlslsm)
        k_in  = MAX(MIN(kmx * EXP( -(MIN(psi_in, 0.0_real_jlslsm)/bs)**cs ),   &
                        fr * kmx), k_floor * kmx)
        k_out = MAX(MIN(kmx * EXP( -(MIN(psi_out, 0.0_real_jlslsm)/bs)**cs ),  &
                        fr * kmx), k_floor * kmx)
      ELSE IF ( l_store_q_path .AND. e_s < 0.0 ) THEN
        ! A flow back on an intact segment (as above).
        psi_out = seg_lut_psi_back(pft, iseg, psi_in,                          &
                                   e_s / MAX(kmx, TINY(1.0_real_jlslsm)))
        k_in  = MAX(kmx * EXP( -(MIN(psi_in, 0.0_real_jlslsm)/bs)**cs ),       &
                    k_floor * kmx)
        k_out = MAX(kmx * EXP( -(MIN(psi_out, 0.0_real_jlslsm)/bs)**cs ),      &
                    k_floor * kmx)
      ELSE IF ( l_cap ) THEN
        ! Exact inversion of the capped integral (see above): the straight
        ! capped part down to psi_cap, then the intact table from psi_cap.
        psi_cap = bs * (-LOG(fr))**(1.0 / cs)
        e_pos   = MAX(e_s, 0.0)
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
                              e_s / MAX(kmx, TINY(1.0_real_jlslsm)))
        k_in  = MAX(kmx * EXP( -(ABS(psi_in /bs))**cs ), k_floor * kmx)
        k_out = MAX(kmx * EXP( -(ABS(psi_out/bs))**cs ), k_floor * kmx)
      END IF

      dpsi_de = ( k_in * dpsi_de - 1.0 ) / k_out

      ! The intact chain at the same water potentials (hydraulic cost).
      IF ( PRESENT(leaf_k_intact) )                                            &
        dpsi_de_int = ( MAX(kmx * EXP( -(MIN(psi_in, 0.0_real_jlslsm)          &
                                         / bs)**cs ), k_floor * kmx)           &
                        * dpsi_de_int - 1.0 )                                  &
                      / MAX(kmx * EXP( -(MIN(psi_out, 0.0_real_jlslsm)         &
                                       / bs)**cs ), k_floor * kmx)

      psi_in  = psi_out
      IF ( iseg == 2 .AND. PRESENT(stem_psi) ) stem_psi(i,j) = psi_out
    END DO

    ! The leaf segment's leaf-side psi is the leaf psi. The whole-plant
    ! conductance is -dE/dpsi_leaf = -1/dpsi_de; it is zero if the chain gave no
    ! usable slope (e.g. a segment at its conductance floor).
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_segments_lut_jls

! *****************************************************************************
! Critical conductance of a leaf path (mol m-2 s-1 Pa-1, the units of
! kmax_ref): the PFT's critical fraction kcrit/kmax_pft of the path's
! intact maximum conductance kmax_ref. One rule for every canopy (big leaf,
! two-leaf, multilayer and the cuticular floor), which wrote it as
! kcrit * LAI or kmax_ref * (kcrit / kmax_pft): equal, but not to the last
! bit.
!
! MDK, 23/10/26
! *****************************************************************************
ELEMENTAL FUNCTION kcrit_of_path( ft, kmax_ref ) RESULT( kc )

USE pftparm, ONLY: kcrit, kmax_pft

IMPLICIT NONE

! Plant functional type.
INTEGER, INTENT(IN) :: ft
! Intact maximum conductance of the path.
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax_ref
! Critical conductance of the path.
REAL(KIND=real_jlslsm) :: kc

kc = kmax_ref * (kcrit(ft) / kmax_pft(ft))

END FUNCTION kcrit_of_path

END MODULE plant_hydraulic_path_mod
