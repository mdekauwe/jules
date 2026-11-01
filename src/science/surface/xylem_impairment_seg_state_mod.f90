! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_seg_state_mod

! Xylem impairment, memory model: the damage state of each part of the
! plant when the plant is split into root, stem and leaf segments
! (l_ximpair_seg_memory). Each segment can lose conductance on its own.
!
! The stem's damage is the plant's damaged kmax (kmax_impaired, leaf area
! basis), the same state the memory model uses without segments. The
! leaf's damage is held here as a cap on its conductance, a fraction of
! the undamaged leaf segment's maximum (1 = undamaged).
!
! This module is separate from xylem_impairment_memory_mod so the segment
! solvers can read it without the two modules using each other.

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: ximpair_kcap_leaf, ximpair_psi_stem, ximpair_seg_state_alloc,        &
          ximpair_seg_frac

! Leaf-segment cap, k_cap / kmax of the leaf segment (1 = intact), of each
! PFT at each land point. Not allocated until the first update (none during
! spin-up) or a dump read, and read as intact until then. Written to and
! read from dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: ximpair_kcap_leaf(:,:)

! Stem water potential (Pa; leaf-side end of the stem segment) at the last
! update, for the PLC_pft diagnostic. Not dumped; 0 until the first update.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: ximpair_psi_stem(:,:)

CHARACTER(LEN=*), PARAMETER :: ModuleName='XYLEM_IMPAIRMENT_SEG_STATE_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Xylem impairment, memory model: create the leaf segment's cap (set to
! 1, undamaged) and the stored stem water potential (set to 0) for every
! land point and PFT, if they do not exist yet.
! ---------------------------------------------------------------------
SUBROUTINE ximpair_seg_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XIMPAIR_SEG_STATE_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(ximpair_kcap_leaf)) THEN
  ALLOCATE(ximpair_kcap_leaf(n_land_pts, npft))
  ximpair_kcap_leaf(:,:) = 1.0
END IF

IF (.NOT. ALLOCATED(ximpair_psi_stem)) THEN
  ALLOCATE(ximpair_psi_stem(n_land_pts, npft))
  ximpair_psi_stem(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ximpair_seg_state_alloc

! ---------------------------------------------------------------------
! Return how much of its undamaged maximum conductance segment iseg
! (1 root, 2 stem, 3 leaf) can still carry at land point l after drought
! damage (0-1, memory model of xylem impairment), which the segment
! solvers use to cap each segment.
!
! kcap_frac is the plant's damaged kmax over its undamaged kmax.
!
!   - The root is undamaged (1), except with l_som_seg_whole (segments on
!     the PFT's whole-plant curve, used with plant capacitance), where it
!     takes kcap_frac.
!   - The stem takes kcap_frac.
!   - The leaf takes kcap_frac, or its own cap with l_ximpair_seg_memory
!     (1 until the cap exists).
!
! MDK, 26/10/26
! ---------------------------------------------------------------------
FUNCTION ximpair_seg_frac( pft, l, iseg, kcap_frac ) RESULT( fr )

USE jules_vegetation_mod, ONLY: l_ximpair_seg_memory, l_som_seg_whole

INTEGER, INTENT(IN) :: pft, l, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: kcap_frac
REAL(KIND=real_jlslsm) :: fr

SELECT CASE (iseg)
CASE (1)
  fr = 1.0
  IF ( l_som_seg_whole ) fr = kcap_frac
CASE (2)
  fr = kcap_frac
CASE DEFAULT
  fr = kcap_frac
  IF ( l_ximpair_seg_memory ) THEN
    fr = 1.0
    IF ( ALLOCATED(ximpair_kcap_leaf) ) fr = ximpair_kcap_leaf(l,pft)
  END IF
END SELECT

END FUNCTION ximpair_seg_frac

END MODULE xylem_impairment_seg_state_mod
