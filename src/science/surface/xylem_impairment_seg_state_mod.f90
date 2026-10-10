! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_seg_state_mod

! Per-segment memory (l_ximpair_seg_memory) for the memory impairment model
! with root / stem / leaf segments. The stem segment's cap is the impaired
! kmax state (k_max / kmax_impaired, leaf basis), as for the single cap;
! the leaf segment's cap is held here as a fraction of the intact segment
! maximum. Kept apart from xylem_impairment_memory_mod so the segment
! solvers can read it without a circular module dependency.

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: ximpair_kcap_leaf, ximpair_seg_state_alloc, ximpair_seg_frac

! Leaf-segment cap, k_cap / kmax of the leaf segment (1 = intact), of each
! PFT at each land point. Not allocated until the first update (none during
! spin-up) or a dump read, and read as intact until then. Written to and
! read from dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: ximpair_kcap_leaf(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the leaf-segment cap (n_land_pts, npft), intact, if not done.
! ---------------------------------------------------------------------
SUBROUTINE ximpair_seg_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(ximpair_kcap_leaf)) THEN
  ALLOCATE(ximpair_kcap_leaf(n_land_pts, npft))
  ximpair_kcap_leaf(:,:) = 1.0
END IF

END SUBROUTINE ximpair_seg_state_alloc

! ---------------------------------------------------------------------
! Cap fraction of segment iseg (1 root, 2 stem, 3 leaf) at land point l:
! the root is never capped; the stem takes kcap_frac (the impaired kmax
! state over the intact kmax); the leaf takes kcap_frac too, or its own
! cap with l_ximpair_seg_memory.
! ---------------------------------------------------------------------
FUNCTION ximpair_seg_frac( pft, l, iseg, kcap_frac ) RESULT( fr )

USE jules_vegetation_mod, ONLY: l_ximpair_seg_memory

INTEGER, INTENT(IN) :: pft, l, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: kcap_frac
REAL(KIND=real_jlslsm) :: fr

SELECT CASE (iseg)
CASE (1)
  fr = 1.0
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
