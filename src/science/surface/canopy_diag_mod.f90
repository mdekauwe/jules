! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE canopy_diag_mod

!-----------------------------------------------------------------------------
! Diagnostics of the canopy that no prognostic or output field holds:
!
! - gc_stomatal_pft: the stomatal conductance the optimiser chose (sum of
!   the leaf classes), before the transpiration-equivalent gc
!   (l_leaf_temp_gc_eq) replaces gc, and before the cuticular leak or the
!   store change it (m s-1, ground).
! - gsoil_pft: the soil evaporation conductance beneath the PFT's canopy as
!   added to the tile's surface conductance (Or or trunk; m s-1).
! - t_can_air_pft, q_can_air_pft: the canopy air temperature (K) and
!   specific humidity (kg kg-1) of the leaf energy balance (l_leaf_temp).
! - h_leaf_pft, le_leaf_pft: the leaves' sensible and latent heat in that
!   balance (W m-2, ground).
!
! Allocated on first use; 0 where not set. Written by sf_stom (and physiol
! for gsoil_pft), read by extract_var.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: gc_stomatal_pft, gsoil_pft, t_can_air_pft, q_can_air_pft,           &
          h_leaf_pft, le_leaf_pft, h_leaf_pts, le_leaf_pts, canopy_diag_alloc

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  gc_stomatal_pft(:,:), gsoil_pft(:,:), t_can_air_pft(:,:),                   &
  q_can_air_pft(:,:), h_leaf_pft(:,:), le_leaf_pft(:,:)
    ! (land_pts, npft): see the header.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: h_leaf_pts(:), le_leaf_pts(:)
    ! (land_pts): the leaf balance's H and LE of the current PFT, passed
    ! from leaf_temp_update_n to the caller (W m-2, ground).

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CANOPY_DIAG_MOD'

CONTAINS

SUBROUTINE canopy_diag_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CANOPY_DIAG_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ALLOCATED(gc_stomatal_pft) ) THEN
  ALLOCATE( gc_stomatal_pft(n_land_pts,npft), gsoil_pft(n_land_pts,npft),     &
            t_can_air_pft(n_land_pts,npft), q_can_air_pft(n_land_pts,npft),   &
            h_leaf_pft(n_land_pts,npft), le_leaf_pft(n_land_pts,npft),        &
            h_leaf_pts(n_land_pts), le_leaf_pts(n_land_pts) )
  gc_stomatal_pft(:,:) = 0.0
  gsoil_pft(:,:)       = 0.0
  t_can_air_pft(:,:)   = 0.0
  q_can_air_pft(:,:)   = 0.0
  h_leaf_pft(:,:)      = 0.0
  le_leaf_pft(:,:)     = 0.0
  h_leaf_pts(:)        = 0.0
  le_leaf_pts(:)       = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE canopy_diag_alloc

END MODULE canopy_diag_mod
