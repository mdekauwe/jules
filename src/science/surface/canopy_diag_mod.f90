! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE canopy_diag_mod

!-----------------------------------------------------------------------------
! Canopy diagnostics that no other model or output variable holds. Values
! marked "ground" are per square metre of ground; "leaf" per square metre
! of leaf.
!
! - gc_stomatal_pft: the canopy stomatal conductance the optimisation chose
!   (summed over sunlit and shaded leaves or canopy layers), before it is
!   replaced by the conductance that gives the same transpiration at the
!   leaf temperature (l_leaf_temp_gc_eq), and before the cuticular leak or
!   the stem store change it (m s-1, ground).
! - gsoil_pft: the soil evaporation conductance under the PFT's canopy, as
!   added to the tile's surface conductance (m s-1), from the Or et al.
!   scheme (l_soil_evap_or) or the standard JULES one.
! - t_can_air_pft, q_can_air_pft: the canopy air temperature (K) and
!   specific humidity (kg kg-1) from the leaf energy balance (l_leaf_temp).
! - h_leaf_pft, le_leaf_pft: the sensible and latent heat from the leaves
!   in that energy balance (W m-2, ground).
! - gs_sun_pft, gs_shd_pft: the stomatal conductance the optimisation chose
!   for the sunlit and the shaded leaves, per unit area of those leaves
!   (m s-1), and lai_sun_pft, the sunlit leaf area (m2 m-2). Two-leaf
!   canopy (can_rad_mod = 7) only. These compare with leaf gas exchange
!   measurements.
! - vcmax_pft, jmax_pft: Vcmax and Jmax of a leaf at the top of the canopy
!   at the canopy temperature tstar (mol m-2 s-1, leaf), with any thermal
!   acclimation and the soil-water factor (l_som_vcmax_psi). With
!   Johnson-Berry photosynthesis, jmax_pft holds Vqmax. They are at tstar,
!   not at the sunlit or shaded leaf temperatures (l_leaf_temp).
! - anetc_pft: net canopy photosynthesis (mol CO2 m-2 s-1, ground).
! - anetl_top_pft, anetl_bot_pft: net leaf photosynthesis of the top and
!   the bottom canopy layer (mol CO2 m-2 s-1, leaf). Layered canopies
!   (can_rad_mod = 4, 5, 6) only.
!   The five above, from vcmax_pft on, are worked out only when one of
!   them is output (l_photo_diag).
! - psi_leaf_sun_pft, psi_leaf_shd_pft: the water potential of the sunlit
!   and the shaded leaves at the conductance the optimisation chose (Pa),
!   before any cuticular leak or stem store balance. Two-leaf canopy
!   (can_rad_mod = 7) only.
!
! The arrays are allocated on first use and are 0 where not set. sf_stom
! writes them (physiol writes gsoil_pft) and extract_var reads them for
! output.
!
! MDK, 30/10/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: gc_stomatal_pft, gsoil_pft, t_can_air_pft, q_can_air_pft,           &
          h_leaf_pft, le_leaf_pft, h_leaf_pts, le_leaf_pts, canopy_diag_alloc, &
          gs_sun_pft, gs_shd_pft, lai_sun_pft, vcmax_pft, jmax_pft,           &
          anetc_pft, anetl_top_pft, anetl_bot_pft, psi_leaf_sun_pft,          &
          psi_leaf_shd_pft

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  ! (land_pts, npft): see the header.
  gc_stomatal_pft(:,:), gsoil_pft(:,:), t_can_air_pft(:,:),                   &
  q_can_air_pft(:,:), h_leaf_pft(:,:), le_leaf_pft(:,:),                     &
  gs_sun_pft(:,:), gs_shd_pft(:,:), lai_sun_pft(:,:), vcmax_pft(:,:),       &
  jmax_pft(:,:), anetc_pft(:,:), anetl_top_pft(:,:), anetl_bot_pft(:,:),     &
  psi_leaf_sun_pft(:,:), psi_leaf_shd_pft(:,:)
! (land_pts): the leaf balance's H and LE of the current PFT, passed
! from leaf_temp_update_n to the caller (W m-2, ground).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: h_leaf_pts(:), le_leaf_pts(:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CANOPY_DIAG_MOD'

CONTAINS

! Allocate the canopy diagnostics, filled with 0, on first use.
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
            h_leaf_pts(n_land_pts), le_leaf_pts(n_land_pts),                   &
            gs_sun_pft(n_land_pts,npft), gs_shd_pft(n_land_pts,npft),         &
            lai_sun_pft(n_land_pts,npft), vcmax_pft(n_land_pts,npft),         &
            jmax_pft(n_land_pts,npft), anetc_pft(n_land_pts,npft),            &
            anetl_top_pft(n_land_pts,npft), anetl_bot_pft(n_land_pts,npft),   &
            psi_leaf_sun_pft(n_land_pts,npft),                                 &
            psi_leaf_shd_pft(n_land_pts,npft) )
  gc_stomatal_pft(:,:) = 0.0
  gsoil_pft(:,:)       = 0.0
  t_can_air_pft(:,:)   = 0.0
  q_can_air_pft(:,:)   = 0.0
  h_leaf_pft(:,:)      = 0.0
  le_leaf_pft(:,:)     = 0.0
  h_leaf_pts(:)        = 0.0
  le_leaf_pts(:)       = 0.0
  gs_sun_pft(:,:)      = 0.0
  gs_shd_pft(:,:)      = 0.0
  lai_sun_pft(:,:)     = 0.0
  vcmax_pft(:,:)       = 0.0
  jmax_pft(:,:)        = 0.0
  anetc_pft(:,:)       = 0.0
  anetl_top_pft(:,:)   = 0.0
  anetl_bot_pft(:,:)   = 0.0
  psi_leaf_sun_pft(:,:) = 0.0
  psi_leaf_shd_pft(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE canopy_diag_alloc

END MODULE canopy_diag_mod
