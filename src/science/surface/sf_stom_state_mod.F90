! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE
MODULE sf_stom_state_mod
! Values that sf_stom keeps from one timestep to the next. They live here,
! outside sf_stom, so they can be written to the dump and a run restarted
! from a dump gives the same answers as one that ran on.
!
!   - f_vcmax_psi_state: the factor that lowers Vcmax as the soil dries
!     (l_som_vcmax_psi), as it was after its lag in recovering when the soil
!     wets up again (0 to 1; dump name fvcmax_psi_state).
!   - cpl_state: with coupled transpiration (l_som_coupled_e) and the
!     two-leaf canopy, the last stomatal and leaf boundary layer
!     conductances of the sunlit and shaded leaves (m/s), the first guess
!     on the first pass of the next timestep (dump names cpl_gl_sun,
!     cpl_gl_shd, cpl_gbw_sun, cpl_gbw_shd). Its fifth slot is the
!     aerodynamic resistance 1 / (ch * |V|) of the last surface flux step
!     (s/m), which without leaf temperatures is the stomatal
!     optimisation's resistance to the first atmospheric level (dump name
!     cpl_ra_ca).
!
! The module also holds what sf_stom hands to the surface flux step within
! a timestep (e_flux_tile, gc_cap_tile) and physiol's root extraction
! weights (wt_root_tile). These are set each timestep and are not in the
! dump.
!
! f_vcmax_psi_state and cpl_state are negative until first set, meaning
! "not yet set". A dump without them gives the same.
!
! MDK, 26/11/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: f_vcmax_psi_state, cpl_state, sf_stom_state_alloc, n_cpl_state,      &
          shared_e_shd_prev, shared_amax_shd_prev, store_e_prev,           &
          e_flux_tile, gc_cap_tile, wt_root_tile, sf_stom_flux_alloc

! gl_sun, gl_shd, gbw_sun, gbw_shd, the flux step's resistance
INTEGER, PARAMETER :: n_cpl_state = 5

! (land_pts, npft)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: f_vcmax_psi_state(:,:)
! (land_pts, npft, n_cpl_state)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cpl_state(:,:,:)

! (land_pts, npft) The shared root and stem (l_som_shared_stem),
! two-leaf: the shaded leaves' transpiration of the last timestep
! (mol H2O m-2 s-1), the flow the sunlit leaves' search puts on the
! shared root and stem in a timestep's first pass. Allocated and set to
! 0 by sf_stom or the dump reading; in the dump, so a restart in daylight
! starts from the same first guess as a continuous run.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shared_e_shd_prev(:,:)
! (land_pts, npft) The same for the shaded leaves' A at their hydraulic
! edge (mol CO2 m-2 s-1), for the whole-plant gain (l_som_shared_gain).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shared_amax_shd_prev(:,:)
! (land_pts, npft) With the capacitance store coupled to the two-leaf
! search: the last timestep's stomatal transpiration of the sunlit and
! shaded leaves together (mol H2O m-2 s-1, >= 0; 0 where the canopy was
! shut). Where it is 0, the first pass does not use the last step's store
! refill. Allocated and set to 0 by sf_stom or the dump reading; in the
! dump with the store's flow (store_q).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_e_prev(:,:)
! (land_pts, npft) The transpiration of this timestep as sf_stom leaves it
! (kg m-2 s-1; negative: none), which the surface flux step delivers.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: e_flux_tile(:,:)
! (land_pts, npft) The most the flux step may set the canopy stomatal
! conductance to in delivering it (m s-1).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: gc_cap_tile(:,:)
! (land_pts, sm_levels, npft) The root extraction weights of the
! transpiration, before soil_evap adds the soil evaporation to layer 1.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: wt_root_tile(:,:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate f_vcmax_psi_state and cpl_state, filled with -1 ("not yet set"),
! unless they are already allocated.
! ---------------------------------------------------------------------
SUBROUTINE sf_stom_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF ( .NOT. ALLOCATED(f_vcmax_psi_state) ) THEN
  ALLOCATE( f_vcmax_psi_state(n_land_pts, npft) )
  f_vcmax_psi_state(:,:) = -1.0
END IF

IF ( .NOT. ALLOCATED(cpl_state) ) THEN
  ALLOCATE( cpl_state(n_land_pts, npft, n_cpl_state) )
  cpl_state(:,:,:) = -1.0
END IF

END SUBROUTINE sf_stom_state_alloc

! ---------------------------------------------------------------------
! Allocate the arrays handed to the surface flux step, with no
! transpiration to deliver, unless they are already allocated.
!
! MDK, 26/11/26
! ---------------------------------------------------------------------
SUBROUTINE sf_stom_flux_alloc( n_land_pts, n_levels )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts, n_levels

!$OMP CRITICAL (sf_stom_flux_alloc_crit)
IF ( .NOT. ALLOCATED(e_flux_tile) ) THEN
  ALLOCATE( gc_cap_tile(n_land_pts, npft) )
  ALLOCATE( wt_root_tile(n_land_pts, n_levels, npft) )
  gc_cap_tile(:,:)    = 0.0
  wt_root_tile(:,:,:) = 0.0
  ALLOCATE( e_flux_tile(n_land_pts, npft) )
  e_flux_tile(:,:) = -1.0
END IF
!$OMP END CRITICAL (sf_stom_flux_alloc_crit)

END SUBROUTINE sf_stom_flux_alloc

END MODULE sf_stom_state_mod
