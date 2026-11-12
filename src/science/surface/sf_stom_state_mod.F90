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
!     conductances of the sunlit and shaded leaves (m/s) and the canopy
!     air resistance (s/m). They are the first guess on the first pass of
!     the next timestep (dump names cpl_gl_sun, cpl_gl_shd, cpl_gbw_sun,
!     cpl_gbw_shd, cpl_ra_ca).
!
! Both are negative until first set, meaning "not yet set". A dump without
! them gives the same.
!
! MDK, 06/11/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: f_vcmax_psi_state, cpl_state, sf_stom_state_alloc, n_cpl_state,      &
          shared_e_shd_prev, shared_amax_shd_prev, store_e_prev

! gl_sun, gl_shd, gbw_sun, gbw_shd, ra_ca
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

END MODULE sf_stom_state_mod
