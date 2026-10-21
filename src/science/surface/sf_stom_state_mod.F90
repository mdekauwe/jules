! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE
MODULE sf_stom_state_mod
! State that sf_stom carries from one timestep to the next, kept here (not
! private to sf_stom) so it can go in the dump and a restart reproduces a
! continuous run:
!   f_vcmax_state : soil-water capacity factor of l_som_vcmax_psi after the
!                   recovery lag (dump: fvcmax_psi_state);
!   cpl_state     : coupled transpiration (l_som_coupled_e), two-leaf: the
!                   last sunlit and shaded stomatal conductances and leaf
!                   boundary layer conductances (m/s) and canopy air
!                   resistance (s/m), the first guess of a timestep's first
!                   pass (dump: cpl_gl_sun, cpl_gl_shd, cpl_gbw_sun,
!                   cpl_gbw_shd, cpl_ra_ca).
! Both are < 0 until first set ("not yet set"), which is also what a dump
! without them gives.
! MDK, 22/10/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: f_vcmax_state, cpl_state, sf_stom_state_alloc, n_cpl_state

INTEGER, PARAMETER :: n_cpl_state = 5
    ! gl_sun, gl_shd, gbw_sun, gbw_shd, ra_ca

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: f_vcmax_state(:,:)
    ! (land_pts, npft)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cpl_state(:,:,:)
    ! (land_pts, npft, n_cpl_state)

CONTAINS

! ---------------------------------------------------------------------
! Allocate both, "not yet set" (-1), if not done.
! ---------------------------------------------------------------------
SUBROUTINE sf_stom_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF ( .NOT. ALLOCATED(f_vcmax_state) ) THEN
  ALLOCATE( f_vcmax_state(n_land_pts, npft) )
  f_vcmax_state(:,:) = -1.0
END IF
IF ( .NOT. ALLOCATED(cpl_state) ) THEN
  ALLOCATE( cpl_state(n_land_pts, npft, n_cpl_state) )
  cpl_state(:,:,:) = -1.0
END IF

END SUBROUTINE sf_stom_state_alloc

END MODULE sf_stom_state_mod
