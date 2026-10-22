! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE cuticular_leak_state_mod
! Plant water potentials carried between timesteps for the cuticular loss
! under the stomatal optimisation (l_som_cuticular_loss). Where the leak
! binds, the psi the impairment memory sees relaxes towards the leak's
! steady state through the stem capacitance (cap_stem), instead of jumping
! to it in one step; where it does not bind (or at night), the state is
! reset to the step's psi. As every night resets it, it is not dumped.
!
! MDK, 22/10/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: cut_psi_leaf, cut_psi_stem, cut_state_alloc, cut_psi_unset

REAL(KIND=real_jlslsm), PARAMETER :: cut_psi_unset = 1.0
      ! Value of a state not yet set (psi is <= 0).

! Leaf and stem water potential (Pa) of each PFT at each land point at the
! end of the last timestep.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_psi_leaf(:,:),              &
                                             cut_psi_stem(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state (n_land_pts, npft), unset, if not done.
! ---------------------------------------------------------------------
SUBROUTINE cut_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(cut_psi_leaf)) THEN
  ALLOCATE(cut_psi_leaf(n_land_pts, npft), cut_psi_stem(n_land_pts, npft))
  cut_psi_leaf(:,:) = cut_psi_unset
  cut_psi_stem(:,:) = cut_psi_unset
END IF

END SUBROUTINE cut_state_alloc

END MODULE cuticular_leak_state_mod
