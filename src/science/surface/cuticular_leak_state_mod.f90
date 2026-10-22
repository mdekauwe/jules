! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE cuticular_leak_state_mod
! Plant water store of the cuticular loss under the stomatal optimisation
! (l_som_cuticular_loss): the stem capacitance cap_stem * LAI, with its
! water potentials carried between timesteps. Where the leak binds, psi
! relaxes towards the leak's steady state and the store pays part of the
! leak; on every other step (and at night) the store relaxes towards the
! step's psi through the root path, refilling (or giving up) water. The
! store's flow, cut_q_store, is taken off the soil extraction in sf_evap,
! so the store's water is conserved. The state is not dumped: a restart
! starts it full (at the step's psi).
!
! MDK, 22/10/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: cut_psi_leaf, cut_psi_stem, cut_q_store, cut_c_store,              &
          cut_state_alloc, cut_psi_unset, cut_m_h2o, cut_store_return

REAL(KIND=real_jlslsm), PARAMETER :: cut_psi_unset = 1.0
      ! Value of a state not yet set (psi is <= 0).
REAL(KIND=real_jlslsm), PARAMETER :: cut_m_h2o = 0.018015
      ! Molar mass of water (kg mol-1).

! Leaf and stem water potential (Pa) of each PFT at each land point at the
! end of the last timestep.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_psi_leaf(:,:),              &
                                             cut_psi_stem(:,:)

! Flow out of the store of each PFT this timestep (kg m-2 s-1): > 0 the
! store gives water to the transpiration, < 0 it refills from the soil.
! The soil gives the transpiration less this.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_q_store(:,:)

! Capacitance of the store of each PFT (mol H2O m-2 Pa-1), cap_stem * LAI.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_c_store(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state (n_land_pts, npft), unset, if not done.
! ---------------------------------------------------------------------
SUBROUTINE cut_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(cut_psi_leaf)) THEN
  ALLOCATE(cut_psi_leaf(n_land_pts, npft), cut_psi_stem(n_land_pts, npft),    &
           cut_q_store(n_land_pts, npft), cut_c_store(n_land_pts, npft))
  cut_psi_leaf(:,:) = cut_psi_unset
  cut_psi_stem(:,:) = cut_psi_unset
  cut_q_store(:,:)  = 0.0
  cut_c_store(:,:)  = 0.0
END IF

END SUBROUTINE cut_state_alloc

! ---------------------------------------------------------------------
! sf_evap: the soil took dq (kg m-2 s-1) less than the store's flow out
! implied at point l, PFT ft (dq > 0: the store gave less, as the
! transpiration was smaller; dq < 0: the soil could not supply the
! refill). The store keeps that water, so its potentials move by
! dq dt / C.
! ---------------------------------------------------------------------
SUBROUTINE cut_store_return( l, ft, timestep, dq )

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, dq

REAL(KIND=real_jlslsm) :: dpsi

IF ( cut_c_store(l,ft) <= 0.0 .OR. cut_psi_leaf(l,ft) > 0.0 ) RETURN
dpsi = dq / cut_m_h2o * timestep / cut_c_store(l,ft)
cut_psi_leaf(l,ft) = MIN(cut_psi_leaf(l,ft) + dpsi, 0.0)
IF ( cut_psi_stem(l,ft) <= 0.0 )                                               &
  cut_psi_stem(l,ft) = MIN(cut_psi_stem(l,ft) + dpsi, 0.0)
cut_q_store(l,ft) = cut_q_store(l,ft) - dq

END SUBROUTINE cut_store_return

END MODULE cuticular_leak_state_mod
