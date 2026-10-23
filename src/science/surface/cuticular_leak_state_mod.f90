! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE cuticular_leak_state_mod
! Plant water store of the cuticular loss under the stomatal optimisation
! (l_som_cuticular_loss): the stem capacitance C = cap_stem * LAI, with its
! water carried between timesteps. Where the leak binds, psi relaxes
! towards the leak's steady state and the store pays part of the leak; on
! every other step (and at night) the store relaxes towards the step's psi
! through the root path, refilling (or giving up) water, and at lit points
! psi_leaf is the store's psi. The store's flow, cut_q_store, is taken off
! the soil extraction in sf_evap, so the store's water is conserved. The
! water, not psi, is carried, so a change of LAI (and so of C) changes the
! store's psi, not its water. C is cap_stem times the leaf area that the
! plant conductance is scaled by (LAI, or the reference LAI with
! kplant_model = 1), so the store's time constant C / k does not change
! with phenology. When C grows the new tissue is born full (the water is
! kept); when it shrinks (leaf loss) psi is kept and the shed tissue takes
! its share of the deficit (cut_w_shed), so shedding does not dry the plant.
! The state (W, psi_stem and C) is dumped.
!
! MDK, 22/10/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: cut_w_store, cut_psi_stem, cut_q_store, cut_c_store, cut_w_shed,   &
          cut_state_alloc, cut_psi_unset, cut_m_h2o, cut_store_return

REAL(KIND=real_jlslsm), PARAMETER :: cut_psi_unset = 1.0
      ! Value of a state not yet set (psi and the store's water are <= 0).
REAL(KIND=real_jlslsm), PARAMETER :: cut_m_h2o = 0.018015
      ! Molar mass of water (kg mol-1).

! Water in the store of each PFT at each land point at the end of the last
! timestep, relative to the store at psi = 0 (mol H2O m-2, <= 0): C times
! the store's (leaf) psi.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_store(:,:)

! Stem water potential (Pa) of each PFT at each land point at the end of
! the last timestep (with plant segments).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_psi_stem(:,:)

! Flow out of the store of each PFT this timestep (kg m-2 s-1): > 0 the
! store gives water to the transpiration, < 0 it refills from the soil.
! The soil gives the transpiration less this.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_q_store(:,:)

! Capacitance of the store of each PFT (mol H2O m-2 Pa-1), cap_stem times
! the reference leaf area of the plant conductance; carried between
! timesteps (dumped) to apply the rule for a change of C.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_c_store(:,:)

! Water deficit carried away by shed tissue this timestep (kg m-2, >= 0):
! the store's water rises by this with no flow, as C falls at fixed psi.
! A diagnostic for the water budget.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_shed(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state (n_land_pts, npft), unset, if not done.
! ---------------------------------------------------------------------
SUBROUTINE cut_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(cut_w_store)) THEN
  ALLOCATE(cut_w_store(n_land_pts, npft), cut_psi_stem(n_land_pts, npft),     &
           cut_q_store(n_land_pts, npft), cut_c_store(n_land_pts, npft),      &
           cut_w_shed(n_land_pts, npft))
  cut_w_store(:,:)  = cut_psi_unset
  cut_psi_stem(:,:) = cut_psi_unset
  cut_q_store(:,:)  = 0.0
  cut_c_store(:,:)  = 0.0
  cut_w_shed(:,:)   = 0.0
END IF

END SUBROUTINE cut_state_alloc

! ---------------------------------------------------------------------
! sf_evap: the soil took dq (kg m-2 s-1) less than the store's flow out
! implied at point l, PFT ft (dq > 0: the store gave less, as the
! transpiration was smaller; dq < 0: the soil could not supply the
! refill). The store keeps that water, dq dt (no more than fills it to
! psi = 0), and its stem psi moves by dq dt / C.
!
! MDK, 22/10/26
! ---------------------------------------------------------------------
SUBROUTINE cut_store_return( l, ft, timestep, dq )

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, dq

REAL(KIND=real_jlslsm) :: dw, dpsi

IF ( cut_c_store(l,ft) <= 0.0 .OR. cut_w_store(l,ft) > 0.0 ) RETURN
dw   = dq / cut_m_h2o * timestep
dpsi = dw / cut_c_store(l,ft)
cut_w_store(l,ft) = MIN(cut_w_store(l,ft) + dw, 0.0)
IF ( cut_psi_stem(l,ft) <= 0.0 )                                               &
  cut_psi_stem(l,ft) = MIN(cut_psi_stem(l,ft) + dpsi, 0.0)
cut_q_store(l,ft) = cut_q_store(l,ft) - dq

END SUBROUTINE cut_store_return

END MODULE cuticular_leak_state_mod
