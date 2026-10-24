! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE cuticular_leak_state_mod
! State of the stomatal optimisation's plant capacitance (l_som_capacitance,
! or l_som_cuticular_loss): a plant water store at the stem node.
!
! The store's water, relative to the full store (psi = 0),
!
!   W = C psi_x      (mol H2O m-2 ground, <= 0),
!
! is carried between timesteps and dumped, with its capacitance
!
!   C = cap_stem L_ref      (mol H2O m-2 ground Pa-1),
!
! (L_ref: LAI, or the reference LAI with kplant_model = 1, so that the
! store's time constant C / k does not change with phenology) and psi_stem.
!
! Each step (sf_stom) the store sets the transpiration, and the leak drains
! it where it binds; it refills from the soil through the root side. Its
! flow out, cut_q_store, is taken off the soil extraction in sf_evap, so
! soil + plant + atmosphere water is conserved. Because W, not psi, is
! carried, a change of C moves psi, not the water:
!
!   C rises (growth):    the new tissue is born full (W kept);
!   C falls (leaf loss): psi is kept and the shed tissue takes its share of
!                        the deficit (cut_w_shed), so shedding does not dry
!                        the plant.
!
! (The module and dump names keep "cut", from the cuticular loss, for
! which the store was first written.)
!
! MDK, 23/10/26

USE um_types, ONLY: real_jlslsm
USE water_constants_mod, ONLY: Water_kg_per_mol   ! MDK, 24/10/26

IMPLICIT NONE

PRIVATE
PUBLIC :: cut_w_store, cut_psi_stem, cut_q_store, cut_c_store, cut_w_shed,     &
          cut_state_alloc, cut_psi_unset, cut_m_h2o, cut_store_return

REAL(KIND=real_jlslsm), PARAMETER :: cut_psi_unset = 1.0
      ! Value of a state not yet set (psi and the store's water are <= 0)
      ! (Pa, or mol m-2).
REAL(KIND=real_jlslsm), PARAMETER :: cut_m_h2o = Water_kg_per_mol
      ! Molar mass of water (kg mol-1), JULES's value.
      !
      ! MDK, 24/10/26

CHARACTER(LEN=*), PARAMETER :: ModuleName='CUTICULAR_LEAK_STATE_MOD'

! Water in the store of each PFT at each land point at the end of the last
! timestep, relative to the full store at psi = 0 (mol H2O m-2 ground,
! <= 0; > 0: not yet set):
!
!   W = C psi_x.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_store(:,:)

! Stem water potential of each PFT at each land point at the end of the
! last timestep, with plant segments (Pa; > 0: not yet set).
!
! MDK, 23/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_psi_stem(:,:)

! Flow out of the store of each PFT this timestep (kg H2O m-2 ground s-1):
!
!   > 0  the store gives water to the transpiration,
!   < 0  it refills from the soil.
!
! The soil gives the transpiration less this.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_q_store(:,:)

! Capacitance of the store of each PFT (mol H2O m-2 ground Pa-1; 0: not yet
! set):
!
!   C = cap_stem L_ref,
!
! carried between timesteps (dumped) to apply the rule for a change of C.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_c_store(:,:)

! Water deficit carried away by shed tissue this timestep (kg H2O m-2
! ground, >= 0): the store's water rises by this with no flow, as C falls
! at fixed psi. A diagnostic for the water budget.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_shed(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state (n_land_pts, npft), not yet set, if not done.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
SUBROUTINE cut_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts
      ! Number of land points.
      !
      ! MDK, 23/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CUT_STATE_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(cut_w_store)) THEN
  ALLOCATE(cut_w_store(n_land_pts, npft), cut_psi_stem(n_land_pts, npft),      &
           cut_q_store(n_land_pts, npft), cut_c_store(n_land_pts, npft),       &
           cut_w_shed(n_land_pts, npft))
  cut_w_store(:,:)  = cut_psi_unset
  cut_psi_stem(:,:) = cut_psi_unset
  cut_q_store(:,:)  = 0.0
  cut_c_store(:,:)  = 0.0
  cut_w_shed(:,:)   = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cut_state_alloc

! ---------------------------------------------------------------------
! sf_evap: the store's flow out at point l, PFT ft, could not all be
! applied, by dq (kg H2O m-2 s-1):
!
!   dq > 0  the store gave less (the transpiration was smaller than its
!           flow out);
!   dq < 0  the soil could not supply the refill.
!
! The store keeps that water,
!
!   dW = dq dt / M_w      (mol m-2; no more than fills it to psi = 0),
!
! its stem psi moves by dW / C, and its flow out is reduced by dq.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
SUBROUTINE cut_store_return( l, ft, timestep, dq )

INTEGER, INTENT(IN) :: l
      ! Land point.
INTEGER, INTENT(IN) :: ft
      ! Plant functional type.
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep
      ! Timestep (s).
REAL(KIND=real_jlslsm), INTENT(IN) :: dq
      ! Flow out that could not be applied (kg H2O m-2 s-1).

REAL(KIND=real_jlslsm) :: dw
      ! Water the store keeps (mol H2O m-2).
REAL(KIND=real_jlslsm) :: dpsi
      ! The change of psi that gives (Pa).
      !
      ! MDK, 23/10/26

! Nothing to do where the store is not set.
IF ( cut_c_store(l,ft) <= 0.0 .OR. cut_w_store(l,ft) > 0.0 ) RETURN
dw   = dq / cut_m_h2o * timestep
dpsi = dw / cut_c_store(l,ft)
! The store's water and stem psi, no wetter than full (psi = 0).
cut_w_store(l,ft) = MIN(cut_w_store(l,ft) + dw, 0.0)
IF ( cut_psi_stem(l,ft) <= 0.0 )                                               &
  cut_psi_stem(l,ft) = MIN(cut_psi_stem(l,ft) + dpsi, 0.0)
cut_q_store(l,ft) = cut_q_store(l,ft) - dq

END SUBROUTINE cut_store_return

END MODULE cuticular_leak_state_mod
