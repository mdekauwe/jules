! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE cuticular_leak_state_mod
! Values kept between timesteps for the plant's stem water store, used by
! the stomatal optimisation with plant capacitance (l_som_capacitance) or
! with water loss through the leaf cuticle (l_som_cuticular_loss).
!
! The store holds water in the stem. Its water content is counted from the
! full store (water potential 0), so it is 0 or negative:
!
!   W = C psi_x      (mol H2O m-2 ground)
!
! where psi_x is the stem's water potential (Pa, <= 0) and C the store's
! capacitance, the water it gives per unit fall in psi_x:
!
!   C = cap_stem L_ref      (mol H2O m-2 ground Pa-1)
!
! cap_stem is the capacitance per unit leaf area and L_ref the leaf area
! index, or the reference leaf area index with kplant_model = 1, so the
! time the store takes to drain through the plant's conductance (C / k)
! does not change as leaves grow and fall. W, C and the stem's water
! potential are carried between timesteps and written to the dump.
!
! Each timestep sf_stom works out how much of the transpiration the store
! gives. Water lost through the cuticle when the stomata are shut also
! drains it. The store refills from the soil through the roots. Its flow
! out, cut_q_store, is taken off the water drawn from the soil in sf_evap,
! so the water in soil, plant and air adds up.
!
! Because W is carried and psi_x is worked out from it, a change in C moves
! psi_x and leaves the water unchanged:
!
!   - when C rises (growth), the new tissue starts full, W is kept;
!   - when C falls (leaves are lost), psi_x is kept, and the lost tissue
!     takes its share of the store's shortfall from full (cut_w_shed), so
!     shedding leaves does not dry the plant.
!
! The "cut" in the module and dump names comes from cuticular.
!
! MDK, 25/10/26

USE um_types, ONLY: real_jlslsm
USE water_constants_mod, ONLY: Water_kg_per_mol

IMPLICIT NONE

PRIVATE
PUBLIC :: cut_w_store, cut_psi_stem, cut_q_store, cut_c_store, cut_w_shed,     &
          cut_e_leak, cut_state_alloc, cut_psi_unset, cut_m_h2o,               &
          cut_store_return

! Value of a state not yet set (psi and the store's water are <= 0)
! (Pa, or mol m-2).
REAL(KIND=real_jlslsm), PARAMETER :: cut_psi_unset = 1.0
! Molar mass of water (kg mol-1), JULES's value.
REAL(KIND=real_jlslsm), PARAMETER :: cut_m_h2o = Water_kg_per_mol

CHARACTER(LEN=*), PARAMETER :: ModuleName='CUTICULAR_LEAK_STATE_MOD'

! Water in the store of each PFT at each land point at the end of the last
! timestep, relative to the full store at psi = 0 (mol H2O m-2 ground,
! <= 0; > 0: not yet set):
!
!   W = C psi_x.
!
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_store(:,:)

! Stem water potential of each PFT at each land point at the end of the
! last timestep, with plant segments (Pa; > 0: not yet set).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_psi_stem(:,:)

! Flow out of the store of each PFT this timestep (kg H2O m-2 ground s-1):
!
!   > 0  the store gives water to the transpiration,
!   < 0  it refills from the soil.
!
! The soil gives the transpiration less this.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_q_store(:,:)

! Capacitance of the store of each PFT (mol H2O m-2 ground Pa-1; 0: not yet
! set):
!
!   C = cap_stem L_ref,
!
! carried between timesteps (dumped) to apply the rule for a change of C.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_c_store(:,:)

! Water deficit carried away by shed tissue this timestep (kg H2O m-2
! ground, >= 0): the store's water rises by this with no flow, as C falls
! at fixed psi. A diagnostic for the water budget.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_w_shed(:,:)

! Transpiration through the leak this timestep, as sf_stom estimates it
! (kg H2O m-2 ground s-1): the canopy flux at points where the leak sets the
! conductance (the stomata shut, l_som_cuticular_loss), after any cut that
! keeps the leaf above the critical-loss psi; 0 elsewhere. The output
! cut_e_leak_pft is the tile's final transpiration at these points (sf_evap
! corrects the flux for the soil supply and the canopy coupling).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: cut_e_leak(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the store's arrays (n_land_pts, npft) unless they are already
! allocated. The water and stem water potential start as "not yet set"
! (cut_psi_unset); the flows and capacitance start at 0.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
SUBROUTINE cut_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points.
INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CUT_STATE_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(cut_w_store)) THEN
  ALLOCATE(cut_w_store(n_land_pts, npft), cut_psi_stem(n_land_pts, npft),      &
           cut_q_store(n_land_pts, npft), cut_c_store(n_land_pts, npft),       &
           cut_w_shed(n_land_pts, npft), cut_e_leak(n_land_pts, npft))
  cut_w_store(:,:)  = cut_psi_unset
  cut_psi_stem(:,:) = cut_psi_unset
  cut_q_store(:,:)  = 0.0
  cut_c_store(:,:)  = 0.0
  cut_w_shed(:,:)   = 0.0
  cut_e_leak(:,:)   = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cut_state_alloc

! ---------------------------------------------------------------------
! Put back into the plant capacitance store the part of its flow that
! could not be used this timestep, so the store's water stays in balance.
! sf_evap calls it at land point l and PFT ft. dq (kg H2O m-2 s-1) is the
! part not used:
!
!   dq > 0: the transpiration was smaller than the store's flow out, so
!           the store gave less;
!   dq < 0: the soil could not supply all of the store's refill.
!
! The store keeps that water:
!
!   dW = dq dt / M_w      (mol H2O m-2 ground)
!
! where dt is the timestep (s) and M_w the molar mass of water
! (kg mol-1). The store's water, and the stem's water potential, which
! moves by dW / C, are kept no wetter than full (0). The flow out
! (cut_q_store) is reduced by dq. Nothing is done where the store is not
! yet set.
!
! MDK, 23/10/26
! ---------------------------------------------------------------------
SUBROUTINE cut_store_return( l, ft, timestep, dq )

! Land point.
INTEGER, INTENT(IN) :: l
! Plant functional type.
INTEGER, INTENT(IN) :: ft
! Timestep (s).
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep
! Flow out that could not be applied (kg H2O m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: dq

! Water the store keeps (mol H2O m-2).
REAL(KIND=real_jlslsm) :: dw
! The change of psi that gives (Pa).
REAL(KIND=real_jlslsm) :: dpsi

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
