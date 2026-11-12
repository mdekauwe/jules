! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE plant_water_store_mod
! Values kept between timesteps for the plant capacitance store: the
! water held in the stem, used by the stomatal optimisation with plant
! capacitance (l_som_capacitance) or with water loss through the leaf
! cuticle (l_som_cuticular_loss). This header describes the whole
! mechanism; sf_stom (store_cpl_q, store_cpl_eval, the leak block) and
! sf_evap point to it.
!
! The store's water is counted from the full store (water potential 0),
! so it is 0 or negative:
!
!   W = C psi_x      (mol H2O m-2 ground)
!
! where psi_x is the stem node's water potential (Pa, <= 0) and C the
! store's capacitance, the water it gives per unit fall in psi_x:
!
!   C = cap_stem L_ref      (mol H2O m-2 ground Pa-1)
!
! cap_stem is the capacitance per unit leaf area and L_ref the leaf area
! index, or the reference leaf area index with kplant_model = 1.
!
! Balance over a timestep dt, solved at the end of the step (backward
! Euler, the "_bal" names in sf_stom):
!
!   q = C (psi_x0 - psi_x) / dt      (mol H2O m-2 ground s-1)
!   E - q = F(psi_s, psi_x)
!   W_new = W - q dt
!
! with q the store's flow out, psi_x0 the stem node's psi at the start of
! the step, E the transpiration, psi_s the soil water potential and F the
! flow the soil-to-stem path (soil, root and stem segments) carries
! between them (store_cpl_q).
!
! What the store assumes:
!
!   1. It has no conductance of its own. It sits at the stem node, ends
!      each step at that node's psi, and drains and refills through the
!      soil-to-stem path, so it relaxes with the time constant
!      tau = C / k_x, with k_x the conductance of that path.
!   2. C does not depend on psi: the store gives the same water per MPa
!      down to the stem's critical-loss psi (psi_floor_stem).
!   3. The optimiser sees the store's flow q as fixed within its search,
!      so stored water costs nothing at the margin; its cost appears
!      through the lower psi_x of later steps.
!   4. The store never sends water back to the soil: q <= E.
!   5. The leaf segment is solved on its own curve from the stem node's
!      psi.
!
! Each timestep sf_stom sets q. Water lost through the cuticle when the
! stomata are shut also drains the store. The flow out, store_q, is taken
! off the water drawn from the soil in sf_evap, so the water in soil,
! plant and air adds up.
!
! Because W is carried and psi_x is worked out from it, a change in C moves
! psi_x and leaves the water unchanged:
!
!   - when C rises (growth), the new tissue starts full, W is kept;
!   - when C falls (leaves are lost), psi_x is kept, and the lost tissue
!     takes its share of the store's shortfall from full (store_w_shed), so
!     shedding leaves does not dry the plant.
!
! W, C and the flow are carried between timesteps and written to the
! dump.
!
! MDK, 12/11/26

USE um_types, ONLY: real_jlslsm
USE water_constants_mod, ONLY: Water_kg_per_mol

IMPLICIT NONE

PRIVATE
PUBLIC :: store_w, store_q, store_c, store_w_shed, leak_e, store_q_gap,        &
          store_w_xyl, store_q_xyl,                                            &
          store_state_alloc, store_unset, store_return, store_wt_root,         &
          store_w_soil, store_w_root_avail

! Value of the store's water not yet set (the water is <= 0) (mol m-2).
REAL(KIND=real_jlslsm), PARAMETER :: store_unset = 1.0

CHARACTER(LEN=*), PARAMETER :: ModuleName='PLANT_WATER_STORE_MOD'

! Water in the store of each PFT at each land point at the end of the last
! timestep, relative to the full store at psi = 0 (mol H2O m-2 ground,
! <= 0; > 0: not yet set):
!
!   W = C psi_x.
!
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_w(:,:)

! Flow out of the store of each PFT this timestep (kg H2O m-2 ground s-1):
!
!   > 0  the store gives water to the transpiration,
!   < 0  it refills from the soil.
!
! The soil gives the transpiration less this.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_q(:,:)

! Capacitance of the store of each PFT (mol H2O m-2 ground Pa-1; 0: not yet
! set):
!
!   C = cap_stem L_ref,
!
! carried between timesteps (dumped) to apply the rule for a change of C.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_c(:,:)
! Water in the stem's conducting xylem of each PFT (kg H2O m-2 ground; < 0:
! not yet set), with the water released by cavitation
! (l_ximpair_cav_water). It is part of the plant's water: the conduits
! that embolise give their water to the store, and the conduits that are
! refilled (repair, new sapwood) take it back from the soil. sf_stom keeps
! it at
!
!   W_xyl = V_w (1 - L)
!
! with V_w the water of the sapwood's conduits when all are full and L the
! stem's locked-in loss, except for water released while the store was
! full, which is held here until the store has room.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_w_xyl(:,:)
! Water drawn from the soil this timestep to fill new xylem (kg H2O m-2
! ground s-1, >= 0; l_ximpair_cav_water): new rings and new sapwood start
! with their conduits full. sf_evap takes it from the soil with the
! store's refill; what the soil cannot give is taken from the store
! (store_return). Set each timestep by sf_stom and not carried over.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_q_xyl(:,:)

! Water deficit carried away by shed tissue this timestep (kg H2O m-2
! ground, >= 0): the store's water rises by this with no flow, as C falls
! at fixed psi. A diagnostic for the water budget.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_w_shed(:,:)

! Transpiration through the cuticular leak this timestep, as sf_stom
! estimates it (kg H2O m-2 ground s-1): the canopy flux at points where the
! leak sets the conductance (the stomata shut, l_som_cuticular_loss), after
! any cut that keeps the leaf above the critical-loss psi; 0 elsewhere. The
! output leak_e_pft is the tile's final transpiration at these points
! (sf_evap corrects the flux for the soil supply and the canopy coupling).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: leak_e(:,:)

! How far the store's flow the optimiser saw in its last pass (held fixed
! there, assumption 3) is from the flow at the committed E this timestep,
! relative to that E (-; 0 where the canopy does not transpire):
!
!   |q_commit - q_last pass| / E.
!
! A diagnostic of the passes' fixed point (output store_q_gap_pft).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_q_gap(:,:)

! Root extraction weights of each PFT at each land point this timestep
! (fraction of the root uptake from each soil layer, summing to 1), set in
! physiol before the soil evaporation is added to layer 1. sf_evap takes
! the store's flow (store_q) from, or gives it back to, the soil layers
! with these weights. Dimensions (land points, soil layers, PFTs).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_wt_root(:,:,:)

! Unfrozen water in each soil layer at each land point this timestep
! (kg H2O m-2), set in physiol: sf_evap does not let the store's refill
! take more than this from a layer. Dimensions (land points, soil layers).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_w_soil(:,:)

! Water the roots of each PFT can take this timestep (kg H2O m-2): the
! water above the closure point (v_close) in each soil layer, weighted by
! the root extraction weights (store_wt_root), set in physiol. sf_evap does
! not let the transpiration and the store's refill together take more than
! this. Dimensions (land points, PFTs).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: store_w_root_avail(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the store's arrays (n_land_pts, npft) unless they are already
! allocated. The water starts as "not yet set" (store_unset); the flows
! and capacitance start at 0. The root weights start with all uptake from
! layer 1 and the soil water at 0; physiol sets them each timestep before
! sf_evap uses them.
!
! MDK, 12/11/26
! ---------------------------------------------------------------------
SUBROUTINE store_state_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft
USE jules_soil_mod, ONLY: sm_levels

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points.
INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STORE_STATE_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(store_w)) THEN
  ALLOCATE(store_w(n_land_pts, npft), store_q(n_land_pts, npft),               &
           store_c(n_land_pts, npft), store_w_shed(n_land_pts, npft),          &
           leak_e(n_land_pts, npft), store_q_gap(n_land_pts, npft),            &
           store_w_xyl(n_land_pts, npft), store_q_xyl(n_land_pts, npft))
  store_w(:,:)      = store_unset
  store_q(:,:)      = 0.0
  store_c(:,:)      = 0.0
  store_w_shed(:,:) = 0.0
  leak_e(:,:)       = 0.0
  store_q_gap(:,:)  = 0.0
  store_w_xyl(:,:)  = -1.0
  store_q_xyl(:,:)  = 0.0
END IF

IF (.NOT. ALLOCATED(store_wt_root)) THEN
  ALLOCATE(store_wt_root(n_land_pts, sm_levels, npft),                         &
           store_w_soil(n_land_pts, sm_levels),                                &
           store_w_root_avail(n_land_pts, npft))
  store_wt_root(:,:,:)    = 0.0
  store_wt_root(:,1,:)    = 1.0
  store_w_soil(:,:)       = 0.0
  store_w_root_avail(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE store_state_alloc

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
! (kg mol-1). The store's water is kept no wetter than full (0); with the
! water of the stem's conduits (store_w_xyl set) the water above full is
! kept there, as water released by the conduits that the store has no
! room for yet, so none is lost. The flow out (store_q) is reduced by dq.
! Nothing is done where the store is not yet set.
!
! MDK, 12/11/26
! ---------------------------------------------------------------------
SUBROUTINE store_return( l, ft, timestep, dq )

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

! Nothing to do where the store is not set.
IF ( store_c(l,ft) <= 0.0 .OR. store_w(l,ft) > 0.0 ) RETURN

dw = dq / Water_kg_per_mol * timestep

! The water above full goes to the stem's conduits when they are tracked.
IF ( ALLOCATED(store_w_xyl) ) THEN
  IF ( store_w_xyl(l,ft) >= 0.0 .AND. store_w(l,ft) + dw > 0.0 )               &
    store_w_xyl(l,ft) = store_w_xyl(l,ft)                                      &
                        + ( store_w(l,ft) + dw ) * Water_kg_per_mol
END IF

! The store's water, no wetter than full (psi = 0).
store_w(l,ft) = MIN(store_w(l,ft) + dw, 0.0)
store_q(l,ft) = store_q(l,ft) - dq

END SUBROUTINE store_return

END MODULE plant_water_store_mod
