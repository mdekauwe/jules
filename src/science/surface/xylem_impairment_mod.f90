! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_MOD'

PRIVATE

PUBLIC :: leaf_conductance_impaired_jls,                                       &
          psi_leaf_impaired_jls,                                               &
          update_xylem_impairment

CONTAINS

! *********************************************************************
! Xylem impairment: drought damage to the plant's water-conducting tissue
! (xylem) that lasts beyond the drought. These routines pick the
! impairment model chosen for each PFT (pft_xylem_impairment_model) and
! pass the work to it:
!
!   0  none: the plant always uses its undamaged vulnerability curve;
!   1  kmax-impaired: damage lowers the plant's kmax, and the curve's b and
!      c are refitted to the lowered kmax (xylem_impairment_kmax_mod);
!   2  trunk-impaired: each part of the xylem keeps the conductance of the
!      driest water potential it has met;
!   3  memory-impaired: the conductance lost in drought is held as a cap on the
!      curve and recovers through new wood and new leaves.
!
! During spin-up no damage is applied or recorded.
!
! The vulnerability curve gives the xylem conductance k (mol m-2 s-1
! Pa-1) at water potential psi (Pa, below 0, more negative when drier).
! It falls from kmax in wet soil, with the shape set by conductance_b (Pa,
! below 0) and conductance_c (no units). kcrit is the lowest conductance
! the model allows. The arguments are:
!
!   kmax_ref: the undamaged kmax at each land point, the PFT value with
!     any canopy scaling the caller applies (for example times LAI);
!   kmax: the kmax after damage, with the same scaling;
!   kcrit: on the same scaling as kmax_ref.
!
! The whole trunk and memory models leave conductance_b and conductance_c
! as the PFT's; the kmax-impaired model gives each point refitted values.
!
! MDK, 21/11/26
! *********************************************************************

! ---------------------------------------------------------------------
! Work out the conductance (mol m-2 s-1 Pa-1) at a given leaf water
! potential of a plant whose xylem may have been damaged by drought, using
! the PFT's impairment model. The stomatal optimiser uses it as the
! reference for the hydraulic cost.
!
! With no impairment model, or during spin-up, this is the undamaged
! curve with kmax_ref. Otherwise the whole trunk or memory model gives
! the damaged value.
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_impaired_jls( pft,                                 &
                                          land_pnts,                           &
                                          water_potential,                     &
                                          kmax_ref,                            &
                                          kmax,                                &
                                          kcrit,                               &
                                          conductance_b,                       &
                                          conductance_c,                       &
                                          psi_leaf_extreme,                    &
                                        ! INTENT OUT
                                          leaf_conductance                     &
  )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_kmax,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls
USE xylem_impairment_whole_trunk_mod,                                          &
    ONLY: leaf_conductance_impaired_whole_trunk_jls
USE xylem_impairment_memory_mod, ONLY: leaf_conductance_impaired_memory_jls

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, land_pnts
                            ! Number of land points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(land_pnts)                                                   &
                            ! Water potential for each land point (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance for each
                            ! land point (mol m-2 s-1 Pa-1)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_conductance(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: errcode

REAL(KIND=real_jlslsm) :: b_ref(land_pnts), c_ref(land_pnts)
                            ! Unimpaired PFT curve parameters.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CONDUCTANCE_IMPAIRED_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and use the default
! conductance model when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  b_ref(:) = conductance_b_pft(pft)
  c_ref(:) = conductance_c_pft(pft)

  CALL leaf_conductance_jls( pft,                                              &
                             land_pnts,                                        &
                             water_potential,                                  &
                             kmax_ref,                                         &
                             kcrit,                                            &
                             b_ref,                                            &
                             c_ref,                                            &
                           ! INTENT OUT
                             leaf_conductance                                  &
                             )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE ( pft_xylem_impairment_model(pft) )

CASE ( xylem_impairment_kmax )
  ! The impairment is in the point's kmax, b and c: the default conductance
  ! model on the impaired curve.
  CALL leaf_conductance_jls( pft,                                              &
                             land_pnts,                                        &
                             water_potential,                                  &
                             kmax,                                             &
                             kcrit,                                            &
                             conductance_b,                                    &
                             conductance_c,                                    &
                           ! INTENT OUT
                             leaf_conductance                                  &
                             )

CASE ( xylem_impairment_whole_trunk )
  CALL leaf_conductance_impaired_whole_trunk_jls( pft,                         &
                                                 land_pnts,                    &
                                                 water_potential,              &
                                                 kmax,                         &
                                                 kcrit,                        &
                                                 conductance_b,                &
                                                 conductance_c,                &
                                                 psi_leaf_extreme,             &
                                               ! INTENT OUT
                                                 leaf_conductance              &
                                                 )

CASE ( xylem_impairment_memory )
  CALL leaf_conductance_impaired_memory_jls( pft,                              &
                                             land_pnts,                        &
                                             water_potential,                  &
                                             kmax_ref,                         &
                                             kmax,                             &
                                             kcrit,                            &
                                             conductance_b,                    &
                                             conductance_c,                    &
                                           ! INTENT OUT
                                             leaf_conductance                  &
                                             )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
    'pft_xylem_impairment_model should be none (0), kmax (1), ' //             &
    'whole_trunk (2) or memory (3)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_jls

! ---------------------------------------------------------------------
! For each candidate transpiration rate E of the stomatal optimiser, work out
! the leaf water potential (Pa) and the xylem conductance at the leaf
! (mol m-2 s-1 Pa-1) of a plant whose xylem may have been damaged by
! drought, at points where the stomata are open. The optimiser tries
! several values of E per point.
!
! With no impairment model, or during spin-up, the undamaged plant is
! solved (psi_leaf_jls). Otherwise the whole trunk or memory model solves
! the path from the root zone to the leaf.
!
! With l_som_rhizo_series (always on where this is called), a soil-to-root
! conductance K_s (mol m-2 s-1 Pa-1) lies in series ahead of the damaged
! plant. Water first crosses the soil, so the water potential where it
! enters the roots is
!
!   psi_in = psi_soil - E / K_s,
!
! where psi_soil is the root-zone water potential (Pa), and the damaged
! plant is solved from psi_in up to the leaf. The soil is not damaged.
! The conductance returned is the change in E for a change in leaf water
! potential along the whole path,
!
!   k_path = k_leaf / (1 + k_in / K_s),
!
! where k_leaf is the damaged plant's conductance at the leaf and k_in is
! its conductance where water enters the roots.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_impaired_jls( pft,                                         &
                                  n_e_leaf,                                    &
                                  land_pnts,                                   &
                                  open_pnts,                                   &
                                  veg_index,                                   &
                                  open_index,                                  &
                                  e_leaf,                                      &
                                  root_zone_psi,                               &
                                  kmax_ref,                                    &
                                  kmax,                                        &
                                  kcrit,                                       &
                                  conductance_b,                               &
                                  conductance_c,                               &
                                  psi_leaf_extreme,                            &
                                  psi_root_extreme,                            &
                               ! INTENT OUT
                                  leaf_psi,                                    &
                                  leaf_k                                       &
  )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_kmax,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls, leaf_conductance_jls
USE xylem_impairment_whole_trunk_mod, ONLY: psi_leaf_impaired_whole_trunk
USE xylem_impairment_memory_mod, ONLY: psi_leaf_impaired_memory
USE plant_hydraulic_path_mod, ONLY: ksr_path, som_psi_in_min

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_e_leaf                                                                     &
                            ! Number of transpirations per land point
, land_pnts                                                                    &
                            ! Number of land points
, open_pnts                                                                    &
                            ! Number of open stomata
, veg_index(land_pnts)                                                         &
                            ! Index of vegetation points on the land grid
, open_index(land_pnts)
                            ! Index of open stomata for each land point

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
, root_zone_psi(land_pnts)                                                     &
                             ! Water potential in the root zone (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance for each
                            ! land point (mol m-2 s-1 Pa-1)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point
                            ! (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)                                                  &
                            ! Minimum leaf water potential for each land point
                            ! (Pa)
, psi_root_extreme(land_pnts)
                            ! Minimum root water potential for each land point
                            ! (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: errcode, i, j, l

REAL(KIND=real_jlslsm) :: b_ref(land_pnts), c_ref(land_pnts)
                            ! Unimpaired PFT curve parameters.
REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(land_pnts), k_inlet(land_pnts)                                        &
                            ! l_som_rhizo_series: root water potential (Pa)
                            ! and impaired conductance there, one sample.
, e1(1, open_pnts), psi1(1, open_pnts), k1(1, open_pnts)                       &
                            ! One sample per open point.
, k_s                                                                          &
                            ! Soil-to-root conductance of the path.
, psi_k(land_pnts)
                            ! Whole trunk: root-end potential the root-end
                            ! conductance is taken at (Pa).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_IMPAIRED_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and just call the base
! conductance models when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  b_ref(:) = conductance_b_pft(pft)
  c_ref(:) = conductance_c_pft(pft)

  CALL psi_leaf_jls( pft,                                                      &
                     n_e_leaf,                                                 &
                     land_pnts,                                                &
                     open_pnts,                                                &
                     veg_index,                                                &
                     open_index,                                               &
                     e_leaf,                                                   &
                     root_zone_psi,                                            &
                     kmax_ref,                                                 &
                     kcrit,                                                    &
                     b_ref,                                                    &
                     c_ref,                                                    &
                   ! INTENT OUT
                     leaf_psi,                                                 &
                     leaf_k                                                    &
            )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE ( pft_xylem_impairment_model(pft) )

  CASE (xylem_impairment_kmax)
    ! The impairment is in the point's kmax, b and c: the default solver
    ! (with its soil-to-root link) on the impaired curve.
    CALL psi_leaf_jls( pft,                                                    &
                       n_e_leaf,                                               &
                       land_pnts,                                              &
                       open_pnts,                                              &
                       veg_index,                                              &
                       open_index,                                             &
                       e_leaf,                                                 &
                       root_zone_psi,                                          &
                       kmax,                                                   &
                       kcrit,                                                  &
                       conductance_b,                                          &
                       conductance_c,                                          &
                     ! INTENT OUT
                       leaf_psi,                                               &
                       leaf_k                                                  &
              )

  CASE (xylem_impairment_whole_trunk, xylem_impairment_memory)
    !-----------------------------------------------------------------------
    ! Soil-to-root conductance in series ahead of the impaired plant, as
    ! psi_leaf_jls: the root water potential of each sample is
    !
    !   psi_in = psi_src - E/K_s,
    !
    ! the impaired plant path runs from there, and the whole-path marginal
    ! conductance is
    !
    !   k_p / (1 + k(psi_in) / K_s),
    !
    ! with k(psi_in) this model's impaired conductance at the root. The soil
    ! is not damaged.
    !-----------------------------------------------------------------------
    DO i = 1, n_e_leaf
      psi_in(:) = root_zone_psi(:)
      DO j = 1, open_pnts
        l = veg_index(open_index(j))
        k_s = ksr_path(l)
        IF ( ksr_path(l) < 0.0 ) THEN
          ! No soil link for this PFT (fsmc_mod /= 2).
        ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
          psi_in(l) = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s,                &
                          som_psi_in_min)
        ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
          psi_in(l) = som_psi_in_min
        END IF
        e1(1,j) = e_leaf(i,j)
      END DO
      CALL impaired_path( 1, e1, psi_in, psi1, k1 )
      IF ( pft_xylem_impairment_model(pft) ==                                  &
           xylem_impairment_whole_trunk ) THEN
        ! Whole trunk: the root end of the path is where water enters, so its
        ! conductance is the unimpaired curve at the root's historic
        ! minimum (as psi_leaf_impaired_whole_trunk uses there), not the
        ! leaf's, which leaf_conductance_impaired_jls would use.
        psi_k(:) = MIN(psi_in(:), psi_root_extreme(:))
        CALL leaf_conductance_jls( pft, land_pnts, psi_k, kmax, kcrit,         &
                                   conductance_b, conductance_c, k_inlet )
      ELSE
        CALL leaf_conductance_impaired_jls( pft, land_pnts, psi_in,            &
                                            kmax_ref, kmax, kcrit,             &
                                            conductance_b, conductance_c,      &
                                            psi_leaf_extreme, k_inlet )
      END IF
      DO j = 1, open_pnts
        l = veg_index(open_index(j))
        k_s = ksr_path(l)
        leaf_psi(i,j) = psi1(1,j)
        IF ( ksr_path(l) < 0.0 ) THEN
          leaf_k(i,j) = k1(1,j)
        ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
          leaf_k(i,j) = k1(1,j) / (1.0 + k_inlet(l) / k_s)
        ELSE
          leaf_k(i,j) = 0.0
        END IF
      END DO
    END DO
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), kmax (1), ' //            &
     'whole_trunk (2) or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! Xylem impairment: solve the damaged plant from the given root water
! potential psi_root (Pa) to the leaf with the whole trunk or memory
! model, giving the leaf water potential and conductance for each E.
SUBROUTINE impaired_path( n_e, e_in, psi_root, psi_out, k_out )

INTEGER, INTENT(IN) :: n_e
REAL(KIND=real_jlslsm), INTENT(IN)  :: e_in(n_e, open_pnts),                   &
                                       psi_root(land_pnts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi_out(n_e, open_pnts),                &
                                       k_out(n_e, open_pnts)

SELECT CASE ( pft_xylem_impairment_model(pft) )
CASE (xylem_impairment_whole_trunk)
    CALL psi_leaf_impaired_whole_trunk( pft,                                   &
                                        n_e,                                   &
                                        land_pnts,                             &
                                        open_pnts,                             &
                                        veg_index,                             &
                                        open_index,                            &
                                        e_in,                                  &
                                        psi_root,                              &
                                        kmax,                                  &
                                        kcrit,                                 &
                                        conductance_b,                         &
                                        conductance_c,                         &
                                        psi_leaf_extreme,                      &
                                        psi_root_extreme,                      &
                                      ! INTENT OUT
                                        psi_out,                               &
                                        k_out                                  &
                                        )

CASE (xylem_impairment_memory)
    CALL psi_leaf_impaired_memory( pft,                                        &
                                   n_e,                                        &
                                   land_pnts,                                  &
                                   open_pnts,                                  &
                                   veg_index,                                  &
                                   open_index,                                 &
                                   e_in,                                       &
                                   psi_root,                                   &
                                   kmax_ref,                                   &
                                   kmax,                                       &
                                   kcrit,                                      &
                                   conductance_b,                              &
                                   conductance_c,                              &
                                 ! INTENT OUT
                                   psi_out,                                    &
                                   k_out                                       &
                                   )

END SELECT

END SUBROUTINE impaired_path

END SUBROUTINE psi_leaf_impaired_jls

! ---------------------------------------------------------------------
! Xylem impairment: record this step's drought damage with the PFT's
! impairment model, once the step's water potentials are known.
!
! The whole trunk model updates the driest leaf and root water potentials
! met so far. The memory model updates its conductance cap (recovery, then
! new damage); with per-segment memory and the two-leaf canopy's leaf
! classes (psi_leaf_cls, w_cls) the leaf segment is damaged by the
! classes' weighted mean conductance. Nothing is recorded during spin-up
! or with no impairment model. The shape of the vulnerability curve is
! chosen elsewhere.
!
! With l_ximpair_leafless_root_psi, a PFT with the budburst gate (phen_bb_doy
! or phen_bb_gdd set) is dormant (phen_dormant = 1) from the step the cold
! leaf drop takes its canopy to the phenology floor, after this year's
! budburst, until the next budburst. While it is dormant the root water
! potential drives both models, and the leaf classes of the two-leaf
! canopy, in place of the leaf and stem, so a leafless tree does not
! embolise from the potential of the few leaves left at the floor. A
! dormant PFT on a frozen root zone takes no damage, as any other, from
! the frozen-soil hold of the memory model (l_ximpair_frozen_hold, which
! l_ximpair_leafless_root_psi needs). The whole trunk model has no frozen test:
! with it the dormant root water potential replaces the leaf potential only
! and its root extreme is unchanged. Recovery is not affected.
! The gate alone does not make a PFT dormant, a PFT whose canopy never
! falls to the floor is never dormant, and leaves shed by drought or heat
! do not make a PFT dormant. The stem water potential kept for output
! (ximpair_psi_stem, psi_stem_pft and the per-segment PLC) is the stem's
! own, also while dormant.
!
! MDK, 17/11/26
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment ( n_land_pts                                &
,                                    n_open_pts                                &
,                                    open_index                                &
,                                    pft                                       &
,                                    psi_leaf                                  &
,                                    psi_root                                  &
,                                    psi_stem                                  &
,                                    impaired_k_max                            &
,                                    psi_leaf_extreme                          &
,                                    psi_root_extreme                          &
,                                    lai                                       &
,                                    canht                                     &
,                                    anetc                                     &
,                                    psi_leaf_cls                              &
,                                    w_cls                                     &
                                   )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_kmax,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model, ximpair_lai_min
USE jules_vegetation_mod, ONLY: l_ximpair_leafless_root_psi
USE phen_budburst_mod, ONLY: phen_bb_on, phen_dormant
USE xylem_impairment_seg_state_mod, ONLY: ximpair_psi_stem

USE xylem_impairment_kmax_mod, ONLY: update_xylem_impairment_kmax
USE xylem_impairment_whole_trunk_mod, ONLY: update_xylem_impairment_whole_trunk
USE xylem_impairment_memory_mod, ONLY: update_xylem_impairment_memory

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  n_land_pts                                                                   &
                            ! Number of land points
, n_open_pts                                                                   &
                            ! Number of land points with open stomata
, open_index(n_open_pts)                                                       &
                            ! Land point of each point with open stomata
, pft
                            ! Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN), TARGET ::                                  &
  psi_leaf(n_land_pts)                                                         &
                            ! Leaf water potential (Pa)
, psi_root(n_land_pts)                                                         &
                            ! Root water potential (Pa)
, psi_stem(n_land_pts)
                            ! Stem water potential (Pa; memory model driver 4)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
  impaired_k_max(n_land_pts)                                                   &
                            ! Maximum xylem conductance (mol m-2 s-1 Pa-1)
,  psi_leaf_extreme(n_land_pts)                                                &
                            ! Minimum historic leaf water potential (Pa)
, psi_root_extreme(n_land_pts)
                            ! Minimum historic root water potential (Pa)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  lai(n_land_pts)                                                              &
                            ! Leaf area index
, canht(n_land_pts)                                                            &
                            ! Canopy height (m)
, anetc(n_land_pts)
                            ! Canopy net photosynthesis (mol CO2 m-2 s-1)
                            ! (canht and anetc are not used; lai only with
                            ! ximpair_lai_min.)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL, TARGET ::                        &
  psi_leaf_cls(n_land_pts,2)                                                   &
                            ! Leaf water potential of the sunlit (1) and
                            ! shaded (2) leaves of the two-leaf canopy (Pa;
                            ! memory model, leaf segment only)
, w_cls(n_land_pts,2)
                            ! Their shares of the leaf area (sum 1)

INTEGER :: errcode

! Open point counter and land point.
INTEGER :: j, l

! The water potentials that drive the impairment (Pa). They point to the
! arguments, or to copies holding the root water potential while the PFT
! is dormant.
REAL(KIND=real_jlslsm), POINTER :: psi_leaf_d(:), psi_stem_d(:), psi_cls_d(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, TARGET :: leaf_dormant(:),                &
                                               stem_dormant(:),                &
                                               cls_dormant(:,:)

! This PFT's xylem is at the root water potential while dormant.
LOGICAL :: l_dormant_root

! The PFT has a leaf area threshold, and this point is leafless.
LOGICAL :: l_lai_gate, l_leafless

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not update the xylem impairment model during spinup
IF (is_spinup) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! A dormant PFT (l_ximpair_leafless_root_psi, budburst gate, canopy fallen to
! the floor): its xylem is taken to be at the root water potential.
l_dormant_root = l_ximpair_leafless_root_psi .AND. phen_bb_on(pft)
IF (l_dormant_root) l_dormant_root = ALLOCATED(phen_dormant)

! A PFT with too little leaf area (ximpair_lai_min) is treated the same.
l_lai_gate = ximpair_lai_min(pft) > 0.0
l_dormant_root = l_dormant_root .OR. l_lai_gate

NULLIFY(psi_cls_d)

IF (l_dormant_root) THEN
  ALLOCATE(leaf_dormant(n_land_pts), stem_dormant(n_land_pts))
  leaf_dormant(:) = psi_leaf(:)
  stem_dormant(:) = psi_stem(:)

  IF (PRESENT(psi_leaf_cls)) THEN
    ALLOCATE(cls_dormant(n_land_pts,2))
    cls_dormant(:,:) = psi_leaf_cls(:,:)
  END IF

  DO j = 1,n_open_pts
    l = open_index(j)

    l_leafless = l_lai_gate
    IF (l_leafless) l_leafless = lai(l) < ximpair_lai_min(pft)
    IF (.NOT. l_leafless .AND. ALLOCATED(phen_dormant) .AND.                   &
        l_ximpair_leafless_root_psi .AND. phen_bb_on(pft))                     &
      l_leafless = phen_dormant(l,pft) > 0.5

    IF (l_leafless) THEN
      leaf_dormant(l) = psi_root(l)
      stem_dormant(l) = psi_root(l)

      IF (ALLOCATED(cls_dormant)) cls_dormant(l,:) = psi_root(l)
    END IF

  END DO

  psi_leaf_d => leaf_dormant
  psi_stem_d => stem_dormant

  IF (ALLOCATED(cls_dormant)) psi_cls_d => cls_dormant
ELSE
  psi_leaf_d => psi_leaf
  psi_stem_d => psi_stem

  IF (PRESENT(psi_leaf_cls)) psi_cls_d => psi_leaf_cls
END IF

SELECT CASE (pft_xylem_impairment_model(pft))
  CASE(xylem_impairment_none)
    ! No xylem impairment

  CASE(xylem_impairment_kmax)
    CALL update_xylem_impairment_kmax( n_land_pts,                             &
                                       n_open_pts,                             &
                                       open_index,                             &
                                       pft,                                    &
                                       psi_leaf_d,                             &
                                       psi_root,                               &
                                       impaired_k_max                          &
                                     )

  CASE(xylem_impairment_whole_trunk)
    CALL update_xylem_impairment_whole_trunk( n_land_pts,                      &
                                              n_open_pts,                      &
                                              open_index,                      &
                                              pft,                             &
                                              psi_leaf_d,                      &
                                              psi_root,                        &
                                              psi_leaf_extreme,                &
                                              psi_root_extreme                 &
                                            )

  CASE(xylem_impairment_memory)
    IF (PRESENT(psi_leaf_cls)) THEN
      CALL update_xylem_impairment_memory( n_land_pts,                         &
                                           n_open_pts,                         &
                                           open_index,                         &
                                           pft,                                &
                                           psi_leaf_d,                         &
                                           psi_root,                           &
                                           psi_stem_d,                         &
                                           impaired_k_max,                     &
                                           psi_leaf_cls = psi_cls_d,           &
                                           w_cls = w_cls                       &
                                         )
    ELSE
      CALL update_xylem_impairment_memory( n_land_pts,                         &
                                           n_open_pts,                         &
                                           open_index,                         &
                                           pft,                                &
                                           psi_leaf_d,                         &
                                           psi_root,                           &
                                           psi_stem_d,                         &
                                           impaired_k_max,                     &
                                           w_cls = w_cls                       &
                                         )
    END IF

    ! The stem water potential kept for output is the stem's own.
    IF (l_dormant_root .AND. ALLOCATED(ximpair_psi_stem))                      &
      ximpair_psi_stem(:,pft) = psi_stem(:)

  CASE DEFAULT
    errcode = 1
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), kmax (1), ' //            &
     'whole_trunk (2) or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment

END MODULE xylem_impairment_mod
