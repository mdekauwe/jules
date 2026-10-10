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

PUBLIC :: leaf_conductance_impaired_jls,                                       &
          xylem_conductance_impaired_stom_opt_jls,                             &
          leaf_psi_impaired_jls,                                               &
          update_xylem_impairment,                                             &
          canopy_impaired_psi_jls

CONTAINS

! *********************************************************************
! Contains routines used to switch between different xylem conductance
! models.
!
! NOTE on the curve parameter arguments passed to the routines below:
!   kmax_ref      : unimpaired maximum xylem conductance for each land
!                   point, i.e. kmax_pft(pft) with whatever canopy
!                   scaling the caller applies (e.g. * fpar for the
!                   big-leaf stomatal optimisation). Used, together with
!                   conductance_b_pft/conductance_c_pft, whenever no
!                   impairment is applied (spinup, xylem_impairment_none).
!   kmax, conductance_b, conductance_c
!                 : the (possibly impaired) per-point curve, with the same
!                   scaling as kmax_ref.
!   kcrit         : critical xylem conductance, with the same scaling as
!                   kmax_ref.
! *********************************************************************

! ---------------------------------------------------------------------
! Function to manage leaf conductance calculations for different
! impairment models.
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
                            ! land point (m/s)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_leaf_extreme(land_pnts)
                            ! Historic minimum leaf water potential for each
                            ! land point (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_conductance(land_pnts)
                            ! Leaf conductance for each land point (m/s)

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

CASE ( xylem_impairment_whole_trunk )
  CALL leaf_conductance_impaired_whole_trunk_jls( pft,                        &
                                                 land_pnts,                  &
                                                 water_potential,            &
                                                 kmax,                       &
                                                 kcrit,                      &
                                                 conductance_b,              &
                                                 conductance_c,              &
                                                 psi_leaf_extreme,           &
                                               ! INTENT OUT
                                                 leaf_conductance            &
                                                 )

CASE ( xylem_impairment_memory )
  CALL leaf_conductance_impaired_memory_jls( pft,                             &
                                             land_pnts,                       &
                                             water_potential,                 &
                                             kmax_ref,                        &
                                             kmax,                            &
                                             kcrit,                           &
                                             conductance_b,                   &
                                             conductance_c,                   &
                                           ! INTENT OUT
                                             leaf_conductance                 &
                                             )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
    'pft_xylem_impairment_model should be none (0), whole_trunk (2) or memory (3)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_jls

! ---------------------------------------------------------------------
! Function to manage xylem conductance calculations for different
! impairment models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It asumes that the input and output arrays are only large
!       enough to contain values for points with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE xylem_conductance_impaired_stom_opt_jls( pft,                       &
                                           n_water_potentials,                 &
                                           open_pnts,                          &
                                           open_index,                         &
                                           veg_index,                          &
                                           land_pnts,                          &
                                           water_potential,                    &
                                           kmax_ref,                           &
                                           kmax,                               &
                                           kcrit,                              &
                                           conductance_b,                      &
                                           conductance_c,                      &
                                           psi_extreme,                        &
                                         ! INTENT OUT
                                           xylem_conductance                   &
  )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls
USE xylem_impairment_whole_trunk_mod, ONLY: xylem_conductance_impaired_whole_trunk_stom_opt_jls
USE xylem_impairment_memory_mod,                                               &
    ONLY: xylem_conductance_impaired_memory_stom_opt_jls

USE um_types, ONLY: real_jlslsm

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, n_water_potentials                                                           &
                            ! Number of water potentials per open point
, open_pnts                                                                    &
                            ! Number of open stomata
, open_index(open_pnts)                                                        &
                            ! Index of open stomata into veg_index
, land_pnts                                                                    &
                            ! Number of land points
, veg_index(land_pnts)
                            ! Index of vegetated points on the land grid

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(n_water_potentials, open_pnts)                               &
                            ! Water potentials for each open point (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance for each
                            ! land point (m/s)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b for each land point (Pa)
, conductance_c(land_pnts)                                                     &
                            ! Conductance parameter c for each land point
, psi_extreme(land_pnts)
                            ! Historic minimum water potential at a given point
                            ! in the tree for each land point (Pa).
                            ! To calculate the leaf (root) conductance this
                            ! value is the historic minm=imum leaf (root) water
                            ! potential.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  xylem_conductance(n_water_potentials, open_pnts)
                            ! Xylem conductance for each open point (m/s)

! Local variables
INTEGER :: errcode, i, l

REAL(KIND=real_jlslsm) ::                                                      &
  kmax_open(open_pnts), kcrit_open(open_pnts), b_open(open_pnts),             &
  c_open(open_pnts)
                            ! Unimpaired curve parameters gathered onto the
                            ! open-point index.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_CONDUCTANCE_IMPAIRED_STOM_OPT_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and use the default
! conductance model when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  DO i = 1, open_pnts
    l = veg_index(open_index(i))
    kmax_open(i)  = kmax_ref(l)
    kcrit_open(i) = kcrit(l)
  END DO
  b_open(:) = conductance_b_pft(pft)
  c_open(:) = conductance_c_pft(pft)

  CALL xylem_conductance_jls( pft,                                             &
                              n_water_potentials,                              &
                              open_pnts,                                       &
                              water_potential,                                 &
                              kmax_open,                                       &
                              kcrit_open,                                      &
                              b_open,                                          &
                              c_open,                                          &
                            ! INTENT OUT
                              xylem_conductance                                &
                              )

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE ( pft_xylem_impairment_model(pft) )

CASE ( xylem_impairment_whole_trunk )
  CALL xylem_conductance_impaired_whole_trunk_stom_opt_jls( pft,               &
                                               n_water_potentials,            &
                                               open_pnts,                     &
                                               open_index,                    &
                                               veg_index,                     &
                                               land_pnts,                     &
                                               water_potential,               &
                                               kmax,                          &
                                               kcrit,                         &
                                               conductance_b,                 &
                                               conductance_c,                 &
                                               psi_extreme,                   &
                                             ! INTENT OUT
                                               xylem_conductance              &
                                               )

CASE ( xylem_impairment_memory )
  CALL xylem_conductance_impaired_memory_stom_opt_jls( pft,                   &
                                               n_water_potentials,            &
                                               open_pnts,                     &
                                               open_index,                    &
                                               veg_index,                     &
                                               land_pnts,                     &
                                               water_potential,               &
                                               kmax_ref,                      &
                                               kmax,                          &
                                               kcrit,                         &
                                               conductance_b,                 &
                                               conductance_c,                 &
                                             ! INTENT OUT
                                               xylem_conductance              &
                                               )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
   'pft_xylem_impairment_model should be none (0), whole_trunk (2) or memory (3)')

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE xylem_conductance_impaired_stom_opt_jls

! ---------------------------------------------------------------------
! Subroutine to manage leaf water potential calculations for different
! conductance models.
! NOTE: This function is designed for the stomatal optimisation model.
!       It returns an output array that only contain values for points
!       with open stomata.
! ---------------------------------------------------------------------
SUBROUTINE leaf_psi_impaired_jls( pft,                                         &
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
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model,                                 &
                   conductance_c_pft, conductance_b_pft

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_impairment_whole_trunk_mod, ONLY: leaf_psi_impaired_whole_trunk
USE xylem_impairment_memory_mod, ONLY: leaf_psi_impaired_memory
USE jules_vegetation_mod, ONLY: l_som_rhizo_series
USE xylem_hydraulics_CW_jls_mod, ONLY: ksr_path, som_psi_in_min

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
                            ! land point (m/s)
, kmax(land_pnts)                                                              &
                            ! Maximum xylem conductance for each land point
                            ! (m/s)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance for each land point
                            ! (m/s)
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
                            ! Leaf water potential for each open point (m/s)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point (m/s)

! Local variables
INTEGER :: errcode, i, j, l

REAL(KIND=real_jlslsm) :: b_ref(land_pnts), c_ref(land_pnts)
                            ! Unimpaired PFT curve parameters.
REAL(KIND=real_jlslsm) ::                                                      &
  psi_in(land_pnts), k_inlet(land_pnts)                                        &
                            ! l_som_rhizo_series: root inlet potential (Pa)
                            ! and impaired conductance there, one sample.
, e1(1, open_pnts), psi1(1, open_pnts), k1(1, open_pnts)                       &
                            ! One sample per open point.
, k_s
                            ! Soil-to-root conductance of the path.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_PSI_IMPAIRED_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not use xylem impairment model during spinup, and just call the base
! conductance models when no impairment model is selected.
IF (is_spinup .OR.                                                             &
    pft_xylem_impairment_model(pft) == xylem_impairment_none) THEN
  b_ref(:) = conductance_b_pft(pft)
  c_ref(:) = conductance_c_pft(pft)

  CALL leaf_psi_jls( pft,                                                      &
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

  CASE (xylem_impairment_whole_trunk, xylem_impairment_memory)
    IF ( .NOT. l_som_rhizo_series ) THEN
      CALL impaired_path( n_e_leaf, e_leaf, root_zone_psi, leaf_psi, leaf_k )
    ELSE
      !-----------------------------------------------------------------------
      ! Soil-to-root conductance in series ahead of the impaired plant, as
      ! leaf_psi_jls: the root inlet of each sample is psi_in = psi_src - E/K_s,
      ! the impaired plant path runs from there, and the whole-path marginal
      ! conductance is k_p / (1 + k(psi_in) / K_s), k(psi_in) this model's
      ! impaired conductance at the inlet. The soil is not damaged.
      !-----------------------------------------------------------------------
      DO i = 1, n_e_leaf
        psi_in(:) = root_zone_psi(:)
        DO j = 1, open_pnts
          l = veg_index(open_index(j))
          k_s = ksr_path(l)
          IF ( ksr_path(l) < 0.0 ) THEN
            ! No soil link for this PFT (fsmc_mod /= 2).
          ELSE IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
            psi_in(l) = MAX(root_zone_psi(l) - e_leaf(i,j) / k_s, som_psi_in_min)
          ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
            psi_in(l) = som_psi_in_min
          END IF
          e1(1,j) = e_leaf(i,j)
        END DO
        CALL impaired_path( 1, e1, psi_in, psi1, k1 )
        CALL leaf_conductance_impaired_jls( pft, land_pnts, psi_in, kmax_ref,  &
                                            kmax, kcrit, conductance_b,        &
                                            conductance_c, psi_leaf_extreme,   &
                                            k_inlet )
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
    END IF
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), whole_trunk (2) or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! The impaired plant path (root inlet to leaf), whole_trunk or memory model.
SUBROUTINE impaired_path( n_e, e_in, psi_root, psi_out, k_out )

INTEGER, INTENT(IN) :: n_e
REAL(KIND=real_jlslsm), INTENT(IN)  :: e_in(n_e, open_pnts),                  &
                                       psi_root(land_pnts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi_out(n_e, open_pnts),               &
                                       k_out(n_e, open_pnts)

SELECT CASE ( pft_xylem_impairment_model(pft) )
CASE (xylem_impairment_whole_trunk)
    CALL leaf_psi_impaired_whole_trunk( pft,                                   &
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
    CALL leaf_psi_impaired_memory( pft,                                        &
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

END SUBROUTINE leaf_psi_impaired_jls

! ---------------------------------------------------------------------
! Function to manage which xylem imparement model is applied.
! Note the differnet types of conductance models are no handled here.
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment ( n_land_pts                                &
,                                    n_open_pts                                &
,                                    open_index                                &
,                                    pft                                       &
,                                    psi_leaf                                  &
,                                    psi_root                                  &
,                                    psi_stem                                  &
,                                    leaf_conductance                          &
,                                    root_conductance                          &
,                                    impaired_k_max                            &
,                                    impaired_conductance_b                    &
,                                    impaired_conductance_c                    &
,                                    psi_leaf_extreme                          &
,                                    psi_root_extreme                          &
,                                    lai                                       &
,                                    canht                                     &
,                                    anetc                                     &
                                   )

USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: xylem_impairment_none,                         &
                                xylem_impairment_whole_trunk,                  &
                                xylem_impairment_memory
USE pftparm, ONLY: pft_xylem_impairment_model, kmax_pft

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
                            ! Number of open plant functional types
, open_index(n_open_pts)                                                       &
                            ! Index of open plant functional types
, pft
                            ! Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  psi_leaf(n_land_pts)                                                         &
                            ! Leaf water potential (Pa)
, psi_root(n_land_pts)                                                         &
                            ! Root water potential (Pa)
, psi_stem(n_land_pts)                                                         &
                            ! Stem water potential (Pa; memory model driver 4)
, leaf_conductance(n_land_pts)                                                 &
                            ! Leaf conductance (m/s)
, root_conductance(n_land_pts)
                            ! Root conductance (m/s)

REAL(KIND=real_jlslsm), INTENT(INOUT) ::                                       &
  impaired_k_max(n_land_pts)                                                   &
                            ! Maximum xylem conductance (m/s)
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

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  impaired_conductance_b(n_land_pts)                                           &
                            ! Conductance model parameter b (Pa). Only
                            ! changed by the kmax impairment model.
, impaired_conductance_c(n_land_pts)
                            ! Conductance model parameter c (unitless)



INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Do not update the xylem impairment model during spinup
IF(is_spinup) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

SELECT CASE (pft_xylem_impairment_model(pft))
  CASE(xylem_impairment_none)
    ! No xylem impairment

  CASE(xylem_impairment_whole_trunk)
    CALL update_xylem_impairment_whole_trunk( n_land_pts,                      &
                                              n_open_pts,                      &
                                              open_index,                      &
                                              pft,                             &
                                              psi_leaf,                        &
                                              psi_root,                        &
                                              psi_leaf_extreme,                &
                                              psi_root_extreme                 &
                                            )

  CASE(xylem_impairment_memory)
    CALL update_xylem_impairment_memory( n_land_pts,                           &
                                         n_open_pts,                           &
                                         open_index,                           &
                                         pft,                                  &
                                         psi_leaf,                             &
                                         psi_root,                             &
                                         psi_stem,                             &
                                         lai,                                  &
                                         canht,                                &
                                         anetc,                                &
                                         impaired_k_max                        &
                                       )

  CASE DEFAULT
    errcode = 1
    CALL ereport(RoutineName, errcode,                                         &
     'pft_xylem_impairment_model should be none (0), whole_trunk (2) or memory (3)')
END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment

! ---------------------------------------------------------------------
! Subroutine to calculate canopy water potential from canoopy
! transpiration and root zone water potential for different conductance
! models.
! NOTE: This function is designed for use outside the stomatal
!       optimisation model. It returns an output array that contains
!       values for all land points.
! ---------------------------------------------------------------------

SUBROUTINE canopy_impaired_psi_jls( land_pts                                   &
,                                   n_surface_pts                              &
,                                   surface_index                              &
,                                   e_canopy                                   &
,                                   root_zone_psi                              &
,                                   conductance_b                              &
,                                   conductance_c                              &
,                                   kmax_pft                                   &
,                                   psi_leaf_extreme                           &
,                                   psi_root_extreme                           &
                                  ! INTENT OUT
,                                   canopy_psi                                 &
,                                   canopy_k                                   &
                                  ! TMP TESTING PURPOSES ONLY
,                                   psi_leaf                                   &
,                                   e_leaf                                     &
,                                   gs_leaf                                     &
                                  )

USE jules_surface_types_mod,  ONLY: npft
USE ancil_info,               ONLY: nsurft
USE pftparm, ONLY: kmax_pft_ref => kmax_pft, kcrit
USE xylem_hydraulics_CW_jls_mod, ONLY: ksr_path

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
                            ! Number of land points
, n_surface_pts(nsurft)                                                      &
                            ! Number of surface points for each surface type
, surface_index(land_pts, nsurft)
                            ! Index for each surface type contained by each land
                            ! point. Starts with vegitation types.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_canopy(land_pts, nsurft)                                                     &
                            ! Transpiration rates for each open point
                            ! (kg/m2/s)
, root_zone_psi(land_pts, npft)                                                &
                             ! Water potential in the root zone (Pa)
, conductance_b(land_pts, npft)                                                &
                            ! Conductance parameter b for each plant
                            ! functional type (Pa)
, conductance_c(land_pts, npft)                                                &
                            ! Conductance parameter c for each plant
                            ! functional type (Pa)
, kmax_pft(land_pts, npft)                                                     &
                            ! Maximum xylem conductance for each plant
                            ! functional type (m/s)
, psi_leaf_extreme(land_pts, npft)                                             &
                            ! Minimum leaf water potential for each land point
                            ! (Pa)
, psi_root_extreme(land_pts, npft)
                            ! Minimum root water potential for each land point
                            ! (Pa)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  canopy_psi(land_pts, npft)                                                 &
                            ! Canopy water potential for each open point
                            ! (Pa)
, canopy_k(land_pts, npft)
                            ! Canopy conductance for each open point
                            ! (mol m-2 s-1)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  psi_leaf(land_pts, npft)                                                     &
                            ! Leaf water potential for each open point (Pa)
, e_leaf(land_pts, npft)                                                      &
                            ! Transpiration rates for each open point (kg/m2/s)
, gs_leaf(land_pts, npft)
                            ! Stomatal conductance for each open point (m s-1)

! Local variables
INTEGER :: i, n, l, pft

INTEGER ::                                                                     &
  open_pnts                                                                    &
                            ! Number of open stomata
, open_index(land_pts)
                            ! Index of open stomata for each land point

REAL(KIND=real_jlslsm) ::                                                      &
  e_canopy_open(land_pts)                                                      &
                            ! Transpiration rates for each open point
                            ! (kg/m2/s)
, canopy_psi_open(land_pts)                                                    &
                            ! Canopy water potential for each open point
                            ! (m/s)
, canopy_k_open(land_pts)                                                      &
                            ! Canopy conductance for each open point (m/s)
, kmax_ref(land_pts)                                                           &
                            ! Unimpaired maximum xylem conductance (m/s)
, kcrit_pts(land_pts)
                            ! Critical xylem conductance (m/s)

! Local variables
INTEGER :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CANOPY_PSI_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! No soil-to-root link here: this runs after the timestep, when ksr/ksr_path
! hold the last PFT's and leaf path's values, so the diagnostic stays plant
! only, as before (sf_stom sets ksr_path again before every solver call).
IF ( ALLOCATED(ksr_path) ) ksr_path(:) = -1.0

DO pft = 1, npft

  canopy_psi(:, pft) = root_zone_psi(:, pft)

  ! Identify points with open stomata
  open_pnts = 0
  DO i = 1, n_surface_pts(pft)
    l = surface_index(i, pft)

    ! If there is transpiration throught the vegitation tile note this point as having open stomata
    IF (e_canopy(l, pft) > TINY(0.0) .AND. e_canopy(l, pft) < 1e-2) THEN
      open_pnts = open_pnts + 1
      open_index(open_pnts) = i
    END IF
  END DO

  ! If there are no points with open stomata skip to the next pft
  IF (open_pnts == 0) CYCLE


  DO n = 1, open_pnts
    l = surface_index(open_index(n), pft)

    ! Initialise outputs
    canopy_psi_open(n) = root_zone_psi(l, pft)
    canopy_k_open(n) = kmax_pft(l, pft)

    ! Construct array of transpiration rates for points with open stomata
    e_canopy_open(n) = e_canopy(l, pft)
  END DO

  ! Leaf-basis (unscaled) curve, as this diagnostic has always used.
  kmax_ref(:)  = kmax_pft_ref(pft)
  kcrit_pts(:) = kcrit(pft)

    ! Calculate canopy water potential for points with open stomata
  CALL leaf_psi_impaired_jls( pft,                                             &
                              1,                                               &
                              land_pts,                                        &
                              open_pnts,                                       &
                              surface_index(:, pft),                           &
                              open_index,                                      &
                              e_canopy_open,                                   &
                              root_zone_psi(:,pft),                            &
                              kmax_ref,                                        &
                              kmax_pft(:,pft),                                 &
                              kcrit_pts,                                       &
                              conductance_b(:,pft),                            &
                              conductance_c(:,pft),                            &
                              psi_leaf_extreme(:,pft),                         &
                              psi_root_extreme(:,pft),                         &
                           ! INTENT OUT
                              canopy_psi_open,                                 &
                              canopy_k_open                                    &
                              )

  ! Insert canopy water potential and conductance values for points with open
  ! stomata into output arrays
  DO n = 1, open_pnts
    l = surface_index(open_index(n), pft)
    canopy_psi(l, pft) = canopy_psi_open(n)
    canopy_k(l, pft) = canopy_k_open(n)
  END DO

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE canopy_impaired_psi_jls


END MODULE xylem_impairment_mod
