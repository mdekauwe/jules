! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE smc_ext_mod

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SMC_EXT_MOD'

CONTAINS

SUBROUTINE smc_ext (npnts,nshyd,surft_pts,surft_index,ft                       &
,                   f_root,k_sat,sthu,v_open,v_sat,v_close                     &
,                   bexp,sathh                                                 &
,                   psi,soil_k,soil_to_root_k,wt_ext,fsmc,psi_root_zone)
!---------------------------------------------------------------------

! Description:
!     Calculates the soil moisture availability factor and
!     the fraction of the transpiration which is extracted from each
!     soil layer.

! Documentation : UM Documentation Paper 25
!---------------------------------------------------------------------
USE pftparm, ONLY: calc_rz_psi, fsmc_mod, root_psi_crit
USE hyd_psi_mod, ONLY: psi_from_sthu, bound_soil_psi
USE jules_vegetation_mod, ONLY: fsmc_shape, leaf_flux_mod, leaf_flux_stom_opt, &
                                stomata_model, stomata_sox, stomata_desica
USE hyd_con_ic_mod, ONLY: hyd_con_ic
USE jules_soil_mod, ONLY: l_bound_soil_wp, ds_psi, dzsoil

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) ::                                                         &
 npnts                                                                         &
                      ! Number of gridpoints.
,nshyd                                                                         &
                      ! Number of soil moisture layers.
,surft_pts                                                                     &
                      ! Number of points containing the
!                     !    given surface type.
,surft_index(npnts)                                                            &
                      ! Indices on the land grid of the
!                     !    points containing the given
!                     !    surface type
,ft                   ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 f_root(npnts,nshyd)                                                           &
                      ! Fraction of roots in each soil
!                     !    layer.
,k_sat(npnts,0:nshyd)                                                          &
!                     ! The saturated hydraulic conductivity (kg/m2/s)
,sthu(npnts,nshyd)                                                             &
                      ! Unfrozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation.
,v_sat(npnts,nshyd)                                                            &
                      ! Volumetric soil moisture
!                     !    concentration at saturation
!                     !    (m3 H2O/m3 soil) from JULES_SOIL_PROPS namelist
,v_open(npnts,nshyd)                                                           &
!                     ! Volumetric soil moisture
!                     !    concentration above which stomatal aperture
!                     !    is not limited by soil water (m3 H2O/m3 soil)
!                     ! When l_use_pft_psi=F, v_open
!                     !    is set from smvccl_soilt - fsmc_p0 *
!                     !    (smvccl_soilt - smvcwt_soilt) (in physiol).
!                     ! When l_use_pft_psi=T, v_open
!                     !    is calculated from psi_open.
,v_close(npnts,nshyd)                                                          &
                      ! Volumetric soil moisture
!                     !    concentration below which
!                     !    stomata close (m3 H2O/m3 soil).
!                     ! When l_use_pft_psi=F, v_close
!                     !    is set from smvcwt_soilt (in physiol).
!                     ! When l_use_pft_psi=T, v_close
!                     !    is calculated from psi_close.
,bexp(npnts,nshyd)                                                             &
                      ! Exponent in soil hydraulic characteristics.
,sathh(npnts,nshyd)                                                            
                      ! If l_vg=False, absolute value of the soil matric
                      ! suction at saturation in m
                      ! If l_vg=True, sathh = 1 / alpha, where alpha
                      ! (in m-1) is a parameter in the van Genuchten model

! psi/soil_k/soil_to_root_k are IN OUT, not OUT: physiol calls smc_ext once
! per PFT with the same (soil-tile) arrays, and they are output as
! diagnostics. Only this PFT's surft points are reset/written below, so a
! later call for a PFT with no points no longer wipes the values set by an
! earlier one (previously every call zeroed all npnts, so the diagnostics
! always came out as the last PFT's - usually all zeros).
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 psi(npnts,nshyd)                                                              &
!                     ! Negative soil water potential in each soil layer
!                     ! (Pa)
,soil_k(npnts,nshyd)                                                           &
!                     ! Soil layer conductivity (kg m-2 s-1)
,soil_to_root_k(npnts,nshyd)
!                     ! Soil to root conductance (kg m-2 s-1)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 wt_ext(npnts,nshyd)
!                     ! Cumulative fraction of transpiration
!                     !    extracted from each soil layer
!                     !    (kg/m2/s).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 psi_root_zone(npnts)                                                          &
!                     ! (negative) soil water potential (in Pa) in root zone
,fsmc(npnts)
!                     ! Soil moisture availability factor.

! work
INTEGER ::                                                                     &
 i,j,n                ! Loop counters

REAL(KIND=real_jlslsm) ::                                                      &
 fsmc_l(npnts,nshyd)                                                           &
!                     ! Soil moisture availability
!                     !    factor for each soil layer.
,psi_open(npnts,nshyd)                                                         &
!                     ! (negative) soil water potential (in Pa)
!                     !    of layers above which stomatal
!                     !    aperture is not limited by soil
!                     !    water.
,psi_close(npnts,nshyd)                                                        &
!                     ! (negative) soil water potential (in Pa)
!                     !    of layers above which stomatal
!                     !    aperture is not limited by soil
!                     !    water.
,soil_dk_dthk(npnts,nshyd)                                                     &
!                     ! The rate of change of soil_k with sthu
!                     !    (kg/m2/s) for each soil layer.
,sthu_open(npnts,nshyd)                                                        &
!                     ! Unfrozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation, above which stomatal
!                     !    aperture is not limited by soil
!                     !    water.
,sthu_close(npnts,nshyd)                                                       &
!                     ! Unfrozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation, below which stomata
!                     !    are closed.
,v_layer(npnts,nshyd)
!                     ! Volumetric soil moisture
!                     !    concentration of layer
!                     !    (m3 H2O/m3 soil).


REAL(KIND=real_jlslsm) ::                                                      &
v_root_zone(npnts)                                                             &
!                     ! Volumetric soil moisture
!                     !    concentration of root zone
!                     !    (m3 H2O/m3 soil)
,v_close_root_zone(npnts)                                                      &
!                     ! Volumetric soil moisture
!                     !    concentration of root zone below which
!                     !    stomata close (m3 H2O/m3 soil)
,v_open_root_zone(npnts)                                                       &
!                     ! Volumetric soil moisture
!                     !    concentration above which stomatal aperture is
!                     !    not limited by soil water.
!                     !    (m3 H2O/m3 soil).
,soil_weights(nshyd)
!                     ! Weightings for soil layer height.

REAL(KIND=real_jlslsm), PARAMETER :: sthu_min = 0.01 ! required by psi_from_sthu

REAL(KIND=real_jlslsm) :: ones(npnts)
REAL(KIND=real_jlslsm) :: soil_to_root_k_ft(npnts,nshyd)
                      ! soil_to_root_conductance result for this PFT, copied
                      ! into soil_to_root_k at this PFT's surft points only
                      ! (see the note on the psi/soil_k/soil_to_root_k
                      ! declaration).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SMC_EXT'

!----------------------------------------------------------------------
! Initialisations
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(n,i,j)                                                           &
!$OMP SHARED(nshyd,npnts,psi,v_layer,sthu,v_sat,wt_ext,psi_root_zone,fsmc,ones,soil_k, soil_to_root_k, &
!$OMP        surft_pts,surft_index)
DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
  DO i = 1,npnts
    v_layer(i,n) = sthu(i,n) * v_sat(i,n)
    wt_ext(i,n)  = 0.0
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
  DO j = 1,surft_pts
    i = surft_index(j)
    psi(i,n)            = 0.0
    soil_k(i,n)         = 0.0
    soil_to_root_k(i,n) = 0.0
  END DO
!$OMP END DO NOWAIT
END DO

!$OMP DO SCHEDULE(STATIC)
DO i = 1,npnts
  psi_root_zone(i)  = 0.0
  fsmc(i)           = 0.0
  ones(i)           = 1.0
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF ( fsmc_mod(ft) == 1 ) THEN

  v_root_zone       = calc_weighted_mean(npnts, nshyd, surft_pts,              &
                                         surft_index, v_layer, f_root)
  v_close_root_zone = calc_weighted_mean(npnts, nshyd, surft_pts,              &
                                         surft_index, v_close, f_root)
  v_open_root_zone  = calc_weighted_mean(npnts, nshyd, surft_pts,              &
                                         surft_index, v_open, f_root)

  ! SOX (stomata_model = stomata_sox) also needs psi_root_zone.
  IF ( fsmc_shape == 1 .OR. leaf_flux_mod == leaf_flux_stom_opt                &
       .OR. calc_rz_psi(ft) .OR. stomata_model == stomata_sox                  &
       .OR. stomata_model == stomata_desica ) THEN
!$OMP PARALLEL DO                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i)                                                             &
!$OMP SHARED(surft_pts,surft_index,psi_root_zone,v_root_zone,v_sat,sathh,bexp)
    DO j = 1,surft_pts
      i = surft_index(j)
      ! Have already checked for soil_props_const_z in init_pftparm in
      ! standalone jules, so just using soil props of top layer. UM doesn't
      ! allow this option (fsmc_shape=1) yet anyway
      psi_root_zone(i) = psi_from_sthu(v_root_zone(i) / v_sat(i,1),            &
                                       sathh(i,1), bexp(i,1), sthu_min)
    END DO
!$OMP END PARALLEL DO
  END IF

  IF ( stomata_model == stomata_sox ) THEN
    ! fsmc_shape is not used by SOX - fsmc is always 1.0 (psi_root_zone is
    ! set above).
    fsmc(:) = 1.0

  ELSE
    fsmc = fsmc_layer(npnts, surft_pts, surft_index, ft, v_root_zone,          &
                      v_close_root_zone, v_open_root_zone, ones, psi_root_zone)
  END IF !stomata_model

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,wt_ext,f_root,v_layer,v_close,nshyd)
  DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      ! Calculate the fraction of transpiration extracted from each soil layer
      ! based on the water in that layer above the point where stomata close.
      wt_ext(i,n) = MAX(f_root(i,n) * (v_layer(i,n) - v_close(i,n)), 0.0)
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  wt_ext = calc_norm_weights(npnts, nshyd, surft_pts, surft_index,             &
                             wt_ext)

ELSE IF (fsmc_mod(ft) == 2) THEN
  ! ----------------------------------------------------------------------
  ! Calculate the soil moisture availability factor for each layer from
  ! the soil to root resistivity.
  !
  ! Refrence G. B. Bonan et al 2014
  ! Note: This implimentation dose not use the root to stem conductance
  !       (equation A24) when calculating the below ground resitivity
  !       (equation A25) only the soil to root conductance (equation A23).
  !
  ! Jake Baguley 16/11/23
  ! ----------------------------------------------------------------------

  ! Calculat the soil moisture concentration limits
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,sthu_open,V_open,V_sat,sthu_close,V_close,nshyd)
  DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      sthu_open(i,n)  = V_open(i,n)  / V_sat(i,n)
      sthu_close(i,n) = V_close(i,n) / V_sat(i,n)
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  ! Calculate the soil water potential for each layer
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,psi,sthu,sathh,bexp,sthu_min,nshyd)
  DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      psi(i,n) = psi_from_sthu(sthu(i,n), sathh(i,n), bexp(i,n), sthu_min)
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  ! Apply boundry conditions to soil water potentials
  IF (l_bound_soil_wp) THEN
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,sathh,bexp,sthu_min,psi_open,sthu_open,     &
!$OMP        psi_close,sthu_close,nshyd)
    DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
      DO j = 1,surft_pts
      i = surft_index(j)
      psi_open(i,n)  = psi_from_sthu(sthu_open(i,n),  sathh(i,n),              &
                                     bexp(i,n), sthu_min)
      psi_close(i,n) = psi_from_sthu(sthu_close(i,n), sathh(i,n),              &
                                     bexp(i,n), sthu_min)
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL
    CALL bound_soil_psi(npnts,nshyd,surft_pts,surft_index,ft,                  &
                        psi_close, psi_open, psi,                              &
                        sthu, sthu_close)
  END IF

  ! NOTE: Inside hyd_con_ic it is designed to iterate over soil_points. By
  !       substituting surft_pts for soil_pts, and surft_index for soil_index
  !       this call makes it iterate over surface points.
  DO n = 1,nshyd
    CALL hyd_con_ic( npnts, surft_pts, surft_index,                            &
                     bexp(:,n), k_sat(:,n), sthu(:,n),                         &
                     soil_k(:,n), soil_dk_dthk(:,n))
  END DO

  soil_to_root_k_ft = soil_to_root_conductance(npnts,nshyd,surft_pts,       &
                                              surft_index,ft,f_root,soil_k)
  DO n = 1,nshyd
    DO j = 1,surft_pts
      i = surft_index(j)
      soil_to_root_k(i,n) = soil_to_root_k_ft(i,n)
    END DO
  END DO

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(nshyd,surft_pts,surft_index,wt_ext,soil_to_root_k,psi,            &
!$OMP        root_psi_crit)
  DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      ! Calculate the transpiration extracted from each soil layer asuming
      ! the root zone water potential is at the critical value
      ! (root_psi_crit).
      wt_ext(i,n) = MAX(soil_to_root_k(i,n) * (psi(i,n) - root_psi_crit(ft)),  &
                        1.0e-9)
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  ! If every layer is at or below root_psi_crit, every raw weight above sits
  ! at the 1.0e-9 floor, and normalising would give each layer an equal
  ! share regardless of thickness. That both makes psi_root_zone a plain
  ! per-layer mean (dominated by the thin, driest top layers) and extracts
  ! the same absolute amount from a 5 cm layer as a 21 cm one, drying the
  ! thin layers out first. The psi_root_zone == 0.0 fallback below can
  ! never fire in that case (the floor keeps the weights non-zero), so
  ! fall back to thickness weights here instead.
  DO j = 1,surft_pts
    i = surft_index(j)
    IF (MAXVAL(wt_ext(i,:)) <= 1.0e-9) THEN
      wt_ext(i,:) = dzsoil(1:nshyd)
    END IF
  END DO

  ! Calculate fraction of water uptake from each layer
  wt_ext = calc_norm_weights(npnts, nshyd, surft_pts, surft_index,             &
                             wt_ext)

  !----------------------------------------------------------------------
  ! Calculate the soil moisture availability factor for each layer and
  ! weight with the fractional water uptake from each layer to calculate
  ! the total availability factor.
  !----------------------------------------------------------------------
  DO n = 1,nshyd
    fsmc_l(:,n) = fsmc_layer(npnts,surft_pts,surft_index,                      &
                             ft, sthu(:,n), v_close(:,n),                      &
                             v_open(:,n), v_sat(:,n), psi(:,n))
  END DO

  fsmc = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,              &
                            fsmc_l, wt_ext)


  !---------------------------------------------------------------------
  ! If needed calculate the root zone water potential. This is the
  ! average of the water potential in all the layers weighted by the
  ! fraction of the transpiration taken from each layer.
  !---------------------------------------------------------------------

  IF ( leaf_flux_mod == leaf_flux_stom_opt .OR. calc_rz_psi(ft)               &
       .OR. stomata_model == stomata_desica ) THEN
    psi_root_zone = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,   &
                                       psi, wt_ext)

    ! If the root zone water potential is zero then realculate it weighted by
    ! the layer height. This occurs when fsmc, and subsiquently wt_ext, is zero
    !for all layers.
    soil_weights(:) = dzsoil(:) / SUM(dzsoil(:))
    DO i = 1,npnts
      IF (psi_root_zone(i) == 0.0) THEN
        psi_root_zone(i) = SUM(psi(i,:) * soil_weights(:))
      END IF
    END DO
  END IF

ELSE

  IF ( fsmc_shape == 1 .OR. leaf_flux_mod == leaf_flux_stom_opt                &
       .OR. calc_rz_psi(ft) .OR. stomata_model == stomata_desica ) THEN
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,psi,sthu,sathh,bexp,nshyd)
    DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
      DO j = 1,surft_pts
        i = surft_index(j)
        psi(i,n) = psi_from_sthu(sthu(i,n), sathh(i,n), bexp(i,n), sthu_min)
      END DO
!$OMP END DO NOWAIT
    END DO
!$OMP END PARALLEL
  END IF
  !----------------------------------------------------------------------
  ! Calculate the soil moisture availability factor for each layer and
  ! weight with the root fraction to calculate the total availability
  ! factor.
  !----------------------------------------------------------------------
  DO n = 1,nshyd
    fsmc_l(:,n) = fsmc_layer(npnts,surft_pts,surft_index,                      &
                             ft, sthu(:,n), v_close(:,n),                      &
                             v_open(:,n), v_sat(:,n), psi(:,n))
  END DO

  fsmc = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,              &
                            fsmc_l, f_root)

  !----------------------------------------------------------------------
  ! Calculate the fraction of the transpiration which is extracted from
  ! each soil layer.
  !----------------------------------------------------------------------

!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,fsmc,wt_ext,f_root,fsmc_l,nshyd)
  DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      IF (fsmc(i) > 0.0)                                                       &
        wt_ext(i,n) = f_root(i,n) * fsmc_l(i,n) / fsmc(i)
    END DO
!$OMP END DO NOWAIT
  END DO
!$OMP END PARALLEL

  !---------------------------------------------------------------------
  ! If needed calculate the root zone water potential. This is the
  ! average of the water potential in all the layers weighted by the
  ! fraction of the transpiration taken from each layer.
  !---------------------------------------------------------------------

  IF ( leaf_flux_mod == leaf_flux_stom_opt .OR. calc_rz_psi(ft)               &
       .OR. stomata_model == stomata_desica ) THEN
    psi_root_zone = calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,   &
                                       psi, wt_ext)

    ! If the root zone water potential is zero then realculate it weighted by
    ! the layer height. This occurs when fsmc, and subsiquently wt_ext, is zero
    !for all layers.
    soil_weights(:) = dzsoil(:) / SUM(dzsoil(:))
    DO i = 1,npnts
      IF (psi_root_zone(i) == 0.0) THEN
        psi_root_zone(i) = SUM(psi(i,:) * soil_weights(:))
      END IF
    END DO
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE smc_ext

FUNCTION calc_norm_weights(npnts, nshyd, surft_pts, surft_index,               &
                           weights) RESULT (norm_weights)
!-----------------------------------------------------------------------------
! Description:
!   Normalises an array of weights over the soil levels dimension.
!   The weights should be >= 0.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

! function arguments
INTEGER, INTENT(IN) :: npnts              ! Number of gridpoints.
INTEGER, INTENT(IN) :: nshyd              ! Number of soil moisture layers.
INTEGER, INTENT(IN) :: surft_pts          ! Number of points containing the
!                                           !    given surface type.
INTEGER, INTENT(IN) :: surft_index(npnts) ! Indices on the land grid of the
!                                           !    points containing the given
!                                           !    surface type.

REAL(KIND=real_jlslsm), INTENT(IN) :: weights(npnts,nshyd)
                                          ! Unnormalised weights.
!                                           ! Should all be >= 0.

! work
REAL(KIND=real_jlslsm) :: norm_factor     ! Factor used in the normalisation
INTEGER :: i,j,n                          ! Loop counters

! returns
REAL(KIND=real_jlslsm) :: norm_weights(npnts,nshyd)         ! Normalised weights

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
INTEGER :: errcode

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_NORM_WEIGHTS'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

norm_weights(:, :) = 0.0
norm_factor = 0.0

DO n = 1,nshyd
  DO j = 1,surft_pts
    i = surft_index(j)

    norm_factor = SUM(weights(i,:))

    IF (norm_factor > 0.0) THEN
      norm_weights(i,n) = MAX(weights(i,n) / norm_factor, 1.0e-9)
    ELSE
      WRITE(*,*) 'ERROR: Normalisation factor is zero in calc_norm_weights'
      WRITE(*,*) 'npnts = ', npnts, ' nshyd = ', nshyd
      WRITE(*,*) 'surft_pts = ', surft_pts
      WRITE(*,*) 'surft_index = ', surft_index
      WRITE(*,*) 'weights = ', weights(i,:)

      errcode = 1
      CALL ereport(RoutineName, errcode,                                         &
                   'Normalisation factor is zero in calc_norm_weights')

    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION calc_norm_weights

FUNCTION calc_weighted_mean(npnts, nshyd, surft_pts, surft_index,              &
                            var, norm_weights) RESULT (weighted_mean)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the weighted mean of a variable over soil levels
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

! function arguments
INTEGER, INTENT(IN) :: npnts              ! Number of gridpoints.
INTEGER, INTENT(IN) :: nshyd              ! Number of soil moisture layers.
INTEGER, INTENT(IN) :: surft_pts          ! Number of points containing the
!                                           !    given surface type.
INTEGER, INTENT(IN) :: surft_index(npnts) ! Indices on the land grid of the
!                                           !    points containing the given
!                                           !    surface type.

REAL(KIND=real_jlslsm), INTENT(IN) :: var(npnts,nshyd)
                                          ! Variable to take the weighted mean of
REAL(KIND=real_jlslsm), INTENT(IN) :: norm_weights(npnts,nshyd)
!                                           ! Normalised weights to use in the mean

! work
INTEGER :: i,j,n                          ! Loop counters

! returns
REAL(KIND=real_jlslsm) :: weighted_mean(npnts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_WEIGHTED_MEAN'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

weighted_mean(:)=0.0

DO n = 1,nshyd
  !CDIR NODEP
  DO j = 1,surft_pts
    i = surft_index(j)

    weighted_mean(i) = weighted_mean(i) + norm_weights(i,n) * var(i,n)

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION calc_weighted_mean

FUNCTION fsmc_layer(npnts,surft_pts,surft_index,                               &
                    pft, v, v_close,                                           &
                    v_open, extra_factor, psi) RESULT (fsmc_l)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the soil water availability factor of a soil layer
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE jules_vegetation_mod, ONLY: fsmc_shape
USE pftparm, ONLY: psi_open, psi_close, fsmc_q

IMPLICIT NONE

INTEGER, INTENT(IN) :: npnts
INTEGER, INTENT(IN) :: surft_pts
INTEGER, INTENT(IN) :: surft_index(npnts)
INTEGER, INTENT(IN) :: pft  ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) :: v(npnts), extra_factor(npnts)
                              ! v*extra_factor is the volumetric soil moisture
!                             !    concentration in layer
!                             !    (m3 H2O/m3 soil).
REAL(KIND=real_jlslsm), INTENT(IN) :: v_close(npnts)
                              ! Volumetric soil moisture
!                             !    concentration below which
!                             !    stomata close (m3 H2O/m3 soil).
REAL(KIND=real_jlslsm), INTENT(IN) :: v_open(npnts)
                              ! Volumetric soil moisture
!                             !    concentration above which
!                             !    stomatal aperture is not limited by soil water.
!                             !    (m3 H2O/m3 soil).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi(npnts)
!                             ! (negative) soil water potential (in Pa)

! internal
REAL(KIND=real_jlslsm) :: x, x_open, x_close

! returns
REAL(KIND=real_jlslsm) :: fsmc_l(npnts)
                              ! soil water availability factor for this soil
!                             !    layer.

INTEGER :: i,j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FSMC_LAYER'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,x_open,x_close,x)                                            &
!$OMP SHARED(surft_pts,surft_index,fsmc_shape,fsmc_l,extra_factor,psi_open,    &
!$OMP        psi_close,psi,v_open,v_close,v,pft)
DO j = 1,surft_pts
  i = surft_index(j)

  IF ( fsmc_shape == 0 ) THEN
    x_open  = v_open(i)
    x_close = v_close(i)
    x       = v(i)

    IF ( ABS(x_open - x_close) > 0.0 ) THEN
      fsmc_l(i) = (x * extra_factor(i) - x_close) / (x_open - x_close)
    ELSE
      fsmc_l(i) = 0.0
    END IF
  ELSE
    x_open  = psi_open(pft)
    x_close = psi_close(pft)
    x       = psi(i)

    IF ( ABS(x_open - x_close) > 0.0 ) THEN
      fsmc_l(i) = (x - x_close) / (x_open - x_close)
    ELSE
      fsmc_l(i) = 0.0
    END IF
  END IF

  fsmc_l(i) = MAX(fsmc_l(i),0.0)
  fsmc_l(i) = MIN(fsmc_l(i),1.0)
  ! Optional curvature of the stress factor (fsmc_q; 1 = linear).
  IF ( fsmc_q(pft) /= 1.0 ) fsmc_l(i) = fsmc_l(i)**fsmc_q(pft)
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION fsmc_layer

FUNCTION soil_to_root_conductance(npnts,nshyd,surft_pts,surft_index,ft         &
,                                 f_root,soil_k)                     &
                           RESULT(soil_to_root_k)
!---------------------------------------------------------------------

! Description:
!     Calculates the conductance to water between the soil and root
!     per unit area of soil (kg m-3 s-1).

! Refrence: G.B.Bonan et al, 2014, Geoscientific Model Development
!---------------------------------------------------------------------
! TODO: Could add OMP to loops
USE pftparm, ONLY: min_rootc_pft, root_radi_pft, rootc_density_pft, rmass

USE conversions_mod, ONLY: pi
USE jules_soil_mod, ONLY: dzsoil

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) ::                                                         &
 npnts                                                                         &
                      ! Number of gridpoints.
,nshyd                                                                         &
                      ! Number of soil moisture layers.
,surft_pts                                                                     &
                      ! Number of points containing the
!                     !    given surface type.
,surft_index(npnts)                                                            &
                      ! Indices on the land grid of the
!                     !    points containing the given
!                     !    surface type
,ft                   ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 f_root(npnts,nshyd)                                                           &
                      ! Fraction of roots in each soil
!                     !    layer.
,soil_k(npnts,nshyd)
                      ! Soil conductivity in each soil
!                     !    layer (kg m-2 s-1).

! Internal variables
INTEGER ::                                                                     &
 i,j,n                ! Loop counters

REAL(KIND=real_jlslsm) ::                                                      &
 root_l(npnts,nshyd)                                                           &
                      ! Root length per unit soil volume
!                     ! (m m-3).
,rootc_pft_local(npnts)                                                        &
                      ! Holds local copy of root mass per unit soil area
!                     ! (kg m-2).
,rootc_per_unit_length
                      ! Holds local calculaion of the root mass per unit root
!                     ! length (kg m-1).

! Return variable
REAL(KIND=real_jlslsm) ::                                                      &
 soil_to_root_k(npnts,nshyd)
                      ! Soil to root water conductance per unit area of soil
!                     ! (kg m-3 s-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='soil_to_root_conductance'

!----------------------------------------------------------------------
! Initialisations
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Init local arrays
rootc_pft_local(:) = min_rootc_pft(ft)
root_l(:,:) = 0.0
soil_to_root_k(:,:) = 0.0

! Init local copy of root mass with inforced minimum values
! Convert to root mass by dividing by root carbon to root mass ratio (rmass)
!DO j = 1,surft_pts
!  i = surft_index(j)
!  rootc_pft_local(i) = MAX(rootc_pft(i), min_rootc_pft(ft))
!  rootc_pft_local(i) = rootc_pft_local(i) / rmass(ft)
!END DO

! The following two loops calcualte equation A23 in Bonan et al 2014
! The first calculates the soil to root conductance per unit root length
! (kg m-2 s-1).
! Note r_{sj}/r{r} is substituted by:
!
!  r_{s,j}     f_{r,i} m_{r}
! -------- = (--------------) ^ (-1/2)
!  r_{r}       rho_r dz_{j}
!
! Where f_{r,j} is the root fraction in the layer, m_{r} is the total root mass (kg)
! and rho_{r} is the root density (kg m-3).
DO n = 1,nshyd
  DO j = 1,surft_pts
    i = surft_index(j)
    soil_to_root_k(i,n) = 4*pi*soil_k(i,n) /                                       &
                          LOG(rootc_density_pft(ft) * dzsoil(n) /                  &
                              (f_root(i,n) * rootc_pft_local(i)))
  END DO
END DO

! To convert root conductance in each layer from per unit root length to per
! unit soil area we need to multiply by root length per unit area of soil (m m-2).
rootc_per_unit_length = rootc_density_pft(ft) * pi * root_radi_pft(ft)**2
DO n = 1,nshyd
  DO j = 1,surft_pts
    i = surft_index(j)
    soil_to_root_k(i,n) = soil_to_root_k(i,n) *(f_root(i,n) * rootc_pft_local(i))/ &
                                                rootc_per_unit_length
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION soil_to_root_conductance

END MODULE smc_ext_mod
