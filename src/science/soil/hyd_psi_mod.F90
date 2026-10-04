! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE hyd_psi_mod

USE water_constants_mod, ONLY: rho_water  !  Density of pure water (kg m-3).
USE planet_constants_mod, ONLY: g
  !  Mean acceleration due to gravity at earth's surface (m s-2)   .

USE jules_soil_mod, ONLY: l_vg_soil
! If l_vg_soil=False, uses Brooks and Corey relation.
! If l_vg_soil=True, uses van Genuchten relation.

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HYD_PSI_MOD'

CONTAINS

FUNCTION psi_from_sthu(sthu, sathh, b, sthu_min) RESULT (psi)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the (negative) soil water potential (in Pa) from the volumetric
!   soil water content as a fraction of saturation.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sthu,                                                                        &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.
  sathh,                                                                       &
    ! If l_vg_soil=False, absolute value of the soil matric suction at
    ! saturation (m).
    ! If l_vg_soil=True, sathh = 1 / alpha, where alpha (in m-1) is a
    ! parameter in the van Genuchten model.
  b,                                                                           &
    ! Exponent in soil hydraulic characteristics.
  sthu_min
    ! Minimum value of sthu for use in the calculation.

!-----------------------------------------------------------------------------
! Function result:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  psi  ! (Negative) soil water potential (Pa).

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  hh  ! Absolute soil water suction (m).

hh = hh_from_sthu(sthu, sathh, b, sthu_min)
psi = - hh * rho_water * g

END FUNCTION psi_from_sthu

!#############################################################################

FUNCTION sthu_from_psi(psi, sathh, b, sthu_min) RESULT (sthu)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the volumetric soil water content as
!   a fraction of saturation from the (negative) soil water potential (in Pa).
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  psi,                                                                         &
    ! (Negative) soil water potential (Pa).
  sathh,                                                                       &
    ! If l_vg_soil=False, absolute value of the soil matric suction at
    ! saturation (m).
    ! If l_vg_soil=True, sathh = 1 / alpha, where alpha (m-1) is a parameter
    ! in the van Genuchten model.
  b,                                                                           &
    ! Exponent in soil hydraulic characteristics.
  sthu_min
    ! Minimum value of sthu for use in the calculation.

!-----------------------------------------------------------------------------
! Function result:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 sthu
   ! Unfrozen soil moisture content of each layer as a fraction of
   ! saturation.

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  hh ! Absolute soil water suction (m).

hh = - psi / ( rho_water * g )
sthu = sthu_from_hh(hh, sathh, b, sthu_min)

END FUNCTION sthu_from_psi

!#############################################################################

FUNCTION hh_from_sthu(sthu, sathh, b, sthu_min) RESULT (hh)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the absolute soil water suction (in m) from the volumetric
!   soil water content as a fraction of saturation.
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  sthu,                                                                        &
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.
  sathh,                                                                       &
    ! If l_vg_soil=False, absolute value of the soil matric
    ! suction at saturation (m).
    ! If l_vg_soil=True, sathh = 1 / alpha, where alpha
    ! (m-1) is a parameter in the van Genuchten model.
  b,                                                                           &
    ! Exponent in soil hydraulic characteristics.
  sthu_min
    ! Minimum value of sthu for use in the calculation.

!-----------------------------------------------------------------------------
! Function result:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  hh  ! Absolute soil water suction (m).

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  bracket,                                                                     &
    ! Intermediate step in calculation when l_vg=True.
  sthu_local
    ! Local value of sthu (after making sure it is not below sthu_min).

IF (sthu < sthu_min) THEN
  sthu_local = sthu_min
ELSE
  sthu_local = sthu
END IF

IF ( l_vg_soil ) THEN
  bracket = -1.0 + sthu_local ** ( -b - 1.0 )
  hh = sathh * bracket ** ( b / ( b + 1.0 ) )
ELSE
  hh  =  sathh / (sthu_local ** b)
END IF

END FUNCTION hh_from_sthu

!#############################################################################

FUNCTION sthu_from_hh(hh, sathh, b, sthu_min) RESULT (sthu)
!-----------------------------------------------------------------------------
! Description:
!   Calculates the volumetric soil water content as
!   a fraction of saturation from the absolute soil water suction (in m).
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN):
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  hh,                                                                          &
    ! Absolute soil water suction (m).
  sathh,                                                                       &
    ! If l_vg_soil=False, absolute value of the soil matric
    ! suction at saturation (m).
    ! If l_vg_soil=True, sathh = 1 / alpha, where alpha
    ! (m-1) is a parameter in the van Genuchten model.
  b,                                                                           &
    ! Exponent in soil hydraulic characteristics.
  sthu_min
    ! Minimum value of sthu for use in the calculation .

!-----------------------------------------------------------------------------
! Function result:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  sthu
    ! Unfrozen soil moisture content of each layer as a fraction of
    ! saturation.

!-----------------------------------------------------------------------------
! Local variables:
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
  bracket
    ! Intermediate step in calculation when l_vg=True.

IF ( l_vg_soil ) THEN
  bracket = 1.0 + ( hh / sathh )** ( ( b + 1.0 ) / b )
  sthu    = ( 1.0 / bracket )** ( 1.0 / (1.0 + b ) )
ELSE
  sthu    = ( hh / sathh )** ( -1.0 / b )
END IF

IF (sthu < sthu_min) THEN
  sthu = sthu_min
END IF

END FUNCTION sthu_from_hh

SUBROUTINE bound_soil_psi(npnts,nshyd,surft_pts,surft_index,ft,                &
                          min_psi, max_psi, psi,                               &
                          sthu, sthu_at_min_psi)
!-----------------------------------------------------------------------------
! Description:
!   Applies boundry condtions to soil water potential.
!
!   When l_ds_correction = false
!     Limits the minimum water potential to min_psi.
!   When l_ds_correction = true
!     Applies the dry soil aproximaion from S.W. Webb 2000.
!
!   Dry soil aproximation:
!     The dry soil approximation follows that outlined in S.W.
!     Webb 2000. For the purpose of computation efficiency the
!     implementation asumes that the soil moisture content
!     of the matching point is known (sthu_at_min_psi, where the
!     retention curve reaches min_psi = psi_close). Below it,
!     log10(-psi) is a straight line in S_l between the two points
!     (0, log_10(-ds_psi)) and (sthu_at_min_psi, log_10(-min_psi)):
!
!      d log_10(P_cap)      log_10(-ds_psi) - log_10(-min_psi)
!     ----------------  = -------------------------------------
!          d S_l                    sthu_at_min_psi
!
!     i.e. psi = ds_psi * (min_psi/ds_psi)**(S_l/sthu_at_min_psi),
!     continuous at the matching point and falling to ds_psi (oven
!     dry) at S_l = 0.
!
! Code Description:
!   Language: Fortran 90.
!-----------------------------------------------------------------------------

USE jules_soil_mod, ONLY: l_ds_correction, ds_psi

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
 min_psi(npnts,nshyd)                                                          &
                      ! Minimum water potential boundry (Pa)
!                     ! When l_ds_correction = false
!                     !     Minimum possible water potential
!                     ! When l_ds_correction = true
!                     !     Water potential at the matching point
!                     !     above which the dry soil approximation
!                     !     is applied. Refrence S.W. Webb 2000.
,max_psi(npnts,nshyd)                                                          &
                      ! Maximum water potential (Pa)
,sthu(npnts,nshyd)                                                             &
                      ! Unfrozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation. Only used
!                     !    if l_ds_correction is True
,sthu_at_min_psi(npnts,nshyd)
                      ! Unfrozen soil moisture content of
!                     !    each layer as a fraction of
!                     !    saturation below which the dry soil
!                     !    approximation is used: the matching
!                     !    point, where psi = min_psi (psi_close).
!                     !    Only used when l_ds_correction is True

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 psi(npnts,nshyd)     ! Negative soil water potential in each soil layer
!                     ! (Pa)

! work
INTEGER ::                                                                     &
 i,j,n                ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BOUND_SOIL_PSI'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Apply maximum value limit to psi
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,psi,max_psi,nshyd)
DO n = 1,nshyd
!$OMP DO SCHEDULE(STATIC)
  DO j = 1,surft_pts
    i = surft_index(j)
    psi(i,n) = MIN(psi(i,n),max_psi(i,n))
  END DO
!$OMP END DO NOWAIT
END DO
!$OMP END PARALLEL

! Modify psi if bellow psi_close
!$OMP PARALLEL                                                                 &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(j,i,n)                                                           &
!$OMP SHARED(surft_pts,surft_index,psi,min_psi,sthu,sthu_at_min_psi,          &
!$OMP        l_ds_correction,nshyd)
DO n = 1,nshyd
  ! Apply dry soil corection to psi
  ! Refrence; Webb 2000
  IF (l_ds_correction(n)) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      if(psi(i,n) < min_psi(i,n)) THEN
        psi(i,n) = ds_psi * (min_psi(i,n) / ds_psi)**(sthu(i,n)/sthu_at_min_psi(i,n))
      end if
    END DO
!$OMP END DO NOWAIT
  ! Apply minimum cut to psi
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j = 1,surft_pts
      i = surft_index(j)
      psi(i,n) = MAX(psi(i,n),min_psi(i,n))
    END DO
!$OMP END DO NOWAIT
  END IF
END DO
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bound_soil_psi



END MODULE hyd_psi_mod
