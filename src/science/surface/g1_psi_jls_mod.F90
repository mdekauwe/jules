! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE g1_psi_jls_mod

! *********************************************************************
! Medlyn et al. (2011) stomatal model with g1 reduced by the pre-dawn water
! potential (stomata_model = 6), following Zhou et al. (2013) and Eqn. 3 of
! De Kauwe et al. (2015):
!
!   g1 = g1_stomata exp(g1b_stomata psi_pd)
!
! with psi_pd (MPa) the pre-dawn water potential and g1_stomata the
! well-watered g1. As De Kauwe et al. (2015), psi_pd is taken as the root
! zone soil water potential, assuming that leaf and soil equilibrate
! overnight (no night-time transpiration).
!
! Pre-dawn: psi_soil_pd follows psi_root_zone at every dark timestep (no
! incident PAR) and is held through the day, so in daylight it is the soil
! water potential of the last dark step before sunrise. This needs no clock
! time, so it works for any timestep length and longitude. If there has been
! no dark step for a day (polar summer) it is updated anyway.
!
! There is no soil moisture factor (fsmc) on photosynthesis: water stress
! acts through g1 only. sf_stom passes the factor exp(g1b_stomata psi_pd)
! to leaf_limits in place of fsmc, as for DESICA.
!
! psi_soil_pd is dumped ('psi_soil_pd', MPa). A value >= 0 (e.g. an initial
! condition of 0) means "not set": it then starts at psi_root_zone.
!
! References:
! Medlyn et al. (2011) Glob. Change Biol. 17: 2134-2144.
! Zhou et al. (2013) Agric. For. Meteorol. 182-183: 204-214.
! De Kauwe et al. (2015) Biogeosciences 12: 7503-7518,
!   doi:10.5194/bg-12-7503-2015.
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: g1_psi_alloc, g1_psi_factor, psi_soil_pd

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  psi_soil_pd(:,:),                                                            &
                            ! Pre-dawn water potential (MPa), (land_pts, npft).
  time_since_dark(:,:)
                            ! Time since psi_soil_pd was last updated (s).

REAL(KIND=real_jlslsm), PARAMETER :: max_time_since_dark = 86400.0
                            ! Update psi_soil_pd at least this often (s).
REAL(KIND=real_jlslsm), PARAMETER :: pa_to_mpa = 1.0e-6

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='G1_PSI_JLS_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! Allocate the state (land_pts, npft), zero = not set, if not done yet.
!-----------------------------------------------------------------------------
SUBROUTINE g1_psi_alloc( land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: land_pts

IF ( .NOT. ALLOCATED(psi_soil_pd) ) THEN
  ALLOCATE( psi_soil_pd(land_pts,npft), time_since_dark(land_pts,npft) )
  psi_soil_pd(:,:)     = 0.0
  time_since_dark(:,:) = 0.0
END IF

END SUBROUTINE g1_psi_alloc

!-----------------------------------------------------------------------------
! Update the pre-dawn water potential and return the g1 factor
! exp(g1b_stomata psi_pd) (0-1) that multiplies g1_stomata in leaf_limits.
!-----------------------------------------------------------------------------
SUBROUTINE g1_psi_factor( ft, land_pts, veg_pts, veg_index, ipar,             &
                          psi_root_zone, fg1 )

USE pftparm, ONLY: g1b_stomata
USE timestep_mod, ONLY: timestep

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: ipar(land_pts)
                            ! Incident PAR (W m-2).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root_zone(land_pts)
                            ! Root zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(OUT) :: fg1(land_pts)
                            ! g1 / g1_stomata (0-1).

INTEGER :: l, m

CALL g1_psi_alloc( land_pts )

fg1(:) = 1.0
DO m = 1,veg_pts
  l = veg_index(m)
  ! Dark (night), not set yet, or no night for a day: take the soil water
  ! potential now. Otherwise (daylight) hold the last night value.
  IF ( ipar(l) <= 0.0 .OR. psi_soil_pd(l,ft) >= 0.0 .OR.                      &
       time_since_dark(l,ft) >= max_time_since_dark ) THEN
    psi_soil_pd(l,ft)     = MIN(psi_root_zone(l), 0.0) * pa_to_mpa
    time_since_dark(l,ft) = 0.0
  ELSE
    time_since_dark(l,ft) = time_since_dark(l,ft) + timestep
  END IF
  fg1(l) = EXP(g1b_stomata(ft) * psi_soil_pd(l,ft))
END DO

END SUBROUTINE g1_psi_factor

END MODULE g1_psi_jls_mod
