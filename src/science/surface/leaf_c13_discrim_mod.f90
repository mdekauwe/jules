! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE leaf_c13_discrim_mod

! Diagnostic leaf 13C discrimination (permil) from the ci of whichever
! stomatal model is in use (discrimination_model). Ported from the
! equations of the vn5.6_jules_Cisotopes branch (A. Lavergne; used by
! Lavergne et al. 2022, GCB, and Palmer et al. 2022, JGR-B), without its
! Prentice (2014) ci, so that the discrimination follows the model's own
! stomatal response to soil water. Diagnostic only: no flux depends on it.
!
! The canopy value is the assimilation-weighted mean over leaf classes
! (sun/shade with can_rad_mod = 7), since the discrimination recorded in
! the assimilates is weighted by the carbon fixed, not by conductance.
! Time means should also be weighted by GPP (offline, with gpp_pft).

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: d13c_disc_pft, d13c_disc_alloc, d13c_disc_leaf

! Canopy leaf 13C discrimination of each PFT at each land point (permil);
! 0 where there is no assimilation. Not dumped.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: d13c_disc_pft(:,:)

! Fractionations (permil): diffusion through stomata, Rubisco (simple and
! full forms), photorespiration.
REAL(KIND=real_jlslsm), PARAMETER :: frac_a        = 4.4
REAL(KIND=real_jlslsm), PARAMETER :: frac_b_simple = 27.0
REAL(KIND=real_jlslsm), PARAMETER :: frac_b        = 28.0
REAL(KIND=real_jlslsm), PARAMETER :: frac_f        = 12.0

CONTAINS

! ---------------------------------------------------------------------
! Allocate the diagnostic (n_land_pts, npft), zeroed, if not done.
! ---------------------------------------------------------------------
SUBROUTINE d13c_disc_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

IF (.NOT. ALLOCATED(d13c_disc_pft)) THEN
  ALLOCATE(d13c_disc_pft(n_land_pts, npft))
  d13c_disc_pft(:,:) = 0.0
END IF

END SUBROUTINE d13c_disc_alloc

! ---------------------------------------------------------------------
! Leaf discrimination (permil) from ci, ca and the CO2 compensation
! point without dark respiration, ccp (all Pa). C4 leaves take a.
! ---------------------------------------------------------------------
FUNCTION d13c_disc_leaf( c3, ci, ca, ccp ) RESULT( disc )

USE jules_vegetation_mod, ONLY: discrimination_model,                         &
    discrimination_simple, discrimination_photoresp

LOGICAL, INTENT(IN) :: c3
REAL(KIND=real_jlslsm), INTENT(IN) :: ci, ca, ccp
REAL(KIND=real_jlslsm) :: disc

disc = 0.0
IF ( ca <= TINY(ca) ) RETURN

IF ( .NOT. c3 ) THEN
  disc = frac_a
  RETURN
END IF

SELECT CASE ( discrimination_model )
CASE ( discrimination_simple )
  disc = frac_a + ( frac_b_simple - frac_a ) * ci / ca
CASE ( discrimination_photoresp )
  disc = frac_a + ( frac_b - frac_a ) * ci / ca - frac_f * ccp / ca
END SELECT

END FUNCTION d13c_disc_leaf

END MODULE leaf_c13_discrim_mod
