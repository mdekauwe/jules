! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE leaf_c13_discrim_mod

! Leaf carbon isotope discrimination against 13C (permil), worked out from
! the leaf's internal CO2 (ci) given by whichever stomatal model is in use.
! It is a diagnostic: no flux depends on it.
!
! The equations come from the vn5.6_jules_Cisotopes branch (A. Lavergne;
! used by Lavergne et al. 2022, GCB, and Palmer et al. 2022, JGR-B). That
! branch also had its own ci from Prentice et al. (2014); this module
! leaves it out and uses the model's ci, so the discrimination follows the
! model's own stomatal response to soil water. The form is chosen by
! discrimination_model.
!
! In the two-leaf canopy (can_rad_mod = 7), the canopy value is the mean
! over the sunlit and shaded leaves weighted by each one's gross
! photosynthesis, since the discrimination stored in the sugars is set by
! the carbon each leaf fixes. Other canopies use the canopy ci. Means over
! time should also be weighted by GPP (afterwards, with the output
! gpp_pft).

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
! Allocate d13c_disc_pft (n_land_pts, npft), filled with 0, unless it is
! already allocated.
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
! Return how strongly one leaf discriminates against 13C in
! photosynthesis (permil), for the leaf d13C output. ci is the leaf's
! internal CO2, ca the CO2 in the air around the leaf and ccp the CO2
! compensation point without dark respiration (all partial pressures,
! Pa). For C3 leaves:
!
!   simple form:           D = a + (b - a) ci / ca
!   with photorespiration: D = a + (b - a) ci / ca - f ccp / ca
!
! where a = 4.4 is the fractionation in diffusion through the stomata,
! b = 27 (simple) or 28 is that by Rubisco, and f = 12 is that in
! photorespiration. C4 leaves take D = a. D is 0 when ca is 0.
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
