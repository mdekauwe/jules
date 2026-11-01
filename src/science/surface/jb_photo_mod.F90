! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE jb_photo_mod

! *********************************************************************
! Electron transport after Johnson and Berry (2021), used with
! photo_model = photo_johnson.
!
! Photosynthesis is limited either by the enzyme Rubisco or by the supply
! of electrons from the light reactions. This module changes only how the
! electron supply responds to light. The Rubisco-limited rate, the CO2
! compensation point and all the temperature responses stay those of the
! Farquhar model. In the Farquhar model the electron flux J follows an
! empirical curve in light with a maximum Jmax. In the Johnson-Berry model
! it is set by the flow of electrons through the cytochrome b6f complex
! (their Eqns 15c, 30a and 30c):
!
!   JP700 = Vqmax I1 / (I1 + Vqmax)
!   J     = JP700 / eta
!   eta   = 1 - nL/nC + (3 + 7 G*/Cc) / ((4 + 8 G*/Cc) nC)
!
! where JP700 is the electron flux through photosystem I and J that
! through photosystem II (both mol e- m-2 s-1), Vqmax the maximum activity
! of cytochrome b6f (mol e- m-2 s-1), I1 the light used by photosystem I
! (as an electron flux, mol m-2 s-1), G* the CO2 compensation point and Cc
! the CO2 at the chloroplast (both Pa), and nL and nC the ATP made per
! electron in linear and cyclic electron flow. J then sets the
! light-limited rate as in the Farquhar model:
!
!   Aj = J/4 (Cc - G*) / (Cc + 2 G*)
!
! The parameters follow Lamour et al. (2026), who matched the
! Johnson-Berry model to the Farquhar parameters used in vegetation
! models:
!
!   1. I1 is set so that at low light J rises with light at the same rate
!      as in the Farquhar model, at Cc = jb_ci_ref (their Eqn 10):
!        I1 = eta_ref alpha_elec acr
!      where eta_ref is eta at jb_ci_ref, alpha_elec the quantum
!      efficiency of electron transport and acr the absorbed light.
!   2. Vqmax at 25 degC is set so the two models give the same J at the
!      light and CO2 of a measured A-Ci curve (their Eqn 11: incident
!      light 1800 umol m-2 s-1, leaf absorptance 0.85, Ci = 800 ppm), from
!      the PFT's Jmax at 25 degC. See jb_vqmax25.
!   3. Vqmax has the temperature response of Jmax, as in ELM-FATES.
!
! So the existing Vcmax and Jmax parameters can be used as they are. A
! fit of the Johnson-Berry model to A-Ci curves would give better Vqmax.
!
! In JULES the array jmax holds Vqmax, and je holds JP700 / eta_ref, which
! is J at Cc = jb_ci_ref. The code that splits je between sunlit and
! shaded leaves or canopy layers and passes it around is then unchanged.
! The effect of CO2 on J is added where the light-limited rate (wlite) is
! formed, as je times jb_eta_scale. Mesophyll conductance is taken as
! infinite, so Cc = ci, the leaf's internal CO2.
!
! References:
! Johnson & Berry, 2021, Photosynth. Res., 148: 101--136,
!   https://doi.org/10.1007/s11120-021-00840-4 (code: doi 10.5281/zenodo.4759246)
! Lamour et al., 2026, bioRxiv, https://doi.org/10.64898/2026.05.08.723728
! *********************************************************************

USE um_types, ONLY: real_jlslsm
USE jules_vegetation_mod, ONLY: light_curvature_fvcb

IMPLICIT NONE

PRIVATE
PUBLIC :: jb_eta_scale, jb_vqmax25, jb_electron_flux

! The Farquhar curvature used to convert Jmax to Vqmax is
! light_curvature_fvcb (jules_vegetation namelist), as in calc_electron_flux.
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  jb_nl = 0.75,                                                                &
    ! Coupling efficiency of linear electron flow (mol ATP mol-1 e-).
  jb_nc = 1.0,                                                                 &
    ! Coupling efficiency of cyclic electron flow (mol ATP mol-1 e-).
  jb_qsat = 1800.0e-6,                                                         &
    ! Incident PAR at which Jmax is taken to have been measured, i.e. of
    ! a standard A-Ci curve (mol photons m-2 s-1).
  jb_leaf_abs = 0.85,                                                          &
    ! Leaf absorptance of PAR used with jb_qsat.
  jb_ci_ref = 800.0e-6,                                                        &
    ! Cc (mol mol-1) at which JB matches the Farquhar J (Lamour et al.).
  jb_ccp25 = 4.73078,                                                          &
    ! CO2 compensation point at 25 degC (Pa), as in sf_stom (Bernacchi).
  jb_p_ref = 101325.0,                                                         &
    ! Pressure (Pa) used to express jb_ci_ref as a partial pressure.
  jb_x_ref = jb_ccp25 / ( jb_ci_ref * jb_p_ref ),                              &
    ! G*/Cc at the reference state.
  jb_eta_ref = 1.0 - jb_nl / jb_nc                                             &
               + ( 3.0 + 7.0 * jb_x_ref )                                      &
                 / ( ( 4.0 + 8.0 * jb_x_ref ) * jb_nc )
    ! eta (JB Eqn 15c) at the reference state.

CONTAINS

!-----------------------------------------------------------------------------
! Return the factor that adjusts the Johnson-Berry electron transport rate
! to the leaf's internal CO2: it takes je (J at Cc = jb_ci_ref) to J at
! the leaf's internal CO2 ci,
!
!   scale = eta_ref / eta(ci)
!
! with eta as in the module header and ci no lower than the CO2
! compensation point ccp (both Pa). It is 1 at ci = jb_ci_ref, rises
! towards eta_ref at high ci and falls to eta_ref / 1.083 at ci = ccp.
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION jb_eta_scale( ccp, ci ) RESULT( scale )

REAL(KIND=real_jlslsm), INTENT(IN) :: ccp, ci
    ! CO2 compensation point and internal CO2 pressure (Pa).
REAL(KIND=real_jlslsm) :: scale

REAL(KIND=real_jlslsm) :: x, eta

x     = ccp / MAX( ci, ccp, TINY(ci) )
eta   = 1.0 - jb_nl / jb_nc + ( 3.0 + 7.0 * x ) / ( ( 4.0 + 8.0 * x ) * jb_nc )
scale = jb_eta_ref / eta

END FUNCTION jb_eta_scale

!-----------------------------------------------------------------------------
! Return Vqmax at 25 degC, the Johnson-Berry model's capacity for
! electron transport, from the Farquhar model's Jmax at 25 degC, so the
! existing Jmax parameters can be used (Lamour et al., 2026, Eqn 11).
!
! At the light of a measured A-Ci curve, the light reaching photosystem II
! as an electron flux is
!
!   i2sat = alpha_elec 0.85 1800e-6      (mol m-2 s-1)
!
! (quantum efficiency times leaf absorptance times incident light), and
! the Farquhar model gives Jsat from Jmax25 and the curvature
! light_curvature_fvcb. Vqmax25 is the value that gives the same J in the
! Johnson-Berry model:
!
!   Vqmax25 = Jsat eta_ref i2sat / (i2sat - Jsat)      (mol m-2 s-1)
!
! It is 0 if Jsat is not positive.
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION jb_vqmax25( jmax25, alpha_elec ) RESULT( vqmax25 )

REAL(KIND=real_jlslsm), INTENT(IN) :: jmax25, alpha_elec
    ! Jmax at 25 degC (mol m-2 s-1) and the quantum efficiency of
    ! electron transport (mol e- mol-1 absorbed photons).
REAL(KIND=real_jlslsm) :: vqmax25

REAL(KIND=real_jlslsm) :: i2sat, jsat

! Light to PS II, and the Farquhar J, at the A-Ci curve light level.
i2sat = alpha_elec * jb_leaf_abs * jb_qsat
jsat  = ( i2sat + jmax25                                                       &
          - SQRT( ( i2sat + jmax25 )**2                                        &
                  - 4.0 * light_curvature_fvcb * i2sat * jmax25 ) )            &
        / ( 2.0 * light_curvature_fvcb )

! jsat < i2sat always (it is the smaller root), so this is finite.
IF ( jsat > 0.0 ) THEN
  vqmax25 = jsat * jb_eta_ref * i2sat / ( i2sat - jsat )
ELSE
  vqmax25 = 0.0
END IF

END FUNCTION jb_vqmax25

!-----------------------------------------------------------------------------
! Work out the electron flux that drives photosynthesis in the
! Johnson-Berry model, at Cc = jb_ci_ref (mol m-2 s-1), from the light
! reaching photosystem II as an electron flux, i2, and Vqmax:
!
!   I1 = eta_ref i2
!   je = Vqmax I1 / (I1 + Vqmax) / eta_ref
!
! je is 0 where I1 + Vqmax is not positive.
!-----------------------------------------------------------------------------
SUBROUTINE jb_electron_flux( land_pts, veg_pts, veg_index, i2, vqmax, je )

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  i2(land_pts),                                                                &
    ! Light to PS II as an electron flux, alpha_elec * acr (mol m-2 s-1).
  vqmax(land_pts)
    ! Maximum Cyt b6f activity (mol e- m-2 s-1).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  je(land_pts)
    ! Electron transport rate at Cc = jb_ci_ref (mol m-2 s-1).

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: i1

DO m = 1,veg_pts
  l = veg_index(m)
  ! Light used for PS I photochemistry, matching the Farquhar quantum yield.
  i1 = jb_eta_ref * i2(l)
  IF ( i1 + vqmax(l) > 0.0 ) THEN
    je(l) = vqmax(l) * i1 / ( i1 + vqmax(l) ) / jb_eta_ref
  ELSE
    je(l) = 0.0
  END IF
END DO

END SUBROUTINE jb_electron_flux

END MODULE jb_photo_mod
