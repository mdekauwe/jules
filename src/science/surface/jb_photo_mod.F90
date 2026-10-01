! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE jb_photo_mod

! *********************************************************************
! Johnson & Berry (2021) electron transport (photo_model = photo_johnson).
!
! The Rubisco-limited rate, the CO2 compensation point and all the
! temperature responses are those of the Farquhar model; only the light
! response of electron transport changes. The empirical non-rectangular
! hyperbola in Jmax and a curvature is replaced by electron flow through
! Cyt b6f (JB Eqns 15c, 30a, 30c):
!   JP700 = Vqmax I1 / (I1 + Vqmax)          PS I electron flux
!   J     = JP700 / eta(Cc)                  PS II (linear) electron flux
!   eta   = 1 - nL/nC + (3 + 7 G*/Cc) / ((4 + 8 G*/Cc) nC)
! where Vqmax is the maximum Cyt b6f activity, I1 the light used for PS I
! photochemistry and G* the CO2 compensation point. J then enters
! Aj = J/4 (Cc - G*)/(Cc + 2 G*) as in Farquhar.
!
! Parameters follow Lamour et al. (2026, bioRxiv 2026.05.08.723728), who
! harmonised JB with the FvCB parameters used in vegetation models:
!   (i)  I1 is set so the initial quantum yield of J equals that of the
!        Farquhar model (their Eqn 10) at Cc = jb_ci_ref:
!        I1 = eta_ref alpha_elec acr, i.e. J = i2 eta_ref / eta at low light.
!   (ii) Vqmax at 25 degC is set so J of JB equals J of FvCB at the light
!        and CO2 of an A-Ci curve (their Eqn 11, Qsat = 1800 umol m-2 s-1,
!        leaf absorptance 0.85, Ci = 800 ppm; supplement and ELM-FATES
!        set-up), from the PFT's Jmax at 25 degC:
!        Vqmax25 = Jsat eta_ref i2sat / (i2sat - Jsat).
!   (iii) Vqmax has the temperature response of Jmax (as in ELM-FATES).
! This lets the existing Vcmax/Jmax parameters be used directly; Vqmax
! should ultimately be recalibrated (e.g. by fitting JB to A-Ci curves).
!
! In JULES the jmax array carries Vqmax and je carries JP700 / eta_ref
! (= J at Cc = jb_ci_ref), so that the sun/shade and layer ratios of je,
! and the code that passes je around, are unchanged. The Cc dependence is
! applied where wlite is formed: wlite uses je * jb_eta_scale(ccp, ci).
! Mesophyll conductance is infinite, so Cc = ci.
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
! eta_ref / eta(Cc): the factor taking je (J at Cc = jb_ci_ref) to J at ci.
! Lies between 1 and eta_ref / 1.125.
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
! Vqmax at 25 degC from Jmax at 25 degC (Lamour et al. Eqn 11).
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
! je = JP700 / eta_ref, the JB electron flux at Cc = jb_ci_ref.
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
