! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE soil_evap_or_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: gsoil_or, w1_or

! Water in soil layer 1 available to soil evaporation (kg m-2): liquid water
! above half the wilting point, as CABLE's lower limit on soil evaporation
! (cable_canopy.F90, flower_limit = wb(1) - swilt/2). Set in physiol each
! timestep with l_soil_evap_or and used by sf_evap to cap the soil
! evaporation alone (the Or resistance does not stop it in dry soil).
REAL(KIND=real_jlslsm), ALLOCATABLE :: w1_or(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOIL_EVAP_OR_MOD'

!-----------------------------------------------------------------------------
! Constants, as in CABLE cable_psm.F90 (Decker et al., 2017).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  rt_dff = 2.5e-5,                                                             &
      ! Diffusivity of water vapour in air (m2 s-1).
  lm = 1.73e-5,                                                                &
      ! Unit conversion constant in the liquid-supply resistance (m).
  rtevap_max = 10000.0,                                                        &
      ! Upper limit on the soil evaporation resistance (s m-1).
  dv_litt = 3.1415841138194147e-05,                                            &
      ! Vapour diffusivity through litter (m2 s-1), CABLE canopy%DvLitt
      ! (Matthews 2006; u = 1 m s-1, bulk litter density 63.5 kg m-3).
  sublayer_dz_guess = 0.005,                                                   &
      ! Height (m) at which the canopy wind profile is evaluated to set the
      ! eddy shape. CABLE resets canopy%sublayer_dz to this every call.
  pi_or = 3.14159,                                                             &
      ! CABLE's pi (pi_r_2), kept for bit-comparability of the formulae.
  ! CABLE roughness constants (cable_data.F90) used for the within-canopy
  ! wind extinction coefficient (rough%coexp).
  csd = 0.003,                                                                 &
      ! Substrate drag coefficient.
  crd = 0.3,                                                                   &
      ! Element drag coefficient.
  ccd = 15.0,                                                                  &
      ! Constant in d/h equation.
  ccw_c = 2.0,                                                                 &
      ! ccw = (zw-d)/(h-d).
  usuhm = 0.3,                                                                 &
      ! Maximum of us/uh.
  vonk_or = 0.40,                                                              &
      ! von Karman constant (CABLE value).
  grav_or = 9.81,                                                              &
      ! Gravitational acceleration (m s-2), CABLE value.
  ! CABLE constants for the below-canopy turbulent resistance rt0
  ! (cable_roughness.F90, cable_canopy.F90).
  csw = 0.50,                                                                  &
      ! Canopy sw decay (Weil theory).
  a33 = 1.25,                                                                  &
      ! Inertial sublayer sw/us.
  ctl = 0.40,                                                                  &
      ! Lagrangian timescale constant.
  zdlin = 1.0,                                                                 &
      ! Height fraction of d below which TL is linear.
  lai_thresh_or = 0.001,                                                       &
      ! CABLE LAI_THRESH: below it the surface is bare soil.
  rt_min_or = 5.0
      ! CABLE minimum turbulent resistance (s m-1).

CONTAINS

!-----------------------------------------------------------------------------
! Description:
!   Soil surface conductance for evaporation (m s-1) from the pore-scale
!   model of Or and co-workers (Haghighi et al., 2013; Haghighi & Or, 2015),
!   as coded in CABLE by Decker et al. (2017, JAMES) - cable_psm.F90,
!   SUBROUTINE or_soil_evap_resistance, used when cable_user%or_evap = T.
!
!   The resistance is the sum of a liquid-supply term (lm / 4K) and vapour
!   diffusion across a viscous sublayer plus the pore-scale boundary layer,
!   the sublayer depth set from the turbulent eddy spectrum at the soil
!   surface. Under a canopy CABLE's below-canopy turbulent resistance rt0
!   (cable_roughness.F90 rt0us, eq. 3.54 of the CSIRO SCAM manual;
!   cable_canopy.F90 ssnow%rtsoil = rt0) is added in series; over bare soil
!   JULES's own aerodynamic resistance (ch vshr) already plays CABLE's
!   rt0 + rt1 role. It returns the conductance of the unsaturated part of the
!   surface (CABLE's rtevap_unsat), i.e. a saturated fraction of zero, as
!   CABLE's groundwater configuration with real slopes. (Offline CABLE
!   without a groundwater elevation file has satfrac ~ 0.99, so nearly the
!   whole surface uses the moisture-independent rtevap_sat there.)
!
!   Differences from CABLE:
!   - The sublayer depth uses visc / u*_surface, dividing by the canopy
!     wind-profile factor; CABLE multiplies by it (see below).
!   - u* is supplied by the caller (JULES computes the tile u* after physiol,
!     so a neutral-stability estimate is used).
!   - K(theta) uses CABLE's Campbell form, K = Ks S**(2b+3), with the JULES
!     b and Ks, even when l_vg_soil = T. The JULES van Genuchten K for the
!     same b is far smaller in the top layer (e.g. ~50x at S = 0.7 for
!     FR-Pue, b = 6.7), which pins the resistance at rtevap_max almost
!     permanently.
!   - No snow adjustments (JULES handles snow separately); litter is a
!     series resistance (see litter_dz).
!   - Residual water content watr = 0 (CABLE default 0.05).
!   - The sublayer depth is evaluated each call from a fixed first guess
!     (CABLE carries canopy%sublayer_dz over from the previous step).
!   - rt0 uses CABLE's default (non-SLI) formulation, neutral stability, and
!     the same z0soil (scaled by z0soil_fac) as the Or resistance.
!   - The soil surface relative humidity (CABLE rh_srf) is not applied: with
!     CABLE's -10 m cap on the suction it is >= 0.9992, i.e. ~1.
!   - The eddy_mod gamma-function ratio is evaluated in log space (same
!     maths, avoids overflow for large eddy shapes).
!-----------------------------------------------------------------------------
FUNCTION gsoil_or(ua, ustar, hc, lai, t_k, sathh, bexp, satcon,          &
                  theta_liq, theta_sat, litter_dz, z0soil_fac)

USE jules_hydrology_mod, ONLY: or_k_min

IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ua,                                                                          &
      ! Wind speed at the reference height (m s-1).
  ustar,                                                                       &
      ! Friction velocity above the surface (m s-1).
  hc,                                                                          &
      ! Canopy height (m); 0 for bare soil.
  lai,                                                                         &
      ! Leaf area index above the soil (m2 m-2).
  t_k,                                                                         &
      ! Air/surface temperature for the kinematic viscosity (K).
  sathh,                                                                       &
      ! Absolute saturated soil water suction of the top layer (m).
  bexp,                                                                        &
      ! Clapp-Hornberger exponent of the top layer.
  satcon,                                                                      &
      ! Saturated hydraulic conductivity of the top layer (kg m-2 s-1).
  theta_liq,                                                                   &
      ! Volumetric liquid water content of the top layer (m3 m-3).
  theta_sat,                                                                   &
      ! Volumetric water content at saturation of the top layer (m3 m-3).
  litter_dz,                                                                   &
      ! Litter layer depth (m); its diffusion resistance litter_dz / dv_litt
      ! is added in series, as CABLE's default-scheme litter resistance
      ! (relitt). CABLE's Or scheme instead adds the litter depth to the
      ! viscous sublayer, which, through the z0soil / sublayer_dz factor,
      ! lowers the resistance; not followed here.
  z0soil_fac
      ! Multiplier on z0soil (1 = CABLE).

REAL(KIND=real_jlslsm) :: gsoil_or
      ! Soil surface conductance (m s-1).

REAL(KIND=real_jlslsm) ::                                                      &
  hruff, usuh, xx, dh, coexp, visc, us, us_surf, z0soil, eddy_shape, log_em,            &
  eddy_mod, sublayer_dz, wb_liq, rel_s, hk_zero, soil_moisture_mod, pore_radius,      &
  pore_size, rtevap, disp, term2, term3, term5, rt0us, rt0

INTEGER :: k, int_eddy_shape

!-----------------------------------------------------------------------------
! Within-canopy wind extinction coefficient (CABLE cable_roughness.F90).
!-----------------------------------------------------------------------------
hruff = MAX(1.0e-6, hc)
usuh  = MIN(SQRT(csd + crd * (lai * 0.5)), usuhm)
xx    = SQRT(ccd * MAX(lai * 0.5, 0.0005))
dh    = 1.0 - (1.0 - EXP(-xx)) / xx
coexp = usuh / (vonk_or * ccw_c * (1.0 - dh))

! Kinematic viscosity of air (CABLE cable_air.F90).
visc = 1.0e-5 * MAX(1.0, 1.35 + 0.0092 * (t_k - 273.15))

us = MAX(1.0e-3, ustar)

! Soil roughness length used by CABLE when or_evap = T.
z0soil = z0soil_fac * ( 0.01 * MIN(1.0, lai) + 0.02 * MIN(us**2 / grav_or, 1.0) )

!-----------------------------------------------------------------------------
! Viscous sublayer depth from the eddy spectrum at the surface.
!-----------------------------------------------------------------------------
! Friction velocity at the soil surface, from the exponential within-canopy
! wind profile evaluated sublayer_dz_guess above the ground.
us_surf = us * EXP(-coexp * (1.0 - sublayer_dz_guess / MAX(1.0e-2, hruff)))

eddy_shape = 0.3 * ua / MAX(1.0e-4, us_surf)
int_eddy_shape = FLOOR(eddy_shape)

log_em = LOG(2.2 * SQRT(112.0 * pi_or)) - (eddy_shape + 1.0) * LOG(2.0)       &
         - 0.5 * LOG(eddy_shape + 1.0)
IF (int_eddy_shape > 0) THEN
  log_em = log_em - LOG_GAMMA(eddy_shape + 1.0) + LOG(2.0 * eddy_shape + 1.0)
  DO k = 1, int_eddy_shape
    log_em = log_em + LOG(2.0 * (eddy_shape - REAL(k)) + 1.0)
  END DO
END IF
eddy_mod = EXP(log_em)

! Sublayer depth scales with visc / u* at the soil surface. CABLE
! (cable_psm.F90) multiplies visc/u* by the wind-profile factor, i.e. uses
! u* / factor, which makes the sublayer thinner under a canopy instead of
! thicker; here it divides, consistent with eddy_shape above.
sublayer_dz = MIN(0.05, MAX(eddy_mod * visc / MAX(1.0e-4, us_surf), 1.0e-7))

!-----------------------------------------------------------------------------
! Pore-scale resistance for the unsaturated surface.
!-----------------------------------------------------------------------------
wb_liq = MAX(0.0001, MIN(pi_or / 4.0, theta_liq))

! Relative saturation (watr = 0; CABLE's default is 0.05) and Campbell
! conductivity.
! kg m-2 s-1 == mm s-1, so 0.001 converts to m s-1 as CABLE does for hyds.
rel_s   = MAX(wb_liq, 0.0) / theta_sat
hk_zero = MAX(0.001 * satcon * (MIN(MAX(rel_s, 0.001), 1.0)                  &
              **(2.0 * bexp + 3.0)), or_k_min)

soil_moisture_mod = 1.0 / pi_or / SQRT(wb_liq)                                 &
                    * (SQRT(pi_or / (4.0 * wb_liq)) - 1.0)

! Pore radius from the air-entry suction (0.148 = 2 * surface tension).
pore_radius = 0.148 / (1000.0 * grav_or * ABS(sathh))
pore_size   = pore_radius * SQRT(pi_or)

IF (sublayer_dz >= 1.0e-7) THEN
  rtevap = MIN(rtevap_max, z0soil / sublayer_dz * (lm / (4.0 * hk_zero)       &
           + (sublayer_dz + pore_size * soil_moisture_mod) / rt_dff))
ELSE
  rtevap = MIN(rtevap_max, lm / (4.0 * hk_zero)                                &
           + (sublayer_dz + pore_size * soil_moisture_mod) / rt_dff)
END IF

!-----------------------------------------------------------------------------
! Below-canopy turbulent resistance (CABLE default scheme, neutral):
! rt0 = max(rt_min, rt0us / u*), rt0us from cable_roughness.F90 (eq. 3.54,
! CSIRO SCAM manual). Only under a canopy (lai >= LAI_THRESH and the canopy
! taller than z0soil); over bare soil JULES's ch vshr is the full path.
!-----------------------------------------------------------------------------
rt0 = 0.0
IF ( lai >= lai_thresh_or .AND. hruff >= MAX(z0soil, 1.0e-7) ) THEN
  disp  = dh * hruff
  term2 = EXP(2.0 * csw * lai * (1.0 - disp / hruff))
  term3 = a33**2 * ctl * 2.0 * csw * lai
  term5 = MAX((2.0 / 3.0) * hruff / disp, 1.0)
  rt0us = term5 * (zdlin * LOG(zdlin * disp / MAX(z0soil, 1.0e-7))           &
                   + (1.0 - zdlin))                                            &
          * (EXP(2.0 * csw * lai) - term2) / term3
  rt0   = MAX(rt_min_or, rt0us / us)
END IF

gsoil_or = 1.0 / (rt0 + rtevap + MAX(litter_dz, 0.0) / dv_litt)

END FUNCTION gsoil_or

END MODULE soil_evap_or_mod
