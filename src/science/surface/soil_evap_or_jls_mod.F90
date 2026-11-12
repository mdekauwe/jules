! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE soil_evap_or_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: gsoil_or, w_layer1_evap_or

! Water in soil layer 1 available to soil evaporation (kg m-2): liquid water
! above half the wilting point, as CABLE's lower limit on soil evaporation
! (cable_canopy.F90, flower_limit = wb(1) - swilt/2). Set in physiol each
! timestep with l_soil_evap_or and used by sf_evap to cap the soil
! evaporation alone (the Or resistance does not stop it in dry soil).
!
! MDK, 12/11/26
REAL(KIND=real_jlslsm), ALLOCATABLE :: w_layer1_evap_or(:,:)

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
! Return the conductance of the soil surface to evaporation (m s-1),
! which limits soil evaporation as the surface dries, from the
! pore-scale model of Or and co-workers (Haghighi et al., 2013; Haghighi
! and Or, 2015) as coded in CABLE by Decker et al. (2017, JAMES):
! cable_psm.F90, or_soil_evap_resistance, used when cable_user%or_evap = T.
!
! Water leaving the soil meets three resistances in series (s m-1):
!
!   1. Inside the soil (rtevap). Liquid water must reach the surface,
!      lm / (4 K), where K is the hydraulic conductivity of the top layer
!      (m s-1) and lm a constant (m). Then vapour diffuses out of the
!      surface pores and across a thin layer of still air (the viscous
!      sublayer) just above the soil. The sublayer's depth is set by the
!      size of the turbulent eddies that reach the soil surface, which
!      depends on the wind and the friction velocity there. Both parts are
!      scaled by z0soil / (sublayer depth), z0soil being the soil's
!      roughness length, and the sum is capped at rtevap_max.
!   2. Through the air under the canopy (rt0), CABLE's below-canopy
!      turbulent resistance (cable_roughness.F90 rt0us, eq. 3.54 of the
!      CSIRO SCAM manual), only where there is a canopy. Over bare soil,
!      JULES's own aerodynamic resistance already covers this path.
!   3. Through the litter, litter_dz / dv_litt, where litter_dz is the
!      litter depth (m) and dv_litt the vapour diffusivity through litter
!      (m2 s-1).
!
!   gsoil_or = 1 / (rt0 + rtevap + litter_dz / dv_litt)
!
! It is the conductance of the unsaturated part of the soil surface
! (CABLE's rtevap_unsat), as if none of the surface were saturated, as in
! CABLE's groundwater set-up with real slopes. Offline CABLE without a
! groundwater elevation file treats about 99 % of the surface as saturated
! and uses a resistance that does not depend on soil moisture there.
!
! Differences from CABLE:
!   - The sublayer depth scales with the air's viscosity divided by the
!     friction velocity at the soil surface, which is lower under a
!     canopy, so the sublayer is thicker there (Decker et al., 2017,
!     eq. 17 and 21, and CABLE CMIP6-GM2). CABLE trunk and Decker's later
!     branch multiply by the canopy wind factor instead of dividing.
!   - The friction velocity comes from the caller, as a neutral-stability
!     estimate, since JULES works out the tile's value after physiol.
!   - K uses CABLE's Campbell form, K = Ks S**(2b+3), with JULES's b and
!     saturated conductivity Ks and relative saturation S, even when
!     l_vg_soil = T. The JULES van Genuchten K for the same b is much
!     smaller in the top layer, which would hold the resistance at
!     rtevap_max nearly all the time. K is kept no lower than or_k_min.
!   - Snow is left to the rest of JULES. Litter is a separate resistance
!     in series (see litter_dz).
!   - The residual water content is 0 (CABLE's default is 0.05).
!   - The sublayer depth is worked out each call from a fixed first guess
!     (CABLE carries it over from the previous step).
!   - rt0 uses CABLE's default (non-SLI) form, neutral stability, and the
!     same z0soil (times z0soil_fac) as the soil resistance.
!   - The relative humidity at the soil surface (CABLE rh_srf) is taken as
!     1: with CABLE's -10 m cap on the suction it is at least 0.9992.
!   - The eddy-size factor's gamma-function ratio is worked out with
!     logarithms, which gives the same result without overflow for large
!     eddy shapes.
!
! MDK, 12/11/26
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

! The names of the locals follow CABLE's cable_psm.F90 and
! cable_roughness.F90, for line-by-line comparison with CABLE.
! Canopy height with a floor of 1e-6 m (m; CABLE name).
REAL(KIND=real_jlslsm) :: hruff
! u* / u at the canopy top, capped at usuhm (CABLE name).
REAL(KIND=real_jlslsm) :: usuh
! Argument of the d/h equation, SQRT(ccd LAI / 2) (CABLE name).
REAL(KIND=real_jlslsm) :: xx
! Displacement height over canopy height, d/h (CABLE name).
REAL(KIND=real_jlslsm) :: dh
! Within-canopy wind extinction coefficient (CABLE name).
REAL(KIND=real_jlslsm) :: coexp
! Kinematic viscosity of air (m2 s-1; CABLE name).
REAL(KIND=real_jlslsm) :: visc
! Friction velocity with a floor of 1e-3 (m s-1; CABLE name).
REAL(KIND=real_jlslsm) :: us
! Friction velocity at the soil surface, from the exponential
! within-canopy wind profile (m s-1).
REAL(KIND=real_jlslsm) :: us_surf
! Soil roughness length (m; CABLE name).
REAL(KIND=real_jlslsm) :: z0soil
! Shape parameter of the eddy-size distribution at the soil
! surface (CABLE name).
REAL(KIND=real_jlslsm) :: eddy_shape
! Logarithm of the eddy-size factor's gamma-function ratio
! (CABLE's eddy_mod, worked out in logs).
REAL(KIND=real_jlslsm) :: log_em
! Eddy-size factor on the viscous sublayer depth (CABLE name).
REAL(KIND=real_jlslsm) :: eddy_mod
! Depth of the viscous sublayer above the soil (m; CABLE
! name).
REAL(KIND=real_jlslsm) :: sublayer_dz
! Liquid water content of the top layer, bounded to (1e-4, pi/4)
! (m3 m-3; CABLE name).
REAL(KIND=real_jlslsm) :: wb_liq
! Relative saturation of the top layer (-).
REAL(KIND=real_jlslsm) :: rel_s
! Campbell hydraulic conductivity of the top layer, floored at
! or_k_min (m s-1; CABLE name).
REAL(KIND=real_jlslsm) :: hk_zero
! Moisture factor on the pore size in rtevap (CABLE name).
REAL(KIND=real_jlslsm) :: soil_moisture_mod
! Pore radius from the air-entry suction (m; CABLE name).
REAL(KIND=real_jlslsm) :: pore_radius
! Pore size, pore_radius SQRT(pi) (m; CABLE name).
REAL(KIND=real_jlslsm) :: pore_size
! Resistance inside the soil and its viscous sublayer, capped at
! rtevap_max (s m-1; CABLE name).
REAL(KIND=real_jlslsm) :: rtevap
! Displacement height, dh hruff (m; CABLE name).
REAL(KIND=real_jlslsm) :: disp
! EXP(2 csw LAI (1 - d / h)), a term of rt0us (CABLE name).
REAL(KIND=real_jlslsm) :: term2
! a33**2 ctl 2 csw LAI, the denominator of rt0us (CABLE name).
REAL(KIND=real_jlslsm) :: term3
! Height factor MAX(2 h / (3 d), 1) of rt0us (CABLE name).
REAL(KIND=real_jlslsm) :: term5
! Below-canopy turbulent resistance times u* (m s-1 s m-1;
! CABLE name, eq. 3.54 of the CSIRO SCAM manual).
REAL(KIND=real_jlslsm) :: rt0us
! Below-canopy turbulent resistance (s m-1; CABLE name).
REAL(KIND=real_jlslsm) :: rt0
! Loop counter over the integer part of eddy_shape.
INTEGER :: k
! Integer part of eddy_shape.
INTEGER :: int_eddy_shape

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

! Sublayer depth scales with visc / u* at the soil surface (Decker et al.
! 2017 eq. 17, with u* scaled down beneath the canopy, eq. 21), so the
! sublayer is thicker under a canopy, consistent with eddy_shape above.
! (CABLE trunk multiplies visc/u* by the wind-profile factor instead.)
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
