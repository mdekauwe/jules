! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE leaf_energy_balance_mod

!-----------------------------------------------------------------------------
! Leaf energy balance of the sunlit and shaded leaf classes (l_leaf_temp):
! the canopy aerodynamics (leaf_aero_geom), the terms fixed over the passes
! of the stomatal optimisation (leaf_temp_setup) and the leaf and canopy air
! temperature update (leaf_temp_update, leaf_temp_update_n). Called by
! sf_stom.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEAF_ENERGY_BALANCE_MOD'

PRIVATE
PUBLIC :: leaf_aero_geom, leaf_temp_setup, leaf_temp_update,                  &
          leaf_temp_update_n

CONTAINS

SUBROUTINE leaf_aero_geom( land_pts, veg_pts, veg_index, lai, canht, z1,       &
                           z0m, zref, rt1ab, zrd, usc, usuh, coexp )

! CABLE canopy geometry (cable_roughness; Raupach, 1994, BLM 71: 211-216;
! Raupach et al., 1997, CSIRO SCAM) for the leaf energy balance
! (l_leaf_temp, leaf_aero_model = 1): displacement d and roughness z0m from
! the leaf area, the reference height zref above d, u*/u_h (usuh), the wind
! extinction coefficient coexp, the roughness-sublayer terms of rt1 (rt1ab =
! rt1usa + rt1usb, times 1/u*), and zrd = zruffs - d and usc (1 if zref is
! above the roughness sublayer) for rt1usc. z1 is taken as height above the
! ground. Points without a canopy get zref = 0 (and are skipped).
! coexp is returned per unit cumulative leaf area: as CABLE's gbhu uses it,
! or with l_leaf_coexp_lai as derived, per normalised height (divided by
! the leaf area).
!
! MDK, 24/10/26

USE jules_vegetation_mod, ONLY: l_leaf_coexp_lai

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), canht(land_pts),          &
                                      z1(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: z0m(land_pts), zref(land_pts),          &
                                       rt1ab(land_pts), zrd(land_pts),         &
                                       usc(land_pts), usuh(land_pts),          &
                                       coexp(land_pts)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  ! CABLE physical constants (cable_data).
  vonk = 0.40, a33 = 1.25, csw = 0.50, ctl = 0.40, crd = 0.3, csd = 0.003,     &
  ccd = 15.0, ccw_c = 2.0, usuhm = 0.3

REAL(KIND=real_jlslsm) :: xx, dh, disp, h, lai_h, term2, term3, term5, zruffs
INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_AERO_GEOM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

z0m(:)   = 1.0
zref(:)  = 0.0
rt1ab(:) = 0.0
zrd(:)   = 1.0
usc(:)   = 0.0
usuh(:)  = usuhm
coexp(:) = 0.0

DO m = 1,veg_pts
  l = veg_index(m)
  h = canht(l)
  IF ( lai(l) < 0.01 .OR. h <= 0.0 ) CYCLE
  lai_h = 0.5 * lai(l)

  usuh(l)  = MIN( SQRT( csd + crd * lai_h ), usuhm )
  xx       = SQRT( ccd * MAX( lai_h, 0.0005 ) )
  dh       = 1.0 - ( 1.0 - EXP( -xx ) ) / xx
  disp     = dh * h
  z0m(l)   = ( (1.0 - dh) * EXP( LOG( ccw_c ) - 1.0 + 1.0 / ccw_c              &
                                 - vonk / usuh(l) ) ) * h
  zref(l)  = MAX( 3.5 + z0m(l), z1(l) - disp, h - disp )
  coexp(l) = usuh(l) / ( vonk * ccw_c * ( 1.0 - dh ) )
  IF ( l_leaf_coexp_lai ) coexp(l) = coexp(l) / lai(l)

  term2    = EXP( 2.0 * csw * lai(l) * ( 1.0 - disp / h ) )
  term3    = a33**2 * ctl * 2.0 * csw * lai(l)
  term5    = MAX( (2.0 / 3.0) * h / disp, 1.0 )
  zruffs   = disp + h * a33**2 * ctl / vonk / term5
  rt1ab(l) = term5 * ( term2 - 1.0 ) / term3                                   &
             + MAX( term5 * ( MIN( zref(l) + disp, zruffs ) - h )              &
                    / ( a33**2 * ctl * h ), 0.0 )
  zrd(l)   = MAX( zruffs - disp, z0m(l) )
  IF ( zref(l) + disp > zruffs ) usc(l) = 1.0
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_aero_geom

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,    &
                            dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,    &
                            wind_ext,                                          &
                            apar_sun, apar_shd,                                &
                            sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,        &
                            wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,      &
                            fapar_sun_vis, fapar_shd_vis, lai_c, sw_c, wlw_c,  &
                            wlws_c, frad_c, ff_c )

! Leaf energy balance (l_leaf_temp; two-leaf, and the multilayer canopy's
! classes): the terms that are fixed over
! the passes of the stomatal optimisation, summed over the canopy layers
! with the sunlit fraction of each layer (as the absorbed PAR is):
! - sw: absorbed PAR + near-IR (W m-2 ground). The incident near-IR is
!   taken as equal to the incident PAR (standalone sets PAR = sw_down/2).
! - wlw, wlws: layer weights of the isothermal long-wave (Wang & Leuning,
!   1998, B18/B19): the sky deficit absorbed from above, kd exp(-kd xi), and
!   the part reflected by the soil, kd exp(-kd (L - xi)), integrated over
!   each layer.
! - frad: the radiative conductance scalar, kd [exp(-kd xi) +
!   exp(-kd (L - xi))] integrated over the class (2 x leaf area for a thin
!   canopy).
! - ff: leaf area weighted by the relative wind speed for forced
!   convection, exp(-wind_ext xi / 2) (Wang & Leuning, 1998; CABLE);
!   the class leaf area when wind_ext = 0.
! With the optional arguments (the multilayer canopy with ml_leaf_temp_model
! = 2), also the same terms of every layer's sunlit and shaded leaf, as
! classes 2n-1 (sunlit) and 2n (shaded) of layer n, the PAR part from the
! layer's own absorbed PAR (fapar_*_vis); the class leaf area too. The
! layer terms sum to the canopy ones.
!
! MDK, 24/10/26

USE pftparm, ONLY: orient

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts), ilayers

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  lai(land_pts), dlai(land_pts), fsun(land_pts,ilayers),                       &
  fapar_sun_nir(land_pts,ilayers), fapar_shd_nir(land_pts,ilayers),            &
  ipar(land_pts), apar_sun(land_pts), apar_shd(land_pts), wind_ext(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  sw_sun(land_pts), sw_shd(land_pts), wlw_sun(land_pts), wlw_shd(land_pts),    &
  wlws_sun(land_pts), wlws_shd(land_pts), frad_sun(land_pts),                  &
  frad_shd(land_pts), ff_sun(land_pts), ff_shd(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  fapar_sun_vis(land_pts,ilayers), fapar_shd_vis(land_pts,ilayers)
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  lai_c(land_pts,2*ilayers), sw_c(land_pts,2*ilayers),                         &
  wlw_c(land_pts,2*ilayers), wlws_c(land_pts,2*ilayers),                       &
  frad_c(land_pts,2*ilayers), ff_c(land_pts,2*ilayers)

REAL(KIND=real_jlslsm) :: kd, xi0, xi1, dtop, dbot, dwind
INTEGER :: l, m, n
LOGICAL :: l_lyr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_SETUP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Diffuse (long-wave) extinction coefficient: spherical leaves (Wang &
! Leuning, 1998), or horizontal.
IF ( orient(ft) == 1 ) THEN
  kd = 1.0
ELSE
  kd = 0.78
END IF

sw_sun(:)   = 0.0
sw_shd(:)   = 0.0
wlw_sun(:)  = 0.0
wlw_shd(:)  = 0.0
wlws_sun(:) = 0.0
wlws_shd(:) = 0.0
frad_sun(:) = 0.0
frad_shd(:) = 0.0
ff_sun(:)   = 0.0
ff_shd(:)   = 0.0
l_lyr = PRESENT(lai_c)
IF ( l_lyr ) THEN
  lai_c(:,:)  = 0.0
  sw_c(:,:)   = 0.0
  wlw_c(:,:)  = 0.0
  wlws_c(:,:) = 0.0
  frad_c(:,:) = 0.0
  ff_c(:,:)   = 0.0
END IF

DO m = 1,veg_pts
  l = veg_index(m)
  sw_sun(l) = apar_sun(l)
  sw_shd(l) = apar_shd(l)
  DO n = 1,ilayers
    xi0  = REAL(n-1) * dlai(l)
    xi1  = REAL(n) * dlai(l)
    dtop = EXP(-kd * xi0) - EXP(-kd * xi1)
    dbot = EXP(-kd * (lai(l) - xi1)) - EXP(-kd * (lai(l) - xi0))

    sw_sun(l)   = sw_sun(l)                                                    &
                  + fapar_sun_nir(l,n) * fsun(l,n) * dlai(l) * ipar(l)
    sw_shd(l)   = sw_shd(l)                                                    &
                  + fapar_shd_nir(l,n) * (1.0 - fsun(l,n)) * dlai(l) * ipar(l)
    wlw_sun(l)  = wlw_sun(l)  + fsun(l,n) * dtop
    wlw_shd(l)  = wlw_shd(l)  + (1.0 - fsun(l,n)) * dtop
    wlws_sun(l) = wlws_sun(l) + fsun(l,n) * dbot
    wlws_shd(l) = wlws_shd(l) + (1.0 - fsun(l,n)) * dbot
    frad_sun(l) = frad_sun(l) + fsun(l,n) * (dtop + dbot)
    frad_shd(l) = frad_shd(l) + (1.0 - fsun(l,n)) * (dtop + dbot)
    IF ( wind_ext(l) > 0.0 ) THEN
      dwind = 2.0 / wind_ext(l) * ( EXP(-0.5 * wind_ext(l) * xi0)              &
                                    - EXP(-0.5 * wind_ext(l) * xi1) )
    ELSE
      dwind = dlai(l)
    END IF
    ff_sun(l)   = ff_sun(l)   + fsun(l,n) * dwind
    ff_shd(l)   = ff_shd(l)   + (1.0 - fsun(l,n)) * dwind
    IF ( l_lyr ) THEN
      lai_c(l,2*n-1)  = fsun(l,n) * dlai(l)
      lai_c(l,2*n)    = (1.0 - fsun(l,n)) * dlai(l)
      sw_c(l,2*n-1)   = ( fapar_sun_vis(l,n) + fapar_sun_nir(l,n) )            &
                        * fsun(l,n) * dlai(l) * ipar(l)
      sw_c(l,2*n)     = ( fapar_shd_vis(l,n) + fapar_shd_nir(l,n) )            &
                        * (1.0 - fsun(l,n)) * dlai(l) * ipar(l)
      wlw_c(l,2*n-1)  = fsun(l,n) * dtop
      wlw_c(l,2*n)    = (1.0 - fsun(l,n)) * dtop
      wlws_c(l,2*n-1) = fsun(l,n) * dbot
      wlws_c(l,2*n)   = (1.0 - fsun(l,n)) * dbot
      frad_c(l,2*n-1) = fsun(l,n) * (dtop + dbot)
      frad_c(l,2*n)   = (1.0 - fsun(l,n)) * (dtop + dbot)
      ff_c(l,2*n-1)   = fsun(l,n) * dwind
      ff_c(l,2*n)     = (1.0 - fsun(l,n)) * dwind
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_temp_setup

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             ra, lw_down, u_wind, lai_sun, lai_shd,            &
                             sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,       &
                             wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,     &
                             z0m_a, zref_a, rt1ab_a, zrd_a, usc_a, usuh_a,     &
                             gl_sun, gl_shd, g_soil, qs_soil, zeta_a, t_sun,   &
                             t_shd, t_c, q_c, gbw_sun, gbw_shd, ra_out,        &
                             n_clamp, lt_a, lt_b )

! Two-leaf leaf energy balance (l_leaf_temp): sunlit and shaded leaf
! temperatures, and canopy air temperature and humidity, for the stomatal
! conductances of the last pass. The balance of leaf_temp_update_n for
! the two classes (sunlit, shaded).
!
! MDK, 24/10/26

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tair(land_pts), q1(land_pts), pstar(land_pts), ra(land_pts),                 &
  lw_down(land_pts), u_wind(land_pts), lai_sun(land_pts), lai_shd(land_pts),   &
  sw_sun(land_pts), sw_shd(land_pts), wlw_sun(land_pts), wlw_shd(land_pts),    &
  wlws_sun(land_pts), wlws_shd(land_pts), frad_sun(land_pts),                  &
  frad_shd(land_pts), ff_sun(land_pts), ff_shd(land_pts), gl_sun(land_pts),    &
  gl_shd(land_pts), z0m_a(land_pts), zref_a(land_pts), rt1ab_a(land_pts),      &
  zrd_a(land_pts), usc_a(land_pts), usuh_a(land_pts), g_soil(land_pts),        &
  qs_soil(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  zeta_a(land_pts), t_sun(land_pts), t_shd(land_pts), t_c(land_pts),           &
  q_c(land_pts), gbw_sun(land_pts), gbw_shd(land_pts), ra_out(land_pts),       &
  n_clamp(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT), OPTIONAL ::                            &
  ! The sunlit (1) and shaded (2) leaves' temperature as a linear
  ! function of their transpiration (see leaf_temp_update_n;
  ! l_leaf_temp_ci).
  lt_a(land_pts,2), lt_b(land_pts,2)

REAL(KIND=real_jlslsm) ::                                                      &
  lai_x(land_pts,2), sw_x(land_pts,2), wlw_x(land_pts,2), wlws_x(land_pts,2),  &
  frad_x(land_pts,2), ff_x(land_pts,2), gl_x(land_pts,2), t_x(land_pts,2),     &
  gbw_x(land_pts,2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lai_x(:,1)  = lai_sun(:)
lai_x(:,2)  = lai_shd(:)
sw_x(:,1)   = sw_sun(:)
sw_x(:,2)   = sw_shd(:)
wlw_x(:,1)  = wlw_sun(:)
wlw_x(:,2)  = wlw_shd(:)
wlws_x(:,1) = wlws_sun(:)
wlws_x(:,2) = wlws_shd(:)
frad_x(:,1) = frad_sun(:)
frad_x(:,2) = frad_shd(:)
ff_x(:,1)   = ff_sun(:)
ff_x(:,2)   = ff_shd(:)
gl_x(:,1)   = gl_sun(:)
gl_x(:,2)   = gl_shd(:)
t_x(:,1)    = t_sun(:)
t_x(:,2)    = t_shd(:)
gbw_x(:,1)  = gbw_sun(:)
gbw_x(:,2)  = gbw_shd(:)

CALL leaf_temp_update_n( land_pts, veg_pts, veg_index, 2, tair, q1, pstar,     &
                         ra, lw_down, u_wind, lai_x, sw_x, wlw_x, wlws_x,      &
                         frad_x, ff_x, z0m_a, zref_a, rt1ab_a, zrd_a, usc_a,   &
                         usuh_a, gl_x, g_soil, qs_soil, zeta_a, t_x, t_c,      &
                         q_c, gbw_x, ra_out, n_clamp, lt_a, lt_b )

t_sun(:)   = t_x(:,1)
t_shd(:)   = t_x(:,2)
gbw_sun(:) = gbw_x(:,1)
gbw_shd(:) = gbw_x(:,2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_temp_update

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_update_n( land_pts, veg_pts, veg_index, nc, tair, q1,     &
                               pstar, ra, lw_down, u_wind, lai_c, sw_c, wlw_c, &
                               wlws_c, frad_c, ff_c, z0m_a, zref_a, rt1ab_a,   &
                               zrd_a, usc_a, usuh_a, gl_c, g_soil, qs_soil,    &
                               zeta_a, t_cls, t_c, q_c, gbw_c, ra_out,         &
                               n_clamp, lt_a_c, lt_b_c )

! Leaf energy balance (l_leaf_temp) of nc leaf classes sharing one canopy
! air: two-leaf (nc = 2, sunlit and shaded), the aggregate multilayer canopy
! (the same) or the multilayer canopy layer by layer (nc = 2 x ilayers, the
! sunlit and shaded leaf of every layer; ml_leaf_temp_model = 2, as
! MAESPA). Leaf temperatures, and canopy air temperature and humidity, for
! the stomatal conductances of the last pass.
!
! Each leaf class exchanges with the canopy air (t_c, q_c) through its
! boundary layer, as a big leaf: Penman-Monteith with the isothermal net
! radiation, linearised about t_c (Leuning et al., 1995, App. E; Wang &
! Leuning, 1998), as in the two-leaf model of De Kauwe (two_leaf_at_WTC):
!
!   LE = (s Rn + D gh cp) / (s + gamma gh / gw),
!   T  = t_c + (Rn - LE) / (cp gh),   gh = 2 gbH + grn,
!
! with forced convection for gbH (per unit leaf area at the canopy-top wind
! u_wind, divided by leaf_shelter, times the class leaf area weighted by
! the wind profile in the canopy, ff) and free convection (times the class
! leaf area), and the radiative conductance grn.
!
! The canopy air is coupled to level 1 through ra (from leaf_aero_model:
! 0, so t_c = tair and q_c = q1; 1, CABLE's rt1; 2, JULES's ra):
!
!   q_c = q1 + ra (LE_sun + LE_shd + LE_soil) / (rho L),
!   t_c = tair + ra (H_sun + H_shd) / (rho cp),
!
! with L the latent heat at t_c and LE_soil = rho L g_soil (qs_soil - q_c)
! the soil evaporation into the canopy air, as in CABLE (soil sensible heat
! is not included: there is no soil surface temperature separate from tstar
! on vegetated tiles).
! The balance is solved for these conductances: LE is linear in q_c, so q_c
! is explicit for a given t_c, and t_c is found by regula falsi within
! dt_can of tair (at once when ra = 0). Leaf temperatures are kept
! within dt_leaf of t_c.
!
! MDK, 24/10/26

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE csigma, ONLY: sbcon
USE canopy_diag_mod, ONLY: h_leaf_pts, le_leaf_pts
USE planet_constants_mod, ONLY: cp, r, repsilon
USE water_constants_mod, ONLY: lc, Water_kg_per_mol
USE jules_vegetation_mod, ONLY: leaf_width, leaf_shelter, leaf_aero_model
USE qsat_mod, ONLY: qsat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! nc: number of leaf classes.
INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts), nc

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Of each class (see leaf_temp_setup; for two-leaf the sunlit and
  ! shaded leaf's lai_sun, sw_sun, ... gl_sun and lai_shd, ... gl_shd):
  ! leaf area, absorbed short-wave, long-wave weights, radiative
  ! conductance scalar, wind-weighted leaf area and stomatal
  ! conductance (m s-1, ground).
  tair(land_pts), q1(land_pts), pstar(land_pts), ra(land_pts),                 &
  lw_down(land_pts), u_wind(land_pts),                                         &
  lai_c(land_pts,nc), sw_c(land_pts,nc), wlw_c(land_pts,nc),                   &
  wlws_c(land_pts,nc), frad_c(land_pts,nc), ff_c(land_pts,nc),                 &
  gl_c(land_pts,nc),                                                           &
  ! CABLE canopy geometry (leaf_aero_model = 1; see leaf_aero_geom).
  z0m_a(land_pts), zref_a(land_pts), rt1ab_a(land_pts), zrd_a(land_pts),       &
  usc_a(land_pts), usuh_a(land_pts),                                           &
  ! Soil evaporation conductance into the canopy air (m s-1, ground
  ! area) and the saturated humidity it evaporates at (kg kg-1).
  g_soil(land_pts), qs_soil(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ! zref/L (leaf_aero_model = 1): the start of the Monin-Obukhov
  ! iteration, and the last iterate on output, so successive passes in
  ! a timestep continue the iteration (0 for the first).
  zeta_a(land_pts),                                                            &
  ! t_cls: leaf temperature of each class.
  t_cls(land_pts,nc), t_c(land_pts), q_c(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ! Leaf boundary layer conductance for water of each class (m/s, class
  ! basis; 0 without leaf area) and the canopy air resistance (s/m) of
  ! the balance found, for l_som_coupled_e. Points not solved keep
  ! their values.
  gbw_c(land_pts,nc), ra_out(land_pts)

! Number of leaf classes (0-nc) whose temperature was held at dt_leaf
! from the canopy air in the balance found.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_clamp(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN OUT), OPTIONAL ::                            &
  ! Each class's leaf temperature as a linear function of its
  ! transpiration E (mol m-2 s-1, ground) at the canopy air of the
  ! balance found, T_leaf = lt_a - lt_b E, i.e. lt_a = t_c + Rn / (cp gh)
  ! (K) and lt_b = lambda / (cp gh) (K per mol m-2 s-1); lt_b = 0
  ! without leaf area (l_leaf_temp_ci). Points not solved keep their
  ! values.
  lt_a_c(land_pts,nc), lt_b_c(land_pts,nc)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  ! CABLE constants: von Karman, gravity, limits on zref/L, minimum rt1.
  vonk = 0.40, grav = 9.81, zetneg = -15.0, zetpos = 1.0, rt_min = 5.0,        &
  ! CABLE constants: minimum wind at zref for u* (m s-1), Pohlhausen
  ! coefficient (single-sided plate), Prandtl number, minimum
  ! boundary-layer conductance of the top leaf (mol m-2 s-1).
  umin_cable = 1.0, apol = 0.70, prandt = 0.71, gbv_min = 0.05,                &
  ! Leaf emissivity (Wang & Leuning, 1998, Table 3).
  eps_leaf = 0.96,                                                             &
  ! Soil emissivity (Wang & Leuning, 1998, Table 3).
  eps_soil = 0.94,                                                             &
  ! Molecular diffusivity of heat (m2 s-1).
  dheat = 21.5e-6,                                                             &
  ! Ratio of boundary-layer conductances for water vapour and heat.
  gbh_2_gbw = 1.075,                                                           &
  ! Molar mass of water (kg mol-1).
  m_h2o = Water_kg_per_mol,                                                    &
  ! Minimum wind speed for forced convection (m s-1).
  u_min = 0.1,                                                                 &
  ! Limits on leaf - canopy air and canopy air - tair (K).
  dt_leaf = 15.0, dt_can = 10.0

! Maximum evaluations of the canopy balance in the t_c search.
INTEGER, PARAMETER :: n_solve = 8
! Tolerance on t_c (K).
REAL(KIND=real_jlslsm), PARAMETER :: tol_tc = 0.01

! Monin-Obukhov iterations for leaf_aero_model = 1 (CABLE niter).
INTEGER, PARAMETER :: n_mo_cable = 4

REAL(KIND=real_jlslsm) ::                                                      &
  ! ncl_last: classes held at dt_leaf in the last canopy_balance call.
  ! tx_v, gbw_last: of each class.
  t_lo(veg_pts), t_hi(veg_pts), r_lo(veg_pts), r_hi(veg_pts), t_try(veg_pts),  &
  t_try1(veg_pts), p_v(veg_pts), qs0_v(veg_pts), qs1_v(veg_pts),               &
  tx_v(nc,veg_pts), qc_v(veg_pts), h_v(veg_pts), le_v(veg_pts),                &
  hb_v(veg_pts), leb_v(veg_pts), res, ra_m(land_pts),                          &
  u_m(land_pts), us, rt1usc, zr, gbw_last(nc), ncl_last, lta_last(nc),       &
  ltb_last(nc)
LOGICAL :: l_done(veg_pts)
INTEGER :: l, m, it, side(veg_pts), it_mo, n_mo

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_UPDATE_N'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( leaf_aero_model == 1 ) THEN
  n_mo = n_mo_cable
ELSE
  n_mo = 1
END IF
ra_m(:) = ra(:)
u_m(:)  = u_wind(:)

DO it_mo = 1,n_mo

IF ( leaf_aero_model == 1 ) THEN
  !---------------------------------------------------------------------------
  ! CABLE (comp_friction_vel, cable_canopy): u* from the wind at zref with
  ! the stability functions, the canopy-top wind u_h = u* / (u*/u_h), and
  ! rt1 = (rt1usa + rt1usb + rt1usc) / u*, rt1usc being the stability-
  ! corrected log profile above the roughness sublayer.
  !---------------------------------------------------------------------------
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( zref_a(l) <= 0.0 ) CYCLE          ! no canopy: keep ra, u_wind
    us = vonk * MAX(u_wind(l), umin_cable)                                     &
         / ( LOG(zref_a(l) / z0m_a(l)) - psim(zeta_a(l))                       &
             + psim(zeta_a(l) * z0m_a(l) / zref_a(l)) )
    us = MIN(MAX(1.0e-6, us), 10.0)
    zr = zrd_a(l)
    rt1usc = usc_a(l) * ( LOG(zref_a(l) / zr) - psis(zeta_a(l))                &
                          + psis(zeta_a(l) * zr / zref_a(l)) ) / vonk
    ra_m(l) = MAX(rt_min, ( rt1ab_a(l) + rt1usc ) / us)
    u_m(l)  = us / usuh_a(l)
  END DO
END IF

!-----------------------------------------------------------------------------
! Find t_c with res(t_c) = 0 by regula falsi (Illinois) in the bracket
! [tair - dt_can, tair + dt_can]; res increases with t_c. Start from the
! bracket ends; points whose root is outside it take the nearer end.
!-----------------------------------------------------------------------------
DO m = 1,veg_pts
  l = veg_index(m)
  p_v(m)  = pstar(l)
  t_lo(m) = tair(l) - dt_can
  t_hi(m) = tair(l) + dt_can
  side(m) = 0
END DO

DO it = 1,2
  DO m = 1,veg_pts
    IF ( it == 1 ) THEN
      t_try(m) = t_lo(m)
    ELSE
      t_try(m) = t_hi(m)
    END IF
    t_try1(m) = t_try(m) + 0.1
  END DO
  CALL qsat(qs0_v, t_try, p_v, veg_pts)
  CALL qsat(qs1_v, t_try1, p_v, veg_pts)
  DO m = 1,veg_pts
    l = veg_index(m)
    CALL canopy_balance( l, t_try(m), qs0_v(m), qs1_v(m), res, tx_v(:,m),      &
                         qc_v(m), hb_v(m), leb_v(m) )
    IF ( it == 1 ) THEN
      r_lo(m) = res
      ! Root below the bracket: keep the lower end.
      l_done(m) = ( res >= 0.0 )
      IF ( l_done(m) ) CALL set_out( l, m, t_try(m), qc_v(m) )
    ELSE
      r_hi(m) = res
      IF ( .NOT. l_done(m) .AND. res <= 0.0 ) THEN
        ! Root above the bracket: keep the upper end.
        l_done(m) = .TRUE.
        CALL set_out( l, m, t_try(m), qc_v(m) )
      END IF
    END IF
  END DO
END DO

DO it = 1,n_solve
  IF ( ALL(l_done) ) EXIT
  DO m = 1,veg_pts
    IF ( l_done(m) ) THEN
      t_try(m) = t_lo(m)
    ELSE
      t_try(m) = t_hi(m) - r_hi(m) * ( t_hi(m) - t_lo(m) )                     &
                 / ( r_hi(m) - r_lo(m) )
    END IF
    t_try1(m) = t_try(m) + 0.1
  END DO
  CALL qsat(qs0_v, t_try, p_v, veg_pts)
  CALL qsat(qs1_v, t_try1, p_v, veg_pts)
  DO m = 1,veg_pts
    IF ( l_done(m) ) CYCLE
    l = veg_index(m)
    CALL canopy_balance( l, t_try(m), qs0_v(m), qs1_v(m), res, tx_v(:,m),      &
                         qc_v(m), hb_v(m), leb_v(m) )
    IF ( ABS(res) < tol_tc .OR. it == n_solve ) THEN
      l_done(m) = .TRUE.
      CALL set_out( l, m, t_try(m), qc_v(m) )
    ELSE IF ( res > 0.0 ) THEN
      t_hi(m) = t_try(m)
      r_hi(m) = res
      IF ( side(m) == 1 ) r_lo(m) = 0.5 * r_lo(m)
      side(m) = 1
    ELSE
      t_lo(m) = t_try(m)
      r_lo(m) = res
      IF ( side(m) == -1 ) r_hi(m) = 0.5 * r_hi(m)
      side(m) = -1
    END IF
  END DO
END DO

! Stability for the next iteration (CABLE): zref/L from the canopy sensible
! heat and the buoyancy of the latent heat,
!
!   zref/L = -k g zref (H + 0.07 LE) / (rho cp T u*^3).
!
! Only the leaves' fluxes: the soil's sensible heat is not known here (nor
! is its latent heat added). Also after the last iteration, for the next
! pass to start from.
IF ( leaf_aero_model == 1 ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( zref_a(l) <= 0.0 ) CYCLE
    us = u_m(l) * usuh_a(l)
    zeta_a(l) = -vonk * grav * zref_a(l) * ( h_v(m) + 0.07 * le_v(m) )         &
              / ( pstar(l) / ( r * tair(l) ) * cp * tair(l) * us**3 )
    zeta_a(l) = MIN(MAX(zeta_a(l), zetneg), zetpos)
  END DO
END IF

END DO   ! it_mo

DO m = 1,veg_pts
  l = veg_index(m)
  ra_out(l) = ra_m(l)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

SUBROUTINE set_out( l, m, tc, qc )
! The outputs at point l (the class temperatures from tx_v).
INTEGER, INTENT(IN) :: l, m
REAL(KIND=real_jlslsm), INTENT(IN) :: tc, qc
t_c(l)   = tc
t_cls(l,:) = tx_v(:,m)
q_c(l)   = qc
! (gbw_last is from the canopy_balance call just made for this point.)
gbw_c(l,:) = gbw_last(:)
h_v(m)   = hb_v(m)
le_v(m)  = leb_v(m)

! For the canopy diagnostics (canopy_diag_mod).
IF ( ALLOCATED(h_leaf_pts) ) THEN
  h_leaf_pts(l)  = hb_v(m)
  le_leaf_pts(l) = leb_v(m)
END IF

! ncl_last is from that call too, and so are lta_last and ltb_last
! (l_leaf_temp_ci).
n_clamp(l) = ncl_last
IF ( PRESENT(lt_a_c) ) lt_a_c(l,:) = lta_last(:)
IF ( PRESENT(lt_b_c) ) lt_b_c(l,:) = ltb_last(:)
END SUBROUTINE set_out

ELEMENTAL FUNCTION psim( zeta_in ) RESULT( r_psi )
! CABLE psim (cable_canopy): stability function for momentum.
REAL(KIND=real_jlslsm), INTENT(IN) :: zeta_in
REAL(KIND=real_jlslsm) :: r_psi, x
REAL(KIND=real_jlslsm), PARAMETER :: gu = 16.0, a = 1.0, b = 0.667, xc = 5.0,  &
                                     d = 0.35, pi_c = 3.14159265
IF ( zeta_in >= 0.0 ) THEN
  r_psi = -a * zeta_in - b * (zeta_in - xc / d) * EXP(-d * zeta_in) - b * xc / d
ELSE
  x     = (1.0 + gu * ABS(zeta_in))**0.25
  r_psi = LOG((1.0 + x * x) * (1.0 + x)**2 / 8.0) - 2.0 * ATAN(x) + pi_c * 0.5
END IF
END FUNCTION psim

ELEMENTAL FUNCTION psis( zeta_in ) RESULT( r_psi )
! CABLE psis (cable_canopy): stability function for scalars.
REAL(KIND=real_jlslsm), INTENT(IN) :: zeta_in
REAL(KIND=real_jlslsm) :: r_psi, y
REAL(KIND=real_jlslsm), PARAMETER :: gu = 16.0, a = 1.0, b = 0.667, c = 5.0,   &
                                     d = 0.35
IF ( zeta_in >= 0.0 ) THEN
  r_psi = -(1.0 + 2.0 / 3.0 * a * zeta_in)**1.5                                &
          - b * (zeta_in - c / d) * EXP(-d * zeta_in) - b * c / d + 1.0
ELSE
  y     = (1.0 + gu * ABS(zeta_in))**0.5
  r_psi = 2.0 * LOG((1.0 + y) * 0.5)
END IF
END FUNCTION psis

SUBROUTINE canopy_balance( l, tc, qs0, qs1, res, t_new, qc_new, h_out, le_out )

! Leaf temperatures (of the nc classes) and canopy humidity for canopy air
! temperature tc, and the residual of the canopy air sensible heat balance,
!
!   res = tc - tair - ra (sum of the classes' H) / (rho cp).
!

INTEGER, INTENT(IN) :: l
REAL(KIND=real_jlslsm), INTENT(IN) :: tc, qs0, qs1
REAL(KIND=real_jlslsm), INTENT(OUT) :: res, t_new(nc), qc_new, h_out, le_out

REAL(KIND=real_jlslsm) ::                                                      &
  cp_mol, cmolar, lhv, slope, gamma_p, sig_t4, eps_air, grn1, gbhw, rho, c_q,  &
  rn(nc), gh(nc), gbh(nc), a(nc), b(nc), le, dt, h_tot, lai_x(nc), gl_x(nc),   &
  t_old(nc), frad(nc), ff(nc), gbhf, gsw, gbw, gw, den, grashof, ra_x, visc,   &
  sum_a, sum_b
INTEGER :: ic, it_q
LOGICAL :: l_act(nc), l_neg

cp_mol  = cp * rmol / r
ra_x    = ra_m(l)
cmolar  = pstar(l) / ( rmol * tc )
lhv     = ( lc - 2.365e3 * (tc - zerodegc) ) * m_h2o
slope   = ( qs1 - qs0 ) / 0.1 * pstar(l) / repsilon
gamma_p = cp_mol * pstar(l) / lhv
sig_t4  = sbcon * tc**4
eps_air = lw_down(l) / sig_t4
grn1    = 4.0 * eps_leaf * sbcon * tc**3 / cp_mol
IF ( leaf_aero_model == 1 ) THEN
  ! CABLE's gbvtop (Raupach et al., 1997, eq. 3.12), with the kinematic
  ! viscosity of air (cable_air) and the floor CABLE applies.
  visc  = 1.0e-5 * MAX(1.0, 1.35 + 0.0092 * (tc - zerodegc))
  gbhw  = cmolar * apol * visc / prandt / leaf_width                           &
          * SQRT(MAX(u_m(l), u_min) * leaf_width / visc)                       &
          * prandt**(1.0 / 3.0) / leaf_shelter
  gbhw  = MAX(gbv_min, gbhw)
ELSE
  gbhw  = 0.003 * SQRT(MAX(u_m(l), u_min) / leaf_width) * cmolar               &
          / leaf_shelter
END IF
rho     = pstar(l) / ( r * tair(l) )

lai_x = lai_c(l,:)
gl_x  = gl_c(l,:)
t_old = t_cls(l,:)
frad  = frad_c(l,:)
ff    = ff_c(l,:)
DO ic = 1,nc
  rn(ic) = sw_c(l,ic) + eps_leaf * (lw_down(l) - sig_t4) * wlw_c(l,ic)         &
           - sig_t4 * (1.0 - eps_soil) * (eps_leaf - eps_air) * wlws_c(l,ic)
END DO

! LE of each class is a + b (qs(tc) - q_c) (W m-2).
DO ic = 1,nc
  a(ic)   = 0.0
  b(ic)   = 0.0
  gh(ic)  = 0.0
  gbh(ic) = 0.0
  gbw_last(ic) = 0.0
  IF ( lai_x(ic) > EPSILON(0.0) ) THEN
    ! Free convection from the last leaf - canopy air difference.
    IF ( ABS(t_old(ic) - tc) > 1.0e-6 ) THEN
      grashof = MAX(1.0e-6, 1.6e8 * ABS(t_old(ic) - tc) * leaf_width**3)
      gbhf    = 0.5 * dheat * grashof**0.25 / leaf_width * cmolar
    ELSE
      gbhf    = 0.0
    END IF
    gbh(ic) = gbhw * ff(ic) + gbhf * lai_x(ic)
    gh(ic)  = 2.0 * gbh(ic) + grn1 * frad(ic)
    gsw     = MAX(gl_x(ic), 0.0) * cmolar
    gbw     = gbh_2_gbw * gbh(ic)
    gbw_last(ic) = gbw / cmolar
    IF ( gsw > 0.0 .AND. gbw > 0.0 ) THEN
      gw    = gbw * gsw / ( gbw + gsw )
      den   = slope + gamma_p * gh(ic) / gw
      a(ic) = slope * rn(ic) / den
      b(ic) = gh(ic) * cp_mol * pstar(l) / repsilon / den
    END IF
  END IF
END DO

! Canopy humidity:
!
!   q_c = q1 + c_q sum(LE),   c_q = ra / (rho L),
!
! with L the latent heat at tc that the leaves' LE uses (lhv / m_h2o,
! J kg-1; the constant lc is 2-3 % off at 25-35 C). The soil's vapour flux
! is rho g_soil (qs_soil - q_c), so ra multiplies g_soil.
c_q    = ra_x * m_h2o / ( rho * lhv )

! Only the classes whose LE is positive at the solved q_c take part, as the
! leaf heat balance below clips each class's LE at 0: a class that would
! condense (at night) adds no vapour to, and removes none from, the canopy
! air. Classes are dropped and q_c re-solved until none of those left is
! negative (at most nc passes; with none negative, a single solve).
l_act(:) = .TRUE.
DO it_q = 1,nc
  sum_a  = 0.0
  sum_b  = 0.0
  DO ic = 1,nc
    IF ( l_act(ic) ) THEN
      sum_a = sum_a + a(ic)
      sum_b = sum_b + b(ic)
    END IF
  END DO
  qc_new = ( q1(l) + c_q * ( sum_a + sum_b * qs0 )                             &
             + ra_x * g_soil(l) * qs_soil(l) )                                 &
           / ( 1.0 + c_q * sum_b + ra_x * g_soil(l) )
  l_neg = .FALSE.
  DO ic = 1,nc
    IF ( l_act(ic) .AND. a(ic) + b(ic) * (qs0 - qc_new) < 0.0 ) THEN
      l_act(ic) = .FALSE.
      l_neg = .TRUE.
    END IF
  END DO
  IF ( .NOT. l_neg ) EXIT
END DO

h_tot  = 0.0
le_out = 0.0
ncl_last = 0.0
DO ic = 1,nc
  ! The leaf temperature as a function of the class's transpiration at
  ! this tc (l_leaf_temp_ci),
  !
  !   T = tc + (rn - lhv E) / (cp_mol gh).
  !
  lta_last(ic) = tc
  ltb_last(ic) = 0.0
  IF ( lai_x(ic) > EPSILON(0.0) .AND. gh(ic) > 0.0 ) THEN
    lta_last(ic) = tc + rn(ic) / ( cp_mol * gh(ic) )
    ltb_last(ic) = lhv / ( cp_mol * gh(ic) )
  END IF
  IF ( lai_x(ic) > EPSILON(0.0) ) THEN
    le     = MAX(0.0, a(ic) + b(ic) * (qs0 - qc_new))
    le_out = le_out + le
    dt     = ( rn(ic) - le ) / ( cp_mol * gh(ic) )
    ! Count the classes held at the limit (t_leaf_clamp).
    IF ( ABS(dt) > dt_leaf ) ncl_last = ncl_last + 1.0
    dt     = MIN(MAX(dt, -dt_leaf), dt_leaf)
    t_new(ic) = tc + dt
    h_tot  = h_tot + 2.0 * gbh(ic) * cp_mol * dt
  ELSE
    t_new(ic) = tc
  END IF
END DO

res     = tc - tair(l) - ra_x * h_tot / ( rho * cp )
h_out   = h_tot

END SUBROUTINE canopy_balance

END SUBROUTINE leaf_temp_update_n

END MODULE leaf_energy_balance_mod
