! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT****************************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE leaf_energy_balance_mod

!-----------------------------------------------------------------------------
! Leaf energy balance (l_leaf_temp): the temperature of the sunlit and
! shaded leaves from the light they absorb, the heat they lose and the
! water they transpire. The leaves exchange heat and water with the air
! of the first atmospheric level, as in MAESPA. sf_stom calls these
! routines:
!
!   - leaf_temp_setup: the radiation and wind terms of each set of leaves
!     that stay the same on every pass of the stomatal optimisation;
!   - leaf_temp_update and leaf_temp_update_n: the leaf temperatures for
!     the stomatal conductances of the last pass.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEAF_ENERGY_BALANCE_MOD'

PRIVATE
PUBLIC :: leaf_temp_setup, leaf_temp_update, leaf_temp_update_n

CONTAINS

! Work out the radiation and wind terms of the leaf energy balance
! (l_leaf_temp) once per timestep, before the stomatal optimisation, since
! they stay the same on every pass of it. They are for the sunlit and the
! shaded leaves (two-leaf canopy, and the multilayer canopy taken as one
! set of sunlit and one of shaded leaves). Each is summed
! over the canopy layers using the sunlit fraction of each layer, as the
! absorbed PAR is. xi is the leaf area above a point in the canopy, lai
! the canopy's leaf area and kd the extinction coefficient for diffuse light
! (1 for horizontal leaves, 0.78 otherwise).
!
!   - sw: absorbed short-wave, PAR plus near-infrared (W m-2 ground). The
!     incident near-infrared is taken as equal to the incident PAR, as
!     standalone JULES sets PAR to half the incoming short-wave.
!   - wlw, wlws: weights for the long-wave radiation the leaves absorb
!     (Wang and Leuning, 1998, eqs B18 and B19): from the sky, kd
!     exp(-kd xi), and from the soil, kd exp(-kd (lai - xi)), each
!     integrated over the layer.
!   - frad: the sum of the two, which scales the leaves' long-wave loss
!     (about 2 kd times the leaf area for a thin canopy).
!   - ff: leaf area weighted by the relative wind speed, exp(-wind_ext
!     xi / 2), for heat loss to the wind (Wang and Leuning, 1998; CABLE).
!     It is the leaf area when wind_ext = 0.
!
! With the optional arguments (multilayer canopy), the same terms are
! also given for the sunlit and shaded leaves of every layer: the sunlit
! leaves of layer n in position 2n-1, the shaded
! in 2n. Their PAR comes from each layer's own absorbed PAR
! (fapar_*_vis), and lai_c gives their leaf area. The layer terms add up
! to the canopy ones.
!
! MDK, 24/10/26
SUBROUTINE leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,    &
                            dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,    &
                            wind_ext,                                          &
                            apar_sun, apar_shd,                                &
                            sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,        &
                            wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,      &
                            fapar_sun_vis, fapar_shd_vis, lai_c, sw_c, wlw_c,  &
                            wlws_c, frad_c, ff_c )

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

! Work out the sunlit and shaded leaf temperatures of the two-leaf
! canopy from the leaf energy balance (l_leaf_temp) for the stomatal
! conductances of the last pass. It packs the sunlit and shaded leaves
! into arrays of two and calls leaf_temp_update_n.
!
! MDK, 04/11/26
SUBROUTINE leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             lw_down, u_wind, lai_sun, lai_shd,                &
                             sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,       &
                             wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,     &
                             gl_sun, gl_shd, t_sun, t_shd, gbw_sun, gbw_shd,   &
                             n_clamp, lt_a, lt_b )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tair(land_pts), q1(land_pts), pstar(land_pts),                               &
  lw_down(land_pts), u_wind(land_pts), lai_sun(land_pts), lai_shd(land_pts),   &
  sw_sun(land_pts), sw_shd(land_pts), wlw_sun(land_pts), wlw_shd(land_pts),    &
  wlws_sun(land_pts), wlws_shd(land_pts), frad_sun(land_pts),                  &
  frad_shd(land_pts), ff_sun(land_pts), ff_shd(land_pts), gl_sun(land_pts),    &
  gl_shd(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  t_sun(land_pts), t_shd(land_pts), gbw_sun(land_pts), gbw_shd(land_pts),      &
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
                         lw_down, u_wind, lai_x, sw_x, wlw_x, wlws_x,          &
                         frad_x, ff_x, gl_x, t_x, gbw_x, n_clamp, lt_a, lt_b )

t_sun(:)   = t_x(:,1)
t_shd(:)   = t_x(:,2)
gbw_sun(:) = gbw_x(:,1)
gbw_shd(:) = gbw_x(:,2)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_temp_update

!#############################################################################
!#############################################################################

! Work out the leaf temperatures of nc sets of leaves by solving the leaf
! energy balance (l_leaf_temp) for the stomatal conductances of the last
! pass. The sets are the sunlit and shaded leaves (nc = 2) of the
! two-leaf canopy, or the sunlit and shaded leaves of every layer of the
! multilayer canopy (nc = 2 x ilayers, as in MAESPA).
!
! Each set of leaves exchanges heat and water through its leaf boundary
! layer with the air of the first atmospheric level (temperature tair,
! humidity q1), as one big leaf, as in MAESPA and MAESTRA. Its latent
! heat follows Penman-Monteith with the net radiation it would have at
! the air temperature (Leuning et al., 1995, App. E; Wang and Leuning,
! 1998), as in De Kauwe's two-leaf model (two_leaf_at_WTC):
!
!   LE = (s Rn + D gh cp) / (s + gamma gh / gw)
!   T  = tair + (Rn - LE) / (cp gh)
!   gh = 2 gbH + grn
!
! where LE is the latent heat (W m-2), Rn the net radiation, s the slope
! of saturation humidity with temperature, D the humidity deficit of the
! air, gamma the psychrometric constant, cp the heat capacity of air, gw
! the conductance for water (stomata and boundary layer in series) and T
! the leaf temperature. gh is the conductance for heat: twice the
! boundary layer conductance gbH (both sides of the leaf) plus the
! radiative conductance grn. gbH has a part from the wind (per unit leaf
! area at the level-1 wind, divided by leaf_shelter, times the
! wind-weighted leaf area ff) and a part from rising warm air (times the
! leaf area). Each leaf temperature is kept within dt_leaf of tair.
!
! MDK, 04/11/26
SUBROUTINE leaf_temp_update_n( land_pts, veg_pts, veg_index, nc, tair, q1,     &
                               pstar, lw_down, u_wind, lai_c, sw_c, wlw_c,     &
                               wlws_c, frad_c, ff_c, gl_c, t_cls, gbw_c,       &
                               n_clamp, lt_a_c, lt_b_c )

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE csigma, ONLY: sbcon
USE canopy_diag_mod, ONLY: h_leaf_pts, le_leaf_pts
USE planet_constants_mod, ONLY: cp, r, repsilon
USE water_constants_mod, ONLY: lc, Water_kg_per_mol
USE jules_vegetation_mod, ONLY: leaf_width, leaf_shelter
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
  tair(land_pts), q1(land_pts), pstar(land_pts),                               &
  lw_down(land_pts), u_wind(land_pts),                                         &
  lai_c(land_pts,nc), sw_c(land_pts,nc), wlw_c(land_pts,nc),                   &
  wlws_c(land_pts,nc), frad_c(land_pts,nc), ff_c(land_pts,nc),                 &
  gl_c(land_pts,nc)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ! t_cls: leaf temperature of each class.
  t_cls(land_pts,nc),                                                          &
  ! Leaf boundary layer conductance for water of each class (m/s, class
  ! basis; 0 without leaf area).
  gbw_c(land_pts,nc)

! Number of leaf classes (0-nc) whose temperature was held at dt_leaf
! from the air.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_clamp(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN OUT), OPTIONAL ::                            &
  ! Each class's leaf temperature as a linear function of its
  ! transpiration E (mol m-2 s-1, ground), T_leaf = lt_a - lt_b E, i.e.
  ! lt_a = tair + Rn / (cp gh) (K) and lt_b = lambda / (cp gh) (K per
  ! mol m-2 s-1); lt_b = 0 without leaf area (l_leaf_temp_ci).
  lt_a_c(land_pts,nc), lt_b_c(land_pts,nc)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
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
  ! Limit on leaf - air temperature (K).
  dt_leaf = 15.0

REAL(KIND=real_jlslsm) ::                                                      &
  ! Air temperature and 0.1 K above it, and the saturated humidity at
  ! each (kg kg-1), for the slope s.
  t0_v(veg_pts), t1_v(veg_pts), p_v(veg_pts), qs0_v(veg_pts), qs1_v(veg_pts)
INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_UPDATE_N'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO m = 1,veg_pts
  l = veg_index(m)
  p_v(m)  = pstar(l)
  t0_v(m) = tair(l)
  t1_v(m) = t0_v(m) + 0.1
END DO
CALL qsat(qs0_v, t0_v, p_v, veg_pts)
CALL qsat(qs1_v, t1_v, p_v, veg_pts)

DO m = 1,veg_pts
  l = veg_index(m)
  CALL leaf_balance( l, qs0_v(m), qs1_v(m) )
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! Solve the energy balance of every class at point l and store the leaf
! temperatures, boundary layer conductances, lt_a and lt_b, the number
! of classes held at dt_leaf and the leaves' total sensible and latent
! heat (canopy_diag_mod). qs0 and qs1 are the saturation humidity
! (kg kg-1) at tair and tair + 0.1 K, which give the slope s.
!
! Sets of leaves whose latent heat would be negative (dew, at night) have
! their latent heat set to 0.
SUBROUTINE leaf_balance( l, qs0, qs1 )

INTEGER, INTENT(IN) :: l
REAL(KIND=real_jlslsm), INTENT(IN) :: qs0, qs1

REAL(KIND=real_jlslsm) ::                                                      &
  tc, cp_mol, cmolar, lhv, slope, gamma_p, sig_t4, eps_air, grn1, gbhw,        &
  rn, gh, gbh, a, b, le, dt, h_tot, le_tot, gbhf, gsw, gbw, gw, den, grashof,  &
  n_cl
INTEGER :: ic

tc      = tair(l)
cp_mol  = cp * rmol / r
cmolar  = pstar(l) / ( rmol * tc )
lhv     = ( lc - 2.365e3 * (tc - zerodegc) ) * m_h2o
slope   = ( qs1 - qs0 ) / 0.1 * pstar(l) / repsilon
gamma_p = cp_mol * pstar(l) / lhv
sig_t4  = sbcon * tc**4
eps_air = lw_down(l) / sig_t4
grn1    = 4.0 * eps_leaf * sbcon * tc**3 / cp_mol
gbhw    = 0.003 * SQRT(MAX(u_wind(l), u_min) / leaf_width) * cmolar            &
          / leaf_shelter

h_tot  = 0.0
le_tot = 0.0
n_cl   = 0.0
DO ic = 1,nc
  rn = sw_c(l,ic) + eps_leaf * (lw_down(l) - sig_t4) * wlw_c(l,ic)             &
       - sig_t4 * (1.0 - eps_soil) * (eps_leaf - eps_air) * wlws_c(l,ic)

  ! LE of the class is a + b (qs(tair) - q1) (W m-2).
  a   = 0.0
  b   = 0.0
  gh  = 0.0
  gbh = 0.0
  gbw_c(l,ic) = 0.0
  IF ( lai_c(l,ic) > EPSILON(0.0) ) THEN
    ! Free convection from the last leaf - air difference.
    IF ( ABS(t_cls(l,ic) - tc) > 1.0e-6 ) THEN
      grashof = MAX(1.0e-6, 1.6e8 * ABS(t_cls(l,ic) - tc) * leaf_width**3)
      gbhf    = 0.5 * dheat * grashof**0.25 / leaf_width * cmolar
    ELSE
      gbhf    = 0.0
    END IF
    gbh = gbhw * ff_c(l,ic) + gbhf * lai_c(l,ic)
    gh  = 2.0 * gbh + grn1 * frad_c(l,ic)
    gsw = MAX(gl_c(l,ic), 0.0) * cmolar
    gbw = gbh_2_gbw * gbh
    gbw_c(l,ic) = gbw / cmolar
    IF ( gsw > 0.0 .AND. gbw > 0.0 ) THEN
      gw  = gbw * gsw / ( gbw + gsw )
      den = slope + gamma_p * gh / gw
      a   = slope * rn / den
      b   = gh * cp_mol * pstar(l) / repsilon / den
    END IF
  END IF

  ! The leaf temperature as a function of the class's transpiration
  ! (l_leaf_temp_ci),
  !
  !   T = tair + (rn - lhv E) / (cp_mol gh).
  !
  IF ( PRESENT(lt_a_c) ) THEN
    lt_a_c(l,ic) = tc
    lt_b_c(l,ic) = 0.0
    IF ( lai_c(l,ic) > EPSILON(0.0) .AND. gh > 0.0 ) THEN
      lt_a_c(l,ic) = tc + rn / ( cp_mol * gh )
      lt_b_c(l,ic) = lhv / ( cp_mol * gh )
    END IF
  END IF

  IF ( lai_c(l,ic) > EPSILON(0.0) ) THEN
    le     = MAX(0.0, a + b * (qs0 - q1(l)))
    le_tot = le_tot + le
    dt     = ( rn - le ) / ( cp_mol * gh )
    ! Count the classes held at the limit (t_leaf_clamp).
    IF ( ABS(dt) > dt_leaf ) n_cl = n_cl + 1.0
    dt     = MIN(MAX(dt, -dt_leaf), dt_leaf)
    t_cls(l,ic) = tc + dt
    h_tot  = h_tot + 2.0 * gbh * cp_mol * dt
  ELSE
    t_cls(l,ic) = tc
  END IF
END DO

n_clamp(l) = n_cl

! For the canopy diagnostics (canopy_diag_mod).
IF ( ALLOCATED(h_leaf_pts) ) THEN
  h_leaf_pts(l)  = h_tot
  le_leaf_pts(l) = le_tot
END IF

END SUBROUTINE leaf_balance

END SUBROUTINE leaf_temp_update_n

END MODULE leaf_energy_balance_mod
