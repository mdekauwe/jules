!******************************COPYRIGHT********************************************
! (c) University of Exeter 2022
! All rights reserved.
!
! This routine has been licensed to the other JULES partners for use and
! distribution under the JULES collaboration agreement, subject to the terms and
! conditions set out therein.
!
! [Met Office Ref SC0237]
!******************************COPYRIGHT********************************************

MODULE sugar_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SUGAR_MOD'

PRIVATE
PUBLIC sugar, sugar_sink_demand

CONTAINS
! *********************************************************************
! Routines to calculate the utilisation of non-structural carbohydrate
! by respiration and growth for each PFT from the SUGAR model
!
! References:
!   Jones et al., 2020, Biogeosciences, 17, 3589-3612,
!        https://doi.org/10.5194/bg-17-3589-2020
! *********************************************************************

SUBROUTINE sugar(ft, leafc, leafc_bal, woodc, rootc, f_nsc, resp_l,            &
                 resp_w, resp_r, resp_p_m, resp_p_g,                           &
                 resp_p, growth, tstar, gpp, f_turgor, flush)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: sugar_model, sugar_mm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                         &
 ft
    ! IN Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tstar                                                                         &
    ! IN Surface temperature (K)
,gpp                                                                           &
    ! IN Gross Primary Productivity (kg C/m2/s)
,leafc                                                                         &
    ! IN (Phenological) Leaf carbon (kg C/m2)
,leafc_bal                                                                     &
    ! IN (Balanced) Leaf carbon (kg C/m2)
,woodc                                                                         &
    ! IN Wood carbon (kg C/m2)
,rootc                                                                         &
    ! IN Root carbon (kg C/m2).
,f_turgor                                                                      &
    ! IN Turgor limit on structural growth (0-1; 1 = none, l_sugar_turgor=F)
,flush
    ! IN Leaf-flush construction, new structural leaf carbon (kg C/m2/s;
    !    0 = none, l_sugar_leaf_flush=F). sugar_model = 2 only.

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
f_nsc
    ! INOUT NSC mass fraction (kg kg-1)

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 resp_w                                                                        &
    ! OUT Wood maintenance respiration rate (kg C/m2/s).
,resp_r                                                                        &
    ! OUT Root maintenance respiration rate (kg C/m2/s)
,resp_l                                                                        &
    ! OUT Leaf maintenance respiration rate (kg C/m2/s)
,resp_p_m                                                                      &
    ! OUT Plant maintenance respiration (kg C/m2/s)
,resp_p_g                                                                      &
    ! OUT Plant growth respiration (kg C/m2/s)
,resp_p                                                                        &
    ! OUT Plant respiration (kg C/m2/s)
,growth
    ! OUT Plant growth rate (kg C/m2/s)

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 df_nsc
    ! WORK Increment in the NSC mass fraction

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SUGAR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Jones et al. (2020) form (sugar_model = 2).
!------------------------------------------------------------------------------
IF ( sugar_model == sugar_mm ) THEN
  CALL sugar_mm_step(ft, leafc, leafc_bal, woodc, rootc, f_nsc, resp_l,        &
                     resp_w, resp_r, resp_p_m, resp_p_g, resp_p, growth,       &
                     tstar, gpp, f_turgor, flush)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!------------------------------------------------------------------------------
! Calculate explicit fluxes using previous timestep f_nsc value
!------------------------------------------------------------------------------
CALL calculate_sug_fluxes(ft, leafc, woodc, rootc, f_nsc, resp_l,              &
                          resp_w, resp_r, resp_p_m, resp_p_g,                  &
                          resp_p, growth, tstar, f_turgor)

!------------------------------------------------------------------------------
! Update NSC mass fraction (f_nsc)
!------------------------------------------------------------------------------
CALL update_nsc_pool(ft, gpp, resp_p, growth, f_nsc, leafc, leafc_bal,         &
                     woodc, rootc, df_nsc, tstar, f_turgor)

!------------------------------------------------------------------------------
! Perform flux corrections to ensure they match the f_nsc increment
! This is only needed if a forward timestep weighting is used
! (which it is by default) but will still be performed if the explicit
! scheme (forw_gamma=0) is used (the corrections will just be zero).
!------------------------------------------------------------------------------
CALL sugar_flux_correction(ft, leafc_bal, leafc, woodc, rootc, df_nsc, resp_l, &
                           resp_w, resp_r, resp_p_m, resp_p_g,                 &
                           resp_p, growth, tstar, gpp, f_turgor)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sugar

SUBROUTINE calculate_sug_fluxes(ft, leafc, woodc, rootc, f_nsc, resp_l,        &
                                resp_w, resp_r, resp_p_m, resp_p_g,            &
                                resp_p, growth, tstar, f_turgor)

USE pftparm, ONLY:                                                             &
! imported parameters
    sug_grec, sug_g0, sug_yg, q10_leaf

USE conversions_mod, ONLY: zerodegc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE
!------------------------------------------------------------------------------
! Arguments with INTENT(IN).
!------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
 ft
    ! IN Plant functional type

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tstar                                                                         &
    ! IN Surface temperature (K)
,f_nsc                                                                         &
    ! INOUT NSC mass fraction (kg kg-1)
,leafc                                                                         &
    ! IN (Phenological) Leaf carbon (kg C/m2)
,woodc                                                                         &
    ! IN Wood carbon (kg C/m2)
,rootc                                                                         &
    ! IN Root carbon (kg C/m2).
,f_turgor
    ! IN Turgor limit on structural growth (0-1)

!------------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 resp_w                                                                        &
    ! OUT Wood maintenance respiration rate (kg C/m2/s).
,resp_r                                                                        &
    ! OUT Root maintenance respiration rate (kg C/m2/s)
,resp_l                                                                        &
    ! OUT Leaf maintenance respiration rate (kg C/m2/s)
,resp_p_m                                                                      &
    ! OUT Plant maintenance respiration (kg C/m2/s)
,resp_p_g                                                                      &
    ! OUT Plant growth respiration (kg C/m2/s)
,resp_p                                                                        &
    ! OUT Plant respiration (kg C/m2/s)
,growth
    ! OUT Plant growth rate (kg C/m2/s)

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 r0                                                                            &
    ! WORK Specific rate of NSC utilisation by respiration
,g0                                                                            &
    ! WORK Specific rate of NSC utilisation by structural growth
,grec0                                                                         &
    ! WORK Specific rate of structural C decomposition
,f_temp                                                                        &
    ! WORK Temperature function used for respiration and growth rates
,growth_l                                                                      &
    ! WORK Leaf growth rate (kg C/m2/s)
,growth_w                                                                      &
    ! WORK Wood growth rate (kg C/m2/s)
,growth_r                                                                      &
    ! WORK Root growth rate (kg C/m2/s)
,grecl                                                                         &
    ! WORK Decay rate of structural leaf carbon into NSC (leaf maintenance)
,grecw                                                                         &
    ! WORK Decay rate of structural stem carbon into NSC (stem maintenance)
,grecr
    ! WORK Decay rate of structural root carbon into NSC (root maintenance)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALCULATE_SUG_FLUXES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Calculate temperature function
!------------------------------------------------------------------------------
f_temp = q10_leaf(ft)**( 0.1 * (tstar - zerodegc - 25.0) )

!------------------------------------------------------------------------------
! Calculate temperature dependent coefficients
!------------------------------------------------------------------------------
! The turgor limit acts on structural production (g0) only: respiration
! (r0) keeps its NSC and temperature dependence.
r0    = sug_g0(ft)   * f_temp * ( 1.0 - sug_yg(ft) ) / sug_yg(ft)
g0    = sug_g0(ft)   * f_temp * f_turgor
grec0 = sug_grec(ft) * f_temp

!------------------------------------------------------------------------------
! Calculate fluxes
!------------------------------------------------------------------------------
resp_l   = r0 * f_nsc * leafc
resp_w   = r0 * f_nsc * woodc
resp_r   = r0 * f_nsc * rootc

grecl    = grec0 * leafc * (1.0 - f_nsc)
grecw    = grec0 * woodc * (1.0 - f_nsc)
grecr    = grec0 * rootc * (1.0 - f_nsc)

growth_l = g0 * f_nsc * leafc - grecl
growth_w = g0 * f_nsc * woodc - grecw
growth_r = g0 * f_nsc * rootc - grecr

growth   = growth_l + growth_w + growth_r
resp_p   = resp_l + resp_w + resp_r
resp_p_m = (1.0 - sug_yg(ft)) * (grecl + grecw + grecr) / sug_yg(ft)
resp_p_g = resp_p - resp_p_m

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE calculate_sug_fluxes

SUBROUTINE update_NSC_pool(ft, gpp, resp_p, growth, f_nsc, leafc,              &
                           leafc_bal, woodc, rootc, df_nsc_corr, tstar,        &
                           f_turgor)

USE timestep_mod, ONLY: timestep

USE pftparm, ONLY:                                                             &
! imported parameters
    sug_grec, sug_g0, sug_yg, q10_leaf

USE conversions_mod, ONLY: zerodegc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
 ft
    ! IN Plant functional type

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 gpp                                                                           &
    ! IN Gross primary productivity (kg C/m2/s)
,tstar                                                                         &
    ! IN Surface temperature (K)
,leafc                                                                         &
    ! IN Phenological leaf carbon (kg/m2)
,leafc_bal                                                                     &
    ! IN Balanced leaf carbon (kg/m2)
,woodc                                                                         &
    ! IN Wood carbon (kg/m2)
,rootc                                                                         &
    ! IN Root carbon (kg/m2)
,resp_p                                                                        &
    ! IN Total plant respiration rate (kg C/m2/s)
,growth                                                                        &
    ! IN Total plant growth rate (kg C/m2/s)
,f_turgor
    ! IN Turgor limit on structural growth (0-1)

!-----------------------------------------------------------------------------
! Arguments with INTENT(INOUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
f_nsc
    ! INOUT NSC mass fraction (kg kg-1)

!------------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!------------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 df_nsc_corr
    ! OUT Corrected increment to NSC mass fraction

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 r0                                                                            &
    ! WORK Specific rate of NSC utilisation by respiration
,g0                                                                            &
    ! WORK Specific rate of NSC utilisation by structural growth
,grec0                                                                         &
    ! WORK Specific rate of structural C decomposition
,f_temp                                                                        &
    ! WORK Temperature function used for respiration and growth rates
,c_veg                                                                         &
    ! WORK Total phenological vegetation carbon (kg C/m2)
    !      (leafC+woodC+rootC)
,c_veg_bal                                                                     &
    ! WORK Total balanced vegetation carbon (kg C/m2)
    !      (leafC_bal+woodC+rootC)
,df_nsc_expl                                                                   &
    ! WORK Explicit increment to NSC mass fractioon
,j_fnsc
    ! WORK Jacobian of f_nsc rate equation evaluated using previous time-step
    !      f_nsc value.

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  forw_gamma = 1.0
    ! Forward time-step weighting. (0 = explicit), (1 = approx implicit)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_NSC_POOL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Calculate temperature function
!------------------------------------------------------------------------------
f_temp = q10_leaf(ft)**( 0.1 * (tstar - zerodegc - 25.0) )

!------------------------------------------------------------------------------
! Calculate temperature dependent coefficients
!------------------------------------------------------------------------------
! The turgor limit acts on structural production (g0) only: respiration
! (r0) keeps its NSC and temperature dependence.
r0    = sug_g0(ft)   * f_temp * ( 1.0 - sug_yg(ft) ) / sug_yg(ft)
g0    = sug_g0(ft)   * f_temp * f_turgor
grec0 = sug_grec(ft) * f_temp

!-----------------------------------------------------------------------------
! Calculate balanced and phenological biomass
!-----------------------------------------------------------------------------
c_veg      = leafc     + woodc + rootc
c_veg_bal  = leafc_bal + woodc + rootc

!-----------------------------------------------------------------------------
! Calculate increment to NSC mass fraction
!-----------------------------------------------------------------------------
j_fnsc      = (2.0* r0 * f_nsc - (r0 + g0 + grec0))* c_veg / c_veg_bal         &
              - gpp / c_veg_bal

df_nsc_expl = timestep * (1.0 / c_veg_bal)                                     &
              * (( 1.0 - f_nsc ) * (gpp - resp_p) - growth)
df_nsc_corr = df_nsc_expl / (1.0 - forw_gamma * j_fnsc)
df_nsc_corr = MIN(1.0 - f_nsc, df_nsc_corr)
df_nsc_corr = MAX(-f_nsc,df_nsc_corr)

!-----------------------------------------------------------------------------
! Increment f_nsc to next timestep value
!-----------------------------------------------------------------------------
f_nsc = f_nsc + df_nsc_corr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_NSC_pool

SUBROUTINE sugar_flux_correction(ft, leafc_bal, leafc, woodc, rootc, df_nsc,   &
                                 resp_l, resp_w, resp_r, resp_p_m, resp_p_g,   &
                                 resp_p, growth, tstar, gpp, f_turgor)

USE pftparm, ONLY:                                                             &
! imported parameters
    sug_grec, sug_g0, sug_yg, q10_leaf

USE conversions_mod, ONLY: zerodegc

USE timestep_mod, ONLY: timestep

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                         &
 ft
    ! IN Plant functional type

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tstar                                                                         &
    ! IN Surface temperature (K)
,gpp                                                                           &
    ! IN Gross primary productivity (kg C/m2/s)
,leafc                                                                         &
    ! IN Phenological leaf carbon (kg/m2)
,leafc_bal                                                                     &
    ! IN Balanced leaf carbon (kg/m2)
,woodc                                                                         &
    ! IN Wood carbon (kg/m2)
,rootc                                                                         &
    ! IN Root carbon (kg/m2)
,df_nsc                                                                        &
    ! IN Increment change in NSC mass fraction
,f_turgor
    ! IN Turgor limit on structural growth (0-1)

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 resp_w                                                                        &
    ! OUT Wood maintenance respiration rate (kg C/m2/s).
,resp_r                                                                        &
    ! OUT Root maintenance respiration rate (kg C/m2/s)
,resp_l                                                                        &
    ! OUT Leaf maintenance respiration rate (kg C/m2/s)
,resp_p_m                                                                      &
    ! OUT Plant maintenance respiration (kg C/m2/s)
,resp_p_g                                                                      &
    ! OUT Plant growth respiration (kg C/m2/s)
,resp_p                                                                        &
    ! OUT Plant respiration (kg C/m2/s)
,growth
    ! OUT Plant growth rate (kg C/m2/s)

!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 r0                                                                            &
    ! WORK Specific rate of NSC utilisation by respiration
,g0                                                                            &
    ! WORK Specific rate of NSC utilisation by structural growth
,grec0                                                                         &
    ! WORK Specific rate of structural C decomposition
,c_veg                                                                         &
    ! WORK Total phenological vegetation carbon (kg C/m2)
    !      (leafC+woodC+rootC)
,c_veg_bal                                                                     &
    ! WORK Total balanced vegetation carbon (kg C/m2)
    !      (leafC_bal+woodC+rootC)
,f_temp                                                                        &
    ! WORK Temperature function used for respiration and growth rates
,f_nsc_av                                                                      &
    ! WORK Average NSC mass fraction over time-step that balances rate of
    !      change equation
,a                                                                             &
    ! WORK Coefficient of second order term in quadratic function for
    !      average f_nsc over time-step
,b                                                                             &
    ! WORK Coefficient of first order term in quadratic function for
    !      average f_nsc over time-step
,c
    ! WORK Constant in quadratic function for average f_nsc over time-step

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SUGAR_FLUX_CORRECTION'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Calculate temperature function
!------------------------------------------------------------------------------
f_temp = q10_leaf(ft)**( 0.1 * (tstar - zerodegc - 25.0) )

!------------------------------------------------------------------------------
! Calculate temperature dependent coefficients
!------------------------------------------------------------------------------
! The turgor limit acts on structural production (g0) only: respiration
! (r0) keeps its NSC and temperature dependence.
r0    = sug_g0(ft)   * f_temp * ( 1.0 - sug_yg(ft) ) / sug_yg(ft)
g0    = sug_g0(ft)   * f_temp * f_turgor
grec0 = sug_grec(ft) * f_temp

!-----------------------------------------------------------------------------
! Calculate balanced and phenological biomass
!-----------------------------------------------------------------------------
c_veg      = leafc     + woodc + rootc
c_veg_bal  = leafc_bal + woodc + rootc

!-----------------------------------------------------------------------------
! Flux correction
! - need to recalculate fluxes to match corrected f_nsc increment
!-----------------------------------------------------------------------------
! Calculate average f_nsc value over timestep that corresponds to increment.
! This is the root of the rate equation of update_NSC_pool,
!   df_nsc/timestep = ((1 - f)(gpp - r0 f c_veg) - (g0 f - grec0 (1 - f)) c_veg)
!                     / c_veg_bal,
! so gpp is divided by c_veg_bal, as there (it was divided by c_veg, which
! differs whenever lai /= lai_bal).
a = r0 * c_veg / c_veg_bal
b = -(r0 + g0 + grec0) * c_veg / c_veg_bal - gpp / c_veg_bal
c = gpp / c_veg_bal + grec0 * c_veg / c_veg_bal - df_nsc / timestep
f_nsc_av = ( -b - ( b**2.0 - 4.0 * a * c )**0.5 ) / ( 2.0 * a )

! Recalculate fluxes with average f_nsc
CALL calculate_sug_fluxes(ft, leafc, woodc, rootc, f_nsc_av, resp_l,           &
                          resp_w, resp_r, resp_p_m, resp_p_g,                  &
                          resp_p, growth, tstar, f_turgor)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sugar_flux_correction

!------------------------------------------------------------------------------
! Sink demand for the SUGAR_stress cap (som_nsc_feedback = 2 or 3): the gross
! photosynthesis that moves f_nsc towards f_full over tau_fill. For
! sugar_model = 1, from the rate equation of update_NSC_pool,
!   c_veg_bal df_nsc/dt = (1 - f)(gpp - resp_p) - growth,
! setting df_nsc/dt = (f_full - f)/tau_fill gives
!   gpp = resp_p + (growth + c_veg_bal (f_full - f)/tau_fill) / (1 - f).
! Carbon beyond this has no use (the pool is full and the sinks are
! saturated). Floored at 0. Fluxes are from the f_nsc passed in (the
! previous time-step value when called before sugar).
!------------------------------------------------------------------------------
FUNCTION sugar_sink_demand(ft, leafc, leafc_bal, woodc, rootc, f_nsc, tstar,   &
                           f_turgor, f_full, tau_fill, flush) RESULT( demand )

USE jules_vegetation_mod, ONLY: sugar_model, sugar_mm
USE pftparm, ONLY: sug_yg

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft
    ! Plant functional type
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  leafc, leafc_bal, woodc, rootc,                                              &
    ! Phenological and balanced leaf, wood and root carbon (kg C/m2)
  f_nsc,                                                                       &
    ! NSC mass fraction (kg kg-1)
  tstar,                                                                       &
    ! Surface temperature (K)
  f_turgor,                                                                    &
    ! Turgor limit on structural growth (0-1)
  f_full,                                                                      &
    ! f_nsc at which the pool is full (kg kg-1)
  tau_fill,                                                                    &
    ! Time to fill the pool (s)
  flush
    ! Leaf-flush construction (kg C/m2/s), sugar_model = 2
REAL(KIND=real_jlslsm) :: demand
    ! Gross photosynthesis the plant can use (kg C/m2/s)

REAL(KIND=real_jlslsm) :: resp_l, resp_w, resp_r, resp_p_m, resp_p_g,          &
                          resp_p, growth, c_veg_bal, u_max, sat

c_veg_bal = leafc_bal + woodc + rootc

IF ( sugar_model == sugar_mm ) THEN
  ! Jones et al. (2020): c_veg_bal dW/dt = gpp - flush/Yg - U(W), so the gpp
  ! that moves W towards f_full over tau_fill is
  ! U(W) + flush/Yg + c_veg_bal (f_full - W)/tau_fill.
  CALL sugar_mm_rates(ft, leafc, woodc, rootc, f_nsc, tstar, f_turgor,         &
                      u_max, sat)
  demand = u_max * sat + flush / sug_yg(ft)                                    &
           + c_veg_bal * ( f_full - f_nsc ) / tau_fill
ELSE
  CALL calculate_sug_fluxes(ft, leafc, woodc, rootc, f_nsc, resp_l,            &
                            resp_w, resp_r, resp_p_m, resp_p_g,                &
                            resp_p, growth, tstar, f_turgor)
  demand = resp_p + ( growth + c_veg_bal * ( f_full - f_nsc ) / tau_fill )     &
                    / MAX( 1.0 - f_nsc, 1.0e-6_real_jlslsm )
END IF
demand = MAX( demand, 0.0_real_jlslsm )

END FUNCTION sugar_sink_demand

!------------------------------------------------------------------------------
! Jones et al. (2020) SUGAR (sugar_model = 2): rate at pool saturation and the
! saturation term. With W = f_nsc = C_NSC / c_veg_bal,
!   G  = G0  FQ(T) f_turgor c_veg W / (W + Km)          growth, Eqn 4
!   Rm = Rm0 FQ(T)          c_veg W / (W + Km)          maintenance, Eqn 8
!   Rg = (1 - Yg) / Yg G                                growth resp., Eqn 7
!   U  = G + Rm + Rg = u_max W / (W + Km)               utilisation, Eqn 11
! so u_max = (Rm0 + f_turgor G0 / Yg) FQ c_veg. The turgor limit acts on
! growth and its respiration only; maintenance continues.
!------------------------------------------------------------------------------
SUBROUTINE sugar_mm_rates(ft, leafc, woodc, rootc, f_nsc, tstar, f_turgor,     &
                          u_max, sat)

USE pftparm, ONLY: sug_g0, sug_rm0, sug_km, sug_yg, q10_leaf
USE conversions_mod, ONLY: zerodegc

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  leafc, woodc, rootc, f_nsc, tstar, f_turgor
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  u_max,                                                                       &
    ! Utilisation at pool saturation (kg C/m2/s)
  sat
    ! W / (W + Km)

REAL(KIND=real_jlslsm) :: f_temp

f_temp = q10_leaf(ft)**( 0.1 * (tstar - zerodegc - 25.0) )
u_max  = ( sug_rm0(ft) + f_turgor * sug_g0(ft) / sug_yg(ft) ) * f_temp         &
         * ( leafc + woodc + rootc )
sat    = MAX( f_nsc, 0.0 ) / ( MAX( f_nsc, 0.0 ) + sug_km(ft) )

END SUBROUTINE sugar_mm_rates

!------------------------------------------------------------------------------
! One time step of the Jones et al. (2020) SUGAR (sugar_model = 2):
!   c_veg_bal dW/dt = gpp - flush/Yg - U(W),
! linearised implicitly in W, as update_NSC_pool. The leaf flush (new leaf
! carbon, l_sugar_leaf_flush) is built from the pool with growth
! respiration. The fluxes are then those of the utilisation the increment
! implies, U_eff = gpp - flush/Yg - c_veg_bal dW/dt, split by the fixed
! ratios G : Rm : Rg of the time step, so the carbon balance closes exactly.
! Plant respiration is shared between leaf, wood and root by their carbon.
!------------------------------------------------------------------------------
SUBROUTINE sugar_mm_step(ft, leafc, leafc_bal, woodc, rootc, f_nsc, resp_l,    &
                         resp_w, resp_r, resp_p_m, resp_p_g, resp_p, growth,   &
                         tstar, gpp, f_turgor, flush)

USE pftparm, ONLY: sug_g0, sug_rm0, sug_km, sug_yg, q10_leaf
USE conversions_mod, ONLY: zerodegc
USE timestep_mod, ONLY: timestep

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  leafc, leafc_bal, woodc, rootc, tstar, gpp, f_turgor, flush
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: f_nsc
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  resp_l, resp_w, resp_r, resp_p_m, resp_p_g, resp_p, growth

REAL(KIND=real_jlslsm) ::                                                      &
  u_max, sat, c_veg, c_veg_bal, f_temp, g_max, rm_max, dudf, df_nsc, u_eff,  &
  src
    ! gpp less the leaf-flush cost, flush/Yg (kg C/m2/s)

CALL sugar_mm_rates(ft, leafc, woodc, rootc, f_nsc, tstar, f_turgor,           &
                    u_max, sat)

c_veg     = leafc + woodc + rootc
c_veg_bal = leafc_bal + woodc + rootc

resp_l = 0.0; resp_w = 0.0; resp_r = 0.0
resp_p_m = 0.0; resp_p_g = 0.0; resp_p = 0.0; growth = 0.0
IF ( c_veg_bal <= 0.0 .OR. u_max <= 0.0 ) RETURN
src = gpp - flush / sug_yg(ft)

! Implicit increment of W.
dudf   = u_max * sug_km(ft) / ( MAX( f_nsc, 0.0 ) + sug_km(ft) )**2
df_nsc = timestep / c_veg_bal * ( src - u_max * sat )                          &
         / ( 1.0 + timestep / c_veg_bal * dudf )
df_nsc = MAX( df_nsc, -f_nsc )
f_nsc  = f_nsc + df_nsc

! Fluxes consistent with the increment.
u_eff  = src - c_veg_bal * df_nsc / timestep
f_temp = q10_leaf(ft)**( 0.1 * (tstar - zerodegc - 25.0) )
g_max  = f_turgor * sug_g0(ft) * f_temp * c_veg
rm_max = sug_rm0(ft) * f_temp * c_veg
growth   = u_eff * g_max / u_max
resp_p_m = u_eff * rm_max / u_max
! Leaf flush: structural growth with its growth respiration.
growth   = growth + flush
resp_p_g = growth * ( 1.0 - sug_yg(ft) ) / sug_yg(ft)
resp_p   = resp_p_m + resp_p_g
IF ( c_veg > 0.0 ) THEN
  resp_l = resp_p * leafc / c_veg
  resp_w = resp_p * woodc / c_veg
  resp_r = resp_p * rootc / c_veg
END IF

END SUBROUTINE sugar_mm_step

END MODULE sugar_mod
