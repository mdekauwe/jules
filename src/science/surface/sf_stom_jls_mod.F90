! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE sf_stom_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_STOM_MOD'

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE :: f_vcmax_state(:,:)
    ! Soil-water capacity factor of l_som_vcmax_psi after the recovery lag
    ! (land_pts, npft); < 0 until first set.

PRIVATE
PUBLIC sf_stom

CONTAINS

! *********************************************************************
! Routines to calculate the bulk stomatal resistance and the canopy
! CO2 fluxes.
!
! References:
!   Bernacchi et al., 2001, Plant Cell and Environment, 24, 253--259,
!      https://doi.org/10.1111/j.1365-3040.2001.00668.x.
!   Medlyn et al., 2002, Plant, Cell and Environment, 25: 1167--1179,
!     https://doi.org/10.1046/j.1365-3040.2002.00891.x.
!   Mercado et al., 2018, New Phytologist, 218: 1462--1477,
!     https://doi.org/10.1111/nph.15100.
!
! *********************************************************************

SUBROUTINE sf_stom  (land_pts,land_index                                       &
,                    veg_pts,veg_index                                         &
,                    ft,co2,co2_3d,co2_dim_len                                 &
,                    co2_dim_row,l_co2_interactive                             &
,                    fsmc_in,veg_state,ht,ipar,lai                             &
,                    canht,pstar                                               &
,                    q1,ra,tstar,o3,t_home_gb,t_growth_gb,psi_root_zone        &
,                    e_supply                                                  &
,                    can_rad_mod,ilayers,leaf_flux_mod,som_base_parm,faparv    &
,                    lwp_c                                                     &
,                    el,gpp,npp,resp_p,resp_l,resp_r,resp_w                    &
,                    growth_sug, f_nsc                                         &
,                    n_leaf,n_root,n_stem,lai_bal,gc                           &
,                    fapar_sun,fapar_shd,fsun                                  &
,                    fapar_sun_nir,fapar_shd_nir                               &
,                    flux_o3,fo3,fapar_diag,apar_diag,psi_leaf,cica_ratio      &
,                    leaf_k                                                    &
,                    tair,lw_down,u_wind,z1_wind                                &
,                    t_leaf_sun,t_leaf_shd,t_can                               &
,                    isoprene,terpene,methanol,acetone                         &
,                    open_index,open_pts                                       &
,                    carbon_gain,hydraulic_cost                                &
                    !New arguments replacing USE statements
                    !crop_vars_mod (IN)
,                   dvi_cpft,rootc_cpft)

USE leaf_mod, ONLY: leaf
USE leaf_limits_mod, ONLY: leaf_limits
USE leaf_processes_sox_mod, ONLY: leaf_processes_sox
USE bvoc_emissions_mod, ONLY: bvoc_emissions

USE conversions_mod, ONLY: zerodegc
USE theta_field_sizes, ONLY: t_i_length

USE jules_surface_types_mod, ONLY: nnpft, ncpft

USE pftparm, ONLY:                                                             &
        kmax_pft, conductance_b, conductance_c, kcrit, gcuticular, min_gl_pft
USE jules_vegetation_mod, ONLY:                                                &
! imported model ids. JBaguley
    leaf_flux_fsmc, leaf_flux_stom_opt,                                        &
! imported parameters
    photo_collatz, photo_farquhar, photo_sox_collatz, photo_johnson,           &
    stomata_medlyn, stomata_sox, stomata_desica, stomata_profit_max,           &
    stomata_sox_profit,                                                        &
    photo_adapt, photo_acclim, photo_adapt_acclim,                             &
    photo_act_model, photo_act_pft, photo_act_gb, n_photo_coef,                &
! imported scalars that are not changed
    dsj_coef, dsv_coef, jv25_coef, act_j_coef, act_v_coef,                     &
    l_bvoc_emis, l_fapar_diag, l_trait_phys, l_stem_resp_fix, l_o3_damage,     &
    l_scale_resp_pm, photo_acclim_model, photo_model, stomata_model, l_sugar,  &
    som_gl_max, l_som_supply_limit,                                            &
    l_som_cuticular_floor, l_som_gravity, l_red, l_som_rhizo_series,           &
    l_som_vcmax_psi,                                                           &
    l_leaf_temp, leaf_temp_iter, l_leaf_temp_gc_eq, leaf_aero_model

USE CN_utils_mod, ONLY:                                                        &
! imported procedures
    get_can_ave_fac, nleaf_from_lai

USE pftparm, ONLY:                                                             &
! imported arrays that are not changed
    a_wl, a_ws, act_jmax, act_vcmax, alpha_elec, b_wl, c3, deact_jmax,         &
    deact_vcmax, ds_jmax, ds_vcmax, eta_sl, kpar, nl0, nr_nl, ns_nl, omega,    &
    r_grow, sigl, lma, nmass, kn, knl, tupp, tlow,  q10_leaf, nsw, nr, hw_sw,  &
    jv25_ratio

USE ccarbon, ONLY:                                                             &
! imported scalar parameters
   epco2,epo2

USE c_rmol, ONLY: rmol

USE jules_surface_mod, ONLY:                                                   &
! imported scalar parameters
   iter,o2,cmass

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE crop_utils_mod, ONLY:                                                      &
   stemc_from_prognostics,                                                     &
   lma_from_prognostics
USE cropparm, ONLY: cfrac_l

USE qsat_mod, ONLY: qsat

USE ereport_mod, ONLY: ereport

USE veg3_field_mod, ONLY: veg_state_type

USE sugar_mod, ONLY: sugar


USE stom_opt_jls_mod, ONLY: stom_opt_mod

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: set_ksr_path, ksr

USE planet_constants_mod, ONLY: repsilon, g, r
USE water_constants_mod, ONLY: rho_water
USE desica_jls_mod, ONLY: desica_fw, desica_hydraulics, tuzet_fw,             &
                          desica_store_inputs
USE timestep_mod, ONLY: timestep



IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(in).
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                         &
 land_pts                                                                      &
                            ! IN Number of land points to be
!                                 !    processed.
,land_index(land_pts)                                                          &
                            ! IN Index of land points on the
!                                 !    P-grid.
,veg_pts                                                                       &
                            ! IN Number of vegetated points.
,veg_index(land_pts)                                                           &
                            ! IN Index of vegetated points
!                                 !    on the land grid.
,co2_dim_len                                                                   &
                            ! IN Length of a CO2 field row.
,co2_dim_row                ! IN Number of CO2 field rows.

INTEGER, INTENT(IN) ::                                                         &
 ft                         ! IN Plant functional type.

LOGICAL, INTENT(IN) :: l_co2_interactive   ! switch for 3D CO2 field

INTEGER, INTENT(IN) ::                                                         &
  can_rad_mod                                                                  &
!                           !Switch for canopy radiation model
 ,ilayers                                                                      &
!                           !No of layers in canopy radiation model
 ,leaf_flux_mod                                                                &
!                           !Switch for leaf flux model
 ,som_base_parm
!                           !Switch for stomatal optimisation model base
!                           !   physical parameter.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 co2                                                                           &
                            ! IN Atmospheric CO2 concentration
,co2_3d(co2_dim_len,co2_dim_row)                                               &
!                                 ! IN 3D atmos CO2 concentration
!                                 !    (kg CO2/kg air).
,fsmc_in(land_pts)                                                             &
                            ! IN Soil water factor.
,ht(land_pts)                                                                  &
                            ! IN Canopy height (m).
,ipar(land_pts)                                                                &
                            ! IN Incident PAR (W/m2).
,lai(land_pts)                                                                 &
                            ! IN Leaf area index.
,canht(land_pts)                                                               &
                            ! IN Canopy Height
,pstar(land_pts)                                                               &
                            ! IN Surface pressure (Pa).
,faparv(land_pts,ilayers)                                                      &
                            ! IN Profile of absorbed PAR.
,fapar_shd(land_pts,ilayers)                                                   &
                            ! IN Profile of absorbed DIFF_PAR.
,fapar_sun(land_pts,ilayers)                                                   &
                            ! IN Profile of absorbed DIR_PAR.
,fsun(land_pts,ilayers)                                                        &
                            ! IN fraction of sunlit leaves
,fapar_sun_nir(land_pts,ilayers)                                               &
,fapar_shd_nir(land_pts,ilayers)                                               &
                            ! IN as fapar_sun, fapar_shd for near-IR
                            ! (l_leaf_temp)
,q1(land_pts)                                                                  &
                            ! IN Specific humidity at level 1
,ra(land_pts)                                                                  &
                            ! IN Aerodynamic resistance (s/m).
,tstar(land_pts)                                                               &
                            ! IN Surface temperature (K).
,o3(land_pts)                                                                  &
                            ! IN Surface ozone concentration (ppb).
,t_home_gb(land_pts)                                                           &
                            ! IN Static (home) temperature for adaptation of
                            ! photosynthesis (K).
,t_growth_gb(land_pts)                                                         &
                            ! IN Running mean (growth) temperature for
                            ! acclimation of photosynthesis (K).
,psi_root_zone(land_pts)                                                       &
                            ! IN Root zone water potential (Pa). Used by the
                            !    stomatal optimisation and by SOX
                            !    (stomata_model = stomata_sox).
,e_supply(land_pts)
                            ! IN Transpiration the soil can supply this
                            !    timestep (kg m-2 s-1); negative = no limit
                            !    (l_som_supply_limit, see physiol)

TYPE(veg_state_type), INTENT(IN OUT) :: veg_state

!-----------------------------------------------------------------------------
! Arguments with INTENT(out).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 el(land_pts)                                                                  &
                            ! OUT Transpiration rate (mol m-2 s-1) from the
                            !      stomatal optimisation model.
,gpp(land_pts)                                                                 &
                            ! OUT Gross Primary Productivity
!                                 !     (kg C/m2/s).
,npp(land_pts)                                                                 &
                            ! OUT Net Primary Productivity
!                                 !     (kg C/m2/s).
,resp_p(land_pts)                                                              &
                            ! OUT Plant respiration rate
!                                 !     (kg C/m2/sec).
,resp_r(land_pts)                                                              &
                            ! OUT Root respiration rate
!                                 !     (kg C/m2/sec).
,resp_l(land_pts)                                                              &
                            ! OUT Leaf maintanence respiration rate
!                                 !     (kg C/m2/sec).
,resp_w(land_pts)                                                              &
                            ! OUT Wood respiration rate
!                                 !     (kg C/m2/sec).
,growth_sug(land_pts)                                                          &
                            ! OUT Net structural C growth rate (SUGAR only)
                                  !     (kg C/m2/sec).
,flux_o3(land_pts)                                                             &
                            ! OUT Flux of O3 to stomata (nmol O3/m2/s).
,fo3(land_pts)                                                                 &
                            ! OUT Ozone exposure factor.
,fapar_diag(land_pts)                                                          &
                            ! OUT FAPAR diagnostic
,apar_diag(land_pts)                                                           &
                            ! OUT APAR diagnostic
,psi_leaf(land_pts)                                                            &
                            ! OUT Leaf water potential (Pa)
,cica_ratio(land_pts)                                                          &
                            ! OUT Ratio of intecellular and atmospheric CO2
,leaf_k(land_pts)                                                              &
                            ! OUT Xylem conductance at leaf level (m/s).
,lwp_c(land_pts)
                            ! OUT Canopy leaf water potential (MPa), SOX

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 t_leaf_sun(land_pts)                                                          &
,t_leaf_shd(land_pts)                                                          &
,t_can(land_pts)
                            ! INOUT Sunlit, shaded and canopy (leaf-area
                            ! weighted) leaf temperature used by
                            ! photosynthesis (K): tstar unless l_leaf_temp.
                            ! The input (last timestep's) is the first
                            ! guess with l_leaf_temp.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 tair(land_pts)                                                                &
                            ! IN Level-1 air temperature (K).
,lw_down(land_pts)                                                             &
                            ! IN Downward longwave radiation (W m-2).
,u_wind(land_pts)                                                              &
                            ! IN Level-1 wind speed (m s-1).
,z1_wind(land_pts)
                            ! IN Height of the level-1 wind (m).


REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 gc(land_pts)                                                                  &
                            ! INOUT Canopy resistance to H2O (m/s).
                            ! The input value is only used if can_rad_mod=1.
,f_nsc(land_pts)
                            ! INOUT Non-structural carbohydrate mass fraction
                            !      (kgC/kgC)

! BVOC variables
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 isoprene(land_pts)                                                            &
                   ! OUT Isoprene Emission Flux (kgC/m2/s)
,terpene(land_pts)                                                             &
                   ! OUT (Mono-)Terpene Emission Flux (kgC/m2/s)
,methanol(land_pts)                                                            &
                   ! OUT Methanol Emission Flux (kgC/m2/s)
,acetone(land_pts)
                   ! OUT Acetone Emission Flux (kgC/m2/s)

INTEGER, INTENT(OUT) ::                                                        &
open_index(land_pts)                                                           &
                            ! OUT Index of land points
!                                 !      with open stomata.
,open_pts                   ! OUT Number of land points
!                                 !      with open stomata.

! TEMPORARY: output variables for testing
REAL(KIND=real_jlslsm) ::                                                      &
 carbon_gain(land_pts)                                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(land_pts)
                            ! Hydraulic cost for each leaf state

!New arguments replacing USE statements
!crop_vars_mod (IN)
REAL(KIND=real_jlslsm), INTENT(IN) :: dvi_cpft(land_pts,ncpft)
REAL(KIND=real_jlslsm), INTENT(IN) :: rootc_cpft(land_pts,ncpft)

!-----------------------------------------------------------------------------
! Local parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  cconu = 12.0e-3,                                                             &
    ! kg C in 1 mol CO2.
  conpar = 2.19e5,                                                             &
    ! Conversion from mol s-1 to W for PAR (J/mol photons).
  t_ref = zerodegc + 25.0,                                                     &
    ! Reference temperature (K).
  tref_rmol = t_ref * rmol
    ! The product of t_ref and rmol (J mol-1).

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 i,j,k,l,m,n                                                                   &
                            ! WORK Loop counters.
,clos_pts                                                                      &
                            ! WORK Number of land points
!                                 !      with closed stomata.
,errcode                                                                       &
                            ! Error code to pass to ereport.
,pft_photo_model
                            ! Indicates which photosynthesis model to use for
                            ! the current PFT.

REAL(KIND=real_jlslsm) ::                                                      &
 dq_min                                                                        &
   ! Minimum-allowed specific humidity deficit (kg H20 vapour/kg air).
,expkn                                                                         &
   ! Decay term.
,fstem                                                                         &
   ! Ratio of respiring stem wood to total wood.
,lma_tmp                                                                       &
   ! Temporary leaf mass per area for crops (kg leaf per m2 leaf area).
,stem_resp_scaling                                                             &
   ! Scaling factor to reduce stem respiration
,stemc                                                                         &
   ! Stem carbon (kg m-2)
,sun_term                                                                      &
   ! Conversion from PAR to electron flux (mol electrons J-1).
,th_degc, tg_degc
   ! Temperatures t_home_gb and t_growth_gb in degrees Celsius.

!-----------------------------------------------------------------------------
! Local array variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 clos_index(land_pts)
                            ! WORK Index of land points
!                                 !      with closed stomata.

REAL(KIND=real_jlslsm) ::                                                      &
 anetc(land_pts)                                                               &
                            ! WORK Net canopy photosynthesis
!                                 !     (mol CO2/m2/s).
,co2c(land_pts)                                                                &
                            ! WORK Canopy level CO2 concentration
!                                 !      (kg CO2/kg air).
,ci(land_pts)                                                                  &
                            ! WORK Internal CO2 pressure (Pa).
,dq(land_pts)                                                                  &
                            ! WORK Specific humidity deficit
!                                 !      (kg H2O/kg air).
,dqc(land_pts)                                                                 &
                            ! WORK Canopy level specific humidity
!                                 !      deficit (kg H2O/kg air).
,fpar(land_pts)                                                                &
                            ! WORK PAR absorption factor.
,i2_shd (land_pts)                                                             &
                            ! Radiation that goes to Photosystem II for
                            ! shaded leaves, expressed as an electron flux
                            ! (mol electrons m-2 s-1).
,i2_sun(land_pts)                                                              &
                            ! Radiation that goes to Photosystem II for
                            ! sunlit leaves, expressed as an electron flux
                            ! (mol electrons m-2 s-1).
,lai_bal(land_pts)                                                             &
                            ! WORK Leaf area index in balanced
!                                 !      growth state.
,nleaf_top(land_pts)                                                           &
                            ! WORK Nitrogen concentration of top leaf.
!                           ! if(l_trait_phys)= g N/m2
                            ! else= kg N/kg C
,n_leaf(land_pts)                                                              &
                            ! WORK Nitrogen contents of the leaf,
,n_root(land_pts)                                                              &
                            !      root,
,n_stem(land_pts)                                                              &
                            !      and stem (kg N/m2).
,leafc(land_pts)                                                               &
                            ! WORK Leaf Carbon (kg C/m2)
,leafc_bal(land_pts)                                                           &
                            ! WORK Carbon of balanced LAI (kg C/m2)
,woodc(land_pts)                                                               &
                            ! WORK Wood Carbon (kg C/m2)
,rootc(land_pts)                                                               &
                            ! WORK Root Carbon (kg C/m2).
,qs(land_pts)                                                                  &
                            ! WORK Saturated specific humidity
!                                 !      (kg H2O/kg air).
,ra_rc(land_pts)                                                               &
                            ! WORK Ratio of aerodynamic resistance
!                                 !      to canopy resistance.
,gc_ml_prev(land_pts)                                                          &
                            ! WORK Multilayer profit max: canopy
!                                 !      conductance (m s-1) from the previous
!                                 !      pass of the humidity-deficit
!                                 !      iteration.
,rdc(land_pts)                                                                 &
                            ! WORK Canopy dark respiration,
!                                 !      without soil water dependence
!                                 !      (mol CO2/m2/s).
,rdmean(land_pts)                                                              &
                            ! WORK Mean dark respiration
!                                 !      per unit leaf area
!                                 !      over canopy
!                                 !      without soil water dependence
!                                 !   (mol CO2/s per m2 leaf area).
,nlmean(land_pts)                                                              &
                            ! WORK Mean nitrogen per unit leaf area
!                                 !      over canopy
!                                 !      without soil water dependence
!                                 !      (kg N per m2 leaf area).
!                                 !      Used in place of n_leaf when LAI
!                                 !      is very small to avoid blowing up
!                                 !      the calculation of respiration.
,resp_p_g(land_pts)                                                            &
                            ! WORK Plant growth respiration rate
!                                 !      (kg C/m2/sec).
,resp_p_m(land_pts)                                                            &
                            ! WORK Plant maintenance respiration
!                                 !      rate (kg C/m2/sec).
,root(land_pts)                                                                &
                            ! WORK Root carbon (kg C/m2).
,faparv_layer(land_pts,ilayers)                                                &
                            ! WORK absorbed par(layers)
,flux_o3_l(land_pts)                                                           &
                            ! WORK Flux of O3 to stomata (nmol O3/m2/s).
,flux_o3_l_sun(land_pts)                                                       &
                            ! WORK Flux of O3 to stomata
                            !      for sunlit leaves
!                             !      (for can_rad_mod=5)
                            !      (nmol O3/m2/s).
,flux_o3_l_shd(land_pts)                                                       &
                            ! WORK Flux of O3 to stomata
                            !      for shaded leaves
                            !      (for can_rad_mod=5)
                            !      (nmol O3/m2/s).
,fo3_l(land_pts)                                                               &
                            ! WORK Ozone exposure factor.
,fo3_l_sun(land_pts)                                                           &
                            ! WORK Ozone exposure factor
                            !      for sunlit leaves
                            !      (for can_rad_mod=5)
,fo3_l_shd(land_pts)                                                           &
                            ! WORK Ozone exposure factor
                            !      for shaded leaves
                            !      (for can_rad_mod=5)
,o3mol(land_pts)                                                               &
                            ! WORK Surface ozone concentration (moles).
,fsmc_scale(land_pts)                                                          &
                            ! WORK Scaling of stem and root plant maintenance
                            ! respiration.
,denom(land_pts)                                                               &
   ! Denominator in temperature-dependency of Vcmax (Collatz model only).
,dsj(land_pts)                                                                 &
   ! Entropy factor for Jmax, including any acclimation (J mol-1 K-1).
,dsv(land_pts)                                                                 &
   ! Entropy factor for Vcmax, including any acclimation (J mol-1 K-1).
,actj(land_pts)                                                                &
   ! Activation energy for Jmax, including any acclimation (J mol-1).
,actv(land_pts)                                                                &
   ! Activation energy for Vcmax, including any acclimation (J mol-1).
,ccp(land_pts)                                                                 &
   ! Photorespiratory compensatory point (Pa). This is zero for C4 plants.
,i2(land_pts)                                                                  &
   ! Radiation that goes to Photosystem II, expressed as an electron flux
   ! (mol electrons m-2 s-1).
,je(land_pts)                                                                  &
   ! Electron transport rate (mol m-2 s-1).
,je_shd(land_pts)                                                              &
   ! Electron transport rate for shaded leaves (mol m-2 s-1).
,je_shd_ratio(land_pts)                                                        &
   ! Ratio je_shd : je.
,je_sun(land_pts)                                                              &
   ! Electron transport rate for sunlit leaves (mol m-2 s-1).
,je_sun_ratio(land_pts)                                                        &
   ! Ratio je_sun : je.
,jmax(land_pts)                                                                &
   ! Maximum rate of electron transport (mol CO2 m-2 s-1).
,jv25(land_pts)                                                                &
   ! Ratio of Jmax to Vcmax at 25 degC, including any acclimation.
,kc(land_pts)                                                                  &
   ! Michaelis-Menten constant for CO2 (Pa).
,km(land_pts)                                                                  &
   ! A combination of Michaelis-Menten and other terms.
,ko(land_pts)                                                                  &
   ! Michaelis-Menten constant for O2 (Pa).
,jmax_temp(land_pts)                                                           &
   ! Factor expressing the effect of temperature on Jmax.
,qtenf_term(land_pts)                                                          &
   ! Q10 temperature term used for Vcmax.
,vcmax_temp(land_pts)                                                          &
   ! Factor expressing the effect of temperature on Vcmax.
,vcmax(land_pts)                                                               &
   ! Maximum rate of carboxylation of Rubisco (mol CO2/m2/s).
,vcmaxc(land_pts)                                                              &
   ! Canopy-scaled vcmax for the can_rad_mod=1 stomatal optimisation branch
   ! (vcmax(l) * fpar(l)). Kept separate from vcmax so that repeated passes
   ! of the humidity-deficit iteration (DO k = 1,iter) rescale from the
   ! unscaled leaf-level vcmax each time, rather than compounding fpar into
   ! vcmax itself.
,jmaxc(land_pts)                                                               &
   ! Canopy-scaled jmax, analogous to vcmaxc above.
,acrc(land_pts)                                                                &
,i2c(land_pts)                                                                 &
   ! Canopy-scaled acr and i2 (x fpar), analogous to vcmaxc above: the
   ! big leaf's absorbed light scales with fpar like its capacity.
,anetl(land_pts)                                                               &
                            ! WORK Net leaf photosynthesis
!                                 !      (mol CO2/m2/s/LAI).
,anetl_sun(land_pts)                                                           &
!                                 ! WORK Net leaf photosynthesis of
!                                 !      sunlit leaves
!                                 !      (mol CO2/m2/s/LAI)
,anetl_shd(land_pts)                                                           &
!                                 ! WORK Net leaf photosynthesis of
!                                 !      shaded leaves
!                                 !      (mol CO2/m2/s/LAI
,apar(land_pts)                                                                &
!                                 ! WORK PAR absorbed by the top leaf
!                                 !      (W/m2).
,acr(land_pts)                                                                 &
                            ! WORK Absorbed PAR
!                                 !      (mol photons/m2/s).
,ca(land_pts)                                                                  &
                            ! WORK Canopy level CO2 pressure
!                                 !      (Pa).
,gl(land_pts)                                                                  &
                            ! WORK Leaf conductance for H2O
!                                 !      (m/s).
,gl_sun(land_pts)                                                              &
                            ! WORK Leaf conductance for H2O of
!                                 !      sunlit leaves (m/s).
,gl_shd(land_pts)                                                              &
                            ! WORK Leaf conductance for H2O of
!                                 !      shaded leaves (m/s).
,el_sun(land_pts)                                                              &
                            ! WORK Leaf transpiration
!                                 !      shaded leaves (mol m-2 s-1)
,el_shd(land_pts)                                                              &
                            ! WORK Leaf transpiration
!                                 !      shaded leaves (mol m-2 s-1)
,leaf_k_sun(land_pts)                                                          &
                            ! WORK Xylem conductance at leaf level (m/s).
,leaf_k_shd(land_pts)                                                          &
                            ! WORK Xylem conductance at leaf level (m/s).
,HC_sun(land_pts)                                                            &
                            ! WORK hydraulic cost for each leaf state
,HC_shd(land_pts)                                                            &
                            ! WORK hydraulic cost for each leaf state
,CG_sun(land_pts)                                                            &
                            ! WORK Carbon gain for each leaf state
,CG_shd(land_pts)                                                            &
                            ! WORK Carbon gain for each leaf state
,psi_leaf_sun(land_pts)                                                        &
                            ! WORK Leaf water potential for
!                                 !      sunlit leaves (Pa).
,psi_leaf_shd(land_pts)                                                        &
                            ! WORK Leaf water potential for
!                                 !      shaded leaves (Pa).
,icr(land_pts)                                                                 &
                            ! WORK Incident PAR (mol photons/m2/s).
,oa(land_pts)                                                                  &
                            ! WORK Atmospheric O2 pressure
!                                 !      (Pa).
,rd(land_pts)                                                                  &
                            ! WORK Dark respiration, including any effect of
                            !      light inhibition (mol CO2/m2/s).
,rd_dark(land_pts)                                                             &
                            ! WORK Dark respiration, excluding effect of
                            !      light inhibition (mol CO2/m2/s).
,rd_sun(land_pts)                                                              &
                            ! WORK Dark respiration of sunlit leaves
!                                 !      (mol CO2/m2/s).
,rd_shd(land_pts)                                                              &
                            ! WORK Dark respiration of shaded leaves
!                                 !      (mol CO2/m2/s).

,wcarb(land_pts)                                                               &
                            ! WORK Carboxylation, ...
,wlite(land_pts)                                                               &
                            !      ... Light, and ...
,wexpt(land_pts)                                                               &
                            !      ... export limited gross ...
!                                 !      ... photosynthetic rates ...
!                                 !      ... (mol CO2/m2/s).
,wlitev(land_pts)                                                              &
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for each layer
!                                 !      (mol CO2/m2/s).
,wlitev_sun(land_pts)                                                          &
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for sunlit leaves
!                                 !      (mol CO2/m2/s).
,wlitev_shd(land_pts)                                                          &
!                                 ! WORK Light limited gross
!                                 !      photosynthetic rates
!                                 !      for shaded leaves
!                                 !      (mol CO2/m2/s).
,dlai(land_pts)                                                                &
                            ! WORK LAI Increment.
,nleaf_layer(land_pts)                                                         &
                            ! WORK Leaf nitrogen concentration in a layer.
                            ! kgN/kgC if l_trait_phys=F
                            ! gN/m2   if l_trait_phys=T.
,can_averaging_fac(land_pts)
                            ! WORK factor to convert top of canopy
                            ! value to canopy average.

! 29 Apr, MGDK
REAL(KIND=real_jlslsm) :: je_dummy(land_pts), fapar_dummy(land_pts)
REAL(KIND=real_jlslsm) :: kmax_per_lyr(land_pts), kcrit_per_lyr(land_pts)
                            ! kcrit_per_lyr scales kcrit(ft) by the same
                            ! fraction of whole-plant kmax that each layer
                            ! received, so the critical-conductance
                            ! threshold stays at (1 - kcrit_fractional_loss)
                            ! of THAT layer's own capacity.
REAL(KIND=real_jlslsm) :: kmax_leaf_lyr(land_pts), kcrit_leaf_lyr(land_pts)
REAL(KIND=real_jlslsm) :: kleaf_prof, kleaf_mean, kleaf_a
                            ! Conductance profile through the canopy
                            ! (multilayer): layer value, canopy mean and
                            ! per-layer decay exponent (see kmax_leaf_lyr).
                            ! kmax_leaf_lyr is the conductance per unit leaf
                            ! area of this layer's whole soil-to-leaf path
                            ! in the multilayer stomatal optimisation:
                            ! kmax_per_lyr redistributed by the leaf-N
                            ! profile (canopy mean kmax_pft), each layer on
                            ! its own parallel path from psi_src.
REAL(KIND=real_jlslsm) :: fsmc_leaf_resp(land_pts)
                            ! fsmc as applied to leaf dark respiration in the
                            ! GPP/respiration diagnostics below: fsmc for
                            ! leaf_flux_fsmc (where leaf() already scaled rd
                            ! by fsmc inside anetl), 1.0 for
                            ! leaf_flux_stom_opt (where stom_opt_mod uses
                            ! unscaled rd, al = wl - rd).
REAL(KIND=real_jlslsm) :: fsmc_unity(land_pts)
REAL(KIND=real_jlslsm) :: fsmc(land_pts)
                            ! Soil water factor applied to the leaf fluxes:
                            ! fsmc_in, or 1.0 for stomata_desica (stress acts
                            ! through psi_leaf only).
REAL(KIND=real_jlslsm) :: fsmc_lim(land_pts)
                            ! fsmc passed to leaf_limits: fsmc, or the Tuzet
                            ! factor fw for stomata_desica (leaf_limits sets
                            ! ci from it).
REAL(KIND=real_jlslsm) :: gl_cut_ds
                            ! DESICA canopy cuticular conductance (m s-1).
INTEGER, PARAMETER :: n_fw_bisect = 12
                            ! DESICA bisection steps on fw (to 2.4e-4).
INTEGER :: n_pass, i_pass
                            ! Passes of the big-leaf flux calculation.
REAL(KIND=real_jlslsm) :: ksr_ds(land_pts)
                            ! DESICA: soil-to-root conductance in series
                            ! (ksr, l_som_rhizo_series), < 0: none.
REAL(KIND=real_jlslsm) :: fw_lo(land_pts), fw_hi(land_pts),                    &
                          el_try(land_pts), psi_try(land_pts), k_try(land_pts),&
                          el_hyd(land_pts)
                            ! DESICA bisection bracket on fw, and the trial
                            ! transpiration (mol m-2 s-1), end-of-step
                            ! psi_leaf (Pa) and plant conductance, and the
                            ! transpiration the plant can deliver.
REAL(KIND=real_jlslsm) :: f_vc(land_pts)
                            ! Soil-water capacity factor (l_som_vcmax_psi).
REAL(KIND=real_jlslsm) :: psi_src(land_pts)
                            ! Water potential at the base of the plant path
                            ! seen by the stomatal optimisation (Pa):
                            ! psi_root_zone, less the gravitational drop
                            ! rho_water g ht when l_som_gravity.
REAL(KIND=real_jlslsm) :: gl_max_lf(land_pts), gl_max_bigleaf(land_pts)
                            ! som_gl_max on the basis each stom_opt_mod call
                            ! works on: per leaf area for the multilayer
                            ! calls, canopy (x LAI) for big-leaf.
REAL(KIND=real_jlslsm) :: ci_sun_ml(land_pts), ci_shd_ml(land_pts),            &
                          ci_gl_ml(land_pts), ci_lai_ml(land_pts)
                            ! Multilayer stomatal optimisation: ci of the
                            ! layer's sunlit / shaded leaf (Pa), and the
                            ! canopy sums of ci weighted by gl*dlai and by
                            ! dlai, giving the gl-weighted canopy ci (as
                            ! two-leaf), or the LAI mean when all are closed.
                            ! fsmc passed to leaf_limits in the multilayer
                            ! stomatal optimisation path: 1.0 everywhere, so
                            ! leaf_limits' fsmc == 0 closure test never
                            ! fires there - see the note at that call.
REAL(KIND=real_jlslsm) :: kmax_bigleaf(land_pts), kcrit_bigleaf(land_pts)
REAL(KIND=real_jlslsm) :: gl_max_eff(land_pts), share_sup(land_pts)
REAL(KIND=real_jlslsm) :: gl_cut(land_pts), gl_cut_eff(land_pts),              &
                          kmax_cut(land_pts), kcrit_cut(land_pts)
                            ! Cuticular floor (l_som_cuticular_floor):
                            ! canopy floor conductance (m s-1) before and
                            ! after the soil-supply cap, and the canopy
                            ! kmax/kcrit used to re-solve psi_leaf.
INTEGER, PARAMETER :: n_cut = 20
                            ! Floor fluxes tried, E_floor * i / n_cut, to
                            ! find the largest one the xylem can carry.
REAL(KIND=real_jlslsm) :: e_cut(n_cut, land_pts), psi_cut(n_cut, land_pts),    &
                          k_cut(n_cut, land_pts), el_cut(land_pts)
LOGICAL :: l_cut(land_pts)
INTEGER :: i_cut
                            ! gl_max for a stom_opt_mod call with the soil
                            ! supply cap applied (apply_supply_limit), and the
                            ! share of e_supply that call may use.

! Two-leaf (sunlit/shaded big-leaf) canopy, can_rad_mod = 7. Canopy (per m2
! ground) totals for each leaf class, aggregated from the layered radiation
! profile.
REAL(KIND=real_jlslsm) ::                                                      &
  lai_sun_2l(land_pts), lai_shd_2l(land_pts),                                  &
      ! Sunlit / shaded leaf area (m2 leaf m-2 ground).
  nw_sun_2l(land_pts), nw_shd_2l(land_pts),                                    &
      ! Leaf area weighted by the relative leaf-N profile exp(-kpar*L),
      ! i.e. the sunlit / shaded parts of the big-leaf fpar.
  apar_sun_2l(land_pts), apar_shd_2l(land_pts),                                &
      ! Absorbed PAR by the sunlit / shaded leaves (W m-2 ground).
  acr_sun_2l(land_pts), acr_shd_2l(land_pts),                                  &
      ! As apar_*_2l, in mol photons m-2 s-1.
  vcmax_sun_2l(land_pts), vcmax_shd_2l(land_pts),                              &
  jmax_sun_2l(land_pts), jmax_shd_2l(land_pts),                                &
      ! Canopy Vcmax / Jmax of each leaf class (mol m-2 ground s-1).
  kmax_sun_2l(land_pts), kmax_shd_2l(land_pts),                                &
  kcrit_sun_2l(land_pts), kcrit_shd_2l(land_pts),                              &
      ! Hydraulic conductance of each leaf class (mol m-2 s-1 Pa-1).
  gl_max_sun_2l(land_pts), gl_max_shd_2l(land_pts),                            &
      ! Maximum stomatal conductance of each leaf class (m s-1).
  ci_sun_2l(land_pts), ci_shd_2l(land_pts),                                    &
      ! Internal CO2 of each leaf class (Pa).
  f_sun_2l,                                                                    &
      ! Sunlit fraction of LAI.
  dnw_2l
      ! Layer integral of exp(-kpar*L) over dlai.

! Two-leaf leaf temperatures (l_leaf_temp).
REAL(KIND=real_jlslsm) ::                                                      &
  t_sun_lt(land_pts), t_shd_lt(land_pts),                                      &
      ! Sunlit / shaded leaf temperature of the current pass (K).
  t_c_lt(land_pts), q_c_lt(land_pts),                                          &
      ! Canopy air temperature (K) and specific humidity (kg kg-1).
  sw_sun_lt(land_pts), sw_shd_lt(land_pts),                                    &
      ! Absorbed PAR + near-IR of each class (W m-2 ground).
  wlw_sun_lt(land_pts), wlw_shd_lt(land_pts),                                  &
  wlws_sun_lt(land_pts), wlws_shd_lt(land_pts),                                &
      ! Layer weights of the isothermal long-wave (sky, soil-reflected).
  frad_sun_lt(land_pts), frad_shd_lt(land_pts),                                &
      ! Radiative conductance scalars of each class.
  ff_sun_lt(land_pts), ff_shd_lt(land_pts),                                    &
      ! Wind-weighted leaf area of each class, for forced convection.
  ccp_sun_lt(land_pts), ccp_shd_lt(land_pts),                                  &
  kc_sun_lt(land_pts), kc_shd_lt(land_pts),                                    &
  ko_sun_lt(land_pts), ko_shd_lt(land_pts),                                    &
  km_sun_lt(land_pts), km_shd_lt(land_pts),                                    &
      ! Temperature-dependent photosynthesis terms of each class.
  qs_sun_lt(land_pts), qs_shd_lt(land_pts),                                    &
  dq_sun_lt(land_pts), dq_shd_lt(land_pts),                                    &
      ! Saturated humidity at, and humidity deficit from, each leaf to the
      ! canopy air (kg kg-1).
  rdc_tstar(land_pts),                                                         &
      ! Canopy dark respiration at tstar, for the stem and root respiration.
  ra_lt(land_pts),                                                             &
      ! Canopy air - level 1 resistance for the leaf energy balance (s m-1).
  t_g_lt(land_pts), qs_g_lt(land_pts),                                         &
      ! First-guess canopy leaf temperature (K) and its saturated humidity
      ! (kg kg-1), for the first-guess canopy air humidity.
  zeta_lt(land_pts),                                                           &
      ! zref/L of the canopy (leaf_aero_model = 1), carried between passes.
  wind_ext(land_pts),                                                          &
      ! Extinction coefficient of the wind with leaf area in the canopy.
  z0m_a(land_pts), zref_a(land_pts), rt1ab_a(land_pts), zrd_a(land_pts),       &
  usc_a(land_pts), usuh_a(land_pts),                                           &
      ! CABLE canopy geometry (leaf_aero_model = 1).
  gc_cap, r_eq
      ! Cap and resistance for the transpiration-equivalent gc
      ! (l_leaf_temp_gc_eq).
REAL(KIND=real_jlslsm), PARAMETER :: m_h2o_lt = 0.018015
      ! Molar mass of water (kg mol-1).
INTEGER :: n_pass_2l
      ! Number of passes of the two-leaf optimisation.
                            ! kmax_pft(ft)/kcrit(ft) broadcast onto a
                            ! land_pts array for the big-leaf (l_multilayer
                            ! = .FALSE.) call to stom_opt_mod, whose kmax/
                            ! kcrit dummy arguments are land_pts-sized.
                            ! NOTE: passing the bare scalars kmax_pft(ft)/
                            ! kcrit(ft) here previously (pre-existing code)
                            ! relied on undefined Fortran sequence
                            ! association and could read out of bounds.
LOGICAL :: l_multilayer
! veg_pts_index is the trivial (1,2,3,...) map into veg_index that makes
! leaf_psi_jls's veg_index(open_index(j)) indexing select veg_index(j)
! directly, i.e. every vegetated point.
INTEGER :: veg_pts_index(land_pts)
REAL(KIND=real_jlslsm) :: k_eff(land_pts)



REAL(KIND=real_jlslsm) ::                                                      &
 act_j_tmp(n_photo_coef)                                                       &
   ! Coefficients governigng the acclimation of activation energy for Jmax.
,act_v_tmp(n_photo_coef)
   ! Coefficients governing the acclimation of activation energy for Vcmax.



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_STOM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Gravity (l_som_gravity). Lifting water to the top of the canopy costs
! rho_water g ht (0.01 MPa per m) on top of the frictional drop, so the
! stomatal optimisation draws on psi_root_zone - rho_water g ht, as in
! DESICA's psi_h (Xu et al. 2016, Eqn S2a). The whole canopy is put at the
! canopy height, the conservative end for the multilayer scheme. Applied as
! a shift of the source, the vulnerability curve is evaluated on the
! shifted potentials, i.e. at the pressure in the upper xylem.
!-----------------------------------------------------------------------------
psi_src(:) = psi_root_zone(:)
IF ( l_som_gravity .AND. leaf_flux_mod == leaf_flux_stom_opt ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    psi_src(l) = psi_root_zone(l) - rho_water * g * MAX(ht(l), 0.0)
  END DO
END IF

! Soil-water capacity factor of Vcmax/Jmax (l_som_vcmax_psi), once per
! timestep, so the recovery lag advances once whatever the number of leaf
! classes and leaf temperature iterations.
f_vc(:) = 1.0
IF ( l_som_vcmax_psi ) CALL vcmax_psi_factor( ft, land_pts, veg_pts,          &
                                              veg_index, psi_root_zone, f_vc )

!-----------------------------------------------------------------------------
! Initialisation.
!
! Note that the initialisation is appropriate for the existing code! If an
! existing variable is used differently (e.g. it is made available as a
! diagnostic) it might need to be treated differently - e.g. initialised.
!
! In general we only initialise output variables (diagnostics and/or those
! required elsewhere - so that these have a defined value at all location)
! and any local variables that need to be initialised (e.g. so as to allow
! accumulation over canopy layers). Many other variables are not initialised
! because they are always calculated at all locations for which a value is
! required (e.g. all vegetated points, or all points with open stomata).
!-----------------------------------------------------------------------------

!$OMP PARALLEL DO IF(land_pts > 1)                                             &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(land_pts,l_o3_damage,fo3,flux_o3,isoprene,                        &
!$OMP        terpene,methanol,acetone,fsmc_scale,                              &
!$OMP        stomata_model,lwp_c)                                              &
!$OMP PRIVATE(l)                                                               &
!$OMP SCHEDULE(STATIC)
DO l = 1, land_pts

  ! Ozone variables that are output.
  flux_o3(l)    = 0.0
  IF ( l_o3_damage ) THEN
    ! Initialise to zero to allow accumulation.
    fo3(l)      = 0.0
  ELSE
    ! Initialise to 1 to indicate no effect.
    fo3(l)      = 1.0
  END IF

  ! BVOC fluxes: values used at non-vegetated points or if fluxes are not
  ! calculated (l_bvoc_emis=F).
  isoprene(l)   = 0.0
  terpene(l)    = 0.0
  methanol(l)   = 0.0
  acetone(l)    = 0.0

  fsmc_scale(l) = 1.0  !  Value used if l_scale_resp_pm =.FALSE..
  ! leaf water potential variables for lwp_c diagnostic if SOX is used
  IF ( stomata_model == stomata_sox ) THEN
    ! initialise to zero to allow accumulation
    lwp_c(l) = 0.0
  END IF
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------------
! Initialise variables that are accumulated over layers.
! Note that some or all of these variables are also used by the big-leaf
! model (can_rad_mod=1) but do not need to be initialised in that case.
!-----------------------------------------------------------------------------
SELECT CASE ( can_rad_mod )
CASE ( 4,5,6 )
!$OMP PARALLEL DO IF(land_pts > 1)                                             &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(land_pts,anetc,gc,rdc,rdmean,psi_leaf,el,leaf_k,hydraulic_cost,   &
!$OMP        carbon_gain,kmax_per_lyr,kcrit_per_lyr,je_dummy,fapar_dummy)      &
!$OMP PRIVATE(l)                                                               &
!$OMP SCHEDULE(STATIC)
  DO l = 1, land_pts
    anetc(l)  = 0.0
    gc(l)     = 0.0
    rdc(l)    = 0.0
    rdmean(l) = 0.0

    ! gs_opt multi-layer, mgdk 9th June, 2026
    psi_leaf(l)  = 0.0
    el(l) = 0.0
    leaf_k(l) = 0.0
    hydraulic_cost(l) = 0.0
    carbon_gain(l) = 0.0
    kmax_per_lyr(l) = 0.0
    kcrit_per_lyr(l) = 0.0
    je_dummy(l) = 1.0
    fapar_dummy(l) = 1.0

  END DO
!$OMP END PARALLEL DO
END SELECT

!-----------------------------------------------------------------------------
! Decide which photosynthesis model will be used for this PFT.
!-----------------------------------------------------------------------------
SELECT CASE ( photo_model )
CASE ( photo_collatz )
  ! C3 and C4 both use a Collatz model.
  pft_photo_model = photo_collatz
CASE ( photo_farquhar )
  ! C3 uses Farquhar, C4 uses Collatz.
  IF ( c3(ft) == 1 ) THEN
    pft_photo_model = photo_farquhar
  ELSE
    pft_photo_model = photo_collatz
  END IF
CASE ( photo_sox_collatz )
  pft_photo_model = photo_collatz
  ! For the purposes of initialisation these are both photo_collatz
CASE ( photo_johnson )
  ! C3 uses Farquhar with Johnson-Berry electron transport, C4 uses Collatz.
  ! Johnson-Berry runs through the Farquhar code throughout: it differs only
  ! in Jmax (which holds Vqmax), je and wlite, which test photo_model
  ! (calc_photo_parameters, calc_electron_flux, leaf_limits, stom_opt).
  IF ( c3(ft) == 1 ) THEN
    pft_photo_model = photo_farquhar
  ELSE
    pft_photo_model = photo_collatz
  END IF
END SELECT

!-----------------------------------------------------------------------------
! Initialise ratios of electron fluxes.
!-----------------------------------------------------------------------------
IF ( pft_photo_model == photo_farquhar ) THEN
  SELECT CASE ( can_rad_mod )
  CASE ( 5, 6 )
!$OMP PARALLEL DO IF(land_pts > 1)                                             &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(land_pts, je_shd_ratio, je_sun_ratio)                             &
!$OMP PRIVATE(l)                                                               &
!$OMP SCHEDULE(STATIC)
    DO l = 1, land_pts
      je_shd_ratio(l)  = 0.0
      je_sun_ratio(l)  = 0.0
    END DO
!$OMP END PARALLEL DO
  END SELECT
END IF

!-----------------------------------------------------------------------------
! Set the CO2 compensation point to zero for C4 plants.
!-----------------------------------------------------------------------------
IF ( c3(ft) == 0 ) THEN
  ccp(:) = 0.0
END IF

!-----------------------------------------------------------------------------
! Calculate humidity at saturation.
!-----------------------------------------------------------------------------
CALL qsat(qs,tstar,pstar,land_pts)

! Set the minimum-allowed humidity deficit.
! (The stomatal optimisation keeps the dq_min it had when selected with
! leaf_flux_mod = 2 and stomata_model = 2.)
IF ( ( stomata_model == stomata_medlyn ) .OR. ( stomata_model == stomata_sox ) &
     .OR. ( stomata_model == stomata_profit_max )                              &
     .OR. ( stomata_model == stomata_sox_profit ) ) THEN
  ! Avoid dq=0 as this would cause the model to blow up.
  dq_min = 0.0001
ELSE
  dq_min = 0.0
END IF

!-----------------------------------------------------------------------------
! Set the canopy CO2 concentration.
!-----------------------------------------------------------------------------
!$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(i, j, l, m, lma_tmp)      &
!$OMP SHARED(l_co2_interactive, l_o3_damage, acr, veg_pts, veg_index,          &
!$OMP        land_index, t_i_length, co2c, co2_3d, co2, dq, qs, q1,            &
!$OMP        can_rad_mod, fpar, icr, kpar, lai, ft, apar, omega, ipar,         &
!$OMP        l_trait_phys, nnpft, dvi_cpft, nleaf_top, nmass, nl0,             &
!$OMP        dlai, ilayers, o3mol, o3, pstar, tstar, lma, ca, stomata_model,   &
!$OMP        dq_min, c3, oa)
IF ( l_co2_interactive ) THEN
  !       Use full 3D CO2 field.
!$OMP  DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)
    j = (land_index(l) - 1) / t_i_length + 1
    i = land_index(l) - (j-1) * t_i_length
    co2c(l) = co2_3d(i,j)
  END DO
!$OMP END DO NOWAIT
ELSE
  !       Use single CO2_MMR value.
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)
    co2c(l) = co2
  END DO
!$OMP END DO NOWAIT
END IF

!-----------------------------------------------------------------------------
! Calculate the surface to level 1 humidity deficit.
!-----------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO m = 1,veg_pts
  l = veg_index(m)
  dq(l) = MAX(dq_min,(qs(l) - q1(l)))
END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------------
! Calculate the PAR absorption factor
!-----------------------------------------------------------------------------
IF ( can_rad_mod == 1 ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)
    fpar(l) = (1.0 - EXP(-kpar(ft) * lai(l))) / kpar(ft)
  END DO
!$OMP END DO NOWAIT
END IF

!-----------------------------------------------------------------------------
! Calculate the PAR absorbed by the top leaf and set top leaf N value.
!-----------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO m  = 1,veg_pts
  l = veg_index(m)
  apar(l) = (1.0 - omega(ft)) * ipar(l)
  acr(l)  = apar(l) / conpar
  IF ( l_trait_phys ) THEN
    IF ( ft > nnpft ) THEN
      lma_tmp = lma_from_prognostics(ft - nnpft, dvi_cpft(l,ft - nnpft))
      nleaf_top(l) = nmass(ft) * lma_tmp * 1000.0
        ! gN/gleaf * gleaf/m2 = gN/m2
    ELSE
      nleaf_top(l) = nmass(ft) * lma(ft) * 1000.0
        ! gN/gleaf * gleaf/m2 = gN/m2
    END IF
  ELSE
    nleaf_top(l) = nl0(ft)
  END IF
END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------------
! Calculate incoming PAR (mol photons m-2 s-1).
!-----------------------------------------------------------------------------
SELECT CASE ( can_rad_mod )
CASE ( 5:6 )
!$OMP DO SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    icr(l) = ipar(l) / conpar
  END DO
!$OMP END DO NOWAIT
END SELECT

!-----------------------------------------------------------------------------
! Calculate the LAI in each canopy layer.
!-----------------------------------------------------------------------------
SELECT CASE ( can_rad_mod )
CASE ( 4:7 )
!$OMP DO SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    dlai(l) = lai(l) / REAL(ilayers)

    ! Kmax of the layer - uniform split
    ! kmax_per_lyr(l) = kmax_pft(ft) / REAL(ilayers)
  END DO
!$OMP END DO NOWAIT
END SELECT


!-----------------------------------------------------------------------------
! Convert O3 concentration from ppb to moles.
!-----------------------------------------------------------------------------
IF ( l_o3_damage ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    o3mol(l) = o3(l) * pstar(l) / (rmol * tstar(l))
  END DO
!$OMP END DO NOWAIT
END IF

!-----------------------------------------------------------------------------
! Calculate partial pressure of oxygen.
!-----------------------------------------------------------------------------
IF ( c3(ft) == 1 ) THEN
  ! Calculate partial pressure of oxygen.
!$OMP DO SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    ! Convert O2 from mass fraction to partial pressure.
    oa(l)  = o2 / epo2 * pstar(l)
  END DO
!$OMP END DO
END IF  !  C3

!-----------------------------------------------------------------------------
! Convert CO2 from mass fraction to partial pressure.
! We ignore the (small) difference between the canopy- and
! reference-level CO2 concentrations.
!-----------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO m  = 1,veg_pts
  l = veg_index(m)
  ca(l) = co2c(l) / epco2 * pstar(l)
END DO
!$OMP END DO
!$OMP END PARALLEL

!-----------------------------------------------------------------------------
! Pre-calculate some expensive terms that do not change between iterations
! and layers.
! tstar is for the current pft, so can't be moved up another level.
!-----------------------------------------------------------------------------

SELECT CASE ( pft_photo_model )

CASE ( photo_collatz )
  ! Use the Collatz model (for C3 or C4 plants).
  CALL leaf_temp_responses( ft, land_pts, veg_pts, veg_index,                  &
                            pft_photo_model, tstar, oa, acr, actj, actv,       &
                            dsj, dsv, ccp, kc, ko, km, denom, qtenf_term,      &
                            vcmax_temp, jmax_temp, i2 )

CASE ( photo_farquhar )
  ! Use the Farquhar model (for C3 plants).
  ! Calculate a constant.
  sun_term = alpha_elec(ft) / conpar

  ! Load parameter values, depending on options.
  SELECT CASE ( photo_acclim_model )
  CASE ( 0 )
    ! No acclimation.
    ! Copy the PFT parameters, including fixed J:V.
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l,m)                   &
!$OMP SHARED(ds_jmax, ds_vcmax, dsj, dsv, ft, jv25, jv25_ratio,                &
!$OMP        actj, act_jmax, actv, act_vcmax, veg_index, veg_pts)              &
!$OMP SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      dsj(l)  = ds_jmax(ft)
      dsv(l)  = ds_vcmax(ft)
      jv25(l) = jv25_ratio(ft)
      actj(l) = act_jmax(ft)
      actv(l) = act_vcmax(ft)
    END DO
!$OMP END PARALLEL DO

  CASE ( photo_adapt, photo_acclim, photo_adapt_acclim )
    ! These use the same forms but t_growth_gb will generally be
    ! different. Although there is no dependency on PFT here (meaning
    ! this could be moved up and out of a PFT loop), we leave it here
    ! so that these parameters are calculated here regardless of the
    ! acclimation model selected.

    ! Decide whether the activation energies are subject to acclimation.  If
    ! they are, then the energies vary by gridbox but not by PFT, otherwise
    ! the energies vary by PFT but not by gridbox.
    SELECT CASE ( photo_act_model )
    CASE ( photo_act_pft )
      act_j_tmp(:) = [act_jmax(ft), 0.0, 0.0]
      act_v_tmp(:) = [act_vcmax(ft), 0.0, 0.0]
    CASE ( photo_act_gb )
      act_j_tmp(:) = act_j_coef(:)
      act_v_tmp(:) = act_v_coef(:)
    CASE DEFAULT
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
                   'photo_act_model should be photo_act_pft or photo_act_gb')
    END SELECT

!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l,m,th_degc,tg_degc)   &
!$OMP SHARED(dsj, dsj_coef, dsv, dsv_coef, jv25, jv25_coef,                    &
!$OMP        actj, act_j_tmp, actv, act_v_tmp,                                 &
!$OMP        t_home_gb, t_growth_gb, veg_index, veg_pts)                       &
!$OMP SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      th_degc = t_home_gb(l) - zerodegc
      tg_degc = t_growth_gb(l) - zerodegc
      dsj(l)  = dsj_coef(1) + dsj_coef(2) * th_degc + dsj_coef(3) * tg_degc
      dsv(l)  = dsv_coef(1) + dsv_coef(2) * th_degc + dsv_coef(3) * tg_degc
      jv25(l) = jv25_coef(1) + jv25_coef(2) * th_degc + jv25_coef(3) * tg_degc
      actj(l) = act_j_tmp(1) + act_j_tmp(2) * th_degc + act_j_tmp(3) * tg_degc
      actv(l) = act_v_tmp(1) + act_v_tmp(2) * th_degc + act_v_tmp(3) * tg_degc
    END DO
!$OMP END PARALLEL DO

  END SELECT  !  photo_acclim_model

  CALL leaf_temp_responses( ft, land_pts, veg_pts, veg_index,                  &
                            pft_photo_model, tstar, oa, acr, actj, actv,       &
                            dsj, dsv, ccp, kc, ko, km, denom, qtenf_term,      &
                            vcmax_temp, jmax_temp, i2,                  &
                            f_vc = f_vc )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'pft_photo_model should be photo_collatz or photo_farquhar')

END SELECT  !  pft_photo_model

!-----------------------------------------------------------------------------
! Soil water factor for the leaf fluxes. DESICA: no fsmc, the stomata close
! on psi_leaf through the Tuzet factor, from psi_leaf of the last timestep.
!-----------------------------------------------------------------------------
IF ( stomata_model == stomata_desica ) THEN
  CALL desica_fw( ft, land_pts, veg_pts, veg_index, psi_root_zone, fsmc_lim )
  fsmc(:) = 1.0
ELSE
  fsmc(:)     = fsmc_in(:)
  fsmc_lim(:) = fsmc_in(:)
END IF

!-----------------------------------------------------------------------------
! Leaf temperatures output: tstar unless l_leaf_temp (set in can_rad_mod 7).
!-----------------------------------------------------------------------------
IF ( .NOT. l_leaf_temp ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    t_leaf_sun(l) = tstar(l)
    t_leaf_shd(l) = tstar(l)
    t_can(l)      = tstar(l)
  END DO
END IF

!-----------------------------------------------------------------------------
! Calculate fluxes.
!-----------------------------------------------------------------------------

SELECT CASE ( can_rad_mod )

CASE ( 4 )

  !---------------------------------------------------------------------------
  ! Varying N model+altered leaf respiration
  ! Multiple canopy layers
  ! N varies through canopy as exponential function of layers.
  !---------------------------------------------------------------------------

  DO n = 1,ilayers
    !-------------------------------------------------------------------------
    ! Initialise GL and calculate the PAR absorbed in this layer.
    !-------------------------------------------------------------------------
    expkn = EXP( REAL(n-1) / REAL(ilayers) * (-kn(ft)) )
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                  &
!$OMP             SHARED(dlai, expkn, faparv, faparv_layer, gl, nleaf_layer, n,&
!$OMP                    nleaf_top, veg_index, veg_pts, kmax_per_lyr,          &
!$OMP                    kmax_pft, ft, ilayers)                                &
!$OMP             SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      gl(l)             = 0.0
      nleaf_layer(l)    = nleaf_top(l) * expkn
      faparv_layer(l,n) = faparv(l,n) * dlai(l)

      ! Kmax is assumed to co-vary with leaf nitrogen, so it is scaled using
      ! the same vertical decay profile as nleaf_layer.
      !
      ! Two equivalent formulations are used depending on canopy scheme:
      ! (i) CASE 4 uses pre-computed expkn
      ! (ii) CASE 5/6 uses explicit layer-normalised exponential form
      !
      ! Note: Kmax is intentionally NOT split by sunlit/shaded fractions (fsun),
      ! as this represents radiative partitioning rather than structural
      ! limitation.
      kmax_per_lyr(l) = (kmax_pft(ft) / REAL(ilayers)) * expkn

    END DO
!$OMP END PARALLEL DO

    !-------------------------------------------------------------------------
    ! Calculate photosynthetic parameters.
    !-------------------------------------------------------------------------
    CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,        &
                                veg_index, denom, jmax_temp, jv25,             &
                                nleaf_layer, qtenf_term, vcmax_temp,           &
                                jmax, rd_dark, vcmax )

    !-------------------------------------------------------------------------
    ! Iterate to ensure that the canopy humidity deficit is consistent with
    ! the H2O flux.
    !-------------------------------------------------------------------------

    DO k = 1,iter
      !-----------------------------------------------------------------------
      ! Diagnose the canopy-level humidity deficit.
      ! Initialise dark respiration with the uninhibited rate.
      !-----------------------------------------------------------------------
!$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                     &
!$OMP          SHARED(dq, dqc, gl, ra, ra_rc, rd, rd_dark, veg_index, veg_pts)
!$OMP DO SCHEDULE(STATIC)
      DO m = 1,veg_pts
        l = veg_index(m)
        ra_rc(l) = ra(l) * gl(l)
        dqc(l)   = dq(l) / (1.0 + ra_rc(l))
        rd(l)    = rd_dark(l)
      END DO
!$OMP END DO
!$OMP END PARALLEL

      !-----------------------------------------------------------------------
      ! Calculate the limiting factors for leaf photosynthesis
      !-----------------------------------------------------------------------
      CALL leaf_limits (ft, land_pts, pft_photo_model ,veg_pts, veg_index      &
,                       acr, apar, ca, ccp, dqc, fsmc_lim, je, kc, km, ko, oa  &
,                       pstar, vcmax                                           &
,                       clos_pts, open_pts, clos_index, open_index             &
,                       ci, wcarb, wexpt, wlite )

!$OMP PARALLEL DO IF(open_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                 &
!$OMP             SHARED(acr, apar, faparv, faparv_layer, ipar, land_index,    &
!$OMP                    n, open_index, open_pts, rd, t_i_length,              &
!$OMP                    wlite, wlitev, veg_index) SCHEDULE(STATIC)
      DO m = 1,open_pts
        l = veg_index(open_index(m))
        wlitev(l) = wlite(l) / apar(l) * faparv(l,n) * ipar(l)
        ! Calculate light inhibition of dark respiration.
        ! This does not change between iterations (though open_index might).
        IF (acr(l) * 1.0e6 * faparv_layer(l,n) >  10.0) THEN
          rd(l) = ( 0.5 - 0.05 * LOG(acr(l) * faparv_layer(l,n) * 1.0e6) )     &
                  * rd(l)
        END IF
      END DO
!$OMP END PARALLEL DO

      !-----------------------------------------------------------------------
      ! Calculate leaf-level fluxes.
      !-----------------------------------------------------------------------
      CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model, veg_pts    &
,                clos_index, open_index, veg_index                             &
,                ca, ci, fsmc, o3mol, ra, tstar, wcarb, wexpt, wlitev, rd      &
,                anetl, flux_o3_l, fo3_l, gl)

    END DO                 ! K-ITER

    !-------------------------------------------------------------------------
    ! Add to canopy-level values.
    !-------------------------------------------------------------------------
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                  &
!$OMP             SHARED(anetc, anetl, dlai, flux_o3, flux_o3_l, fo3, fo3_l,   &
!$OMP                    rdmean,ilayers,gc, gl, rd, rdc, veg_index, veg_pts,   &
!$OMP                    LAI,l_o3_damage                                     ) &
!$OMP             SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      anetc(l) = anetc(l) + anetl(l) * dlai(l)
      gc(l)    = gc(l)    + gl(l) * dlai(l)
      rdc(l)   = rdc(l)   + rd(l) * dlai(l)

      rdmean(l) = rdmean(l) + rd(l) / REAL(ilayers)

      IF (l_o3_damage) THEN
        flux_o3(l) = flux_o3(l) + flux_o3_l(l) * dlai(l)
        fo3(l)     = fo3(l) + fo3_l(l) * dlai(l) / lai(l)
      END IF
    END DO
!$OMP END PARALLEL DO

  END DO                   ! N LAYERS


CASE ( 5, 6 )
!---------------------------------------------------------------------------
! Sunlit and shaded leaves treated separately
! Multiple canopy layers
! N varies through canopy as exponential
!---------------------------------------------------------------------------

  ! ---------------------
  ! Stomatal optimisation
  !
  ! Each sunlit and shaded leaf in each layer draws from psi_src on its own
  ! parallel soil-to-leaf path, as in the big-leaf and two-leaf schemes, so
  ! every leaf is charged the hydraulic cost of the whole path. There is no
  ! shared canopy node: an earlier version solved the leaves from a shared
  ! node psi (outer bracketed iteration on the summed transpiration) and
  ! charged each leaf only its own leaf segment, som_leaf_resist_frac of the
  ! path resistance. That left the root-to-canopy drop uncharged, so each
  ! leaf saw about half the marginal hydraulic cost in moist soil, and E and
  ! gs came out high relative to big-leaf and two-leaf (see
  ! runs/roses/RED/notes/multilayer_profit_max_design.md).
  !
  ! ---------------------
  IF (leaf_flux_mod == leaf_flux_stom_opt) THEN

    fsmc_unity(:) = 1.0
    gl_max_lf(:)  = som_gl_max

        ! The humidity-deficit iteration (DO k = 1,iter) runs over the whole
        ! canopy, as in the big-leaf and two-leaf schemes: the deficit at the
        ! leaves is dq / (1 + ra * gc) with gc the canopy conductance (per m2
        ! ground) from the previous pass, the conductance the surface energy
        ! balance (sf_evap) sees. It used to iterate within each layer on that
        ! layer's per-leaf-area gl, so ra * gl ~ 0 and every leaf saw nearly
        ! the full above-canopy deficit. Same number of stom_opt_mod calls;
        ! only the last pass is kept.
        gc_ml_prev(:) = 0.0

        DO k = 1,iter

        el(:) = 0.0
        anetl(:) = 0.0
        anetc(:) = 0.0
        gc(:) = 0.0
        rdc(:) = 0.0
        psi_leaf(:) = 0.0
        leaf_k(:) = 0.0
        carbon_gain(:) = 0.0
        hydraulic_cost(:) = 0.0
        ci_gl_ml(:) = 0.0
        ci_lai_ml(:) = 0.0
        rdmean(:) = 0.0
        flux_o3(:) = 0.0
        IF ( l_o3_damage ) fo3(:) = 0.0

        DO n = 1,ilayers

            !-------------------------------------------------------------------
            ! Initialise GL for this layer.
            ! We could initialise to gl(n-1) here, but simpler to use zero
            ! and seems to converge pretty quickly anyway.
            !-------------------------------------------------------------------
        !$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE)                        &
        !$OMP             PRIVATE(l, m, kleaf_prof, kleaf_mean, kleaf_a)        &
        !$OMP             SHARED(ft, gl, ilayers, kn, knl, n, nleaf_top, nleaf_layer,  &
        !$OMP   veg_index,dlai,can_rad_mod, veg_pts, kmax_pft, kcrit,          &
        !$OMP   kmax_per_lyr, kcrit_per_lyr, kmax_leaf_lyr, kcrit_leaf_lyr,    &
        !$OMP   share_sup, lai) SCHEDULE(STATIC)
            DO m = 1,veg_pts
            l = veg_index(m)
            gl(l) = 0.0
            IF ( can_rad_mod == 6 ) THEN
                nleaf_layer(l) = nleaf_top(l) * EXP((n-1) * dlai(l) *          &
                                    (-knl(ft)))
            ELSE
                nleaf_layer(l) = nleaf_top(l) * EXP((n-1) /                    &
                                    REAL(ilayers) * (-kn(ft)))
            END IF

            ! Kmax is NOT scaled with leaf nitrogen (it used to follow the
            ! nleaf_layer decay profile): kmax_pft is the leaf-area-basis
            ! whole-plant conductance, the same in every layer.
            !
            ! Note: Kmax is intentionally NOT split by sunlit/shaded fractions
            ! (fsun), as this represents radiative partitioning rather than
            ! structural limitation.
            ! NOTE: no /REAL(ilayers) here, matching nleaf_layer above
            ! exactly (same decay profile assumption, same per-leaf-area
            ! basis). This used to divide by ilayers here despite
            ! nleaf_layer having no such factor - kmax_per_lyr is used as
            ! the hydraulic-capacity reference for el_sun/el_shd inside
            ! this layer's own stom_opt_mod call below, and those are
            ! per-leaf-area quantities (not yet dlai-integrated - see the
            ! "ground-area integrated fluxes" comment at their
            ! accumulation below), the same basis as vcmax/nleaf_layer.
            ! The spurious /ilayers made kmax_per_lyr ilayers-times too
            ! small relative to that, i.e. the per-layer hydraulic-cost
            ! check saw far less conductance capacity than the layer
            ! actually has, forcing stomata far more closed than warranted
            ! - very likely the dominant cause of multilayer's GPP running
            ! well below big-leaf's for the same site/forcing.
            ! kmax_pft is the leaf-area-basis whole-plant conductance,
            ! uniform through the canopy (see the big-leaf note below): the
            ! same value per unit leaf area in every layer, so the canopy
            ! total (sum of kmax_per_lyr * dlai) is kmax_pft * LAI. It used
            ! to decay with leaf N like nleaf_layer, which made kmax_pft a
            ! top-of-canopy value.
            kmax_per_lyr(l) = kmax_pft(ft)

            ! Scale kcrit by the same fraction of whole-plant kmax that this
            ! layer received, so the critical-conductance threshold stays
            ! calibrated at (1 - kcrit_fractional_loss) of THIS layer's own
            ! capacity (Brodribb & Cochard, 2009; Sabot et al., 2020).
            kcrit_per_lyr(l) = kmax_per_lyr(l) * (kcrit(ft) / kmax_pft(ft))

            ! This layer's whole soil-to-leaf path, on its own parallel path
            ! from psi_src (no shared canopy node). kcrit keeps the same
            ! kcrit/kmax fraction, so the kl > kcrit feasibility test is
            ! unchanged in terms of psi.
            ! The path conductance follows the leaf-N profile (sun leaves
            ! have higher leaf hydraulic conductance than shade leaves and it
            ! co-varies with photosynthetic capacity; Sack & Holbrook 2006,
            ! Brodribb et al. 2007), normalised to a canopy mean of 1 so the
            ! canopy total stays kmax_pft * LAI: conductance is
            ! redistributed from shaded to sunlit layers, not changed in
            ! total. This is the same convention as the two-leaf scheme,
            ! which splits kmax_pft * LAI by N-weighted leaf area.
            ! kleaf_mean is the closed-form mean of the same discrete
            ! profile as nleaf_layer over the ilayers layers (equal dlai).
            IF ( can_rad_mod == 6 ) THEN
                kleaf_a = knl(ft) * dlai(l)
            ELSE
                kleaf_a = kn(ft) / REAL(ilayers)
            END IF
            kleaf_prof = EXP(-REAL(n-1) * kleaf_a)
            IF ( kleaf_a > 1.0e-6 ) THEN
                kleaf_mean = (1.0 - EXP(-REAL(ilayers) * kleaf_a)) /          &
                             (REAL(ilayers) * (1.0 - EXP(-kleaf_a)))
            ELSE
                kleaf_mean = 1.0
            END IF
            kmax_leaf_lyr(l)  = kmax_per_lyr(l)  * (kleaf_prof / kleaf_mean)
            kcrit_leaf_lyr(l) = kcrit_per_lyr(l) * (kleaf_prof / kleaf_mean)

            ! Soil supply cap per unit leaf area (el is per leaf area in
            ! each layer): each leaf may use its kmax share of the supply,
            ! as in the two-leaf scheme. Summed over the canopy (dlai
            ! weights) the shares add to 1. Sunlit and shaded leaves of a
            ! layer have the same kmax per leaf area, so the same share.
            share_sup(l) = (kleaf_prof / kleaf_mean)                         &
                           / MAX(lai(l), EPSILON(1.0))

            END DO
        !$OMP END PARALLEL DO

            !-------------------------------------------------------------------
            ! Calculate photosynthetic parameters.
            !-------------------------------------------------------------------
            CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,&
                                        veg_index, denom, jmax_temp, jv25,     &
                                        nleaf_layer, qtenf_term, vcmax_temp,   &
                                        jmax, rd_dark, vcmax )

            IF ( pft_photo_model == photo_farquhar ) THEN
            !-------------------------------------------------------------------
            ! Calculate sunlit and shaded radiation terms.
            !-------------------------------------------------------------------
        !$OMP PARALLEL DO IF(veg_pts > 1)                                              &
        !$OMP DEFAULT(NONE)                                                            &
        !$OMP SHARED(n, veg_pts, veg_index, fapar_sun, fapar_shd, i2_shd, i2_sun,      &
        !$OMP        ipar, sun_term)                                                   &
        !$OMP PRIVATE(l,m)                                                             &
        !$OMP SCHEDULE(STATIC)
            DO m = 1,veg_pts
                l = veg_index(m)
                i2_sun(l) = sun_term * fapar_sun(l,n) * ipar(l)
                i2_shd(l) = sun_term * fapar_shd(l,n) * ipar(l)
            END DO
        !$OMP END PARALLEL DO

            !-------------------------------------------------------------------
            ! Calculate the electron fluxes.
            ! Although we ultimately only need je_sun and je_shd, we also calculate
            ! je because we use that to calculate a single value of wlite, from which
            ! sunlit and shaded terms are calculated using je_sun and je_shd.
            !-----------------------------------------------------------------------
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je)
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_sun,     &
                                     jmax, je_sun)
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_shd,     &
                                     jmax, je_shd)

            ! Calculate ratios of electron flux terms.
            ! These are used to scale the light-limited photosynthesis, which is
            ! linearly related to the electron flux.
        !$OMP PARALLEL DO IF(veg_pts > 1)                                              &
        !$OMP DEFAULT(NONE)                                                            &
        !$OMP SHARED(veg_pts, veg_index, je, je_shd, je_shd_ratio, je_sun,             &
        !$OMP        je_sun_ratio)                                                     &
        !$OMP PRIVATE(l,m)                                                             &
        !$OMP SCHEDULE(STATIC)
            DO m = 1,veg_pts
                l = veg_index(m)
                IF ( je(l) > TINY(je(l)) ) THEN
                je_sun_ratio(l) = je_sun(l) / je(l)
                je_shd_ratio(l) = je_shd(l) / je(l)
                END IF
            END DO
        !$OMP END PARALLEL DO

            END IF  !  pft_photo_model

            !-----------------------------------------------------------------------
            ! Diagnose the canopy-level humidity deficit from the canopy
            ! conductance of the previous pass (see DO k above).
            ! Initialise the sunlit and shaded respiration rates with the
            ! uninhibited, sunlit respiration rate.
            !-----------------------------------------------------------------------
        !$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                     &
        !$OMP          SHARED(dq, dqc, gc_ml_prev, ra, ra_rc, rd_dark, rd_shd, rd_sun, &
        !$OMP                 veg_index, veg_pts)
        !$OMP DO SCHEDULE(STATIC)
            DO m = 1,veg_pts
                l = veg_index(m)
                ra_rc(l)  = ra(l) * gc_ml_prev(l)
                dqc(l)    = dq(l) / (1.0 + ra_rc(l))
                rd_sun(l) = rd_dark(l)
                rd_shd(l) = rd_dark(l)
            END DO
        !$OMP END DO
        !$OMP END PARALLEL

            !-----------------------------------------------------------------------
            ! Calculate the limiting factors for leaf photosynthesis
            !-----------------------------------------------------------------------
            ! fsmc_unity, not fsmc: leaf_limits closes stomata wherever
            ! fsmc == 0 (soil moisture at/below wilting), which is the
            ! standard JULES empirical soil-moisture stress. The stomatal
            ! optimisation handles water stress itself through psi_root_zone
            ! and the hydraulic cost, and must not also have JULES' fsmc
            ! switch shut it off - otherwise every layer is forced closed
            ! whenever fsmc hits 0, regardless of what the hydraulics say
            ! (this is what collapsed multilayer GPP to ~0 through the
            ! FR-Pue summer 2003 drought while big-leaf, which classifies
            ! open/closed on apar alone - see the leaf_flux_stom_opt branch
            ! of the big-leaf path - kept photosynthesising). With fsmc_unity
            ! the Medlyn closure test reduces to apar == 0, the same as
            ! big-leaf. (The Jacobs model's dq >= dqcrit closure still
            ! applies here.)
            CALL leaf_limits (ft, land_pts, pft_photo_model, veg_pts           &
        ,                     veg_index, acr, apar, ca, ccp, dqc, fsmc_unity   &
        ,                     je, kc                                           &
        ,                     km, ko, oa, pstar, vcmax                         &
        ,                     clos_pts, open_pts, clos_index, open_index       &
        ,                     ci, wcarb, wexpt, wlite)

        !$OMP PARALLEL IF(veg_pts > 1)                                                 &
        !$OMP DEFAULT(NONE)                                                            &
        !$OMP PRIVATE(m,l)                                                             &
        !$OMP SHARED(veg_pts,veg_index,rd_sun,rd_shd,rd,open_pts,pft_photo_model,      &
        !$OMP        open_index,wlitev_sun,wlite,apar,fapar_sun,ipar,je_shd_ratio,     &
        !$OMP        je_sun_ratio,wlitev_shd,fapar_shd,icr,n,fsun,dlai)

            IF ( pft_photo_model == photo_collatz ) THEN
        !$OMP DO SCHEDULE(STATIC)
                DO m = 1,open_pts
                l = veg_index(open_index(m))
                wlitev_sun(l) = wlite(l) / apar(l) * fapar_sun(l,n) * ipar(l)
                wlitev_shd(l) = wlite(l) / apar(l) * fapar_shd(l,n) * ipar(l)

                END DO
        !$OMP END DO NOWAIT
            ELSE
        !$OMP DO SCHEDULE(STATIC)
                DO m = 1,open_pts
                l = veg_index(open_index(m))
                wlitev_sun(l) = wlite(l) * je_sun_ratio(l)
                wlitev_shd(l) = wlite(l) * je_shd_ratio(l)


                ! Kmax of the layer for sunlit and shaded
                !dKmax_sun(l,n) = kmax_per_lyr(l) * fsun(l,n)
                !dKmax_shd(l,n) = kmax_per_lyr(l) * (1.0 - fsun(l,n))

                END DO

        !$OMP END DO NOWAIT
            END IF  !  photo_model


            !-----------------------------------------------------------------------
            ! Introducing inhibition of leaf respiration in the light for sunlit
            ! and shaded leaves, from papers by Atkin et al. This is an
            ! improvement over the description used for can_rad_mod=4.
            ! This does not change between iterations (though open_index might).
            !-----------------------------------------------------------------------
        !$OMP DO SCHEDULE(STATIC)
            DO m = 1,open_pts
                l = veg_index(open_index(m))
                IF ( fapar_sun(l,n) * icr(l) * fsun(l,n) * dlai(l) *           &
                    1.0e6 >  10.0 ) rd_sun(l) = 0.7 * rd_sun(l)
                IF ( fapar_shd(l,n) * icr(l) * (1.0 - fsun(l,n)) * dlai(l) *   &
                    1.0e6 >  10.0 ) rd_shd(l) = 0.7 * rd_shd(l)
            END DO
        !$OMP END DO
        !$OMP END PARALLEL

            ! Added CASE block to do gs_opt on multi-layer canopy, 29 Apr, MGDK
            SELECT CASE ( leaf_flux_mod)
                CASE (leaf_flux_fsmc)
                !-----------------------------------------------------------------------
                ! Calculate leaf-level fluxes separately for sunlit and shaded leaves.
                !-----------------------------------------------------------------------
                CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model   &
            ,              veg_pts, clos_index, open_index, veg_index          &
            ,              ca, ci, fsmc, o3mol, ra, tstar                      &
            ,              wcarb, wexpt, wlitev_sun, rd_sun                    &
            ,              anetl_sun, flux_o3_l_sun, fo3_l_sun, gl_sun)

                CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model   &
            ,              veg_pts, clos_index, open_index, veg_index          &
            ,              ca, ci, fsmc, o3mol, ra, tstar                      &
            ,              wcarb, wexpt, wlitev_shd, rd_shd                    &
            ,              anetl_shd, flux_o3_l_shd, fo3_l_shd, gl_shd)

            CASE (leaf_flux_stom_opt)

                l_multilayer = .TRUE.
                ! share_sup: this layer's kmax share, set with kmax_leaf_lyr.
                CALL apply_supply_limit( land_pts, veg_pts, veg_index,         &
                                         e_supply, share_sup, dqc, tstar,      &
                                         pstar, gl_max_lf, gl_max_eff )
                CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
                ! Added gs opt call for sunlit, 29 Apr, MGDK
                CALL stom_opt_mod (                                            &
                    ! IN
                    land_pts, som_base_parm, ft, open_pts, open_index,         &
                    pft_photo_model, veg_pts, veg_index,                       &
                    ca, psi_src, acr, apar, oa, vcmax, kc, ko, ccp, pstar,     &
                    km, dqc, qs, je, tstar, je_sun_ratio, fapar_sun(:,n),      &
                    kmax_leaf_lyr, kcrit_leaf_lyr, gl_max_eff, ipar,           &
                    l_multilayer,                                              &
                    ! IN OUT
                    rd_sun,                                                    &
                    ! OUT
                    ci_sun_ml, anetl_sun, el_sun, flux_o3_l_sun, fo3_l_sun,    &
                    gl_sun, psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun           &
                )

                ! Added gs opt call for shaded , 29 Apr, MGDK
                CALL stom_opt_mod (                                            &
                    ! IN
                    land_pts, som_base_parm, ft, open_pts, open_index,         &
                    pft_photo_model, veg_pts, veg_index,                       &
                    ca, psi_src, acr, apar, oa, vcmax, kc, ko, ccp, pstar,     &
                    km, dqc, qs, je, tstar,  je_shd_ratio, fapar_shd(:,n),     &
                    kmax_leaf_lyr, kcrit_leaf_lyr, gl_max_eff, ipar,           &
                    l_multilayer,                                              &
                    ! IN OUT
                    rd_shd,                                                    &
                    ! OUT
                    ci_shd_ml, anetl_shd, el_shd, flux_o3_l_shd, fo3_l_shd,    &
                    gl_shd, psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd           &
                )

                ! Closed leaves get min_gl_pft from stom_opt_mod, a canopy
                ! value: divide by LAI so the canopy total (sum of gl * dlai
                ! over the layers) is min_gl_pft, as in big leaf and two-leaf
                ! (review #11).
                DO i = 1,clos_pts
                    l = veg_index(clos_index(i))
                    IF ( lai(l) > EPSILON(0.0) ) THEN
                        gl_sun(l) = gl_sun(l) / lai(l)
                        gl_shd(l) = gl_shd(l) / lai(l)
                    END IF
                END DO

            CASE DEFAULT
                errcode = 101  !  a hard error
                CALL ereport(RoutineName, errcode,                                     &
                    'leaf_flux_mod should be 1 fsmc regulated or 2 stomatal optimisation')
            END SELECT ! leaf_flux_model
            !-----------------------------------------------------------------------
            ! Update layer conductance.
            !-----------------------------------------------------------------------
        !$OMP PARALLEL DO                                                              &
        !$OMP SCHEDULE(STATIC)                                                         &
        !$OMP DEFAULT(NONE)                                                            &
        !$OMP PRIVATE(m,l)                                                             &
        !$OMP SHARED(veg_pts,veg_index,gl,fsun,n,gl_sun,gl_shd,rd,rd_sun,rd_shd)
            DO m = 1,veg_pts
                l = veg_index(m)
                gl(l) = fsun(l,n) * gl_sun(l) + (1.0 - fsun(l,n)) * gl_shd(l)
                rd(l) = fsun(l,n) * rd_sun(l) + (1.0 - fsun(l,n)) * rd_shd(l)
            END DO
        !$OMP END PARALLEL DO

            !-------------------------------------------------------------------------
            ! Add to canopy-level values.
            !-------------------------------------------------------------------------
        !$OMP PARALLEL DO                                                      &
        !$OMP SCHEDULE(STATIC)                                                 &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP PRIVATE(m,l)                                                     &
        !$OMP SHARED(veg_pts,veg_index,anetl,fsun,anetl_sun,anetl_shd,anetc,dlai,gc,   &
        !$OMP        gl,rdc,rd,rdmean,ilayers,l_o3_damage,flux_o3_l,flux_o3_l_sun,     &
        !$OMP        flux_o3_l_shd,fo3_l_sun,fo3_l_shd,flux_o3,fo3,lai,n,fo3_l,        &
        !$OMP        el,el_sun,el_shd,hydraulic_cost,HC_sun,HC_shd,carbon_gain,        &
        !$OMP        CG_sun,CG_shd,psi_leaf,psi_leaf_sun,psi_leaf_shd,leaf_k,          &
        !$OMP        leaf_k_sun,leaf_k_shd,ci_gl_ml,ci_lai_ml,gl_sun,gl_shd,           &
        !$OMP        ci_sun_ml,ci_shd_ml)
            DO m = 1,veg_pts
            l = veg_index(m)

            anetl(l)     = fsun(l,n) * anetl_sun(l)                                  &
                            + (1.0 - fsun(l,n)) * anetl_shd(l)
            anetc(l)     = anetc(l) + anetl(l) * dlai(l)

            gc(l)        = gc(l)  + gl(l) * dlai(l)
            rdc(l)       = rdc(l) + rd(l) * dlai(l)

            rdmean(l)    = rdmean(l) + rd(l) / REAL(ilayers)

            ! These are ground-area integrated fluxes, mdgk 9th June 20226
            el(l) = el(l) +                                                    &
                        (fsun(l,n) * el_sun(l) +                               &
                        (1.0 - fsun(l,n)) * el_shd(l)) *                       &
                        dlai(l)



            hydraulic_cost(l) = hydraulic_cost(l) +                            &
                                    (fsun(l,n) * HC_sun(l) +                   &
                                    (1.0 - fsun(l,n)) * HC_shd(l)) *           &
                                    dlai(l)

            carbon_gain(l) = carbon_gain(l) +                                  &
                                    (fsun(l,n) * CG_sun(l) +                   &
                                    (1.0 - fsun(l,n)) * CG_shd(l)) *           &
                                    dlai(l)

            ! Calculate LAI weighted integral over the canopy
            ! This is turned into a canopy weighted mean below...
            ! Probably better called sum or something until the mean is calculated
            psi_leaf(l) = psi_leaf(l) +                                        &
                            (fsun(l,n) * psi_leaf_sun(l) +                     &
                            (1.0 - fsun(l,n)) * psi_leaf_shd(l)) *             &
                            dlai(l)

            leaf_k(l) = leaf_k(l) + (fsun(l,n) * leaf_k_sun(l) +               &
                                    (1.0 - fsun(l,n)) * leaf_k_shd(l)) *       &
                                    dlai(l)

            ci_gl_ml(l) = ci_gl_ml(l) +                                        &
                          (fsun(l,n) * gl_sun(l) * ci_sun_ml(l) +              &
                          (1.0 - fsun(l,n)) * gl_shd(l) * ci_shd_ml(l)) *      &
                          dlai(l)
            ci_lai_ml(l) = ci_lai_ml(l) +                                      &
                           (fsun(l,n) * ci_sun_ml(l) +                         &
                           (1.0 - fsun(l,n)) * ci_shd_ml(l)) * dlai(l)





            IF (l_o3_damage) THEN
                flux_o3_l(l) = fsun(l,n) * flux_o3_l_sun(l)                    &
                            + (1.0 - fsun(l,n)) * flux_o3_l_shd(l)
                fo3_l(l)     = fsun(l,n) * fo3_l_sun(l)                        &
                            + (1.0 - fsun(l,n)) * fo3_l_shd(l)

                flux_o3(l)   = flux_o3(l) + flux_o3_l(l) * dlai(l)
                fo3(l)       = fo3(l)     + fo3_l(l) * dlai(l) / lai(l)
            END IF
            END DO
        !$OMP END PARALLEL DO

        END DO                   ! N LAYERS

        ! Canopy conductance of this pass, for the next pass's deficit.
        gc_ml_prev(:) = gc(:)

        END DO                   ! K-ITER (whole canopy)

        ! Calculate canopy mean properties
        DO m = 1, veg_pts
            l = veg_index(m)

            ! LAI-weighted canopy means. Carbon gain and hydraulic cost are
            ! per-leaf 0-1 quantities, as in the big-leaf and two-leaf
            ! schemes; their dlai-weighted sums ranged 0..LAI.
            IF (lai(l) > 0.0) THEN
                psi_leaf(l)       = psi_leaf(l)       / lai(l)
                leaf_k(l)         = leaf_k(l)         / lai(l)
                carbon_gain(l)    = carbon_gain(l)    / lai(l)
                hydraulic_cost(l) = hydraulic_cost(l) / lai(l)
            END IF

            ! gl-weighted canopy ci (as two-leaf); the LAI mean when no
            ! stomata are open (closed leaves have ci = ca).
            IF ( gc(l) > TINY(gc(l)) ) THEN
                ci(l) = ci_gl_ml(l) / gc(l)
            ELSE IF ( lai(l) > 0.0 ) THEN
                ci(l) = ci_lai_ml(l) / lai(l)
            ELSE
                ci(l) = ca(l)
            END IF

        END DO
ELSE

    !-----------------------------------------------
    ! No stomatal optimisation
    !-----------------------------------------------

    DO n = 1,ilayers

      !-------------------------------------------------------------------------
      ! Initialise GL for this layer.
      ! We could initialise to gl(n-1) here, but simpler to use zero
      ! and seems to converge pretty quickly anyway.
      !-------------------------------------------------------------------------
  !$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                  &
  !$OMP             SHARED(ft, gl, ilayers, kn, knl, n, nleaf_top, nleaf_layer,  &
  !$OMP   veg_index,dlai,can_rad_mod, veg_pts) SCHEDULE(STATIC)
      DO m = 1,veg_pts
      l = veg_index(m)
      gl(l) = 0.0
      IF ( can_rad_mod == 6 ) THEN
          nleaf_layer(l) = nleaf_top(l) * EXP((n-1) * dlai(l) * (-knl(ft)))
      ELSE
          nleaf_layer(l) = nleaf_top(l) * EXP((n-1) / REAL(ilayers) * (-kn(ft)))
      END IF

      END DO
  !$OMP END PARALLEL DO

      !-------------------------------------------------------------------------
      ! Calculate photosynthetic parameters.
      !-------------------------------------------------------------------------
      CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,        &
                                  veg_index, denom, jmax_temp, jv25,             &
                                  nleaf_layer, qtenf_term, vcmax_temp,           &
                                  jmax, rd_dark, vcmax )

      IF ( pft_photo_model == photo_farquhar ) THEN
      !-----------------------------------------------------------------------
      ! Calculate sunlit and shaded radiation terms.
      !-----------------------------------------------------------------------
  !$OMP PARALLEL DO IF(veg_pts > 1)                                              &
  !$OMP DEFAULT(NONE)                                                            &
  !$OMP SHARED(n, veg_pts, veg_index, fapar_sun, fapar_shd, i2_shd, i2_sun,      &
  !$OMP        ipar, sun_term)                                                   &
  !$OMP PRIVATE(l,m)                                                             &
  !$OMP SCHEDULE(STATIC)
      DO m = 1,veg_pts
          l = veg_index(m)
          i2_sun(l) = sun_term * fapar_sun(l,n) * ipar(l)
          i2_shd(l) = sun_term * fapar_shd(l,n) * ipar(l)
      END DO
  !$OMP END PARALLEL DO

      !-----------------------------------------------------------------------
      ! Calculate the electron fluxes.
      ! Although we ultimately only need je_sun and je_shd, we also calculate
      ! je because we use that to calculate a single value of wlite, from which
      ! sunlit and shaded terms are calculated using je_sun and je_shd.
      !-----------------------------------------------------------------------
      CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je)
      CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_sun, jmax,     &
                              je_sun)
      CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_shd, jmax,     &
                              je_shd)

      ! Calculate ratios of electron flux terms.
      ! These are used to scale the light-limited photosynthesis, which is
      ! linearly related to the electron flux.
  !$OMP PARALLEL DO IF(veg_pts > 1)                                              &
  !$OMP DEFAULT(NONE)                                                            &
  !$OMP SHARED(veg_pts, veg_index, je, je_shd, je_shd_ratio, je_sun,             &
  !$OMP        je_sun_ratio)                                                     &
  !$OMP PRIVATE(l,m)                                                             &
  !$OMP SCHEDULE(STATIC)
      DO m = 1,veg_pts
          l = veg_index(m)
          IF ( je(l) > TINY(je(l)) ) THEN
          je_sun_ratio(l) = je_sun(l) / je(l)
          je_shd_ratio(l) = je_shd(l) / je(l)
          END IF
      END DO
  !$OMP END PARALLEL DO

      END IF  !  pft_photo_model

      !-------------------------------------------------------------------------
      ! Iterate to ensure that the canopy humidity deficit is consistent with
      ! the H2O flux.
      !-------------------------------------------------------------------------

      DO k = 1,iter

      !-----------------------------------------------------------------------
      ! Diagnose the canopy-level humidity deficit.
      ! Initialise the sunlit and shaded respiration rates with the
      ! uninhibited, sunlit respiration rate.
      !-----------------------------------------------------------------------
  !$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                     &
  !$OMP          SHARED(dq, dqc, gl, ra, ra_rc, rd_dark, rd_shd, rd_sun,         &
  !$OMP                 veg_index, veg_pts)
  !$OMP DO SCHEDULE(STATIC)
      DO m = 1,veg_pts
          l = veg_index(m)
          ra_rc(l)  = ra(l) * gl(l)
          dqc(l)    = dq(l) / (1.0 + ra_rc(l))
          rd_sun(l) = rd_dark(l)
          rd_shd(l) = rd_dark(l)
      END DO
  !$OMP END DO
  !$OMP END PARALLEL

      !-----------------------------------------------------------------------
      ! Calculate the limiting factors for leaf photosynthesis
      !-----------------------------------------------------------------------
      CALL leaf_limits (ft, land_pts, pft_photo_model, veg_pts, veg_index      &
  ,                       acr, apar, ca, ccp, dqc, fsmc_lim, je, kc, km, ko, oa  &
  ,                       pstar, vcmax                                           &
  ,                       clos_pts, open_pts, clos_index, open_index             &
  ,                       ci, wcarb, wexpt, wlite)

  !$OMP PARALLEL IF(veg_pts > 1)                                                 &
  !$OMP DEFAULT(NONE)                                                            &
  !$OMP PRIVATE(m,l)                                                             &
  !$OMP SHARED(veg_pts,veg_index,rd_sun,rd_shd,rd,open_pts,pft_photo_model,      &
  !$OMP        open_index,wlitev_sun,wlite,apar,fapar_sun,ipar,je_shd_ratio,     &
  !$OMP        je_sun_ratio,wlitev_shd,fapar_shd,icr,n,fsun,dlai)

      IF ( pft_photo_model == photo_collatz ) THEN
  !$OMP DO SCHEDULE(STATIC)
          DO m = 1,open_pts
          l = veg_index(open_index(m))
          wlitev_sun(l) = wlite(l) / apar(l) * fapar_sun(l,n) * ipar(l)
          wlitev_shd(l) = wlite(l) / apar(l) * fapar_shd(l,n) * ipar(l)


          END DO
  !$OMP END DO NOWAIT
      ELSE
  !$OMP DO SCHEDULE(STATIC)
          DO m = 1,open_pts
          l = veg_index(open_index(m))
          wlitev_sun(l) = wlite(l) * je_sun_ratio(l)
          wlitev_shd(l) = wlite(l) * je_shd_ratio(l)


          END DO

  !$OMP END DO NOWAIT
      END IF  !  photo_model


      !-----------------------------------------------------------------------
      ! Introducing inhibition of leaf respiration in the light for sunlit
      ! and shaded leaves, from papers by Atkin et al. This is an
      ! improvement over the description used for can_rad_mod=4.
      ! This does not change between iterations (though open_index might).
      !-----------------------------------------------------------------------
  !$OMP DO SCHEDULE(STATIC)
      DO m = 1,open_pts
          l = veg_index(open_index(m))
          IF ( fapar_sun(l,n) * icr(l) * fsun(l,n) * dlai(l) *                   &
              1.0e6 >  10.0 ) rd_sun(l) = 0.7 * rd_sun(l)
          IF ( fapar_shd(l,n) * icr(l) * (1.0 - fsun(l,n)) * dlai(l) *           &
              1.0e6 >  10.0 ) rd_shd(l) = 0.7 * rd_shd(l)
      END DO
  !$OMP END DO
  !$OMP END PARALLEL

      ! Added CASE block to do gs_opt on multi-layer canopy, 29 Apr, MGDK
      SELECT CASE ( leaf_flux_mod)
          CASE (leaf_flux_fsmc)
          !-----------------------------------------------------------------------
          ! Calculate leaf-level fluxes separately for sunlit and shaded leaves.
          !-----------------------------------------------------------------------
          CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model, veg_pts &
      ,                clos_index, open_index, veg_index                         &
      ,                ca, ci, fsmc, o3mol, ra, tstar                            &
      ,                wcarb, wexpt, wlitev_sun, rd_sun                          &
      ,                anetl_sun, flux_o3_l_sun, fo3_l_sun, gl_sun)

          CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model, veg_pts &
      ,                clos_index, open_index, veg_index                         &
      ,                ca, ci, fsmc, o3mol, ra, tstar                            &
      ,                wcarb, wexpt, wlitev_shd, rd_shd                          &
      ,                anetl_shd, flux_o3_l_shd, fo3_l_shd, gl_shd)



      CASE DEFAULT
          errcode = 101  !  a hard error
          CALL ereport(RoutineName, errcode,                                     &
              'leaf_flux_mod should be 1 fsmc regulated or 2 stomatal optimisation')
      END SELECT ! leaf_flux_model
      !-----------------------------------------------------------------------
      ! Update layer conductance.
      !-----------------------------------------------------------------------
  !$OMP PARALLEL DO                                                              &
  !$OMP SCHEDULE(STATIC)                                                         &
  !$OMP DEFAULT(NONE)                                                            &
  !$OMP PRIVATE(m,l)                                                             &
  !$OMP SHARED(veg_pts,veg_index,gl,fsun,n,gl_sun,gl_shd,rd,rd_sun,rd_shd)
      DO m = 1,veg_pts
          l = veg_index(m)
          gl(l) = fsun(l,n) * gl_sun(l) + (1.0 - fsun(l,n)) * gl_shd(l)
          rd(l) = fsun(l,n) * rd_sun(l) + (1.0 - fsun(l,n)) * rd_shd(l)
      END DO
  !$OMP END PARALLEL DO

      END DO                 ! K-ITER

      !-------------------------------------------------------------------------
      ! Add to canopy-level values.
      !-------------------------------------------------------------------------
  !$OMP PARALLEL DO                                                              &
  !$OMP SCHEDULE(STATIC)                                                         &
  !$OMP DEFAULT(NONE)                                                            &
  !$OMP PRIVATE(m,l)                                                             &
  !$OMP SHARED(veg_pts,veg_index,anetl,fsun,anetl_sun,anetl_shd,anetc,dlai,gc,   &
  !$OMP        gl,rdc,rd,rdmean,ilayers,l_o3_damage,flux_o3_l,flux_o3_l_sun,     &
  !$OMP        flux_o3_l_shd,fo3_l_sun,fo3_l_shd,flux_o3,fo3,lai,n,fo3_l)
      DO m = 1,veg_pts
      l = veg_index(m)

      anetl(l)     = fsun(l,n) * anetl_sun(l)                                  &
                      + (1.0 - fsun(l,n)) * anetl_shd(l)
      anetc(l)     = anetc(l) + anetl(l) * dlai(l)

      gc(l)        = gc(l)  + gl(l) * dlai(l)
      rdc(l)       = rdc(l) + rd(l) * dlai(l)

      rdmean(l)    = rdmean(l) + rd(l) / REAL(ilayers)




      IF (l_o3_damage) THEN
          flux_o3_l(l) = fsun(l,n) * flux_o3_l_sun(l)                            &
                      + (1.0 - fsun(l,n)) * flux_o3_l_shd(l)
          fo3_l(l)     = fsun(l,n) * fo3_l_sun(l)                                &
                      + (1.0 - fsun(l,n)) * fo3_l_shd(l)

          flux_o3(l)   = flux_o3(l) + flux_o3_l(l) * dlai(l)
          fo3(l)       = fo3(l)     + fo3_l(l) * dlai(l) / lai(l)
      END IF
      END DO
  !$OMP END PARALLEL DO

  END DO                   ! N LAYERS


END IF

CASE ( 1 )

  !---------------------------------------------------------------------------
  ! "Big leaf" model.
  ! N varies through canopy according to Beers Law
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Calculate photosynthetic parameters.
  ! There is no light limitation of dark respiration in this case.
  !---------------------------------------------------------------------------
  CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,          &
                              veg_index, denom, jmax_temp, jv25,               &
                              nleaf_top, qtenf_term, vcmax_temp,               &
                              jmax, rd, vcmax )

  IF ( pft_photo_model == photo_farquhar ) THEN
    ! Calculate the electron flux.
    CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je)
  END IF

  !---------------------------------------------------------------------------
  ! DESICA: solve fw = Tuzet(psi_leaf at the end of the step, for the E that
  ! fw gives) by bisection on fw in [0, 1], re-running the big-leaf fluxes
  ! below for each trial fw (see desica_jls_mod). Other models: one pass.
  !---------------------------------------------------------------------------
  n_pass = 1
  IF ( stomata_model == stomata_desica ) THEN
    n_pass = n_fw_bisect + 1
    fw_lo(:) = 0.0
    fw_hi(:) = 1.0
    ksr_ds(:) = -1.0
    IF ( l_som_rhizo_series ) ksr_ds(:) = ksr(:)
  END IF

  DO i_pass = 1,n_pass

  IF ( stomata_model == stomata_desica ) THEN
    fsmc_lim(:) = 0.5 * (fw_lo(:) + fw_hi(:))
  END IF

  !---------------------------------------------------------------------------
  ! Iterate to ensure that the canopy humidity deficit is consistent with the
  ! H2O flux. The first estimate of the canopy humidity deficit uses the
  ! stomatal conductance from the previous timestep.
  !---------------------------------------------------------------------------
  DO k = 1,iter

    !-------------------------------------------------------------------------
    ! Diagnose the canopy-level humidity deficit.
    !-------------------------------------------------------------------------
    DO m = 1,veg_pts
      l = veg_index(m)
      ra_rc(l) = ra(l) * gc(l)
      dqc(l)   = dq(l) / (1.0 + ra_rc(l))
    END DO

    SELECT CASE ( leaf_flux_mod)
      CASE (leaf_flux_fsmc)
        IF (stomata_model == stomata_sox) THEN
          !-------------------------------------------------------------------
          ! SOX (Eller et al. 2020, semi-analytical; trunk vn7.9).
          !-------------------------------------------------------------------
          CALL leaf_processes_sox(ft, land_pts, veg_pts, veg_index             &
,                                 acr, apar, ca, ccp, dqc, kc, ko, oa          &
,                                 pstar, vcmax, tstar, ht, psi_root_zone, rd   &
,                                 clos_pts, open_pts, clos_index, open_index   &
,                                 ci, o3mol, ra, anetl, flux_o3_l              &
,                                 fo3_l, gl, lwp_c, pft_photo_model)

        ELSE
          !-------------------------------------------------------------------
          ! Calculate the limiting factors for leaf photosynthesis.
          !-------------------------------------------------------------------
          CALL leaf_limits (ft, land_pts, pft_photo_model, veg_pts, veg_index  &
          ,                 acr, apar, ca, ccp, dqc, fsmc_lim, je, kc, km, ko  &
          ,                 oa, pstar, vcmax                                   &
          ,                 clos_pts, open_pts, clos_index, open_index         &
          ,                 ci, wcarb, wexpt, wlite)

          !-------------------------------------------------------------------
          ! Calculate leaf-level fluxes.
          !-------------------------------------------------------------------
          CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model         &
          ,          veg_pts, clos_index, open_index, veg_index                &
          ,          ca, ci, fsmc, o3mol, ra, tstar, wcarb, wexpt, wlite, rd   &
          ,          anetl, flux_o3_l, fo3_l, gl)
        END IF ! stomata_model

        !-------------------------------------------------------------------------
        ! Scale to canopy level.
        !-------------------------------------------------------------------------
        DO m = 1,veg_pts
          l = veg_index(m)
          gc(l) = fpar(l) * gl(l)
        END DO

      CASE (leaf_flux_stom_opt)

        ! Identify locations with no par and hence closed stomata.
        clos_pts = 0
        open_pts = 0
        DO i = 1,veg_pts
          l = veg_index(i)
          ! NOTE: fsmc is not used by the stomatal optimisation model, hence
          !       no need to check it here.
          IF ( apar(l) == 0.0 ) THEN
            clos_pts = clos_pts + 1
            clos_index(clos_pts) = i
          ELSE
            open_pts = open_pts + 1
            open_index(open_pts) = i
          END IF
        END DO

        !-------------------------------------------------------------------------
        ! Scale to canopy level.
        !-------------------------------------------------------------------------
        DO m = 1,veg_pts
          l = veg_index(m)
          rdc(l) = rd(l) * fpar(l)
          ! Scale into vcmaxc/jmaxc, not vcmax/jmax in place: this block runs
          ! once per pass of the enclosing DO k = 1,iter loop, and vcmax/jmax
          ! are only set once (by calc_photo_parameters) before that loop, so
          ! overwriting them here would compound fpar into them on every
          ! iteration (vcmax*fpar, then *fpar again, ...).
          vcmaxc(l) = vcmax(l) * fpar(l)
          jmaxc(l) = jmax(l) * fpar(l)
          ! Light is scaled by fpar too. The fsmc big leaf computes the top
          ! leaf and multiplies by fpar, and the Farquhar/Collatz rates are
          ! homogeneous of degree one in (light, vcmax, jmax, rd), so this
          ! is the same canopy. Without it, the canopy had fpar times the
          ! capacity but only the top leaf's light, understating the
          ! light-limited rate by up to fpar.
          acrc(l) = acr(l) * fpar(l)

          ! Recalculate je now that jmaxc has been scaled up.
          ! This is done here to keep the canopy scaling in a similar place to
          ! the scaling in the leaf_flux_fsmc model.
          IF ( pft_photo_model == photo_farquhar ) THEN
            ! Calculate the electron flux (i2 is only set for Farquhar).
            i2c(l) = i2(l) * fpar(l)
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2c, jmaxc, je)
          END IF
        END DO

       ! Passing two dummy vars to make it consistent with multilyr, 29 Apr, MGDK
       l_multilayer = .FALSE.

       ! stom_opt_mod's kmax/kcrit dummy arguments are land_pts-sized;
       ! broadcast the whole-plant scalars onto land_pts arrays here rather
       ! than passing the bare scalars kmax_pft(ft)/kcrit(ft) directly (the
       ! latter previously relied on undefined Fortran sequence association
       ! between a scalar actual argument and an array dummy argument, and
       ! could read out of bounds of the kmax_pft/kcrit arrays).
       !
       ! kmax_pft/kcrit are the whole-plant hydraulic conductance per unit
       ! leaf area (mol m-2 s-1 Pa-1), as estimated from site data
       ! (E / (LAI * (psi_predawn - psi_midday)) on well-watered days) and
       ! as used by CABLE. e_leaf passed into stom_opt_mod is canopy scale
       ! (per m2 ground), so the conductance must be too, or leaf_psi_jls
       ! solves psi_leaf = psi_root - e_leaf/k against too small a k. All
       ! leaves draw on the plant's conductance in parallel, so the canopy
       ! value is kmax_pft * LAI. This used to be kmax_pft * fpar, the
       ! leaf-N (vcmax-like) scaling, which made kmax_pft a top-of-canopy
       ! value and understated canopy conductance by LAI / fpar (1.3-1.7
       ! for LAI 1.5-2.3, kpar 0.5) for a data-derived kmax_pft.
       kmax_bigleaf(:) = kmax_pft(ft) * lai(:)
       kcrit_bigleaf(:) = kcrit(ft) * lai(:)
       ! som_gl_max is per unit leaf area, so the canopy cap is x LAI, as
       ! kmax_bigleaf (and the multilayer and two-leaf totals).
       gl_max_bigleaf(:) = som_gl_max * lai(:)
       share_sup(:) = 1.0
       CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,        &
                                share_sup, dqc, tstar, pstar,                  &
                                gl_max_bigleaf, gl_max_eff )
       CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

       CALL stom_opt_mod (                                                  &
              ! IN
                land_pts, som_base_parm, ft, open_pts, open_index,           &
                pft_photo_model, veg_pts, veg_index,                         &
                ca, psi_src, acrc, apar, oa, vcmaxc, kc, ko, ccp, pstar,      &
                km, dqc, qs, je, tstar, je_dummy, fapar_dummy,               &
                kmax_bigleaf, kcrit_bigleaf, gl_max_eff, ipar,               &
                l_multilayer,                                                &
              ! IN OUT
                rdc,                                                           &
              ! OUT
                ci, anetc, el, flux_o3, fo3, gc, psi_leaf,                     &
                carbon_gain, hydraulic_cost, leaf_k                            &
        )

        ! leaf_k from stom_opt_mod is canopy scale (see kmax_bigleaf). Return
        ! it per unit leaf area, as the multilayer path does, so that it is
        ! comparable with kmax_pft (e.g. the PLC_pft diagnostic).
        leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))

      CASE DEFAULT
        errcode = 101  !  a hard error
        CALL ereport(RoutineName, errcode,                                     &
               'leaf_flux_mod should be 1 fsmc regulated or 2 stomatal optimisation')
    END SELECT ! leaf_flux_model

  END DO   ! End of iteration loop

  !---------------------------------------------------------------------------
  ! Calculate canopy-level fluxes.
  !---------------------------------------------------------------------------
  SELECT CASE (leaf_flux_mod)
    CASE (leaf_flux_fsmc)
      DO m = 1,veg_pts
        l = veg_index(m)

        anetc(l) = anetl(l) * fpar(l)
        rdc(l)   = rd(l) * fpar(l)

        IF ( lai(l) > EPSILON(0.0) ) THEN
          rdmean(l) = rd(l) * fpar(l) / lai(l)
        ELSE
          rdmean(l) = rd(l)
        END IF

        IF (l_o3_damage) THEN
          flux_o3(l) = flux_o3_l(l) * fpar(l)
          fo3(l)     = fo3_l(l)
        END IF

      END DO

    ! The stomatal optimisation model is already scaled to canopy.
    ! Only need to calculate rdmean from rdc
    CASE (leaf_flux_stom_opt)
      DO m = 1,veg_pts
        l = veg_index(m)

        IF ( lai(l) > EPSILON(0.0) ) THEN
          rdmean(l) = rdc(l)/ lai(l)
        ELSE
          rdmean(l) = rdc(l)
        END IF
      END DO
  END SELECT

  IF ( stomata_model == stomata_desica .AND. i_pass < n_pass ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      ! Trial E for this fw, with the cuticular floor if it is on (here
      ! without its supply/xylem bounds, which only act near closure).
      gl_cut_ds = 0.0
      IF ( l_som_cuticular_floor ) gl_cut_ds = gcuticular(ft) * 1.0e-3 * rmol  &
                                               * tstar(l) / pstar(l) * lai(l)
      el_try(l) = MAX(dqc(l), 0.0) * pstar(l) / repsilon                     &
                  * MAX(gc(l), gl_cut_ds) / (rmol * tstar(l))
    END DO
    CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,       &
                            lai, ht, psi_root_zone, ksr_ds, el_try, .FALSE., &
                            psi_try, k_try, el_hyd )
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( fsmc_lim(l) > tuzet_fw(ft, psi_try(l)) ) THEN
        fw_hi(l) = fsmc_lim(l)
      ELSE
        fw_lo(l) = fsmc_lim(l)
      END IF
    END DO
  END IF

  END DO  ! i_pass (DESICA bisection)


CASE ( 7 )

  !---------------------------------------------------------------------------
  ! Two-leaf model: one sunlit and one shaded big leaf (de Pury & Farquhar,
  ! 1997; Wang & Leuning, 1998, as in CABLE), both at the surface
  ! temperature tstar, or each at its own leaf temperature with l_leaf_temp.
  !
  ! The sunlit fraction and the absorbed PAR per unit sunlit/shaded leaf area
  ! come from the JULES two-stream profile over ilayers (as for
  ! can_rad_mod = 5/6), and are summed over the layers into two canopy-scale
  ! leaves. Leaf N (and so Vcmax, Jmax, Rd) follows the big-leaf profile
  ! exp(-kpar*L), so the canopy totals of Vcmax, Jmax, Rd and gl_max match
  ! can_rad_mod = 1. The absorbed light is the two-stream canopy total
  ! (<= incident PAR), whereas big-leaf uses acr*fpar, (1-omega)*ipar*
  ! (1-exp(-kpar*L))/kpar, about 1/kpar times more; so expect lower
  ! light-limited GPP than big-leaf, beyond the effect of the sun/shade split. Kmax is per unit
  ! leaf area as in big-leaf, and the canopy kmax_pft*LAI is shared between
  ! the classes by N-weighted leaf area; each class draws from psi_root_zone
  ! on its own parallel path.
  !
  ! l_leaf_temp: each leaf has its own temperature from a leaf energy
  ! balance (Penman-Monteith; Leuning et al., 1995; Wang & Leuning, 1998),
  ! solved with its absorbed PAR, near-IR and isothermal long-wave, its
  ! boundary-layer and radiative conductances and its stomatal conductance
  ! from the last pass. The leaves exchange with a canopy air space that is
  ! coupled to level 1 through ra. Each pass runs the photosynthesis
  ! temperature responses, the humidity deficit and the stomatal
  ! optimisation of each leaf at its own temperature, then updates the leaf
  ! temperatures: leaf_temp_iter passes, so no more stom_opt_mod calls than
  ! without l_leaf_temp (for leaf_temp_iter = iter).
  !---------------------------------------------------------------------------
  IF ( leaf_flux_mod /= leaf_flux_stom_opt ) THEN
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
                 'can_rad_mod = 7 is only coded for leaf_flux_mod = 2')
  END IF

  ! Top-leaf photosynthetic parameters.
  CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,          &
                              veg_index, denom, jmax_temp, jv25,               &
                              nleaf_top, qtenf_term, vcmax_temp,               &
                              jmax, rd, vcmax )

  !---------------------------------------------------------------------------
  ! Aggregate the layer profile into the two leaf classes.
  !---------------------------------------------------------------------------
  DO m = 1,veg_pts
    l = veg_index(m)
    lai_sun_2l(l)  = 0.0
    lai_shd_2l(l)  = 0.0
    nw_sun_2l(l)   = 0.0
    nw_shd_2l(l)   = 0.0
    apar_sun_2l(l) = 0.0
    apar_shd_2l(l) = 0.0
    DO n = 1,ilayers
      IF ( kpar(ft) > EPSILON(0.0) ) THEN
        dnw_2l = ( EXP(-kpar(ft) * REAL(n-1) * dlai(l))                        &
                   - EXP(-kpar(ft) * REAL(n) * dlai(l)) ) / kpar(ft)
      ELSE
        dnw_2l = dlai(l)
      END IF
      lai_sun_2l(l)  = lai_sun_2l(l)  + fsun(l,n) * dlai(l)
      lai_shd_2l(l)  = lai_shd_2l(l)  + (1.0 - fsun(l,n)) * dlai(l)
      nw_sun_2l(l)   = nw_sun_2l(l)   + fsun(l,n) * dnw_2l
      nw_shd_2l(l)   = nw_shd_2l(l)   + (1.0 - fsun(l,n)) * dnw_2l
      apar_sun_2l(l) = apar_sun_2l(l)                                          &
                       + fapar_sun(l,n) * fsun(l,n) * dlai(l) * ipar(l)
      apar_shd_2l(l) = apar_shd_2l(l)                                          &
                       + fapar_shd(l,n) * (1.0 - fsun(l,n)) * dlai(l) * ipar(l)
    END DO

    acr_sun_2l(l) = apar_sun_2l(l) / conpar
    acr_shd_2l(l) = apar_shd_2l(l) / conpar

    vcmax_sun_2l(l) = vcmax(l) * nw_sun_2l(l)
    vcmax_shd_2l(l) = vcmax(l) * nw_shd_2l(l)
    jmax_sun_2l(l)  = jmax(l)  * nw_sun_2l(l)
    jmax_shd_2l(l)  = jmax(l)  * nw_shd_2l(l)

    ! Canopy conductance kmax_pft*LAI (per unit leaf area, as big-leaf),
    ! shared between the classes by their N-weighted leaf area (CABLE's
    ! scalex), so sun leaves get more conductance per unit leaf area than
    ! shade leaves without changing the canopy total (cf. the multilayer
    ! leaf-segment profile). CABLE itself uses kmax_pft * scalex, which
    ! lowers the canopy total by fpar/LAI.
    kmax_sun_2l(l) = kmax_pft(ft) * lai(l) * nw_sun_2l(l)                      &
                     / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    kmax_shd_2l(l) = kmax_pft(ft) * lai(l) * nw_shd_2l(l)                      &
                     / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    kcrit_sun_2l(l) = kmax_sun_2l(l) * (kcrit(ft) / kmax_pft(ft))
    kcrit_shd_2l(l) = kmax_shd_2l(l) * (kcrit(ft) / kmax_pft(ft))
    ! som_gl_max is per unit leaf area: x the class leaf area.
    gl_max_sun_2l(l) = som_gl_max * lai_sun_2l(l)
    gl_max_shd_2l(l) = som_gl_max * lai_shd_2l(l)

    ! Radiation to photosystem II of each class (cf. i2 = alpha_elec*acr).
    i2_sun(l) = alpha_elec(ft) * acr_sun_2l(l)
    i2_shd(l) = alpha_elec(ft) * acr_shd_2l(l)
  END DO

  ! The non-rectangular hyperbola is homogeneous of degree one in (I, Jmax),
  ! so applying it to the class totals is the same as summing identical
  ! leaves.
  IF ( pft_photo_model == photo_farquhar ) THEN
    CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_sun,             &
                             jmax_sun_2l, je_sun )
    CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_shd,             &
                             jmax_shd_2l, je_shd )
  END IF

  l_multilayer = .FALSE.

  !---------------------------------------------------------------------------
  ! Leaf temperature: radiation and conductance terms that do not change
  ! between passes, and the first guess of the leaf temperatures (last
  ! timestep's, or tstar on the first timestep).
  !---------------------------------------------------------------------------
  IF ( l_leaf_temp ) THEN
    wind_ext(:) = 0.0
    zref_a(:)   = 0.0
    IF ( leaf_aero_model == 1 ) THEN
      CALL leaf_aero_geom( land_pts, veg_pts, veg_index, lai, canht, z1_wind,  &
                           z0m_a, zref_a, rt1ab_a, zrd_a, usc_a, usuh_a,       &
                           wind_ext )
    END IF
    CALL leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,      &
                          dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,      &
                          wind_ext,                                            &
                          apar_sun_2l, apar_shd_2l,                            &
                          sw_sun_lt, sw_shd_lt, wlw_sun_lt, wlw_shd_lt,        &
                          wlws_sun_lt, wlws_shd_lt, frad_sun_lt, frad_shd_lt,  &
                          ff_sun_lt, ff_shd_lt )
    n_pass_2l = leaf_temp_iter
    t_g_lt(:)  = tstar(:)
    zeta_lt(:) = 0.0
    ! Store the tstar-based canopy respiration for the stem and root
    ! respiration, which stay at tstar.
    DO m = 1,veg_pts
      l = veg_index(m)
      rdc_tstar(l) = rd(l) * ( nw_sun_2l(l) + nw_shd_2l(l) )
      IF ( t_leaf_sun(l) < 150.0 .OR. t_leaf_shd(l) < 150.0 ) THEN
        t_sun_lt(l) = tstar(l)
        t_shd_lt(l) = tstar(l)
      ELSE
        t_sun_lt(l) = t_leaf_sun(l)
        t_shd_lt(l) = t_leaf_shd(l)
      END IF
      IF ( lai(l) > EPSILON(0.0) ) THEN
        t_g_lt(l) = ( lai_sun_2l(l) * t_sun_lt(l)                              &
                      + lai_shd_2l(l) * t_shd_lt(l) ) / lai(l)
      END IF
      ! Canopy air - level 1 resistance (leaf_aero_model).
      SELECT CASE ( leaf_aero_model )
      CASE ( 1 )
        ! CABLE's rt1, computed in leaf_temp_update (ra where there is no
        ! canopy); here neutral, for the first guess of q_c.
        IF ( zref_a(l) > 0.0 ) THEN
          ra_lt(l) = MAX( 5.0, rt1ab_a(l) * LOG(zref_a(l) / z0m_a(l))          &
                               / ( 0.4 * MAX(u_wind(l), 1.0) ) )
        ELSE
          ra_lt(l) = ra(l)
        END IF
      CASE ( 2 )
        ! JULES's ra.
        ra_lt(l) = ra(l)
      CASE DEFAULT
        ra_lt(l) = 0.0
      END SELECT
      t_c_lt(l) = tair(l)
    END DO
    ! Canopy air humidity from last timestep's conductance and leaf
    ! temperatures (q1 for leaf_aero_model = 0).
    CALL qsat(qs_g_lt, t_g_lt, pstar, land_pts)
    DO m = 1,veg_pts
      l = veg_index(m)
      q_c_lt(l) = ( q1(l) + ra_lt(l) * gc(l) * qs_g_lt(l) )                    &
                  / ( 1.0 + ra_lt(l) * gc(l) )
    END DO
    IF ( leaf_aero_model == 1 ) THEN
      ! The neutral rt1 above is for the first guess only: leaf_temp_update
      ! takes ra where there is no canopy.
      DO m = 1,veg_pts
        l = veg_index(m)
        ra_lt(l) = ra(l)
      END DO
    END IF
    ! Respiration below uses rd and the classes start from it.
    ccp_sun_lt(:) = ccp(:)
    ccp_shd_lt(:) = ccp(:)
    kc_sun_lt(:)  = kc(:)
    kc_shd_lt(:)  = kc(:)
    ko_sun_lt(:)  = ko(:)
    ko_shd_lt(:)  = ko(:)
    km_sun_lt(:)  = km(:)
    km_shd_lt(:)  = km(:)
  ELSE
    n_pass_2l = iter
  END IF

  !---------------------------------------------------------------------------
  ! Iterate for the canopy humidity deficit, starting from the previous
  ! timestep's conductance (as for big-leaf), and with l_leaf_temp for the
  ! leaf temperatures.
  !---------------------------------------------------------------------------
  DO k = 1,n_pass_2l

    DO m = 1,veg_pts
      l = veg_index(m)
      ra_rc(l) = ra(l) * gc(l)
      dqc(l)   = dq(l) / (1.0 + ra_rc(l))
      ! No light inhibition of Rd, as for big-leaf (rd is reset each pass
      ! because stom_opt_mod treats it as IN OUT).
      rd_sun(l) = rd(l) * nw_sun_2l(l)
      rd_shd(l) = rd(l) * nw_shd_2l(l)
    END DO

    IF ( l_leaf_temp ) THEN
      ! Photosynthesis parameters and humidity deficit at each leaf's
      ! temperature. The deficit is from the leaf to the canopy air,
      ! qs(T_leaf) - q_c (with T_leaf = tstar this is dqc).
      CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_sun_lt, oa, actj, actv,      &
                               dsj, dsv, jv25, nleaf_top, nw_sun_2l, i2_sun,   &
                               pstar, q_c_lt, dq_min,                          &
                               ccp_sun_lt, kc_sun_lt, ko_sun_lt, km_sun_lt,    &
                               vcmax_sun_2l, je_sun, rd_sun, qs_sun_lt,        &
                               dq_sun_lt, f_vc = f_vc )
      CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_shd_lt, oa, actj, actv,      &
                               dsj, dsv, jv25, nleaf_top, nw_shd_2l, i2_shd,   &
                               pstar, q_c_lt, dq_min,                          &
                               ccp_shd_lt, kc_shd_lt, ko_shd_lt, km_shd_lt,    &
                               vcmax_shd_2l, je_shd, rd_shd, qs_shd_lt,        &
                               dq_shd_lt, f_vc = f_vc )
    END IF

    !-------------------------------------------------------------------------
    ! Sunlit leaf.
    !-------------------------------------------------------------------------
    clos_pts = 0
    open_pts = 0
    DO i = 1,veg_pts
      l = veg_index(i)
      IF ( apar_sun_2l(l) > 0.0 .AND. lai_sun_2l(l) > EPSILON(0.0) ) THEN
        open_pts = open_pts + 1
        open_index(open_pts) = i
      ELSE
        clos_pts = clos_pts + 1
        clos_index(clos_pts) = i
      END IF
    END DO

    ! Soil supply cap: this class may use its kmax share of the supply.
    DO m = 1,veg_pts
      l = veg_index(m)
      share_sup(l) = kmax_sun_2l(l) /                                          &
                     MAX(kmax_sun_2l(l) + kmax_shd_2l(l), TINY(1.0))
    END DO

    IF ( l_leaf_temp ) THEN
      CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,         &
                               share_sup, dq_sun_lt, t_sun_lt, pstar,          &
                               gl_max_sun_2l, gl_max_eff )
      CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

      CALL stom_opt_mod (                                                      &
              ! IN
                land_pts, som_base_parm, ft, open_pts, open_index,             &
                pft_photo_model, veg_pts, veg_index,                           &
                ca, psi_src, acr_sun_2l, apar_sun_2l, oa, vcmax_sun_2l,        &
                kc_sun_lt, ko_sun_lt, ccp_sun_lt, pstar,                       &
                km_sun_lt, dq_sun_lt, qs_sun_lt, je_sun, t_sun_lt, je_dummy,   &
                fapar_dummy,                                                   &
                kmax_sun_2l, kcrit_sun_2l, gl_max_eff, ipar,                   &
                l_multilayer,                                                  &
              ! IN OUT
                rd_sun,                                                        &
              ! OUT
                ci_sun_2l, anetl_sun, el_sun, flux_o3_l_sun, fo3_l_sun,        &
                gl_sun, psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun               &
        )
    ELSE
    CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,           &
                             share_sup, dqc, tstar, pstar, gl_max_sun_2l,      &
                             gl_max_eff )
    CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

    CALL stom_opt_mod (                                                        &
            ! IN
              land_pts, som_base_parm, ft, open_pts, open_index,               &
              pft_photo_model, veg_pts, veg_index,                             &
              ca, psi_src, acr_sun_2l, apar_sun_2l, oa, vcmax_sun_2l,          &
              kc, ko, ccp, pstar,                                              &
              km, dqc, qs, je_sun, tstar, je_dummy, fapar_dummy,               &
              kmax_sun_2l, kcrit_sun_2l, gl_max_eff, ipar,                     &
              l_multilayer,                                                    &
            ! IN OUT
              rd_sun,                                                          &
            ! OUT
              ci_sun_2l, anetl_sun, el_sun, flux_o3_l_sun, fo3_l_sun, gl_sun,  &
              psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun                         &
      )
    END IF

    ! Closed leaves get min_gl_pft from stom_opt_mod as a canopy value;
    ! share it by leaf area so the two classes together give big-leaf's
    ! (half each when there is no leaf area).
    DO i = 1,clos_pts
      l = veg_index(clos_index(i))
      IF ( lai(l) > EPSILON(0.0) ) THEN
        gl_sun(l) = gl_sun(l) * lai_sun_2l(l) / lai(l)
      ELSE
        gl_sun(l) = 0.5 * gl_sun(l)
      END IF
    END DO

    !-------------------------------------------------------------------------
    ! Shaded leaf.
    !-------------------------------------------------------------------------
    clos_pts = 0
    open_pts = 0
    DO i = 1,veg_pts
      l = veg_index(i)
      IF ( apar_shd_2l(l) > 0.0 .AND. lai_shd_2l(l) > EPSILON(0.0) ) THEN
        open_pts = open_pts + 1
        open_index(open_pts) = i
      ELSE
        clos_pts = clos_pts + 1
        clos_index(clos_pts) = i
      END IF
    END DO

    ! Soil supply cap: this class may use its kmax share of the supply.
    DO m = 1,veg_pts
      l = veg_index(m)
      share_sup(l) = kmax_shd_2l(l) /                                          &
                     MAX(kmax_sun_2l(l) + kmax_shd_2l(l), TINY(1.0))
    END DO

    IF ( l_leaf_temp ) THEN
      CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,         &
                               share_sup, dq_shd_lt, t_shd_lt, pstar,          &
                               gl_max_shd_2l, gl_max_eff )
      CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

      CALL stom_opt_mod (                                                      &
              ! IN
                land_pts, som_base_parm, ft, open_pts, open_index,             &
                pft_photo_model, veg_pts, veg_index,                           &
                ca, psi_src, acr_shd_2l, apar_shd_2l, oa, vcmax_shd_2l,        &
                kc_shd_lt, ko_shd_lt, ccp_shd_lt, pstar,                       &
                km_shd_lt, dq_shd_lt, qs_shd_lt, je_shd, t_shd_lt, je_dummy,   &
                fapar_dummy,                                                   &
                kmax_shd_2l, kcrit_shd_2l, gl_max_eff, ipar,                   &
                l_multilayer,                                                  &
              ! IN OUT
                rd_shd,                                                        &
              ! OUT
                ci_shd_2l, anetl_shd, el_shd, flux_o3_l_shd, fo3_l_shd,        &
                gl_shd, psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd               &
        )
    ELSE
    CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,           &
                             share_sup, dqc, tstar, pstar, gl_max_shd_2l,      &
                             gl_max_eff )
    CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

    CALL stom_opt_mod (                                                        &
            ! IN
              land_pts, som_base_parm, ft, open_pts, open_index,               &
              pft_photo_model, veg_pts, veg_index,                             &
              ca, psi_src, acr_shd_2l, apar_shd_2l, oa, vcmax_shd_2l,          &
              kc, ko, ccp, pstar,                                              &
              km, dqc, qs, je_shd, tstar, je_dummy, fapar_dummy,               &
              kmax_shd_2l, kcrit_shd_2l, gl_max_eff, ipar,                     &
              l_multilayer,                                                    &
            ! IN OUT
              rd_shd,                                                          &
            ! OUT
              ci_shd_2l, anetl_shd, el_shd, flux_o3_l_shd, fo3_l_shd, gl_shd,  &
              psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd                         &
      )
    END IF

    DO i = 1,clos_pts
      l = veg_index(clos_index(i))
      IF ( lai(l) > EPSILON(0.0) ) THEN
        gl_shd(l) = gl_shd(l) * lai_shd_2l(l) / lai(l)
      ELSE
        gl_shd(l) = 0.5 * gl_shd(l)
      END IF
    END DO

    DO m = 1,veg_pts
      l = veg_index(m)
      gc(l) = gl_sun(l) + gl_shd(l)
    END DO

    ! New leaf and canopy air temperatures for the next pass (none after
    ! the last pass, so the temperatures output are those used).
    IF ( l_leaf_temp .AND. k < n_pass_2l ) THEN
      CALL leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             ra_lt, lw_down, u_wind, lai_sun_2l, lai_shd_2l,      &
                             sw_sun_lt, sw_shd_lt, wlw_sun_lt, wlw_shd_lt,     &
                             wlws_sun_lt, wlws_shd_lt, frad_sun_lt,            &
                             frad_shd_lt, ff_sun_lt, ff_shd_lt, z0m_a, zref_a, &
                             rt1ab_a, zrd_a, usc_a, usuh_a, gl_sun,            &
                             gl_shd, zeta_lt, t_sun_lt, t_shd_lt, t_c_lt,      &
                             q_c_lt )
    END IF

  END DO   ! End of iteration loop

  !---------------------------------------------------------------------------
  ! Canopy totals and leaf-area-weighted means.
  !---------------------------------------------------------------------------
  DO m = 1,veg_pts
    l = veg_index(m)

    IF ( lai(l) > EPSILON(0.0) ) THEN
      f_sun_2l = lai_sun_2l(l) / lai(l)
    ELSE
      f_sun_2l = 0.5
    END IF

    anetc(l) = anetl_sun(l) + anetl_shd(l)
    rdc(l)   = rd_sun(l) + rd_shd(l)
    el(l)    = el_sun(l) + el_shd(l)
    ! Sun + shade is the canopy conductance; per unit leaf area, as the
    ! big-leaf path.
    leaf_k(l) = ( leaf_k_sun(l) + leaf_k_shd(l) )                              &
                / MAX(lai(l), TINY(1.0_real_jlslsm))

    psi_leaf(l)       = f_sun_2l * psi_leaf_sun(l)                             &
                        + (1.0 - f_sun_2l) * psi_leaf_shd(l)
    carbon_gain(l)    = f_sun_2l * CG_sun(l) + (1.0 - f_sun_2l) * CG_shd(l)
    hydraulic_cost(l) = f_sun_2l * HC_sun(l) + (1.0 - f_sun_2l) * HC_shd(l)

    IF ( gc(l) > TINY(gc(l)) ) THEN
      ci(l) = ( gl_sun(l) * ci_sun_2l(l) + gl_shd(l) * ci_shd_2l(l) ) / gc(l)
    ELSE
      ci(l) = f_sun_2l * ci_sun_2l(l) + (1.0 - f_sun_2l) * ci_shd_2l(l)
    END IF

    IF ( lai(l) > EPSILON(0.0) ) THEN
      rdmean(l) = rdc(l) / lai(l)
    ELSE
      rdmean(l) = rdc(l)
    END IF

    IF (l_o3_damage) THEN
      flux_o3(l) = flux_o3_l_sun(l) + flux_o3_l_shd(l)
      fo3(l)     = f_sun_2l * fo3_l_sun(l) + (1.0 - f_sun_2l) * fo3_l_shd(l)
    END IF

    IF ( l_leaf_temp ) THEN
      ! Leaf temperatures used. A class with no leaf area takes the other's.
      IF ( lai_sun_2l(l) > EPSILON(0.0) ) THEN
        t_leaf_sun(l) = t_sun_lt(l)
      ELSE
        t_leaf_sun(l) = t_shd_lt(l)
      END IF
      IF ( lai_shd_2l(l) > EPSILON(0.0) ) THEN
        t_leaf_shd(l) = t_shd_lt(l)
      ELSE
        t_leaf_shd(l) = t_sun_lt(l)
      END IF
      t_can(l) = f_sun_2l * t_leaf_sun(l) + (1.0 - f_sun_2l) * t_leaf_shd(l)

      ! The surface energy balance applies gc with qs(tstar) - q1 through
      ! ra. Return the gc that gives the transpiration the optimisation
      ! chose (el, at the leaf temperatures) there:
      !   rho (qs(tstar) - q1) / (ra + 1/gc) = el * m_h2o,
      ! capped at the canopy som_gl_max.
      IF ( l_leaf_temp_gc_eq ) THEN
        gc_cap = gl_max_sun_2l(l) + gl_max_shd_2l(l)
        IF ( el(l) * m_h2o_lt > TINY(1.0) ) THEN
          r_eq = pstar(l) / ( r * tstar(l) ) * dq(l) / ( el(l) * m_h2o_lt )    &
                 - ra(l)
          IF ( r_eq > 0.0 ) THEN
            gc(l) = 1.0 / r_eq
            IF ( gc_cap > 0.0 ) gc(l) = MIN(gc(l), gc_cap)
          ELSE IF ( gc_cap > 0.0 .AND. dq(l) > 0.0 ) THEN
            ! Demand above what any gc gives: the cap. (With dq <= 0 the
            ! tile balance has condensation at any gc; keep the optimised
            ! gc rather than condensing through stomata at the cap.)
            gc(l) = gc_cap
          END IF
        END IF
      END IF
    END IF
  END DO

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'can_rad_mod should be 1, 4, 5, 6 or 7')

END SELECT  ! can_rad_mod

!-----------------------------------------------------------------------------
! DESICA: as the stomatal optimisation, no stomatal loss with no light
! (gc = min_gl_pft, as stom_opt_mod's closed points), and the canopy
! transpiration (as in stom_opt_mod_ci) for the cuticular floor below and
! the plant hydraulics after it. psi_leaf/leaf_k here are only the
! steady-state values the floor block expects; desica_hydraulics replaces
! them.
!-----------------------------------------------------------------------------
IF ( stomata_model == stomata_desica ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( apar(l) == 0.0 ) gc(l) = min_gl_pft(ft)
    el(l)       = MAX(dqc(l), 0.0) * pstar(l) / repsilon * gc(l)               &
                  / (rmol * tstar(l))
    psi_leaf(l) = psi_root_zone(l)
    leaf_k(l)   = kmax_pft(ft)   ! per unit leaf area, as the other paths
  END DO
END IF

!-----------------------------------------------------------------------------
! Cuticular floor (l_som_cuticular_floor). Water still leaks through the
! cuticle once the stomata have (nearly) shut, whatever the carbon: the
! canopy conductance is not allowed below gcuticular * LAI. This is applied after
! the optimisation, so it is not traded against carbon (A is unchanged), but
! it is charged to the plant water: the soil-supply cap still bounds it and
! psi_leaf / leaf_k are re-solved for the total flux. Applies at night too
! (closed points), as cuticular loss is not stomatal.
!-----------------------------------------------------------------------------
IF ( l_som_cuticular_floor .AND. ( leaf_flux_mod == leaf_flux_stom_opt .OR.  &
                                   stomata_model == stomata_desica ) ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    ! mmol m-2 leaf s-1 -> m s-1, times LAI for the canopy.
    gl_cut(l) = gcuticular(ft) * 1.0e-3 * rmol * tstar(l) / pstar(l) * lai(l)
    share_sup(l) = 1.0
    veg_pts_index(m) = m
  END DO
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,             &
                           share_sup, dqc, tstar, pstar, gl_cut, gl_cut_eff )
  CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
  ! Without plant water storage the leak must also pass the xylem: try
  ! E_floor * i / n_cut and keep the largest flux with k > kcrit (the
  ! feasibility rule of stom_opt_mod), so psi_leaf cannot run away in a
  ! drought once the soil is near root_psi_crit.
  DO m = 1,veg_pts
    l = veg_index(m)
    kmax_cut(l)  = MAX(kmax_pft(ft) * lai(l), TINY(1.0_real_jlslsm))
    kcrit_cut(l) = kcrit(ft) * lai(l)
    l_cut(l) = gl_cut(l) > 0.0 .AND. gc(l) < gl_cut_eff(l)
    ! Canopy transpiration (mol H2O m-2 s-1) at the floor, as in
    ! stom_opt_mod_ci.
    el_cut(l) = 0.0
    IF ( l_cut(l) ) el_cut(l) = MAX(dqc(l), 0.0) * pstar(l) / repsilon         &
                                * gl_cut_eff(l) / (rmol * tstar(l))
    DO i_cut = 1,n_cut
      e_cut(i_cut,m) = el_cut(l) * REAL(i_cut) / REAL(n_cut)
    END DO
  END DO
  IF ( ANY(l_cut(veg_index(1:veg_pts))) ) THEN
    CALL leaf_psi_jls( ft, n_cut, land_pts, veg_pts, veg_index,                &
                       veg_pts_index, e_cut, psi_src,                          &
                       kmax_cut, kcrit_cut,                                    &
                     ! OUT
                       psi_cut, k_cut )
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( .NOT. l_cut(l) ) CYCLE
      DO i_cut = n_cut,1,-1
        IF ( k_cut(i_cut,m) > kcrit_cut(l) .AND.                               &
             psi_cut(i_cut,m) <= psi_src(l) ) EXIT
      END DO
      ! Only raise the flux: the stomatal optimum itself was feasible.
      IF ( i_cut >= 1 ) THEN
        IF ( gl_cut_eff(l) * REAL(i_cut) / REAL(n_cut) > gc(l) ) THEN
          gc(l)       = gl_cut_eff(l) * REAL(i_cut) / REAL(n_cut)
          el(l)       = e_cut(i_cut,m)
          psi_leaf(l) = psi_cut(i_cut,m)
          ! k_cut is canopy scale (kmax_cut); leaf_k is per unit leaf area.
          leaf_k(l)   = k_cut(i_cut,m) / MAX(lai(l), TINY(1.0_real_jlslsm))
        END IF
      END IF
    END DO
  END IF
END IF

!-----------------------------------------------------------------------------
! DESICA: project psi_leaf and psi_stem for this timestep's transpiration
! (including any cuticular floor). Where the plant cannot deliver it
! without psi_stem or psi_leaf passing their lower bounds, gc is cut to the
! transpiration it can deliver (A is not re-solved: by then the Tuzet
! factor has shut the stomata). The state is advanced in sf_evap, with the
! actual transpiration; here the step's inputs are stored for that.
!-----------------------------------------------------------------------------
IF ( stomata_model == stomata_desica ) THEN
  CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,         &
                          lai, ht, psi_root_zone, ksr_ds, el, .FALSE.,       &
                          psi_try, k_try, el_hyd )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( el_hyd(l) < el(l) ) THEN
      gc(l) = gc(l) * el_hyd(l) / el(l)
      el(l) = el_hyd(l)
    END IF
  END DO
  CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,         &
                          lai, ht, psi_root_zone, ksr_ds, el, .FALSE.,       &
                          psi_leaf, leaf_k, el_hyd )
  ! desica_hydraulics returns Kplant (soil-to-leaf) per m2 ground; per unit
  ! leaf area, as the other paths (kplant_pft).
  leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))
  CALL desica_store_inputs( ft, land_pts, veg_pts, veg_index, lai, ht,        &
                            psi_root_zone, ksr_ds )
END IF

!-----------------------------------------------------------------------------
! Calculate ratio of intercellular CO2 to canopy level atmospheric CO2
!  concentration. Only used as output.
!  Note: TINY is used to avoid division by zero.
!-----------------------------------------------------------------------------
!TODO: Add OMP
DO m = 1,veg_pts
  l = veg_index(m)
  cica_ratio(l) = ci(l) / (ca(l) + TINY(ca(l)))
END DO

!-----------------------------------------------------------------------------
! Calculate plant-level respiration, NPP and GPP.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Calculate the conversion from top-leaf to canopy-average nitrogen.
!-----------------------------------------------------------------------------
can_averaging_fac(:) =                                                         &
    get_can_ave_fac( ft, land_pts, veg_pts, veg_index, lai )

!$OMP PARALLEL DO                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(SHARED)                                                          &
!$OMP PRIVATE(m,l,stemc,lma_tmp,fstem,stem_resp_scaling)
DO m = 1,veg_pts
  l = veg_index(m)

  !---------------------------------------------------------------------------
  ! Calculate the actual and balanced mean leaf nitrogen concentration
  ! assuming perfect light acclimation, then
  ! calculate the total nitrogen content of the leaf, root and stem
  ! Assume that root biomass is equal to balanced growth leaf biomass
  !---------------------------------------------------------------------------
  lai_bal(l) = (a_ws(ft) * eta_sl(ft) * ht(l) / a_wl(ft))                      &
               **(1.0 / (b_wl(ft) - 1.0))

  IF (l_red) lai_bal(l) = veg_state%lai_bal(l,ft)

  !---------------------------------------------------------------------------
  ! Calculate the total nitrogen content of the leaf, root and stem
  !---------------------------------------------------------------------------
  IF ( ft > nnpft ) THEN
    !Crop PFTs
    stemc   = stemc_from_prognostics(ft - nnpft, canht(l) )
    root(l) = rootc_cpft(l,ft - nnpft)

    IF ( l_trait_phys ) THEN

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft, can_averaging_fac(l) )
      n_root(l) = nr_nl(ft) * nmass(ft) * (1.0 / cfrac_l(ft - nnpft))          &
                  * root(l)
      n_stem(l) = ns_nl(ft) * nmass(ft) * (1.0 / cfrac_l(ft - nnpft)) * stemc
        ! nmass(ft)/cfrac_l(ft-nnpft) is equivalent to nl0(ft)

      lma_tmp   = lma_from_prognostics(ft - nnpft, dvi_cpft(l,ft - nnpft))
      nlmean(l) = nmass(ft) * lma_tmp * can_averaging_fac(l)

    ELSE

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft, can_averaging_fac(l) )
      n_root(l) = nr_nl(ft) * nl0(ft) * root(l)
      n_stem(l) = ns_nl(ft) * nl0(ft) * stemc

      nlmean(l) = nl0(ft) * can_averaging_fac(l) * cfrac_l(ft - nnpft)         &
                  * lma_from_prognostics(ft - nnpft, dvi_cpft(l,ft - nnpft))
    END IF

  ELSE

    !Non-crop PFTs
    IF ( l_trait_phys ) THEN
      root(l) = lma(ft) * lai_bal(l)
      !Note new units of temporary variable root:
      !kg root/m2= gleaf/m2 * kg/g

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft, can_averaging_fac(l))
      nlmean(l) = nmass(ft) * lma(ft) * can_averaging_fac(l)

      n_root(l) = nr(ft) * root(l) * cmass

      !Initial calculation of N content in respiring stem wood
      n_stem(l) = eta_sl(ft) * ht(l) * lai_bal(l) * nsw(ft)
      !USE veg3 allometry
      IF (l_red) n_stem(l) = veg_state%woodC(l,ft) / a_ws(ft) * nsw(ft)

      !Reduce n_stem for consistency with non-trait n_stem for now
      !This must be done to achieve realistic respiration rates.
      fstem            = 1.0 / a_ws(ft)
      stem_resp_scaling = fstem + (1.0 - fstem) * hw_sw(ft)
      n_stem(l)         = n_stem(l) * stem_resp_scaling

    ELSE

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft, can_averaging_fac(l))
      root(l)   = sigl(ft) * lai_bal(l)
      n_root(l) = nr_nl(ft) * nl0(ft) * root(l)

      IF ( l_stem_resp_fix ) THEN
        n_stem(l) = ns_nl(ft) * nl0(ft) * eta_sl(ft) * ht(l) * lai_bal(l)
      ELSE
        n_stem(l) = ns_nl(ft) * nl0(ft) * eta_sl(ft) * ht(l) * lai(l)
      END IF

      nlmean(l) = nl0(ft) * can_averaging_fac(l) * sigl(ft)

    END IF

  END IF

  !---------------------------------------------------------------------------
  ! Calculate the Gross Primary Productivity, the plant maintenance
  ! respiration rate, and the wood maintenance respiration rate
  ! in kg C/m2/sec
  !---------------------------------------------------------------------------
  ! GPP = net + leaf dark respiration, using the same rd scaling that went
  ! into anetc. leaf() (leaf_flux_fsmc) computes al = (wl - rd) * fsmc, so
  ! rd * fsmc is added back; stom_opt_mod (leaf_flux_stom_opt) computes
  ! al = wl - rd with no fsmc (water stress acts through the hydraulics
  ! instead), so rd must be added back unscaled. Using fsmc for both left
  ! stom_opt's GPP short by rdc * (1 - fsmc) - and negative at night,
  ! where closed-stomata anetc = -rdc.
  IF (leaf_flux_mod == leaf_flux_stom_opt) THEN
    fsmc_leaf_resp(l) = 1.0
  ELSE
    fsmc_leaf_resp(l) = fsmc(l)
  END IF

  gpp(l) = cconu * (anetc(l) + rdc(l) * fsmc_leaf_resp(l))

  IF (l_scale_resp_pm) THEN
    fsmc_scale(l) = fsmc(l)
  END IF

  IF ( l_sugar ) THEN
    !-------------------------------------------------------------------------
    ! Calculate carbon contents of leaf root and wood
    !-------------------------------------------------------------------------
    IF ( l_trait_phys ) THEN

      leafc(l)     = lma(ft) * lai(l) * cmass
      leafc_bal(l) = lma(ft) * lai_bal(l) * cmass

      IF ( l_red ) THEN
        rootc(l) = veg_state%rootC(l,ft)
        woodc(l) = veg_state%woodC(l,ft)
      ELSE
        rootc(l) = lma(ft) * lai_bal(l) * cmass
        woodc(l) = a_ws(ft) * eta_sl(ft) * ht(l) * lai_bal(l)
      END IF ! ( l_red )

    ELSE

      leafc(l)     = sigl(ft) * lai(l)
      leafc_bal(l) = sigl(ft) * lai_bal(l)
      rootc(l) = lma(ft) * lai_bal(l) * cmass
      woodc(l) = a_ws(ft) * eta_sl(ft) * ht(l) * lai_bal(l)

    END IF ! ( l_trait_phys )
    !---------------------------------------------------------------------------
    ! Calculate respration using SUGAR model
    !---------------------------------------------------------------------------
    CALL sugar(ft, leafc(l), leafc_bal(l) , woodc(l), rootc(l), f_nsc(l),      &
               resp_l(l), resp_w(l), resp_r(l), resp_p_m(l), resp_p_g(l),      &
               resp_p(l), growth_sug(l), tstar(l), gpp(l))

  ELSE
    IF ( lai(l) > EPSILON(0.0) .AND. l_leaf_temp ) THEN
      ! Leaf respiration at the leaf temperatures; stem and root respiration
      ! scale with the canopy Rd at tstar, as without l_leaf_temp.
      resp_w(l) = cconu * rdc_tstar(l) * n_stem(l) * fsmc_scale(l) / n_leaf(l)
      resp_r(l) = cconu * rdc_tstar(l) * n_root(l) * fsmc_scale(l) / n_leaf(l)
      resp_l(l) = cconu * rdc(l) * fsmc_leaf_resp(l)
      resp_p_m(l) = resp_l(l) + resp_w(l) + resp_r(l)
    ELSE IF ( lai(l) > EPSILON(0.0) ) THEN
      resp_p_m(l) = cconu * rdc(l)                                             &
           * (n_leaf(l) * fsmc_leaf_resp(l) + n_stem(l) * fsmc_scale(l) +      &
              n_root(l) * fsmc_scale(l)) / n_leaf(l)
      resp_w(l) = cconu * rdc(l) * n_stem(l) * fsmc_scale(l) / n_leaf(l)
      resp_r(l) = cconu * rdc(l) * n_root(l) * fsmc_scale(l) / n_leaf(l)
      resp_l(l) = cconu * rdc(l) * fsmc_leaf_resp(l)
    ELSE
      resp_w(l) = cconu * rdmean(l) * n_stem(l) * fsmc_scale(l) / nlmean(l)
      resp_r(l) = cconu * rdmean(l) * n_root(l) * fsmc_scale(l) / nlmean(l)
      resp_l(l) = cconu * rdc(l) * fsmc_leaf_resp(l)
      resp_p_m(l) = resp_w(l) + resp_r(l) + resp_l(l)
    END IF
    !-------------------------------------------------------------------------
    ! Calculate the total plant respiration
    !-------------------------------------------------------------------------
    resp_p_g(l) = r_grow(ft) * (gpp(l) - resp_p_m(l))
    resp_p(l)   = resp_p_m(l) + resp_p_g(l)
  END IF !l_sugar

  !---------------------------------------------------------------------------
  ! Calculate Net Primary Productivity
  !---------------------------------------------------------------------------
  npp(l)      = gpp(l) - resp_p(l)

END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------------
! Calculate BVOC emissions
!-----------------------------------------------------------------------------
IF ( l_bvoc_emis ) THEN
  CALL bvoc_emissions(land_pts,ft,veg_index,                                   &
                      open_pts,open_index,clos_pts,clos_index,                 &
                      lai,ci,gpp,tstar,                                        &
                      isoprene,terpene,methanol,acetone,                       &
                      !New arguments replacing USE statements
                      !crop_vars_mod (IN)
                      dvi_cpft)
END IF

!-----------------------------------------------------------------------------
! If this functional type is a crop that has not emerged in a particular grid
! box, set all the output variables to zero.
!-----------------------------------------------------------------------------
IF ( ft > nnpft ) THEN

!$OMP PARALLEL DO                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(veg_pts,veg_index,dvi_cpft,nnpft,gpp,resp_w,resp_p,npp,gc,fo3,    &
!$OMP         flux_o3,isoprene,terpene,methanol,acetone,ft)                    &
!$OMP PRIVATE(m,l)
  DO m = 1,veg_pts
    l = veg_index(m)

    IF ( dvi_cpft(l,ft - nnpft) < 0.0 ) THEN
      gpp(l)      = 0.0
      resp_w(l)   = 0.0
      resp_p(l)   = 0.0
      npp(l)      = 0.0
      gc(l)       = 0.0
      fo3(l)      = 0.0
      flux_o3(l)  = 0.0
      isoprene(l) = 0.0
      terpene(l)  = 0.0
      methanol(l) = 0.0
      acetone(l)  = 0.0
    END IF
  END DO
!$OMP END PARALLEL DO
END IF

!-----------------------------------------------------------------------------
! Calculate FAPAR diagnostic (fraction of absorbed photosynthetically
! active radiation). N.b. this is not per LAI.
! Calculate APAR diagnostic (absorbed photosynthetically
! active radiation).
!-----------------------------------------------------------------------------

fapar_diag(:) = 0.0
apar_diag(:) = 0.0

IF ( l_fapar_diag ) THEN
  SELECT CASE ( can_rad_mod )
  CASE ( 1 )
    DO m = 1,veg_pts
      l = veg_index(m)
      fapar_diag(l) = (1.0 - omega(ft)) * fpar(l) * kpar(ft)
    END DO
  CASE ( 4 )
    DO m = 1,veg_pts
      l = veg_index(m)
      fapar_diag(l) = SUM(faparv(l,:)) * dlai(l)
    END DO
  CASE ( 5, 6, 7 )
    DO m = 1,veg_pts
      l = veg_index(m)
      fapar_diag(l) = SUM( fsun(l,:) * fapar_sun(l,:) +                        &
                           ( 1.0 - fsun(l,:) ) * fapar_shd(l,:)                &
                         ) * dlai(l)
    END DO
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
                 'can_rad_mod should be 1, 4, 5, 6 or 7 (l_fapar_diag)')
  END SELECT
  apar_diag(:) = fapar_diag(:) * ipar(:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE sf_stom

!#############################################################################
!#############################################################################

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

USE jules_vegetation_mod, ONLY: l_leaf_coexp_lai

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), canht(land_pts),          &
                                      z1(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: z0m(land_pts), zref(land_pts),          &
                                       rt1ab(land_pts), zrd(land_pts),         &
                                       usc(land_pts), usuh(land_pts),          &
                                       coexp(land_pts)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  vonk = 0.40, a33 = 1.25, csw = 0.50, ctl = 0.40, crd = 0.3, csd = 0.003,     &
  ccd = 15.0, ccw_c = 2.0, usuhm = 0.3
      ! CABLE physical constants (cable_data).

REAL(KIND=real_jlslsm) :: xx, dh, disp, h, lai_h, term2, term3, term5, zruffs
INTEGER :: l, m

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
  z0m(l)   = ( (1.0 - dh) * EXP( LOG( ccw_c ) - 1.0 + 1.0 / ccw_c             &
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

END SUBROUTINE leaf_aero_geom

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,    &
                            dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,    &
                            wind_ext,                                          &
                            apar_sun, apar_shd,                                &
                            sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,        &
                            wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd )

! Two-leaf leaf energy balance (l_leaf_temp): the terms that are fixed over
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

USE pftparm, ONLY: orient

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

REAL(KIND=real_jlslsm) :: kd, xi0, xi1, dtop, dbot, dwind
INTEGER :: l, m, n

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
  END DO
END DO

END SUBROUTINE leaf_temp_setup

!#############################################################################
!#############################################################################

SUBROUTINE leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_leaf, oa, actj, actv, dsj,   &
                               dsv, jv25, nleaf_top, nw, i2c, pstar, q_c,      &
                               dq_min, ccp, kc, ko, km, vcmax_c, je_c, rd_c,   &
                               qs_c, dq_c, f_vc )

! Two-leaf leaf energy balance (l_leaf_temp): the photosynthesis parameters
! of one leaf class at its temperature t_leaf - the temperature responses,
! the class totals of Vcmax and Rd (N-weighted leaf area nw), the electron
! flux - and its humidity deficit to the canopy air, qs(t_leaf) - q_c.

USE jules_vegetation_mod, ONLY: photo_farquhar
USE qsat_mod, ONLY: qsat

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts),             &
                       pft_photo_model

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  t_leaf(land_pts), oa(land_pts), actj(land_pts), actv(land_pts),              &
  dsj(land_pts), dsv(land_pts), jv25(land_pts), nleaf_top(land_pts),           &
  nw(land_pts), i2c(land_pts), pstar(land_pts), q_c(land_pts), dq_min

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ccp(land_pts), kc(land_pts), ko(land_pts), km(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  vcmax_c(land_pts), je_c(land_pts), rd_c(land_pts), qs_c(land_pts),           &
  dq_c(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: f_vc(land_pts)
    ! Soil-water capacity factor (l_som_vcmax_psi).

REAL(KIND=real_jlslsm) ::                                                      &
  acr0(land_pts), denom_t(land_pts), qtenf_t(land_pts), vtemp(land_pts),       &
  jtemp(land_pts), i2_t(land_pts), jmax_t(land_pts), rd_t(land_pts),           &
  vcmax_t(land_pts), jmax_c(land_pts)
INTEGER :: l, m

acr0(:)    = 0.0
denom_t(:) = 1.0
qtenf_t(:) = 1.0
vtemp(:)   = 0.0
jtemp(:)   = 0.0
i2_t(:)    = 0.0
jmax_t(:)  = 0.0
rd_t(:)    = 0.0
vcmax_t(:) = 0.0
vcmax_c(:) = 0.0
jmax_c(:)  = 0.0
je_c(:)    = 0.0
rd_c(:)    = 0.0
dq_c(:)    = 0.0

CALL leaf_temp_responses( ft, land_pts, veg_pts, veg_index,                    &
                          pft_photo_model, t_leaf, oa, acr0, actj, actv,       &
                          dsj, dsv, ccp, kc, ko, km, denom_t, qtenf_t,         &
                          vtemp, jtemp, i2_t, f_vc = f_vc )

CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,            &
                            veg_index, denom_t, jtemp, jv25,                   &
                            nleaf_top, qtenf_t, vtemp,                         &
                            jmax_t, rd_t, vcmax_t )

DO m = 1,veg_pts
  l = veg_index(m)
  vcmax_c(l) = vcmax_t(l) * nw(l)
  jmax_c(l)  = jmax_t(l)  * nw(l)
  rd_c(l)    = rd_t(l)    * nw(l)
END DO

IF ( pft_photo_model == photo_farquhar ) THEN
  CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2c, jmax_c, je_c )
END IF

CALL qsat(qs_c, t_leaf, pstar, land_pts)
DO m = 1,veg_pts
  l = veg_index(m)
  dq_c(l) = MAX(dq_min, qs_c(l) - q_c(l))
END DO

END SUBROUTINE leaf_class_at_temp

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             ra, lw_down, u_wind, lai_sun, lai_shd,            &
                             sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,       &
                             wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,     &
                             z0m_a, zref_a, rt1ab_a, zrd_a, usc_a, usuh_a,     &
                             gl_sun, gl_shd, zeta_a, t_sun, t_shd, t_c, q_c )

! Two-leaf leaf energy balance (l_leaf_temp): sunlit and shaded leaf
! temperatures, and canopy air temperature and humidity, for the stomatal
! conductances of the last pass.
!
! Each leaf class exchanges with the canopy air (t_c, q_c) through its
! boundary layer, as a big leaf: Penman-Monteith with the isothermal net
! radiation, linearised about t_c (Leuning et al., 1995, App. E; Wang &
! Leuning, 1998), as in the two-leaf model of De Kauwe (two_leaf_at_WTC):
!   LE = (s Rn + D gh cp) / (s + gamma gh / gw),
!   T  = t_c + (Rn - LE) / (cp gh),   gh = 2 gbH + grn,
! with forced convection for gbH (per unit leaf area at the canopy-top wind
! u_wind, divided by leaf_shelter, times the class leaf area weighted by
! the wind profile in the canopy, ff) and free convection (times the class
! leaf area), and the radiative conductance grn.
!
! The canopy air is coupled to level 1 through ra (from leaf_aero_model:
! 0, so t_c = tair and q_c = q1; 1, CABLE's rt1; 2, JULES's ra):
!   q_c = q1 + ra (LE_sun + LE_shd) / (rho lc),
!   t_c = tair + ra (H_sun + H_shd) / (rho cp).
! The balance is solved for these conductances: LE is linear in q_c, so q_c
! is explicit for a given t_c, and t_c is found by regula falsi within
! dt_can of tair (at once when ra = 0). Leaf temperatures are kept
! within dt_leaf of t_c.

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE csigma, ONLY: sbcon
USE planet_constants_mod, ONLY: cp, r, repsilon
USE water_constants_mod, ONLY: lc
USE jules_vegetation_mod, ONLY: leaf_width, leaf_shelter, leaf_aero_model
USE qsat_mod, ONLY: qsat

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tair(land_pts), q1(land_pts), pstar(land_pts), ra(land_pts),                 &
  lw_down(land_pts), u_wind(land_pts), lai_sun(land_pts), lai_shd(land_pts),   &
  sw_sun(land_pts), sw_shd(land_pts), wlw_sun(land_pts), wlw_shd(land_pts),    &
  wlws_sun(land_pts), wlws_shd(land_pts), frad_sun(land_pts),                  &
  frad_shd(land_pts), ff_sun(land_pts), ff_shd(land_pts), gl_sun(land_pts),    &
  gl_shd(land_pts),                                                            &
  z0m_a(land_pts), zref_a(land_pts), rt1ab_a(land_pts), zrd_a(land_pts),       &
  usc_a(land_pts), usuh_a(land_pts)
      ! CABLE canopy geometry (leaf_aero_model = 1; see leaf_aero_geom).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  zeta_a(land_pts),                                                            &
      ! zref/L (leaf_aero_model = 1): the start of the Monin-Obukhov
      ! iteration, and the last iterate on output, so successive passes in
      ! a timestep continue the iteration (0 for the first).
  t_sun(land_pts), t_shd(land_pts), t_c(land_pts), q_c(land_pts)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  vonk = 0.40, grav = 9.81, zetneg = -15.0, zetpos = 1.0, rt_min = 5.0,        &
      ! CABLE constants: von Karman, gravity, limits on zref/L, minimum rt1.
  umin_cable = 1.0, apol = 0.70, prandt = 0.71, gbv_min = 0.05,                &
      ! CABLE constants: minimum wind at zref for u* (m s-1), Pohlhausen
      ! coefficient (single-sided plate), Prandtl number, minimum
      ! boundary-layer conductance of the top leaf (mol m-2 s-1).
  eps_leaf = 0.96,                                                             &
      ! Leaf emissivity (Wang & Leuning, 1998, Table 3).
  eps_soil = 0.94,                                                             &
      ! Soil emissivity (Wang & Leuning, 1998, Table 3).
  dheat = 21.5e-6,                                                             &
      ! Molecular diffusivity of heat (m2 s-1).
  gbh_2_gbw = 1.075,                                                           &
      ! Ratio of boundary-layer conductances for water vapour and heat.
  m_h2o = 0.018015,                                                            &
      ! Molar mass of water (kg mol-1).
  u_min = 0.1,                                                                 &
      ! Minimum wind speed for forced convection (m s-1).
  dt_leaf = 15.0, dt_can = 10.0
      ! Limits on leaf - canopy air and canopy air - tair (K).

INTEGER, PARAMETER :: n_solve = 8
      ! Maximum evaluations of the canopy balance in the t_c search.
REAL(KIND=real_jlslsm), PARAMETER :: tol_tc = 0.01
      ! Tolerance on t_c (K).

INTEGER, PARAMETER :: n_mo_cable = 4
      ! Monin-Obukhov iterations for leaf_aero_model = 1 (CABLE niter).

REAL(KIND=real_jlslsm) ::                                                      &
  t_lo(veg_pts), t_hi(veg_pts), r_lo(veg_pts), r_hi(veg_pts), t_try(veg_pts),  &
  t_try1(veg_pts), p_v(veg_pts), qs0_v(veg_pts), qs1_v(veg_pts),               &
  ts_v(veg_pts), tsh_v(veg_pts), qc_v(veg_pts), h_v(veg_pts), le_v(veg_pts),   &
  hb_v(veg_pts), leb_v(veg_pts), res, ra_m(land_pts),                          &
  u_m(land_pts), us, rt1usc, zr
LOGICAL :: l_done(veg_pts)
INTEGER :: l, m, it, side(veg_pts), it_mo, n_mo

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
    CALL canopy_balance( l, t_try(m), qs0_v(m), qs1_v(m), res, ts_v(m),        &
                         tsh_v(m), qc_v(m), hb_v(m), leb_v(m) )
    IF ( it == 1 ) THEN
      r_lo(m) = res
      ! Root below the bracket: keep the lower end.
      l_done(m) = ( res >= 0.0 )
      IF ( l_done(m) ) CALL set_out( l, m, t_try(m), ts_v(m), tsh_v(m), qc_v(m) )
    ELSE
      r_hi(m) = res
      IF ( .NOT. l_done(m) .AND. res <= 0.0 ) THEN
        ! Root above the bracket: keep the upper end.
        l_done(m) = .TRUE.
        CALL set_out( l, m, t_try(m), ts_v(m), tsh_v(m), qc_v(m) )
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
    CALL canopy_balance( l, t_try(m), qs0_v(m), qs1_v(m), res, ts_v(m),        &
                         tsh_v(m), qc_v(m), hb_v(m), leb_v(m) )
    IF ( ABS(res) < tol_tc .OR. it == n_solve ) THEN
      l_done(m) = .TRUE.
      CALL set_out( l, m, t_try(m), ts_v(m), tsh_v(m), qc_v(m) )
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
! heat and the buoyancy of the latent heat, -k g zref (H + 0.07 LE) /
! (rho cp T u*^3). Only the leaves' fluxes: the soil's are not known here.
! Also after the last iteration, for the next pass to start from.
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

CONTAINS

SUBROUTINE set_out( l, m, tc, ts, tsh, qc )
INTEGER, INTENT(IN) :: l, m
REAL(KIND=real_jlslsm), INTENT(IN) :: tc, ts, tsh, qc
t_c(l)   = tc
t_sun(l) = ts
t_shd(l) = tsh
q_c(l)   = qc
h_v(m)   = hb_v(m)
le_v(m)  = leb_v(m)
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

SUBROUTINE canopy_balance( l, tc, qs0, qs1, res, ts_new, tsh_new, qc_new,      &
                           h_out, le_out )

! Leaf temperatures and canopy humidity for canopy air temperature tc, and
! the residual of the canopy air sensible heat balance,
!   res = tc - tair - ra (H_sun + H_shd) / (rho cp).

INTEGER, INTENT(IN) :: l
REAL(KIND=real_jlslsm), INTENT(IN) :: tc, qs0, qs1
REAL(KIND=real_jlslsm), INTENT(OUT) :: res, ts_new, tsh_new, qc_new, h_out,  &
                                       le_out

REAL(KIND=real_jlslsm) ::                                                      &
  cp_mol, cmolar, lhv, slope, gamma_p, sig_t4, eps_air, grn1, gbhw, rho, c_q,  &
  rn(2), gh(2), gbh(2), a(2), b(2), le, dt, t_x(2), h_tot, lai_x(2), gl_x(2),  &
  t_old(2), frad(2), ff(2), gbhf, gsw, gbw, gw, den, grashof, ra_x, visc
INTEGER :: ic

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

lai_x = [ lai_sun(l), lai_shd(l) ]
gl_x  = [ gl_sun(l), gl_shd(l) ]
t_old = [ t_sun(l), t_shd(l) ]
frad  = [ frad_sun(l), frad_shd(l) ]
ff    = [ ff_sun(l), ff_shd(l) ]
rn(1) = sw_sun(l) + eps_leaf * (lw_down(l) - sig_t4) * wlw_sun(l)              &
        - sig_t4 * (1.0 - eps_soil) * (eps_leaf - eps_air) * wlws_sun(l)
rn(2) = sw_shd(l) + eps_leaf * (lw_down(l) - sig_t4) * wlw_shd(l)              &
        - sig_t4 * (1.0 - eps_soil) * (eps_leaf - eps_air) * wlws_shd(l)

! LE of each class is a + b (qs(tc) - q_c) (W m-2).
DO ic = 1,2
  a(ic)   = 0.0
  b(ic)   = 0.0
  gh(ic)  = 0.0
  gbh(ic) = 0.0
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
    IF ( gsw > 0.0 .AND. gbw > 0.0 ) THEN
      gw    = gbw * gsw / ( gbw + gsw )
      den   = slope + gamma_p * gh(ic) / gw
      a(ic) = slope * rn(ic) / den
      b(ic) = gh(ic) * cp_mol * pstar(l) / repsilon / den
    END IF
  END IF
END DO

! Canopy humidity: q_c = q1 + c_q sum(LE), c_q = ra / (rho lc).
c_q    = ra_x / ( rho * lc )
qc_new = ( q1(l) + c_q * ( a(1) + a(2) + ( b(1) + b(2) ) * qs0 ) )             &
         / ( 1.0 + c_q * ( b(1) + b(2) ) )

h_tot  = 0.0
le_out = 0.0
DO ic = 1,2
  IF ( lai_x(ic) > EPSILON(0.0) ) THEN
    le     = MAX(0.0, a(ic) + b(ic) * (qs0 - qc_new))
    le_out = le_out + le
    dt     = ( rn(ic) - le ) / ( cp_mol * gh(ic) )
    dt     = MIN(MAX(dt, -dt_leaf), dt_leaf)
    t_x(ic) = tc + dt
    h_tot  = h_tot + 2.0 * gbh(ic) * cp_mol * dt
  ELSE
    t_x(ic) = tc
  END IF
END DO

ts_new  = t_x(1)
tsh_new = t_x(2)
res     = tc - tair(l) - ra_x * h_tot / ( rho * cp )
h_out   = h_tot

END SUBROUTINE canopy_balance

END SUBROUTINE leaf_temp_update

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_responses( ft, land_pts, veg_pts, veg_index,              &
                                pft_photo_model, t_leaf, oa, acr, actj, actv,  &
                                dsj, dsv, ccp, kc, ko, km, denom, qtenf_term,  &
                                vcmax_temp, jmax_temp, i2, f_vc )

! Temperature responses of the photosynthesis parameters at leaf
! temperature t_leaf: denom, qtenf_term, ccp, kc, ko (Collatz; C3 only for
! ccp, kc, ko) or ccp, km, i2, vcmax_temp, jmax_temp (Farquhar). Called with
! tstar for the canopy, and with each leaf's temperature when l_leaf_temp.

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE jules_vegetation_mod, ONLY: photo_collatz, photo_farquhar
USE pftparm, ONLY: alpha_elec, c3, deact_jmax, deact_vcmax, q10_leaf, tlow,    &
                   tupp

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts),             &
                       pft_photo_model

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  t_leaf(land_pts),                                                            &
      ! Leaf temperature (K).
  oa(land_pts), acr(land_pts),                                                 &
      ! O2 partial pressure (Pa); absorbed PAR (mol m-2 s-1).
  actj(land_pts), actv(land_pts), dsj(land_pts), dsv(land_pts)
      ! Activation energies and entropy terms (Farquhar only).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ccp(land_pts), kc(land_pts), ko(land_pts), km(land_pts),                     &
  denom(land_pts), qtenf_term(land_pts), vcmax_temp(land_pts),                 &
  jmax_temp(land_pts), i2(land_pts)
      ! IN OUT: each model sets only its own terms (and ccp stays zero for
      ! C4 with Collatz).

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: f_vc(land_pts)
      ! Soil-water capacity factor (l_som_vcmax_psi, vcmax_psi_factor):
      ! Vcmax and Jmax (Farquhar) are scaled by it.

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  t_ref = zerodegc + 25.0,                                                     &
    ! Reference temperature (K).
  tref_rmol = t_ref * rmol
    ! The product of t_ref and rmol (J mol-1).

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: jmax_numerator, kc_val, ko_val, power, tau, tdegc,  &
                          t_minus_ref, t_term, vcmax_numerator

SELECT CASE ( pft_photo_model )

CASE ( photo_collatz )
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l,m,power,tau,tdegc)   &
!$OMP SHARED(c3, veg_pts, veg_index, ccp, denom, ft, kc, ko, oa, tlow,         &
!$OMP        q10_leaf,  qtenf_term, t_leaf, tupp) SCHEDULE(STATIC)
  DO m  = 1,veg_pts
    l = veg_index(m)
    tdegc         = t_leaf(l) - zerodegc
    power         = 0.1 * (tdegc- 25.0)
    denom(l)      = (1.0 + EXP (0.3 * (tdegc - tupp(ft)))) *                   &
                    (1.0 + EXP (0.3 * (tlow(ft) - tdegc)))
    qtenf_term(l) = q10_leaf(ft)** power

    IF ( c3(ft) == 1 ) THEN
      ! Calculate terms that are only needed for C3 plants.
      ! Although oa, kc and ko are always used together we keep them separate
      ! to maintain bit comparability.
      tau    = 2600.0  * (0.57 ** power)
      ccp(l) = 0.5 * oa(l) / tau
      kc(l)  = 30.0    * (2.1 ** power)
      ko(l)  = 30000.0 * (1.2 ** power)
    END IF

  END DO
!$OMP END PARALLEL DO

CASE ( photo_farquhar )
!$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE)                                &
!$OMP PRIVATE(l, m, jmax_numerator, kc_val, ko_val, t_minus_ref, t_term,       &
!$OMP         vcmax_numerator)                                                 &
!$OMP SHARED(veg_pts, veg_index, acr, actj, actv, alpha_elec,                  &
!$OMP        ccp, deact_jmax, deact_vcmax, dsj, dsv, ft, i2, jmax_temp, km,    &
!$OMP        oa, t_leaf, vcmax_temp) SCHEDULE(STATIC)
  DO m = 1,veg_pts

    l = veg_index(m)
    ! Temperature responses of carboxylation, oxygenation,and CO2 compensation
    ! point, from Bernacchi et al. (2001).
    t_minus_ref = t_leaf(l) - t_ref
    t_term      = t_minus_ref / ( tref_rmol * t_leaf(l) )
    ccp(l)      = 4.73078 * EXP( 37830.0 * t_term )
    ! For the Farquhar model we combine oa, kc and ko into km.
    kc_val      = 44.8    * EXP( 79430.0 * t_term )
    ko_val      = 30808.2 * EXP( 36380.0 * t_term )
    km(l)       = kc_val * ( 1.0 + oa(l) / ko_val )
    ! Radiation that goes to Photosystem II.
    i2(l)       = alpha_elec(ft) * acr(l)

    ! Calculate the temperature response of Vcmax and Jmax, Eq.17 of
    ! Medlyn et al. (2002).
    vcmax_numerator = EXP( actv(l) * t_minus_ref                               &
                           / ( tref_rmol * t_leaf(l) ) )                       &
                      * ( 1.0 + EXP( ( t_ref * dsv(l) - deact_vcmax(ft) )      &
                                     / tref_rmol ) )
    vcmax_temp(l)   = vcmax_numerator                                          &
                      / ( 1.0 + EXP( ( t_leaf(l) * dsv(l) - deact_vcmax(ft) )  &
                                     / ( t_leaf(l) * rmol ) ) )

    jmax_numerator = EXP( actj(l) * t_minus_ref                                &
                           / ( tref_rmol * t_leaf(l) ) )                       &
                      * ( 1.0 + EXP( ( t_ref * dsj(l) - deact_jmax(ft) )       &
                                     / tref_rmol ) )
    jmax_temp(l)   = jmax_numerator                                            &
                     / ( 1.0 + EXP( ( t_leaf(l) * dsj(l) - deact_jmax(ft) )    &
                                    / ( t_leaf(l) * rmol ) ) )

  END DO
!$OMP END PARALLEL DO

  ! Soil-water down-regulation of photosynthetic capacity (l_som_vcmax_psi;
  ! see vcmax_psi_factor). Applied to Vcmax and Jmax, so Rd follows.
  IF ( PRESENT(f_vc) ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      vcmax_temp(l) = vcmax_temp(l) * f_vc(l)
      jmax_temp(l)  = jmax_temp(l)  * f_vc(l)
    END DO
  END IF

END SELECT

END SUBROUTINE leaf_temp_responses

!#############################################################################
!#############################################################################

SUBROUTINE vcmax_psi_factor( ft, land_pts, veg_pts, veg_index, psi_rz, f_vc )

! Soil-water down-regulation factor of photosynthetic capacity
! (l_som_vcmax_psi): the Zhou et al. (2013) form on the root-zone (predawn)
! water potential, as used by De Kauwe et al. (2015, Biogeosciences 12:
! 7503-7518, Eq. 4),
!   f = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi))), psi in MPa,
! f = 1 in wet soil and about 0.5 at psi = psi_f.
! Adjustments to De Kauwe et al. (2015), which had neither (both added
! Oct 2026 after FR-Pue tests, where f -> 0 drove drought GPP to ~0 and
! capacity recovered within the timestep of the first rain):
! (1) an optional floor,
!   f_vc = fmin + (1 - fmin) f   (psi_vcmax_fmin, default 0),
! so the stomatal (hydraulic) limitation, not capacity, sets photosynthesis
! in the driest soil.
! (2) a recovery lag: the factor follows f at once as the soil dries
! (psi_rz itself changes slowly) but relaxes back towards it with a fixed
! 5-day e-folding
! time (90 % back in about 11 days) as the soil rewets, so capacity recovers
! over days rather than within the timestep of the rain. 5 days sits between
! fast reversible downregulation (Rubisco activation, mesophyll conductance:
! hours to days) and protein resynthesis (a week or more), weighted towards
! the slow end because FR-Pue droughts are long and severe. 3 days changes
! little (FR-Pue 2005-07 GPP RMSE 0.86 vs 0.85; instant 0.93), so this is a
! fixed choice, not a parameter. The state is not in the dump (it starts at
! the instant value).

USE pftparm, ONLY: psi_vcmax_f, sf_vcmax, psi_vcmax_fmin
USE jules_surface_types_mod, ONLY: npft
USE timestep_mod, ONLY: timestep

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_rz(land_pts)
    ! Root-zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: f_vc(land_pts)
    ! Capacity factor (-); only veg points are set.

REAL(KIND=real_jlslsm), PARAMETER :: tau_rec = 5.0 * 86400.0
    ! Recovery e-folding time (s).

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: psi_mpa, f_psi, f_now, w_rec

IF ( .NOT. ALLOCATED(f_vcmax_state) ) THEN
  ALLOCATE( f_vcmax_state(land_pts, npft) )
  f_vcmax_state(:,:) = -1.0
END IF
w_rec = 1.0 - EXP( -timestep / tau_rec )

DO m = 1,veg_pts
  l = veg_index(m)
  psi_mpa = MIN(psi_rz(l), 0.0) * 1.0e-6
  f_psi = ( 1.0 + EXP( sf_vcmax(ft) * psi_vcmax_f(ft) * 1.0e-6 ) )            &
          / ( 1.0 + EXP( sf_vcmax(ft) * ( psi_vcmax_f(ft) * 1.0e-6             &
                                          - psi_mpa ) ) )
  f_psi = MIN(MAX(f_psi, 0.0_real_jlslsm), 1.0_real_jlslsm)
  f_now = psi_vcmax_fmin(ft) + (1.0 - psi_vcmax_fmin(ft)) * f_psi
  IF ( f_vcmax_state(l,ft) < 0.0 .OR. f_now <= f_vcmax_state(l,ft) ) THEN
    f_vcmax_state(l,ft) = f_now
  ELSE
    f_vcmax_state(l,ft) = f_vcmax_state(l,ft)                                  &
                          + w_rec * ( f_now - f_vcmax_state(l,ft) )
  END IF
  f_vc(l) = f_vcmax_state(l,ft)
END DO

END SUBROUTINE vcmax_psi_factor

!#############################################################################
!#############################################################################

SUBROUTINE calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,      &
                                  veg_index, denom, jmax_temp, jv25,           &
                                  nleaf, qtenf_term, vcmax_temp,               &
                                  jmax, rd_dark, vcmax )

! Calculate the maximum rates of carboxylation of Rubisco and electron
! transport, and dark respiration without light inhibition.

USE jules_vegetation_mod, ONLY:                                                &
! imported parameters
    jv_ntotal, jv_scale, photo_collatz, photo_farquhar, photo_johnson,         &
! imported scalars that are not changed
    n_alloc_jmax, n_alloc_vcmax, l_trait_phys, photo_jv_model, photo_model

USE pftparm, ONLY: alpha_elec, fd, jv25_ratio, neff, vint, vsl

USE jb_photo_mod, ONLY: jb_vqmax25

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::                                                          &
  ft,                                                                          &
    ! Index of plant functional type.
  land_pts,                                                                    &
    ! Number of land points.
  pft_photo_model,                                                             &
    ! Indicates which photosynthesis model to use for the current PFT.
  veg_pts,                                                                     &
    ! Number of vegetated points.
  veg_index(land_pts)
    ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  denom(land_pts),                                                             &
    ! Denominator in equation for Vcmax with the Collatz model.
  jmax_temp(land_pts),                                                         &
    ! Factor expressing the effect of temperature on Jmax.
    ! Only used with the Farquhar model.
  jv25(land_pts),                                                              &
    ! Ratio of Jmax to Vcmax at 25 degC, including any acclimation.
    ! Only used with the Farquhar model.
  nleaf(land_pts),                                                             &
    ! Leaf nitrogen concentration.
    ! If l_trait_phys = (kg N m-2),  else = (kgN [kgC]-1).
  qtenf_term(land_pts),                                                        &
   ! Q10 temperature term used for Vcmax with the Collatz model.
  vcmax_temp(land_pts)
    ! Factor expressing the effect of temperature on Vcmax.
    ! Only used with the Farquhar model.

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  jmax(land_pts),                                                              &
    ! Maximum rate of electron transport (mol CO2 m-2 s-1).
    ! Only calculated with the Farquhar model. With photo_johnson, the
    ! maximum Cyt b6f activity Vqmax (mol e- m-2 s-1).
  rd_dark(land_pts),                                                           &
    ! Dark respiration before light inhibition (mol CO2/m2/s).
  vcmax(land_pts)
    ! Maximum rate of carboxylation of Rubisco (mol CO2/m2/s).

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  l, m
    ! Indices.

REAL(KIND=real_jlslsm) ::                                                      &
  n_total,                                                                     &
    ! Total N allocated to photosynthetic components (kg m-2).
  recip_j,                                                                     &
    ! Reciprocal of n_alloc_jmax (kg m-2 of N [mol CO2 m-2 s-1]).
  recip_v
    ! Reciprocal of n_alloc_vcmax (kg m-2 of N [mol CO2 m-2 s-1]).

!-----------------------------------------------------------------------------
! Local array variables.
!-----------------------------------------------------------------------------
REAL ::                                                                        &
  vcmax_ref(land_pts)
    ! Maximum rate of carboxylation of Rubisco at the reference temperature,
    ! ignoring the effects of acclimation and N allocation (mol CO2/m2/s).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_PHOTO_PARAMETERS'

!-----------------------------------------------------------------------------
!end of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Calculate some constants.
!-----------------------------------------------------------------------------
IF ( photo_jv_model == jv_ntotal ) THEN
  recip_j  = 1.0 / n_alloc_jmax
  recip_v  = 1.0 / n_alloc_vcmax
END IF

!-----------------------------------------------------------------------------
! Calculate Vcmax at the reference temperature, without any acclimation.
!-----------------------------------------------------------------------------
!$OMP PARALLEL IF(veg_pts > 1)  DEFAULT(NONE)                                  &
!$OMP PRIVATE(l, m, n_total)                                                   &
!$OMP SHARED(ft, pft_photo_model, photo_jv_model, veg_index, veg_pts,          &
!$OMP        photo_model, alpha_elec,                                          &
!$OMP        denom, fd, jmax, jmax_temp, jv25, jv25_ratio, neff, nleaf,        &
!$OMP        qtenf_term, rd_dark, recip_j, recip_v, vcmax, vcmax_ref,          &
!$OMP        vcmax_temp, vint, vsl, l_trait_phys )

!$OMP DO SCHEDULE(STATIC)
DO m = 1,veg_pts

  l = veg_index(m)

  IF (l_trait_phys) THEN
    vcmax_ref(l) = (vsl(ft) * nleaf(l) + vint(ft)) * 1.0e-6  ! Kattge 2009
  ELSE
    vcmax_ref(l) = neff(ft) * nleaf(l)
  END IF

END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------------
! Calculate Vcmax and Jmax.
!-----------------------------------------------------------------------------
SELECT CASE ( pft_photo_model )
CASE ( photo_collatz )

  !---------------------------------------------------------------------------
  ! Use the Collatz model.
  !---------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)
    ! Using brackets here to recreate existing results.
    vcmax(l) = ( vcmax_ref(l) * qtenf_term(l) ) / denom(l)
  END DO
!$OMP END DO NOWAIT

CASE ( photo_farquhar )

  !---------------------------------------------------------------------------
  ! Use the Farquhar model (for C3 plants).
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! Calculate values at the reference temperature, including any acclimation
  ! of jv25 (but excluding other acclimation terms).
  !---------------------------------------------------------------------------
  SELECT CASE ( photo_jv_model )

  CASE ( jv_scale )
    ! Find J25 by scaling V25.
!$OMP DO SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      vcmax(l) = vcmax_ref(l)
      jmax(l)  = vcmax_ref(l) * jv25(l)
    END DO
!$OMP END DO NOWAIT

  CASE ( jv_ntotal )
    ! Assume the total N allocated to photosynthetic capacity is constant.
!$OMP DO SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      ! Calculate total N allocated to photosynthetic capacity, using the
      ! prescribed parameters at the reference temperature.
      ! This is Eq.5 of Mercado et al. (2018).
      n_total = vcmax_ref(l) * recip_v                                         &
                + vcmax_ref(l) * jv25_ratio(ft) * recip_j
      ! Calculate Vcmax and Jmax at 25degC, including temperature acclimation
      ! of J:V.
      vcmax(l) = n_total / ( recip_v + jv25(l) * recip_j )
      jmax(l)  = n_total / ( recip_v / jv25(l) + recip_j )
    END DO
!$OMP END DO NOWAIT

  END SELECT  !  photo_jv_model

  !---------------------------------------------------------------------------
  ! Johnson-Berry: jmax holds the maximum Cyt b6f activity, Vqmax, found from
  ! Jmax at 25 degC (Lamour et al. 2026, Eqn 11). It takes the temperature
  ! response of Jmax below.
  !---------------------------------------------------------------------------
  IF ( photo_model == photo_johnson ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      jmax(l) = jb_vqmax25( jmax(l), alpha_elec(ft) )
    END DO
!$OMP END DO
  END IF

  !---------------------------------------------------------------------------
  ! Calculate final values, including temperature effect.
  !---------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO m = 1,veg_pts
    l = veg_index(m)

    ! Calculate rates according to acclimated ratio and N allocation to
    ! photosynthesis, and including temperature term.
    ! At present neither acclimation nor N allocation are represented.
    vcmax(l) = vcmax(l) * vcmax_temp(l)
    jmax(l)  = jmax(l)  * jmax_temp(l)

  END DO
!$OMP END DO NOWAIT

END SELECT  !  pft_photo_model

!-----------------------------------------------------------------------------
! Calculate dark respiration. Any effect of light inhibition is added later.
!-----------------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO m = 1,veg_pts
  l = veg_index(m)
  rd_dark(l) = fd(ft) * vcmax(l)
END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_photo_parameters

!#############################################################################
!#############################################################################

SUBROUTINE calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je )

! Calculate the electron flux for the Farquhar model, or for Johnson-Berry
! (photo_model = photo_johnson; jmax is then Vqmax, see jb_photo_mod).

USE jules_vegetation_mod, ONLY: photo_johnson, photo_model,                    &
                                light_curvature => light_curvature_fvcb
USE jb_photo_mod, ONLY: jb_electron_flux
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER,INTENT(IN) ::                                                          &
  land_pts,                                                                    &
    ! Number of land points.
  veg_pts,                                                                     &
    ! Number of vegetated points.
  veg_index(land_pts)
    ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  i2(land_pts),                                                                &
    ! Radiation that goes to Photosystem II, expressed as an electron flux
    ! (mol m-2 s-1).
  jmax(land_pts)
    ! Maximum rate of electron transport (mol CO2 m-2 s-1).

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  je(land_pts)
    ! Electron transport rate (mol m-2 s-1).

!-----------------------------------------------------------------------------
! Local variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
  l, m
    ! Indices.

REAL(KIND=real_jlslsm) ::                                                      &
 recip_denom
   ! The reciprocal of the denominator.


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CALC_ELECTRON_FLUX'

!-----------------------------------------------------------------------------
!end of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( photo_model == photo_johnson ) THEN
  CALL jb_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je )
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!-----------------------------------------------------------------------------
! Calculate a constant.
!-----------------------------------------------------------------------------
recip_denom = 1.0 / ( 2.0 * light_curvature )

!-----------------------------------------------------------------------------
! Calculate electron flux by finding a root of a quadratic equation.
! This is the solution of Eq.4 of Medlyn et al. (2002).
!-----------------------------------------------------------------------------
DO m = 1,veg_pts
  l = veg_index(m)
  je(l)  = ( i2(l) + jmax(l)                                                   &
                    - SQRT( ( i2(l) + jmax(l) )**2                             &
                            - 4.0 * light_curvature * i2(l) * jmax(l) )        &
           ) * recip_denom
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE calc_electron_flux

!-----------------------------------------------------------------------------
! Soil supply limit on transpiration (l_som_supply_limit).
!
! stom_opt_mod costs E_sample = dq * gl * pstar / (repsilon * rmol * T)
! (mol H2O m-2 s-1), so E_sample <= E_supply is the same condition as
! gl <= gl_supply = E_supply * repsilon * rmol * T / (dq * pstar). That
! is applied through gl_max, which every Ci search path already checks,
! so gs, A and E stay consistent: the optimiser picks the most profitable
! state the soil can supply. gl_max <= 0 means "no cap" in stom_opt_mod,
! so a zero supply is passed as TINY (no sample feasible: stomata close).
!-----------------------------------------------------------------------------
SUBROUTINE apply_supply_limit( land_pts, veg_pts, veg_index, e_supply, share,  &
                               dq, tstar, pstar, gl_max_in, gl_max_out )

USE c_rmol, ONLY: rmol
USE planet_constants_mod, ONLY: repsilon
USE jules_vegetation_mod, ONLY: l_som_supply_limit

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts
INTEGER, INTENT(IN) :: veg_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_supply(land_pts),                                                          &
      ! Transpiration the soil can supply (kg m-2 s-1); < 0: no limit.
  share(land_pts),                                                             &
      ! Share of e_supply for this call (1 for the big leaf).
  dq(land_pts),                                                                &
      ! Humidity deficit passed to stom_opt_mod (kg kg-1).
  tstar(land_pts),                                                             &
      ! Leaf temperature passed to stom_opt_mod (K).
  pstar(land_pts),                                                             &
      ! Surface pressure (Pa).
  gl_max_in(land_pts)
      ! gl_max without the supply cap (m s-1); <= 0: no cap.

REAL(KIND=real_jlslsm), INTENT(OUT) :: gl_max_out(land_pts)

REAL(KIND=real_jlslsm), PARAMETER :: m_h2o = 0.018015
      ! Molar mass of water (kg mol-1).

REAL(KIND=real_jlslsm) :: gl_sup
INTEGER :: l, m

gl_max_out(:) = gl_max_in(:)
IF ( .NOT. l_som_supply_limit ) RETURN

DO m = 1,veg_pts
  l = veg_index(m)
  IF ( e_supply(l) >= 0.0 .AND. dq(l) > 0.0 ) THEN
    gl_sup = e_supply(l) * share(l) / m_h2o                                    &
             * repsilon * rmol * tstar(l) / ( dq(l) * pstar(l) )
    IF ( gl_max_in(l) > 0.0 ) gl_sup = MIN(gl_sup, gl_max_in(l))
    gl_max_out(l) = MAX(gl_sup, TINY(1.0_real_jlslsm))
  END IF
END DO

END SUBROUTINE apply_supply_limit

END MODULE sf_stom_mod
