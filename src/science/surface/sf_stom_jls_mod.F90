! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE sf_stom_mod

USE um_types, ONLY: real_jlslsm
USE sf_stom_state_mod, ONLY: f_vcmax_state, cpl_state, sf_stom_state_alloc
                            ! MDK, 22/10/26

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_STOM_MOD'

! f_vcmax_state (l_som_vcmax_psi) and cpl_state (l_som_coupled_e, two-leaf)
! are in sf_stom_state_mod, so they go in the dump.
! MDK, 22/10/26

INTEGER, SAVE, PRIVATE :: n_lt_clamp_warn = 0
INTEGER, PARAMETER, PRIVATE :: max_lt_clamp_warn = 20
    ! l_leaf_temp: warnings issued so far, and the most issued, for leaf
    ! temperatures held at the dt_leaf limit (see t_leaf_clamp).
    ! MDK, 18/10/26

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
,                    e_supply,gsoil_can                                        &
,                    can_rad_mod,ilayers,leaf_flux_mod,faparv                  &
,                    lwp_c                                                     &
,                    el,gpp,npp,resp_p,resp_l,resp_r,resp_w                    &
,                    growth_sug, f_nsc                                         &
,                    n_leaf,n_root,n_stem,lai_bal,gc                           &
,                    fapar_sun,fapar_shd,fsun                                  &
,                    fapar_sun_nir,fapar_shd_nir                               &
,                    flux_o3,fo3,fapar_diag,apar_diag,psi_leaf,cica_ratio      &
,                    leaf_k                                                    &
,                    tair,lw_down,u_wind,z1_wind                               &
,                    t_leaf_sun,t_leaf_shd,t_can,t_leaf_clamp                  &
,                    isoprene,terpene,methanol,acetone                         &
,                    open_index,open_pts                                       &
,                    carbon_gain,hydraulic_cost,gc_som                         &
,                    k_max,conductance_b,conductance_c                         &
,                    psi_leaf_extreme, psi_root_extreme                        &
                     !New arguments replacing USE statements
                     !crop_vars_mod (IN)
,                    dvi_cpft,rootc_cpft)

USE leaf_mod, ONLY: leaf
USE leaf_limits_mod, ONLY: leaf_limits
USE leaf_processes_sox_mod, ONLY: leaf_processes_sox
USE bvoc_emissions_mod, ONLY: bvoc_emissions

USE conversions_mod, ONLY: zerodegc
USE theta_field_sizes, ONLY: t_i_length

USE jules_surface_types_mod, ONLY: nnpft, ncpft, npft

USE pftparm, ONLY:                                                             &
        kmax_pft, conductance_b_pft, conductance_c_pft, gcuticular,            &
        min_gl_pft, lai_ref_k   ! lai_ref_k: MDK, 21/10/26
USE jules_vegetation_mod, ONLY:                                                &
! imported model ids. JBaguley
    leaf_flux_fsmc, leaf_flux_stom_opt,                                        &
! imported parameters
    photo_collatz, photo_farquhar, photo_sox_collatz, photo_johnson,           &
    stomata_medlyn, stomata_sox_analytical, stomata_desica,                    &
    stomata_profit_max,                                                        &
    stomata_sox_opt, stomata_cmax, stomata_cgain,                              &
    stomata_supply_loss, stomata_cap, stomata_pm2,                             &
    photo_adapt, photo_acclim, photo_adapt_acclim,                             &
    photo_act_model, photo_act_pft, photo_act_gb, n_photo_coef,                &
! imported scalars that are not changed
    dsj_coef, dsv_coef, jv25_coef, act_j_coef, act_v_coef,                     &
    l_bvoc_emis, l_fapar_diag, l_trait_phys, l_stem_resp_fix, l_o3_damage,     &
    l_scale_resp_pm, photo_acclim_model, photo_model, stomata_model, l_sugar,  &
    som_gl_max, l_som_supply_limit,                                            &
    l_som_cuticular_loss, l_som_gravity, l_red, l_som_rhizo_series,            &
    l_som_coupled_e, l_som_capacitance_on, som_leaf_resist_frac,               &
    l_som_leak_capacitance_only,                                               &
    l_som_vcmax_psi, kplant_model, kplant_model_ground,                        &
    l_leaf_temp, leaf_temp_iter, l_leaf_temp_gc_eq, leaf_aero_model,           &
    leaf_wind_ext, ml_leaf_temp_model, ml_leaf_temp_layers

USE CN_utils_mod, ONLY:                                                        &
! imported procedures
    get_can_ave_fac, nleaf_from_lai

USE calc_c_comps_triffid_mod, ONLY: lai_bal_from_ht, live_stemwood
                            ! MDK, 20/10/26

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


USE stom_opt_jls_mod, ONLY: stom_opt_mod, som_lambda_eff_pft, stom_opt_at_e,   &
                            e_of_g
                            ! MDK, 23/10/26
USE cuticular_leak_state_mod, ONLY: cut_w_store, cut_psi_stem,                 &
                                    cut_q_store, cut_c_store, cut_m_h2o,       &
                                    cut_state_alloc, cut_psi_unset,            &
                                    cut_w_shed
USE pftparm, ONLY: cap_stem, kcrit_fractional_loss, conductance_b_seg,         &
                   conductance_c_seg
USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: l_som_plant_segments
USE pftparm, ONLY: pft_xylem_impairment_model
USE jules_vegetation_mod, ONLY: xylem_impairment_none, SOX_conductance
USE pftparm, ONLY: pft_conductance_model

USE xylem_impairment_mod, ONLY: update_xylem_impairment
! MDK, 17/10/26
USE jules_vegetation_mod, ONLY: l_triffid, l_triffid_nsc, l_phen_psi
USE triffid_nsc_mod, ONLY: triffid_nsc_daily, triffid_nsc_turgor,              &
                           phen_psi_update
USE jules_vegetation_mod, ONLY: discrimination_model, discrimination_off
USE leaf_c13_discrim_mod, ONLY: d13c_disc_pft, d13c_disc_alloc, d13c_disc_leaf
USE canopy_diag_mod, ONLY: canopy_diag_alloc, gc_stomatal_pft, t_can_air_pft,  &
                           q_can_air_pft, h_leaf_pft, le_leaf_pft,            &
                           h_leaf_pts, le_leaf_pts   ! MDK, 24/10/26


USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls
USE plant_hydraulic_path_mod, ONLY: k_path_zero_flow, kcrit_of_path
USE plant_hydraulic_path_mod, ONLY: set_ksr_path, ksr, set_k_ground_basis,     &
                                    kplant_scale
                                       ! MDK, 21/10/26
USE plant_hydraulic_path_mod, ONLY: set_store_q_path, clear_store_q_path,      &
                                    psi_leaf_segments_jls
USE jules_vegetation_mod, ONLY: som_store_coupling
                                       ! MDK, 24/10/26

USE planet_constants_mod, ONLY: g, r
USE water_constants_mod, ONLY: rho_water, Water_kg_per_mol
USE desica_jls_mod, ONLY: desica_fw, desica_hydraulics, tuzet_fw,              &
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
 ,leaf_flux_mod
!                           !Switch for leaf flux model

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
                            !    (stomata_model = stomata_sox_analytical).
,e_supply(land_pts)                                                            &
                            ! IN Transpiration the soil can supply this
                            !    timestep (kg m-2 s-1); negative = no limit
                            !    (l_som_supply_limit, see physiol)
,gsoil_can(land_pts)
                            ! IN Soil evaporation conductance beneath the
                            !    canopy (m s-1, ground area, as soil_evap
                            !    adds it to gc; see physiol)

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
,t_can(land_pts)                                                               &
                            ! INOUT Sunlit, shaded and canopy (leaf-area
                            ! weighted) leaf temperature used by
                            ! photosynthesis (K): tstar unless l_leaf_temp.
                            ! The input (last timestep's) is the first
                            ! guess with l_leaf_temp.
,t_leaf_clamp(land_pts)
                            ! INOUT l_leaf_temp: fraction of the leaf
                            ! classes (with leaf area) whose temperature was
                            ! held at the dt_leaf limit from the canopy air
                            ! this timestep (0 without l_leaf_temp).
                            ! MDK, 18/10/26

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
,f_nsc(land_pts)                                                               &
                            ! INOUT Non-structural carbohydrate mass fraction
                            !      (kgC/kgC)
,k_max(land_pts)                                                               &
                            ! INOUT Maximum stomatal conductance (m/s).
,conductance_b(land_pts)                                                       &
                            ! INOUT Conductance curve parameter b for each land
                            ! point (Pa): the PFT's (never changed here).
,conductance_c(land_pts)                                                       &
                            ! INOUT Conductance curve parameter c for each land
                            ! point: the PFT's.
,psi_leaf_extreme(land_pts)                                                    &
                            ! INOUT Extreme leaf water potential (Pa)
,psi_root_extreme(land_pts)
                            ! INOUT Extreme root water potential (Pa)

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

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 carbon_gain(land_pts)                                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(land_pts)                                                      &
                            ! Hydraulic cost for each leaf state
,gc_som(land_pts)
                            ! Canopy resistance to H2O (m/s) for the stomatal
                            ! optimisation model.

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
,HC_sun(land_pts)                                                              &
                            ! WORK hydraulic cost for each leaf state
,HC_shd(land_pts)                                                              &
                            ! WORK hydraulic cost for each leaf state
,CG_sun(land_pts)                                                              &
                            ! WORK Carbon gain for each leaf state
,CG_shd(land_pts)                                                              &
                            ! WORK Carbon gain for each leaf state
,psi_stem(land_pts)                                                            &
                            ! WORK Stem water potential (Pa): the stem
                            !      segment's leaf-side end at the optimum (xylem
                            !      impairment driver 4), else psi_src.
,psi_stem_sun(land_pts), psi_stem_shd(land_pts)                                &
                            ! WORK psi_stem of the sunlit / shaded leaves.
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

! MDK, 29/04/26
REAL(KIND=real_jlslsm) :: je_dummy(land_pts), fapar_dummy(land_pts)
REAL(KIND=real_jlslsm) :: lam_sun(land_pts), lam_shd(land_pts)
                            ! Two-leaf: dA/dE at the optimum of the sunlit
                            ! and shaded leaves (lambda_eff of stom_opt_mod),
                            ! for the diagnostic som_lambda_eff_pft.
                            !
                            ! MDK, 23/10/26
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
                            ! fsmc passed to leaf_limits in the multilayer
                            ! stomatal optimisation path: 1.0 everywhere, so
                            ! leaf_limits' fsmc == 0 closure test never
                            ! fires there - see the note at that call.
                            !
                            ! MDK, 24/10/26
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
REAL(KIND=real_jlslsm) :: gl_leaf_cap
                            ! som_gl_max, or 0 (no cap) for the supply-loss
                            ! stomata (stomata_model = 9), which have none.
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
REAL(KIND=real_jlslsm) :: vcmax_sun_ml(land_pts), vcmax_shd_ml(land_pts),      &
                          je_sun_ml(land_pts), je_shd_ml(land_pts),            &
                          rd_sun_ml(land_pts), rd_shd_ml(land_pts),            &
                          gl_sun_cls(land_pts), gl_shd_cls(land_pts),          &
                          ones_c(land_pts)
                            ! Multilayer with l_leaf_temp: the layer's
                            ! Vcmax, electron flux and dark respiration
                            ! (per unit leaf area) of the sunlit / shaded
                            ! leaves at their temperatures; the canopy
                            ! stomatal conductance of each class (m s-1,
                            ! ground), for the leaf energy balance; 1.
                            ! MDK, 18/10/26
INTEGER :: n_pass_ml
                            ! Number of passes of the multilayer
                            ! optimisation: iter, or leaf_temp_iter with
                            ! l_leaf_temp.
                            ! MDK, 18/10/26
REAL(KIND=real_jlslsm) :: n_clamp_lt(land_pts)
                            ! l_leaf_temp: number of leaf classes whose
                            ! temperature was held at dt_leaf from the
                            ! canopy air in the last leaf temperature
                            ! update of the timestep.
                            ! MDK, 18/10/26
REAL(KIND=real_jlslsm) :: b_pft_pts(land_pts), c_pft_pts(land_pts)
                            ! conductance_b_pft(ft)/conductance_c_pft(ft)
                            ! broadcast onto land_pts arrays, for the
                            ! multilayer stomatal optimisation, which does
                            ! not apply xylem impairment.
REAL(KIND=real_jlslsm) :: kmax_ref_bigleaf(land_pts)
                            ! Unimpaired kmax_pft(ft) * LAI for the big-leaf
                            ! stomatal optimisation with xylem impairment.
REAL(KIND=real_jlslsm) :: kmax_ref_sun_2l(land_pts), kmax_ref_shd_2l(land_pts)
                            ! Two-leaf: each class's share of the unimpaired
                            ! canopy kmax (reference for the hydraulic cost
                            ! under xylem impairment, as kmax_ref_bigleaf).
                            !
                            ! MDK, 24/10/26
INTEGER :: open_land_index(land_pts)
                            ! Land point of each open point (i.e.
                            ! veg_index(open_index(:))), for the impairment
                            ! update.
REAL(KIND=real_jlslsm) :: kmax_bigleaf(land_pts), kcrit_bigleaf(land_pts)
                            ! Big leaf: the canopy kmax (k_max * LAI, times
                            ! kplant_scale) and its kcrit, on land_pts for
                            ! stom_opt_mod.
REAL(KIND=real_jlslsm) :: gl_max_eff(land_pts), share_sup(land_pts)
                            ! gl_max for a stom_opt_mod call with the soil
                            ! supply cap applied (apply_supply_limit), and the
                            ! share of e_supply that call may use.
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: gl_cut(land_pts), gl_cut_eff(land_pts),              &
                          kmax_cut(land_pts), kcrit_cut(land_pts)
REAL(KIND=real_jlslsm) :: kmax_ref_cut(land_pts)
                            ! The intact canopy kmax of the store's path
                            ! (mol H2O m-2 ground s-1 Pa-1), the optimiser's
                            ! kmax_ref (cut_k_path).
                            !
                            ! MDK, 23/10/26
                            ! Cuticular loss (l_som_cuticular_loss):
                            ! canopy floor conductance (m s-1) before and
                            ! after the soil-supply cap, and the canopy
                            ! kmax/kcrit (DESICA: to re-solve psi_leaf; the
                            ! stomatal optimisation: kmax_cut is the kmax of
                            ! the store's path, cut_k_path).
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm) :: lai_ref_kg(land_pts)
                            ! kplant_model = 1: reference LAI of the root
                            ! and stem conductances (lai_ref_k, or the
                            ! balanced LAI from canopy height).
                            ! MDK, 21/10/26
INTEGER, PARAMETER :: n_cut = 20
                            ! Floor fluxes tried, E_floor * i / n_cut, to
                            ! find the largest one the xylem can carry
                            ! (DESICA).
LOGICAL :: l_store_cpl
                            ! som_store_coupling = 1 (with l_cap_store): the
                            ! optimiser sees the store, through its flow out
                            ! q_f held fixed in each pass (store_cpl_q).
REAL(KIND=real_jlslsm) :: q_f(land_pts)
                            ! som_store_coupling = 1: the store's flow out
                            ! at the stem node the optimiser sees (mol H2O
                            ! m-2 ground s-1; > 0: it gives, < 0: it refills).
REAL(KIND=real_jlslsm) :: e_sup_use(land_pts)
                            ! The soil supply the two-leaf optimiser is
                            ! capped by (kg m-2 s-1; < 0: no cap): e_supply,
                            ! or with som_store_coupling = 1 e_supply + q_f,
                            ! as the soil supplies only E - q_f.
REAL(KIND=real_jlslsm) :: e_be_2l(land_pts)
                            ! som_store_coupling = 1: the canopy E the
                            ! store's balance is taken at (mol H2O m-2 s-1).
REAL(KIND=real_jlslsm) :: q_be(land_pts), psi_st_be(land_pts),                 &
                          psi_l_be(land_pts)
                            ! som_store_coupling = 1: the store's balance
                            ! over the step at the canopy's E (store_cpl_q):
                            ! its flow out (mol H2O m-2 s-1), the stem node
                            ! psi at the end of the step and the leaf psi
                            ! (Pa).
                            !
                            ! MDK, 24/10/26
LOGICAL :: l_cap_store
                            ! The stomatal optimisation's plant capacitance
                            ! is on (l_som_capacitance_on: l_som_capacitance
                            ! or l_som_cuticular_loss): the store sets the
                            ! transpiration, the leak (if on) drains it, and
                            ! the xylem impairment update is made after the
                            ! leak block, with the store's water potentials.
                            ! (Was l_cut_ximpair.)
                            !
                            ! MDK, 23/10/26
LOGICAL :: lit_cut(land_pts)
                            ! Points where the leak may apply: lit points
                            ! (those of the impairment update) with
                            ! apar >= apar_cut_min.
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm), PARAMETER :: apar_cut_min = 2.0
                            ! Least absorbed PAR (W m-2) for the leak,
                            ! about 5 W m-2 of shortwave: not apar > 0,
                            ! as forcing noise at night (SW up to a few
                            ! W m-2) would turn the leak on in the dark.
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm), PARAMETER :: dq_ratio_cut_max = 3.0
                            ! Most the leak's conductance is scaled up by
                            ! the ratio of the leaf-to-air deficit to dqc,
                            ! so that a dqc near 0 (tstar near the dew
                            ! point) cannot give a burst of transpiration.
                            !
                            ! MDK, 22/10/26
! Plant capacitance and cuticular loss (the store at the stem node).
REAL(KIND=real_jlslsm) :: r_crit
                            ! Fraction of the conductance left at the
                            ! critical loss, 1 - kcrit_fractional_loss (-).
REAL(KIND=real_jlslsm) :: psi_fl_leaf
                            ! Floor on psi_leaf: psi at the critical loss on
                            ! the leaf segment (or PFT) curve (Pa).
REAL(KIND=real_jlslsm) :: psi_fl_stem
                            ! Floor on psi_stem, as psi_fl_leaf on the stem
                            ! segment curve (Pa).
REAL(KIND=real_jlslsm) :: c_cut
                            ! The store's capacitance C (mol H2O m-2 ground
                            ! Pa-1).
REAL(KIND=real_jlslsm) :: psi0_store
                            ! The store's (stem node) psi at the start of the
                            ! step, W / C (Pa).
REAL(KIND=real_jlslsm) :: psi0_stem
                            ! Last step's psi_stem (Pa).
LOGICAL :: l_seg_cut
                            ! The cuticular loss's psi_stem is solved (plant
                            ! segments under the impairment).
                            !
                            ! MDK, 24/10/26
LOGICAL :: l_ml_ximp
                            ! Multilayer stomatal optimisation with xylem
                            ! impairment (not in spin-up).
REAL(KIND=real_jlslsm) :: ml_kfrac(land_pts), kmax_imp_lyr(land_pts)
                            ! Multilayer: the impaired fraction of the
                            ! conductance, k_max / kmax_pft (1 without
                            ! impairment), and this layer's impaired path
                            ! conductance per unit leaf area,
                            !
                            !   kmax_imp_lyr = kmax_leaf_lyr ml_kfrac
                            !
                            ! (mol H2O m-2 leaf s-1 Pa-1).
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm) :: qs_air_cut(land_pts)
                            ! Saturated humidity at the air temperature
                            ! (kg kg-1).
REAL(KIND=real_jlslsm) :: dq_cut(land_pts)
                            ! The leaf-to-air deficit driving the leak
                            ! (kg kg-1).
REAL(KIND=real_jlslsm) :: k_cut_eff
                            ! The whole-plant conductance of the store's path
                            ! (cut_k_path, ksr in series), and in the
                            ! relaxation its root side (mol H2O m-2 ground
                            ! s-1 Pa-1).
REAL(KIND=real_jlslsm) :: q_st
                            ! The store's flow out (mol H2O m-2 ground s-1;
                            ! < 0 it refills).
REAL(KIND=real_jlslsm) :: e_fl
                            ! The largest E that leaves the leaf at
                            ! psi_fl_leaf at the end of the step (mol H2O
                            ! m-2 ground s-1).
REAL(KIND=real_jlslsm) :: psi_tgt
                            ! The step's psi, the optimiser's psi_leaf (Pa).
REAL(KIND=real_jlslsm) :: r_st
                            ! Fraction of the store's departure from its
                            ! steady state left after the step (-).
REAL(KIND=real_jlslsm) :: psi_st
                            ! The store's psi at the end of the step (Pa).
REAL(KIND=real_jlslsm) :: w_new_cut
                            ! The store's water after a fall of C (mol H2O
                            ! m-2 ground).
LOGICAL :: lk_cut(land_pts)
                            ! Points where the leak binds.
LOGICAL :: tgt_on(land_pts)
                            ! Points where the store sets the transpiration
                            ! this step.
LOGICAL :: lit_tgt(land_pts)
                            ! The lit points the store is tried at.
LOGICAL :: ml_tank_pass
                            ! The multilayer canopy's store pass.
INTEGER :: n_tgt
                            ! Number of open points where the store sets the
                            ! transpiration.
INTEGER :: idx_tgt(land_pts)
                            ! Those points, as indices of veg_index.
INTEGER :: n_tank_ml
                            ! Number of store passes of the multilayer canopy
                            ! (0 or 1).
REAL(KIND=real_jlslsm) :: tgt_kr(land_pts)
                            ! Root-side conductance of the store's node (mol
                            ! H2O m-2 ground s-1 Pa-1).
REAL(KIND=real_jlslsm) :: tgt_kl(land_pts)
                            ! Leaf-side conductance of the store's node (mol
                            ! H2O m-2 ground s-1 Pa-1).
REAL(KIND=real_jlslsm) :: e_tgt_sun(land_pts), e_tgt_shd(land_pts)
                            ! Target transpiration of the sunlit and shaded
                            ! classes (mol H2O m-2 ground s-1).
REAL(KIND=real_jlslsm) :: e_star_tgt(land_pts)
                            ! The optimiser's canopy transpiration E* (mol
                            ! H2O m-2 ground s-1).
REAL(KIND=real_jlslsm) :: psi_t_tgt(land_pts)
                            ! The optimiser's (canopy) psi_leaf, the target
                            ! (Pa).
REAL(KIND=real_jlslsm) :: e_new_tgt(land_pts)
                            ! The store's transpiration E' (mol H2O m-2
                            ! ground s-1).
REAL(KIND=real_jlslsm) :: r_tgt(land_pts)
                            ! E' / E* (-).
REAL(KIND=real_jlslsm) :: f_leaf_store
                            ! Share of the plant resistance on the leaf side
                            ! of the store: som_leaf_resist_frac kept within
                            ! [0.01, 0.99], so that neither side has an
                            ! infinite conductance (-). Set once, used by
                            ! store_target_e and the carry loop.
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm), PARAMETER :: dpsi_tgt_min = 1.0e4
                            ! Least drop from the soil to the optimiser's
                            ! leaf psi (Pa, 0.01 MPa) for its secant plant
                            ! conductance E* / (psi_s - psi_t) to be used;
                            ! below it (stomata about shut) the secant is
                            ! round-off (one such step gave a leaf 2.6 MPa
                            ! below its stem under the leak).
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm), PARAMETER :: lai_cut_min = 0.01
                            ! Least LAI (m2 m-2) for the store: below it the
                            ! store is dropped, as psi = W / C has no
                            ! meaning as C goes to 0.
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm) :: e_cut(n_cut, land_pts), psi_cut(n_cut, land_pts),    &
                          k_cut(n_cut, land_pts), el_cut(land_pts)
                            ! Cuticular loss (DESICA): the floor fluxes
                            ! tried and their psi_leaf and k, and the flux
                            ! chosen (mol H2O m-2 s-1).
LOGICAL :: l_cut(land_pts)
                            ! Points where the floor applies.
                            !
                            ! MDK, 24/10/26
REAL(KIND=real_jlslsm) :: dq_opt(land_pts), r_bl_c(land_pts),                  &
                          r_ca_c(land_pts), g_oth_c(land_pts), zero_c(land_pts)
                            ! Coupled transpiration (l_som_coupled_e): the
                            ! humidity deficit, leaf boundary layer and canopy
                            ! air resistances (s/m) and the other leaf class's
                            ! effective conductance (m/s) passed to
                            ! stom_opt_mod (dq_opt the deficit used by the
                            ! optimisation also without it); zeros.
REAL(KIND=real_jlslsm) :: gbw_sun_lt(land_pts), gbw_shd_lt(land_pts),          &
                          ra_ca_lt(land_pts), ra_ca_prev(land_pts)
                            ! Two-leaf, coupled: the leaf boundary layer
                            ! conductance for water of each class (m/s, class
                            ! basis) and the canopy air resistance (s/m), from
                            ! the last leaf energy balance; the last
                            ! timestep's canopy air resistance (-1 if none).
REAL(KIND=real_jlslsm) :: g_soil_c(land_pts)
                            ! Soil evaporation conductance beneath the
                            ! canopy (m/s), which shares the canopy air and
                            ! ra with the leaves.
INTEGER :: i_cut
                            ! Loop counter over the floor fluxes (n_cut).
                            !
                            ! MDK, 24/10/26

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
  f_sun_2l, n_cls_lt,                                                          &
      ! Sunlit fraction of LAI; with ml_leaf_temp_model = 2 the number of
      ! leaf classes with leaf area (MDK, 19/10/26).
  dnw_2l,                                                                      &
      ! Layer integral of exp(-kn_2l*L) over dlai.
  kn_2l
      ! N extinction per unit LAI of the two-leaf classes (= knl).

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
REAL(KIND=real_jlslsm), PARAMETER :: m_h2o_lt = Water_kg_per_mol
      ! Molar mass of water (kg mol-1).
INTEGER :: n_pass_2l
      ! Number of passes of the two-leaf optimisation.
LOGICAL :: l_lt_lyr
      ! l_leaf_temp in the multilayer canopy with ml_leaf_temp_model = 2:
      ! a leaf energy balance for the sunlit and shaded leaf of every layer.
      ! MDK, 19/10/26
INTEGER :: nc_lt
      ! With l_lt_lyr, the number of leaf classes, 2 x ilayers (class 2n-1
      ! the sunlit and 2n the shaded leaf of layer n). MDK, 19/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
  lai_lt_c(:,:), sw_lt_c(:,:), wlw_lt_c(:,:), wlws_lt_c(:,:), frad_lt_c(:,:),  &
  ff_lt_c(:,:), gl_lt_c(:,:), t_lt_c(:,:), gbw_lt_c(:,:)
      ! With l_lt_lyr (land_pts, nc_lt): of each class, the leaf area, the
      ! terms of the leaf energy balance fixed over the passes (see
      ! leaf_temp_setup), the stomatal conductance (m s-1, ground), the leaf
      ! temperature (K) and the boundary-layer conductance for water.
      ! MDK, 19/10/26
LOGICAL :: l_multilayer
                            ! .TRUE. for the multilayer canopy's calls to
                            ! stom_opt_mod (it then uses je_dummy and
                            ! fapar_dummy).
                            !
                            ! MDK, 24/10/26
! veg_pts_index is the trivial (1,2,3,...) map into veg_index that makes
! psi_leaf_jls's veg_index(open_index(j)) indexing select veg_index(j)
! directly, i.e. every vegetated point.
INTEGER :: veg_pts_index(land_pts)



REAL(KIND=real_jlslsm) ::                                                      &
 act_j_tmp(n_photo_coef)                                                       &
   ! Coefficients governigng the acclimation of activation energy for Jmax.
,act_v_tmp(n_photo_coef)
   ! Coefficients governing the acclimation of activation energy for Vcmax.

REAL(KIND=real_jlslsm) ::                                                      &
  d13c_canopy(land_pts),                                                       &
    ! Canopy leaf 13C discrimination (permil), assimilation-weighted over
    ! leaf classes (discrimination_model).
  wa_sun, wa_shd
    ! Gross assimilation of the sun and shade classes, as weights.
LOGICAL :: l_d13c_classes
    ! The discrimination was weighted over leaf classes (can_rad_mod = 7);
    ! otherwise it is taken from the canopy ci.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SF_STOM'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set by the two-leaf canopy (can_rad_mod = 7) at points with assimilation.
d13c_canopy(:) = 0.0

gl_leaf_cap = som_gl_max
zero_c(:) = 0.0
IF ( stomata_model == stomata_supply_loss ) gl_leaf_cap = 0.0
! (Initialised here so that every path, not only two-leaf, has them set.)
gbw_sun_lt(:) = 0.0
gbw_shd_lt(:) = 0.0
ra_ca_lt(:)   = 0.0

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
psi_stem(:) = psi_root_zone(:)
! Plant capacitance (the plant water store) of the stomatal optimisation
! (l_som_capacitance_on: l_som_capacitance or the cuticular loss) moves the
! xylem impairment update
! after it (see below).
!
! MDK, 23/10/26
l_cap_store = l_som_capacitance_on
! The store coupled to the optimiser (som_store_coupling = 1). MDK, 24/10/26
l_store_cpl = l_cap_store .AND. som_store_coupling == 1
psi_stem_sun(:) = psi_root_zone(:)
psi_stem_shd(:) = psi_root_zone(:)
IF ( l_som_gravity .AND. leaf_flux_mod == leaf_flux_stom_opt ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    psi_src(l) = psi_root_zone(l) - rho_water * g * MAX(ht(l), 0.0)
  END DO
  psi_stem(:) = psi_src(:)
  psi_stem_sun(:) = psi_src(:)
  psi_stem_shd(:) = psi_src(:)
END IF

! Soil-water capacity factor of Vcmax/Jmax (l_som_vcmax_psi), once per
! timestep, so the recovery lag advances once whatever the number of leaf
! classes and leaf temperature iterations.
f_vc(:) = 1.0
IF ( l_som_vcmax_psi ) CALL vcmax_psi_factor( ft, land_pts, veg_pts,           &
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
  IF ( stomata_model == stomata_sox_analytical ) THEN
    ! initialise to zero to allow accumulation
    lwp_c(l) = 0.0
  END IF
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------------
! Initialise variables that are accumulated over layers.
! Note that some or all of these variables are also used by the big-leaf
! model (can_rad_mod=1) but do not need to be initialised in that case.
!
! je_dummy and fapar_dummy are passed to stom_opt_mod by every canopy, but
! read only by multilayer; set for all, so no canopy passes them undefined.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
je_dummy(:) = 1.0
fapar_dummy(:) = 1.0

! The stomatal optimisation outputs are not set by every canopy (the fsmc
! paths and, for gc_som, the multilayer canopy), so start them at 0 here.
! gc_som is then set from gc after the canopy switch.
!
! MDK, 24/10/26
el(:)             = 0.0
psi_leaf(:)       = 0.0
leaf_k(:)         = 0.0
carbon_gain(:)    = 0.0
hydraulic_cost(:) = 0.0
gc_som(:)         = 0.0
! Canopy diagnostics of this PFT (canopy_diag_mod). MDK, 24/10/26
CALL canopy_diag_alloc( land_pts )
gc_stomatal_pft(:,ft) = 0.0
t_can_air_pft(:,ft)   = 0.0
q_can_air_pft(:,ft)   = 0.0
h_leaf_pft(:,ft)      = 0.0
le_leaf_pft(:,ft)     = 0.0

SELECT CASE ( can_rad_mod )
CASE ( 4,5,6 )
!$OMP PARALLEL DO IF(land_pts > 1)                                             &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(land_pts,anetc,gc,rdc,rdmean,psi_leaf,el,leaf_k,hydraulic_cost,   &
!$OMP        carbon_gain,kmax_per_lyr,kcrit_per_lyr)                           &
!$OMP PRIVATE(l)                                                               &
!$OMP SCHEDULE(STATIC)
  DO l = 1, land_pts
    anetc(l)  = 0.0
    gc(l)     = 0.0
    rdc(l)    = 0.0
    rdmean(l) = 0.0

    ! gs_opt multi-layer. MDK, 09/06/26
    psi_leaf(l)  = 0.0
    el(l) = 0.0
    leaf_k(l) = 0.0
    hydraulic_cost(l) = 0.0
    carbon_gain(l) = 0.0
    kmax_per_lyr(l) = 0.0
    kcrit_per_lyr(l) = 0.0

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

! Soil evaporation beneath the canopy, at qs(tstar), as the tile flux.
g_soil_c(:) = MAX(gsoil_can(:), 0.0)

! Set the minimum-allowed humidity deficit.
! (The stomatal optimisation keeps the dq_min it had when selected with
! leaf_flux_mod = 2 and stomata_model = 2.)
IF ( ( stomata_model == stomata_medlyn )                                       &
     .OR. ( stomata_model == stomata_sox_analytical )                          &
     .OR. ( stomata_model == stomata_profit_max )                              &
     .OR. ( stomata_model == stomata_sox_opt )                                 &
     .OR. ( stomata_model == stomata_cmax )                                    &
     .OR. ( stomata_model == stomata_cgain )                                   &
     .OR. ( stomata_model == stomata_supply_loss )                             &
     .OR. ( stomata_model == stomata_cap )                                     &
     .OR. ( stomata_model == stomata_pm2 ) ) THEN
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
                            vcmax_temp, jmax_temp, i2,                         &
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
! Leaf temperatures output: tstar unless l_leaf_temp (set in can_rad_mod 5,
! 6 and 7).
! MDK, 18/10/26
!-----------------------------------------------------------------------------
n_clamp_lt(:) = 0.0
ones_c(:)     = 1.0
! MDK, 19/10/26
l_lt_lyr = l_leaf_temp .AND. ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) .AND.  &
           ml_leaf_temp_model == ml_leaf_temp_layers .AND.                     &
           leaf_flux_mod == leaf_flux_stom_opt
nc_lt = 2 * ilayers
DO m = 1,veg_pts
  l = veg_index(m)
  t_leaf_clamp(l) = 0.0
END DO
IF ( .NOT. l_leaf_temp ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    t_leaf_sun(l) = tstar(l)
    t_leaf_shd(l) = tstar(l)
    t_can(l)      = tstar(l)
  END DO
END IF

!-----------------------------------------------------------------------------
! Root and stem conductance per unit ground area (kplant_model = 1): the
! reference LAI is lai_ref_k, or where that is unset the balanced LAI from
! canopy height (as lai_bal below, so it follows TRIFFID growth).
! set_k_ground_basis gives the whole-plant factor kplant_scale on kmax,
! its intact reference and kcrit (applied to the paths below) and the
! segment shares used by the segment solvers. With kplant_model = 0 it
! sets 1.
! MDK, 21/10/26
!-----------------------------------------------------------------------------
IF ( kplant_model == kplant_model_ground ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( lai_ref_k(ft) > 0.0 ) THEN
      lai_ref_kg(l) = lai_ref_k(ft)
    ELSE IF ( l_red ) THEN
      lai_ref_kg(l) = veg_state%lai_bal(l,ft)
    ELSE
      lai_ref_kg(l) = lai_bal_from_ht(ft, ht(l))
    END IF
  END DO
ELSE
  lai_ref_kg(:) = lai(:)
END IF
CALL set_k_ground_basis( ft, land_pts, veg_pts, veg_index, lai, lai_ref_kg )

! Plant capacitance: this step's capacitance of the store,
!
!   C = cap_stem L_ref      (mol H2O m-2 ground Pa-1),
!
! with cap_stem per unit leaf area (mol m-2 leaf Pa-1) and L_ref the leaf
! area the plant conductance is scaled by (lai_ref_kg: LAI, or the
! reference LAI with kplant_model = 1), so that C / k does not change with
! phenology. Where C changes, with the store's water W = C psi:
!
!   C falls (leaf loss): psi is kept, W_new = W C_new / C_old; the shed
!                        tissue takes the rest of the deficit (cut_w_shed),
!                        so shedding does not make psi more negative;
!   C rises (growth):    W is kept: the new tissue is born full.
!
! Set before the canopy fluxes, which depend on the store.
!
! MDK, 23/10/26
! No point has its transpiration set by the store until store_target_e.
tgt_on(:) = .FALSE.
IF ( l_cap_store ) THEN
  CALL cut_state_alloc( land_pts )
  f_leaf_store = MIN(MAX(som_leaf_resist_frac, 0.01), 0.99)
  cut_w_shed(:,ft) = 0.0
  DO m = 1,veg_pts
    l = veg_index(m)
    ! This step's C; last step's is in cut_c_store (0 if not yet set).
    c_cut = MAX(cap_stem(ft) * lai_ref_kg(l), TINY(1.0_real_jlslsm))
    ! C has fallen and the store is set: keep psi, shed the difference.
    IF ( cut_w_store(l,ft) <= 0.0 .AND. cut_c_store(l,ft) > c_cut ) THEN
      w_new_cut = cut_w_store(l,ft) * c_cut / cut_c_store(l,ft)
      cut_w_shed(l,ft)  = (w_new_cut - cut_w_store(l,ft)) * cut_m_h2o
      cut_w_store(l,ft) = w_new_cut
    END IF
    cut_c_store(l,ft) = c_cut
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
!$OMP                    nleaf_top, veg_index, veg_pts)                        &
!$OMP             SCHEDULE(STATIC)
    DO m = 1,veg_pts
      l = veg_index(m)
      gl(l)             = 0.0
      nleaf_layer(l)    = nleaf_top(l) * expkn
      faparv_layer(l,n) = faparv(l,n) * dlai(l)
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

    ! Vulnerability curve: the PFT's, or with xylem impairment the
    ! point's impaired curve (conductance_b, conductance_c), and the
    ! impaired fraction of the plant conductance, k_max / kmax_pft, which
    ! scales each layer's path (kmax_imp_lyr) while the intact path
    ! (kmax_leaf_lyr) is the reference for the hydraulic cost, as the
    ! two-leaf classes. With plant segments each leaf's parallel path is
    ! the root, stem and leaf segments in series (stom_opt_mod), with the
    ! segment memory under impairment model 3, as big-leaf and two-leaf.
    !
    ! MDK, 23/10/26
    l_ml_ximp = .NOT. is_spinup .AND.                                          &
                pft_xylem_impairment_model(ft) /= xylem_impairment_none
    ml_kfrac(:) = 1.0
    IF ( l_ml_ximp ) THEN
      b_pft_pts(:) = conductance_b(:)
      c_pft_pts(:) = conductance_c(:)
      DO m = 1,veg_pts
        l = veg_index(m)
        ml_kfrac(l) = k_max(l) / kmax_pft(ft)
      END DO
    ELSE
      b_pft_pts(:) = conductance_b_pft(ft)
      c_pft_pts(:) = conductance_c_pft(ft)
    END IF

    fsmc_unity(:) = 1.0
    gl_max_lf(:)  = gl_leaf_cap

        ! The humidity-deficit iteration (DO k = 1,iter) runs over the whole
        ! canopy, as in the big-leaf and two-leaf schemes: the deficit at the
        ! leaves is dq / (1 + ra * gc) with gc the canopy conductance (per m2
        ! ground) from the previous pass, the conductance the surface energy
        ! balance (sf_evap) sees. It used to iterate within each layer on that
        ! layer's per-leaf-area gl, so ra * gl ~ 0 and every leaf saw nearly
        ! the full above-canopy deficit. Same number of stom_opt_mod calls;
        ! only the last pass is kept.
        gc_ml_prev(:) = 0.0

        ! l_leaf_temp, ml_leaf_temp_model = 1: a sunlit and a shaded leaf
        ! temperature for the whole canopy, from the two-leaf leaf energy
        ! balance (lt_init, leaf_temp_update) with the multilayer canopy's
        ! sunlit / shaded leaf area, absorbed radiation and stomatal
        ! conductance summed over the layers. Every layer's sunlit (shaded)
        ! leaves are at T_sun (T_shd): their photosynthesis temperature
        ! terms, respiration and humidity deficit qs(T_leaf) - q_c.
        !
        ! ml_leaf_temp_model = 2 (l_lt_lyr): the same balance for the sunlit
        ! and the shaded leaf of every layer (leaf_temp_update_n with 2 x
        ! ilayers classes), each from its own absorbed short- and long-wave
        ! radiation, stomatal conductance and forced-convection boundary
        ! layer at the wind of its depth, all with one canopy air; each
        ! layer's classes use their own temperatures.
        !
        ! Either way the temperatures are updated after each pass but the
        ! last (leaf_temp_iter passes), so there are no more stom_opt_mod
        ! calls than without l_leaf_temp for leaf_temp_iter = iter.
        !
        ! MDK, 19/10/26
        n_pass_ml = iter
        IF ( l_leaf_temp ) THEN
          DO m = 1,veg_pts
            l = veg_index(m)
            lai_sun_2l(l)  = 0.0
            lai_shd_2l(l)  = 0.0
            apar_sun_2l(l) = 0.0
            apar_shd_2l(l) = 0.0
            DO n = 1,ilayers
              lai_sun_2l(l)  = lai_sun_2l(l)  + fsun(l,n) * dlai(l)
              lai_shd_2l(l)  = lai_shd_2l(l)  + (1.0 - fsun(l,n)) * dlai(l)
              apar_sun_2l(l) = apar_sun_2l(l)                                  &
                               + fapar_sun(l,n) * fsun(l,n) * dlai(l) * ipar(l)
              apar_shd_2l(l) = apar_shd_2l(l)                                  &
                               + fapar_shd(l,n) * (1.0 - fsun(l,n)) * dlai(l)  &
                               * ipar(l)
            END DO
          END DO
          IF ( l_lt_lyr ) THEN
            ! ml_leaf_temp_model = 2: the classes of every layer.
            ! MDK, 19/10/26
            ALLOCATE( lai_lt_c(land_pts,nc_lt), sw_lt_c(land_pts,nc_lt),       &
                      wlw_lt_c(land_pts,nc_lt), wlws_lt_c(land_pts,nc_lt),     &
                      frad_lt_c(land_pts,nc_lt), ff_lt_c(land_pts,nc_lt),      &
                      gl_lt_c(land_pts,nc_lt), t_lt_c(land_pts,nc_lt),         &
                      gbw_lt_c(land_pts,nc_lt) )
            gl_lt_c(:,:)  = 0.0
            gbw_lt_c(:,:) = 0.0
          END IF
          CALL lt_init()
          n_pass_ml = leaf_temp_iter
          IF ( l_lt_lyr ) THEN
            ! First guess of every layer: last timestep's sunlit / shaded
            ! canopy temperature (or tstar), as lt_init set.
            ! MDK, 19/10/26
            DO n = 1,ilayers
              t_lt_c(:,2*n-1) = t_sun_lt(:)
              t_lt_c(:,2*n)   = t_shd_lt(:)
            END DO
          END IF
        END IF

        ! Plant capacitance, multilayer canopy: one more pass over the
        ! canopy after the optimiser's passes (n_tank_ml = 1), in which
        !
        !   1. E' comes from the canopy E* (el) and the leaf-area-weighted
        !      psi_leaf of the last pass (store_target_e);
        !   2. each sunlit and shaded leaf of each layer is re-solved at
        !
        !        E_leaf' = E_leaf E' / E*
        !
        !      (ml_class, with ml_tank_pass).
        !
        ! MDK, 23/10/26
        n_tank_ml = 0
        IF ( l_cap_store ) n_tank_ml = 1
        ml_tank_pass = .FALSE.
        DO k = 1,n_pass_ml + n_tank_ml

        ! The store's pass: E' from the totals of the last pass (before
        ! they are reset below).
        IF ( k > n_pass_ml ) THEN
          DO m = 1,veg_pts
            l = veg_index(m)
            e_star_tgt(l) = el(l)
            psi_t_tgt(l)  = psi_src(l)
            IF ( lai(l) > 0.0 ) psi_t_tgt(l) = psi_leaf(l) / lai(l)
            lit_tgt(l)    = apar(l) > 0.0
          END DO
          CALL store_target_e()
          ml_tank_pass = ANY(tgt_on(:))
        END IF

        IF ( l_leaf_temp ) THEN
          rdc_tstar(:)  = 0.0
          gl_sun_cls(:) = 0.0
          gl_shd_cls(:) = 0.0
        END IF

        el(:) = 0.0
        anetl(:) = 0.0
        anetc(:) = 0.0
        gc(:) = 0.0
        rdc(:) = 0.0
        psi_leaf(:) = 0.0
        ! With plant segments, the leaf-area-weighted stem psi of the
        ! leaves' paths (MDK, 23/10/26).
        IF ( l_som_plant_segments ) psi_stem(:) = 0.0
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
        !$OMP             PRIVATE(l, m, kleaf_prof, kleaf_mean, kleaf_a)       &
        !$OMP             SHARED(ft, gl, ilayers, kn, knl, n, nleaf_top,       &
        !$OMP                    nleaf_layer,                                  &
        !$OMP   veg_index,dlai,can_rad_mod, veg_pts, kmax_pft,                 &
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
            !
            ! kmax_per_lyr is per unit leaf area, as vcmax and nleaf_layer
            ! are, so the canopy total (sum of kmax_per_lyr * dlai) is
            ! kmax_pft * LAI. (It was once divided by ilayers and decayed
            ! with leaf N; both removed.)
            !
            ! MDK, 23/10/26
            kmax_per_lyr(l) = kmax_pft(ft)

            ! Scale kcrit by the same fraction of whole-plant kmax that this
            ! layer received, so the critical-conductance threshold stays
            ! calibrated at (1 - kcrit_fractional_loss) of THIS layer's own
            ! capacity (Brodribb & Cochard, 2009; Sabot et al., 2020).
            kcrit_per_lyr(l) = kcrit_of_path(ft, kmax_per_lyr(l))

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
                kleaf_mean = (1.0 - EXP(-REAL(ilayers) * kleaf_a)) /           &
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
            share_sup(l) = (kleaf_prof / kleaf_mean)                           &
                           / MAX(lai(l), EPSILON(1.0))

            END DO
        !$OMP END PARALLEL DO

            ! This layer's impaired path (MDK, 23/10/26).
            DO m = 1,veg_pts
              l = veg_index(m)
              kmax_imp_lyr(l) = kmax_leaf_lyr(l) * ml_kfrac(l)
            END DO

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
        !$OMP PARALLEL DO IF(veg_pts > 1)                                      &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP SHARED(n, veg_pts, veg_index, fapar_sun, fapar_shd, i2_shd,      &
        !$OMP        i2_sun,                                                   &
        !$OMP        ipar, sun_term)                                           &
        !$OMP PRIVATE(l,m)                                                     &
        !$OMP SCHEDULE(STATIC)
            DO m = 1,veg_pts
                l = veg_index(m)
                i2_sun(l) = sun_term * fapar_sun(l,n) * ipar(l)
                i2_shd(l) = sun_term * fapar_shd(l,n) * ipar(l)
            END DO
        !$OMP END PARALLEL DO

            !-------------------------------------------------------------------
            ! Calculate the electron fluxes.
            ! Although we ultimately only need je_sun and je_shd, we also
            ! calculate je because we use that to calculate a single value of
            ! wlite, from which sunlit and shaded terms are calculated using
            ! je_sun and je_shd.
            !-------------------------------------------------------------------
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je)
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_sun,     &
                                     jmax, je_sun)
            CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2_shd,     &
                                     jmax, je_shd)

            ! Calculate ratios of electron flux terms.
            ! These are used to scale the light-limited photosynthesis, which is
            ! linearly related to the electron flux.
        !$OMP PARALLEL DO IF(veg_pts > 1)                                      &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP SHARED(veg_pts, veg_index, je, je_shd, je_shd_ratio, je_sun,     &
        !$OMP        je_sun_ratio)                                             &
        !$OMP PRIVATE(l,m)                                                     &
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

            !-------------------------------------------------------------------
            ! Diagnose the canopy-level humidity deficit from the canopy
            ! conductance of the previous pass (see DO k above).
            ! Initialise the sunlit and shaded respiration rates with the
            ! uninhibited, sunlit respiration rate.
            !-------------------------------------------------------------------
        !$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)             &
        !$OMP          SHARED(dq, dqc, gc_ml_prev, ra, ra_rc, rd_dark, rd_shd, &
        !$OMP                 rd_sun,                                          &
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

            !-------------------------------------------------------------------
            ! Calculate the limiting factors for leaf photosynthesis
            !-------------------------------------------------------------------
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

        !$OMP PARALLEL IF(veg_pts > 1)                                         &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP PRIVATE(m,l)                                                     &
        !$OMP SHARED(veg_pts,veg_index,rd_sun,rd_shd,rd,open_pts,              &
        !$OMP        pft_photo_model,                                          &
        !$OMP        open_index,wlitev_sun,wlite,apar,fapar_sun,ipar,          &
        !$OMP        je_shd_ratio,                                             &
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


            !-------------------------------------------------------------------
            ! Introducing inhibition of leaf respiration in the light for sunlit
            ! and shaded leaves, from papers by Atkin et al. This is an
            ! improvement over the description used for can_rad_mod=4.
            ! This does not change between iterations (though open_index might).
            !-------------------------------------------------------------------
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

            ! Multilayer stomatal optimisation (this branch only runs with
            ! leaf_flux_mod = leaf_flux_stom_opt; the fsmc case is handled
            ! below). (Was a SELECT CASE on leaf_flux_mod whose fsmc and default
            ! branches could not be reached; MDK, 23/10/26.)

                l_multilayer = .TRUE.
                IF ( l_leaf_temp ) THEN
                  ! Each class at its own leaf temperature (see DO k):
                  ! Vcmax, Jmax (and so the electron flux), Rd, the
                  ! Michaelis-Menten terms and the compensation point of
                  ! this layer's leaf N, and the deficit qs(T_leaf) - q_c.
                  ! The electron flux is the class's own, so je_ratio = 1.
                  ! Stem and root respiration stay with Rd at tstar (as
                  ! two-leaf): sum the layer's tstar Rd here.
                  ! MDK, 18/10/26
                  DO m = 1,veg_pts
                    l = veg_index(m)
                    rdc_tstar(l) = rdc_tstar(l)                                &
                                   + ( fsun(l,n) * rd_sun(l)                   &
                                       + (1.0 - fsun(l,n)) * rd_shd(l) )       &
                                   * dlai(l)
                  END DO
                  ! ml_leaf_temp_model = 2: this layer's own sunlit and
                  ! shaded leaf temperatures.
                  ! MDK, 19/10/26
                  IF ( l_lt_lyr ) THEN
                    t_sun_lt(:) = t_lt_c(:,2*n-1)
                    t_shd_lt(:) = t_lt_c(:,2*n)
                  END IF
                  CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,   &
                         pft_photo_model, t_sun_lt, oa, actj, actv, dsj, dsv,  &
                         jv25, nleaf_layer, ones_c, i2_sun, pstar, q_c_lt,     &
                         dq_min, ccp_sun_lt, kc_sun_lt, ko_sun_lt, km_sun_lt,  &
                         vcmax_sun_ml, je_sun_ml, rd_sun_ml, qs_sun_lt,        &
                         dq_sun_lt, f_vc = f_vc )
                  CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,   &
                         pft_photo_model, t_shd_lt, oa, actj, actv, dsj, dsv,  &
                         jv25, nleaf_layer, ones_c, i2_shd, pstar, q_c_lt,     &
                         dq_min, ccp_shd_lt, kc_shd_lt, ko_shd_lt, km_shd_lt,  &
                         vcmax_shd_ml, je_shd_ml, rd_shd_ml, qs_shd_lt,        &
                         dq_shd_lt, f_vc = f_vc )
                  ! Rd at the leaf temperatures, with the light inhibition
                  ! applied above.
                  DO m = 1,veg_pts
                    l = veg_index(m)
                    rd_sun(l) = rd_sun_ml(l)
                    rd_shd(l) = rd_shd_ml(l)
                  END DO
                  DO m = 1,open_pts
                    l = veg_index(open_index(m))
                    IF ( fapar_sun(l,n) * icr(l) * fsun(l,n) * dlai(l) *       &
                         1.0e6 >  10.0 ) rd_sun(l) = 0.7 * rd_sun(l)
                    IF ( fapar_shd(l,n) * icr(l) * (1.0 - fsun(l,n)) *         &
                         dlai(l) * 1.0e6 >  10.0 ) rd_shd(l) = 0.7 * rd_shd(l)
                  END DO

                  ! Sunlit, then shaded leaf of the layer, each at its own
                  ! leaf temperature (MDK, 21/10/26).
                  CALL ml_class( vcmax_sun_ml, kc_sun_lt, ko_sun_lt,           &
                                 ccp_sun_lt, km_sun_lt, dq_sun_lt,             &
                                 je_sun_ml, t_sun_lt, ones_c, fapar_sun(:,n),  &
                                 rd_sun, ci_sun_ml, anetl_sun, el_sun,         &
                                 flux_o3_l_sun, fo3_l_sun, gl_sun,             &
                                 psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun,     &
                                 psi_stem_sun )
                  CALL ml_class( vcmax_shd_ml, kc_shd_lt, ko_shd_lt,           &
                                 ccp_shd_lt, km_shd_lt, dq_shd_lt,             &
                                 je_shd_ml, t_shd_lt, ones_c, fapar_shd(:,n),  &
                                 rd_shd, ci_shd_ml, anetl_shd, el_shd,         &
                                 flux_o3_l_shd, fo3_l_shd, gl_shd,             &
                                 psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd,     &
                                 psi_stem_shd )
                ELSE
                ! share_sup: this layer's kmax share, set with kmax_leaf_lyr.
                ! Added gs opt calls for sunlit and shaded. MDK, 29/04/26
                CALL ml_class( vcmax, kc, ko, ccp, km, dqc, je, tstar,         &
                               je_sun_ratio, fapar_sun(:,n),                   &
                               rd_sun, ci_sun_ml, anetl_sun, el_sun,           &
                               flux_o3_l_sun, fo3_l_sun, gl_sun,               &
                               psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun,       &
                               psi_stem_sun )
                CALL ml_class( vcmax, kc, ko, ccp, km, dqc, je, tstar,         &
                               je_shd_ratio, fapar_shd(:,n),                   &
                               rd_shd, ci_shd_ml, anetl_shd, el_shd,           &
                               flux_o3_l_shd, fo3_l_shd, gl_shd,               &
                               psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd,       &
                               psi_stem_shd )
                END IF

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
                ! The same for open leaves the optimiser closed (no flux,
                ! gl = min_gl_pft from stom_opt_mod).
                !
                ! MDK, 21/10/26
                DO i = 1,open_pts
                    l = veg_index(open_index(i))
                    IF ( lai(l) > EPSILON(0.0) ) THEN
                        IF ( el_sun(l) <= 0.0 .AND.                            &
                             gl_sun(l) <= min_gl_pft(ft) )                     &
                          gl_sun(l) = gl_sun(l) / lai(l)
                        IF ( el_shd(l) <= 0.0 .AND.                            &
                             gl_shd(l) <= min_gl_pft(ft) )                     &
                          gl_shd(l) = gl_shd(l) / lai(l)
                    END IF
                END DO

            !-------------------------------------------------------------------
            ! Update layer conductance.
            !-------------------------------------------------------------------
        !$OMP PARALLEL DO                                                      &
        !$OMP SCHEDULE(STATIC)                                                 &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP PRIVATE(m,l)                                                     &
        !$OMP SHARED(veg_pts,veg_index,gl,fsun,n,gl_sun,gl_shd,rd,rd_sun,rd_shd)
            DO m = 1,veg_pts
                l = veg_index(m)
                gl(l) = fsun(l,n) * gl_sun(l) + (1.0 - fsun(l,n)) * gl_shd(l)
                rd(l) = fsun(l,n) * rd_sun(l) + (1.0 - fsun(l,n)) * rd_shd(l)
            END DO
        !$OMP END PARALLEL DO

            !-------------------------------------------------------------------
            ! Add to canopy-level values.
            !-------------------------------------------------------------------
        !$OMP PARALLEL DO                                                      &
        !$OMP SCHEDULE(STATIC)                                                 &
        !$OMP DEFAULT(NONE)                                                    &
        !$OMP PRIVATE(m,l)                                                     &
        !$OMP SHARED(veg_pts,veg_index,anetl,fsun,anetl_sun,anetl_shd,anetc,   &
        !$OMP        dlai,gc,                                                  &
        !$OMP        gl,rdc,rd,rdmean,ilayers,l_o3_damage,flux_o3_l,           &
        !$OMP        flux_o3_l_sun,                                            &
        !$OMP        flux_o3_l_shd,fo3_l_sun,fo3_l_shd,flux_o3,fo3,lai,n,      &
        !$OMP        fo3_l,                                                    &
        !$OMP        el,el_sun,el_shd,hydraulic_cost,HC_sun,HC_shd,            &
        !$OMP        carbon_gain,                                              &
        !$OMP        CG_sun,CG_shd,psi_leaf,psi_leaf_sun,psi_leaf_shd,leaf_k,  &
        !$OMP        leaf_k_sun,leaf_k_shd,ci_gl_ml,ci_lai_ml,gl_sun,gl_shd,   &
        !$OMP        ci_sun_ml,ci_shd_ml,l_som_plant_segments,psi_stem,        &
        !$OMP        psi_stem_sun,psi_stem_shd)
            DO m = 1,veg_pts
            l = veg_index(m)

            anetl(l)     = fsun(l,n) * anetl_sun(l)                            &
                            + (1.0 - fsun(l,n)) * anetl_shd(l)
            anetc(l)     = anetc(l) + anetl(l) * dlai(l)

            gc(l)        = gc(l)  + gl(l) * dlai(l)
            rdc(l)       = rdc(l) + rd(l) * dlai(l)

            rdmean(l)    = rdmean(l) + rd(l) / REAL(ilayers)

            ! These are ground-area integrated fluxes. MDK, 09/06/26
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
            ! (a sum until it is turned into the canopy weighted mean below).
            psi_leaf(l) = psi_leaf(l) +                                        &
                            (fsun(l,n) * psi_leaf_sun(l) +                     &
                            (1.0 - fsun(l,n)) * psi_leaf_shd(l)) *             &
                            dlai(l)

            leaf_k(l) = leaf_k(l) + (fsun(l,n) * leaf_k_sun(l) +               &
                                    (1.0 - fsun(l,n)) * leaf_k_shd(l)) *       &
                                    dlai(l)

            ! MDK, 23/10/26
            IF ( l_som_plant_segments )                                        &
              psi_stem(l) = psi_stem(l) +                                      &
                            (fsun(l,n) * psi_stem_sun(l) +                     &
                            (1.0 - fsun(l,n)) * psi_stem_shd(l)) * dlai(l)

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

            ! l_leaf_temp: the class conductances (m s-1, ground) for the
            ! leaf energy balance.
            ! MDK, 18/10/26
            IF ( l_leaf_temp ) THEN
              DO m = 1,veg_pts
                l = veg_index(m)
                gl_sun_cls(l) = gl_sun_cls(l) + fsun(l,n) * gl_sun(l) * dlai(l)
                gl_shd_cls(l) = gl_shd_cls(l)                                  &
                                + (1.0 - fsun(l,n)) * gl_shd(l) * dlai(l)
              END DO
              ! (and of each layer's classes; MDK, 19/10/26)
              IF ( l_lt_lyr ) THEN
                DO m = 1,veg_pts
                  l = veg_index(m)
                  gl_lt_c(l,2*n-1) = fsun(l,n) * gl_sun(l) * dlai(l)
                  gl_lt_c(l,2*n)   = (1.0 - fsun(l,n)) * gl_shd(l) * dlai(l)
                END DO
              END IF
            END IF

        END DO                   ! N LAYERS

        ! Canopy conductance of this pass, for the next pass's deficit.
        ! Not after the optimiser's last pass when the store's pass
        ! follows, so that pass solves the same humidity state, and
        ! E_leaf' = E_leaf E' / E* scales the E_leaf behind E*.
        !
        ! MDK, 23/10/26
        IF ( k /= n_pass_ml .OR. n_tank_ml == 0 ) gc_ml_prev(:) = gc(:)

        ! New leaf and canopy air temperatures for the next pass (none after
        ! the last pass, so the temperatures output are those used), as
        ! two-leaf.
        ! MDK, 18/10/26
        ! ml_leaf_temp_model = 2: the balance of every layer's sunlit and
        ! shaded leaf (MDK, 19/10/26).
        IF ( l_lt_lyr .AND. k < n_pass_ml ) THEN
          CALL leaf_temp_update_n( land_pts, veg_pts, veg_index, nc_lt, tair,  &
                                   q1, pstar, ra_lt, lw_down, u_wind,          &
                                   lai_lt_c, sw_lt_c, wlw_lt_c, wlws_lt_c,     &
                                   frad_lt_c, ff_lt_c, z0m_a, zref_a, rt1ab_a, &
                                   zrd_a, usc_a, usuh_a, gl_lt_c, g_soil_c,    &
                                   qs, zeta_lt, t_lt_c, t_c_lt, q_c_lt,        &
                                   gbw_lt_c, ra_ca_lt, n_clamp_lt )
        ELSE IF ( l_leaf_temp .AND. k < n_pass_ml ) THEN
          CALL leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1,       &
                                 pstar, ra_lt, lw_down, u_wind, lai_sun_2l,    &
                                 lai_shd_2l, sw_sun_lt, sw_shd_lt, wlw_sun_lt, &
                                 wlw_shd_lt, wlws_sun_lt, wlws_shd_lt,         &
                                 frad_sun_lt, frad_shd_lt, ff_sun_lt,          &
                                 ff_shd_lt, z0m_a, zref_a, rt1ab_a, zrd_a,     &
                                 usc_a, usuh_a, gl_sun_cls, gl_shd_cls,        &
                                 g_soil_c, qs, zeta_lt, t_sun_lt, t_shd_lt,    &
                                 t_c_lt, q_c_lt, gbw_sun_lt, gbw_shd_lt,       &
                                 ra_ca_lt, n_clamp_lt )
        END IF

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
                IF ( l_som_plant_segments )                                    &
                  psi_stem(l) = psi_stem(l) / lai(l)   ! MDK, 23/10/26
            ELSE IF ( l_som_plant_segments ) THEN
                psi_stem(l) = psi_src(l)               ! MDK, 23/10/26
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

            ! l_leaf_temp: leaf temperatures used and the transpiration-
            ! equivalent gc (as two-leaf), capped at som_gl_max x LAI.
            ! MDK, 18/10/26
            IF ( l_leaf_temp ) THEN
              IF ( lai(l) > EPSILON(0.0) ) THEN
                f_sun_2l = lai_sun_2l(l) / lai(l)
              ELSE
                f_sun_2l = 0.5
              END IF
              gc_cap = gl_leaf_cap * lai(l)
              IF ( l_lt_lyr ) THEN
                ! ml_leaf_temp_model = 2: the sunlit (shaded) temperature
                ! output is the leaf-area weighted mean of the layers'
                ! sunlit (shaded) leaves (so t_can is that of all leaves);
                ! a class without leaf area takes the canopy air's (and
                ! lt_finish then the other class's). The clamp fraction is
                ! over all the classes with leaf area.
                ! MDK, 19/10/26
                t_sun_lt(l) = 0.0
                t_shd_lt(l) = 0.0
                n_cls_lt    = 0.0
                DO n = 1,ilayers
                  t_sun_lt(l) = t_sun_lt(l)                                    &
                                + lai_lt_c(l,2*n-1) * t_lt_c(l,2*n-1)
                  t_shd_lt(l) = t_shd_lt(l)                                    &
                                + lai_lt_c(l,2*n) * t_lt_c(l,2*n)
                  IF ( lai_lt_c(l,2*n-1) > EPSILON(0.0) )                      &
                    n_cls_lt = n_cls_lt + 1.0
                  IF ( lai_lt_c(l,2*n) > EPSILON(0.0) )                        &
                    n_cls_lt = n_cls_lt + 1.0
                END DO
                IF ( lai_sun_2l(l) > EPSILON(0.0) ) THEN
                  t_sun_lt(l) = t_sun_lt(l) / lai_sun_2l(l)
                ELSE
                  t_sun_lt(l) = t_c_lt(l)
                END IF
                IF ( lai_shd_2l(l) > EPSILON(0.0) ) THEN
                  t_shd_lt(l) = t_shd_lt(l) / lai_shd_2l(l)
                ELSE
                  t_shd_lt(l) = t_c_lt(l)
                END IF
                CALL lt_finish( l, f_sun_2l, n_cls_lt )
              ELSE
              CALL lt_finish( l, f_sun_2l )
              END IF
            END IF

        END DO
        ! MDK, 19/10/26
        IF ( l_lt_lyr ) THEN
          DEALLOCATE( lai_lt_c, sw_lt_c, wlw_lt_c, wlws_lt_c, frad_lt_c,       &
                      ff_lt_c, gl_lt_c, t_lt_c, gbw_lt_c )
        END IF

        ! Xylem impairment state update, once per timestep (after the
        ! humidity iteration), as two-leaf: the canopy (leaf-area-weighted)
        ! psi_leaf and, with segments, psi_stem drive one damage state per
        ! point, at the lit points. With the plant capacitance the update
        ! is made after the store, below.
        !
        ! MDK, 23/10/26
        IF ( l_ml_ximp .AND. .NOT. l_cap_store ) THEN
          open_pts = 0
          DO m = 1,veg_pts
            l = veg_index(m)
            IF ( apar(l) > 0.0 ) THEN
              open_pts = open_pts + 1
              open_land_index(open_pts) = l
            END IF
          END DO
          CALL update_xylem_impairment(                                        &
                ! IN
                  land_pts, open_pts, open_land_index, ft, psi_leaf, psi_src,  &
                  psi_stem,                                                    &
                ! IN OUT
                  k_max, psi_leaf_extreme, psi_root_extreme,                   &
                ! IN (recovery)
                  lai, canht, anetc )
        END IF
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
  !$OMP PARALLEL DO IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                &
  !$OMP             SHARED(ft, gl, ilayers, kn, knl, n, nleaf_top,             &
  !$OMP                    nleaf_layer,                                        &
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
      CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,      &
                                  veg_index, denom, jmax_temp, jv25,           &
                                  nleaf_layer, qtenf_term, vcmax_temp,         &
                                  jmax, rd_dark, vcmax )

      IF ( pft_photo_model == photo_farquhar ) THEN
      !-----------------------------------------------------------------------
      ! Calculate sunlit and shaded radiation terms.
      !-----------------------------------------------------------------------
  !$OMP PARALLEL DO IF(veg_pts > 1)                                            &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP SHARED(n, veg_pts, veg_index, fapar_sun, fapar_shd, i2_shd, i2_sun,    &
  !$OMP        ipar, sun_term)                                                 &
  !$OMP PRIVATE(l,m)                                                           &
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
  !$OMP PARALLEL DO IF(veg_pts > 1)                                            &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP SHARED(veg_pts, veg_index, je, je_shd, je_shd_ratio, je_sun,           &
  !$OMP        je_sun_ratio)                                                   &
  !$OMP PRIVATE(l,m)                                                           &
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
  !$OMP PARALLEL IF(veg_pts > 1) DEFAULT(NONE) PRIVATE(l, m)                   &
  !$OMP          SHARED(dq, dqc, gl, ra, ra_rc, rd_dark, rd_shd, rd_sun,       &
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
  ,                       acr, apar, ca, ccp, dqc, fsmc_lim, je, kc, km, ko    &
  ,                       oa, pstar, vcmax                                     &
  ,                       clos_pts, open_pts, clos_index, open_index           &
  ,                       ci, wcarb, wexpt, wlite)

  !$OMP PARALLEL IF(veg_pts > 1)                                               &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP PRIVATE(m,l)                                                           &
  !$OMP SHARED(veg_pts,veg_index,rd_sun,rd_shd,rd,open_pts,pft_photo_model,    &
  !$OMP        open_index,wlitev_sun,wlite,apar,fapar_sun,ipar,je_shd_ratio,   &
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
          IF ( fapar_sun(l,n) * icr(l) * fsun(l,n) * dlai(l) *                 &
              1.0e6 >  10.0 ) rd_sun(l) = 0.7 * rd_sun(l)
          IF ( fapar_shd(l,n) * icr(l) * (1.0 - fsun(l,n)) * dlai(l) *         &
              1.0e6 >  10.0 ) rd_shd(l) = 0.7 * rd_shd(l)
      END DO
  !$OMP END DO
  !$OMP END PARALLEL

      ! Added CASE block to do gs_opt on multi-layer canopy. MDK, 29/04/26
      SELECT CASE ( leaf_flux_mod)
          CASE (leaf_flux_fsmc)
          !---------------------------------------------------------------------
          ! Calculate leaf-level fluxes separately for sunlit and shaded leaves.
          !---------------------------------------------------------------------
          CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model, veg_pts&
      ,                clos_index, open_index, veg_index                       &
      ,                ca, ci, fsmc, o3mol, ra, tstar                          &
      ,                wcarb, wexpt, wlitev_sun, rd_sun                        &
      ,                anetl_sun, flux_o3_l_sun, fo3_l_sun, gl_sun)

          CALL leaf (clos_pts, ft, land_pts, open_pts, pft_photo_model, veg_pts&
      ,                clos_index, open_index, veg_index                       &
      ,                ca, ci, fsmc, o3mol, ra, tstar                          &
      ,                wcarb, wexpt, wlitev_shd, rd_shd                        &
      ,                anetl_shd, flux_o3_l_shd, fo3_l_shd, gl_shd)



      CASE DEFAULT
          errcode = 101  !  a hard error
          CALL ereport(RoutineName, errcode,                                   &
              'leaf_flux_mod should be 1 fsmc regulated or 2 stomatal ' //     &
              'optimisation')
      END SELECT ! leaf_flux_model
      !-----------------------------------------------------------------------
      ! Update layer conductance.
      !-----------------------------------------------------------------------
  !$OMP PARALLEL DO                                                            &
  !$OMP SCHEDULE(STATIC)                                                       &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP PRIVATE(m,l)                                                           &
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
  !$OMP PARALLEL DO                                                            &
  !$OMP SCHEDULE(STATIC)                                                       &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP PRIVATE(m,l)                                                           &
  !$OMP SHARED(veg_pts,veg_index,anetl,fsun,anetl_sun,anetl_shd,anetc,dlai,gc, &
  !$OMP        gl,rdc,rd,rdmean,ilayers,l_o3_damage,flux_o3_l,flux_o3_l_sun,   &
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
          flux_o3_l(l) = fsun(l,n) * flux_o3_l_sun(l)                          &
                      + (1.0 - fsun(l,n)) * flux_o3_l_shd(l)
          fo3_l(l)     = fsun(l,n) * fo3_l_sun(l)                              &
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
        IF (stomata_model == stomata_sox_analytical) THEN
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

        !-----------------------------------------------------------------------
        ! Scale to canopy level.
        !-----------------------------------------------------------------------
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

        !-----------------------------------------------------------------------
        ! Scale to canopy level.
        !-----------------------------------------------------------------------
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

          ! Canopy light for the electron flux (i2 is only set for Farquhar).
          IF ( pft_photo_model == photo_farquhar ) i2c(l) = i2(l) * fpar(l)
        END DO

        ! Recalculate je now that jmaxc has been scaled up, once for all
        ! points: calc_electron_flux loops over every vegetated point, so it
        ! is called after the loop, not inside it. This is done here to keep
        ! the canopy scaling in a similar place to the scaling in the
        ! leaf_flux_fsmc model.
        !
        ! MDK, 24/10/26
        IF ( pft_photo_model == photo_farquhar ) THEN
          CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2c, jmaxc,   &
                                   je )
        END IF

       ! Big leaf: not the multilayer canopy (stom_opt_mod uses je_dummy and
       ! fapar_dummy only for that). MDK, 29/04/26
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
       ! (per m2 ground), so the conductance must be too, or psi_leaf_jls
       ! solves psi_leaf = psi_root - e_leaf/k against too small a k. All
       ! leaves draw on the plant's conductance in parallel, so the canopy
       ! value is kmax_pft * LAI. This used to be kmax_pft * fpar, the
       ! leaf-N (vcmax-like) scaling, which made kmax_pft a top-of-canopy
       ! value and understated canopy conductance by LAI / fpar (1.3-1.7
       ! for LAI 1.5-2.3, kpar 0.5) for a data-derived kmax_pft.
       kmax_bigleaf(:) = kmax_pft(ft) * lai(:)
       kcrit_bigleaf(:) = kcrit_of_path(ft, kmax_bigleaf(:))
       ! som_gl_max is per unit leaf area, so the canopy cap is x LAI, as
       ! kmax_bigleaf (and the multilayer and two-leaf totals).
       gl_max_bigleaf(:) = gl_leaf_cap * lai(:)
       share_sup(:) = 1.0
       ! Coupled transpiration: the deficit to the level-1 air through ra
       ! (shared with the soil evaporation, as in the tile flux).
       dq_opt(:) = dqc(:)
       r_ca_c(:) = 0.0
       IF ( l_som_coupled_e ) THEN
         dq_opt(:) = dq(:)
         r_ca_c(:) = ra(:)
       END IF
       CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,        &
                                share_sup, dq_opt, tstar, pstar,               &
                                gl_max_bigleaf, gl_max_eff,                    &
                                zero_c, r_ca_c, g_soil_c )
       CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

       ! Xylem impairment (JBaguley): k_max/conductance_b/conductance_c hold
       ! the (leaf-basis) impaired vulnerability curve for each point, so
       ! scale the impaired k_max by LAI in exactly the same way as
       ! kmax_pft above, keeping kmax_pft(ft) * LAI as the unimpaired
       ! reference for the hydraulic cost.
       kmax_ref_bigleaf(:) = kmax_bigleaf(:)
       kmax_bigleaf(:) = k_max(:) * lai(:)
       ! In spin-up, or without an impairment model, the intact kmax: the
       ! impairment is off there (stom_opt_mod), so a damaged k_max from a
       ! dump must not act as if it were intact. MDK, 24/10/26
       IF ( is_spinup .OR. pft_xylem_impairment_model(ft) ==                   &
                           xylem_impairment_none )                             &
         kmax_bigleaf(:) = kmax_ref_bigleaf(:)
       ! kplant_model = 1: root and stem per unit ground area, the whole
       ! plant in series. The same factor on the impaired and intact plant
       ! and kcrit, so the damage fraction (k_max / kmax_pft) is unchanged;
       ! the impairment update below uses the unscaled leaf-basis values.
       ! MDK, 21/10/26
       IF ( kplant_model == kplant_model_ground ) THEN
         kmax_bigleaf(:)     = kmax_bigleaf(:)     * kplant_scale(:)
         kmax_ref_bigleaf(:) = kmax_ref_bigleaf(:) * kplant_scale(:)
         kcrit_bigleaf(:)    = kcrit_bigleaf(:)    * kplant_scale(:)
       END IF

       CALL stom_opt_mod (                                                     &
              ! IN
                land_pts, ft, open_pts, open_index,                            &
                pft_photo_model, veg_pts, veg_index,                           &
                ca, psi_src, acrc, apar, oa, vcmaxc, kc, ko, ccp, pstar,       &
                km, dq_opt, je, tstar, je_dummy, fapar_dummy,                  &
                kmax_bigleaf, kcrit_bigleaf, gl_max_eff, ipar,                 &
                l_multilayer, zero_c, r_ca_c, g_soil_c,                        &
                kmax_ref_bigleaf, conductance_b, conductance_c,                &
                psi_leaf_extreme, psi_root_extreme, .TRUE.,                    &
              ! IN OUT
                rdc,                                                           &
              ! OUT
                ci, anetc, el, flux_o3, fo3, gc, psi_leaf,                     &
                carbon_gain, hydraulic_cost, leaf_k,                           &
                psi_stem = psi_stem                                            &
        )

        ! Plant capacitance, big leaf: on the final pass of the humidity
        ! iteration,
        !
        !   1. E' from the optimiser's E* (el) and psi_leaf at the open
        !      points (store_target_e);
        !   2. ci, A (anetc), gc and E re-solved at E' at the points where
        !      the store applies (stom_opt_at_e).
        !
        ! MDK, 23/10/26
        IF ( l_cap_store .AND. k == iter ) THEN
          ! The optimiser's canopy E* and psi_leaf, and the open points.
          DO m = 1,veg_pts
            l = veg_index(m)
            e_star_tgt(l) = el(l)
            psi_t_tgt(l)  = psi_leaf(l)
            lit_tgt(l)    = .FALSE.
          END DO
          DO m = 1,open_pts
            lit_tgt(veg_index(open_index(m))) = .TRUE.
          END DO
          CALL store_target_e()
          ! The open points where the store sets E (indices of veg_index).
          n_tgt = 0
          DO m = 1,open_pts
            IF ( tgt_on(veg_index(open_index(m))) ) THEN
              n_tgt = n_tgt + 1
              idx_tgt(n_tgt) = open_index(m)
            END IF
          END DO
          ! Re-solve the big leaf there at E'.
          IF ( n_tgt > 0 )                                                     &
            CALL stom_opt_at_e( land_pts, ft, n_tgt, idx_tgt, pft_photo_model, &
                                veg_index, rdc, ca, psi_src, acrc, apar, oa,   &
                                vcmaxc, kc, ko, ccp, pstar, km, dq_opt,        &
                                je, tstar, je_dummy, fapar_dummy, ipar,        &
                                kmax_bigleaf, kcrit_bigleaf, gl_max_eff,       &
                                l_multilayer, zero_c, r_ca_c, g_soil_c,        &
                                kmax_ref_bigleaf, conductance_b,               &
                                conductance_c, psi_leaf_extreme,               &
                                psi_root_extreme, .TRUE., e_new_tgt, ci,       &
                                anetc, gc, el )
        END IF

        ! el stays in mol H2O m-2 s-1 (ground), as in every path: physiol uses
        ! it for the layer uptake at the actual root potential
        ! (l_som_rhizo_series). (It used to be converted to kg here, a
        ! leftover with no other consumer.)

        ! leaf_k from stom_opt_mod is canopy scale (see kmax_bigleaf). Return
        ! it per unit leaf area, as the multilayer path does (kplant_pft).
        leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))

        ! Update the impairment state only once per timestep, on the final
        ! pass of the humidity-deficit iteration (DO k = 1,iter): earlier
        ! passes use an unconverged dqc, and the kmax model's weighted update
        ! would otherwise be applied iter times per timestep. JBaguley
        IF (k == iter) THEN

          ! The leaf and root-zone conductances once computed here for the
          ! update were never used by it, so they are no longer computed.
          !
          ! MDK, 22/10/26

          ! Apply xylem imapairment to the vulnerability curve.
          ! NOTE: open_index indexes veg_index, while update_xylem_impairment
          !       indexes the land arrays with it directly, so pass it the
          !       land points of the open stomata.
          DO m = 1, open_pts
            open_land_index(m) = veg_index(open_index(m))
          END DO

          ! With the plant store or the cuticular loss (l_cap_store), the
          ! update is made after it, below.
          !
          ! MDK, 24/10/26
          IF ( .NOT. l_cap_store ) THEN
            CALL update_xylem_impairment(                                      &
                  ! IN
                    land_pts,                                                  &
                    open_pts,                                                  &
                    open_land_index,                                           &
                    ft,                                                        &
                    psi_leaf,                                                  &
                    psi_src,                                                   &
                    psi_stem,                                                  &
                  ! IN OUT
                    k_max,                                                     &
                    psi_leaf_extreme,                                          &
                    psi_root_extreme,                                          &
                  ! IN (recovery)
                    lai,                                                       &
                    canht,                                                     &
                    anetc                                                      &
            )
          END IF

        END IF ! k == iter

      CASE DEFAULT
        errcode = 101  !  a hard error
        CALL ereport(RoutineName, errcode,                                     &
               'leaf_flux_mod should be 1 fsmc regulated or 2 stomatal ' //    &
               'optimisation')
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
      IF ( l_som_cuticular_loss ) gl_cut_ds = gcuticular(ft) * 1.0e-3 * rmol   &
                                               * tstar(l) / pstar(l) * lai(l)
      el_try(l) = e_of_g(MAX(gc(l), gl_cut_ds), MAX(dqc(l), 0.0), pstar(l),    &
                         tstar(l))
    END DO
    CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,        &
                            lai, ht, psi_root_zone, ksr_ds, el_try, .FALSE.,   &
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
  ! light-limited GPP than big-leaf, beyond the effect of the sun/shade
  ! split. Kmax is per unit leaf area as in big-leaf, and the canopy
  ! kmax_pft*LAI is shared between
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
    ! N extinction per unit LAI (extkn of Wang & Leuning 1998): knl, as the
    ! multilayer scheme (can_rad_mod = 6). kpar is the PAR extinction only.
    kn_2l = knl(ft)
    DO n = 1,ilayers
      IF ( kn_2l > EPSILON(0.0) ) THEN
        dnw_2l = ( EXP(-kn_2l * REAL(n-1) * dlai(l))                           &
                   - EXP(-kn_2l * REAL(n) * dlai(l)) ) / kn_2l
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
    ! Xylem impairment: both classes draw on the one (impaired) plant, so
    ! the impaired leaf-basis k_max (= kmax_pft when unimpaired) is shared
    ! the same way, and the unimpaired kmax_pft share is the reference for
    ! the hydraulic cost; kcrit stays on the unimpaired curve, as big-leaf.
    kmax_ref_sun_2l(l) = kmax_pft(ft) * lai(l) * nw_sun_2l(l)                  &
                         / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    kmax_ref_shd_2l(l) = kmax_pft(ft) * lai(l) * nw_shd_2l(l)                  &
                         / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    kmax_sun_2l(l) = k_max(l) * lai(l) * nw_sun_2l(l)                          &
                     / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    kmax_shd_2l(l) = k_max(l) * lai(l) * nw_shd_2l(l)                          &
                     / MAX(nw_sun_2l(l) + nw_shd_2l(l), TINY(1.0))
    ! In spin-up, or without an impairment model, the intact kmax (as big
    ! leaf). MDK, 24/10/26
    IF ( is_spinup .OR. pft_xylem_impairment_model(ft) ==                     &
                        xylem_impairment_none ) THEN
      kmax_sun_2l(l) = kmax_ref_sun_2l(l)
      kmax_shd_2l(l) = kmax_ref_shd_2l(l)
    END IF
    kcrit_sun_2l(l) = kcrit_of_path(ft, kmax_ref_sun_2l(l))
    kcrit_shd_2l(l) = kcrit_of_path(ft, kmax_ref_shd_2l(l))
    ! kplant_model = 1: root and stem per unit ground area, each class's
    ! path in series. One factor for both classes (so the shares are
    ! unchanged) and for the impaired, intact and critical conductances
    ! (so the damage fraction is unchanged); the impairment update divides
    ! it back out.
    ! MDK, 21/10/26
    IF ( kplant_model == kplant_model_ground ) THEN
      kmax_sun_2l(l)     = kmax_sun_2l(l)     * kplant_scale(l)
      kmax_shd_2l(l)     = kmax_shd_2l(l)     * kplant_scale(l)
      kmax_ref_sun_2l(l) = kmax_ref_sun_2l(l) * kplant_scale(l)
      kmax_ref_shd_2l(l) = kmax_ref_shd_2l(l) * kplant_scale(l)
      kcrit_sun_2l(l)    = kcrit_sun_2l(l)    * kplant_scale(l)
      kcrit_shd_2l(l)    = kcrit_shd_2l(l)    * kplant_scale(l)
    END IF
    ! som_gl_max is per unit leaf area: x the class leaf area.
    gl_max_sun_2l(l) = gl_leaf_cap * lai_sun_2l(l)
    gl_max_shd_2l(l) = gl_leaf_cap * lai_shd_2l(l)

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
    ! First guesses and fixed terms of the leaf energy balance (lt_init,
    ! shared with the multilayer scheme).
    ! MDK, 18/10/26
    CALL lt_init()
    n_pass_2l = leaf_temp_iter
    ! Store the tstar-based canopy respiration for the stem and root
    ! respiration, which stay at tstar.
    DO m = 1,veg_pts
      l = veg_index(m)
      rdc_tstar(l) = rd(l) * ( nw_sun_2l(l) + nw_shd_2l(l) )
    END DO
  ELSE
    n_pass_2l = iter
    ! Without leaf temperatures both classes are at tstar, with the tile's
    ! photosynthesis terms, so the class routines take the class arrays
    ! either way (and the 13C discrimination reads set values).
    !
    ! MDK, 24/10/26
    t_sun_lt(:)   = tstar(:)
    t_shd_lt(:)   = tstar(:)
    ccp_sun_lt(:) = ccp(:)
    ccp_shd_lt(:) = ccp(:)
    kc_sun_lt(:)  = kc(:)
    kc_shd_lt(:)  = kc(:)
    ko_sun_lt(:)  = ko(:)
    ko_shd_lt(:)  = ko(:)
    km_sun_lt(:)  = km(:)
    km_shd_lt(:)  = km(:)
  END IF

  !---------------------------------------------------------------------------
  ! Coupled transpiration: first guesses of the class conductances and leaf
  ! boundary layers (the last timestep's) and the canopy air resistance.
  !---------------------------------------------------------------------------
  IF ( l_som_coupled_e ) THEN
    CALL cpl_state_load( ft, land_pts, veg_pts, veg_index, lai, lai_sun_2l,    &
                         lai_shd_2l, gc, gl_sun, gl_shd, gbw_sun_lt,           &
                         gbw_shd_lt, ra_ca_prev )
    IF ( l_leaf_temp .AND. leaf_aero_model == 1 ) THEN
      ! CABLE's rt1: the last timestep's where there is a canopy, else the
      ! neutral first guess set above (JULES's ra where there is no canopy).
      ! (ra_lt is JULES's ra by now, which is the canopy air resistance of
      ! the passes only where there is no canopy.)
      DO m = 1,veg_pts
        l = veg_index(m)
        IF ( zref_a(l) > 0.0 .AND. ra_ca_prev(l) >= 0.0 )                      &
          ra_ca_lt(l) = ra_ca_prev(l)
      END DO
    ELSE IF ( l_leaf_temp ) THEN
      ra_ca_lt(:) = ra_lt(:)
    ELSE
      ra_ca_lt(:) = ra(:)
    END IF
  END IF

  !---------------------------------------------------------------------------
  ! The soil supply cap of the optimiser, and with the store coupled to it
  ! (som_store_coupling = 1) the store's flow out it sees in the first pass:
  ! last step's (0 where the store is not set).
  !
  ! cut_q_store is deliberately not in the dump, so it is 0 after a restart.
  ! That changes nothing for a dump written at midnight: q_f is used only at
  ! lit points, and the night steps set cut_q_store before the first of
  ! them (FR-Pue 2022: a restart with and without it in the dump gave
  ! identical output). Only a dump written in daylight would differ, in its
  ! first step's first pass.
  !
  ! MDK, 24/10/26
  !---------------------------------------------------------------------------
  e_sup_use(:) = e_supply(:)
  IF ( l_store_cpl ) THEN
    q_f(:) = 0.0
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( cut_w_store(l,ft) <= 0.0 .AND. lai(l) >= lai_cut_min )              &
        q_f(l) = cut_q_store(l,ft) / cut_m_h2o
    END DO
    CALL store_cpl_supply()
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
    ! Sunlit, then shaded leaf (the shaded leaf sees this pass's sunlit gl
    ! through coupled_inputs).
    !
    ! MDK, 21/10/26
    !-------------------------------------------------------------------------
    CALL two_leaf_class( apar_sun_2l, lai_sun_2l, kmax_ref_sun_2l,             &
                         qs_sun_lt, dq_sun_lt, gbw_sun_lt, gbw_shd_lt,         &
                         gl_shd, t_sun_lt, gl_max_sun_2l, acr_sun_2l,          &
                         vcmax_sun_2l, kc_sun_lt, ko_sun_lt, ccp_sun_lt,       &
                         km_sun_lt, je_sun, kmax_sun_2l, kcrit_sun_2l,         &
                         rd_sun, ci_sun_2l, anetl_sun, el_sun,                 &
                         flux_o3_l_sun, fo3_l_sun, gl_sun, psi_leaf_sun,       &
                         CG_sun, HC_sun, leaf_k_sun, psi_stem_sun, lam_sun )
    CALL two_leaf_class( apar_shd_2l, lai_shd_2l, kmax_ref_shd_2l,             &
                         qs_shd_lt, dq_shd_lt, gbw_shd_lt, gbw_sun_lt,         &
                         gl_sun, t_shd_lt, gl_max_shd_2l, acr_shd_2l,          &
                         vcmax_shd_2l, kc_shd_lt, ko_shd_lt, ccp_shd_lt,       &
                         km_shd_lt, je_shd, kmax_shd_2l, kcrit_shd_2l,         &
                         rd_shd, ci_shd_2l, anetl_shd, el_shd,                 &
                         flux_o3_l_shd, fo3_l_shd, gl_shd, psi_leaf_shd,       &
                         CG_shd, HC_shd, leaf_k_shd, psi_stem_shd, lam_shd )

    DO m = 1,veg_pts
      l = veg_index(m)
      gc(l) = gl_sun(l) + gl_shd(l)
    END DO

    ! The store coupled to the optimiser (som_store_coupling = 1): the
    ! store's flow out for the next pass, from its balance over the step at
    ! this pass's canopy E (a fixed-point iteration on q_f; the last pass's
    ! balance is committed after the leak block).
    !
    ! MDK, 24/10/26
    IF ( l_store_cpl .AND. k < n_pass_2l ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)
        e_be_2l(l) = MAX(el_sun(l) + el_shd(l), 0.0)
      END DO
      CALL store_cpl_q( .FALSE., e_be_2l )
      DO m = 1,veg_pts
        l = veg_index(m)
        IF ( cut_w_store(l,ft) <= 0.0 .AND. lai(l) >= lai_cut_min )            &
          q_f(l) = q_be(l)
      END DO
      CALL store_cpl_supply()
    END IF

    ! New leaf and canopy air temperatures for the next pass (none after
    ! the last pass, so the temperatures output are those used).
    IF ( l_leaf_temp .AND. k < n_pass_2l ) THEN
      CALL leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             ra_lt, lw_down, u_wind, lai_sun_2l, lai_shd_2l,   &
                             sw_sun_lt, sw_shd_lt, wlw_sun_lt, wlw_shd_lt,     &
                             wlws_sun_lt, wlws_shd_lt, frad_sun_lt,            &
                             frad_shd_lt, ff_sun_lt, ff_shd_lt, z0m_a, zref_a, &
                             rt1ab_a, zrd_a, usc_a, usuh_a, gl_sun,            &
                             gl_shd, g_soil_c, qs, zeta_lt, t_sun_lt,          &
                             t_shd_lt, t_c_lt, q_c_lt, gbw_sun_lt, gbw_shd_lt, &
                             ra_ca_lt, n_clamp_lt )
    END IF

  END DO   ! End of iteration loop

  IF ( l_som_coupled_e ) CALL cpl_state_save( ft, land_pts, veg_pts, veg_index,&
                                              gl_sun, gl_shd, gbw_sun_lt,      &
                                              gbw_shd_lt, ra_ca_lt )

  ! Diagnostic marginal water price of the canopy (som_lambda_eff_pft):
  ! the sunlit and shaded leaves' dA/dE weighted by their transpiration.
  !
  ! MDK, 22/10/26
!$OMP CRITICAL (som_lambda_alloc)
  IF ( .NOT. ALLOCATED(som_lambda_eff_pft) ) THEN
    ALLOCATE( som_lambda_eff_pft(land_pts, npft) )
    som_lambda_eff_pft(:,:) = 0.0
  END IF
!$OMP END CRITICAL (som_lambda_alloc)
  DO m = 1,veg_pts
    l = veg_index(m)
    som_lambda_eff_pft(l,ft) = 0.0
    IF ( el_sun(l) + el_shd(l) > TINY(1.0_real_jlslsm) )                       &
      som_lambda_eff_pft(l,ft) = ( el_sun(l) * lam_sun(l)                      &
                                   + el_shd(l) * lam_shd(l) )                  &
                                 / ( el_sun(l) + el_shd(l) )
  END DO

  ! Plant capacitance: the transpiration that brings the leaf to the
  ! optimiser's psi_leaf at the end of the step through the plant water
  ! store, with A, gl and ci re-solved there (store_target_2l). Not with
  ! the store coupled to the optimiser (som_store_coupling = 1), whose
  ! choice already includes the store.
  !
  ! MDK, 24/10/26
  IF ( l_cap_store .AND. .NOT. l_store_cpl ) CALL store_target_2l()

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
    ! Sun + shade is the canopy conductance (per unit leaf area below, after
    ! the class loop, as the big-leaf path).
    leaf_k(l) = leaf_k_sun(l) + leaf_k_shd(l)

    psi_leaf(l)       = f_sun_2l * psi_leaf_sun(l)                             &
                        + (1.0 - f_sun_2l) * psi_leaf_shd(l)
    psi_stem(l)       = f_sun_2l * psi_stem_sun(l)                             &
                        + (1.0 - f_sun_2l) * psi_stem_shd(l)
    carbon_gain(l)    = f_sun_2l * CG_sun(l) + (1.0 - f_sun_2l) * CG_shd(l)
    hydraulic_cost(l) = f_sun_2l * HC_sun(l) + (1.0 - f_sun_2l) * HC_shd(l)

    IF ( gc(l) > TINY(gc(l)) ) THEN
      ci(l) = ( gl_sun(l) * ci_sun_2l(l) + gl_shd(l) * ci_shd_2l(l) ) / gc(l)
    ELSE
      ci(l) = f_sun_2l * ci_sun_2l(l) + (1.0 - f_sun_2l) * ci_shd_2l(l)
    END IF

    ! Leaf 13C discrimination weighted by the gross assimilation of each
    ! class (ci above is conductance-weighted, which is right for the water
    ! flux but not for the isotopic signal of the fixed carbon).
    IF ( discrimination_model /= discrimination_off ) THEN
      wa_sun = MAX( anetl_sun(l) + rd_sun(l), 0.0 )
      wa_shd = MAX( anetl_shd(l) + rd_shd(l), 0.0 )
      IF ( wa_sun + wa_shd > TINY(wa_sun) ) THEN
        d13c_canopy(l) = ( wa_sun * d13c_disc_leaf( c3(ft) == 1,               &
                             ci_sun_2l(l), ca(l), ccp_sun_lt(l) )              &
                         + wa_shd * d13c_disc_leaf( c3(ft) == 1,               &
                             ci_shd_2l(l), ca(l), ccp_shd_lt(l) ) )            &
                         / ( wa_sun + wa_shd )
      ELSE
        d13c_canopy(l) = 0.0
      END IF
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
      ! Leaf temperatures used and the transpiration-equivalent gc
      ! (lt_finish, shared with the multilayer scheme).
      ! MDK, 18/10/26
      gc_cap = gl_max_sun_2l(l) + gl_max_shd_2l(l)
      CALL lt_finish( l, f_sun_2l )
    END IF
  END DO

  ! As the big-leaf path: el in mol H2O m-2 s-1 and leaf_k per unit leaf
  ! area. (gc_som is set after the canopy switch.)
  !
  ! MDK, 24/10/26
  leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))

  ! Xylem impairment state update, once per timestep (after the humidity
  ! iteration, cf. k == iter in the big-leaf path), driven by the canopy
  ! (leaf-area-weighted) psi_leaf: both classes share the one plant, so
  ! there is a single damage state per point, as for the big leaf.
  ! The class and root-zone conductances once computed here for the update
  ! were never used by it, so they are no longer computed.
  !
  ! MDK, 22/10/26
  open_pts = 0
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( apar_sun_2l(l) > 0.0 .OR. apar_shd_2l(l) > 0.0 ) THEN
      open_pts = open_pts + 1
      open_land_index(open_pts) = l
    END IF
  END DO
  ! With the plant store or the cuticular loss (l_cap_store), the update
  ! is made after it, below.
  !
  ! MDK, 24/10/26
  IF ( .NOT. l_cap_store ) THEN
    CALL update_xylem_impairment(                                              &
          ! IN
            land_pts, open_pts, open_land_index, ft, psi_leaf, psi_src,        &
            psi_stem,                                                          &
          ! IN OUT
            k_max, psi_leaf_extreme, psi_root_extreme,                         &
          ! IN (recovery)
            lai, canht, anetc )
  END IF

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'can_rad_mod should be 1, 4, 5, 6 or 7')

END SELECT  ! can_rad_mod

! Pass the stomatal optimisation conductance out, for every canopy. (The
! cuticular leak and the store below reset it where they change gc.)
!
! MDK, 24/10/26
IF ( leaf_flux_mod == leaf_flux_stom_opt ) gc_som(:) = gc(:)
! The stomatal conductance itself: gc here, unless the leaf energy balance
! replaced it by gc_eq (then lt_finish kept it). MDK, 24/10/26
IF ( .NOT. l_leaf_temp ) gc_stomatal_pft(:,ft) = gc(:)

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
    el(l)       = e_of_g(gc(l), MAX(dqc(l), 0.0), pstar(l), tstar(l))
    psi_leaf(l) = psi_root_zone(l)
    leaf_k(l)   = kmax_pft(ft)   ! per unit leaf area, as the other paths
  END DO
END IF

!-----------------------------------------------------------------------------
! Cuticular loss (l_som_cuticular_loss). Water still leaks through the
! cuticle once the stomata have (nearly) shut, whatever the carbon: the
! canopy conductance is not allowed below gcuticular * LAI. This is applied
! after the optimisation, so it is not traded against carbon (A is
! unchanged), but it is charged to the plant water: the soil-supply cap
! still bounds it and
! psi_leaf / leaf_k are re-solved for the total flux.
!
! DESICA (here): applies at night too (closed points), and the flux is cut
! to the largest the xylem can carry with k > kcrit.
!
! The stomatal optimisation: see the next block.
!
! MDK, 22/10/26
!-----------------------------------------------------------------------------
IF ( l_som_cuticular_loss .AND. stomata_model == stomata_desica ) THEN
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
    kmax_cut(l)  = MAX(k_max(l) * lai(l), TINY(1.0_real_jlslsm))
    kcrit_cut(l) = kcrit_of_path(ft, kmax_pft(ft) * lai(l))
    IF ( kplant_model == kplant_model_ground ) THEN   ! MDK, 21/10/26
      kmax_cut(l)  = kmax_cut(l)  * kplant_scale(l)
      kcrit_cut(l) = kcrit_cut(l) * kplant_scale(l)
    END IF
    l_cut(l) = gl_cut(l) > 0.0 .AND. gc(l) < gl_cut_eff(l)
    ! Canopy transpiration (mol H2O m-2 s-1) at the floor, as in
    ! stom_opt_mod_ci.
    el_cut(l) = 0.0
    IF ( l_cut(l) ) el_cut(l) = e_of_g(gl_cut_eff(l), MAX(dqc(l), 0.0),        &
                                       pstar(l), tstar(l))
    DO i_cut = 1,n_cut
      e_cut(i_cut,m) = el_cut(l) * REAL(i_cut) / REAL(n_cut)
    END DO
  END DO
  IF ( ANY(l_cut(veg_index(1:veg_pts))) ) THEN
    CALL psi_leaf_jls( ft, n_cut, land_pts, veg_pts, veg_index,                &
                       veg_pts_index, e_cut, psi_src,                          &
                       kmax_cut, kcrit_cut, conductance_b, conductance_c,      &
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
          gc_som(l)   = gc(l)
          ! el is mol H2O m-2 s-1 and leaf_k per unit leaf area here.
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
! Cuticular loss under the stomatal optimisation (l_som_cuticular_loss).
! What it represents (Choat et al. 2018, Nature 558: 531-539, Fig. 2):
! stomata close in a drought before significant cavitation, but the plant
! goes on losing water through the cuticle, leaky stomata and bark,
! 100-1,000 times more slowly than with the stomata open. gcuticular is
! that residual (minimum) conductance, per unit leaf area. So psi_x goes on
! falling after closure; the fall is buffered by the release of stored
! water (the capacitance store below), while the conductance of the path
! falls (the xylem impairment memory and, in dry soil, the soil-to-root
! conductance), so the soil supplies less and less and the plant lives off
! its store. If the drought persists, psi_x reaches the threshold at which
! emboli propagate (about P50 in conifers, P88 in angiosperms) and the
! hydraulic system fails; here psi is held at the critical loss of
! conductance (kcrit_fractional_loss, e.g. 0.95, beyond P88), as there is
! no mortality.
!
! Where the optimiser has (nearly) shut the stomata, the canopy conductance
! is raised to gcuticular * LAI. There is no carbon trade at that point, so
! A is unchanged; the leak goes into the transpiration (and so ET and LE)
! through gc. Its water comes from the capacitance store and the soil (the
! soil-supply cap applies: no water is taken that is not there). Unlike
! DESICA:
! - only at lit points with apar >= apar_cut_min (about 5 W m-2 of
!   shortwave), as no stomatal conductance is computed at night;
! - the leak is not cut back when the plant nears kcrit, so in a drought
!   it keeps drying the plant;
! - the leak drains the capacitance store at the stem node (the carry loop
!   below), psi_leaf follows it, and the impairment update is made after
!   this block, so the leak drives the loss of conductance.
! psi is not allowed below the psi of the critical loss of conductance
! (kcrit_fractional_loss, e.g. 0.95) on the PFT curve (or the leaf
! and stem segments), so where the path cannot carry the leak (no steady
! state without storage) psi does not fall to the solver's floor and the
! memory does not go past that loss. The transpiration is uncoupled at
! tstar (dqc), as the canopy flux in sf_evap; for a conductance this small
! the canopy air resistance makes ~1 % difference.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
IF ( l_cap_store ) THEN
  CALL cut_state_alloc( land_pts )
  cut_q_store(:,ft) = 0.0
  ! The multilayer canopy has no lit set from its own path: the points
  ! with absorbed PAR. MDK, 23/10/26
  IF ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) THEN
    open_pts = 0
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( apar(l) > 0.0 ) THEN
        open_pts = open_pts + 1
        open_land_index(open_pts) = l
      END IF
    END DO
  END IF
  lk_cut(:)  = .FALSE.
  lit_cut(:) = .FALSE.
  DO m = 1,open_pts
    l = open_land_index(m)
    lit_cut(l) = apar(l) >= apar_cut_min
  END DO
  ! The leak is driven by the leaf-to-air deficit (as gmin * VPD): with
  ! leaf temperatures (two-leaf), qs at the sunlit and shaded leaf
  ! temperatures weighted by leaf area, less q1; otherwise the air deficit,
  ! qs(tair) - q1. Not dqc, at tstar, which in a drought is the hot tile
  ! skin, up to 10 K above the leaves (leak ~2x too large). sf_evap takes
  ! the flux at tstar, so gc is scaled by dq_cut / dqc to give it.
  !
  ! MDK, 22/10/26
  CALL qsat(qs_air_cut, tair, pstar, land_pts)
  DO m = 1,veg_pts
    l = veg_index(m)
    dq_cut(l) = MAX(qs_air_cut(l) - q1(l), 0.0)
    IF ( l_leaf_temp .AND. can_rad_mod == 7 ) THEN
      IF ( lai_sun_2l(l) + lai_shd_2l(l) > EPSILON(0.0) )                      &
        dq_cut(l) = ( lai_sun_2l(l) * MAX(qs_sun_lt(l) - q1(l), 0.0)           &
                      + lai_shd_2l(l) * MAX(qs_shd_lt(l) - q1(l), 0.0) )       &
                    / ( lai_sun_2l(l) + lai_shd_2l(l) )
    END IF
  END DO
  DO m = 1,veg_pts
    l = veg_index(m)
    ! mmol m-2 leaf s-1 -> m s-1, times LAI for the canopy.
    gl_cut(l) = gcuticular(ft) * 1.0e-3 * rmol * tstar(l) / pstar(l) * lai(l)
    ! No leak with capacitance alone (l_som_capacitance). MDK, 23/10/26
    IF ( .NOT. l_som_cuticular_loss ) gl_cut(l) = 0.0
    IF ( dqc(l) > 0.0 ) gl_cut(l) = gl_cut(l)                                  &
                                    * MIN(dq_cut(l) / dqc(l), dq_ratio_cut_max)
    share_sup(l) = 1.0
  END DO
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,             &
                           share_sup, dqc, tstar, pstar, gl_cut, gl_cut_eff )
  CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( .NOT. lit_cut(l) .OR. gl_cut(l) <= 0.0 .OR.                           &
         gc(l) >= gl_cut_eff(l) ) CYCLE
    gc(l)     = gl_cut_eff(l)
    gc_som(l) = gc(l)
    ! Canopy transpiration (mol H2O m-2 s-1), as el elsewhere.
    el(l)     = e_of_g(gc(l), MAX(dqc(l), 0.0), pstar(l), tstar(l))
    lk_cut(l) = .TRUE.
  END DO
  ! The canopy path of the optimisation (cf. kmax_bigleaf).
  DO m = 1,veg_pts
    l = veg_index(m)
    kmax_cut(l)     = k_max(l) * lai(l)
    kmax_ref_cut(l) = kmax_pft(ft) * lai(l)
    kcrit_cut(l)    = kcrit_of_path(ft, kmax_ref_cut(l))
    IF ( kplant_model == kplant_model_ground ) THEN
      kmax_cut(l)     = kmax_cut(l)     * kplant_scale(l)
      kmax_ref_cut(l) = kmax_ref_cut(l) * kplant_scale(l)
      kcrit_cut(l)    = kcrit_cut(l)    * kplant_scale(l)
    END IF
  END DO
  ! Floors on psi: the psi of the critical loss of conductance
  ! (kcrit_fractional_loss), where the memory stops, on the curve of the
  ! leaf and stem segments, or the PFT curve without segments.
  r_crit = 1.0 - kcrit_fractional_loss(ft)
  l_seg_cut = l_som_plant_segments .AND. .NOT. is_spinup .AND.                 &
              pft_xylem_impairment_model(ft) /= xylem_impairment_none
  IF ( l_seg_cut ) THEN
    psi_fl_leaf = -ABS(conductance_b_seg(ft,3))                                &
                  * (-LOG(r_crit))**(1.0 / conductance_c_seg(ft,3))
    psi_fl_stem = -ABS(conductance_b_seg(ft,2))                                &
                  * (-LOG(r_crit))**(1.0 / conductance_c_seg(ft,2))
  ELSE IF ( pft_conductance_model(ft) == SOX_conductance ) THEN
    psi_fl_leaf = -ABS(conductance_b_pft(ft))                                  &
                  * (1.0 / r_crit - 1.0)**(1.0 / conductance_c_pft(ft))
    psi_fl_stem = psi_fl_leaf
  ELSE
    psi_fl_leaf = -ABS(conductance_b_pft(ft))                                  &
                  * (-LOG(r_crit))**(1.0 / conductance_c_pft(ft))
    psi_fl_stem = psi_fl_leaf
  END IF
  ! Advance the capacitance store over the step. Its state is its water,
  !
  !   W = C psi_x      (mol H2O m-2 ground, <= 0; psi_x its psi, Pa),
  !
  ! carried between steps (and dumped). Four cases:
  !
  ! 1. Almost no leaves (LAI < lai_cut_min): the store is dropped (unset).
  !
  ! 2. The store sets E (tgt_on), or the leak binds: the store sits at the
  !    stem node, filled from the soil through the root side k_r and drawn
  !    by the step's transpiration E, held over the step:
  !
  !      psi_x = psi_e + (psi_x0 - psi_e) exp(-dt k_r / C),
  !      psi_e = psi_s - E / k_r      (the steady state for this E),
  !
  !    and the leaf is psi_leaf = psi_x - E / k_l. Where the leak binds,
  !    E is cut so that the leaf ends no lower than psi_fl_leaf, the
  !    critical-loss psi.
  !
  ! 3. The store is not yet set (first step, or a dump without it): it
  !    starts at the step's psi (in balance with it, not full).
  !
  ! 4. Everywhere else (the night, unlit points, stomata shut with no
  !    leak): the store relaxes towards the step's psi, psi_t, through the
  !    root side,
  !
  !      psi_x = psi_t + (psi_x0 - psi_t) exp(-dt k_r / C),
  !
  !    with k_r from the mean of the plant conductance at psi_x0 and psi_t
  !    on the optimiser's own path (cut_k_path: segments with the
  !    vulnerable root, or the PFT curve; impaired kmax; ksr in series;
  !    floored at the critical loss): a cheap stand-in for the exact path
  !    solve, which at every step made the run ~25x slower. It refills from
  !    the soil when drier than psi_t and gives up water to the
  !    transpiration when wetter (at most the transpiration: no flow back
  !    to the soil).
  !
  ! The store's flow out over the step,
  !
  !   q = C (psi_x0 - psi_x) / dt      (> 0: it gives; < 0: it refills),
  !
  ! goes to cut_q_store (kg m-2 s-1), which sf_evap takes off the soil
  ! extraction, so the water is conserved.
  !
  ! With the store coupled to the optimiser (som_store_coupling = 1), cases
  ! 2 and 4 are replaced, wherever the store is set, day and night, by its
  ! balance over the step at the final canopy E (store_cpl_q, exact):
  !
  !   q = C (psi_x0 - psi_x) / dt = E - E_r,   psi_x = psi_st(E_r),
  !
  ! the stem node psi_st where the soil, root and stem carry E_r, and the
  ! leaf psi_leaf where the leaf segment then carries E.
  !
  ! MDK, 24/10/26
  IF ( l_store_cpl ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      e_be_2l(l) = MAX(el(l), 0.0)
    END DO
    CALL store_cpl_q( .TRUE., e_be_2l )
  END IF
  DO m = 1,veg_pts
    l = veg_index(m)
    ! The store's capacitance (mol m-2 Pa-1), the step's psi (the
    ! optimiser's psi_leaf, Pa) and last step's stem psi (Pa).
    c_cut = cut_c_store(l,ft)
    psi_tgt   = psi_leaf(l)
    psi0_stem = cut_psi_stem(l,ft)
    cut_psi_stem(l,ft) = psi_stem(l)

    ! Case 1: almost no leaves, the store is dropped.
    IF ( lai(l) < lai_cut_min ) THEN
      cut_w_store(l,ft)  = cut_psi_unset
      cut_psi_stem(l,ft) = cut_psi_unset
      cut_q_store(l,ft)  = 0.0
      CYCLE
    END IF
    ! som_store_coupling = 1, the store set: its balance at the final E.
    ! (MDK, 24/10/26)
    IF ( l_store_cpl .AND. cut_w_store(l,ft) <= 0.0 ) THEN
      psi0_store = cut_w_store(l,ft) / c_cut
      CALL store_commit( l, c_cut, psi0_store, psi_st_be(l) )
      cut_psi_stem(l,ft) = psi_st_be(l)
      psi_leaf(l) = psi_l_be(l)
      psi_stem(l) = psi_st_be(l)
      CYCLE
    END IF
    ! Case 2: the store sets E, or the leak binds (and the store is set).
    IF ( tgt_on(l) .OR. ( lk_cut(l) .AND.                                      &
                          cut_w_store(l,ft) <= 0.0 ) ) THEN
      ! The store's psi at the start of the step (Pa).
      psi0_store = cut_w_store(l,ft) / c_cut
      ! Where the leak binds at points the optimiser shut, there is no
      ! secant conductance from store_target_e: the root and leaf sides
      ! come from the optimiser's own plant path (cut_k_path: the segments,
      ! with the vulnerable root, or the PFT curve; impaired kmax; ksr in
      ! series), the mean of its conductance at the store's psi and the
      ! source's, split by som_leaf_resist_frac.
      IF ( .NOT. tgt_on(l) ) THEN
        k_cut_eff = 0.5 * ( cut_k_path(l, kmax_cut(l), psi0_store)             &
                            + cut_k_path(l, kmax_cut(l), psi_src(l)) )
        k_cut_eff = MAX(k_cut_eff, TINY(1.0_real_jlslsm))
        tgt_kl(l) = k_cut_eff / f_leaf_store
        tgt_kr(l) = k_cut_eff / (1.0 - f_leaf_store)
      END IF
      ! The fraction of the store's departure from its steady state left
      ! after the step (-).
      r_st = EXP(-timestep * tgt_kr(l) / c_cut)
      ! l_som_leak_capacitance_only: where the leak binds, the store alone
      ! pays it (no soil supply), so
      !
      !   psi_x = psi_x0 - E dt / C,
      !
      ! and E is cut where the leaf would end below psi_fl_leaf.
      !
      ! MDK, 23/10/26
      IF ( l_som_leak_capacitance_only .AND. lk_cut(l) ) THEN
        IF ( el(l) > 0.0 ) THEN
          e_fl = MAX( (psi0_store - psi_fl_leaf)                               &
                      / (timestep / c_cut + 1.0 / tgt_kl(l)), 0.0 )
          IF ( e_fl < el(l) ) THEN
            gc(l)     = gc(l) * e_fl / el(l)
            gc_som(l) = gc(l)
            el(l)     = e_fl
          END IF
        END IF
        psi_st = psi0_store - MAX(el(l), 0.0) * timestep / c_cut
        CALL store_commit( l, c_cut, psi0_store, psi_st )
        cut_psi_stem(l,ft) = psi_st
        psi_leaf(l) = psi_st - MAX(el(l), 0.0) / tgt_kl(l)
        IF ( l_seg_cut ) psi_stem(l) = psi_st
        CYCLE
      END IF
      ! The leak cut to the largest E that leaves the leaf at psi_fl_leaf
      ! at the end of the step (the E' equation of store_target_e with the
      ! target psi_fl_leaf), with gc scaled to it.
      IF ( lk_cut(l) .AND. el(l) > 0.0 ) THEN
        e_fl = e_leaf_end( psi_src(l), psi0_store, r_st, tgt_kr(l),            &
                           tgt_kl(l), psi_fl_leaf )
        e_fl = MAX(e_fl, 0.0)
        IF ( e_fl < el(l) ) THEN
          gc(l)     = gc(l) * e_fl / el(l)
          gc_som(l) = gc(l)
          el(l)     = e_fl
        END IF
      END IF
      ! The stem node over the step, its flow out and its water, and the
      ! leaf at psi_x - E / k_l.
      psi_st = psi_src(l) - MAX(el(l), 0.0) / tgt_kr(l)
      psi_st = psi_st + (psi0_store - psi_st) * r_st
      CALL store_commit( l, c_cut, psi0_store, psi_st )
      cut_psi_stem(l,ft) = psi_st
      psi_leaf(l) = psi_st - MAX(el(l), 0.0) / tgt_kl(l)
      IF ( l_seg_cut ) psi_stem(l) = psi_st
      CYCLE
    END IF
    ! Case 3: the store is not yet set: it starts at the step's psi. (Leak
    ! points with the store set went to case 2, so only the unset test is
    ! needed here.)
    !
    ! MDK, 23/10/26
    IF ( cut_w_store(l,ft) > 0.0 ) THEN
      cut_w_store(l,ft) = c_cut * psi_leaf(l)
      CYCLE
    END IF
    ! Case 4: the store relaxes towards the step's psi through the root
    ! side (unchanged where it is within 1 Pa of it).
    psi0_store = cut_w_store(l,ft) / c_cut
    psi_st     = psi_tgt
    IF ( ABS(psi0_store - psi_tgt) > 1.0 ) THEN
      ! The plant conductance between psi_x0 and psi_t on the optimiser's
      ! own path (cut_k_path: the segments, with the vulnerable root, or
      ! the PFT curve; ksr in series, so a dry rhizosphere slows the
      ! refill too) (mol m-2 s-1 Pa-1).
      k_cut_eff = 0.5 * ( cut_k_path(l, kmax_cut(l), psi0_store)               &
                          + cut_k_path(l, kmax_cut(l), psi_tgt) )
      ! With the store at the stem node, only the root side fills it.
      !
      ! MDK, 23/10/26
      k_cut_eff = k_cut_eff / (1.0 - f_leaf_store)
      ! The relaxation over the step, and the store's flow out (mol m-2
      ! s-1). The store gives at most the transpiration (no flow back to
      ! the soil): where it would give more, the relaxation is slowed.
      r_st = EXP(-timestep * k_cut_eff / c_cut)
      q_st = c_cut * (psi0_store - psi_tgt) * (1.0 - r_st) / timestep
      IF ( q_st > MAX(el(l), 0.0) ) THEN
        q_st = MAX(el(l), 0.0)
        r_st = 1.0 - q_st * timestep / (c_cut * (psi0_store - psi_tgt))
      END IF
      psi_st = psi_tgt + (psi0_store - psi_tgt) * r_st
      IF ( l_seg_cut .AND. psi0_stem <= 0.0 .AND. psi_stem(l) < 0.0 )          &
        cut_psi_stem(l,ft) = psi_stem(l) + (psi0_stem - psi_stem(l)) * r_st
    END IF
    ! The store's flow out (kg m-2 s-1) and its water at the end of the step.
    CALL store_commit( l, c_cut, psi0_store, psi_st )
  END DO

  ! At the other lit points, those of the impairment update, psi_leaf (and
  ! psi_stem with segments) is the store's psi, so that the plant's psi,
  ! and the embolism it drives, does not jump between the lagged store and
  ! the optimiser's steady state. Points where the store sets E (tgt_on)
  ! or the leak binds (lk_cut) keep the leaf psi of the carry loop,
  !
  !   psi_leaf = psi_x - E / k_l,
  !
  ! below the stem node by the leaf-side drop. At night psi_leaf is the
  ! optimiser's.
  !
  ! MDK, 23/10/26
  ! (Not with som_store_coupling = 1, where the balance set them. MDK,
  ! 24/10/26)
  DO m = 1,open_pts
    l = open_land_index(m)
    IF ( l_store_cpl ) EXIT
    IF ( cut_w_store(l,ft) > 0.0 .OR. tgt_on(l) .OR. lk_cut(l) ) CYCLE
    psi_leaf(l) = cut_w_store(l,ft) / cut_c_store(l,ft)
    IF ( l_seg_cut .AND. cut_psi_stem(l,ft) <= 0.0 )                           &
      psi_stem(l) = cut_psi_stem(l,ft)
  END DO

  ! Xylem impairment update, once per timestep, in place of those in the
  ! big-leaf and two-leaf paths: the same lit points, with the water
  ! potentials of the total (stomatal plus cuticular) flux.
  !
  ! MDK, 23/10/26
  CALL update_xylem_impairment(                                                &
          ! IN
            land_pts, open_pts, open_land_index, ft, psi_leaf, psi_src,        &
            psi_stem,                                                          &
          ! IN OUT
            k_max, psi_leaf_extreme, psi_root_extreme,                         &
          ! IN (recovery)
            lai, canht, anetc )
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
  CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,          &
                          lai, ht, psi_root_zone, ksr_ds, el, .FALSE.,         &
                          psi_try, k_try, el_hyd )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( el_hyd(l) < el(l) ) THEN
      gc(l) = gc(l) * el_hyd(l) / el(l)
      el(l) = el_hyd(l)
    END IF
  END DO
  CALL desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,          &
                          lai, ht, psi_root_zone, ksr_ds, el, .FALSE.,         &
                          psi_leaf, leaf_k, el_hyd )
  ! desica_hydraulics returns Kplant (soil-to-leaf) per m2 ground; per unit
  ! leaf area, as the other paths (kplant_pft).
  leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))
  CALL desica_store_inputs( ft, land_pts, veg_pts, veg_index, lai, ht,         &
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
! Leaf 13C discrimination diagnostic (discrimination_model). With the
! two-leaf canopy it was weighted over the sun and shade classes above;
! otherwise it is taken from the canopy ci. Zero with no assimilation.
!-----------------------------------------------------------------------------
IF ( discrimination_model /= discrimination_off ) THEN
  CALL d13c_disc_alloc( land_pts )
  d13c_disc_pft(:,ft) = 0.0
  l_d13c_classes = ( can_rad_mod == 7 )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( l_d13c_classes ) THEN
      d13c_disc_pft(l,ft) = d13c_canopy(l)
    ELSE IF ( anetc(l) + rdc(l) > 0.0 ) THEN
      d13c_disc_pft(l,ft) = d13c_disc_leaf( c3(ft) == 1, ci(l), ca(l),         &
                                            ccp(l) )
    END IF
  END DO
END IF

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
  IF (l_red) THEN
    lai_bal(l) = veg_state%lai_bal(l,ft)

  ELSE
    lai_bal(l) = lai_bal_from_ht(ft, ht(l))   ! MDK, 20/10/26

  END IF

  !---------------------------------------------------------------------------
  ! Calculate the total nitrogen content of the leaf, root and stem
  !---------------------------------------------------------------------------
  IF ( ft > nnpft ) THEN
    !Crop PFTs
    stemc   = stemc_from_prognostics(ft - nnpft, canht(l) )
    root(l) = rootc_cpft(l,ft - nnpft)

    IF ( l_trait_phys ) THEN

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft,                     &
                                  can_averaging_fac(l) )
      n_root(l) = nr_nl(ft) * nmass(ft) * (1.0 / cfrac_l(ft - nnpft))          &
                  * root(l)
      n_stem(l) = ns_nl(ft) * nmass(ft) * (1.0 / cfrac_l(ft - nnpft)) * stemc
        ! nmass(ft)/cfrac_l(ft-nnpft) is equivalent to nl0(ft)

      lma_tmp   = lma_from_prognostics(ft - nnpft, dvi_cpft(l,ft - nnpft))
      nlmean(l) = nmass(ft) * lma_tmp * can_averaging_fac(l)

    ELSE

      n_leaf(l) = nleaf_from_lai( l, ft, lai(l), dvi_cpft,                     &
                                  can_averaging_fac(l) )
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
      IF (l_red) THEN
        n_stem(l) = veg_state%woodC(l,ft) / a_ws(ft) * nsw(ft)
      ELSE
        ! MDK, 20/10/26
        n_stem(l) = live_stemwood(ft, ht(l), lai_bal(l)) * nsw(ft)
      END IF

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

! Turgor-limited sink of the TRIFFID store: water status (and the cold
! limit on the daily mean air temperature) over the store period (returns at
! once when neither is set).
! MDK, 18/10/26
IF (l_triffid_nsc)                                                             &
  CALL triffid_nsc_turgor( land_pts, veg_pts, veg_index, ft, timestep,         &
                           psi_root_zone, tair )

! Drought leaf shedding (l_phen_psi): daily water status of the PFT.
! MDK, 18/10/26
IF (l_phen_psi)                                                                &
  CALL phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,            &
                        psi_root_zone )

! TRIFFID store without TRIFFID (prescribed LAI): daily store step.
! MDK, 17/10/26
IF (l_triffid_nsc .AND. .NOT. l_triffid)                                       &
  CALL triffid_nsc_daily( land_pts, veg_pts, veg_index, ft, timestep, ht,      &
                          npp )

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

CONTAINS

!-----------------------------------------------------------------------------
! Plant capacitance: the canopy plant conductance (mol H2O m-2 ground s-1
! Pa-1) at psi (Pa) on the optimiser's own path,
!
!   k = k_path_zero_flow(psi)
!
! (every element at psi): with plant segments the root, stem and leaf
! segments in series, so the calibrated (vulnerable) root segment limits the
! store's exchange with the soil as it does the optimiser's supply, else
! the PFT curve; maximum kmx (the impaired canopy kmax); the soil-to-root
! link (ksr) in series. psi is floored at the critical loss (psi_fl_leaf)
! and capped at 0. This replaces the PFT curve with ksr in series
! (cut_k_curve), which left out the vulnerable root, so the soil went on
! feeding the leak and the refill through a near-intact root (60-80 % of
! the leak in drought summers at FR-Pue).
!
! With segments and xylem impairment the path is the optimiser's
! (stom_opt_mod): the intact kmax_ref, with the stem and leaf segments
! capped at
!
!   kcap_frac = kmx / kmax_ref,
!
! so the root is not impaired. Without segments the impaired kmx scales
! the PFT curve.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION cut_k_path( l_c, kmx, psi )

INTEGER, INTENT(IN) :: l_c
                            ! Land point.
REAL(KIND=real_jlslsm), INTENT(IN) :: kmx
                            ! Impaired canopy kmax (mol H2O m-2 ground s-1
                            ! Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
                            ! Water potential (Pa).
                            !
                            ! MDK, 23/10/26

IF ( l_seg_cut ) THEN
  cut_k_path = k_path_zero_flow( ft, l_c, kmax_ref_cut(l_c),                   &
                 MIN(MAX(psi, psi_fl_leaf), 0.0),                              &
                 kcap_frac = kmx                                               &
                             / MAX(kmax_ref_cut(l_c), TINY(1.0_real_jlslsm)) )
ELSE
  cut_k_path = k_path_zero_flow( ft, l_c, kmx, MIN(MAX(psi, psi_fl_leaf), 0.0) )
END IF

END FUNCTION cut_k_path

!-----------------------------------------------------------------------------
! Plant capacitance: the transpiration, held over the step, that brings the
! leaf to psi_end at the end of the step, with the store at the stem node
! (store_target_e, where psi_end is the optimiser's psi_leaf, and the leak's
! cut in the carry loop, where it is the critical-loss psi):
!
!   E = [psi_s (1 - e) + psi_x0 e - psi_end] / [(1 - e) / k_r + 1 / k_l],
!
! e the fraction of the store's departure from its steady state left after
! the step, exp(-dt k_r / C). Not floored at 0 here.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION e_leaf_end( psi_s, psi_x0, e_r, k_r, k_l,      &
                                            psi_end )

REAL(KIND=real_jlslsm), INTENT(IN) :: psi_s
                            ! Source (root-zone) psi (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_x0
                            ! The store's psi at the start of the step (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: e_r
                            ! exp(-dt k_r / C) (-).
REAL(KIND=real_jlslsm), INTENT(IN) :: k_r
                            ! Root-side conductance (mol H2O m-2 ground s-1
                            ! Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: k_l
                            ! Leaf-side conductance (mol H2O m-2 ground s-1
                            ! Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_end
                            ! The leaf psi at the end of the step (Pa).
                            !
                            ! MDK, 23/10/26

e_leaf_end = ( psi_s * (1.0 - e_r) + psi_x0 * e_r - psi_end )                  &
             / ( (1.0 - e_r) / k_r + 1.0 / k_l )

END FUNCTION e_leaf_end

!-----------------------------------------------------------------------------
! Plant capacitance: the store at land point l_c ends the step at psi_x1
! from psi_x0, so its flow out over the step and its water are
!
!   cut_q_store = C (psi_x0 - psi_x1) / dt M_w      (kg H2O m-2 s-1),
!   cut_w_store = C psi_x1                           (mol H2O m-2).
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_commit( l_c, c_s, psi_x0, psi_x1 )

INTEGER, INTENT(IN) :: l_c
                            ! Land point.
REAL(KIND=real_jlslsm), INTENT(IN) :: c_s
                            ! The store's capacitance C (mol H2O m-2 ground
                            ! Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_x0
                            ! The store's psi at the start of the step (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_x1
                            ! The store's psi at the end of the step (Pa).
                            !
                            ! MDK, 23/10/26

cut_q_store(l_c,ft) = c_s * (psi_x0 - psi_x1) / timestep * cut_m_h2o
cut_w_store(l_c,ft) = c_s * psi_x1

END SUBROUTINE store_commit

!-----------------------------------------------------------------------------
! som_store_coupling = 1: the soil supply cap of the two-leaf optimiser.
! The soil supplies E - q_f, so the cap on E is
!
!   e_supply + q_f      (kg m-2 s-1, no less than 0),
!
! and no cap (e_supply < 0) stays no cap.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_supply()

IMPLICIT NONE

INTEGER :: m, l

DO m = 1,veg_pts
  l = veg_index(m)
  e_sup_use(l) = e_supply(l)
  IF ( e_supply(l) >= 0.0 )                                                    &
    e_sup_use(l) = MAX(e_supply(l) + q_f(l) * cut_m_h2o, 0.0)
END DO

END SUBROUTINE store_cpl_supply

!-----------------------------------------------------------------------------
! som_store_coupling = 1: the store's balance over the step at the canopy
! transpiration e_c, at the points where the store is set. Backward Euler
! at the stem node:
!
!   q   = C' (psi_x0 - psi_st),   C' = C / dt,
!   E_r = E - q >= 0              (soil, root and stem; no flow back),
!
! with psi_st the stem-node psi where the soil link, root and stem carry
! E_r: the segment chain from psi_src on the canopy's path (the two classes'
! kmax, intact kmax and kcrit summed; ksr in full). The residual
!
!   R(q) = q - C' (psi_x0 - psi_st(E - q))
!
! rises with q (a larger q leaves less E_r, so a wetter psi_st), and as
! psi_st <= psi_src the root lies between
!
!   q_lo = C' (psi_x0 - psi_src)  (R <= 0)   and   q_hi = E  (R > 0).
!
! It is found by regula falsi (Illinois), with a bisection step whenever
! the bracket has shrunk by less than half.
!
! Where q_lo >= E, the store is wetter than the soil by more than the
! step's E: it pays all of E (E_r = 0), psi_st = psi_x0 - E / C', and the
! leaf is psi_st less the leaf segment's drop at E (taken from psi_src).
!
! Sets q_be (= C' (psi_x0 - psi_st), as store_commit takes it), psi_st_be
! and psi_l_be, the leaf psi where the leaf segment carries E from psi_st.
! Elsewhere q_be = 0, psi_st_be = psi_src and psi_l_be = psi_leaf.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_q( l_ex, e_c )

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_ex
                            ! .TRUE.: exact segment solves (the commit);
                            ! .FALSE.: the Newton steps (between passes).
REAL(KIND=real_jlslsm), INTENT(IN) :: e_c(land_pts)
                            ! Canopy transpiration (mol H2O m-2 s-1).

INTEGER, PARAMETER :: n_it_be = 40
                            ! Most regula falsi steps.
REAL(KIND=real_jlslsm), PARAMETER :: rtol_be = 1.0e-6
                            ! Tolerance on the bracket of q, relative to
                            ! its first width.

INTEGER :: m, l, j, it, n_be
INTEGER :: idx_be(land_pts)
                            ! Positions in veg_index of the points solved.
LOGICAL :: l_imp_be, l_all_done
LOGICAL :: l_only(land_pts), l_bis(land_pts), l_done(land_pts)
                            ! The store pays all of E; next step bisects;
                            ! converged.
INTEGER :: side_be(land_pts)
                            ! Illinois: the end replaced last (-1 lo, 1 hi).
REAL(KIND=real_jlslsm) ::                                                      &
  c_p(land_pts), psi_x0(land_pts),                                             &
      ! C / dt (mol m-2 s-1 Pa-1) and the store's psi at the start (Pa).
  q_lo(land_pts), q_hi(land_pts), f_lo(land_pts), f_hi(land_pts),              &
      ! The bracket on q and R at its ends (mol m-2 s-1).
  q_try(land_pts), f_try(land_pts), w_tol(land_pts), w_old(land_pts),          &
      ! Trial q and R there; the converged width and the last width.
  kref_c(land_pts), kmx_c(land_pts), kcr_c(land_pts), kcap_c(land_pts)
      ! The canopy's intact and impaired kmax, kcrit and their ratio.
REAL(KIND=real_jlslsm) :: psi1(1,land_pts), st1(1,land_pts)
                            ! One sample per point: the leaf and stem psi of
                            ! the chain.
                            !
                            ! MDK, 24/10/26

q_be(:)      = 0.0
psi_st_be(:) = psi_src(:)
psi_l_be(:)  = psi_leaf(:)
! (q_try is passed whole to set_store_q_path, so the points not solved
! hold 0, not junk. MDK, 24/10/26)
q_try(:)     = 0.0

! The canopy's path (the two classes together) and the soil link in full.
l_imp_be = .NOT. is_spinup .AND.                                               &
           pft_xylem_impairment_model(ft) /= xylem_impairment_none
DO m = 1,veg_pts
  l = veg_index(m)
  kref_c(l) = kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l)
  kmx_c(l)  = kmax_sun_2l(l) + kmax_shd_2l(l)
  kcr_c(l)  = kcrit_sun_2l(l) + kcrit_shd_2l(l)
  kcap_c(l) = kmx_c(l) / MAX(kref_c(l), TINY(1.0_real_jlslsm))
  share_sup(l) = 1.0
END DO
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

! The points with the store set, their bracket, and those where the store
! pays all of E.
n_be = 0
DO m = 1,veg_pts
  l = veg_index(m)
  l_only(l) = .FALSE.
  IF ( cut_w_store(l,ft) > 0.0 .OR. lai(l) < lai_cut_min ) CYCLE
  c_p(l)    = cut_c_store(l,ft) / timestep
  psi_x0(l) = cut_w_store(l,ft) / cut_c_store(l,ft)
  q_lo(l)   = c_p(l) * (psi_x0(l) - psi_src(l))
  q_hi(l)   = e_c(l)
  n_be = n_be + 1
  idx_be(n_be) = m
  IF ( q_lo(l) >= q_hi(l) ) THEN
    l_only(l)    = .TRUE.
    q_be(l)      = e_c(l)
    psi_st_be(l) = psi_x0(l) - e_c(l) / c_p(l)
  END IF
  l_bis(l)   = .FALSE.
  l_done(l)  = l_only(l)
  side_be(l) = 0
  w_tol(l)   = rtol_be * MAX(q_hi(l) - q_lo(l), TINY(1.0_real_jlslsm))
  w_old(l)   = q_hi(l) - q_lo(l)
END DO
IF ( n_be == 0 ) RETURN

! R at both ends of the bracket.
DO j = 1,n_be
  q_try(veg_index(idx_be(j))) = q_lo(veg_index(idx_be(j)))
END DO
CALL store_cpl_eval( n_be, idx_be, l_imp_be, l_ex, e_c, q_try, c_p, psi_x0,    &
                     kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_be
  l = veg_index(idx_be(j))
  f_lo(l) = f_try(l)
  q_try(l) = q_hi(l)
END DO
CALL store_cpl_eval( n_be, idx_be, l_imp_be, l_ex, e_c, q_try, c_p, psi_x0,    &
                     kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_be
  l = veg_index(idx_be(j))
  f_hi(l) = f_try(l)
END DO

! Regula falsi (Illinois) on q, with bisection where the bracket shrinks
! slowly.
DO it = 1,n_it_be
  l_all_done = .TRUE.
  DO j = 1,n_be
    l = veg_index(idx_be(j))
    IF ( l_done(l) ) THEN
      q_try(l) = 0.5 * (q_lo(l) + q_hi(l))
      CYCLE
    END IF
    l_all_done = .FALSE.
    IF ( l_bis(l) .OR. f_hi(l) - f_lo(l) <= 0.0 ) THEN
      q_try(l) = 0.5 * (q_lo(l) + q_hi(l))
    ELSE
      q_try(l) = (q_lo(l) * f_hi(l) - q_hi(l) * f_lo(l)) / (f_hi(l) - f_lo(l))
      q_try(l) = MIN(MAX(q_try(l), q_lo(l)), q_hi(l))
    END IF
  END DO
  IF ( l_all_done ) EXIT
  CALL store_cpl_eval( n_be, idx_be, l_imp_be, l_ex, e_c, q_try, c_p, psi_x0,  &
                       kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
  DO j = 1,n_be
    l = veg_index(idx_be(j))
    IF ( l_done(l) ) CYCLE
    IF ( f_try(l) < 0.0 ) THEN
      q_lo(l) = q_try(l)
      f_lo(l) = f_try(l)
      IF ( side_be(l) == -1 ) f_hi(l) = 0.5 * f_hi(l)
      side_be(l) = -1
    ELSE
      q_hi(l) = q_try(l)
      f_hi(l) = f_try(l)
      IF ( side_be(l) == 1 ) f_lo(l) = 0.5 * f_lo(l)
      side_be(l) = 1
    END IF
    l_bis(l)  = q_hi(l) - q_lo(l) > 0.5 * w_old(l)
    w_old(l)  = q_hi(l) - q_lo(l)
    l_done(l) = w_old(l) <= w_tol(l)
  END DO
END DO

! The state at the solution (the middle of the bracket; q = E where the
! store pays all of E).
DO j = 1,n_be
  l = veg_index(idx_be(j))
  q_try(l) = 0.5 * (q_lo(l) + q_hi(l))
  IF ( l_only(l) ) q_try(l) = e_c(l)
END DO
CALL store_cpl_eval( n_be, idx_be, l_imp_be, l_ex, e_c, q_try, c_p, psi_x0,    &
                     kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_be
  l = veg_index(idx_be(j))
  IF ( l_only(l) ) THEN
    ! The stem node is the store's; the leaf drop is the leaf segment's at
    ! E (from psi_src, where the chain put the stem node).
    psi_l_be(l) = psi_st_be(l) - (psi_src(l) - psi1(1,j))
  ELSE
    psi_st_be(l) = st1(1,j)
    psi_l_be(l)  = psi1(1,j)
    q_be(l)      = c_p(l) * (psi_x0(l) - psi_st_be(l))
  END IF
  ! The commit: no psi below the critical-loss psi of the stem and leaf
  ! (where E cannot be carried the chain returns its floor).
  IF ( l_ex ) THEN
    psi_st_be(l) = MAX(psi_st_be(l), psi_fl_stem)
    psi_l_be(l)  = MAX(psi_l_be(l), psi_fl_leaf)
    q_be(l)      = c_p(l) * (psi_x0(l) - psi_st_be(l))
  END IF
END DO

END SUBROUTINE store_cpl_q

!-----------------------------------------------------------------------------
! som_store_coupling = 1, store_cpl_q: R(q) = q - C' (psi_x0 - psi_st) at
! the n_be points solved, with the store's flow held at q_try on the
! canopy's path (psi1, st1: the chain's leaf and stem psi, Pa).
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_eval( n_be, idx_be, l_imp_be, l_ex, e_c, q_try, c_p,      &
                           psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1,    &
                           f_try )

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_be
                            ! Number of points solved.
INTEGER, INTENT(IN) :: idx_be(land_pts)
                            ! Their positions in veg_index.
LOGICAL, INTENT(IN) :: l_imp_be, l_ex
                            ! Xylem impairment active; exact segment solves.
REAL(KIND=real_jlslsm), INTENT(IN) :: e_c(land_pts), q_try(land_pts),          &
                                      c_p(land_pts), psi_x0(land_pts)
                            ! Canopy E and trial q (mol m-2 s-1), C / dt
                            ! (mol m-2 s-1 Pa-1), the store's psi (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: kref_c(land_pts), kmx_c(land_pts),       &
                                      kcr_c(land_pts), kcap_c(land_pts)
                            ! The canopy's intact and impaired kmax, kcrit
                            ! and their ratio.
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi1(1,land_pts), st1(1,land_pts)
                            ! Leaf and stem psi at each point solved (Pa).
REAL(KIND=real_jlslsm), INTENT(OUT) :: f_try(land_pts)
                            ! R at q_try (mol m-2 s-1).
                            !
                            ! MDK, 24/10/26

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: e1(1,land_pts), k1(1,land_pts)

DO j = 1,n_be
  e1(1,j) = e_c(veg_index(idx_be(j)))
END DO
CALL set_store_q_path( land_pts, veg_pts, veg_index, share_sup, q_try )
IF ( l_imp_be ) THEN
  CALL psi_leaf_segments_jls( ft, 1, land_pts, n_be, veg_index, idx_be,        &
                              e1(:,1:n_be), psi_src, kref_c, kcr_c,            &
                              psi1(:,1:n_be), k1(:,1:n_be),                    &
                              kcap_frac = kcap_c,                              &
                              stem_psi = st1(:,1:n_be), l_exact = l_ex )
ELSE
  CALL psi_leaf_segments_jls( ft, 1, land_pts, n_be, veg_index, idx_be,        &
                              e1(:,1:n_be), psi_src, kmx_c, kcr_c,             &
                              psi1(:,1:n_be), k1(:,1:n_be),                    &
                              stem_psi = st1(:,1:n_be), l_exact = l_ex )
END IF
CALL clear_store_q_path()
DO j = 1,n_be
  l = veg_index(idx_be(j))
  f_try(l) = q_try(l) - c_p(l) * (psi_x0(l) - st1(1,j))
END DO

END SUBROUTINE store_cpl_eval


!-----------------------------------------------------------------------------
! Multilayer canopy, one leaf class (sunlit or shaded) of the current layer:
! the supply cap and soil-to-root share for the class's deficit dq_x and
! temperature t_x, then the optimisation with the class's photosynthesis
! inputs (_x: at the leaf temperature with l_leaf_temp, else the canopy
! values at tstar). This was four copies of the same call. (Without
! l_leaf_temp both classes see dqc and tstar, so the cap is the same for
! both, as when it was set once for the pair.)
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
SUBROUTINE ml_class( vcmax_x, kc_x, ko_x, ccp_x, km_x, dq_x, je_x,             &
                     t_x, je_ratio_x, fapar_x, rd_c, ci_c, anetl_c, el_c,      &
                     flux_o3_c, fo3_c, gl_c, psi_leaf_c, cg_c, hc_c,           &
                     leaf_k_c, psi_stem_c )
IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  vcmax_x(land_pts), kc_x(land_pts), ko_x(land_pts), ccp_x(land_pts),          &
  km_x(land_pts), dq_x(land_pts), je_x(land_pts),                              &
  t_x(land_pts), je_ratio_x(land_pts), fapar_x(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: rd_c(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ci_c(land_pts), anetl_c(land_pts), el_c(land_pts), flux_o3_c(land_pts),      &
  fo3_c(land_pts), gl_c(land_pts), psi_leaf_c(land_pts), cg_c(land_pts),       &
  hc_c(land_pts), leaf_k_c(land_pts)
                            ! OUT from stom_opt_mod; IN OUT here so that
                            ! points it does not set keep their values.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: psi_stem_c(land_pts)
                            ! The stem psi of this leaf's path (Pa; with
                            ! plant segments).
                            !
                            ! MDK, 23/10/26

INTEGER :: i, l
REAL(KIND=real_jlslsm) :: e_tgt_x(land_pts)
                            ! Capacitance pass: this leaf's target E.
                            !
                            ! MDK, 23/10/26

e_tgt_x(:) = 0.0
CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply, share_sup,    &
                         dq_x, t_x, pstar, gl_max_lf, gl_max_eff )
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
CALL stom_opt_mod (                                                            &
    ! IN
    land_pts, ft, open_pts, open_index,                                        &
    pft_photo_model, veg_pts, veg_index,                                       &
    ca, psi_src, acr, apar, oa, vcmax_x, kc_x, ko_x, ccp_x, pstar,             &
    km_x, dq_x, je_x, t_x, je_ratio_x, fapar_x,                                &
    kmax_imp_lyr, kcrit_leaf_lyr, gl_max_eff, ipar,                            &
    l_multilayer, zero_c, zero_c, zero_c,                                      &
    kmax_leaf_lyr, b_pft_pts, c_pft_pts,                                       &
    psi_leaf_extreme, psi_root_extreme, l_ml_ximp,                             &
    ! IN OUT
    rd_c,                                                                      &
    ! OUT
    ci_c, anetl_c, el_c, flux_o3_c, fo3_c,                                     &
    gl_c, psi_leaf_c, cg_c, hc_c, leaf_k_c,                                    &
    psi_stem = psi_stem_c                                                      &
)

! Plant capacitance (the multilayer store's pass): this leaf at
!
!   E_leaf' = E_leaf E' / E*,
!
! its optimal E scaled by the canopy's E' / E* (r_tgt), with ci, A and gl
! re-solved there (stom_opt_at_e).
!
! MDK, 23/10/26
IF ( ml_tank_pass ) THEN
  ! The open points where the store applies, and their target E.
  n_tgt = 0
  DO i = 1,open_pts
    l = veg_index(open_index(i))
    IF ( tgt_on(l) ) THEN
      n_tgt = n_tgt + 1
      idx_tgt(n_tgt) = open_index(i)
      e_tgt_x(l) = el_c(l) * r_tgt(l)
    END IF
  END DO
  IF ( n_tgt > 0 )                                                             &
    CALL stom_opt_at_e( land_pts, ft, n_tgt, idx_tgt, pft_photo_model,         &
                        veg_index, rd_c, ca, psi_src, acr, apar, oa, vcmax_x,  &
                        kc_x, ko_x, ccp_x, pstar, km_x, dq_x, je_x, t_x,       &
                        je_ratio_x, fapar_x, ipar, kmax_imp_lyr,               &
                        kcrit_leaf_lyr, gl_max_eff, l_multilayer, zero_c,      &
                        zero_c, zero_c, kmax_leaf_lyr, b_pft_pts, c_pft_pts,   &
                        psi_leaf_extreme, psi_root_extreme, l_ml_ximp,         &
                        e_tgt_x,                                               &
                        ci_c, anetl_c, gl_c, el_c )
END IF

END SUBROUTINE ml_class

!-----------------------------------------------------------------------------
! Two-leaf canopy (can_rad_mod = 7), one leaf class (sunlit or shaded) for
! one pass: the open points, this class's share of the supply cap and of the
! soil-to-root conductance, the coupled-transpiration inputs (with the other
! class's gl, gl_o), the optimisation (at the class's leaf temperature with
! l_leaf_temp, else at tstar), and min_gl_pft shared by leaf area at the
! points that end up closed. The arguments are the class's arrays (_o: the
! other class's). This was two copies of the same code, one per class.
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
SUBROUTINE two_leaf_class( apar_c, lai_c, kmax_ref_c, qs_lt_c, dq_lt_c,        &
                           gbw_c, gbw_o, gl_o, t_lt_c, gl_max_c, acr_c,        &
                           vcmax_c, kc_lt_c, ko_lt_c, ccp_lt_c, km_lt_c,       &
                           je_c, kmax_c, kcrit_c, rd_c, ci_c, anetl_c, el_c,   &
                           flux_o3_c, fo3_c, gl_c, psi_leaf_c, cg_c, hc_c,     &
                           leaf_k_c, psi_stem_c, lam_c )
IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  apar_c(land_pts), lai_c(land_pts), kmax_ref_c(land_pts),                     &
  qs_lt_c(land_pts), dq_lt_c(land_pts), gbw_c(land_pts), gbw_o(land_pts),      &
  gl_o(land_pts), t_lt_c(land_pts), gl_max_c(land_pts), acr_c(land_pts),       &
  vcmax_c(land_pts), kc_lt_c(land_pts), ko_lt_c(land_pts),                     &
  ccp_lt_c(land_pts), km_lt_c(land_pts), je_c(land_pts), kmax_c(land_pts),     &
  kcrit_c(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: rd_c(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ci_c(land_pts), anetl_c(land_pts), el_c(land_pts), flux_o3_c(land_pts),      &
  fo3_c(land_pts), gl_c(land_pts), psi_leaf_c(land_pts), cg_c(land_pts),       &
  hc_c(land_pts), leaf_k_c(land_pts), psi_stem_c(land_pts)
                            ! OUT from stom_opt_mod; IN OUT here so that
                            ! points it does not set keep their values, as
                            ! with the arrays passed directly before.
REAL(KIND=real_jlslsm), INTENT(OUT) :: lam_c(land_pts)
                            ! dA/dE at the optimum (lambda_eff of
                            ! stom_opt_mod; 0 where not set).
                            !
                            ! MDK, 23/10/26

INTEGER :: i, l, m

clos_pts = 0
open_pts = 0
DO i = 1,veg_pts
  l = veg_index(i)
  IF ( apar_c(l) > 0.0 .AND. lai_c(l) > EPSILON(0.0) ) THEN
    open_pts = open_pts + 1
    open_index(open_pts) = i
  ELSE
    clos_pts = clos_pts + 1
    clos_index(clos_pts) = i
  END IF
END DO

! Soil supply cap, and this class's share of the soil-to-root conductance
! (set_ksr_path): its share of the undamaged plant conductance (kmax_ref),
! so that xylem impairment does not touch the soil term. lai_cost
! (stom_opt_mod): this class's leaf area, for the CMax / CGain cost per
! unit leaf area (the kmax share is N-weighted).
DO m = 1,veg_pts
  l = veg_index(m)
  share_sup(l) = kmax_ref_c(l) /                                               &
                 MAX(kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l), TINY(1.0))
END DO

CALL coupled_inputs( land_pts, veg_pts, veg_index, l_leaf_temp, dq_min,        &
                     dqc, dq, qs_lt_c, q1, dq_lt_c, gbw_c, gbw_o, ra_ca_lt,    &
                     gl_o, g_soil_c, qs, dq_opt, r_bl_c, r_ca_c, g_oth_c )

! The store coupled to the optimiser (som_store_coupling = 1): this class's
! share of the store's flow out, held fixed for its search (cleared after).
!
! MDK, 24/10/26
IF ( l_store_cpl ) CALL set_store_q_path( land_pts, veg_pts, veg_index,        &
                                          share_sup, q_f )

! At the class's leaf temperature (l_leaf_temp), else tstar: the class
! arrays hold the tile's values then (see the two-leaf set-up).
!
! MDK, 24/10/26
CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_use,              &
                         share_sup, dq_opt, t_lt_c, pstar, gl_max_c,           &
                         gl_max_eff, r_bl_c, r_ca_c, g_oth_c )
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
CALL stom_opt_mod (                                                            &
        ! IN
          land_pts, ft, open_pts, open_index,                                  &
          pft_photo_model, veg_pts, veg_index,                                 &
          ca, psi_src, acr_c, apar_c, oa, vcmax_c,                             &
          kc_lt_c, ko_lt_c, ccp_lt_c, pstar,                                   &
          km_lt_c, dq_opt, je_c, t_lt_c, je_dummy,                             &
          fapar_dummy,                                                         &
          kmax_c, kcrit_c, gl_max_eff, ipar,                                   &
          l_multilayer, r_bl_c, r_ca_c, g_oth_c,                               &
          kmax_ref_c, conductance_b, conductance_c,                            &
          psi_leaf_extreme, psi_root_extreme, .TRUE.,                          &
        ! IN OUT
          rd_c,                                                                &
        ! OUT
          ci_c, anetl_c, el_c, flux_o3_c, fo3_c,                               &
          gl_c, psi_leaf_c, cg_c, hc_c, leaf_k_c,                              &
          psi_stem = psi_stem_c, lambda_eff = lam_c,                           &
          lai_cost = lai_c                                                     &
  )
IF ( l_store_cpl ) CALL clear_store_q_path()   ! MDK, 24/10/26

! Closed leaves get min_gl_pft from stom_opt_mod as a canopy value; share it
! by leaf area so the two classes together give big-leaf's (half each when
! there is no leaf area). The same for open leaves the optimiser closed (no
! flux, gl = min_gl_pft).
DO i = 1,clos_pts
  l = veg_index(clos_index(i))
  IF ( lai(l) > EPSILON(0.0) ) THEN
    gl_c(l) = gl_c(l) * lai_c(l) / lai(l)
  ELSE
    gl_c(l) = 0.5 * gl_c(l)
  END IF
END DO
DO i = 1,open_pts
  l = veg_index(open_index(i))
  IF ( el_c(l) <= 0.0 .AND. gl_c(l) <= min_gl_pft(ft) .AND.                    &
       lai(l) > EPSILON(0.0) )                                                 &
    gl_c(l) = gl_c(l) * lai_c(l) / lai(l)
END DO

END SUBROUTINE two_leaf_class

!-----------------------------------------------------------------------------
! Plant capacitance: the transpiration E' that brings the leaf to the
! optimiser's steady-state psi_leaf, psi_t (psi_t_tgt), by the end of the
! step, with the store at the stem node, at the lit points lit_tgt.
!
! The optimiser's secant plant conductance (E* = e_star_tgt),
!
!   k = E* / (psi_s - psi_t),
!
! is split into a root side and a leaf side (f = som_leaf_resist_frac),
!
!   k_r = k / (1 - f),
!   k_l = k / f,
!
! so that with the store at its steady state E' = E*. With E held over
! the step,
!
!   E' = [psi_s (1 - e) + psi_x0 e - psi_t] / [(1 - e) / k_r + 1 / k_l],
!
! where
!
!   e      = exp(-dt k_r / C),
!   psi_x0 = W / C (the store's psi at the start of the step).
!
! So E' > E* while the store is fuller than its steady state (morning),
! E' < E* while it refills (afternoon), and E' -> E* as dt >> C / k_r.
! The leaf ends the step at psi_t, the optimiser's own risk level, where the
! re-solve reaches E' (stom_opt_at_e: within the supply cap and ci <= ca,
! and E' = 0 gives its smallest sample, not zero) and the leak does not
! raise E above it.
!
! Points where the optimiser has shut the stomata (E* ~ 0, or a drop
! psi_s - psi_t below dpsi_tgt_min, where the secant k is meaningless) are
! left to the optimiser. Sets tgt_on, tgt_kr, tgt_kl, e_new_tgt (E', E*
! elsewhere) and r_tgt (E' / E*, 1 elsewhere).
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_target_e()

IMPLICIT NONE

INTEGER :: m, l
REAL(KIND=real_jlslsm) :: k_t
                            ! The optimiser's whole-plant secant conductance
                            ! (mol H2O m-2 ground s-1 Pa-1).
REAL(KIND=real_jlslsm) :: c_s
                            ! The store's capacitance C (mol H2O m-2 ground
                            ! Pa-1).
REAL(KIND=real_jlslsm) :: psi_x0
                            ! The store's (stem node) psi at the start of the
                            ! step (Pa).
REAL(KIND=real_jlslsm) :: e_x
                            ! exp(-dt k_r / C): the fraction of the store's
                            ! departure from its steady state left after one
                            ! step (-).
REAL(KIND=real_jlslsm) :: e_new
                            ! The store's transpiration E' (mol H2O m-2
                            ! ground s-1).
                            !
                            ! MDK, 23/10/26

! The share of the plant resistance on the leaf side of the store is
! f_leaf_store (som_leaf_resist_frac kept off 0 and 1, set with the C rule).
!
! MDK, 23/10/26
DO m = 1,veg_pts
  l = veg_index(m)

  ! By default the optimiser's transpiration stands (E' = E*).
  e_new_tgt(l) = e_star_tgt(l)
  r_tgt(l)     = 1.0

  ! Only lit points: at night there is no stomatal transpiration to set,
  ! and the store just refills (carry loop).
  IF ( .NOT. lit_tgt(l) ) CYCLE

  ! Not where the store is not yet set (first step, or after a restart from
  ! a dump without it) or there are (almost) no leaves.
  IF ( cut_w_store(l,ft) > 0.0 .OR. lai(l) < lai_cut_min ) CYCLE

  ! Not where the optimiser has (about) shut the stomata: the secant
  ! conductance below would be round-off.
  IF ( e_star_tgt(l) <= TINY(1.0_real_jlslsm) .OR.                             &
       psi_src(l) - psi_t_tgt(l) < dpsi_tgt_min ) CYCLE

  ! The optimiser's whole-plant secant conductance at its optimum (mol m-2
  ! s-1 Pa-1), split into the leaf side and the root side of the store, in
  ! series (1/k = 1/k_l + 1/k_r).
  k_t = e_star_tgt(l) / (psi_src(l) - psi_t_tgt(l))
  tgt_kl(l) = k_t / f_leaf_store
  tgt_kr(l) = k_t / (1.0 - f_leaf_store)

  ! The store: its capacitance C (mol m-2 Pa-1) and its psi at the start of
  ! the step, psi_x0 = W / C.
  c_s    = cut_c_store(l,ft)
  psi_x0 = cut_w_store(l,ft) / c_s

  ! The fraction of the store's departure from its steady state that is
  ! left after one step of refilling through the root side (time constant
  ! C / k_r).
  e_x    = EXP(-timestep * tgt_kr(l) / c_s)

  ! The transpiration that puts the leaf at psi_t at the end of the step:
  ! the stem psi at the end of the step, psi_x0 e + psi_s (1 - e)
  ! - E (1 - e) / k_r, less the drop E / k_l across the leaf side, set to
  ! psi_t and solved for E. No negative transpiration.
  e_new  = e_leaf_end( psi_src(l), psi_x0, e_x, tgt_kr(l), tgt_kl(l),          &
                       psi_t_tgt(l) )
  e_new  = MAX(e_new, 0.0)

  ! The store sets the transpiration here: E', and its ratio to the
  ! optimiser's E*, which scales each leaf class (or layer) before A, gl
  ! and ci are re-solved there.
  tgt_on(l)    = .TRUE.
  e_new_tgt(l) = e_new
  r_tgt(l)     = e_new / e_star_tgt(l)
END DO

END SUBROUTINE store_target_e

!-----------------------------------------------------------------------------
! Plant capacitance, two-leaf canopy.
!
! The plant has one store, so E' is set for the canopy (store_target_e):
!
!   E*    = E_sun + E_shd                  (the optimiser's transpiration),
!   psi_t = f_sun psi_sun + (1 - f_sun) psi_shd
!                                          (leaf-area-weighted psi_leaf),
!
! with f_sun = LAI_sun / LAI. Each class then keeps its share of E*:
!
!   E'_sun = E_sun E' / E*,
!   E'_shd = E_shd E' / E*,
!
! and its A, gl and ci are re-solved at that E (two_leaf_class_at_e).
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE store_target_2l()

IMPLICIT NONE

INTEGER :: m, l
REAL(KIND=real_jlslsm) :: f_s
                            ! Sunlit share of the leaf area, LAI_sun / LAI
                            ! (-).
                            !
                            ! MDK, 23/10/26

! The optimiser's canopy transpiration and leaf psi (the target), and the
! points with light on either class.
DO m = 1,veg_pts
  l = veg_index(m)
  e_star_tgt(l) = el_sun(l) + el_shd(l)
  f_s = 0.5
  IF ( lai(l) > EPSILON(0.0) ) f_s = lai_sun_2l(l) / lai(l)
  psi_t_tgt(l) = f_s * psi_leaf_sun(l) + (1.0 - f_s) * psi_leaf_shd(l)
  lit_tgt(l)   = apar_sun_2l(l) > 0.0 .OR. apar_shd_2l(l) > 0.0
END DO

! The canopy E' from the store (and tgt_on where it applies).
CALL store_target_e()

! Each class's target transpiration: its own optimal E where the store
! does not apply, else its share of E*, scaled to E'.
DO m = 1,veg_pts
  l = veg_index(m)
  e_tgt_sun(l) = el_sun(l)
  e_tgt_shd(l) = el_shd(l)
  IF ( tgt_on(l) ) THEN
    e_tgt_sun(l) = el_sun(l) * e_new_tgt(l) / e_star_tgt(l)
    e_tgt_shd(l) = el_shd(l) * e_new_tgt(l) / e_star_tgt(l)
  END IF
END DO
IF ( .NOT. ANY(tgt_on(:)) ) RETURN

! Re-solve A, gl and ci of each class at its target E: sunlit first, then
! shaded, which sees the new sunlit gl through the canopy coupling (as in
! the optimiser's own pass, two_leaf_class).
CALL two_leaf_class_at_e( apar_sun_2l, lai_sun_2l, kmax_ref_sun_2l,            &
                          qs_sun_lt, dq_sun_lt, gbw_sun_lt, gbw_shd_lt,        &
                          gl_shd, t_sun_lt, gl_max_sun_2l, acr_sun_2l,         &
                          vcmax_sun_2l, kc_sun_lt, ko_sun_lt, ccp_sun_lt,      &
                          km_sun_lt, je_sun, kmax_sun_2l, kcrit_sun_2l,        &
                          rd_sun, e_tgt_sun, ci_sun_2l, anetl_sun, el_sun,     &
                          gl_sun )
CALL two_leaf_class_at_e( apar_shd_2l, lai_shd_2l, kmax_ref_shd_2l,            &
                          qs_shd_lt, dq_shd_lt, gbw_shd_lt, gbw_sun_lt,        &
                          gl_sun, t_shd_lt, gl_max_shd_2l, acr_shd_2l,         &
                          vcmax_shd_2l, kc_shd_lt, ko_shd_lt, ccp_shd_lt,      &
                          km_shd_lt, je_shd, kmax_shd_2l, kcrit_shd_2l,        &
                          rd_shd, e_tgt_shd, ci_shd_2l, anetl_shd, el_shd,     &
                          gl_shd )

! The canopy conductance from the re-solved classes.
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( tgt_on(l) ) gc(l) = gl_sun(l) + gl_shd(l)
END DO

END SUBROUTINE store_target_2l

!-----------------------------------------------------------------------------
! Plant capacitance: one leaf class (sunlit or shaded) re-solved at the
! target transpiration e_tgt_c set by the store.
!
! At the points where the store sets the transpiration (tgt_on), the
! class's
!
!   ci, A (anetl_c), gl and E (el_c)
!
! are those at which it transpires e_tgt_c (stom_opt_at_e), with the same
! canopy coupling, soil-supply cap and photosynthesis inputs (at the leaf
! temperature with l_leaf_temp) as the optimiser's pass, two_leaf_class.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE two_leaf_class_at_e( apar_c, lai_c, kmax_ref_c, qs_lt_c, dq_lt_c,   &
                                gbw_c, gbw_o, gl_o, t_lt_c, gl_max_c, acr_c,   &
                                vcmax_c, kc_lt_c, ko_lt_c, ccp_lt_c, km_lt_c,  &
                                je_c, kmax_c, kcrit_c, rd_c, e_tgt_c, ci_c,    &
                                anetl_c, el_c, gl_c )
IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  apar_c(land_pts), lai_c(land_pts), kmax_ref_c(land_pts),                     &
  qs_lt_c(land_pts), dq_lt_c(land_pts), gbw_c(land_pts), gbw_o(land_pts),      &
  gl_o(land_pts), t_lt_c(land_pts), gl_max_c(land_pts), acr_c(land_pts),       &
  vcmax_c(land_pts), kc_lt_c(land_pts), ko_lt_c(land_pts),                     &
  ccp_lt_c(land_pts), km_lt_c(land_pts), je_c(land_pts), kmax_c(land_pts),     &
  kcrit_c(land_pts), rd_c(land_pts)
                            ! The class's inputs, as two_leaf_class.
REAL(KIND=real_jlslsm), INTENT(IN) :: e_tgt_c(land_pts)
                            ! The class's target transpiration (mol H2O m-2
                            ! ground s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: ci_c(land_pts)
                            ! Intercellular CO2 (Pa).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: anetl_c(land_pts)
                            ! Net assimilation (mol CO2 m-2 ground s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: el_c(land_pts)
                            ! Transpiration (mol H2O m-2 ground s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: gl_c(land_pts)
                            ! Conductance (m s-1). Points the store does not
                            ! set keep the optimiser's values.
                            !
                            ! MDK, 23/10/26

INTEGER :: i, l, m

! The open points of this class where the store sets the transpiration.
! (This reuses the host's open_pts and open_index, as two_leaf_class does:
! the two-leaf lit set for the leak and the impairment is rebuilt after it.)
!
! MDK, 23/10/26
open_pts = 0
DO i = 1,veg_pts
  l = veg_index(i)
  IF ( tgt_on(l) .AND. apar_c(l) > 0.0 .AND. lai_c(l) > EPSILON(0.0) ) THEN
    open_pts = open_pts + 1
    open_index(open_pts) = i
  END IF
END DO
IF ( open_pts == 0 ) RETURN

! This class's share of the soil supply and of the soil-to-root
! conductance: its share of the undamaged plant conductance (as
! two_leaf_class).
DO m = 1,veg_pts
  l = veg_index(m)
  share_sup(l) = kmax_ref_c(l) /                                               &
                 MAX(kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l), TINY(1.0))
END DO

! The deficit and resistances that couple this class's transpiration to
! the canopy air and the other class (gl_o).
CALL coupled_inputs( land_pts, veg_pts, veg_index, l_leaf_temp, dq_min,        &
                     dqc, dq, qs_lt_c, q1, dq_lt_c, gbw_c, gbw_o, ra_ca_lt,    &
                     gl_o, g_soil_c, qs, dq_opt, r_bl_c, r_ca_c, g_oth_c )

! The soil-supply cap on gl, then ci, A, gl and E at the target E: at
! the leaf temperature (l_leaf_temp), else at tstar (the class arrays hold
! the tile's values then).
!
! MDK, 24/10/26
CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,               &
                         share_sup, dq_opt, t_lt_c, pstar, gl_max_c,           &
                         gl_max_eff, r_bl_c, r_ca_c, g_oth_c )
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
CALL stom_opt_at_e( land_pts, ft, open_pts, open_index, pft_photo_model,       &
                    veg_index, rd_c, ca, psi_src, acr_c, apar_c, oa,           &
                    vcmax_c, kc_lt_c, ko_lt_c, ccp_lt_c, pstar, km_lt_c,       &
                    dq_opt, je_c, t_lt_c, je_dummy, fapar_dummy,               &
                    ipar, kmax_c, kcrit_c, gl_max_eff, l_multilayer,           &
                    r_bl_c, r_ca_c, g_oth_c, kmax_ref_c, conductance_b,        &
                    conductance_c, psi_leaf_extreme, psi_root_extreme,         &
                    .TRUE., e_tgt_c, ci_c, anetl_c, gl_c, el_c )

END SUBROUTINE two_leaf_class_at_e

!-----------------------------------------------------------------------------
! Leaf energy balance (l_leaf_temp), two-leaf and multilayer: set-up before
! the passes.
!
! 1. The radiation and conductance terms that do not change between
!    passes, from the class leaf areas and absorbed PAR (lai_sun_2l,
!    lai_shd_2l, apar_sun_2l, apar_shd_2l), for each layer as well with
!    ml_leaf_temp_model = 2.
! 2. First guesses: the leaf temperatures (last timestep's, or tstar on
!    the first timestep), the canopy air temperature (tair) and humidity,
!    and the canopy air resistances.
! 3. The photosynthesis temperature parameters of each class start from the
!    tile's.
!
! (Moved here from the two-leaf code unchanged.)
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE lt_init()
IMPLICIT NONE
INTEGER :: l, m
    ! Its own loop counters, not the host's. MDK, 24/10/26

  ! The wind extinction: leaf_wind_ext (0 by default) for leaf_aero_model
  ! 0 and 2, CABLE's coexp for 1 (MDK, 19/10/26).
  wind_ext(:) = leaf_wind_ext
  zref_a(:)   = 0.0
  IF ( leaf_aero_model == 1 ) THEN
    CALL leaf_aero_geom( land_pts, veg_pts, veg_index, lai, canht, z1_wind,    &
                         z0m_a, zref_a, rt1ab_a, zrd_a, usc_a, usuh_a,         &
                         wind_ext )
  END IF
  IF ( l_lt_lyr ) THEN
    ! Also the terms of every layer's sunlit and shaded leaf.
    ! MDK, 19/10/26
    CALL leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,      &
                          dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,      &
                          wind_ext,                                            &
                          apar_sun_2l, apar_shd_2l,                            &
                          sw_sun_lt, sw_shd_lt, wlw_sun_lt, wlw_shd_lt,        &
                          wlws_sun_lt, wlws_shd_lt, frad_sun_lt,               &
                          frad_shd_lt, ff_sun_lt, ff_shd_lt,                   &
                          fapar_sun_vis = fapar_sun,                           &
                          fapar_shd_vis = fapar_shd, lai_c = lai_lt_c,         &
                          sw_c = sw_lt_c, wlw_c = wlw_lt_c,                    &
                          wlws_c = wlws_lt_c, frad_c = frad_lt_c,              &
                          ff_c = ff_lt_c )
  ELSE
  CALL leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,        &
                        dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,        &
                        wind_ext,                                              &
                        apar_sun_2l, apar_shd_2l,                              &
                        sw_sun_lt, sw_shd_lt, wlw_sun_lt, wlw_shd_lt,          &
                        wlws_sun_lt, wlws_shd_lt, frad_sun_lt, frad_shd_lt,    &
                        ff_sun_lt, ff_shd_lt )
  END IF
  ! First guesses: the leaf-area-weighted leaf temperature (K) and the
  ! stability parameter of the canopy air (-).
  t_g_lt(:)  = tstar(:)
  zeta_lt(:) = 0.0
  ! Defined at every land point (leaf_class_at_temp takes qsat of the
  ! whole array, and the per-layer classes copy it).
  ! MDK, 19/10/26
  t_sun_lt(:) = tstar(:)
  t_shd_lt(:) = tstar(:)
  DO m = 1,veg_pts
    l = veg_index(m)
    ! Last timestep's leaf temperatures (K), or tstar where there are none
    ! yet (below 150 K: unset).
    IF ( t_leaf_sun(l) < 150.0 .OR. t_leaf_shd(l) < 150.0 ) THEN
      t_sun_lt(l) = tstar(l)
      t_shd_lt(l) = tstar(l)
    ELSE
      t_sun_lt(l) = t_leaf_sun(l)
      t_shd_lt(l) = t_leaf_shd(l)
    END IF
    ! Their leaf-area-weighted mean (K).
    IF ( lai(l) > EPSILON(0.0) ) THEN
      t_g_lt(l) = ( lai_sun_2l(l) * t_sun_lt(l)                                &
                    + lai_shd_2l(l) * t_shd_lt(l) ) / lai(l)
    END IF
    ! Canopy air - level 1 resistance (s m-1; leaf_aero_model).
    SELECT CASE ( leaf_aero_model )
    CASE ( 1 )
      ! CABLE's rt1, computed in leaf_temp_update (ra where there is no
      ! canopy); here neutral, for the first guess of q_c.
      IF ( zref_a(l) > 0.0 ) THEN
        ra_lt(l) = MAX( 5.0, rt1ab_a(l) * LOG(zref_a(l) / z0m_a(l))            &
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
    ! Canopy air temperature (K): the air temperature to start.
    t_c_lt(l) = tair(l)
  END DO
  ! Canopy air humidity from last timestep's conductance and leaf
  ! temperatures (q1 for leaf_aero_model = 0). (gc on entry is last
  ! timestep's tile conductance, which already includes the soil's.)
  CALL qsat(qs_g_lt, t_g_lt, pstar, land_pts)
  DO m = 1,veg_pts
    l = veg_index(m)
    q_c_lt(l) = ( q1(l) + ra_lt(l) * gc(l) * qs_g_lt(l) )                      &
                / ( 1.0 + ra_lt(l) * gc(l) )
  END DO
  IF ( leaf_aero_model == 1 ) THEN
    ! The neutral rt1 above is for the first guess only: leaf_temp_update
    ! takes ra where there is no canopy. Keep it as the first guess of
    ! the canopy air resistance of the coupled transpiration.
    DO m = 1,veg_pts
      l = veg_index(m)
      ra_ca_lt(l) = ra_lt(l)
      ra_lt(l) = ra(l)
    END DO
  END IF
  ! Respiration below uses rd and the classes start from it. The CO2
  ! compensation point (Pa) and the Michaelis-Menten constants (Pa) of each
  ! class start from the tile's (they are reset at the leaf temperatures
  ! each pass, leaf_class_at_temp).
  ccp_sun_lt(:) = ccp(:)
  ccp_shd_lt(:) = ccp(:)
  kc_sun_lt(:)  = kc(:)
  kc_shd_lt(:)  = kc(:)
  ko_sun_lt(:)  = ko(:)
  ko_shd_lt(:)  = ko(:)
  km_sun_lt(:)  = km(:)
  km_shd_lt(:)  = km(:)

END SUBROUTINE lt_init

!-----------------------------------------------------------------------------
! Leaf energy balance (l_leaf_temp), two-leaf and multilayer, at point lp
! after the passes: the leaf temperatures used (t_leaf_sun, t_leaf_shd,
! t_can, with f_sun the sunlit fraction of the leaf area), the fraction of
! the classes held at the dt_leaf limit, and with l_leaf_temp_gc_eq the
! transpiration-equivalent gc, capped at gc_cap (set by the caller).
! (Moved here from the two-leaf code unchanged, except t_leaf_clamp.)
! MDK, 18/10/26
!-----------------------------------------------------------------------------
SUBROUTINE lt_finish( lp, f_sun, n_cls_in )
IMPLICIT NONE
INTEGER, INTENT(IN) :: lp
REAL(KIND=real_jlslsm), INTENT(IN) :: f_sun
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: n_cls_in
    ! Number of leaf classes with leaf area, for t_leaf_clamp (2 at most
    ! when absent; with ml_leaf_temp_model = 2 those of all the layers).
    ! MDK, 19/10/26
REAL(KIND=real_jlslsm) :: n_cls
REAL(KIND=real_jlslsm) :: dq_eq, cap_eq
    ! gc_eq: the tile's actual deficit qs(tstar) - q1 (kg kg-1) and the
    ! cap on gc (m s-1). MDK, 24/10/26
CHARACTER(LEN=300) :: clamp_msg

  ! Leaf temperatures used. A class with no leaf area takes the other's.
  IF ( lai_sun_2l(lp) > EPSILON(0.0) ) THEN
    t_leaf_sun(lp) = t_sun_lt(lp)
  ELSE
    t_leaf_sun(lp) = t_shd_lt(lp)
  END IF
  IF ( lai_shd_2l(lp) > EPSILON(0.0) ) THEN
    t_leaf_shd(lp) = t_shd_lt(lp)
  ELSE
    t_leaf_shd(lp) = t_sun_lt(lp)
  END IF
  t_can(lp) = f_sun * t_leaf_sun(lp) + (1.0 - f_sun) * t_leaf_shd(lp)

  ! The surface energy balance applies gc with qs(tstar) - q1 through
  ! ra. Return the gc that gives the transpiration the optimisation
  ! chose (el, at the leaf temperatures) there:
  !   rho (qs(tstar) - q1) / (ra + 1/gc) = el * m_h2o,
  ! capped at the canopy som_gl_max. The soil conductance g_s shares ra
  ! in the tile flux (soil_evap adds it to gc) and the transpiration is
  ! the stomatal share,
  !   rho (qs(tstar) - q1) gc / (1 + ra (gc + g_s)) = el * m_h2o,
  ! so 1/gc is the r_eq above divided by 1 + ra g_s.
  !
  ! The deficit is the actual qs(tstar) - q1, not dq, which the stomatal
  ! optimisation floors at dq_min: with the floor a near-saturated or
  ! dewy step gave a huge gc that JULES then applied to the real
  ! deficit. Below dq_min (or with condensation) the optimised gc is
  ! kept. The cap is always finite: som_gl_max x LAI where the canopy
  ! cap is off (supply-loss).
  !
  ! MDK, 24/10/26
  ! Diagnostics (canopy_diag_mod): the stomatal gc before gc_eq replaces
  ! it, and the canopy air and leaf fluxes of the energy balance.
  ! MDK, 24/10/26
  gc_stomatal_pft(lp,ft) = gc(lp)
  t_can_air_pft(lp,ft)   = t_c_lt(lp)
  q_can_air_pft(lp,ft)   = q_c_lt(lp)
  h_leaf_pft(lp,ft)      = h_leaf_pts(lp)
  le_leaf_pft(lp,ft)     = le_leaf_pts(lp)

  IF ( l_leaf_temp_gc_eq ) THEN
    dq_eq  = qs(lp) - q1(lp)
    cap_eq = gc_cap
    IF ( cap_eq <= 0.0 ) cap_eq = som_gl_max * lai(lp)
    IF ( el(lp) * m_h2o_lt > TINY(1.0) .AND. dq_eq > 1.0e-4 ) THEN
      r_eq = pstar(lp) / ( r * tstar(lp) ) * dq_eq / ( el(lp) * m_h2o_lt )     &
             - ra(lp)
      r_eq = r_eq / ( 1.0 + ra(lp) * g_soil_c(lp) )
      IF ( r_eq > 0.0 ) THEN
        gc(lp) = MIN(1.0 / r_eq, cap_eq)
      ELSE
        ! Demand above what any gc gives: the cap.
        gc(lp) = cap_eq
      END IF
    END IF
  END IF

  ! Fraction of the classes with leaf area held at the dt_leaf limit.
  ! MDK, 18/10/26
  n_cls = 0.0
  IF ( lai_sun_2l(lp) > EPSILON(0.0) ) n_cls = n_cls + 1.0
  IF ( lai_shd_2l(lp) > EPSILON(0.0) ) n_cls = n_cls + 1.0
  IF ( PRESENT(n_cls_in) ) n_cls = n_cls_in
  IF ( n_cls > 0.0 ) t_leaf_clamp(lp) = MIN(n_clamp_lt(lp), n_cls) / n_cls
  IF ( n_clamp_lt(lp) > 0.0 .AND. n_lt_clamp_warn < max_lt_clamp_warn ) THEN
!$OMP CRITICAL (sf_stom_lt_clamp_warn)
    n_lt_clamp_warn = n_lt_clamp_warn + 1
    WRITE(clamp_msg, '(A,I0,A,I0,A,3F8.2,A,I0,A)')                             &
      'l_leaf_temp: leaf temperature held at dt_leaf from the canopy air, '//  &
      'PFT ', ft, ', point ', lp, ', T_sun, T_shd, T_c (K) ',                  &
      t_sun_lt(lp), t_shd_lt(lp), t_c_lt(lp),                                  &
      ' (see t_leaf_clamp_pft; at most ', max_lt_clamp_warn,                   &
      ' such warnings)'
    errcode = -1
    CALL ereport(RoutineName, errcode, TRIM(clamp_msg))
!$OMP END CRITICAL (sf_stom_lt_clamp_warn)
  END IF

END SUBROUTINE lt_finish

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

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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
! MDK, 19/10/26

USE pftparm, ONLY: orient

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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

! MDK, 19/10/26
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  fapar_sun_vis(land_pts,ilayers), fapar_shd_vis(land_pts,ilayers)
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  lai_c(land_pts,2*ilayers), sw_c(land_pts,2*ilayers),                         &
  wlw_c(land_pts,2*ilayers), wlws_c(land_pts,2*ilayers),                       &
  frad_c(land_pts,2*ilayers), ff_c(land_pts,2*ilayers)

REAL(KIND=real_jlslsm) :: kd, xi0, xi1, dtop, dbot, dwind
INTEGER :: l, m, n
LOGICAL :: l_lyr


! Diffuse (long-wave) extinction coefficient: spherical leaves (Wang &
! Leuning, 1998), or horizontal.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_SETUP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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

SUBROUTINE leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_leaf, oa, actj, actv, dsj,   &
                               dsv, jv25, nleaf_top, nw, i2c, pstar, q_c,      &
                               dq_min, ccp, kc, ko, km, vcmax_c, je_c, rd_c,   &
                               qs_c, dq_c, f_vc )

! Leaf energy balance (l_leaf_temp; two-leaf and multilayer): the
! photosynthesis parameters
! of one leaf class at its temperature t_leaf - the temperature responses,
! the class totals of Vcmax and Rd (N-weighted leaf area nw), the electron
! flux - and its humidity deficit to the canopy air, qs(t_leaf) - q_c.

USE jules_vegetation_mod, ONLY: photo_farquhar
USE qsat_mod, ONLY: qsat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_CLASS_AT_TEMP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_class_at_temp

!#############################################################################
!#############################################################################

SUBROUTINE leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             ra, lw_down, u_wind, lai_sun, lai_shd,            &
                             sw_sun, sw_shd, wlw_sun, wlw_shd, wlws_sun,       &
                             wlws_shd, frad_sun, frad_shd, ff_sun, ff_shd,     &
                             z0m_a, zref_a, rt1ab_a, zrd_a, usc_a, usuh_a,     &
                             gl_sun, gl_shd, g_soil, qs_soil, zeta_a, t_sun,   &
                             t_shd, t_c, q_c, gbw_sun, gbw_shd, ra_out,        &
                             n_clamp )

! Two-leaf leaf energy balance (l_leaf_temp): sunlit and shaded leaf
! temperatures, and canopy air temperature and humidity, for the stomatal
! conductances of the last pass. The balance of leaf_temp_update_n for
! the two classes (sunlit, shaded).
! MDK, 19/10/26

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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
                         q_c, gbw_x, ra_out, n_clamp )

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
                               n_clamp )

! Leaf energy balance (l_leaf_temp) of nc leaf classes sharing one canopy
! air: two-leaf (nc = 2, sunlit and shaded), the aggregate multilayer canopy
! (the same) or the multilayer canopy layer by layer (nc = 2 x ilayers, the
! sunlit and shaded leaf of every layer; ml_leaf_temp_model = 2, as
! MAESPA). Leaf temperatures, and canopy air temperature and humidity, for
! the stomatal conductances of the last pass. (Generalised from the
! two-leaf balance to nc classes, MDK, 19/10/26.)
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
!   q_c = q1 + ra (LE_sun + LE_shd + LE_soil) / (rho lc),
!   t_c = tair + ra (H_sun + H_shd) / (rho cp),
! LE_soil = rho lc g_soil (qs_soil - q_c) being the soil evaporation into
! the canopy air, as in CABLE (soil sensible heat is not included: there is
! no soil surface temperature separate from tstar on vegetated tiles).
! The balance is solved for these conductances: LE is linear in q_c, so q_c
! is explicit for a given t_c, and t_c is found by regula falsi within
! dt_can of tair (at once when ra = 0). Leaf temperatures are kept
! within dt_leaf of t_c.

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE csigma, ONLY: sbcon
USE canopy_diag_mod, ONLY: h_leaf_pts, le_leaf_pts   ! MDK, 24/10/26
USE planet_constants_mod, ONLY: cp, r, repsilon
USE water_constants_mod, ONLY: lc, Water_kg_per_mol
USE jules_vegetation_mod, ONLY: leaf_width, leaf_shelter, leaf_aero_model
USE qsat_mod, ONLY: qsat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts), nc
      ! nc: number of leaf classes. MDK, 19/10/26

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  tair(land_pts), q1(land_pts), pstar(land_pts), ra(land_pts),                 &
  lw_down(land_pts), u_wind(land_pts),                                         &
  lai_c(land_pts,nc), sw_c(land_pts,nc), wlw_c(land_pts,nc),                   &
  wlws_c(land_pts,nc), frad_c(land_pts,nc), ff_c(land_pts,nc),                 &
  gl_c(land_pts,nc),                                                           &
      ! Of each class (see leaf_temp_setup; the sunlit and shaded leaf of

      ! two-leaf were lai_sun, sw_sun, ... gl_sun and lai_shd, ... gl_shd):
      ! leaf area, absorbed short-wave, long-wave weights, radiative
      ! conductance scalar, wind-weighted leaf area and stomatal
      ! conductance (m s-1, ground). MDK, 19/10/26
  z0m_a(land_pts), zref_a(land_pts), rt1ab_a(land_pts), zrd_a(land_pts),       &
  usc_a(land_pts), usuh_a(land_pts),                                           &
      ! CABLE canopy geometry (leaf_aero_model = 1; see leaf_aero_geom).
  g_soil(land_pts), qs_soil(land_pts)
      ! Soil evaporation conductance into the canopy air (m s-1, ground
      ! area) and the saturated humidity it evaporates at (kg kg-1).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  zeta_a(land_pts),                                                            &
      ! zref/L (leaf_aero_model = 1): the start of the Monin-Obukhov
      ! iteration, and the last iterate on output, so successive passes in
      ! a timestep continue the iteration (0 for the first).
  t_cls(land_pts,nc), t_c(land_pts), q_c(land_pts)
      ! (t_cls: leaf temperature of each class; MDK, 19/10/26)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  gbw_c(land_pts,nc), ra_out(land_pts)
      ! Leaf boundary layer conductance for water of each class (m/s, class
      ! basis; 0 without leaf area) and the canopy air resistance (s/m) of
      ! the balance found, for l_som_coupled_e. Points not solved keep
      ! their values.

REAL(KIND=real_jlslsm), INTENT(IN OUT) :: n_clamp(land_pts)
      ! Number of leaf classes (0-nc) whose temperature was held at dt_leaf
      ! from the canopy air in the balance found.
      ! MDK, 18/10/26

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
  m_h2o = Water_kg_per_mol,                                                    &
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
  tx_v(nc,veg_pts), qc_v(veg_pts), h_v(veg_pts), le_v(veg_pts),                &
  hb_v(veg_pts), leb_v(veg_pts), res, ra_m(land_pts),                          &
  u_m(land_pts), us, rt1usc, zr, gbw_last(nc), ncl_last
      ! (ncl_last: classes held at dt_leaf in the last canopy_balance call;
      ! MDK, 18/10/26. tx_v, gbw_last of each class; MDK, 19/10/26)
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
! heat and the buoyancy of the latent heat, -k g zref (H + 0.07 LE) /
! (rho cp T u*^3). Only the leaves' fluxes: the soil's sensible heat is not
! known here (nor is its latent heat added).
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

DO m = 1,veg_pts
  l = veg_index(m)
  ra_out(l) = ra_m(l)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

SUBROUTINE set_out( l, m, tc, qc )
! (The class temperatures from tx_v; MDK, 19/10/26.)
INTEGER, INTENT(IN) :: l, m
REAL(KIND=real_jlslsm), INTENT(IN) :: tc, qc
t_c(l)   = tc
t_cls(l,:) = tx_v(:,m)
q_c(l)   = qc
! (gbw_last is from the canopy_balance call just made for this point.)
gbw_c(l,:) = gbw_last(:)
h_v(m)   = hb_v(m)
le_v(m)  = leb_v(m)
! For the canopy diagnostics (canopy_diag_mod). MDK, 24/10/26
IF ( ALLOCATED(h_leaf_pts) ) THEN
  h_leaf_pts(l)  = hb_v(m)
  le_leaf_pts(l) = leb_v(m)
END IF
! (ncl_last too.) MDK, 18/10/26
n_clamp(l) = ncl_last
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

! Leaf temperatures and canopy humidity for canopy air temperature tc, and
! the residual of the canopy air sensible heat balance,
!   res = tc - tair - ra (sum of the classes' H) / (rho cp).
! (nc classes; MDK, 19/10/26)

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

! Canopy humidity: q_c = q1 + c_q sum(LE), c_q = ra / (rho L), with L the
! latent heat at tc that the leaves' LE uses (lhv / m_h2o, J kg-1; it was
! the constant lc, 2-3 % off at 25-35 C). The soil's vapour flux is
! rho g_soil (qs_soil - q_c), so ra multiplies g_soil.
!
! MDK, 24/10/26
c_q    = ra_x * m_h2o / ( rho * lhv )
! Only the classes whose LE is positive at the solved q_c take part, as the
! leaf heat balance below clips each class's LE at 0: a class that would
! condense (at night) adds no vapour to, and removes none from, the canopy
! air. Classes are dropped and q_c re-solved until none of those left is
! negative (at most nc passes; with none negative this is the single solve
! as before).
! MDK, 22/10/26
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
  IF ( lai_x(ic) > EPSILON(0.0) ) THEN
    le     = MAX(0.0, a(ic) + b(ic) * (qs0 - qc_new))
    le_out = le_out + le
    dt     = ( rn(ic) - le ) / ( cp_mol * gh(ic) )
    ! Count the classes held at the limit (t_leaf_clamp). MDK, 18/10/26
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

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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
REAL(KIND=real_jlslsm) :: jmax_numerator, kc_val, ko_val, power, tau, tdegc,   &
                          t_minus_ref, t_term, vcmax_numerator

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_TEMP_RESPONSES'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

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
! fixed choice, not a parameter. The state is dumped (fvcmax_psi_state), so
! a restart continues it; a run without it starts at the instant value.
! MDK, 24/10/26

USE pftparm, ONLY: psi_vcmax_f, sf_vcmax, psi_vcmax_fmin
USE timestep_mod, ONLY: timestep

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='VCMAX_PSI_FACTOR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL sf_stom_state_alloc( land_pts )   ! MDK, 22/10/26
w_rec = 1.0 - EXP( -timestep / tau_rec )

DO m = 1,veg_pts
  l = veg_index(m)
  psi_mpa = MIN(psi_rz(l), 0.0) * 1.0e-6
  f_psi = ( 1.0 + EXP( sf_vcmax(ft) * psi_vcmax_f(ft) * 1.0e-6 ) )             &
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

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
                               dq, tstar, pstar, gl_max_in, gl_max_out,        &
                               r_bl, r_ca, g_o )

USE c_rmol, ONLY: rmol
USE planet_constants_mod, ONLY: repsilon
USE jules_vegetation_mod, ONLY: l_som_supply_limit, l_som_coupled_e
USE stom_opt_jls_mod, ONLY: g_of_geff
USE water_constants_mod, ONLY: Water_kg_per_mol   ! MDK, 24/10/26

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
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

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  r_bl(land_pts), r_ca(land_pts), g_o(land_pts)
      ! l_som_coupled_e: the resistances and other-class conductance of
      ! stom_opt_mod, so that the cap is the g whose coupled E is the supply.

REAL(KIND=real_jlslsm), PARAMETER :: m_h2o = Water_kg_per_mol
      ! Molar mass of water (kg mol-1). MDK, 24/10/26

REAL(KIND=real_jlslsm) :: gl_sup
INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='APPLY_SUPPLY_LIMIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

gl_max_out(:) = gl_max_in(:)
IF ( .NOT. l_som_supply_limit ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO m = 1,veg_pts
  l = veg_index(m)
  IF ( e_supply(l) >= 0.0 .AND. dq(l) > 0.0 ) THEN
    gl_sup = e_supply(l) * share(l) / m_h2o                                    &
             * repsilon * rmol * tstar(l) / ( dq(l) * pstar(l) )
    IF ( l_som_coupled_e .AND. PRESENT(r_ca) )                                 &
      gl_sup = g_of_geff(gl_sup, r_bl(l), r_ca(l), g_o(l))
    IF ( gl_max_in(l) > 0.0 ) gl_sup = MIN(gl_sup, gl_max_in(l))
    gl_max_out(l) = MAX(gl_sup, TINY(1.0_real_jlslsm))
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE apply_supply_limit

!#############################################################################
! Coupled transpiration (l_som_coupled_e), two-leaf: the inputs of one leaf
! class's stom_opt_mod call. Without it, dq_opt is the deficit used before
! (dq_lt with l_leaf_temp, else dqc) and the resistances are 0. With it, the
! deficit is from the leaf to the level-1 air (qs(T_leaf) - q1 with
! l_leaf_temp, else dq), r_bl the class's leaf boundary layer (l_leaf_temp
! only), r_ca the canopy air resistance and g_oth the other class's
! stomatal conductance in series with its boundary layer.
!
! The soil evaporation beneath the canopy (conductance g_s, at
! qs_soil = qs(tstar)) is a further path into the canopy air. Solving the
! canopy air humidity with the leaf, the other class (taken at this leaf's
! qs, as before) and the soil gives
!   E = g' (D + r_ca g_s (qs_lt - qs_soil)) / (1 + r_ca (g' + g'_other + g_s)),
! so g_s is added to g_oth and r_ca g_s (qs_lt - qs_soil) to the deficit,
! which does not depend on g' (0 without l_leaf_temp, where the leaf is at
! tstar).
!#############################################################################
SUBROUTINE coupled_inputs( land_pts, veg_pts, veg_index, l_lt, dq_min, dqc,    &
                           dq, qs_lt, q1, dq_lt, gbw_this, gbw_other, r_ca_in, &
                           gl_other, g_soil, qs_soil, dq_opt, r_bl, r_ca,      &
                           g_oth )

USE jules_vegetation_mod, ONLY: l_som_coupled_e

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
LOGICAL, INTENT(IN) :: l_lt
REAL(KIND=real_jlslsm), INTENT(IN) :: dq_min
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  dqc(land_pts), dq(land_pts), qs_lt(land_pts), q1(land_pts), dq_lt(land_pts), &
  gbw_this(land_pts), gbw_other(land_pts), r_ca_in(land_pts),                  &
  gl_other(land_pts), g_soil(land_pts), qs_soil(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  dq_opt(land_pts), r_bl(land_pts), r_ca(land_pts), g_oth(land_pts)

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: r_bl_o, g_o


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='COUPLED_INPUTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( l_lt ) THEN
  dq_opt(:) = dq_lt(:)
ELSE
  dq_opt(:) = dqc(:)
END IF
r_bl(:)  = 0.0
r_ca(:)  = 0.0
g_oth(:) = 0.0
IF ( .NOT. l_som_coupled_e ) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO m = 1,veg_pts
  l = veg_index(m)
  r_bl_o = 0.0
  IF ( l_lt ) THEN
    dq_opt(l) = MAX(dq_min, qs_lt(l) - q1(l))
    IF ( gbw_this(l) > 0.0 )  r_bl(l) = 1.0 / gbw_this(l)
    IF ( gbw_other(l) > 0.0 ) r_bl_o  = 1.0 / gbw_other(l)
  ELSE
    dq_opt(l) = dq(l)
  END IF
  r_ca(l)  = r_ca_in(l)
  g_o      = MAX(gl_other(l), 0.0)
  g_oth(l) = g_o / ( 1.0 + g_o * r_bl_o )
  g_oth(l) = g_oth(l) + g_soil(l)
  IF ( l_lt ) dq_opt(l) = MAX(dq_min, dq_opt(l)                                &
                              + r_ca(l) * g_soil(l) * (qs_lt(l) - qs_soil(l)))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE coupled_inputs

!#############################################################################
! First guesses for coupled transpiration in a timestep's first two-leaf
! pass: the last timestep's class conductances, boundary layers and canopy
! air resistance, or, before the first, the canopy conductance shared by
! leaf area, no boundary layer and a canopy air resistance of -1 (unknown).
!#############################################################################
SUBROUTINE cpl_state_load( ft, land_pts, veg_pts, veg_index, lai, lai_sun,     &
                           lai_shd, gc, gl_sun, gl_shd, gbw_sun, gbw_shd,      &
                           ra_ca )


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), lai_sun(land_pts),        &
                                      lai_shd(land_pts), gc(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gl_sun(land_pts), gl_shd(land_pts),     &
                                       gbw_sun(land_pts), gbw_shd(land_pts),   &
                                       ra_ca(land_pts)
INTEGER :: l, m
REAL(KIND=real_jlslsm) :: f_sun


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CPL_STATE_LOAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL sf_stom_state_alloc( land_pts )   ! MDK, 22/10/26
gl_sun(:)  = 0.0
gl_shd(:)  = 0.0
gbw_sun(:) = 0.0
gbw_shd(:) = 0.0
ra_ca(:)   = -1.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( cpl_state(l,ft,1) >= 0.0 ) THEN
    gl_sun(l)  = cpl_state(l,ft,1)
    gl_shd(l)  = cpl_state(l,ft,2)
    gbw_sun(l) = cpl_state(l,ft,3)
    gbw_shd(l) = cpl_state(l,ft,4)
    ra_ca(l)   = cpl_state(l,ft,5)
  ELSE
    f_sun = 0.5
    IF ( lai(l) > EPSILON(0.0) ) f_sun = lai_sun(l) / lai(l)
    gl_sun(l) = MAX(gc(l), 0.0) * f_sun
    gl_shd(l) = MAX(gc(l), 0.0) * (1.0 - f_sun)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cpl_state_load

SUBROUTINE cpl_state_save( ft, land_pts, veg_pts, veg_index, gl_sun, gl_shd,   &
                           gbw_sun, gbw_shd, ra_ca )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
! (DrHook: MDK, 24/10/26)
IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: gl_sun(land_pts), gl_shd(land_pts),      &
                                      gbw_sun(land_pts), gbw_shd(land_pts),    &
                                      ra_ca(land_pts)
INTEGER :: l, m


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CPL_STATE_SAVE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO m = 1,veg_pts
  l = veg_index(m)
  cpl_state(l,ft,1) = MAX(gl_sun(l), 0.0)
  cpl_state(l,ft,2) = MAX(gl_shd(l), 0.0)
  cpl_state(l,ft,3) = MAX(gbw_sun(l), 0.0)
  cpl_state(l,ft,4) = MAX(gbw_shd(l), 0.0)
  cpl_state(l,ft,5) = MAX(ra_ca(l), 0.0)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cpl_state_save

END MODULE sf_stom_mod
