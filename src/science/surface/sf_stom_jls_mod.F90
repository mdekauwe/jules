! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE sf_stom_mod

USE um_types, ONLY: real_jlslsm
USE sf_stom_state_mod, ONLY: f_vcmax_psi_state, cpl_state, sf_stom_state_alloc
                            ! f_vcmax_psi_state (l_som_vcmax_psi) and cpl_state
                            ! (l_som_coupled_e, two-leaf) are in
                            ! sf_stom_state_mod, so they go in the dump.
                            !
                            ! MDK, 22/10/26
USE sf_stom_state_mod, ONLY: shared_e_shd_prev, shared_amax_shd_prev
                            ! The shared root and stem (l_som_shared_stem): the
                            ! shaded leaves' values from the last timestep (see
                            ! sf_stom_state_mod).
                            !
                            ! MDK, 30/10/26
USE sf_stom_state_mod, ONLY: store_e_prev
                            ! The capacitance store coupled to the two-leaf
                            ! search: the last timestep's transpiration (see
                            ! sf_stom_state_mod).
                            !
                            ! MDK, 01/11/26
USE leaf_energy_balance_mod, ONLY: leaf_temp_setup,                         &
                                   leaf_temp_update, leaf_temp_update_n
                            ! Leaf energy balance.
                            !
                            ! MDK, 24/10/26

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SF_STOM_MOD'

INTEGER, SAVE, PRIVATE :: n_lt_clamp_warn = 0
INTEGER, PARAMETER, PRIVATE :: max_lt_clamp_warn = 20
    ! l_leaf_temp: warnings issued so far, and the most issued, for leaf
    ! temperatures held at the dt_leaf limit (see t_leaf_clamp).
    !
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
        min_gl_pft, lai_ref_k
                            ! lai_ref_k: reference LAI of the root and
                            ! stem conductance (kplant_model = 1).
                            !
                            ! MDK, 21/10/26
USE jules_vegetation_mod, ONLY:                                                &
! imported model ids.
    leaf_flux_fsmc, leaf_flux_stom_opt,                                        &
! imported parameters
    photo_collatz, photo_farquhar, photo_sox_collatz, photo_johnson,           &
    stomata_medlyn, stomata_sox_analytical,                                    &
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
    l_som_cuticular_loss, l_som_gravity, l_red,                                &
    l_som_cuticular_day_only, l_som_shared_stem, l_som_shared_gain,            &
    l_som_coupled_e, l_som_capacitance_on,                                     &
    l_som_vcmax_psi, kplant_model, kplant_model_ground,                        &
    l_nsc_sink_feedback,                                                       &
    l_leaf_temp, leaf_temp_iter, l_leaf_temp_gc_eq,                            &
    leaf_wind_ext, l_leaf_temp_ci

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


USE stom_opt_jls_mod, ONLY: stom_opt_mod, som_lambda_eff_pft,                  &
                            set_leaf_temp_ci, clear_leaf_temp_ci
                            ! MDK, 04/11/26
USE plant_water_store_mod, ONLY: store_w, store_q, store_c,                    &
                                 store_state_alloc, store_unset,               &
                                 store_w_shed, leak_e, store_q_gap,            &
                                 store_w_xyl, store_q_xyl
                            ! MDK, 05/11/26
USE pftparm, ONLY: cap_stem, kcrit_fractional_loss, conductance_b_seg,         &
                   conductance_c_seg
USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: l_som_plant_segments
USE pftparm, ONLY: pft_xylem_impairment_model
USE jules_vegetation_mod, ONLY: xylem_impairment_none

USE xylem_impairment_mod, ONLY: update_xylem_impairment
! MDK, 17/10/26
USE jules_vegetation_mod, ONLY: l_triffid, l_triffid_nsc, l_phen_psi
USE triffid_nsc_mod, ONLY: triffid_nsc_daily, triffid_nsc_turgor
USE jules_vegetation_mod, ONLY: l_phen_fuse, l_phen_heat
USE leaf_shed_mod, ONLY: phen_psi_update, phen_fuse_update,                   &
                         heat_damage_update
! MDK, 09/11/26
USE jules_vegetation_mod, ONLY: discrimination_model, discrimination_off
USE leaf_c13_discrim_mod, ONLY: d13c_disc_pft, d13c_disc_alloc, d13c_disc_leaf
USE canopy_diag_mod, ONLY: gs_sun_pft, gs_shd_pft, lai_sun_pft
USE canopy_diag_mod, ONLY: psi_leaf_sun_pft, psi_leaf_shd_pft
                            ! MDK, 30/10/26
USE canopy_diag_mod, ONLY: vcmax_pft, jmax_pft, anetc_pft, anetl_top_pft,      &
                           anetl_bot_pft
USE jules_vegetation_mod, ONLY: l_photo_diag
                            ! MDK, 30/10/26
USE canopy_diag_mod, ONLY: canopy_diag_alloc, gc_stomatal_pft, t_can_air_pft,  &
                           q_can_air_pft, h_leaf_pft, le_leaf_pft,            &
                           h_leaf_pts, le_leaf_pts
                            ! MDK, 24/10/26


USE plant_hydraulic_path_mod, ONLY: k_path_zero_flow, kcrit_of_path
USE plant_hydraulic_path_mod, ONLY: set_ksr_path, set_k_ground_basis,          &
                                    kplant_scale
                                       ! MDK, 04/11/26
USE plant_hydraulic_path_mod, ONLY: l_kscale
USE jules_vegetation_mod, ONLY: l_nsc_k_sapwood
                            ! MDK, 31/10/26
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_n
USE triffid_nsc_mod, ONLY: nsc_sapwood
USE jules_vegetation_mod, ONLY: l_ximpair_cav_water
USE xylem_impairment_memory_mod, ONLY: ximpair_lock
USE pftparm, ONLY: wood_density, xyl_lumen_frac
                            ! MDK, 12/11/26
                            ! MDK, 04/11/26
USE triffid_nsc_mod, ONLY: triffid_nsc_alloc, nsc_k_area, nsc_k_area_frac
                           ! MDK, 31/10/26
USE jules_vegetation_mod, ONLY: nsc_damage_response, nsc_damage_trugman,       &
                                k_ground_ratio_max
USE triffid_nsc_mod, ONLY: nsc_k_area_raw, nsc_lai_sapwood, nsc_f_a_min
                           ! MDK, 04/11/26
USE plant_hydraulic_path_mod, ONLY: set_store_q_path, clear_store_q_path,      &
                                    psi_leaf_segments_jls
                                       ! MDK, 24/10/26
USE plant_hydraulic_path_mod, ONLY: psi_leaf_segments_lut_jls
                            ! MDK, 31/10/26
USE plant_hydraulic_path_mod, ONLY: lut_seg_step, seg_kmax_fac
USE plant_hydraulic_path_mod, ONLY: set_grav_drop, clear_grav_drop,            &
                                    grav_node_drop
                            ! MDK, 04/11/26
USE supply_table_mod, ONLY: build_supply_lut_seg
USE xylem_impairment_seg_state_mod, ONLY: ximpair_seg_frac
                                       ! MDK, 04/11/26
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut
                            ! MDK, 31/10/26
USE plant_hydraulic_path_mod, ONLY: set_shared_stem, clear_shared_stem,       &
                                    shared_amax_self
                                       ! MDK, 30/10/26

USE planet_constants_mod, ONLY: g, r
USE water_constants_mod, ONLY: rho_water, Water_kg_per_mol
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
                            ! OUT Canopy transpiration from the stomatal
                            !     optimisation (mol H2O m-2 ground s-1). As
                            !     a work array inside the multilayer layer
                            !     loop it is summed over the layers.
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
                            ! OUT Whole-path (soil-to-leaf) conductance,
                            !     Kplant, per unit leaf area
                            !     (mol H2O m-2 leaf s-1 Pa-1).
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
                            ! held at the dt_leaf limit from the air
                            ! this timestep (0 without l_leaf_temp).
                            !
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
                            ! INOUT Canopy conductance to H2O (m/s).
                            ! The input (last timestep's) value is used with
                            ! can_rad_mod = 1 and by the two-leaf canopy
                            ! (first pass, and the state it falls back on).
,f_nsc(land_pts)                                                               &
                            ! INOUT Non-structural carbohydrate mass fraction
                            !      (kgC/kgC)
,k_max(land_pts)                                                               &
                            ! INOUT Maximum plant conductance with the xylem
                            ! impairment (mol H2O m-2 leaf s-1 Pa-1, per
                            ! unit leaf area).
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
                            ! Canopy conductance to H2O (m/s) of the stomatal
                            ! optimisation: gc after the canopy switch, as
                            ! the cuticular leak and the store set it.

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
                            ! WORK Leaf conductance for H2O (m/s), per unit
                            !      leaf area: of the top leaf in the big
                            !      leaf canopy without the stomatal
                            !      optimisation (gc = fpar gl), and of one
                            !      layer in the multilayer canopies. The
                            !      big leaf optimisation returns the canopy
                            !      value (per unit ground area) in gc, and
                            !      the two-leaf canopy uses gl_sun and
                            !      gl_shd.
,gl_sun(land_pts)                                                              &
                            ! WORK Stomatal conductance for H2O of the
                            !      sunlit leaves (m/s): per unit leaf area
                            !      of a layer in the multilayer canopies;
                            !      the class total per unit ground area in
                            !      the two-leaf canopy.
,gl_shd(land_pts)                                                              &
                            ! WORK As gl_sun, for the shaded leaves.
,el_sun(land_pts)                                                              &
                            ! WORK Transpiration of the sunlit leaves
                            !      (mol H2O m-2 s-1): per unit leaf area of
                            !      a layer in the multilayer canopy; the
                            !      class total per unit ground area in the
                            !      two-leaf canopy.
,el_shd(land_pts)                                                              &
                            ! WORK As el_sun, for the shaded leaves.
,leaf_k_sun(land_pts)                                                          &
                            ! WORK Whole-path conductance of the sunlit
                            !      leaves (mol H2O m-2 s-1 Pa-1, on the
                            !      basis of gl_sun).
,leaf_k_shd(land_pts)                                                          &
                            ! WORK As leaf_k_sun, for the shaded leaves.
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
                            !      impairment driver 4), else psi_zf.
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

REAL(KIND=real_jlslsm) :: je_dummy(land_pts), fapar_dummy(land_pts)
                            ! Placeholders for the big-leaf optimiser call
                            ! (stom_opt_mod uses them only for the
                            ! multilayer canopy).
                            !
                            ! MDK, 29/04/26
REAL(KIND=real_jlslsm) :: lam_sun(land_pts), lam_shd(land_pts)
                            ! Two-leaf: dA/dE at the optimum of the sunlit
                            ! and shaded leaves (lambda_eff of stom_opt_mod),
                            ! for the diagnostic som_lambda_eff_pft.
                            !
                            ! MDK, 23/10/26
REAL(KIND=real_jlslsm) :: kmax_leaf_mean(land_pts)
                            ! Multilayer: kmax_pft, the canopy mean per unit
                            ! leaf area of the soil-to-leaf conductance that
                            ! kmax_leaf_lyr spreads over the layers.
REAL(KIND=real_jlslsm) :: kcrit_leaf_mean(land_pts)
                            ! Its critical conductance (kcrit_of_path), spread
                            ! over the layers as kcrit_leaf_lyr.
REAL(KIND=real_jlslsm) :: kmax_leaf_lyr(land_pts)
                            ! Multilayer: this layer's soil-to-leaf
                            ! conductance per unit leaf area, kmax_leaf_mean
                            ! redistributed by the leaf-N profile (canopy mean
                            ! kmax_pft), each layer on its own parallel path
                            ! from psi_src.
REAL(KIND=real_jlslsm) :: kcrit_leaf_lyr(land_pts)
                            ! The same for kcrit_leaf_mean.
REAL(KIND=real_jlslsm) :: kleaf_lyr_fac
                            ! This layer's value of the conductance profile
                            ! through the canopy, exp(-(n - 1) kleaf_decay).
REAL(KIND=real_jlslsm) :: kleaf_canopy_mean
                            ! The profile's mean over the layers.
REAL(KIND=real_jlslsm) :: kleaf_decay
                            ! The profile's decay exponent per layer.
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
                            ! Soil water factor applied to the leaf fluxes
                            ! and passed to leaf_limits (fsmc_in).
                            !
                            ! MDK, 04/11/26
REAL(KIND=real_jlslsm) :: f_vcmax_down(land_pts)
                            ! Factor on Vcmax and Jmax (0-1): the soil-water
                            ! factor (l_som_vcmax_psi) times the sink
                            ! feedback (l_nsc_sink_feedback).
REAL(KIND=real_jlslsm) :: f_sink_feedback(land_pts)
                            ! Sink feedback factor (l_nsc_sink_feedback;
                            ! 1 = none).
                            !
                            ! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: vcmax_top_d(land_pts), jmax_top_d(land_pts),         &
                          rd_top_d(land_pts)
                            ! Vcmax, Jmax and Rd of a leaf at the top of
                            ! the canopy at tstar (vcmax_pft, jmax_pft).
                            !
                            ! MDK, 30/10/26
REAL(KIND=real_jlslsm) :: psi_src(land_pts)
                            ! Water potential at the base of the plant path
                            ! seen by the stomatal optimisation (Pa):
                            ! psi_root_zone, less the gravitational drop
                            ! rho_water g ht with l_som_gravity on a plant
                            ! without segments.
REAL(KIND=real_jlslsm) :: psi_zf(land_pts)
                            ! Water potential of the stem and leaf at zero
                            ! flow (Pa): psi_root_zone, less rho_water g ht
                            ! with l_som_gravity.
REAL(KIND=real_jlslsm) :: drop_grav(land_pts)
                            ! The gravitational drop rho_water g ht (Pa).
REAL(KIND=real_jlslsm) :: gl_leaf_cap
                            ! som_gl_max, or 0 (no cap) for the supply-loss
                            ! stomata (stomata_model = 9), which have none.
REAL(KIND=real_jlslsm) :: gl_max_lf(land_pts), gl_max_bigleaf(land_pts)
                            ! som_gl_max on the basis each stom_opt_mod call
                            ! works on: per leaf area for the multilayer
                            ! calls, canopy (x LAI) for big-leaf.
REAL(KIND=real_jlslsm) :: ci_sun_ml(land_pts)
                            ! Multilayer stomatal optimisation: ci of the
                            ! layer's sunlit leaves (Pa).
REAL(KIND=real_jlslsm) :: ci_shd_ml(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ci_sum_gl_ml(land_pts)
                            ! Canopy sum of ci weighted by gl dlai, giving the
                            ! gl-weighted canopy ci (as two-leaf).
REAL(KIND=real_jlslsm) :: ci_sum_lai_ml(land_pts)
                            ! Canopy sum of ci weighted by dlai, giving the
                            ! LAI mean when all the leaves are closed.
REAL(KIND=real_jlslsm) :: vcmax_sun_ml(land_pts)
                            ! Multilayer with l_leaf_temp: the layer's Vcmax
                            ! of the sunlit leaves at their temperature (per
                            ! unit leaf area).
REAL(KIND=real_jlslsm) :: vcmax_shd_ml(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: je_sun_ml(land_pts)
                            ! The layer's electron flux of the sunlit leaves
                            ! at their temperature (per unit leaf area).
REAL(KIND=real_jlslsm) :: je_shd_ml(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: rd_sun_ml(land_pts)
                            ! The layer's dark respiration of the sunlit
                            ! leaves at their temperature (per unit leaf
                            ! area).
REAL(KIND=real_jlslsm) :: rd_shd_ml(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ones_pts(land_pts)
                            ! Ones.
REAL(KIND=real_jlslsm) :: r_bl_sun_ml(land_pts)
                            ! Multilayer with l_leaf_temp: the leaf boundary
                            ! layer resistance to water vapour of the layer's
                            ! sunlit leaves (s m-1, leaf area), from the last
                            ! leaf energy balance (0 before the first). The
                            ! optimiser's transpiration then goes through the
                            ! same boundary layer as the latent heat of the
                            ! balance.
REAL(KIND=real_jlslsm) :: r_bl_shd_ml(land_pts)
                            ! The same for the shaded leaves.
                            !
                            ! MDK, 12/11/26
INTEGER :: n_pass_ml
                            ! Number of passes of the multilayer
                            ! optimisation: iter, or leaf_temp_iter with
                            ! l_leaf_temp.
REAL(KIND=real_jlslsm) :: n_clamp_lt(land_pts)
                            ! l_leaf_temp: number of leaf classes whose
                            ! temperature was held at dt_leaf from the
                            ! air in the last leaf temperature
                            ! update of the timestep.
                            !
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
INTEGER :: n_lit
                            ! Number of lit points in open_land_index (for
                            ! the impairment update and the leak), kept apart
                            ! from open_pts, which goes with open_index.
                            !
                            ! MDK, 25/10/26
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
REAL(KIND=real_jlslsm) :: gc_cuticular(land_pts)
                            ! Cuticular loss (l_som_cuticular_loss): the
                            ! canopy's cuticular conductance (m s-1) on the
                            ! deficit dqc, before the soil-supply cap (so
                            ! gc_cuticular dqc is the leak's flux).
REAL(KIND=real_jlslsm) :: gc_cuticular_capped(land_pts)
                            ! The same after the soil-supply cap (m s-1).
REAL(KIND=real_jlslsm) :: e_per_gc_leak
                            ! The leak's tile flux per unit conductance,
                            ! rho (qs(tstar) - q1) (kg m-3).
REAL(KIND=real_jlslsm) :: e_target_leak
                            ! The leak's target flux rho gc_cuticular_capped dqc
                            ! (kg m-2 s-1).
REAL(KIND=real_jlslsm) :: gc_leak
                            ! The canopy conductance that gives the target
                            ! flux through ra (m s-1).
                            !
                            ! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: lai_ref_kg(land_pts)
                            ! kplant_model = 1: reference LAI of the root
                            ! and stem conductances (lai_ref_k, or the
                            ! balanced LAI from canopy height).
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm) :: k_area(land_pts)
                            ! l_nsc_k_sapwood: the sapwood-area factor on
                            ! the plant's kmax (0-1).
                            !
                            ! MDK, 31/10/26
REAL(KIND=real_jlslsm) :: k_raw(land_pts), r_stem(land_pts), lai_sw
                            ! l_nsc_k_sapwood: the sapwood-area ratio before
                            ! the MIN, the stem's factor from its sapwood
                            ! area (Trugman response; < 0 = not set)
                            ! and the leaf area the sapwood can supply.
                            !
                            ! MDK, 01/11/26
REAL(KIND=real_jlslsm) :: q_store(land_pts)
                            ! The capacitance store: the store's flow out
                            ! at the stem node the optimiser sees (mol H2O
                            ! m-2 ground s-1; > 0: it gives, < 0: it refills).
REAL(KIND=real_jlslsm) :: kmax_ref_store(land_pts)
                            ! The capacitance store: the whole canopy's intact
                            ! plant conductance kmax for the store's balance
                            ! (store_cpl_q; mol m-2 ground s-1 Pa-1), set by
                            ! each canopy scheme before its first store solve.
REAL(KIND=real_jlslsm) :: kmax_store(land_pts)
                            ! The same, with the xylem impairment (mol m-2
                            ! ground s-1 Pa-1).
REAL(KIND=real_jlslsm) :: kcrit_store(land_pts)
                            ! The same canopy's kcrit, the conductance at the
                            ! critical loss (mol m-2 ground s-1 Pa-1).
REAL(KIND=real_jlslsm), ALLOCATABLE :: el_sun_lyr(:,:)
                            ! The capacitance store, multilayer: each layer's
                            ! sunlit E of the last pass (mol H2O m-2 leaf s-1;
                            ! land_pts, ilayers).
REAL(KIND=real_jlslsm), ALLOCATABLE :: el_shd_lyr(:,:)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: e_ml_prev(land_pts)
                            ! Multilayer: the canopy E of the last pass (mol
                            ! H2O m-2 ground s-1).
REAL(KIND=real_jlslsm) :: q_store_sun_lyr(land_pts)
                            ! This layer's sunlit part of q_store (mol H2O m-2
                            ! leaf s-1),
                            !
                            !   q_c = q_store E_c / E,
                            !
                            ! from the last pass (by the layer's kmax share,
                            ! share_sup, in the first pass or where E = 0), so
                            ! the parts over the canopy add up to q_store.
REAL(KIND=real_jlslsm) :: q_store_shd_lyr(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: q_store_sun(land_pts)
                            ! The capacitance store, two-leaf canopy: the
                            ! sunlit class's part of q_store (mol H2O m-2
                            ! ground s-1), in proportion to its flow,
                            !
                            !   q_c = q_store E_c / E,
                            !
                            ! from the previous pass's E_c (by the kmax_ref
                            ! shares in the first pass, or where E = 0), so the
                            ! class's root flow E_c (1 - q_store / E) is not
                            ! negative where q_store <= E (no class is searched
                            ! in back-flow, where water is cheap), and with
                            ! q_store = 0 the classes are the independent paths
                            ! of capacitance off. The parts add up to q_store.
REAL(KIND=real_jlslsm) :: q_store_shd(land_pts)
                            ! The same for the shaded class.
REAL(KIND=real_jlslsm) :: share_one(land_pts)
                            ! A share of 1, for set_store_q_path.
                            !
                            ! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: e_oth_2l(land_pts)
                            ! The shared root and stem (l_som_shared_stem):
                            ! the other class's transpiration (mol H2O m-2
                            ! ground s-1) that a class's search puts on the
                            ! shared root and stem: for the sunlit leaves the
                            ! shaded leaves' of the last pass (of the last
                            ! timestep in the first pass), for the shaded
                            ! leaves the sunlit leaves' of this pass.
                            !
                            ! MDK, 30/10/26
REAL(KIND=real_jlslsm) :: a_oth_2l(land_pts), amax_sun_2l(land_pts),          &
                          amax_shd_2l(land_pts)
                            ! The shared root and stem: the other class's A
                            ! at its hydraulic edge (the gain normalisation,
                            ! mol CO2 m-2 ground s-1), taken as e_oth_2l;
                            ! and each class's, from its last search.
                            !
                            ! MDK, 30/10/26
REAL(KIND=real_jlslsm) :: e_sup_use(land_pts)
                            ! The soil supply the optimiser is capped by, in
                            ! every canopy (kg m-2 s-1; < 0: no cap):
                            ! e_supply, or with the capacitance store
                            ! e_supply + q_store, as the soil supplies only
                            ! E - q_store.
REAL(KIND=real_jlslsm) :: e_store_bal(land_pts)
                            ! The capacitance store: the canopy E the
                            ! store's balance is taken at (mol H2O m-2 s-1).
REAL(KIND=real_jlslsm) :: q_store_bal(land_pts)
                            ! The capacitance store's balance over the step at
                            ! the canopy's E (store_cpl_q): the store's flow
                            ! out (mol H2O m-2 s-1).
REAL(KIND=real_jlslsm) :: psi_stem_bal(land_pts)
                            ! The stem node's psi at the end of the step (Pa).
REAL(KIND=real_jlslsm) :: psi_leaf_bal(land_pts)
                            ! The leaf psi at the end of the step (Pa).
                            !
                            ! MDK, 12/11/26
LOGICAL :: l_leaf_floor_hit(land_pts)
                            ! The capacitance store with the leak: the
                            ! committed leaf psi was below psi_floor_leaf before
                            ! its floor (store_cpl_q).
REAL(KIND=real_jlslsm) :: kplant_bal(land_pts)
                            ! The capacitance store: the canopy's plant
                            ! conductance at the committed E on the path
                            ! without the store's shift (the output
                            ! kplant), mol m-2 s-1 Pa-1, as kmax_sun_2l +
                            ! kmax_shd_2l.
                            !
                            ! MDK, 25/10/26
LOGICAL :: l_cap_store
                            ! The stomatal optimisation's plant capacitance
                            ! is on (l_som_capacitance_on: l_som_capacitance
                            ! or l_som_cuticular_loss): the optimiser sets
                            ! the transpiration, seeing the store through
                            ! its flow q_store, the leak (if on) drains it, and
                            ! the xylem impairment update is made after the
                            ! leak block, with the store's water potentials.
                            !
                            ! MDK, 04/11/26
REAL(KIND=real_jlslsm) :: psi_leaf_cls(land_pts,2), w_cls(land_pts,2)
                            ! Two-leaf canopy: the leaf water potential of
                            ! the sunlit (1) and shaded (2) leaves (Pa) and
                            ! their shares of the leaf area, for the leaf
                            ! segment's damage (update_xylem_impairment).
                            !
                            ! MDK, 04/11/26
REAL(KIND=real_jlslsm) :: w_xyl_full
                            ! l_ximpair_cav_water: the water of the stem's
                            ! conduits when all are full, V_w (kg m-2).
REAL(KIND=real_jlslsm) :: w_xyl_hold
                            ! The water the conduits should hold now, those
                            ! not lost, V_w (1 - L) (kg m-2).
REAL(KIND=real_jlslsm) :: w_xyl_excess
                            ! The conduits' water above w_xyl_hold (kg m-2):
                            ! > 0 is released to the store (embolism), < 0
                            ! is drawn from the soil (new xylem).
REAL(KIND=real_jlslsm) :: w_to_store
                            ! The part of w_xyl_excess the store takes this
                            ! step, at most up to full (kg m-2).
REAL(KIND=real_jlslsm), PARAMETER :: c_w_cav = 0.46
                            ! Carbon fraction of wood dry mass (Quercus mean,
                            ! Doraisami et al. 2022).
REAL(KIND=real_jlslsm), PARAMETER :: rho_cw_cav = 1500.0
                            ! Density of the wood substance in water, for the
                            ! wood's lumen volume when xyl_lumen_frac is
                            ! unset (kg m-3; Kellogg & Wangaard 1969).
                            !
                            ! MDK, 12/11/26
LOGICAL :: l_in_impair(land_pts)
                            ! Points already in the impairment update (the
                            ! lit points), so the leak's other points are
                            ! added once.
                            !
                            ! MDK, 26/10/26
LOGICAL :: l_leak_lit(land_pts)
                            ! Points where the daytime-only leak may apply
                            ! (l_som_cuticular_day_only): lit points (those
                            ! of the impairment update) with apar >=
                            ! apar_leak_min.
REAL(KIND=real_jlslsm), PARAMETER :: apar_leak_min = 2.0
                            ! Least absorbed PAR (W m-2) for the
                            ! daytime-only leak, about 5 W m-2 of
                            ! shortwave: not apar > 0, as forcing noise at
                            ! night (SW up to a few W m-2) would turn the
                            ! leak on in the dark.
                            !
                            ! MDK, 29/10/26
REAL(KIND=real_jlslsm), PARAMETER :: dq_ratio_leak_max = 3.0
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
REAL(KIND=real_jlslsm) :: psi_floor_leaf
                            ! Floor on psi_leaf: psi at the critical loss on
                            ! the leaf segment (or PFT) curve (Pa).
REAL(KIND=real_jlslsm) :: psi_floor_stem
                            ! Floor on psi_stem, as psi_floor_leaf on the stem
                            ! segment curve (Pa).
REAL(KIND=real_jlslsm) :: c_store
                            ! The store's capacitance C (mol H2O m-2 ground
                            ! Pa-1).
REAL(KIND=real_jlslsm) :: psi0_store
                            ! The store's (stem node) psi at the start of the
                            ! step, W / C (Pa).
LOGICAL :: l_ml_ximp
                            ! Multilayer stomatal optimisation with xylem
                            ! impairment (not in spin-up).
REAL(KIND=real_jlslsm) :: kfrac_imp_ml(land_pts)
                            ! Multilayer: the impaired fraction of the
                            ! conductance, k_max / kmax_pft (1 without
                            ! impairment).
REAL(KIND=real_jlslsm) :: kmax_imp_lyr(land_pts)
                            ! This layer's impaired path conductance per unit
                            ! leaf area (mol H2O m-2 leaf s-1 Pa-1),
                            !
                            !   kmax_imp_lyr = kmax_leaf_lyr kfrac_imp_ml
                            !
                            ! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: qs_air_leak(land_pts)
                            ! Saturated humidity at the air temperature
                            ! (kg kg-1).
REAL(KIND=real_jlslsm) :: dq_leak(land_pts)
                            ! The leaf-to-air deficit driving the leak
                            ! (kg kg-1).
REAL(KIND=real_jlslsm) :: w_store_new
                            ! The store's water after a fall of C (mol H2O
                            ! m-2 ground).
LOGICAL :: l_leak_on(land_pts)
                            ! Points where the leak binds.
REAL(KIND=real_jlslsm), PARAMETER :: lai_store_min = 0.01
                            ! Least LAI (m2 m-2) for the store: below it the
                            ! store is dropped, as psi = W / C has no
                            ! meaning as C goes to 0.
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm) :: dq_opt(land_pts)
                            ! Humidity deficit the stomatal optimisation
                            ! uses (kg kg-1), with and without coupled
                            ! transpiration (l_som_coupled_e).
REAL(KIND=real_jlslsm) :: r_bl_opt(land_pts)
                            ! Coupled transpiration: the leaf boundary layer
                            ! resistance passed to stom_opt_mod (s m-1).
REAL(KIND=real_jlslsm) :: r_ca_opt(land_pts)
                            ! Coupled transpiration: the canopy air
                            ! resistance passed to stom_opt_mod (s m-1).
REAL(KIND=real_jlslsm) :: g_other_opt(land_pts)
                            ! Coupled transpiration: the other leaf class's
                            ! effective conductance passed to stom_opt_mod
                            ! (m s-1).
REAL(KIND=real_jlslsm) :: zeros_pts(land_pts)
                            ! Zeros, for stom_opt_mod arguments not used.
REAL(KIND=real_jlslsm) :: gbw_sun_lt(land_pts)
                            ! Two-leaf, coupled: the sunlit leaves' boundary
                            ! layer conductance for water (m s-1, class
                            ! basis), from the last leaf energy balance.
REAL(KIND=real_jlslsm) :: gbw_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ra_ca_lt(land_pts)
                            ! Two-leaf, coupled: the canopy air resistance
                            ! (s m-1), from the last leaf energy balance.
REAL(KIND=real_jlslsm) :: g_soil_c(land_pts)
                            ! Soil evaporation conductance beneath the
                            ! canopy (m s-1), which shares the canopy air and
                            ! ra with the leaves.
                            !
                            ! MDK, 12/11/26

! Name suffixes in the canopy code below: _2l the two-leaf canopy
! (can_rad_mod = 7), _lt the leaf energy balance (l_leaf_temp), _ml and _lyr
! the multilayer canopy, _sun and _shd the sunlit and shaded leaves.

! l_leaf_temp_ci (two-leaf): each class's leaf temperature as a linear
! function of its transpiration from the last pass's balance,
!
!   T = tleaf_a_2l - tleaf_b_2l E
!
! (tleaf_b_2l = 0 before the first balance of the step), and the
! temperature slopes of the deficit and the photosynthesis parameters
! about it.
!
! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: tleaf_a_2l(land_pts,2)
                            ! T_leaf at E = 0 of the sunlit (1) and shaded
                            ! (2) leaves (K).
REAL(KIND=real_jlslsm) :: tleaf_b_2l(land_pts,2)
                            ! The fall of T_leaf per unit E of each class
                            ! (K per mol m-2 s-1).
REAL(KIND=real_jlslsm) :: dqs_dt_2l(land_pts)
                            ! Change of the saturation humidity per kelvin,
                            ! dqs / dT (kg kg-1 K-1).
REAL(KIND=real_jlslsm) :: dlnpar_dt_2l(land_pts,7)
                            ! Change per kelvin of the log of each
                            ! photosynthesis parameter (Vcmax, Je, Kc, Ko,
                            ! Gamma*, Km, Rd).

! Two-leaf (sunlit/shaded big-leaf) canopy, can_rad_mod = 7. Canopy (per m2
! ground) totals for each leaf class, aggregated from the layered radiation
! profile.
!
! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: lai_sun_2l(land_pts)
                            ! Sunlit leaf area (m2 leaf m-2 ground).
REAL(KIND=real_jlslsm) :: lai_shd_2l(land_pts)
                            ! Shaded leaf area (m2 leaf m-2 ground).
REAL(KIND=real_jlslsm) :: nweight_sun_2l(land_pts)
                            ! Sunlit leaf area weighted by the relative
                            ! leaf-N profile exp(-knl L), the sunlit N
                            ! weight of the canopy.
REAL(KIND=real_jlslsm) :: nweight_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: apar_sun_2l(land_pts)
                            ! Absorbed PAR by the sunlit leaves (W m-2
                            ! ground).
REAL(KIND=real_jlslsm) :: apar_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: apar_mol_sun_2l(land_pts)
                            ! As apar_sun_2l, in mol photons m-2 s-1.
REAL(KIND=real_jlslsm) :: apar_mol_shd_2l(land_pts)
                            ! As apar_shd_2l, in mol photons m-2 s-1.
REAL(KIND=real_jlslsm) :: vcmax_sun_2l(land_pts)
                            ! Canopy Vcmax of the sunlit leaves (mol m-2
                            ! ground s-1).
REAL(KIND=real_jlslsm) :: vcmax_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: jmax_sun_2l(land_pts)
                            ! Canopy Jmax of the sunlit leaves (mol m-2
                            ! ground s-1).
REAL(KIND=real_jlslsm) :: jmax_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: kmax_sun_2l(land_pts)
                            ! Hydraulic conductance of the sunlit leaves'
                            ! path (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: kmax_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: kcrit_sun_2l(land_pts)
                            ! kcrit of the sunlit leaves' path (mol m-2 s-1
                            ! Pa-1).
REAL(KIND=real_jlslsm) :: kcrit_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: gl_max_sun_2l(land_pts)
                            ! Maximum stomatal conductance of the sunlit
                            ! leaves (m s-1).
REAL(KIND=real_jlslsm) :: gl_max_shd_2l(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ci_sun_2l(land_pts)
                            ! Internal CO2 of the sunlit leaves (Pa).
REAL(KIND=real_jlslsm) :: ci_shd_2l(land_pts)
                            ! Internal CO2 of the shaded leaves (Pa).
REAL(KIND=real_jlslsm) :: f_sun_2l
                            ! Sunlit fraction of LAI.
REAL(KIND=real_jlslsm) :: n_cls_with_lai
                            ! Multilayer canopy: the number of leaf classes
                            ! with leaf area.
REAL(KIND=real_jlslsm) :: dnweight_2l
                            ! Layer integral of exp(-kn_2l L) over dlai.
REAL(KIND=real_jlslsm) :: kn_2l
                            ! N extinction per unit LAI of the two-leaf
                            ! classes (= knl).
REAL(KIND=real_jlslsm) :: nweight_tot_2l
                            ! Canopy N weight N_sun + N_shd (floored at
                            ! TINY); each class gets its share N_class /
                            ! nweight_tot_2l of the canopy kmax.

! Two-leaf leaf temperatures (l_leaf_temp).
!
! MDK, 12/11/26
REAL(KIND=real_jlslsm) :: t_sun_lt(land_pts)
                            ! Sunlit leaf temperature of the current pass
                            ! (K).
REAL(KIND=real_jlslsm) :: t_shd_lt(land_pts)
                            ! Shaded leaf temperature of the current pass
                            ! (K).
REAL(KIND=real_jlslsm) :: t_c_lt(land_pts)
                            ! Temperature of the air the leaves are in, tair
                            ! (K).
REAL(KIND=real_jlslsm) :: q_c_lt(land_pts)
                            ! Specific humidity of that air, q1 (kg kg-1).
REAL(KIND=real_jlslsm) :: sw_sun_lt(land_pts)
                            ! Absorbed PAR + near-IR of the sunlit leaves
                            ! (W m-2 ground).
REAL(KIND=real_jlslsm) :: sw_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: lwwt_sky_sun_lt(land_pts)
                            ! Layer weight of the isothermal long-wave from
                            ! the sky, sunlit leaves.
REAL(KIND=real_jlslsm) :: lwwt_sky_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: lwwt_soil_sun_lt(land_pts)
                            ! Layer weight of the isothermal long-wave from
                            ! the soil, sunlit leaves.
REAL(KIND=real_jlslsm) :: lwwt_soil_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: grad_fac_sun_lt(land_pts)
                            ! Radiative conductance scalar of the sunlit
                            ! leaves.
REAL(KIND=real_jlslsm) :: grad_fac_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: lai_wind_sun_lt(land_pts)
                            ! Wind-weighted leaf area of the sunlit leaves,
                            ! for forced convection.
REAL(KIND=real_jlslsm) :: lai_wind_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ccp_sun_lt(land_pts)
                            ! Gamma* (CO2 compensation point) of the sunlit
                            ! leaves at their temperature.
REAL(KIND=real_jlslsm) :: ccp_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: kc_sun_lt(land_pts)
                            ! Kc of the sunlit leaves at their
                            ! temperature.
REAL(KIND=real_jlslsm) :: kc_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: ko_sun_lt(land_pts)
                            ! Ko of the sunlit leaves at their
                            ! temperature.
REAL(KIND=real_jlslsm) :: ko_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: km_sun_lt(land_pts)
                            ! Km of the sunlit leaves at their
                            ! temperature.
REAL(KIND=real_jlslsm) :: km_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: qs_sun_lt(land_pts)
                            ! Saturated humidity at the sunlit leaves
                            ! (kg kg-1).
REAL(KIND=real_jlslsm) :: qs_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: dq_sun_lt(land_pts)
                            ! Humidity deficit from the sunlit leaves to the
                            ! canopy air (kg kg-1).
REAL(KIND=real_jlslsm) :: dq_shd_lt(land_pts)
                            ! The same for the shaded leaves.
REAL(KIND=real_jlslsm) :: rdc_tstar(land_pts)
                            ! Canopy dark respiration at tstar, for the stem
                            ! and root respiration.
REAL(KIND=real_jlslsm) :: wind_ext(land_pts)
                            ! Extinction coefficient of the wind with leaf
                            ! area in the canopy.
REAL(KIND=real_jlslsm) :: gc_eq_cap
                            ! Cap on the transpiration-equivalent gc
                            ! (l_leaf_temp_gc_eq).
REAL(KIND=real_jlslsm) :: r_eq_gc
                            ! Resistance for the transpiration-equivalent gc
                            ! (l_leaf_temp_gc_eq).
REAL(KIND=real_jlslsm), PARAMETER :: m_h2o_lt = Water_kg_per_mol
                            ! Molar mass of water (kg mol-1).
INTEGER :: n_pass_2l
                            ! Number of passes of the two-leaf optimisation.
LOGICAL :: l_leaf_temp_lyr
                            ! l_leaf_temp in the multilayer canopy: a leaf
                            ! energy balance for the sunlit and shaded leaf
                            ! of every layer.
INTEGER :: n_cls_lt
                            ! With l_leaf_temp_lyr, the number of leaf
                            ! classes, 2 x ilayers (class 2n-1 the sunlit and
                            ! 2n the shaded leaf of layer n).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lai_lt_c(:,:)
                            ! With l_leaf_temp_lyr (land_pts, n_cls_lt): the
                            ! leaf area of each class.
REAL(KIND=real_jlslsm), ALLOCATABLE :: sw_lt_c(:,:)
                            ! Absorbed short-wave of each class (see
                            ! leaf_temp_setup; fixed over the passes).
REAL(KIND=real_jlslsm), ALLOCATABLE :: lwwt_sky_lt_c(:,:)
                            ! Long-wave weight from the sky of each class.
REAL(KIND=real_jlslsm), ALLOCATABLE :: lwwt_soil_lt_c(:,:)
                            ! Long-wave weight from the soil of each class.
REAL(KIND=real_jlslsm), ALLOCATABLE :: grad_fac_lt_c(:,:)
                            ! Radiative conductance scalar of each class.
REAL(KIND=real_jlslsm), ALLOCATABLE :: lai_wind_lt_c(:,:)
                            ! Wind-weighted leaf area of each class.
REAL(KIND=real_jlslsm), ALLOCATABLE :: gl_lt_c(:,:)
                            ! Stomatal conductance of each class (m s-1,
                            ! ground).
REAL(KIND=real_jlslsm), ALLOCATABLE :: t_lt_c(:,:)
                            ! Leaf temperature of each class (K).
REAL(KIND=real_jlslsm), ALLOCATABLE :: gbw_lt_c(:,:)
                            ! Boundary-layer conductance for water of each
                            ! class.
                            !
                            ! MDK, 12/11/26
LOGICAL :: l_multilayer
                            ! .TRUE. for the multilayer canopy's calls to
                            ! stom_opt_mod (it then uses je_dummy and
                            ! fapar_dummy).
                            !
                            ! MDK, 24/10/26



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
zeros_pts(:) = 0.0
IF ( stomata_model == stomata_supply_loss ) gl_leaf_cap = 0.0
! (Initialised here so that every path, not only two-leaf, has them set.)
gbw_sun_lt(:) = 0.0
gbw_shd_lt(:) = 0.0
ra_ca_lt(:)   = 0.0

!-----------------------------------------------------------------------------
! Gravity (l_som_gravity). Lifting water to the top of the canopy costs
!
!   drop = rho_water g ht   (0.01 MPa per m)
!
! on top of the frictional drop (Xu et al. 2016, Eqn S2a). The whole
! canopy is put at the canopy height, the conservative end for the
! multilayer scheme. With segments (l_som_plant_segments) the drop is
! taken at the node between the root and the stem (set_grav_drop): the
! soil-to-root link and the root segment work at psi_root_zone, and the
! stem and leaf segments start from the root's leaf-side water potential
! less the drop, so their curves are evaluated at the pressure in the
! upper xylem. Without segments there is no stem node, so the drop is
! taken at the source (psi_src) and the whole-plant curve is evaluated on
! the shifted potentials. Either way the stem and leaf at zero flow are at
! psi_zf = psi_root_zone - drop.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
CALL clear_grav_drop()
psi_src(:) = psi_root_zone(:)
psi_zf(:) = psi_root_zone(:)
psi_stem(:) = psi_root_zone(:)
! Plant capacitance (the plant water store) of the stomatal optimisation
! (l_som_capacitance_on: l_som_capacitance or the cuticular loss) moves the
! xylem impairment update after it (see below).
!
! MDK, 25/10/26
l_cap_store = l_som_capacitance_on
psi_stem_sun(:) = psi_root_zone(:)
psi_stem_shd(:) = psi_root_zone(:)
IF ( l_som_gravity .AND. leaf_flux_mod == leaf_flux_stom_opt ) THEN
  drop_grav(:) = 0.0
  DO m = 1,veg_pts
    l = veg_index(m)
    drop_grav(l) = rho_water * g * MAX(ht(l), 0.0)
    psi_zf(l) = psi_root_zone(l) - drop_grav(l)
  END DO
  IF ( l_som_plant_segments ) THEN
    CALL set_grav_drop( land_pts, veg_pts, veg_index, drop_grav )
  ELSE
    psi_src(:) = psi_zf(:)
  END IF
  psi_stem(:) = psi_zf(:)
  psi_stem_sun(:) = psi_zf(:)
  psi_stem_shd(:) = psi_zf(:)
END IF

! Soil-water capacity factor of Vcmax/Jmax (l_som_vcmax_psi), once per
! timestep, so the recovery lag advances once whatever the number of leaf
! classes and leaf temperature iterations. The sink feedback
! (l_nsc_sink_feedback) multiplies it: the direct response to drought,
! then a further limit when carbon backs up in the store.
!
! MDK, 31/10/26
f_vcmax_down(:) = 1.0
IF ( l_som_vcmax_psi ) CALL vcmax_psi_factor( ft, land_pts, veg_pts,           &
                                              veg_index, psi_root_zone,        &
                                              f_vcmax_down )
IF ( l_nsc_sink_feedback ) THEN
  f_sink_feedback(:) = 1.0
  CALL nsc_sink_factor( ft, land_pts, veg_pts, veg_index, psi_root_zone,       &
                        f_sink_feedback )
  DO m = 1,veg_pts
    l = veg_index(m)
    f_vcmax_down(l) = f_vcmax_down(l) * f_sink_feedback(l)
  END DO
END IF

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
! Canopy diagnostics of this PFT (canopy_diag_mod).
!
! MDK, 30/10/26
CALL canopy_diag_alloc( land_pts )
gc_stomatal_pft(:,ft) = 0.0
n_lit = 0
gs_sun_pft(:,ft)      = 0.0
gs_shd_pft(:,ft)      = 0.0
lai_sun_pft(:,ft)     = 0.0
vcmax_pft(:,ft)       = 0.0
jmax_pft(:,ft)        = 0.0
anetc_pft(:,ft)       = 0.0
anetl_top_pft(:,ft)   = 0.0
anetl_bot_pft(:,ft)   = 0.0
t_can_air_pft(:,ft)   = 0.0
q_can_air_pft(:,ft)   = 0.0
h_leaf_pft(:,ft)      = 0.0
le_leaf_pft(:,ft)     = 0.0

SELECT CASE ( can_rad_mod )
CASE ( 4,5,6 )
!$OMP PARALLEL DO IF(land_pts > 1)                                             &
!$OMP DEFAULT(NONE)                                                            &
!$OMP SHARED(land_pts,anetc,gc,rdc,rdmean,psi_leaf,el,leaf_k,hydraulic_cost,   &
!$OMP        carbon_gain,kmax_leaf_mean,kcrit_leaf_mean)                       &
!$OMP PRIVATE(l)                                                               &
!$OMP SCHEDULE(STATIC)
  DO l = 1, land_pts
    anetc(l)  = 0.0
    gc(l)     = 0.0
    rdc(l)    = 0.0
    rdmean(l) = 0.0

    ! gs_opt multi-layer.
    !
    ! MDK, 09/06/26
    psi_leaf(l)        = 0.0
    el(l)              = 0.0
    leaf_k(l)          = 0.0
    hydraulic_cost(l)  = 0.0
    carbon_gain(l)     = 0.0
    kmax_leaf_mean(l)  = 0.0
    kcrit_leaf_mean(l) = 0.0

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
                            vcmax_temp, jmax_temp, i2 )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'pft_photo_model should be photo_collatz or photo_farquhar')

END SELECT  !  pft_photo_model

!-----------------------------------------------------------------------------
! Diagnostics (l_photo_diag, only when output): Vcmax and Jmax of a leaf at
! the top of the canopy at tstar, with any thermal acclimation and the
! soil-water factor (canopy_diag_mod).
!
! MDK, 30/10/26
!-----------------------------------------------------------------------------
IF ( l_photo_diag ) THEN
  CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,          &
                              veg_index, denom, jmax_temp, jv25, nleaf_top,    &
                              qtenf_term, vcmax_temp, jmax_top_d, rd_top_d,    &
                              vcmax_top_d, f_vcmax_down = f_vcmax_down )
  DO m = 1,veg_pts
    l = veg_index(m)
    vcmax_pft(l,ft) = vcmax_top_d(l)
    jmax_pft(l,ft)  = jmax_top_d(l)
  END DO
END IF

! Soil water factor for the leaf fluxes.
fsmc(:) = fsmc_in(:)

!-----------------------------------------------------------------------------
! Leaf temperatures output: tstar unless l_leaf_temp (set in can_rad_mod 5,
! 6 and 7).
!
! MDK, 19/10/26
!-----------------------------------------------------------------------------
n_clamp_lt(:)   = 0.0
ones_pts(:)     = 1.0
l_leaf_temp_lyr = l_leaf_temp                                                  &
                  .AND. ( can_rad_mod == 5 .OR. can_rad_mod == 6 )             &
                  .AND. leaf_flux_mod == leaf_flux_stom_opt
n_cls_lt = 2 * ilayers
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
!
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
! Conductance from the sapwood area (l_nsc_k_sapwood): the factor f_A on
! kmax, its intact reference and kcrit (the stem segment's with segments),
! 1 in spin-up. With the Trugman response to damage (nsc_damage_response
! = 1) the stem's factor is the leaf area its sapwood can supply over the
! LAI, L_A / LAI, held between nsc_f_a_min and k_ground_ratio_max (the
! ground basis's bound as LAI -> 0), so the stem's conductance per unit
! ground area follows its sapwood and shedding leaves raises the supply
! per unit leaf area; nsc_k_area_frac is then that factor.
!
! MDK, 04/11/26
IF ( l_nsc_k_sapwood ) THEN
  CALL triffid_nsc_alloc( land_pts )
  k_area(:) = 1.0
  k_raw(:)  = 1.0
  r_stem(:) = -1.0
  IF ( .NOT. is_spinup ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      k_area(l) = nsc_k_area(l, ft, k_raw(l))
      IF ( nsc_damage_response == nsc_damage_trugman ) THEN
        lai_sw = nsc_lai_sapwood(l, ft, ht(l))
        IF ( lai_sw >= 0.0 )                                                   &
          r_stem(l) = MIN(MAX(lai_sw / MAX(lai(l), TINY(1.0_real_jlslsm)),     &
                              nsc_f_a_min), k_ground_ratio_max)
      END IF
    END DO
  END IF
  nsc_k_area_frac(:,ft) = k_area(:)
  nsc_k_area_raw(:,ft)  = k_raw(:)
  IF ( nsc_damage_response == nsc_damage_trugman ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( r_stem(l) > 0.0 ) nsc_k_area_frac(l,ft) = r_stem(l)
    END DO
    CALL set_k_ground_basis( ft, land_pts, veg_pts, veg_index, lai,           &
                             lai_ref_kg, k_area, r_stem )
  ELSE
    CALL set_k_ground_basis( ft, land_pts, veg_pts, veg_index, lai,           &
                             lai_ref_kg, k_area )
  END IF
ELSE
  CALL set_k_ground_basis( ft, land_pts, veg_pts, veg_index, lai, lai_ref_kg )
END IF

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
!                        tissue takes the rest of the deficit (store_w_shed),
!                        so shedding does not make psi more negative;
!   C rises (growth):    W is kept: the new tissue is born full.
!
! Set before the canopy fluxes, which depend on the store.
!
! MDK, 23/10/26
IF ( l_cap_store ) THEN
  CALL store_state_alloc( land_pts )
  store_w_shed(:,ft) = 0.0
  DO m = 1,veg_pts
    l = veg_index(m)
    ! This step's C; last step's is in store_c (0 if not yet set).
    c_store = MAX(cap_stem(ft) * lai_ref_kg(l), TINY(1.0_real_jlslsm))
    ! C has fallen and the store is set: keep psi, shed the difference.
    IF ( store_w(l,ft) <= 0.0 .AND. store_c(l,ft) > c_store ) THEN
      w_store_new = store_w(l,ft) * c_store / store_c(l,ft)
      store_w_shed(l,ft)  = (w_store_new - store_w(l,ft)) * Water_kg_per_mol
      store_w(l,ft) = w_store_new
    END IF
    store_c(l,ft) = c_store
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
                                nleaf_layer, qtenf_term, vcmax_temp, jmax,     &
                                rd_dark, vcmax, f_vcmax_down = f_vcmax_down )

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
,                       acr, apar, ca, ccp, dqc, fsmc, je, kc, km, ko, oa      &
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
    ! Net leaf photosynthesis of the top and the bottom layer
    ! (l_photo_diag, canopy_diag_mod).
    !
    ! MDK, 30/10/26
    IF ( l_photo_diag ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)
        IF ( n == 1 ) anetl_top_pft(l,ft) = anetl(l)
        IF ( n == ilayers ) anetl_bot_pft(l,ft) = anetl(l)
      END DO
    END IF

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
  ! shared canopy node: leaves solved from a shared node psi would each be
  ! charged only their own leaf segment, leaving the root-to-canopy drop
  ! uncharged, so each leaf would see about half the marginal hydraulic cost
  ! in moist soil, and E and gs would be high relative to big-leaf and
  ! two-leaf.
  !
  ! MDK, 30/10/26
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
    l_ml_ximp = .NOT. is_spinup .AND.                                          &
                pft_xylem_impairment_model(ft) /= xylem_impairment_none
    kfrac_imp_ml(:) = 1.0
    IF ( l_ml_ximp ) THEN
      b_pft_pts(:) = conductance_b(:)
      c_pft_pts(:) = conductance_c(:)
      DO m = 1,veg_pts
        l = veg_index(m)
        kfrac_imp_ml(l) = k_max(l) / kmax_pft(ft)
      END DO
    ELSE
      b_pft_pts(:) = conductance_b_pft(ft)
      c_pft_pts(:) = conductance_c_pft(ft)
    END IF

    fsmc_unity(:) = 1.0
    gl_max_lf(:)  = gl_leaf_cap

        ! The humidity-deficit iteration (DO k = 1,iter) runs over the whole
        ! canopy, as in the big-leaf and two-leaf schemes: the deficit at the
        ! leaves is
        !
        !   dq / (1 + ra gc),
        !
        ! with gc the canopy conductance (per m2 ground) from the previous
        ! pass, the conductance the surface energy balance (sf_evap) sees.
        ! Iterating within each layer on that layer's per-leaf-area gl would
        ! give ra gl ~ 0, so every leaf would see nearly the full
        ! above-canopy deficit. Only the last pass is kept.
        gc_ml_prev(:) = 0.0

        ! l_leaf_temp: as in MAESPA, a leaf energy balance for the sunlit
        ! and the shaded leaf of every layer (leaf_temp_update_n with 2 x
        ! ilayers classes), each from its own absorbed short- and long-wave
        ! radiation, stomatal conductance and forced-convection boundary
        ! layer at the wind of its depth, all in the level-1 air. Each
        ! layer's classes use their own temperatures for their
        ! photosynthesis temperature terms, respiration and humidity
        ! deficit qs(T_leaf) - q1.
        !
        ! The temperatures are updated after each pass but the
        ! last (leaf_temp_iter passes), so there are no more stom_opt_mod
        ! calls than without l_leaf_temp for leaf_temp_iter = iter.
        n_pass_ml = iter
        IF ( l_leaf_temp ) THEN
          CALL sun_shade_totals( land_pts, veg_pts, veg_index, ilayers, dlai,  &
                                 fsun, fapar_sun, fapar_shd, ipar,             &
                                 lai_sun_2l, lai_shd_2l, apar_sun_2l,          &
                                 apar_shd_2l )
          ! The classes of every layer.
          ALLOCATE( lai_lt_c(land_pts,n_cls_lt), sw_lt_c(land_pts,n_cls_lt),   &
                    lwwt_sky_lt_c(land_pts,n_cls_lt),                          &
                    lwwt_soil_lt_c(land_pts,n_cls_lt),                         &
                    grad_fac_lt_c(land_pts,n_cls_lt),                          &
                    lai_wind_lt_c(land_pts,n_cls_lt),                          &
                    gl_lt_c(land_pts,n_cls_lt), t_lt_c(land_pts,n_cls_lt),     &
                    gbw_lt_c(land_pts,n_cls_lt) )
          gl_lt_c(:,:)  = 0.0
          gbw_lt_c(:,:) = 0.0
          CALL lt_init()
          n_pass_ml = leaf_temp_iter
          ! First guess of every layer: last timestep's sunlit / shaded
          ! canopy temperature (or tstar), as lt_init set.
          DO n = 1,ilayers
            t_lt_c(:,2*n-1) = t_sun_lt(:)
            t_lt_c(:,2*n)   = t_shd_lt(:)
          END DO
        END IF

        ! The plant capacitance store, multilayer: the store's flow out the
        ! first pass sees (last step's), the soil supply cap with it, and
        ! the canopy's path for its balance (as the big leaf; no
        ! kplant_model = 1 here).
        e_sup_use(:) = e_supply(:)
        IF ( l_cap_store ) THEN
          ! (It needs at least 2 passes: iter and leaf_temp_iter are both
          ! fixed at 3.)
          share_one(:) = 1.0
          e_ml_prev(:) = 0.0
          CALL store_cpl_first_q()
          DO m = 1,veg_pts
            l                 = veg_index(m)
            kmax_ref_store(l) = kmax_pft(ft) * lai(l)
            kmax_store(l)     = kmax_ref_store(l) * kfrac_imp_ml(l)
            kcrit_store(l)    = kcrit_of_path(ft, kmax_ref_store(l))
          END DO
          CALL store_cpl_supply()
          IF ( ALLOCATED(el_sun_lyr) ) DEALLOCATE( el_sun_lyr, el_shd_lyr )
          ALLOCATE( el_sun_lyr(land_pts,ilayers), el_shd_lyr(land_pts,ilayers) )
          el_sun_lyr(:,:) = 0.0
          el_shd_lyr(:,:) = 0.0
        END IF
        DO k = 1,n_pass_ml

        IF ( l_leaf_temp ) THEN
          rdc_tstar(:)  = 0.0
        END IF

        el(:) = 0.0
        anetl(:) = 0.0
        anetc(:) = 0.0
        gc(:) = 0.0
        rdc(:) = 0.0
        psi_leaf(:) = 0.0
        ! With plant segments, the leaf-area-weighted stem psi of the
        ! leaves' paths.
        IF ( l_som_plant_segments ) psi_stem(:) = 0.0
        leaf_k(:) = 0.0
        carbon_gain(:) = 0.0
        hydraulic_cost(:) = 0.0
        ci_sum_gl_ml(:) = 0.0
        ci_sum_lai_ml(:) = 0.0
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
        !$OMP             PRIVATE(l, m, kleaf_lyr_fac, kleaf_canopy_mean,      &
        !$OMP                     kleaf_decay)                                 &
        !$OMP             SHARED(ft, gl, ilayers, kn, knl, n, nleaf_top,       &
        !$OMP                    nleaf_layer,                                  &
        !$OMP   veg_index,dlai,can_rad_mod, veg_pts, kmax_pft,                 &
        !$OMP   kmax_leaf_mean, kcrit_leaf_mean, kmax_leaf_lyr,                &
        !$OMP   kcrit_leaf_lyr,                                                &
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

            ! Kmax is NOT scaled with leaf nitrogen: kmax_pft is the
            ! leaf-area-basis whole-plant conductance, the same in every layer.
            !
            ! Note: Kmax is intentionally NOT split by sunlit/shaded fractions
            ! (fsun), as this represents radiative partitioning rather than
            ! structural limitation.
            !
            ! kmax_leaf_mean is per unit leaf area, as vcmax and nleaf_layer
            ! are, so the canopy total (sum of kmax_leaf_mean * dlai) is
            ! kmax_pft * LAI.
            kmax_leaf_mean(l) = kmax_pft(ft)

            ! Scale kcrit by the same fraction of whole-plant kmax that this
            ! layer received, so the critical-conductance threshold stays
            ! calibrated at (1 - kcrit_fractional_loss) of THIS layer's own
            ! capacity (Brodribb & Cochard, 2009; Sabot et al., 2020).
            kcrit_leaf_mean(l) = kcrit_of_path(ft, kmax_leaf_mean(l))

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
            ! kleaf_canopy_mean is the closed-form mean of the same discrete
            ! profile as nleaf_layer over the ilayers layers (equal dlai).
            IF ( can_rad_mod == 6 ) THEN
                kleaf_decay = knl(ft) * dlai(l)
            ELSE
                kleaf_decay = kn(ft) / REAL(ilayers)
            END IF
            kleaf_lyr_fac = EXP(-REAL(n-1) * kleaf_decay)
            IF ( kleaf_decay > 1.0e-6 ) THEN
                kleaf_canopy_mean = (1.0 - EXP(-REAL(ilayers) * kleaf_decay))  &
                             / (REAL(ilayers) * (1.0 - EXP(-kleaf_decay)))
            ELSE
                kleaf_canopy_mean = 1.0
            END IF
            kmax_leaf_lyr(l) = kmax_leaf_mean(l) * (kleaf_lyr_fac /            &
                               kleaf_canopy_mean)
            kcrit_leaf_lyr(l) = kcrit_leaf_mean(l) * (kleaf_lyr_fac /          &
                                kleaf_canopy_mean)

            ! Soil supply cap per unit leaf area (el is per leaf area in
            ! each layer): each leaf may use its kmax share of the supply,
            ! as in the two-leaf scheme. Summed over the canopy (dlai
            ! weights) the shares add to 1. Sunlit and shaded leaves of a
            ! layer have the same kmax per leaf area, so the same share.
            share_sup(l) = (kleaf_lyr_fac / kleaf_canopy_mean) / MAX(lai(l),   &
                           EPSILON(1.0))

            END DO
        !$OMP END PARALLEL DO

            ! This layer's impaired path.
            DO m = 1,veg_pts
              l = veg_index(m)
              kmax_imp_lyr(l) = kmax_leaf_lyr(l) * kfrac_imp_ml(l)
            END DO

            !-------------------------------------------------------------------
            ! Calculate photosynthetic parameters.
            !-------------------------------------------------------------------
            CALL calc_photo_parameters( ft, land_pts, pft_photo_model,         &
                                        veg_pts, veg_index, denom, jmax_temp,  &
                                        jv25, nleaf_layer, qtenf_term,         &
                                        vcmax_temp, jmax, rd_dark, vcmax,      &
                                        f_vcmax_down = f_vcmax_down )

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
            ! whenever fsmc hits 0, regardless of what the hydraulics say.
            ! With fsmc_unity the Medlyn closure test reduces to apar == 0,
            ! the same as big-leaf (which classifies open/closed on apar
            ! alone; see the leaf_flux_stom_opt branch of the big-leaf path).
            ! (The Jacobs model's dq >= dqcrit closure still applies here.)
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
                !dKmax_sun(l,n) = kmax_leaf_mean(l) * fsun(l,n)
                !dKmax_shd(l,n) = kmax_leaf_mean(l) * (1.0 - fsun(l,n))

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
            ! below).

                l_multilayer = .TRUE.
                IF ( l_leaf_temp ) THEN
                  ! Each class at its own leaf temperature (see DO k):
                  ! Vcmax, Jmax (and so the electron flux), Rd, the
                  ! Michaelis-Menten terms and the compensation point of
                  ! this layer's leaf N, and the deficit qs(T_leaf) - q_c.
                  ! The electron flux is the class's own, so je_ratio = 1.
                  ! Stem and root respiration stay with Rd at tstar (as
                  ! two-leaf): sum the layer's tstar Rd here.
                  DO m = 1,veg_pts
                    l = veg_index(m)
                    rdc_tstar(l) = rdc_tstar(l)                                &
                                   + ( fsun(l,n) * rd_sun(l)                   &
                                       + (1.0 - fsun(l,n)) * rd_shd(l) )       &
                                   * dlai(l)
                  END DO
                  ! This layer's own sunlit and shaded leaf temperatures.
                  t_sun_lt(:) = t_lt_c(:,2*n-1)
                  t_shd_lt(:) = t_lt_c(:,2*n)
                  CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,   &
                         pft_photo_model, t_sun_lt, oa, actj, actv, dsj, dsv,  &
                         jv25, nleaf_layer, ones_pts, i2_sun, pstar, q_c_lt,   &
                         dq_min, ccp_sun_lt, kc_sun_lt, ko_sun_lt, km_sun_lt,  &
                         vcmax_sun_ml, je_sun_ml, rd_sun_ml, qs_sun_lt,        &
                         dq_sun_lt, f_vcmax_down = f_vcmax_down )
                  CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,   &
                         pft_photo_model, t_shd_lt, oa, actj, actv, dsj, dsv,  &
                         jv25, nleaf_layer, ones_pts, i2_shd, pstar, q_c_lt,   &
                         dq_min, ccp_shd_lt, kc_shd_lt, ko_shd_lt, km_shd_lt,  &
                         vcmax_shd_ml, je_shd_ml, rd_shd_ml, qs_shd_lt,        &
                         dq_shd_lt, f_vcmax_down = f_vcmax_down )
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

                  ! The leaf boundary layer of the layer's sunlit and shaded
                  ! leaves per unit leaf area, 1 / gbw: the class's
                  ! boundary layer conductance (class basis) over its leaf
                  ! area.
                  DO m = 1,veg_pts
                    l = veg_index(m)
                    r_bl_sun_ml(l) = 0.0
                    r_bl_shd_ml(l) = 0.0
                    IF ( gbw_lt_c(l,2*n-1) > 0.0 )                             &
                      r_bl_sun_ml(l) = lai_lt_c(l,2*n-1) / gbw_lt_c(l,2*n-1)
                    IF ( gbw_lt_c(l,2*n) > 0.0 )                               &
                      r_bl_shd_ml(l) = lai_lt_c(l,2*n) / gbw_lt_c(l,2*n)
                  END DO

                  IF ( l_cap_store ) CALL store_cpl_ml_q( n, k )
                  ! Sunlit, then shaded leaf of the layer, each at its own
                  ! leaf temperature.
                  CALL ml_class( vcmax_sun_ml, kc_sun_lt, ko_sun_lt,           &
                                 ccp_sun_lt, km_sun_lt, dq_sun_lt, je_sun_ml,  &
                                 t_sun_lt, ones_pts, fapar_sun(:,n),           &
                                 r_bl_sun_ml, rd_sun, ci_sun_ml, anetl_sun,    &
                                 el_sun, flux_o3_l_sun, fo3_l_sun, gl_sun,     &
                                 psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun,     &
                                 psi_stem_sun, q_store_sun_lyr )
                  CALL ml_class( vcmax_shd_ml, kc_shd_lt, ko_shd_lt,           &
                                 ccp_shd_lt, km_shd_lt, dq_shd_lt, je_shd_ml,  &
                                 t_shd_lt, ones_pts, fapar_shd(:,n),           &
                                 r_bl_shd_ml, rd_shd, ci_shd_ml, anetl_shd,    &
                                 el_shd, flux_o3_l_shd, fo3_l_shd, gl_shd,     &
                                 psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd,     &
                                 psi_stem_shd, q_store_shd_lyr )
                  IF ( l_cap_store ) THEN
                    el_sun_lyr(:,n) = el_sun(:)
                    el_shd_lyr(:,n) = el_shd(:)
                  END IF
                ELSE
                ! share_sup: this layer's kmax share, set with kmax_leaf_lyr.
                IF ( l_cap_store ) CALL store_cpl_ml_q( n, k )
                CALL ml_class( vcmax, kc, ko, ccp, km, dqc, je, tstar,         &
                               je_sun_ratio, fapar_sun(:,n), zeros_pts,        &
                               rd_sun, ci_sun_ml, anetl_sun, el_sun,           &
                               flux_o3_l_sun, fo3_l_sun, gl_sun, psi_leaf_sun, &
                               CG_sun, HC_sun, leaf_k_sun, psi_stem_sun,       &
                               q_store_sun_lyr )
                CALL ml_class( vcmax, kc, ko, ccp, km, dqc, je, tstar,         &
                               je_shd_ratio, fapar_shd(:,n), zeros_pts,        &
                               rd_shd, ci_shd_ml, anetl_shd, el_shd,           &
                               flux_o3_l_shd, fo3_l_shd, gl_shd, psi_leaf_shd, &
                               CG_shd, HC_shd, leaf_k_shd, psi_stem_shd,       &
                               q_store_shd_lyr )
                IF ( l_cap_store ) THEN
                  el_sun_lyr(:,n) = el_sun(:)
                  el_shd_lyr(:,n) = el_shd(:)
                END IF
                END IF

                ! Closed leaves get min_gl_pft from stom_opt_mod, a canopy
                ! value: divide by LAI so the canopy total (sum of gl * dlai
                ! over the layers) is min_gl_pft, as in big leaf and two-leaf.
                DO i = 1,clos_pts
                    l = veg_index(clos_index(i))
                    IF ( lai(l) > EPSILON(0.0) ) THEN
                        gl_sun(l) = gl_sun(l) / lai(l)
                        gl_shd(l) = gl_shd(l) / lai(l)
                    END IF
                END DO
                ! The same for open leaves the optimiser closed (no flux,
                ! gl = min_gl_pft from stom_opt_mod).
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
        !$OMP        leaf_k_sun,leaf_k_shd,ci_sum_gl_ml,ci_sum_lai_ml,         &
        !$OMP        gl_sun,gl_shd,                                            &
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

            ! These are ground-area integrated fluxes.
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

            IF ( l_som_plant_segments )                                        &
              psi_stem(l) = psi_stem(l) +                                      &
                            (fsun(l,n) * psi_stem_sun(l) +                     &
                            (1.0 - fsun(l,n)) * psi_stem_shd(l)) * dlai(l)

            ci_sum_gl_ml(l) = ci_sum_gl_ml(l) + (fsun(l,n) * gl_sun(l) *       &
                          ci_sun_ml(l) + (1.0 - fsun(l,n)) * gl_shd(l) *       &
                          ci_shd_ml(l)) * dlai(l)
            ci_sum_lai_ml(l) = ci_sum_lai_ml(l) + (fsun(l,n) * ci_sun_ml(l) +  &
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
            ! Net leaf photosynthesis of the top and the bottom layer
            ! (l_photo_diag, canopy_diag_mod).
            IF ( l_photo_diag ) THEN
              DO m = 1,veg_pts
                l = veg_index(m)
                IF ( n == 1 ) anetl_top_pft(l,ft) = anetl(l)
                IF ( n == ilayers ) anetl_bot_pft(l,ft) = anetl(l)
              END DO
            END IF

            ! l_leaf_temp: the conductances (m s-1, ground) of this layer's
            ! classes for the leaf energy balance.
            IF ( l_leaf_temp ) THEN
              DO m = 1,veg_pts
                l = veg_index(m)
                gl_lt_c(l,2*n-1) = fsun(l,n) * gl_sun(l) * dlai(l)
                gl_lt_c(l,2*n)   = (1.0 - fsun(l,n)) * gl_shd(l) * dlai(l)
              END DO
            END IF

        END DO                   ! N LAYERS

        ! Canopy conductance of this pass, for the next pass's deficit.
        gc_ml_prev(:) = gc(:)

        ! The store coupled to the optimiser: its flow out for the next
        ! pass, from its balance at this pass's canopy E, where some point
        ! transpires (the last pass's balance is committed after the leak
        ! block), and this pass's E for the layers' parts of q_store.
        IF ( l_cap_store .AND. k < n_pass_ml ) THEN
          DO m = 1,veg_pts
            l              = veg_index(m)
            e_store_bal(l) = MAX(el(l), 0.0)
            e_ml_prev(l)   = e_store_bal(l)
          END DO
          IF ( ANY(e_store_bal(veg_index(1:veg_pts)) > 0.0) ) THEN
            CALL store_cpl_q( .FALSE., e_store_bal )
            DO m = 1,veg_pts
              l = veg_index(m)
              ! A point still closed keeps last step's q_store, often a night
              ! refill (< 0): floored at 0, so if it opens in the next pass it
              ! is not searched with its root flow raised by the refill.
              IF ( store_w(l,ft) <= 0.0 .AND. lai(l) >= lai_store_min ) THEN
                IF ( e_store_bal(l) > 0.0 ) THEN
                  q_store(l) = q_store_bal(l)
                ELSE
                  q_store(l) = MAX(q_store(l), 0.0)
                END IF
              END IF
            END DO
            CALL store_cpl_supply()
          END IF
        END IF

        ! New leaf temperatures of every layer's sunlit and shaded leaf for
        ! the next pass (none after the last pass, so the temperatures
        ! output are those used), as two-leaf.
        IF ( l_leaf_temp .AND. k < n_pass_ml ) THEN
          CALL leaf_temp_update_n( land_pts, veg_pts, veg_index, n_cls_lt,     &
                                   tair, q1, pstar, lw_down, u_wind, lai_lt_c, &
                                   sw_lt_c, lwwt_sky_lt_c, lwwt_soil_lt_c,     &
                                   grad_fac_lt_c, lai_wind_lt_c, gl_lt_c,      &
                                   t_lt_c, gbw_lt_c, n_clamp_lt )
        END IF

        END DO                   ! K-ITER (whole canopy)

        ! The canopy's stomatal transpiration for the next step's first
        ! pass.
        IF ( l_cap_store ) CALL store_cpl_e_prev( el )

        ! Calculate canopy mean properties
        DO m = 1, veg_pts
            l = veg_index(m)

            ! LAI-weighted canopy means. Carbon gain and hydraulic cost are
            ! per-leaf 0-1 quantities, as in the big-leaf and two-leaf
            ! schemes; their dlai-weighted sums range 0..LAI.
            IF (lai(l) > 0.0) THEN
                psi_leaf(l)       = psi_leaf(l)       / lai(l)
                leaf_k(l)         = leaf_k(l)         / lai(l)
                carbon_gain(l)    = carbon_gain(l)    / lai(l)
                hydraulic_cost(l) = hydraulic_cost(l) / lai(l)
                IF ( l_som_plant_segments )                                    &
                  psi_stem(l) = psi_stem(l) / lai(l)
            ELSE IF ( l_som_plant_segments ) THEN
                psi_stem(l) = psi_zf(l)
            END IF

            ! gl-weighted canopy ci (as two-leaf); the LAI mean when no
            ! stomata are open (closed leaves have ci = ca).
            IF ( gc(l) > TINY(gc(l)) ) THEN
                ci(l) = ci_sum_gl_ml(l) / gc(l)
            ELSE IF ( lai(l) > 0.0 ) THEN
                ci(l) = ci_sum_lai_ml(l) / lai(l)
            ELSE
                ci(l) = ca(l)
            END IF

            ! l_leaf_temp: leaf temperatures used and the transpiration-
            ! equivalent gc (as two-leaf), capped at som_gl_max x LAI.
            IF ( l_leaf_temp ) THEN
              IF ( lai(l) > EPSILON(0.0) ) THEN
                f_sun_2l = lai_sun_2l(l) / lai(l)
              ELSE
                f_sun_2l = 0.5
              END IF
              gc_eq_cap = gl_leaf_cap * lai(l)
              ! The sunlit (shaded) temperature output is the leaf-area
              ! weighted mean of the layers' sunlit (shaded) leaves (so
              ! t_can is that of all leaves); a class without leaf area
              ! takes the air's (and lt_finish then the other class's). The
              ! clamp fraction is over all the classes with leaf area.
              t_sun_lt(l)    = 0.0
              t_shd_lt(l)    = 0.0
              n_cls_with_lai = 0.0
              DO n = 1,ilayers
                t_sun_lt(l) = t_sun_lt(l)                                    &
                              + lai_lt_c(l,2*n-1) * t_lt_c(l,2*n-1)
                t_shd_lt(l) = t_shd_lt(l)                                    &
                              + lai_lt_c(l,2*n) * t_lt_c(l,2*n)
                IF ( lai_lt_c(l,2*n-1) > EPSILON(0.0) )                      &
                  n_cls_with_lai = n_cls_with_lai + 1.0
                IF ( lai_lt_c(l,2*n) > EPSILON(0.0) )                        &
                  n_cls_with_lai = n_cls_with_lai + 1.0
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
              CALL lt_finish( l, f_sun_2l, n_cls_with_lai )
            END IF

        END DO
        IF ( l_leaf_temp ) THEN
          DEALLOCATE( lai_lt_c, sw_lt_c, lwwt_sky_lt_c, lwwt_soil_lt_c,        &
                      grad_fac_lt_c, lai_wind_lt_c, gl_lt_c, t_lt_c,           &
                      gbw_lt_c )
        END IF

        ! Xylem impairment state update, once per timestep (after the
        ! humidity iteration), as two-leaf: the canopy (leaf-area-weighted)
        ! psi_leaf and, with segments, psi_stem drive one damage state per
        ! point, at the lit points. With the plant capacitance the update
        ! is made after the store, below.
        IF ( l_ml_ximp .AND. .NOT. l_cap_store ) THEN
          n_lit = 0
          DO m = 1,veg_pts
            l = veg_index(m)
            IF ( apar(l) > 0.0 ) THEN
              n_lit = n_lit + 1
              open_land_index(n_lit) = l
            END IF
          END DO
          CALL update_xylem_impairment(                                        &
                ! IN
                  land_pts, n_lit, open_land_index, ft, psi_leaf, psi_src,     &
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
                                  nleaf_layer, qtenf_term, vcmax_temp, jmax,   &
                                  rd_dark, vcmax,                              &
                                  f_vcmax_down = f_vcmax_down )

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
  ,                       acr, apar, ca, ccp, dqc, fsmc, je, kc, km, ko        &
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

      ! Leaf fluxes by leaf_flux_mod: leaf_flux_fsmc here (the stomatal
      ! optimisation is the branch above); any other value stops the run.
      !
      ! MDK, 29/04/26
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
  ! Net leaf photosynthesis of the top and the bottom layer
  ! (l_photo_diag, canopy_diag_mod).
  !
  ! MDK, 30/10/26
  IF ( l_photo_diag ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( n == 1 ) anetl_top_pft(l,ft) = anetl(l)
      IF ( n == ilayers ) anetl_bot_pft(l,ft) = anetl(l)
    END DO
  END IF

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
                              jmax, rd, vcmax, f_vcmax_down = f_vcmax_down )

  IF ( pft_photo_model == photo_farquhar ) THEN
    ! Calculate the electron flux.
    CALL calc_electron_flux( land_pts, veg_pts, veg_index, i2, jmax, je)
  END IF

  ! The store coupled to the optimiser (the capacitance store), big leaf:
  ! the store's flow out the first pass sees, last step's, and the soil
  ! supply cap with it (as the two-leaf canopy). e_sup_use is e_supply
  ! otherwise.
  !
  ! MDK, 25/10/26
  e_sup_use(:) = e_supply(:)
  IF ( l_cap_store ) THEN
    ! (It needs iter >= 2: iter is fixed at 3.)
    CALL store_cpl_first_q()
    CALL store_cpl_supply()
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
          ,                 acr, apar, ca, ccp, dqc, fsmc, je, kc, km, ko      &
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
       ! fapar_dummy only for that).
       !
       ! MDK, 29/04/26
       l_multilayer = .FALSE.

       ! stom_opt_mod's kmax/kcrit dummy arguments are land_pts-sized;
       ! broadcast the whole-plant scalars onto land_pts arrays here rather
       ! than passing the bare scalars kmax_pft(ft)/kcrit(ft) directly, which
       ! would rely on undefined Fortran sequence association between a
       ! scalar actual argument and an array dummy argument, and could read
       ! out of bounds of the kmax_pft/kcrit arrays.
       !
       ! kmax_pft/kcrit are the whole-plant hydraulic conductance per unit
       ! leaf area (mol m-2 s-1 Pa-1), as estimated from site data
       ! (E / (LAI * (psi_predawn - psi_midday)) on well-watered days) and
       ! as used by CABLE. e_leaf passed into stom_opt_mod is canopy scale
       ! (per m2 ground), so the conductance must be too, or psi_leaf_jls
       ! solves
       !
       !   psi_leaf = psi_root - e_leaf / k
       !
       ! against too small a k. All leaves draw on the plant's conductance in
       ! parallel, so the canopy value is kmax_pft * LAI, not kmax_pft * fpar,
       ! the leaf-N (vcmax-like) scaling, which would make kmax_pft a
       ! top-of-canopy value and understate the canopy conductance by
       ! LAI / fpar (1.3-1.7 for LAI 1.5-2.3, kpar 0.5) for a data-derived
       ! kmax_pft.
       kmax_bigleaf(:) = kmax_pft(ft) * lai(:)
       kcrit_bigleaf(:) = kcrit_of_path(ft, kmax_bigleaf(:))
       ! som_gl_max is per unit leaf area, so the canopy cap is x LAI, as
       ! kmax_bigleaf (and the multilayer and two-leaf totals).
       gl_max_bigleaf(:) = gl_leaf_cap * lai(:)
       share_sup(:) = 1.0
       ! Coupled transpiration: the deficit to the level-1 air through ra
       ! (shared with the soil evaporation, as in the tile flux).
       dq_opt(:) = dqc(:)
       r_ca_opt(:) = 0.0
       IF ( l_som_coupled_e ) THEN
         dq_opt(:) = dq(:)
         r_ca_opt(:) = ra(:)
       END IF
       CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_use,       &
                                share_sup, dq_opt, tstar, pstar,               &
                                gl_max_bigleaf, gl_max_eff,                    &
                                zeros_pts, r_ca_opt, g_soil_c )
       CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

       ! Xylem impairment: k_max/conductance_b/conductance_c hold
       ! the (leaf-basis) impaired vulnerability curve for each point, so
       ! scale the impaired k_max by LAI in exactly the same way as
       ! kmax_pft above, keeping kmax_pft(ft) * LAI as the unimpaired
       ! reference for the hydraulic cost.
       kmax_ref_bigleaf(:) = kmax_bigleaf(:)
       kmax_bigleaf(:) = k_max(:) * lai(:)
       ! In spin-up, or without an impairment model, the intact kmax: the
       ! impairment is off there (stom_opt_mod), so a damaged k_max from a
       ! dump must not act as if it were intact.
       !
       ! MDK, 24/10/26
       IF ( is_spinup .OR. pft_xylem_impairment_model(ft) ==                   &
                           xylem_impairment_none )                             &
         kmax_bigleaf(:) = kmax_ref_bigleaf(:)
       ! kplant_model = 1: root and stem per unit ground area, the whole
       ! plant in series. The same factor on the impaired and intact plant
       ! and kcrit, so the damage fraction (k_max / kmax_pft) is unchanged;
       ! the impairment update below uses the unscaled leaf-basis values.
       !
       ! MDK, 21/10/26
       IF ( l_kscale ) THEN
         kmax_bigleaf(:)     = kmax_bigleaf(:)     * kplant_scale(:)
         kmax_ref_bigleaf(:) = kmax_ref_bigleaf(:) * kplant_scale(:)
         kcrit_bigleaf(:)    = kcrit_bigleaf(:)    * kplant_scale(:)
       END IF

       ! The store coupled to the optimiser: the canopy's path for its
       ! balance, and its flow out held fixed for this search (one path,
       ! share 1).
       !
       ! MDK, 25/10/26
       IF ( l_cap_store ) THEN
         kmax_ref_store(:) = kmax_ref_bigleaf(:)
         kmax_store(:)     = kmax_bigleaf(:)
         kcrit_store(:)    = kcrit_bigleaf(:)
         CALL set_store_q_path( land_pts, veg_pts, veg_index, share_sup,       &
                                q_store )
       END IF

       CALL stom_opt_mod (                                                     &
              ! IN
                land_pts, ft, open_pts, open_index,                            &
                pft_photo_model, veg_pts, veg_index,                           &
                ca, psi_src, acrc, apar, oa, vcmaxc, kc, ko, ccp, pstar,       &
                km, dq_opt, je, tstar, je_dummy, fapar_dummy,                  &
                kmax_bigleaf, kcrit_bigleaf, gl_max_eff, ipar,                 &
                l_multilayer, zeros_pts, r_ca_opt, g_soil_c,                   &
                kmax_ref_bigleaf, conductance_b, conductance_c,                &
                psi_leaf_extreme, psi_root_extreme, .TRUE.,                    &
              ! IN OUT
                rdc,                                                           &
              ! OUT
                ci, anetc, el, flux_o3, fo3, gc, psi_leaf,                     &
                carbon_gain, hydraulic_cost, leaf_k,                           &
                psi_stem = psi_stem                                            &
        )
        ! The store's flow held for the call applies no longer.
        !
        ! MDK, 25/10/26
        IF ( l_cap_store ) CALL clear_store_q_path()

        ! The store coupled to the optimiser: its flow out for the next pass,
        ! from its balance at this pass's E, where some point transpires
        ! (the last pass's balance is committed after the leak block).
        !
        ! MDK, 25/10/26
        IF ( l_cap_store .AND. k < iter ) THEN
          DO m = 1,veg_pts
            l = veg_index(m)
            e_store_bal(l) = MAX(el(l), 0.0)
          END DO
          IF ( ANY(e_store_bal(veg_index(1:veg_pts)) > 0.0) ) THEN
            CALL store_cpl_q( .FALSE., e_store_bal )
            DO m = 1,veg_pts
              l = veg_index(m)
              ! A point still closed keeps last step's q_store, often a night
              ! refill (< 0): floored at 0, so if it opens in the next pass it
              ! is not searched with its root flow raised by the refill.
              IF ( store_w(l,ft) <= 0.0 .AND. lai(l) >= lai_store_min ) THEN
                IF ( e_store_bal(l) > 0.0 ) THEN
                  q_store(l) = q_store_bal(l)
                ELSE
                  q_store(l) = MAX(q_store(l), 0.0)
                END IF
              END IF
            END DO
            CALL store_cpl_supply()
          END IF
        END IF

        ! el stays in mol H2O m-2 s-1 (ground), as in every path: physiol uses
        ! it for the layer uptake at the actual root potential
        ! (l_som_rhizo_series).

        ! leaf_k from stom_opt_mod is canopy scale (see kmax_bigleaf). Return
        ! it per unit leaf area, as the multilayer path does (kplant_pft).
        leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))

        ! Update the impairment state once per timestep, on the final pass
        ! of the humidity-deficit iteration (DO k = 1,iter), with the
        ! converged state.
        IF (k == iter) THEN

          ! Apply the xylem impairment to the vulnerability curve.
          ! NOTE: open_index indexes veg_index, while update_xylem_impairment
          !       indexes the land arrays with it directly, so pass it the
          !       land points of the open stomata.
          n_lit = open_pts
          DO m = 1, n_lit
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
                    n_lit,                                                     &
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

  ! The canopy's stomatal transpiration for the next step's first pass.
  IF ( l_cap_store ) CALL store_cpl_e_prev( el )


CASE ( 7 )

  !---------------------------------------------------------------------------
  ! Two-leaf model: one sunlit and one shaded big leaf (de Pury & Farquhar,
  ! 1997; Wang & Leuning, 1998, as in CABLE), both at the surface
  ! temperature tstar, or each at its own leaf temperature with l_leaf_temp.
  !
  ! The sunlit fraction and the absorbed PAR per unit sunlit/shaded leaf area
  ! come from the JULES two-stream profile over ilayers (as for
  ! can_rad_mod = 5/6), and are summed over the layers into two canopy-scale
  ! leaves.
  !
  ! Leaf N, and so Vcmax, Jmax and Rd, declines with cumulative leaf area L
  ! as exp(-knl L), as in the multilayer scheme (can_rad_mod = 6). Layer n,
  ! of leaf area dL, has the N weight
  !
  !   w_n = ( exp(-knl (n-1) dL) - exp(-knl n dL) ) / knl,
  !
  ! and the classes take
  !
  !   N_sun = sum_n f_sun,n w_n,
  !   N_shd = sum_n (1 - f_sun,n) w_n,
  !
  ! with f_sun,n the sunlit fraction of layer n. Each class's Vcmax and Jmax
  ! are the top-leaf values times its N weight, so the canopy total is the
  ! top-leaf value times (1 - exp(-knl LAI)) / knl. Big-leaf (can_rad_mod =
  ! 1) scales by (1 - exp(-kpar LAI)) / kpar, so the canopy totals match
  ! big-leaf only when knl = kpar.
  !
  ! The absorbed light is the two-stream canopy total (<= incident PAR),
  ! whereas big-leaf uses
  !
  !   acr fpar = (1 - omega) ipar (1 - exp(-kpar LAI)) / kpar,
  !
  ! about 1/kpar times more, so expect lower light-limited GPP than big-leaf,
  ! beyond the effect of the sun/shade split. The leaf conductance cap is
  ! gl_leaf_cap times each class's leaf area.
  !
  ! kmax is per unit leaf area as in big-leaf. The canopy conductance
  ! kmax_pft LAI is shared between the classes by their N weights,
  !
  !   kmax_sun = kmax_pft LAI N_sun / (N_sun + N_shd),
  !
  ! and likewise for the shaded leaf, so the canopy total is unchanged. Each
  ! class draws from the source psi (psi_src) on its own path, or with
  ! l_som_shared_stem through a root and stem shared with the other class.
  !
  ! l_leaf_temp: each leaf has its own temperature from a leaf energy
  ! balance (Penman-Monteith; Leuning et al., 1995; Wang & Leuning, 1998),
  ! solved with its absorbed PAR, near-IR and isothermal long-wave, its
  ! boundary-layer and radiative conductances and its stomatal conductance
  ! from the last pass. The leaves exchange with the level-1 air, as in
  ! MAESPA. Each pass runs the photosynthesis
  ! temperature responses, the humidity deficit and the stomatal
  ! optimisation of each leaf at its own temperature, then updates the leaf
  ! temperatures: leaf_temp_iter passes, so no more stom_opt_mod calls than
  ! without l_leaf_temp (for leaf_temp_iter = iter).
  !
  ! MDK, 04/11/26
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
                              jmax, rd, vcmax, f_vcmax_down = f_vcmax_down )

  !---------------------------------------------------------------------------
  ! Aggregate the layer profile into the two leaf classes.
  !---------------------------------------------------------------------------
  CALL sun_shade_totals( land_pts, veg_pts, veg_index, ilayers, dlai, fsun,    &
                         fapar_sun, fapar_shd, ipar, lai_sun_2l, lai_shd_2l,   &
                         apar_sun_2l, apar_shd_2l )
  DO m = 1,veg_pts
    l                 = veg_index(m)
    nweight_sun_2l(l) = 0.0
    nweight_shd_2l(l) = 0.0
    ! N extinction per unit LAI (extkn of Wang & Leuning 1998): knl, as the
    ! multilayer scheme (can_rad_mod = 6). kpar is the PAR extinction only.
    kn_2l = knl(ft)
    DO n = 1,ilayers
      IF ( kn_2l > EPSILON(0.0) ) THEN
        dnweight_2l = ( EXP(-kn_2l * REAL(n-1) * dlai(l))                      &
                        - EXP(-kn_2l * REAL(n) * dlai(l)) ) / kn_2l
      ELSE
        dnweight_2l = dlai(l)
      END IF
      nweight_sun_2l(l) = nweight_sun_2l(l) + fsun(l,n) * dnweight_2l
      nweight_shd_2l(l) = nweight_shd_2l(l) + (1.0 - fsun(l,n)) * dnweight_2l
    END DO

    apar_mol_sun_2l(l) = apar_sun_2l(l) / conpar
    apar_mol_shd_2l(l) = apar_shd_2l(l) / conpar

    vcmax_sun_2l(l) = vcmax(l) * nweight_sun_2l(l)
    vcmax_shd_2l(l) = vcmax(l) * nweight_shd_2l(l)
    jmax_sun_2l(l)  = jmax(l)  * nweight_sun_2l(l)
    jmax_shd_2l(l)  = jmax(l)  * nweight_shd_2l(l)

    ! Canopy conductance kmax_pft*LAI (per unit leaf area, as big-leaf),
    ! shared between the classes by their N-weighted leaf area (CABLE's
    ! scalex), so sun leaves get more conductance per unit leaf area than
    ! shade leaves without changing the canopy total (cf. the multilayer
    ! leaf-segment profile). CABLE itself uses kmax_pft * scalex, which
    ! lowers the canopy total by fpar/LAI.
    !
    ! Xylem impairment: both classes draw on the one (impaired) plant, so
    ! the impaired leaf-basis k_max (= kmax_pft when unimpaired) is shared
    ! the same way, and the unimpaired kmax_pft share is the reference for
    ! the hydraulic cost; kcrit stays on the unimpaired curve, as big-leaf.
    nweight_tot_2l = MAX(nweight_sun_2l(l) + nweight_shd_2l(l), TINY(1.0))
    kmax_ref_sun_2l(l) = kmax_pft(ft) * lai(l) * nweight_sun_2l(l) /           &
                         nweight_tot_2l
    kmax_ref_shd_2l(l) = kmax_pft(ft) * lai(l) * nweight_shd_2l(l) /           &
                         nweight_tot_2l
    kmax_sun_2l(l) = k_max(l) * lai(l) * nweight_sun_2l(l) / nweight_tot_2l
    kmax_shd_2l(l) = k_max(l) * lai(l) * nweight_shd_2l(l) / nweight_tot_2l

    ! In spin-up, or without an impairment model, the intact kmax (as big
    ! leaf).
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
    ! (so the damage fraction is unchanged). The impairment update gets the
    ! unscaled k_max, so it never sees the factor.
    IF ( l_kscale ) THEN
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
    i2_sun(l) = alpha_elec(ft) * apar_mol_sun_2l(l)
    i2_shd(l) = alpha_elec(ft) * apar_mol_shd_2l(l)
  END DO

  ! Each class's electron transport comes from its totals of I and Jmax.
  ! The light response,
  !
  !   J = ( I + Jmax - SQRT((I + Jmax)**2 - 4 theta I Jmax) ) / (2 theta),
  !
  ! with theta = light_curvature (calc_electron_flux), scales with I and
  ! Jmax together, J(c I, c Jmax) = c J(I, Jmax), so J of
  ! the totals equals the sum of J over the class's leaves when they all
  ! have the same I / Jmax. Within a class the ratio varies (light and N
  ! fall at different rates with depth); J is concave, so J of the totals
  ! is then at least the sum over leaves. This is the usual two-leaf
  ! approximation (de Pury & Farquhar, 1997).
  !
  ! MDK, 03/11/26
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
    CALL lt_init()
    n_pass_2l = leaf_temp_iter
    ! Store the tstar-based canopy respiration for the stem and root
    ! respiration, which stay at tstar.
    DO m = 1,veg_pts
      l = veg_index(m)
      rdc_tstar(l) = rd(l) * ( nweight_sun_2l(l) + nweight_shd_2l(l) )
    END DO
  ELSE
    n_pass_2l = iter
    ! Without leaf temperatures both classes are at tstar, with the tile's
    ! photosynthesis terms, so the class routines take the class arrays
    ! either way (and the 13C discrimination reads set values).
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
  ! The store coupled to the optimiser updates q_store between passes, so it
  ! needs at least two: iter and leaf_temp_iter are both fixed at 3.

  !---------------------------------------------------------------------------
  ! Coupled transpiration: first guesses of the class conductances and leaf
  ! boundary layers (the last timestep's) and the canopy air resistance.
  !---------------------------------------------------------------------------
  IF ( l_som_coupled_e ) THEN
    CALL cpl_state_load( ft, land_pts, veg_pts, veg_index, lai, lai_sun_2l,    &
                         lai_shd_2l, gc, gl_sun, gl_shd, gbw_sun_lt,           &
                         gbw_shd_lt )
    IF ( l_leaf_temp ) THEN
      ! The leaves are in the level-1 air.
      ra_ca_lt(:) = 0.0
    ELSE
      ra_ca_lt(:) = ra(:)
    END IF
  END IF

  !---------------------------------------------------------------------------
  ! The soil supply cap of the optimiser, and with the store coupled to it
  ! (the capacitance store) the store's flow out it sees in the first pass
  ! (store_cpl_first_q).
  !
  ! MDK, 04/11/26
  !---------------------------------------------------------------------------
  e_sup_use(:) = e_supply(:)
  IF ( l_cap_store ) THEN
    CALL store_cpl_first_q()
    CALL store_cpl_supply()
  END IF
  share_one(:) = 1.0
  q_store_sun(:)   = 0.0
  q_store_shd(:)   = 0.0
  IF ( l_cap_store ) THEN
    ! The canopy's path for the store's balance: the two classes together.
    DO m = 1,veg_pts
      l                 = veg_index(m)
      kmax_ref_store(l) = kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l)
      kmax_store(l)     = kmax_sun_2l(l) + kmax_shd_2l(l)
      kcrit_store(l)    = kcrit_sun_2l(l) + kcrit_shd_2l(l)
    END DO
    CALL store_cpl_split( .FALSE. )
  END IF

  !---------------------------------------------------------------------------
  ! Iterate for the canopy humidity deficit, starting from the previous
  ! timestep's conductance (as for big-leaf), and with l_leaf_temp for the
  ! leaf temperatures.
  !---------------------------------------------------------------------------
  tleaf_a_2l(:,:) = 0.0
  tleaf_b_2l(:,:) = 0.0

  ! The shared root and stem: the last timestep's shaded flow (0 before the
  ! first).
  e_oth_2l(:) = 0.0
  a_oth_2l(:) = 0.0
  amax_sun_2l(:) = 0.0
  amax_shd_2l(:) = 0.0
  IF ( l_som_shared_stem ) THEN
!$OMP CRITICAL (som_shared_stem_alloc)
    IF ( .NOT. ALLOCATED(shared_e_shd_prev) ) THEN
      ALLOCATE( shared_e_shd_prev(land_pts, npft) )
      ALLOCATE( shared_amax_shd_prev(land_pts, npft) )
      shared_e_shd_prev(:,:) = 0.0
      shared_amax_shd_prev(:,:) = 0.0
    END IF
!$OMP END CRITICAL (som_shared_stem_alloc)
  END IF
  DO k = 1,n_pass_2l

    DO m = 1,veg_pts
      l = veg_index(m)
      ra_rc(l) = ra(l) * gc(l)
      dqc(l)   = dq(l) / (1.0 + ra_rc(l))
      ! No light inhibition of Rd, as for big-leaf (rd is reset each pass
      ! because stom_opt_mod treats it as IN OUT).
      rd_sun(l) = rd(l) * nweight_sun_2l(l)
      rd_shd(l) = rd(l) * nweight_shd_2l(l)
    END DO

    IF ( l_leaf_temp ) THEN
      ! Photosynthesis parameters and humidity deficit at each leaf's
      ! temperature. The deficit is from the leaf to the level-1 air,
      ! qs(T_leaf) - q1.
      CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_sun_lt, oa, actj, actv, dsj, &
                               dsv, jv25, nleaf_top, nweight_sun_2l, i2_sun,   &
                               pstar, q_c_lt, dq_min, ccp_sun_lt, kc_sun_lt,   &
                               ko_sun_lt, km_sun_lt, vcmax_sun_2l, je_sun,     &
                               rd_sun, qs_sun_lt, dq_sun_lt,                   &
                               f_vcmax_down = f_vcmax_down )
      CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_shd_lt, oa, actj, actv, dsj, &
                               dsv, jv25, nleaf_top, nweight_shd_2l, i2_shd,   &
                               pstar, q_c_lt, dq_min, ccp_shd_lt, kc_shd_lt,   &
                               ko_shd_lt, km_shd_lt, vcmax_shd_2l, je_shd,     &
                               rd_shd, qs_shd_lt, dq_shd_lt,                   &
                               f_vcmax_down = f_vcmax_down )
    END IF

    !-------------------------------------------------------------------------
    ! Sunlit, then shaded leaf (the shaded leaf sees this pass's sunlit gl
    ! through coupled_inputs).
    !-------------------------------------------------------------------------
    ! l_leaf_temp_ci: from the second pass, each candidate Ci of the
    ! search at its own leaf temperature (ltc_set).
    IF ( l_leaf_temp_ci .AND. k > 1 ) CALL ltc_set( 1, t_sun_lt,               &
      nweight_sun_2l, i2_sun, ccp_sun_lt, kc_sun_lt, ko_sun_lt, km_sun_lt,     &
      vcmax_sun_2l, je_sun, rd_sun, qs_sun_lt )
    ! The shared root and stem: the shaded leaves' flow on it, of the last
    ! pass, or in the first pass of the last timestep.
    IF ( l_som_shared_stem ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)
        IF ( k == 1 ) THEN
          e_oth_2l(l) = shared_e_shd_prev(l,ft)
          a_oth_2l(l) = shared_amax_shd_prev(l,ft)
        ELSE
          e_oth_2l(l) = MAX(el_shd(l), 0.0)
          a_oth_2l(l) = amax_shd_2l(l)
        END IF
      END DO
    END IF
    CALL two_leaf_class( apar_sun_2l, lai_sun_2l, kmax_ref_sun_2l, qs_sun_lt,  &
                         dq_sun_lt, gbw_sun_lt, gbw_shd_lt, gl_shd, t_sun_lt,  &
                         gl_max_sun_2l, apar_mol_sun_2l, vcmax_sun_2l,         &
                         kc_sun_lt, ko_sun_lt, ccp_sun_lt, km_sun_lt, je_sun,  &
                         kmax_sun_2l, kcrit_sun_2l, rd_sun, ci_sun_2l,         &
                         anetl_sun, el_sun, flux_o3_l_sun, fo3_l_sun, gl_sun,  &
                         psi_leaf_sun, CG_sun, HC_sun, leaf_k_sun,             &
                         psi_stem_sun, lam_sun, q_store_sun, e_oth_2l,         &
                         a_oth_2l, amax_sun_2l )
    IF ( l_leaf_temp_ci .AND. k > 1 ) THEN
      CALL clear_leaf_temp_ci()
      CALL ltc_set( 2, t_shd_lt, nweight_shd_2l, i2_shd, ccp_shd_lt,           &
                    kc_shd_lt, ko_shd_lt, km_shd_lt, vcmax_shd_2l, je_shd,     &
                    rd_shd, qs_shd_lt )
    END IF
    ! The shared root and stem: the sunlit leaves' flow of this pass.
    IF ( l_som_shared_stem ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)
        e_oth_2l(l) = MAX(el_sun(l), 0.0)
        a_oth_2l(l) = amax_sun_2l(l)
      END DO
    END IF
    CALL two_leaf_class( apar_shd_2l, lai_shd_2l, kmax_ref_shd_2l, qs_shd_lt,  &
                         dq_shd_lt, gbw_shd_lt, gbw_sun_lt, gl_sun, t_shd_lt,  &
                         gl_max_shd_2l, apar_mol_shd_2l, vcmax_shd_2l,         &
                         kc_shd_lt, ko_shd_lt, ccp_shd_lt, km_shd_lt, je_shd,  &
                         kmax_shd_2l, kcrit_shd_2l, rd_shd, ci_shd_2l,         &
                         anetl_shd, el_shd, flux_o3_l_shd, fo3_l_shd, gl_shd,  &
                         psi_leaf_shd, CG_shd, HC_shd, leaf_k_shd,             &
                         psi_stem_shd, lam_shd, q_store_shd, e_oth_2l,         &
                         a_oth_2l, amax_shd_2l )
    IF ( l_leaf_temp_ci .AND. k > 1 ) CALL clear_leaf_temp_ci()

    DO m = 1,veg_pts
      l = veg_index(m)
      gc(l) = gl_sun(l) + gl_shd(l)
    END DO

    ! The store coupled to the optimiser (the capacitance store): the
    ! store's flow out for the next pass, from its balance over the step at
    ! this pass's canopy E (a fixed-point iteration on q_store; the last pass's
    ! balance is committed after the leak block). Only where some point
    ! transpires: q_store is used only by the next pass's open points, so
    ! with E = 0 everywhere, e.g. at night, the solve is skipped and q_store
    ! keeps its value. The IFs are nested: Fortran need not short-circuit
    ! .AND., and e_store_bal is set only with the store coupled.
    IF ( l_cap_store .AND. k < n_pass_2l ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)
        e_store_bal(l) = MAX(el_sun(l) + el_shd(l), 0.0)
      END DO
      IF ( ANY(e_store_bal(veg_index(1:veg_pts)) > 0.0) ) THEN
        CALL store_cpl_q( .FALSE., e_store_bal )
        DO m = 1,veg_pts
          l = veg_index(m)
          ! A point still closed keeps last step's q_store, often a night refill
          ! (< 0): floored at 0, so if it opens in the next pass it is
          ! not searched with its root flow raised by the refill.
          IF ( store_w(l,ft) <= 0.0 .AND. lai(l) >= lai_store_min ) THEN
            IF ( e_store_bal(l) > 0.0 ) THEN
              q_store(l) = q_store_bal(l)
            ELSE
              q_store(l) = MAX(q_store(l), 0.0)
            END IF
          END IF
        END DO
        CALL store_cpl_supply()
        CALL store_cpl_split( .TRUE. )
      END IF
    END IF

    ! New leaf temperatures for the next pass (none after
    ! the last pass, so the temperatures output are those used).
    IF ( l_leaf_temp .AND. k < n_pass_2l ) THEN
      CALL leaf_temp_update( land_pts, veg_pts, veg_index, tair, q1, pstar,    &
                             lw_down, u_wind, lai_sun_2l, lai_shd_2l,          &
                             sw_sun_lt, sw_shd_lt, lwwt_sky_sun_lt,            &
                             lwwt_sky_shd_lt, lwwt_soil_sun_lt,                &
                             lwwt_soil_shd_lt, grad_fac_sun_lt,                &
                             grad_fac_shd_lt, lai_wind_sun_lt,                 &
                             lai_wind_shd_lt, gl_sun, gl_shd, t_sun_lt,        &
                             t_shd_lt, gbw_sun_lt, gbw_shd_lt, n_clamp_lt,     &
                             tleaf_a_2l, tleaf_b_2l )
    END IF

  END DO   ! End of iteration loop

  IF ( l_som_coupled_e ) CALL cpl_state_save( ft, land_pts, veg_pts, veg_index,&
                                              gl_sun, gl_shd, gbw_sun_lt,      &
                                              gbw_shd_lt, ra_ca_lt )

  ! The canopy's stomatal transpiration for the next step's first pass.
  !
  ! MDK, 01/11/26
  IF ( l_cap_store ) CALL store_cpl_e_prev( el_sun + el_shd )

  ! The shared root and stem (l_som_shared_stem): the shaded leaves' flow
  ! (and A at their edge) for the next timestep's first pass, and one stem
  ! water potential for both classes. The shaded leaves' search, the last,
  ! put this pass's flow of both classes on the stem, so its stem psi is the
  ! shared stem's; where the shaded leaves do not transpire, the sunlit
  ! leaves' (whose search had the shaded leaves' flow of the last pass on
  ! it).
  IF ( l_som_shared_stem ) THEN
    DO m = 1,veg_pts
      l = veg_index(m)
      shared_e_shd_prev(l,ft) = MAX(el_shd(l), 0.0)
      shared_amax_shd_prev(l,ft) = amax_shd_2l(l)
      IF ( el_shd(l) > 0.0 ) THEN
        psi_stem_sun(l) = psi_stem_shd(l)
      ELSE
        psi_stem_shd(l) = psi_stem_sun(l)
      END IF
    END DO
  END IF

  ! Diagnostic marginal water price of the canopy (som_lambda_eff_pft):
  ! the sunlit and shaded leaves' dA/dE weighted by their transpiration.
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

    ! The sunlit and shaded leaves' stomatal conductance per unit leaf
    ! area of each class, as the optimiser chose it (gl_sun and gl_shd are
    ! class totals per unit ground).
    lai_sun_pft(l,ft) = lai_sun_2l(l)
    IF ( lai_sun_2l(l) > EPSILON(0.0) )                                        &
      gs_sun_pft(l,ft) = gl_sun(l) / lai_sun_2l(l)
    IF ( lai_shd_2l(l) > EPSILON(0.0) )                                        &
      gs_shd_pft(l,ft) = gl_shd(l) / lai_shd_2l(l)
    ! and their leaf water potentials.
    psi_leaf_sun_pft(l,ft) = psi_leaf_sun(l)
    psi_leaf_shd_pft(l,ft) = psi_leaf_shd(l)

    ! The classes' leaf water potentials and shares of the leaf area, for
    ! the leaf segment's damage (the store commit below may reset the
    ! potentials).
    !
    ! MDK, 04/11/26
    psi_leaf_cls(l,1) = psi_leaf_sun(l)
    psi_leaf_cls(l,2) = psi_leaf_shd(l)
    w_cls(l,1)        = f_sun_2l
    w_cls(l,2)        = 1.0 - f_sun_2l

    IF ( l_leaf_temp ) THEN
      ! Leaf temperatures used and the transpiration-equivalent gc
      ! (lt_finish, shared with the multilayer scheme).
      gc_eq_cap = gl_max_sun_2l(l) + gl_max_shd_2l(l)
      CALL lt_finish( l, f_sun_2l )
    END IF
  END DO

  ! The shared root and stem without the store: the plant's conductance
  ! from one solve of the whole plant (kplant_shared_stem). With the store
  ! it is the store's kplant_bal, set below.
  !
  ! MDK, 04/11/26
  IF ( l_som_shared_stem .AND. .NOT. l_cap_store ) CALL kplant_shared_stem()

  ! As the big-leaf path: el in mol H2O m-2 s-1 and leaf_k per unit leaf
  ! area. (gc_som is set after the canopy switch.)
  leaf_k(:) = leaf_k(:) / MAX(lai(:), TINY(1.0_real_jlslsm))

  ! Xylem impairment state update, once per timestep (after the humidity
  ! iteration, cf. k == iter in the big-leaf path), driven by the canopy
  ! (leaf-area-weighted) psi_leaf: both classes share the one plant, so
  ! there is a single damage state per point, as for the big leaf. With
  ! per-segment memory the leaf segment's damage is the leaf-area-weighted
  ! mean of the classes' relative conductances at their own leaf psi
  ! (psi_leaf_cls, w_cls).
  !
  ! MDK, 04/11/26
  n_lit = 0
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( apar_sun_2l(l) > 0.0 .OR. apar_shd_2l(l) > 0.0 ) THEN
      n_lit = n_lit + 1
      open_land_index(n_lit) = l
    END IF
  END DO
  ! With the plant store or the cuticular loss (l_cap_store), the update
  ! is made after it, below.
  IF ( .NOT. l_cap_store ) THEN
    CALL update_xylem_impairment(                                              &
          ! IN
            land_pts, n_lit, open_land_index, ft, psi_leaf, psi_src,           &
            psi_stem,                                                          &
          ! IN OUT
            k_max, psi_leaf_extreme, psi_root_extreme,                         &
          ! IN (recovery)
            lai, canht, anetc,                                                 &
          ! IN (leaf classes)
            psi_leaf_cls = psi_leaf_cls, w_cls = w_cls )
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
! replaced it by gc_eq (then lt_finish kept it).
IF ( .NOT. l_leaf_temp ) gc_stomatal_pft(:,ft) = gc(:)

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
! conductance (kcrit_fractional_loss, e.g. 0.95, beyond P88); mortality, if
! on (l_trif_hydmort), acts through TRIFFID's cover from the accumulated
! loss.
!
! Where the optimiser has (nearly) shut the stomata, the canopy's flux is
! raised to the leak's, gcuticular * LAI times the leaf-to-air humidity
! deficit (dq_leak, at most dq_ratio_leak_max times the tile's deficit); the
! conductance that gives it is solved from the tile's flux. There is no
! carbon trade at that point, so A is unchanged; the leak goes into the
! transpiration (and so ET and LE) through gc. Its water comes from the
! capacitance store and the soil (the soil-supply cap applies: no water is
! taken that is not there). The store is described in
! plant_water_store_mod. In detail:
! - with som_cuticular_leak = 2 the closed canopy (gc = min_gl_pft) is
!   raised to the leak day and night; with som_cuticular_leak = 1 only at
!   lit points with apar >= apar_leak_min;
! - the leak is cut back only where it would take the leaf below the psi
!   of the critical loss of conductance on the leaf (or PFT) curve,
!   psi_floor_leaf (kcrit_fractional_loss, e.g. 0.95; store_cpl_leak_limit),
!   so in a drought it keeps drying the plant down to that psi;
! - the leak drains the capacitance store at the stem node, psi_leaf
!   follows it, and the impairment update is made after this block, so the
!   leak drives the loss of conductance.
! psi is not allowed below the psi of the critical loss of conductance
! (kcrit_fractional_loss, e.g. 0.95) on the PFT curve (or the leaf
! and stem segments), so where the path cannot carry the leak (no steady
! state without storage) psi does not fall to the solver's floor and the
! memory does not go past that loss. The transpiration is the tile's flux
! at tstar, as sf_evap takes it.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
IF ( l_cap_store ) THEN
  CALL store_state_alloc( land_pts )
  store_q(:,ft) = 0.0
  leak_e(:,ft)  = 0.0

  ! The multilayer canopy has no lit set from its own path: the points
  ! with absorbed PAR.
  IF ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) THEN
    n_lit = 0
    DO m = 1,veg_pts
      l = veg_index(m)
      IF ( apar(l) > 0.0 ) THEN
        n_lit = n_lit + 1
        open_land_index(n_lit) = l
      END IF
    END DO
  END IF
  l_leak_on(:)  = .FALSE.
  l_leak_lit(:) = .FALSE.
  DO m = 1,n_lit
    l = open_land_index(m)
    l_leak_lit(l) = apar(l) >= apar_leak_min
  END DO

  ! The leak is driven by the leaf-to-air deficit (as gmin * VPD): with
  ! leaf temperatures (two-leaf), qs at the sunlit and shaded leaf
  ! temperatures weighted by leaf area, less q1; otherwise the air deficit,
  ! qs(tair) - q1. Not dqc, at tstar, which in a drought is the hot tile
  ! skin, up to 10 K above the leaves (leak ~2x too large). sf_evap takes
  ! the flux at tstar, so gc is scaled by dq_leak / dqc. Through ra, shared
  ! with the soil conductance, the tile then delivers slightly less than
  ! gcuticular LAI dq_leak.

  CALL qsat(qs_air_leak, tair, pstar, land_pts)
  DO m = 1,veg_pts
    l = veg_index(m)
    dq_leak(l) = MAX(qs_air_leak(l) - q1(l), 0.0)
    IF ( l_leaf_temp .AND. can_rad_mod == 7 ) THEN
      IF ( lai_sun_2l(l) + lai_shd_2l(l) > EPSILON(0.0) )                      &
        dq_leak(l) = ( lai_sun_2l(l) * MAX(qs_sun_lt(l) - q1(l), 0.0)          &
                     + lai_shd_2l(l) * MAX(qs_shd_lt(l) - q1(l), 0.0) )        &
                   / ( lai_sun_2l(l) + lai_shd_2l(l) )
    END IF
  END DO
  DO m = 1,veg_pts
    l = veg_index(m)
    ! mmol m-2 leaf s-1 -> m s-1, times LAI for the canopy.
    gc_cuticular(l) = gcuticular(ft) * 1.0e-3 * rmol * tstar(l) / pstar(l)     &
                      * lai(l)
    ! No leak with capacitance alone (l_som_capacitance).
    IF ( .NOT. l_som_cuticular_loss ) gc_cuticular(l) = 0.0
    IF ( dqc(l) > 0.0 ) gc_cuticular(l) = gc_cuticular(l)                      &
                        * MIN(dq_leak(l) / dqc(l), dq_ratio_leak_max)
    share_sup(l) = 1.0
  END DO

  ! (Coupled to the optimiser the soil supplies only E - q_store: the cap is
  ! e_sup_use.)
  IF ( l_cap_store ) THEN
    CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_use,          &
                             share_sup, dqc, tstar, pstar, gc_cuticular,       &
                             gc_cuticular_capped )
  ELSE
    CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_supply,           &
                             share_sup, dqc, tstar, pstar, gc_cuticular,       &
                             gc_cuticular_capped )
  END IF
  CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

  ! The leak acts at every point, lit or not, or with
  ! l_som_cuticular_day_only only at the lit points with apar >=
  ! apar_leak_min. Without it, at night and in dim light the closed canopy
  ! (gc = min_gl_pft) is raised to the leak as by day. Its deficit is
  ! dq_leak above (the leaf temperatures, or the air, at night), and its
  ! flux the tile's at tstar, as sf_evap takes it, so where the surface
  ! has cooled below the air (clear nights) the leak is smaller than
  ! gcuticular * LAI * VPD, and nil with dew (qs(tstar) <= q1).
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( ( l_som_cuticular_day_only .AND. .NOT. l_leak_lit(l) )                &
         .OR. gc_cuticular(l) <= 0.0 ) CYCLE
    ! The canopy conductance that gives the leak's flux
    ! X = rho gc_cuticular_capped dqc (gcuticular LAI times the leaf-to-air
    ! deficit, after the supply cap) through the tile's resistances, as
    ! sf_evap takes it (the stomatal share of the flux through ra, shared
    ! with the soil conductance g_s; see gc_eq in lt_finish),
    !
    !   X = A gc / (1 + ra (gc + g_s)),   A = rho (qs(tstar) - q1),
    !   gc = X (1 + ra g_s) / (A - ra X),
    !
    ! or gc_cuticular_capped where A - ra X is near 0 (rounding).
    e_per_gc_leak = pstar(l) / ( r * tstar(l) ) * MAX(qs(l) - q1(l), 0.0)
    e_target_leak = pstar(l) / ( r * tstar(l) ) * gc_cuticular_capped(l)       &
                    * dqc(l)
    gc_leak = gc_cuticular_capped(l)
    IF ( e_per_gc_leak - ra(l) * e_target_leak > TINY(1.0_real_jlslsm) )       &
      gc_leak = e_target_leak * ( 1.0 + ra(l) * g_soil_c(l) )                  &
                / ( e_per_gc_leak - ra(l) * e_target_leak )
    IF ( gc(l) >= gc_leak ) CYCLE
    gc(l)     = gc_leak
    gc_som(l) = gc(l)
    ! Canopy transpiration (mol H2O m-2 s-1), the tile's flux at this gc,
    !
    !   E = rho (qs(tstar) - q1) gc / (1 + ra (gc + g_s)) / M_w,
    !
    ! which is X / M_w.
    el(l)     = pstar(l) / ( r * tstar(l) ) * MAX(qs(l) - q1(l), 0.0)         &
                * gc(l) / ( 1.0 + ra(l) * ( gc(l) + g_soil_c(l) ) )            &
                / m_h2o_lt
    l_leak_on(l) = .TRUE.
  END DO

  ! Lower limits on the leaf and stem psi for the store and the leak. Each
  ! is the psi at which the segment keeps only the share
  !
  !   r = 1 - kcrit_fractional_loss
  !
  ! of its conductance (the critical loss, e.g. 0.95), which is also the
  ! lowest conductance the impairment memory allows. On the Weibull curve,
  ! k / kmax = exp(-(psi / b)**c),
  !
  !   psi_fl = -b (-LOG(r))**(1 / c).
  !
  ! psi_floor_leaf and psi_floor_stem come from the leaf and stem segment curves
  ! (Weibull), which the store's balance and the optimiser use, in spin-up
  ! and without an impairment model too. The leak is cut where the leaf
  ! would end below psi_floor_leaf (store_cpl_leak_limit), and the store's stem
  ! psi is kept at or above psi_floor_stem (store_cpl_q).
  !
  ! MDK, 05/11/26
  r_crit = 1.0 - kcrit_fractional_loss(ft)
  psi_floor_leaf = -ABS(conductance_b_seg(ft,3))                               &
                   * (-LOG(r_crit))**(1.0 / conductance_c_seg(ft,3))
  psi_floor_stem = -ABS(conductance_b_seg(ft,2))                               &
                   * (-LOG(r_crit))**(1.0 / conductance_c_seg(ft,2))

  ! Advance the capacitance store over the step. Its state is its water,
  !
  !   W = C psi_x      (mol H2O m-2 ground, <= 0; psi_x its psi, Pa),
  !
  ! carried between steps (and dumped). Wherever the store is set, day and
  ! night, its balance over the step at the final canopy E (store_cpl_q,
  ! exact):
  !
  !   q = C (psi_x0 - psi_x) / dt = E - E_r,   psi_x = psi_st(E_r),
  !
  ! the stem node psi_st where the soil, root and stem carry E_r, and the
  ! leaf psi_leaf where the leaf segment then carries E. With the leak, the
  ! leak's E is cut where the leaf would end below psi_floor_leaf
  ! (store_cpl_leak_limit). q goes to store_q (kg m-2 s-1), which sf_evap
  ! takes off the soil extraction, so the water is conserved. With almost
  ! no leaves (LAI < lai_store_min) the store is dropped (unset), and the
  ! water it has given and not yet taken back (its deficit, store_w
  ! <= 0) is drawn from the soil over the step, as a refill, so a store set
  ! again later starts full without making water; a store not yet set
  ! (first step, or a dump without it) starts at the step's stem psi.
  DO m = 1,veg_pts
    l = veg_index(m)
    e_store_bal(l) = MAX(el(l), 0.0)
  END DO
  CALL store_cpl_q( .TRUE., e_store_bal )
  ! The passes' fixed point: the store's flow at the committed E against
  ! the flow the last pass held fixed, relative to E (store_q_gap).
  !
  ! MDK, 05/11/26
  DO m = 1,veg_pts
    l = veg_index(m)
    store_q_gap(l,ft) = 0.0
    IF ( e_store_bal(l) > 0.0 )                                                &
      store_q_gap(l,ft) = ABS(q_store_bal(l) - q_store(l)) / e_store_bal(l)
  END DO
  IF ( l_som_cuticular_loss ) CALL store_cpl_leak_limit()
  DO m = 1,veg_pts
    l = veg_index(m)
    ! The store's capacitance (mol m-2 Pa-1).
    c_store = store_c(l,ft)

    ! Almost no leaves: the store is dropped, its deficit refilled from
    ! the soil over the step.
    IF ( lai(l) < lai_store_min ) THEN
      store_q(l,ft)  = 0.0
      IF ( store_w(l,ft) <= 0.0 )                                              &
        store_q(l,ft) = store_w(l,ft) * Water_kg_per_mol / timestep
      store_w(l,ft)  = store_unset
      CYCLE
    END IF
    ! The store set: its balance at the final E, and kplant without the
    ! store's shift, per unit leaf area as leaf_k.
    IF ( store_w(l,ft) <= 0.0 ) THEN
      ! Two-leaf canopy: each class's leaf psi from the committed stem psi
      ! (store_cpl_cls_psi), where the store's balance was solved.
      !
      ! MDK, 04/11/26
      IF ( can_rad_mod == 7 .AND. c_store > 0.0 ) CALL store_cpl_cls_psi( l )
      psi0_store = store_w(l,ft) / c_store
      CALL store_commit( l, c_store, psi0_store, psi_stem_bal(l) )
      psi_leaf(l) = psi_leaf_bal(l)
      psi_stem(l) = psi_stem_bal(l)
      IF ( kplant_bal(l) > 0.0 )                                               &
        leaf_k(l) = kplant_bal(l) / MAX(lai(l), TINY(1.0_real_jlslsm))
      CYCLE
    END IF
    ! The store not yet set: it starts at the step's stem psi.
    store_w(l,ft) = c_store * MIN(psi_stem(l), 0.0)
  END DO

  ! The leak's transpiration (leak_e, kg m-2 s-1): the final canopy flux
  ! where the leak sets the conductance, after its cut.
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( l_leak_on(l) ) leak_e(l,ft) = MAX(el(l), 0.0) * Water_kg_per_mol
  END DO

  ! Xylem impairment update, once per timestep, in place of those in the
  ! big-leaf and two-leaf paths: the same lit points, with the water
  ! potentials of the total (stomatal plus cuticular) flux. Also the unlit
  ! points where the leak acts, so the night leak's psi drives the loss of
  ! conductance.
  l_in_impair(:) = .FALSE.
  DO m = 1,n_lit
    l_in_impair(open_land_index(m)) = .TRUE.
  END DO
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( l_leak_on(l) .AND. .NOT. l_in_impair(l) ) THEN
      n_lit = n_lit + 1
      open_land_index(n_lit) = l
    END IF
  END DO
  ! The two-leaf canopy also gives its classes' leaf psi, for the leaf
  ! segment's damage.
  !
  ! MDK, 04/11/26
  IF ( can_rad_mod == 7 ) THEN
    CALL update_xylem_impairment(                                              &
          ! IN
            land_pts, n_lit, open_land_index, ft, psi_leaf, psi_src,           &
            psi_stem,                                                          &
          ! IN OUT
            k_max, psi_leaf_extreme, psi_root_extreme,                         &
          ! IN (recovery)
            lai, canht, anetc,                                                 &
          ! IN (leaf classes)
            psi_leaf_cls = psi_leaf_cls, w_cls = w_cls )
  ELSE
    CALL update_xylem_impairment(                                              &
          ! IN
            land_pts, n_lit, open_land_index, ft, psi_leaf, psi_src,           &
            psi_stem,                                                          &
          ! IN OUT
            k_max, psi_leaf_extreme, psi_root_extreme,                         &
          ! IN (recovery)
            lai, canht, anetc )
  END IF

  ! Water released by cavitation (l_ximpair_cav_water). The water of the
  ! stem's conduits is part of the plant's water, W_xyl (store_w_xyl), and
  ! each step it is brought to the conduits not lost,
  !
  !   W_xyl = V_w (1 - L),
  !   V_w = rho_w (X_s / (c_w rho)) f   (kg m-2),
  !
  ! with L the stem's locked-in loss, X_s the sapwood carbon (nsc_sapwood),
  ! rho the wood density (wood_density), c_w = 0.46 the carbon fraction of
  ! wood dry mass (the mean of live stem wood of 31 Quercus species,
  ! Doraisami et al. 2022, Sci Data 9: 284) and f the share of the
  ! sapwood's volume in the conduits' lumens (xyl_lumen_frac); unset, f is
  ! all the wood's lumens, 1 - rho / rho_cw, with rho_cw = 1500 kg m-3 the
  ! density of the wood substance in water (hardwoods, Kellogg & Wangaard
  ! 1969, Wood Fiber 1: 180). The difference is exchanged with the plant
  ! store, so the water in soil, store and conduits adds up:
  !
  !   - conduits that embolise (L rises, or sapwood is lost) give their
  !     water to the store, up to full (psi = 0); what the store cannot
  !     take stays in W_xyl until it has room;
  !   - new xylem (a new ring, which also lowers L, or more sapwood)
  !     starts with its conduits full, and that water is drawn from the
  !     soil over the step (store_q_xyl, which sf_evap takes from the soil
  !     with the store's refill).
  !
  ! MDK, 12/11/26
  IF ( l_ximpair_cav_water .AND. ALLOCATED(ximpair_lock) .AND.                 &
       ALLOCATED(nsc_sapwood) ) THEN
    ! The water drawn for new xylem is this step's only.
    store_q_xyl(:,ft) = 0.0

    IF ( wood_density(ft) > 0.0 ) THEN
      DO m = 1,veg_pts
        l = veg_index(m)

        ! Nothing is exchanged while the store is not set or before
        ! TRIFFID has set the sapwood; the conduits' water waits.
        IF ( store_w(l,ft) > 0.0 .OR. nsc_sapwood(l,ft) < 0.0 ) CYCLE

        ! The water of the sapwood's conduits when all are full, V_w
        ! (kg m-2): the sapwood's volume, X_s / (c_w rho), times the share
        ! of it in the conduits' lumens, f, times rho_w,
        !
        !   V_w = rho_w (X_s / (c_w rho)) f,
        !
        ! with f = xyl_lumen_frac, or unset all the wood's lumens,
        ! f = 1 - rho / rho_cw.
        IF ( xyl_lumen_frac(ft) >= 0.0 ) THEN
          w_xyl_full = rho_water * nsc_sapwood(l,ft)                           &
                       / (c_w_cav * wood_density(ft)) * xyl_lumen_frac(ft)
        ELSE
          w_xyl_full = rho_water * nsc_sapwood(l,ft) / c_w_cav                 &
                       * (1.0 / wood_density(ft) - 1.0 / rho_cw_cav)
        END IF

        ! The water the conduits should hold now: the share not lost,
        ! 1 - L.
        w_xyl_hold = w_xyl_full * (1.0 - MIN(MAX(ximpair_lock(l,ft), 0.0), 1.0))

        ! First time here: the conduits not lost start full. This is an
        ! initial condition, so nothing is exchanged.
        IF ( store_w_xyl(l,ft) < 0.0 ) THEN
          store_w_xyl(l,ft) = w_xyl_hold
          CYCLE
        END IF

        ! The conduits' water above what they should now hold (kg m-2):
        ! > 0 when conduits have embolised (or the sapwood has shrunk)
        ! since the last step, plus any water held back while the store
        ! was full; < 0 when new xylem has formed.
        w_xyl_excess = store_w_xyl(l,ft) - w_xyl_hold

        IF ( w_xyl_excess > 0.0 ) THEN
          ! Embolism: the water goes to the store, at most up to full
          ! (store_w = 0); the rest stays in the conduits' water and goes
          ! when the store has room.
          w_to_store = MIN(w_xyl_excess, -store_w(l,ft) * Water_kg_per_mol)
          store_w(l,ft) = store_w(l,ft) + w_to_store / Water_kg_per_mol
          store_w_xyl(l,ft) = store_w_xyl(l,ft) - w_to_store
        ELSE
          ! New xylem starts with its conduits full: that water is drawn
          ! from the soil over this step (sf_evap); what the soil cannot
          ! give comes from the store (store_return).
          store_q_xyl(l,ft) = -w_xyl_excess / timestep
          store_w_xyl(l,ft) = w_xyl_hold
        END IF

      END DO
    END IF
  END IF
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
    ! The balanced LAI from the canopy height (lai_bal_from_ht).
    !
    ! MDK, 20/10/26
    lai_bal(l) = lai_bal_from_ht(ft, ht(l))

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

        !Reduce n_stem for consistency with non-trait n_stem for now
        !This must be done to achieve realistic respiration rates.
        fstem            = 1.0 / a_ws(ft)
        stem_resp_scaling = fstem + (1.0 - fstem) * hw_sw(ft)
        n_stem(l)         = n_stem(l) * stem_resp_scaling
      ELSE
        ! The respiring stem is the live sapwood, with nsw its N per unit
        ! carbon (kg N / kg C): TRIFFID's live stemwood, or the sapwood pool
        ! with l_nsc_sapwood_n once the pool is set. hw_sw (heartwood N)
        ! does not respire.
        ! MDK, 04/11/26
        n_stem(l) = live_stemwood(ft, ht(l), lai_bal(l)) * nsw(ft)
        IF ( l_nsc_sapwood_n ) THEN
          IF ( ALLOCATED(nsc_sapwood) ) THEN
            IF ( nsc_sapwood(l,ft) >= 0.0 )                                    &
              n_stem(l) = nsc_sapwood(l,ft) * nsw(ft)
          END IF
        END IF
      END IF

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
  ! instead), so rd must be added back unscaled. Using fsmc for both would
  ! leave stom_opt's GPP short by rdc * (1 - fsmc) - and negative at night,
  ! where closed-stomata anetc = -rdc.
  IF (leaf_flux_mod == leaf_flux_stom_opt) THEN
    fsmc_leaf_resp(l) = 1.0
  ELSE
    fsmc_leaf_resp(l) = fsmc(l)
  END IF

  gpp(l) = cconu * (anetc(l) + rdc(l) * fsmc_leaf_resp(l))
  ! Net canopy photosynthesis for the photosynthesis diagnostics
  ! (l_photo_diag).
  !
  ! MDK, 30/10/26
  IF ( l_photo_diag ) anetc_pft(l,ft) = anetc(l)

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
!
! MDK, 18/10/26
IF (l_triffid_nsc)                                                             &
  CALL triffid_nsc_turgor( land_pts, veg_pts, veg_index, ft, timestep,         &
                           psi_root_zone, tair )

! Leaf shedding (leaf_shed_mod): the daily drought counts (l_phen_psi),
! the plant's daily water status (l_phen_fuse) and the leaf heat damage
! (l_phen_heat).
!
! MDK, 09/11/26
IF (l_phen_psi)                                                                &
  CALL phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,            &
                        psi_root_zone )
IF (l_phen_fuse)                                                               &
  CALL phen_fuse_update( land_pts, veg_pts, veg_index, ft, timestep,           &
                         psi_root_zone, psi_leaf )
IF (l_phen_heat)                                                               &
  CALL heat_damage_update( land_pts, veg_pts, veg_index, ft, timestep,         &
                           tstar, t_can, t_leaf_sun )

! TRIFFID store without TRIFFID (prescribed LAI): daily store step.
!
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

! The gravity drop at the root-stem node is this timestep's only.
!
! MDK, 04/11/26
CALL clear_grav_drop()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

!-----------------------------------------------------------------------------
! Plant capacitance: save the store's flow and water at the end of the step.
!
! The store is water held in the stem. It provides the water for transpiration
! as the stem dries and refills from the soil as the stem wets. Over the step
! the store's water potential goes from psi_x0 to psi_x1, and C is the water
! it gives up per unit fall in psi.
!
! store_q is the rate at which the store gives water over the step
! (> 0: it gave water to the transpiration; < 0: it refilled from the soil):
!
!   store_q = C (psi_x0 - psi_x1) M_w / dt      (kg H2O m-2 s-1)
!
! store_w is the store's water at the end of the step, counted from a
! full store at psi = 0, so it is 0 when full and negative as the stem dries:
!
!   store_w = C psi_x1                           (mol H2O m-2)
!
! M_w is the molar mass of water and dt the timestep. Both are per m2 of
! ground. The "cut_" prefix is from the module that holds the store,
! plant_water_store_mod.
!
! MDK, 05/11/26
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

store_q(l_c,ft) = c_s * (psi_x0 - psi_x1) / timestep * Water_kg_per_mol
store_w(l_c,ft) = c_s * psi_x1

END SUBROUTINE store_commit

!-----------------------------------------------------------------------------
! The plant capacitance store coupled to the optimiser: the store's flow
! out (q_store, mol H2O m-2 s-1) that the first pass of the step sees, in all
! three canopies. It is last step's (0 where the store is not set). Where
! the canopy was shut on the last step (no stomatal transpiration, e.g. the
! night before the first lit step) that flow is a refill (< 0), floored at
! 0 as between the passes, so leaves that open now are not searched with
! their root flow raised by last night's refill.
!
! store_q and store_e_prev are in the dump, so a restart in daylight
! starts its first pass from the same flow as a continuous run.
!
! MDK, 05/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_first_q()

IMPLICIT NONE

INTEGER :: l, m

!$OMP CRITICAL (som_store_e_prev_alloc)
IF ( .NOT. ALLOCATED(store_e_prev) ) THEN
  ALLOCATE( store_e_prev(land_pts, npft) )
  store_e_prev(:,:) = 0.0
END IF
!$OMP END CRITICAL (som_store_e_prev_alloc)

q_store(:) = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( store_w(l,ft) <= 0.0 .AND. lai(l) >= lai_store_min ) THEN
    q_store(l) = store_q(l,ft) / Water_kg_per_mol
    IF ( store_e_prev(l,ft) <= 0.0 ) q_store(l) = MAX(q_store(l), 0.0)
  END IF
END DO

END SUBROUTINE store_cpl_first_q

!-----------------------------------------------------------------------------
! Keep the canopy's stomatal transpiration e_c (mol H2O m-2 s-1, ground) for
! the next step's first pass (store_cpl_first_q).
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_e_prev( e_c )

IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) :: e_c(land_pts)

INTEGER :: l, m

DO m = 1,veg_pts
  l = veg_index(m)
  store_e_prev(l,ft) = MAX(e_c(l), 0.0)
END DO

END SUBROUTINE store_cpl_e_prev

!-----------------------------------------------------------------------------
! Plant capacitance: the most the canopy may transpire this step is what
! the soil can supply plus what the stem store gives.
!
!   e_sup_use = e_supply + q_store M_w      (kg H2O m-2 s-1, at least 0)
!
! e_supply is the soil's supply (< 0: no cap, which stays no cap), q_store the
! store's flow out (mol H2O m-2 s-1; < 0 when it refills, which lowers the
! cap) and M_w the molar mass of water.
!
! Sets e_sup_use.
!
! MDK, 05/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_supply()

IMPLICIT NONE

INTEGER :: m, l

DO m = 1,veg_pts
  l = veg_index(m)
  e_sup_use(l) = e_supply(l)
  IF ( e_supply(l) >= 0.0 )                                                    &
    e_sup_use(l) = MAX(e_supply(l) + q_store(l) * Water_kg_per_mol, 0.0)
END DO

END SUBROUTINE store_cpl_supply

!-----------------------------------------------------------------------------
! Plant capacitance, two-leaf canopy: split the store's flow q_store between
! the sunlit and shaded leaves (q_store_sun, q_store_shd).
!
! Each leaf gets the share of q_store that it transpired on the last pass
! (el_sun, el_shd). On the first pass (l_pass false), or where the canopy
! did not transpire, each gets its share of the canopy's hydraulic
! conductance (kmax_ref_sun_2l, kmax_ref_shd_2l) instead.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_split( l_pass )

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_pass

INTEGER :: m, l
REAL(KIND=real_jlslsm) :: s_sun
                            ! The sunlit leaves' share of q_store (-).
REAL(KIND=real_jlslsm) :: e_t
                            ! The canopy's transpiration on the last pass
                            ! (mol H2O m-2 ground s-1).

DO m = 1,veg_pts
  l = veg_index(m)

  ! Sunlit leaves' share of the canopy's intact conductance, used on the
  ! first pass or where the canopy did not transpire.
  s_sun = kmax_ref_sun_2l(l) /                                                 &
          MAX(kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l), TINY(1.0_real_jlslsm))

  ! After the first pass: the sunlit leaves' share of the canopy's
  ! transpiration on the last pass instead.
  e_t   = MAX(el_sun(l), 0.0) + MAX(el_shd(l), 0.0)
  IF ( l_pass .AND. e_t > 0.0 ) s_sun = MAX(el_sun(l), 0.0) / e_t

  ! The shaded leaves take the rest, so the two parts add to q_store.
  q_store_sun(l) = s_sun * q_store(l)
  q_store_shd(l) = q_store(l) - q_store_sun(l)
END DO

END SUBROUTINE store_cpl_split

!-----------------------------------------------------------------------------
! Plant capacitance, multilayer canopy: the store's flow q_store given to the
! sunlit and shaded leaves of layer n_l, per unit leaf area (q_store_sun_lyr,
! q_store_shd_lyr), on pass k_p.
!
! Each leaf gets q_store times its transpiration on the last pass (el_sun_lyr,
! el_shd_lyr) over the canopy's (e_ml_prev). On the first pass, or where
! the canopy did not transpire, both get q_store times the layer's share of the
! canopy's hydraulic conductance (share_sup) instead.
!
! MDK, 05/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_ml_q( n_l, k_p )

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_l, k_p
                            ! The layer and the pass.

INTEGER :: m, l

DO m = 1,veg_pts
  l = veg_index(m)

  ! After the first pass: each leaf gets q_store in proportion to its
  ! transpiration per leaf area on the last pass over the canopy's.
  IF ( k_p > 1 .AND. e_ml_prev(l) > 0.0 ) THEN
    q_store_sun_lyr(l) = q_store(l) * MAX(el_sun_lyr(l,n_l), 0.0) / e_ml_prev(l)
    q_store_shd_lyr(l) = q_store(l) * MAX(el_shd_lyr(l,n_l), 0.0) / e_ml_prev(l)

  ! First pass, or the canopy did not transpire: by the layer's share of
  ! the canopy's conductance per leaf area, the same for sunlit and shaded
  ! leaves, as they have the same kmax per leaf area.
  ELSE
    q_store_sun_lyr(l) = q_store(l) * share_sup(l)
    q_store_shd_lyr(l) = q_store(l) * share_sup(l)
  END IF
END DO

END SUBROUTINE store_cpl_ml_q

!-----------------------------------------------------------------------------
! Plant capacitance with the cuticular leak (l_som_cuticular_loss): stop
! the leak drying the leaf past its critical point.
!
! Where the stomata are shut, water still leaks out through the cuticle
! (l_leak_on). If that leak would take the leaf's water potential below
! psi_floor_leaf, where the leaf has lost its critical share of conductance,
! the leak is reduced to the largest transpiration that keeps the leaf at
! psi_floor_leaf at the end of the step.
!
! That transpiration is found by bisection between 0 and the leak's own
! value, solving the store's balance (store_cpl_q) at each guess. 20 steps
! find it to 1e-6 of the leak's value. The canopy conductance gc is then
! set to give that transpiration through the tile's flux equation, and the
! store's balance is committed. Only these points are solved again; the
! others keep the balance already committed.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_leak_limit()

IMPLICIT NONE

INTEGER, PARAMETER :: n_bisect_leak = 20
                            ! Bisection steps on the leak's E.
INTEGER :: m, l, it
LOGICAL :: l_leak_limited(land_pts)
                            ! Points where the leak's E is cut.
REAL(KIND=real_jlslsm) :: e_lo(land_pts)
                            ! Lower end of the bisection bracket on E (mol
                            ! H2O m-2 s-1).
REAL(KIND=real_jlslsm) :: e_hi(land_pts)
                            ! Upper end of that bracket (mol H2O m-2 s-1).
REAL(KIND=real_jlslsm) :: e_per_gc_cut
                            ! rho (qs(tstar) - q1) (kg m-3), for the exact gc.
REAL(KIND=real_jlslsm) :: e_kg_cut
                            ! The cut E (kg m-2 s-1), for the exact gc.

! The points to cut: the leak sets the conductance (l_leak_on), the committed
! leaf psi fell below psi_floor_leaf (l_leaf_floor_hit), the canopy
! transpires, and the store is set up there. The bisection bracket on E runs
! from 0 to the committed E (the leak's own value).
l_leak_limited(:) = .FALSE.
DO m = 1,veg_pts
  l = veg_index(m)
  l_leak_limited(l) = l_leak_on(l) .AND. l_leaf_floor_hit(l)                   &
                      .AND. el(l) > 0.0 .AND. store_w(l,ft) <= 0.0             &
                      .AND. lai(l) >= lai_store_min
  e_lo(l) = 0.0
  e_hi(l) = e_store_bal(l)
END DO
IF ( .NOT. ANY(l_leak_limited) ) RETURN

! Bisection on E. At each guess, solve the store's balance (with the
! search's solver: Newton steps, or the lookup tables with som_psi_solver =
! 3; cut points only) for the leaf psi at the end of the step. A leaf below
! psi_floor_leaf means E is too large, so the guess becomes the upper end;
! otherwise it becomes the lower end. e_lo is always an E the leaf survives.
DO it = 1,n_bisect_leak
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( l_leak_limited(l) ) e_store_bal(l) = 0.5 * (e_lo(l) + e_hi(l))
  END DO
  CALL store_cpl_q( .FALSE., e_store_bal, l_sel = l_leak_limited )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( .NOT. l_leak_limited(l) ) CYCLE
    IF ( psi_leaf_bal(l) < psi_floor_leaf ) THEN
      e_hi(l) = e_store_bal(l)
    ELSE
      e_lo(l) = e_store_bal(l)
    END IF
  END DO
END DO

! Take E = e_lo, and find the canopy conductance gc that gives it. The
! tile's flux is not linear in gc,
!
!   E M = A gc / (1 + ra (gc + g_s)),   A = rho (qs(tstar) - q1),
!
! with M the molar mass of water, ra the aerodynamic resistance and g_s the
! soil's conductance, so
!
!   gc = E M (1 + ra g_s) / (A - ra E M).
!
! The denominator is positive, because the cut E is below the leak's E
! under the same air.
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( .NOT. l_leak_limited(l) ) CYCLE
  e_store_bal(l) = e_lo(l)
  e_per_gc_cut   = pstar(l) / ( r * tstar(l) ) * MAX(qs(l) - q1(l), 0.0)
  e_kg_cut  = e_lo(l) * m_h2o_lt

  ! If the denominator is near 0 (rounding), scale gc by the cut instead.
  IF ( e_per_gc_cut - ra(l) * e_kg_cut > TINY(1.0_real_jlslsm) ) THEN
    gc(l) = e_kg_cut * ( 1.0 + ra(l) * g_soil_c(l) )                           &
            / ( e_per_gc_cut - ra(l) * e_kg_cut )
  ELSE
    gc(l) = gc(l) * e_lo(l) / el(l)
  END IF

  ! No lower than the closed canopy's conductance, min_gl_pft (a canopy
  ! value), when even the source psi is below the floor (e_lo = 0).
  gc(l) = MAX(gc(l), min_gl_pft(ft))
  gc_som(l)  = gc(l)
  el(l)      = e_lo(l)
END DO

! Commit the store's balance at the cut E, with the exact segment solves.
CALL store_cpl_q( .TRUE., e_store_bal, l_sel = l_leak_limited )

END SUBROUTINE store_cpl_leak_limit

!-----------------------------------------------------------------------------
! Plant capacitance: how much of the canopy's transpiration the stem store
! gives this step, and how much comes from the soil.
!
! The canopy transpires E. Part of it, q, comes from the stem store and the
! rest, E - q, is pulled from the soil up through the roots and stem. The
! two must agree:
!
!   - the store gives water as it dries, q = C (psi_x0 - psi_st) / dt,
!     where psi_x0 is its water potential at the start of the step and
!     psi_st at the end;
!   - psi_st is also the stem's water potential when the roots and stem
!     carry E - q from the soil. The more the store gives, the less the
!     soil must supply, so the wetter the stem.
!
! The q where both hold is found by a bracketed search between two limits.
! The lower one is the store's flow if the stem were as wet as the soil
! (negative when the store refills). The upper one is all of E, as the soil
! never takes water back from the plant (E - q >= 0).
!
! If the store is so much wetter than the soil that it could give all of E,
! it does (the soil gives nothing), and the leaf's water potential is
! worked out from the store's.
!
! On the final call of the step (l_ex), the search uses the exact plant
! hydraulics, the stem and leaf are kept no drier than their critical
! points, the store gives no more than E, and the plant's conductance
! (kplant) is recorded as if there were no store. Before that, only points
! that transpire are solved, and with l_sel only the chosen points.
!
! Sets q_store_bal (the store's flow, as store_commit takes it),
! psi_stem_bal (the stem), psi_leaf_bal (the leaf), kplant_bal (kplant) and
! l_leaf_floor_hit (the leaf passed its critical point). Points without a
! store keep the soil-only values. It also sets share_sup to 1; each canopy
! scheme sets it again before using it.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_q( l_ex, e_c, l_sel )

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_ex
                            ! .TRUE.: exact segment solves (the commit);
                            ! .FALSE.: the search's solver (between
                            ! passes): Newton steps, or the lookup tables
                            ! with som_psi_solver = 3.
REAL(KIND=real_jlslsm), INTENT(IN) :: e_c(land_pts)
                            ! Canopy transpiration (mol H2O m-2 s-1).
LOGICAL, INTENT(IN), OPTIONAL :: l_sel(land_pts)
                            ! Only these points are solved; the outputs at
                            ! the others are left as they were (the leak
                            ! cut).

INTEGER, PARAMETER :: n_it_bal = 40
                            ! Most regula falsi steps.
REAL(KIND=real_jlslsm), PARAMETER :: rtol_bal = 1.0e-6
                            ! Tolerance on the bracket of q, relative to
                            ! its first width.

INTEGER :: m, l, j, it, n_bal, n_only
INTEGER :: idx_bal(land_pts)
                            ! Positions in veg_index of the points solved.
INTEGER :: idx_only(land_pts)
                            ! Positions in veg_index of the points where the
                            ! store gives all of E.
LOGICAL :: l_imp_bal
                            ! Xylem impairment active.
LOGICAL :: l_all_done
                            ! Every point solved has converged.
LOGICAL :: l_only(land_pts)
                            ! The store pays all of E.
LOGICAL :: l_bis(land_pts)
                            ! The next step bisects (Illinois fallback).
LOGICAL :: l_done(land_pts)
                            ! Converged.
INTEGER :: side_bal(land_pts)
                            ! Illinois: the end replaced last (-1 lo, 1 hi).
REAL(KIND=real_jlslsm) :: c_p(land_pts)
                            ! C / dt (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: psi_x0(land_pts)
                            ! The store's psi at the start of the step (Pa).
REAL(KIND=real_jlslsm) :: q_lo(land_pts)
                            ! Lower end of the bracket on the store's flow q
                            ! (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: q_hi(land_pts)
                            ! Upper end of that bracket (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: f_lo(land_pts)
                            ! The mismatch R at q_lo (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: f_hi(land_pts)
                            ! The mismatch R at q_hi (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: q_try(land_pts)
                            ! The guessed q (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: f_try(land_pts)
                            ! The mismatch R at q_try (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: w_tol(land_pts)
                            ! The bracket width at which q has converged
                            ! (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: w_old(land_pts)
                            ! The bracket width at the last step (mol m-2
                            ! s-1).
REAL(KIND=real_jlslsm) :: kref_c(land_pts)
                            ! The canopy's intact kmax (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: kmx_c(land_pts)
                            ! The canopy's impaired kmax (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: kcr_c(land_pts)
                            ! The canopy's kcrit (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: kcap_c(land_pts)
                            ! kmx_c / kref_c, the impaired fraction (-).
REAL(KIND=real_jlslsm) :: psi1(1,land_pts)
                            ! One sample per point: the leaf psi of the chain
                            ! (Pa).
REAL(KIND=real_jlslsm) :: st1(1,land_pts)
                            ! The stem psi of the chain (Pa).
REAL(KIND=real_jlslsm) :: k1(1,land_pts)
                            ! The chain's conductance, with l_ex (mol m-2 s-1
                            ! Pa-1).
REAL(KIND=real_jlslsm) :: psi_from(land_pts)
                            ! Where the chain starts (Pa): psi_src, or the
                            ! store's psi where it pays all of E.

! Start every output at its no-store value: no flow from the store, the
! stem at the soil's water potential and the leaf as the canopy scheme left
! it. With l_sel only the chosen points are reset; the others keep what an
! earlier call gave them.
IF ( PRESENT(l_sel) ) THEN
  WHERE ( l_sel(:) )
    q_store_bal(:)      = 0.0
    psi_stem_bal(:)     = psi_zf(:)
    psi_leaf_bal(:)     = psi_leaf(:)
    kplant_bal(:)       = 0.0
    l_leaf_floor_hit(:) = .FALSE.
  END WHERE
ELSE
  q_store_bal(:)      = 0.0
  psi_stem_bal(:)     = psi_zf(:)
  psi_leaf_bal(:)     = psi_leaf(:)
  kplant_bal(:)       = 0.0
  l_leaf_floor_hit(:) = .FALSE.
END IF

! The guessed store flow at every point. q_try is passed whole to
! set_store_q_path, so the points not solved must hold 0, not junk.
q_try(:)     = 0.0

! The plant's hydraulics for the whole canopy, as the canopy scheme set
! them: intact conductance (kmax_ref_store), conductance with today's damage
! (kmax_store), the critical conductance (kcrit_store) and the damaged share of
! the intact one (kcap_c). The soil-to-root link is the whole canopy's
! (share_sup = 1). l_imp_bal: xylem damage is on, so the damaged path is
! used.
l_imp_bal = .NOT. is_spinup .AND.                                              &
           pft_xylem_impairment_model(ft) /= xylem_impairment_none
DO m = 1,veg_pts
  l = veg_index(m)
  kref_c(l) = kmax_ref_store(l)
  kmx_c(l)  = kmax_store(l)
  kcr_c(l)  = kcrit_store(l)
  kcap_c(l) = kmx_c(l) / MAX(kref_c(l), TINY(1.0_real_jlslsm))
  share_sup(l) = 1.0
END DO
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

! Choose the points to solve and set up the search at each. The search
! for the store's flow q runs between q_lo, the flow if the stem were as
! wet as the soil, and q_hi = E, the most the store may give.
n_bal = 0
DO m = 1,veg_pts
  l = veg_index(m)
  l_only(l) = .FALSE.
  ! (No capacity, as from a dump without store_c: the store's psi
  ! W / C is undefined, so it is taken as not set.)
  IF ( store_w(l,ft) > 0.0 .OR. store_c(l,ft) <= 0.0 .OR.                      &
       lai(l) < lai_store_min ) CYCLE
  ! Between passes (not l_ex) only where the point transpires: q_store is used
  ! only by the next pass's open points.
  IF ( .NOT. l_ex .AND. e_c(l) <= 0.0 ) CYCLE
  IF ( PRESENT(l_sel) ) THEN
    IF ( .NOT. l_sel(l) ) CYCLE
  END IF

  ! The store's capacitance per second, and its water potential at the
  ! start of the step from its water, psi = W / C.
  c_p(l)    = store_c(l,ft) / timestep
  psi_x0(l) = store_w(l,ft) / store_c(l,ft)

  ! The two ends of the search.
  q_lo(l)   = c_p(l) * (psi_x0(l) - psi_zf(l))
  q_hi(l)   = e_c(l)
  n_bal = n_bal + 1
  idx_bal(n_bal) = m

  ! The store is so wet that even with the stem as wet as the soil it would
  ! give more than E: it gives all of E, the soil gives nothing, and the
  ! store dries by E over the step. No search is needed.
  IF ( q_lo(l) >= q_hi(l) ) THEN
    l_only(l)       = .TRUE.
    q_store_bal(l)  = e_c(l)
    psi_stem_bal(l) = psi_x0(l) - e_c(l) / c_p(l)
  END IF

  ! Search bookkeeping: no forced halving yet, which end moved last (none),
  ! and the stopping width, a millionth of the starting one.
  l_bis(l)   = .FALSE.
  l_done(l)  = l_only(l)
  side_bal(l) = 0
  w_tol(l)   = rtol_bal * MAX(q_hi(l) - q_lo(l), TINY(1.0_real_jlslsm))
  w_old(l)   = q_hi(l) - q_lo(l)
END DO
IF ( n_bal == 0 ) RETURN

! The mismatch at both ends of the search: the guessed q less the flow the
! store would give at the stem's resulting water potential. It is below 0
! at q_lo and above 0 at q_hi, so the answer lies between.
DO j = 1,n_bal
  q_try(veg_index(idx_bal(j))) = q_lo(veg_index(idx_bal(j)))
END DO
CALL store_cpl_eval( n_bal, idx_bal, l_imp_bal, .FALSE., e_c, q_try, c_p,      &
                     psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_bal
  l = veg_index(idx_bal(j))
  f_lo(l) = f_try(l)
  q_try(l) = q_hi(l)
END DO
CALL store_cpl_eval( n_bal, idx_bal, l_imp_bal, .FALSE., e_c, q_try, c_p,      &
                     psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_bal
  l = veg_index(idx_bal(j))
  f_hi(l) = f_try(l)
END DO

! Narrow the range until it holds the answer to within the stopping width.
! Each new guess is where a straight line drawn between the mismatches at
! the two ends crosses 0. If the last guess shrank the range by less than
! half, the midpoint is used instead, which always halves it. (This is the
! "regula falsi" method, in its "Illinois" form.) Points already done sit
! at their midpoint, so later calls see a settled value.
DO it = 1,n_it_bal
  l_all_done = .TRUE.
  DO j = 1,n_bal
    l = veg_index(idx_bal(j))
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
  CALL store_cpl_eval( n_bal, idx_bal, l_imp_bal, .FALSE., e_c, q_try, c_p,    &
                       psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )

  ! The guess replaces the end whose mismatch has the same sign, so the
  ! answer stays inside the range. If the same end has moved twice running,
  ! the other end is stuck: its mismatch is halved, which pulls the next
  ! guess toward it. Stop once the range is below the stopping width.
  DO j = 1,n_bal
    l = veg_index(idx_bal(j))
    IF ( l_done(l) ) CYCLE
    IF ( f_try(l) < 0.0 ) THEN
      q_lo(l) = q_try(l)
      f_lo(l) = f_try(l)
      IF ( side_bal(l) == -1 ) f_hi(l) = 0.5 * f_hi(l)
      side_bal(l) = -1
    ELSE
      q_hi(l) = q_try(l)
      f_hi(l) = f_try(l)
      IF ( side_bal(l) == 1 ) f_lo(l) = 0.5 * f_lo(l)
      side_bal(l) = 1
    END IF
    l_bis(l)  = q_hi(l) - q_lo(l) > 0.5 * w_old(l)
    w_old(l)  = q_hi(l) - q_lo(l)
    l_done(l) = w_old(l) <= w_tol(l)
  END DO
END DO

! The answer: the middle of the final range (q = E where the store gives
! all of E). Solve the plant once more there, exactly on the final call
! (l_ex), and take the stem and leaf water potentials from it. q_store_bal is
! worked out again from the stem, so the store's water balances exactly.
DO j = 1,n_bal
  l = veg_index(idx_bal(j))
  q_try(l) = 0.5 * (q_lo(l) + q_hi(l))
  IF ( l_only(l) ) q_try(l) = e_c(l)
END DO
CALL store_cpl_eval( n_bal, idx_bal, l_imp_bal, l_ex, e_c, q_try, c_p, psi_x0, &
                     kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try )
DO j = 1,n_bal
  l = veg_index(idx_bal(j))
  IF ( .NOT. l_only(l) ) THEN
    psi_stem_bal(l) = st1(1,j)
    psi_leaf_bal(l) = psi1(1,j)
    q_store_bal(l)  = c_p(l) * (psi_x0(l) - psi_stem_bal(l))
  END IF
END DO

! Where the store gives all of E, the stem is at the store's water
! potential at the end of the step, and the leaf is solved from there:
! the leaf carries E from the stem, and the roots and stem carry nothing.
! Starting from the soil's water potential instead would, in a drought,
! give a far wetter leaf that conducts too well.
! (The root carries nothing, so the chain starts at the stem psi plus the
! gravity drop the root-stem node takes off.)
! Only those points are solved again.
n_only = 0
DO j = 1,n_bal
  IF ( l_only(veg_index(idx_bal(j))) ) THEN
    n_only = n_only + 1
    idx_only(n_only) = idx_bal(j)
  END IF
END DO
IF ( n_only > 0 ) THEN
  psi_from(:) = psi_src(:)
  DO j = 1,n_only
    l = veg_index(idx_only(j))
    psi_from(l) = psi_stem_bal(l) + grav_node_drop(l)
  END DO
  CALL store_cpl_eval( n_only, idx_only, l_imp_bal, l_ex, e_c, q_try, c_p,     &
                       psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try, &
                       psi_from = psi_from )
  DO j = 1,n_only
    l = veg_index(idx_only(j))
    psi_leaf_bal(l) = psi1(1,j)
  END DO
END IF

! Final call only: keep the stem and leaf no drier than their critical
! points (where the plant cannot carry E, the solve returns its floor),
! and the store's flow no more than E (if it would be, the store gives all
! of E), so the store never gives more than the canopy transpires and its
! water balances. Record whether the leaf passed its critical point, for
! the leak cut.
IF ( l_ex ) THEN
  DO j = 1,n_bal
    l                   = veg_index(idx_bal(j))
    l_leaf_floor_hit(l) = psi_leaf_bal(l) < psi_floor_leaf
    psi_stem_bal(l)     = MAX(psi_stem_bal(l), psi_floor_stem)
    q_store_bal(l)      = c_p(l) * (psi_x0(l) - psi_stem_bal(l))
    IF ( q_store_bal(l) > e_c(l) ) THEN
      psi_stem_bal(l) = psi_x0(l) - e_c(l) / c_p(l)
      q_store_bal(l)  = e_c(l)
    END IF

    ! The leaf floor, but never wetter than the stem it draws from (with
    ! a leaf segment more vulnerable than the stem, psi_floor_leaf can lie
    ! above the floored stem).
    psi_leaf_bal(l) = MAX(psi_leaf_bal(l), MIN(psi_floor_leaf, psi_stem_bal(l)))
  END DO
END IF

! Final call only: the plant's conductance for the output kplant, worked
! out at the final E with no store flow (q = 0), as it would be with the
! store off, so kplant means the same with and without capacitance.
IF ( l_ex ) THEN
  q_try(:) = 0.0
  CALL store_cpl_eval( n_bal, idx_bal, l_imp_bal, .FALSE., e_c, q_try, c_p,    &
                       psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1, f_try, &
                       k1_out = k1 )
  DO j = 1,n_bal
    kplant_bal(veg_index(idx_bal(j))) = k1(1,j)
  END DO
END IF

END SUBROUTINE store_cpl_q

!-----------------------------------------------------------------------------
! Plant capacitance, two-leaf canopy: the leaf water potential of the
! sunlit and of the shaded leaves at point l (psi_leaf_cls) once the
! store's balance is committed (store_cpl_q), for the leaf segment's
! damage (update_xylem_impairment).
!
! Each class's leaf segment draws from the committed stem psi psi_stem_bal
! and carries its share of the canopy flow the commit used, e_store_bal (the
! stomatal flow, or the leak's after its cut),
!
!   E_c = s_c e_store_bal,
!   s_sun = E_sun / (E_sun + E_shd),   s_shd = 1 - s_sun,
!
! with E_sun and E_shd the classes' stomatal transpiration (el_sun,
! el_shd), or s_sun = f_sun, the sunlit share of the leaf area, where they
! do not transpire (the night leak). The class's leaf segment has the
! maximum conductance the optimiser gives it, kmax_c seg_kmax_fac, with
! kmax_c the class's share of the canopy conductance by N-weighted leaf
! area (kmax_ref_sun_2l, kmax_ref_shd_2l), and under xylem impairment the
! canopy leaf segment's cap (ximpair_seg_frac), as in store_cpl_eval. It
! is solved by one step of the lookup-table chain (lut_seg_step), with
! either som_psi_solver (the tables are built here when needed), and
! floored as the canopy leaf psi is in store_cpl_q. Where E_c / kmax_c is
! the same for both classes, both get the canopy leaf psi psi_leaf_bal.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_cls_psi( l )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l
                            ! Land point.

INTEGER :: ic
LOGICAL :: l_cap_c
                            ! The leaf segment is capped.
REAL(KIND=real_jlslsm) :: s_c(2)
                            ! Each class's share of e_store_bal (-).
REAL(KIND=real_jlslsm) :: kmx_c(2)
                            ! Each class's leaf segment maximum conductance.
REAL(KIND=real_jlslsm) :: e_t
                            ! The classes' stomatal transpiration (mol m-2
                            ! s-1).
REAL(KIND=real_jlslsm) :: fr_c
                            ! The leaf segment's cap fraction (-).

CALL build_supply_lut_seg( ft )

e_t = MAX(el_sun(l), 0.0) + MAX(el_shd(l), 0.0)
IF ( e_t > 0.0 ) THEN
  s_c(1) = MAX(el_sun(l), 0.0) / e_t
ELSE
  s_c(1) = w_cls(l,1)
END IF
s_c(2) = 1.0 - s_c(1)

! The intact conductance with the cap, or without impairment the
! conductance the optimiser used (as store_cpl_eval).
fr_c    = 1.0
l_cap_c = .FALSE.
IF ( .NOT. is_spinup .AND.                                                     &
     pft_xylem_impairment_model(ft) /= xylem_impairment_none ) THEN
  kmx_c(1) = kmax_ref_sun_2l(l)
  kmx_c(2) = kmax_ref_shd_2l(l)
  fr_c     = MIN(MAX(ximpair_seg_frac(ft, l, 3, kmax_store(l)                  &
                     / MAX(kmax_ref_store(l),                                  &
                           TINY(1.0_real_jlslsm))), 1.0e-12),                  &
                 1.0)
  l_cap_c  = fr_c < 1.0 - 1.0e-6
ELSE
  kmx_c(1) = kmax_sun_2l(l)
  kmx_c(2) = kmax_shd_2l(l)
END IF

DO ic = 1,2
  psi_leaf_cls(l,ic) = lut_seg_step( ft, 3, psi_stem_bal(l),                   &
                                     s_c(ic) * e_store_bal(l),                 &
                                     kmx_c(ic) * seg_kmax_fac(ft, l, 3), fr_c, &
                                     l_cap_c, .FALSE. )
  psi_leaf_cls(l,ic) = MAX(psi_leaf_cls(l,ic),                                 &
                           MIN(psi_floor_leaf, psi_stem_bal(l)))
END DO

END SUBROUTINE store_cpl_cls_psi

!-----------------------------------------------------------------------------
! Two-leaf canopy with the shared root and stem (l_som_shared_stem) and no
! store: the plant's conductance for the output kplant (leaf_k, per m2
! ground here). Each class's path conductance includes the whole root and
! stem, so the sum of the two would count them twice. Instead the whole
! plant (both classes' leaves, the root and stem, and the soil-to-root
! link) is solved once at the canopy's transpiration E_sun + E_shd, as
! store_cpl_q gives kplant with the store (store_cpl_eval with q = 0):
!
!   leaf_k = -dE / dpsi_leaf   of the whole plant.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
SUBROUTINE kplant_shared_stem()

IMPLICIT NONE

INTEGER :: j, l, m, n_k
INTEGER :: idx_k(land_pts)
                            ! The points solved, as positions in veg_index.
LOGICAL :: l_imp_k
                            ! Xylem impairment active.
REAL(KIND=real_jlslsm) :: e_k(land_pts), q0(land_pts), c0(land_pts),         &
                          psi0(land_pts), f0(land_pts)
                            ! Canopy E (mol m-2 s-1); no store flow, no
                            ! capacitance and the store's psi (unused);
                            ! the store's mismatch (unused).
REAL(KIND=real_jlslsm) :: kref_k(land_pts), kmx_k(land_pts),                  &
                          kcr_k(land_pts), kcap_k(land_pts)
                            ! The canopy's intact and impaired kmax, kcrit
                            ! and their ratio.
REAL(KIND=real_jlslsm) :: psi1(1,land_pts), st1(1,land_pts), k1(1,land_pts)
                            ! The leaf and stem psi, and the conductance.

l_imp_k = .NOT. is_spinup .AND.                                                &
          pft_xylem_impairment_model(ft) /= xylem_impairment_none
e_k(:)    = 0.0
q0(:)     = 0.0
c0(:)     = 0.0
psi0(:)   = 0.0
kref_k(:) = 0.0
kmx_k(:)  = 0.0
kcr_k(:)  = 0.0
kcap_k(:) = 1.0
n_k = 0
DO m = 1,veg_pts
  l = veg_index(m)
  share_sup(l) = 1.0
  kref_k(l) = kmax_ref_sun_2l(l) + kmax_ref_shd_2l(l)
  IF ( kref_k(l) <= 0.0 ) CYCLE
  kmx_k(l)  = kmax_sun_2l(l) + kmax_shd_2l(l)
  kcr_k(l)  = kcrit_sun_2l(l) + kcrit_shd_2l(l)
  kcap_k(l) = kmx_k(l) / kref_k(l)
  e_k(l)    = MAX(el_sun(l) + el_shd(l), 0.0)
  n_k = n_k + 1
  idx_k(n_k) = m
END DO
IF ( n_k == 0 ) RETURN

CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
CALL store_cpl_eval( n_k, idx_k, l_imp_k, .FALSE., e_k, q0, c0, psi0,          &
                     kref_k, kmx_k, kcr_k, kcap_k, psi1, st1, f0,              &
                     k1_out = k1 )
DO j = 1,n_k
  leaf_k(veg_index(idx_k(j))) = k1(1,j)
END DO

END SUBROUTINE kplant_shared_stem

!-----------------------------------------------------------------------------
! Plant capacitance: one guess of the search in store_cpl_q. For a guessed
! store flow q_try, solve the plant's water potentials and say how far the
! guess is from what the store would really give.
!
! With the store giving q_try, the roots and stem carry E - q_try from the
! soil, which sets the stem's water potential psi_st (st1) and the leaf's
! (psi1). At that stem water potential the store would give
! C (psi_x0 - psi_st) / dt, so the mismatch is
!
!   f_try = q_try - C (psi_x0 - psi_st) / dt      (mol H2O m-2 s-1)
!
! It is 0 at the answer, below 0 when the guess is too small and above 0
! when it is too large. Only the n_bal points being solved are worked out.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE store_cpl_eval( n_bal, idx_bal, l_imp_bal, l_ex, e_c, q_try, c_p,   &
                           psi_x0, kref_c, kmx_c, kcr_c, kcap_c, psi1, st1,    &
                           f_try, k1_out, psi_from )

IMPLICIT NONE

INTEGER, INTENT(IN) :: n_bal
                            ! Number of points solved.
INTEGER, INTENT(IN) :: idx_bal(land_pts)
                            ! Their positions in veg_index.
LOGICAL, INTENT(IN) :: l_imp_bal
                            ! Xylem impairment active.
LOGICAL, INTENT(IN) :: l_ex
                            ! Exact segment solves.
REAL(KIND=real_jlslsm), INTENT(IN) :: e_c(land_pts)
                            ! Canopy E (mol m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: q_try(land_pts)
                            ! The guessed store flow q (mol m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: c_p(land_pts)
                            ! C / dt (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_x0(land_pts)
                            ! The store's psi at the start of the step (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: kref_c(land_pts)
                            ! The canopy's intact kmax (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: kmx_c(land_pts)
                            ! The canopy's impaired kmax (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: kcr_c(land_pts)
                            ! The canopy's kcrit (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: kcap_c(land_pts)
                            ! kmx_c / kref_c, the impaired fraction (-).
REAL(KIND=real_jlslsm), INTENT(OUT) :: psi1(1,land_pts)
                            ! Leaf psi at each point solved (Pa).
REAL(KIND=real_jlslsm), INTENT(OUT) :: st1(1,land_pts)
                            ! Stem psi at each point solved (Pa).
REAL(KIND=real_jlslsm), INTENT(OUT) :: f_try(land_pts)
                            ! R at q_try (mol m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: k1_out(1,land_pts)
                            ! The chain's conductance (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: psi_from(land_pts)
                            ! Where the chain starts (Pa; psi_src if
                            ! absent).

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: e1(1,land_pts)
                            ! Canopy E at each point solved (mol m-2 s-1).
REAL(KIND=real_jlslsm) :: k1(1,land_pts)
                            ! The chain's conductance (mol m-2 s-1 Pa-1).
REAL(KIND=real_jlslsm) :: psi_in(land_pts)
                            ! Where the chain starts (Pa).

DO j = 1,n_bal
  e1(1,j) = e_c(veg_index(idx_bal(j)))
END DO
psi_in(:) = psi_src(:)
IF ( PRESENT(psi_from) ) psi_in(:) = psi_from(:)
CALL set_store_q_path( land_pts, veg_pts, veg_index, share_sup, q_try )
! som_psi_solver = 3: the segments from their supply tables, also for the
! commit (l_ex), as the table inversion is exact to its interpolation and a
! flow beyond what a segment carries ends at the bottom of its table, below
! the floors applied at the commit.
IF ( som_psi_solver == psi_solver_lut .AND. l_imp_bal ) THEN
  CALL psi_leaf_segments_lut_jls( ft, 1, land_pts, n_bal, veg_index, idx_bal,  &
                                  e1(:,1:n_bal), psi_in, kref_c,               &
                                  psi1(:,1:n_bal), k1(:,1:n_bal),              &
                                  kcap_frac = kcap_c,                          &
                                  stem_psi = st1(:,1:n_bal) )
ELSE IF ( som_psi_solver == psi_solver_lut ) THEN
  CALL psi_leaf_segments_lut_jls( ft, 1, land_pts, n_bal, veg_index, idx_bal,  &
                                  e1(:,1:n_bal), psi_in, kmx_c,                &
                                  psi1(:,1:n_bal), k1(:,1:n_bal),              &
                                  stem_psi = st1(:,1:n_bal) )
ELSE IF ( l_imp_bal ) THEN
  CALL psi_leaf_segments_jls( ft, 1, land_pts, n_bal, veg_index, idx_bal,      &
                              e1(:,1:n_bal), psi_in, kref_c, kcr_c,            &
                              psi1(:,1:n_bal), k1(:,1:n_bal),                  &
                              kcap_frac = kcap_c,                              &
                              stem_psi = st1(:,1:n_bal), l_exact = l_ex )
ELSE
  CALL psi_leaf_segments_jls( ft, 1, land_pts, n_bal, veg_index, idx_bal,      &
                              e1(:,1:n_bal), psi_in, kmx_c, kcr_c,             &
                              psi1(:,1:n_bal), k1(:,1:n_bal),                  &
                              stem_psi = st1(:,1:n_bal), l_exact = l_ex )
END IF
CALL clear_store_q_path()
DO j = 1,n_bal
  l = veg_index(idx_bal(j))
  f_try(l) = q_try(l) - c_p(l) * (psi_x0(l) - st1(1,j))
END DO
IF ( PRESENT(k1_out) ) k1_out(1,1:n_bal) = k1(1,1:n_bal)

END SUBROUTINE store_cpl_eval


!-----------------------------------------------------------------------------
! Multilayer canopy: run the stomatal optimisation for the sunlit or the
! shaded leaves of the current layer.
!
! First the soil's water supply is turned into a cap on the leaves'
! stomatal conductance, at their humidity deficit dq_x and temperature
! t_x, and the leaves are given their share of the soil-to-root
! conductance. With the stem store, the supply is the leaves' share of the
! soil's supply plus their own part of the store's flow out (q_c), and
! that flow is held fixed during the leaves' search. Then the optimisation
! runs with the leaves' photosynthesis inputs (arguments ending _x). With
! l_leaf_temp these are at the leaves' own temperature; otherwise they are
! the canopy values at tstar, and the sunlit and shaded leaves get the
! same cap, as both see the canopy deficit dqc and tstar.
!
! The leaves' transpiration and CO2 uptake go through their boundary
! layer r_bl_x (s m-1, leaf area; 0 without l_leaf_temp), in series with
! the stomata, in the optimisation and in the supply cap. The leaves work
! per unit leaf area, so the CMax and CGain costs take a leaf area of 1
! (lai_cost), not the layer's N-profile share of kmax.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE ml_class( vcmax_x, kc_x, ko_x, ccp_x, km_x, dq_x, je_x,             &
                     t_x, je_ratio_x, fapar_x, r_bl_x, rd_c, ci_c, anetl_c,    &
                     el_c, flux_o3_c, fo3_c, gl_c, psi_leaf_c, cg_c, hc_c,     &
                     leaf_k_c, psi_stem_c, q_c )
IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  vcmax_x(land_pts), kc_x(land_pts), ko_x(land_pts), ccp_x(land_pts),          &
  km_x(land_pts), dq_x(land_pts), je_x(land_pts),                              &
  t_x(land_pts), je_ratio_x(land_pts), fapar_x(land_pts), r_bl_x(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: rd_c(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ci_c(land_pts), anetl_c(land_pts), el_c(land_pts), flux_o3_c(land_pts),      &
  fo3_c(land_pts), gl_c(land_pts), psi_leaf_c(land_pts), cg_c(land_pts),       &
  hc_c(land_pts), leaf_k_c(land_pts)
                            ! OUT from stom_opt_mod, which sets them over
                            ! the whole array.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: psi_stem_c(land_pts)
                            ! The stem psi of this leaf's path (Pa; with
                            ! plant segments).
REAL(KIND=real_jlslsm), INTENT(IN) :: q_c(land_pts)
                            ! The capacitance store: this leaf's part of the
                            ! store's flow out (mol H2O m-2 leaf s-1).

INTEGER :: i, l
REAL(KIND=real_jlslsm) :: e_sup_cls(land_pts)
                            ! With the store: this leaf's supply, its kmax
                            ! share of the soil supply plus its own part of
                            ! the store's flow (kg m-2 s-1).

IF ( l_cap_store ) THEN
  ! This leaf's kmax share of the soil supply plus its own part of the
  ! store's flow (q_c), not share_sup (e_s + q_store).
  e_sup_cls(:) = e_supply(:)
  WHERE ( e_supply(:) >= 0.0 )                                                 &
    e_sup_cls(:) = MAX(share_sup(:) * e_supply(:)                              &
                       + q_c(:) * Water_kg_per_mol, 0.0)
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_cls,            &
                           share_one, dq_x, t_x, pstar, gl_max_lf, gl_max_eff, &
                           r_bl_x, zeros_pts, zeros_pts )
ELSE
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_use,            &
                           share_sup, dq_x, t_x, pstar, gl_max_lf, gl_max_eff, &
                           r_bl_x, zeros_pts, zeros_pts )
END IF
CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )

! The store's flow out held fixed for this leaf's search.
IF ( l_cap_store ) CALL set_store_q_path( land_pts, veg_pts, veg_index,        &
                                          share_one, q_c )
CALL stom_opt_mod (                                                            &
    ! IN
    land_pts, ft, open_pts, open_index,                                        &
    pft_photo_model, veg_pts, veg_index,                                       &
    ca, psi_src, acr, apar, oa, vcmax_x, kc_x, ko_x, ccp_x, pstar,             &
    km_x, dq_x, je_x, t_x, je_ratio_x, fapar_x,                                &
    kmax_imp_lyr, kcrit_leaf_lyr, gl_max_eff, ipar,                            &
    l_multilayer, r_bl_x, zeros_pts, zeros_pts,                                &
    kmax_leaf_lyr, b_pft_pts, c_pft_pts,                                       &
    psi_leaf_extreme, psi_root_extreme, l_ml_ximp,                             &
    ! IN OUT
    rd_c,                                                                      &
    ! OUT
    ci_c, anetl_c, el_c, flux_o3_c, fo3_c, gl_c, psi_leaf_c, cg_c, hc_c,       &
    leaf_k_c, psi_stem = psi_stem_c, lai_cost = ones_pts )
IF ( l_cap_store ) CALL clear_store_q_path()

END SUBROUTINE ml_class

!-----------------------------------------------------------------------------
! Two-leaf canopy (can_rad_mod = 7): run the stomatal optimisation for one
! set of leaves, sunlit or shaded, on one pass. Arguments ending _c are this
! set's; those ending _o are the other set's.
!
! Steps:
!
!   1. Split the points into open (the leaves get light) and closed.
!   2. Give these leaves their share of the soil's water supply and of the
!      soil-to-root conductance, in proportion to their share of the
!      canopy's intact hydraulic conductance (share_sup). With the stem
!      store, add the store's flow to these leaves (q_store_cls).
!   3. Set up the leaves' transpiration, which also depends on how open the
!      other leaves are (gl_o), as both lose water to the same canopy air.
!   4. Run the optimisation, at these leaves' own temperature with
!      l_leaf_temp, or at the canopy temperature tstar otherwise.
!   5. At points that end up closed, share the closed canopy's conductance
!      (min_gl_pft) between the two sets by leaf area.
!
! With a shared root and stem (l_som_shared_stem), the two sets of leaves
! sit on one plant: these leaves' root and stem are the whole plant's, and
! they also carry the other leaves' transpiration (e_o_c). So the whole
! store's flow, the whole soil-to-root conductance and the soil's supply
! less the other leaves' use are seen here, while these leaves keep their
! own share of the leaf conductance. The critical conductance (kcrit) stays
! the same fraction of the plant's intact conductance, now of the whole
! plant's root and stem.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE two_leaf_class( apar_c, lai_c, kmax_ref_c, qs_lt_c, dq_lt_c,        &
                           gbw_c, gbw_o, gl_o, t_lt_c, gl_max_c, acr_c,        &
                           vcmax_c, kc_lt_c, ko_lt_c, ccp_lt_c, km_lt_c,       &
                           je_c, kmax_c, kcrit_c, rd_c, ci_c, anetl_c, el_c,   &
                           flux_o3_c, fo3_c, gl_c, psi_leaf_c, cg_c, hc_c,     &
                           leaf_k_c, psi_stem_c, lam_c, q_store_cls, e_o_c,    &
                           a_o_c, amax_c )
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
                            ! OUT from stom_opt_mod, which sets them over
                            ! the whole array.
REAL(KIND=real_jlslsm), INTENT(OUT) :: lam_c(land_pts)
                            ! dA/dE at the optimum (lambda_eff of
                            ! stom_opt_mod; 0 where not set).
REAL(KIND=real_jlslsm), INTENT(IN) :: q_store_cls(land_pts)
                            ! The capacitance store: this class's part of
                            ! the store's flow out (q_store_sun or q_store_shd).
REAL(KIND=real_jlslsm), INTENT(IN) :: e_o_c(land_pts)
                            ! The shared root and stem: the other class's
                            ! transpiration (mol H2O m-2 ground s-1, >= 0),
                            ! held fixed. Not used otherwise.
REAL(KIND=real_jlslsm), INTENT(IN) :: a_o_c(land_pts)
                            ! The shared root and stem: the other class's A
                            ! at its hydraulic edge (mol CO2 m-2 ground s-1),
                            ! for the whole-plant gain (l_som_shared_gain).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: amax_c(land_pts)
                            ! The shared root and stem: this class's A at its
                            ! hydraulic edge from this search (0 where it was
                            ! not searched); unchanged otherwise.

INTEGER :: i, l, m
REAL(KIND=real_jlslsm) :: kcrit_use(land_pts), k_own(land_pts),               &
                          kcrit_pl(land_pts)
                            ! The shared root and stem: the class's kcrit on
                            ! the path with the whole root and stem, and the
                            ! intact conductance of the class's own path
                            ! (relative to its kmax) it is scaled from.
REAL(KIND=real_jlslsm) :: e_sup_cls(land_pts)
                            ! With the store: the class's supply, its kmax
                            ! share of the soil supply plus its own part of
                            ! the store's flow (kg m-2 s-1).

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

CALL coupled_inputs( land_pts, veg_pts, veg_index, l_leaf_temp, dq_min, dqc,   &
                     dq, qs_lt_c, q1, dq_lt_c, gbw_c, gbw_o, ra_ca_lt, gl_o,   &
                     g_soil_c, qs, dq_opt, r_bl_opt, r_ca_opt, g_other_opt )

! The store coupled to the optimiser (the capacitance store): this class's
! part of the store's flow out (store_cpl_split), held fixed for its search
! (cleared after). With the shared root and stem the whole store's flow,
! which leaves at the shared stem node.
IF ( l_cap_store .AND. l_som_shared_stem ) THEN
  CALL set_store_q_path( land_pts, veg_pts, veg_index, share_one, q_store )
ELSE IF ( l_cap_store ) THEN
  CALL set_store_q_path( land_pts, veg_pts, veg_index, share_one, q_store_cls )
END IF

! The shared root and stem (l_som_shared_stem): kcrit as kcrit_of_path, a
! fixed fraction of the path's intact conductance, which with the whole
! root and stem is k_sh / k_own times the class's own,
!
!   kcrit_sh = kcrit_c k_sh / k_own,
!
! k_sh and k_own the path's intact conductance at psi = 0 with and without
! the whole root and stem (k_path_zero_flow, plant only, per unit kmax).
! Then the class's root and stem are set to the whole plant's, with the
! other class's flow (cleared after the search). The cost is the whole
! plant's (stom_opt_mod), with the whole plant's kcrit, kcrit_of_path of
! its intact conductance: kcrit_c / share_sup in the class's units.
kcrit_use(:) = kcrit_c(:)
IF ( l_som_shared_stem ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    k_own(l) = k_path_zero_flow(ft, l, 1.0_real_jlslsm, 0.0_real_jlslsm,      &
                                l_plant_only = .TRUE.)
    kcrit_pl(l) = kcrit_c(l) / MAX(share_sup(l), TINY(1.0_real_jlslsm))
  END DO
  CALL set_shared_stem( land_pts, veg_pts, veg_index, share_sup, e_o_c,        &
                        kcrit_pl, a_o_c, l_som_shared_gain )
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( k_own(l) > TINY(1.0_real_jlslsm) )                                    &
      kcrit_use(l) = kcrit_c(l)                                                &
                     * k_path_zero_flow(ft, l, 1.0_real_jlslsm,                &
                                        0.0_real_jlslsm,                       &
                                        l_plant_only = .TRUE.) / k_own(l)
  END DO
END IF

! At the class's leaf temperature (l_leaf_temp), else tstar: the class
! arrays hold the tile's values then (see the two-leaf set-up).
IF ( l_som_shared_stem ) THEN
  ! The shared root and stem: the canopy is capped by the soil supply plus
  ! the store's flow (q_store, with the store coupled),
  !
  !   E + E_o <= e_supply + q_store,
  !
  ! so the class by what is left after the other class's E.
  e_sup_cls(:) = e_supply(:)
  DO m = 1,veg_pts
    l = veg_index(m)
    IF ( e_supply(l) >= 0.0 ) THEN
      e_sup_cls(l) = e_supply(l) - e_o_c(l) * Water_kg_per_mol
      IF ( l_cap_store )                                                       &
        e_sup_cls(l) = e_sup_cls(l) + q_store(l) * Water_kg_per_mol
      e_sup_cls(l) = MAX(e_sup_cls(l), 0.0)
    END IF
  END DO
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_cls,            &
                           share_one, dq_opt, t_lt_c, pstar, gl_max_c,         &
                           gl_max_eff, r_bl_opt, r_ca_opt, g_other_opt )
ELSE IF ( l_cap_store ) THEN
  ! The class's kmax share of the soil supply plus its own part of the
  ! store's flow (q_store_cls), not share_sup (e_s + q_store).
  e_sup_cls(:) = e_supply(:)
  WHERE ( e_supply(:) >= 0.0 )                                                 &
    e_sup_cls(:) = MAX(share_sup(:) * e_supply(:)                              &
                       + q_store_cls(:) * Water_kg_per_mol, 0.0)
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_cls,            &
                           share_one, dq_opt, t_lt_c, pstar, gl_max_c,         &
                           gl_max_eff, r_bl_opt, r_ca_opt, g_other_opt )
ELSE
  CALL apply_supply_limit( land_pts, veg_pts, veg_index, e_sup_use,            &
                           share_sup, dq_opt, t_lt_c, pstar, gl_max_c,         &
                           gl_max_eff, r_bl_opt, r_ca_opt, g_other_opt )
END IF

! The shared root and stem: the whole soil-to-root conductance.
IF ( l_som_shared_stem ) THEN
  CALL set_ksr_path( land_pts, veg_pts, veg_index, share_one )
ELSE
  CALL set_ksr_path( land_pts, veg_pts, veg_index, share_sup )
END IF
CALL stom_opt_mod (                                                            &
        ! IN
          land_pts, ft, open_pts, open_index,                                  &
          pft_photo_model, veg_pts, veg_index,                                 &
          ca, psi_src, acr_c, apar_c, oa, vcmax_c,                             &
          kc_lt_c, ko_lt_c, ccp_lt_c, pstar,                                   &
          km_lt_c, dq_opt, je_c, t_lt_c, je_dummy,                             &
          fapar_dummy,                                                         &
          kmax_c, kcrit_use, gl_max_eff, ipar,                                 &
          l_multilayer, r_bl_opt, r_ca_opt, g_other_opt,                       &
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
IF ( l_cap_store ) CALL clear_store_q_path()

! The shared root and stem: this class's A at its edge, then off.
IF ( l_som_shared_stem ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    amax_c(l) = shared_amax_self(l)
  END DO
  CALL clear_shared_stem()
END IF

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
! Prepare the stomatal search to correct photosynthesis for leaf
! temperature (l_leaf_temp_ci), for the sunlit (ic = 1) or shaded
! (ic = 2) leaves.
!
! Opening the stomata more cools the leaf, because it transpires more, and
! a cooler leaf photosynthesises differently and has a smaller humidity
! deficit. So each stomatal opening the search tries should be scored at
! its own leaf temperature. This routine gives the search what it needs
! to do that cheaply, without solving the leaf energy balance for every
! candidate.
!
! The leaf temperature is a straight line in the transpiration E, from
! the last pass's leaf energy balance (tleaf_a_2l, tleaf_b_2l):
!
!   T_leaf = a - b E      (K)
!
! The routine works out the leaves' photosynthesis parameters and
! saturation humidity at their temperature t_cls (from
! leaf_class_at_temp) and at t_cls + 1 K. From these it takes how each
! parameter changes per kelvin, as d ln(p) / dT (0 where a parameter is
! not positive at both temperatures), and how the saturation humidity
! changes per kelvin, dqs / dT (kg kg-1 K-1). It passes these to the
! search (set_leaf_temp_ci).
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE ltc_set( ic, t_cls, nw_c, i2c, ccp_c, kc_c, ko_c, km_c, vcmax_c,    &
                    je_c, rd_c, qs_c )
IMPLICIT NONE
INTEGER, INTENT(IN) :: ic
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  t_cls(land_pts), nw_c(land_pts), i2c(land_pts), ccp_c(land_pts),             &
  kc_c(land_pts), ko_c(land_pts), km_c(land_pts), vcmax_c(land_pts),           &
  je_c(land_pts), rd_c(land_pts), qs_c(land_pts)
REAL(KIND=real_jlslsm) ::                                                      &
  t1(land_pts), ccp1(land_pts), kc1(land_pts), ko1(land_pts), km1(land_pts),   &
  vcmax1(land_pts), je1(land_pts), rd1(land_pts), qs1(land_pts), dq1(land_pts),&
  p0(7), p1(7)
INTEGER :: l, m, ip

t1(:) = t_cls(:) + 1.0
ccp1(:) = ccp_c(:)
kc1(:)  = kc_c(:)
ko1(:)  = ko_c(:)
km1(:)  = km_c(:)
CALL leaf_class_at_temp( ft, land_pts, veg_pts, veg_index, pft_photo_model,    &
                         t1, oa, actj, actv, dsj, dsv, jv25, nleaf_top, nw_c,  &
                         i2c, pstar, q_c_lt, dq_min, ccp1, kc1, ko1, km1,      &
                         vcmax1, je1, rd1, qs1, dq1,                           &
                         f_vcmax_down = f_vcmax_down )
dlnpar_dt_2l(:,:) = 0.0
dqs_dt_2l(:)      = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  p0 = (/ vcmax_c(l), je_c(l), kc_c(l), ko_c(l), ccp_c(l), km_c(l), rd_c(l) /)
  p1 = (/ vcmax1(l), je1(l), kc1(l), ko1(l), ccp1(l), km1(l), rd1(l) /)
  DO ip = 1,7
    IF ( p0(ip) > 0.0 .AND. p1(ip) > 0.0 )                                     &
      dlnpar_dt_2l(l,ip) = LOG(p1(ip) / p0(ip))
  END DO
  dqs_dt_2l(l) = qs1(l) - qs_c(l)
END DO
CALL set_leaf_temp_ci( land_pts, veg_pts, veg_index, tleaf_a_2l(:,ic),         &
                       tleaf_b_2l(:,ic), dqs_dt_2l, dlnpar_dt_2l )

END SUBROUTINE ltc_set

!-----------------------------------------------------------------------------
! Get the leaf energy balance (l_leaf_temp) ready before the passes of the
! stomatal optimisation, in the two-leaf and multilayer canopies. Each
! pass solves the leaves' temperatures, so the terms that stay the same
! on every pass are worked out once here, with first guesses for the rest.
!
!   1. Work out the radiation and wind terms that stay the same on every
!      pass, from the leaf area and absorbed PAR of the sunlit and shaded
!      leaves, and for every layer as well in the multilayer canopy.
!   2. Make first guesses of the leaf temperatures (the last timestep's,
!      or tstar on the first timestep), and set the air the leaves are in
!      to the first atmospheric level's (tair, q1).
!   3. Start the photosynthesis temperature terms of the sunlit and shaded
!      leaves (CO2 compensation point and Michaelis-Menten constants) from
!      the tile's.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE lt_init()
IMPLICIT NONE
INTEGER :: l, m
    ! Its own loop counters, not the host's.

  ! The wind extinction with leaf area (0 by default).
  wind_ext(:) = leaf_wind_ext

  IF ( l_leaf_temp_lyr ) THEN
    ! Also the terms of every layer's sunlit and shaded leaf.
    CALL leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai,      &
                          dlai, fsun, fapar_sun_nir, fapar_shd_nir, ipar,      &
                          wind_ext, apar_sun_2l, apar_shd_2l, sw_sun_lt,       &
                          sw_shd_lt, lwwt_sky_sun_lt, lwwt_sky_shd_lt,         &
                          lwwt_soil_sun_lt, lwwt_soil_shd_lt, grad_fac_sun_lt, &
                          grad_fac_shd_lt, lai_wind_sun_lt, lai_wind_shd_lt,   &
                          fapar_sun_vis = fapar_sun,                           &
                          fapar_shd_vis = fapar_shd,                           &
                          lai_c = lai_lt_c, sw_c = sw_lt_c,                    &
                          lwwt_sky_c = lwwt_sky_lt_c,                          &
                          lwwt_soil_c = lwwt_soil_lt_c,                        &
                          grad_fac_c = grad_fac_lt_c,                          &
                          lai_wind_c = lai_wind_lt_c )
  ELSE
  CALL leaf_temp_setup( ft, land_pts, veg_pts, veg_index, ilayers, lai, dlai,  &
                        fsun, fapar_sun_nir, fapar_shd_nir, ipar, wind_ext,    &
                        apar_sun_2l, apar_shd_2l, sw_sun_lt, sw_shd_lt,        &
                        lwwt_sky_sun_lt, lwwt_sky_shd_lt, lwwt_soil_sun_lt,    &
                        lwwt_soil_shd_lt, grad_fac_sun_lt, grad_fac_shd_lt,    &
                        lai_wind_sun_lt, lai_wind_shd_lt )
  END IF

  ! Defined at every land point (leaf_class_at_temp takes qsat of the
  ! whole array, and the per-layer classes copy it).
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
    ! The air the leaves are in: the first atmospheric level's.
    t_c_lt(l) = tair(l)
    q_c_lt(l) = q1(l)
  END DO

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
! Finish the leaf energy balance (l_leaf_temp) at point lp after the
! passes of the stomatal optimisation, in the two-leaf and multilayer
! canopies: keep the leaf temperatures the passes ended with and hand the
! results to the rest of JULES.
!
!   1. Keep the leaf temperatures used (t_leaf_sun, t_leaf_shd) and their
!      mean t_can, weighted by the sunlit fraction of the leaf area f_sun.
!      Leaves with no leaf area take the other leaves' temperature.
!   2. Fill the canopy diagnostics (canopy_diag_mod).
!   3. With l_leaf_temp_gc_eq, replace the canopy conductance gc with the
!      one that gives the transpiration the optimisation chose when JULES
!      applies gc at the surface temperature tstar:
!
!        rho (qs(tstar) - q1) gc / (1 + ra (gc + g_s)) = E M_w
!
!      where rho is the air density, qs(tstar) - q1 the humidity deficit
!      from the surface to the first atmospheric level (kg kg-1), ra the
!      aerodynamic resistance (s m-1), g_s the soil's evaporation
!      conductance, E the transpiration (mol m-2 s-1) and M_w the molar
!      mass of water. gc is no more than gc_eq_cap, set by the caller, or
!      som_gl_max times LAI where the caller sets none. The optimised gc
!      is kept when the deficit is below 1e-4 kg kg-1.
!   4. Record the fraction of the leaves (sunlit and shaded, or those of
!      every layer, n_cls_in) whose temperature was held at dt_leaf from
!      the air, and warn the first few times it happens.
!
! MDK, 12/11/26
!-----------------------------------------------------------------------------
SUBROUTINE lt_finish( lp, f_sun, n_cls_in )
IMPLICIT NONE
INTEGER, INTENT(IN) :: lp
REAL(KIND=real_jlslsm), INTENT(IN) :: f_sun
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: n_cls_in
    ! Number of leaf classes with leaf area, for t_leaf_clamp (2 at most
    ! when absent; in the multilayer canopy those of all the layers).
REAL(KIND=real_jlslsm) :: n_cls
REAL(KIND=real_jlslsm) :: dq_eq, cap_eq
    ! gc_eq: the tile's actual deficit qs(tstar) - q1 (kg kg-1) and the
    ! cap on gc (m s-1).
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

  ! Diagnostics (canopy_diag_mod): the stomatal gc before gc_eq replaces
  ! it, the air the leaves are in (tair, q1) and the leaf fluxes of the
  ! energy balance.
  gc_stomatal_pft(lp,ft) = gc(lp)
  t_can_air_pft(lp,ft)   = t_c_lt(lp)
  q_can_air_pft(lp,ft)   = q_c_lt(lp)
  h_leaf_pft(lp,ft)      = h_leaf_pts(lp)
  le_leaf_pft(lp,ft)     = le_leaf_pts(lp)

  ! The surface energy balance applies gc with qs(tstar) - q1 through
  ! ra. Return the gc that gives the transpiration the optimisation
  ! chose (el, at the leaf temperatures) there,
  !
  !   rho (qs(tstar) - q1) / (ra + 1/gc) = el m_h2o,
  !
  ! capped at the canopy som_gl_max. The soil conductance g_s shares ra
  ! in the tile flux (soil_evap adds it to gc) and the transpiration is
  ! the stomatal share,
  !
  !   rho (qs(tstar) - q1) gc / (1 + ra (gc + g_s)) = el m_h2o,
  !
  ! so 1/gc is the r_eq_gc above divided by 1 + ra g_s.
  !
  ! The deficit is the actual qs(tstar) - q1, not dq, which the stomatal
  ! optimisation floors at dq_min: with the floor a near-saturated or
  ! dewy step would give a huge gc that JULES then applies to the real
  ! deficit. Below dq_min (or with condensation) the optimised gc is
  ! kept. The cap is always finite: som_gl_max x LAI where the canopy
  ! cap is off (supply-loss).
  IF ( l_leaf_temp_gc_eq ) THEN
    dq_eq  = qs(lp) - q1(lp)
    cap_eq = gc_eq_cap
    IF ( cap_eq <= 0.0 ) cap_eq = som_gl_max * lai(lp)
    IF ( el(lp) * m_h2o_lt > TINY(1.0) .AND. dq_eq > 1.0e-4 ) THEN
      r_eq_gc = pstar(lp) / ( r * tstar(lp) ) * dq_eq / ( el(lp) * m_h2o_lt )  &
             - ra(lp)
      r_eq_gc = r_eq_gc / ( 1.0 + ra(lp) * g_soil_c(lp) )
      IF ( r_eq_gc > 0.0 ) THEN
        gc(lp) = MIN(1.0 / r_eq_gc, cap_eq)
      ELSE
        ! Demand above what any gc gives: the cap.
        gc(lp) = cap_eq
      END IF
    END IF
  END IF

  ! Fraction of the classes with leaf area held at the dt_leaf limit.
  n_cls = 0.0
  IF ( lai_sun_2l(lp) > EPSILON(0.0) ) n_cls = n_cls + 1.0
  IF ( lai_shd_2l(lp) > EPSILON(0.0) ) n_cls = n_cls + 1.0
  IF ( PRESENT(n_cls_in) ) n_cls = n_cls_in
  IF ( n_cls > 0.0 ) t_leaf_clamp(lp) = MIN(n_clamp_lt(lp), n_cls) / n_cls
  IF ( n_clamp_lt(lp) > 0.0 .AND. n_lt_clamp_warn < max_lt_clamp_warn ) THEN
!$OMP CRITICAL (sf_stom_lt_clamp_warn)
    n_lt_clamp_warn = n_lt_clamp_warn + 1
    WRITE(clamp_msg, '(A,I0,A,I0,A,3F8.2,A,I0,A)')                             &
      'l_leaf_temp: leaf temperature held at dt_leaf from the air, '//         &
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

SUBROUTINE sun_shade_totals( land_pts, veg_pts, veg_index, ilayers, dlai,     &
                             fsun, fapar_sun, fapar_shd, ipar, lai_sun,        &
                             lai_shd, apar_sun, apar_shd )

! Add up, over the canopy layers, the leaf area of the sunlit and of the
! shaded leaves (m2 m-2) and the PAR each set absorbs (W m-2 ground). Layer
! n has leaf area dL, sunlit fraction f_n and absorbed PAR per unit sunlit
! and shaded leaf area a_sun,n and a_shd,n (fractions of the incident PAR
! I, from the two-stream profile), so
!
!   L_sun = sum_n f_n dL,          L_shd = sum_n (1 - f_n) dL,
!   A_sun = I sum_n a_sun,n f_n dL,
!   A_shd = I sum_n a_shd,n (1 - f_n) dL.
!
! The two-leaf canopy runs on these totals, and the multilayer canopy's
! leaf energy balance uses them too.
!
! MDK, 03/11/26

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts), ilayers
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  dlai(land_pts),                                                              &
      ! Leaf area of each layer (m2 m-2).
  fsun(land_pts,ilayers),                                                      &
      ! Sunlit fraction of each layer.
  fapar_sun(land_pts,ilayers), fapar_shd(land_pts,ilayers),                    &
      ! PAR absorbed per unit sunlit and shaded leaf area of each layer, as
      ! a fraction of ipar.
  ipar(land_pts)
      ! Incident PAR (W m-2).
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  lai_sun(land_pts), lai_shd(land_pts),                                        &
      ! Sunlit and shaded leaf area (m2 m-2).
  apar_sun(land_pts), apar_shd(land_pts)
      ! PAR absorbed by the sunlit and shaded leaves (W m-2 ground).

INTEGER :: m, l, n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SUN_SHADE_TOTALS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO m = 1,veg_pts
  l = veg_index(m)
  lai_sun(l)  = 0.0
  lai_shd(l)  = 0.0
  apar_sun(l) = 0.0
  apar_shd(l) = 0.0
  DO n = 1,ilayers
    lai_sun(l)  = lai_sun(l)  + fsun(l,n) * dlai(l)
    lai_shd(l)  = lai_shd(l)  + (1.0 - fsun(l,n)) * dlai(l)
    apar_sun(l) = apar_sun(l)                                                  &
                  + fapar_sun(l,n) * fsun(l,n) * dlai(l) * ipar(l)
    apar_shd(l) = apar_shd(l)                                                  &
                  + fapar_shd(l,n) * (1.0 - fsun(l,n)) * dlai(l) * ipar(l)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sun_shade_totals

!#############################################################################
!#############################################################################

SUBROUTINE leaf_class_at_temp( ft, land_pts, veg_pts, veg_index,               &
                               pft_photo_model, t_leaf, oa, actj, actv, dsj,   &
                               dsv, jv25, nleaf_top, nw, i2c, pstar, q_c,      &
                               dq_min, ccp, kc, ko, km, vcmax_c, je_c, rd_c,   &
                               qs_c, dq_c, f_vcmax_down )

! Work out the photosynthesis inputs of one set of leaves (sunlit or
! shaded) at their own temperature t_leaf, for the leaf energy balance
! (l_leaf_temp) in the two-leaf and multilayer canopies. Leaves warmer or
! cooler than the canopy photosynthesise differently, so the stomatal
! optimisation needs these inputs at the leaves' own temperature.
!
! It applies the temperature responses, then gives Vcmax and dark
! respiration (Rd) summed over the leaves, by multiplying the top leaf's
! values by the leaves' nitrogen-weighted leaf area nw, and the electron
! flux. It also gives the leaves' humidity deficit to the canopy air,
! qs(t_leaf) - q_c (kg kg-1), no lower than dq_min.
!
! MDK, 12/11/26

USE jules_vegetation_mod, ONLY: photo_farquhar
USE qsat_mod, ONLY: qsat

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
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

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: f_vcmax_down(land_pts)
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
                          vtemp, jtemp, i2_t )

CALL calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,            &
                            veg_index, denom_t, jtemp, jv25,                   &
                            nleaf_top, qtenf_t, vtemp,                         &
                            jmax_t, rd_t, vcmax_t, f_vcmax_down = f_vcmax_down )

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

SUBROUTINE leaf_temp_responses( ft, land_pts, veg_pts, veg_index,              &
                                pft_photo_model, t_leaf, oa, acr, actj, actv,  &
                                dsj, dsv, ccp, kc, ko, km, denom, qtenf_term,  &
                                vcmax_temp, jmax_temp, i2 )

! Adjust the photosynthesis parameters for the leaf temperature t_leaf.
! It is called with tstar for the canopy, and with each leaf's
! temperature when l_leaf_temp. It sets denom, qtenf_term, ccp, kc, ko
! (Collatz; C3 only for ccp, kc, ko) or ccp, km, i2, vcmax_temp,
! jmax_temp (Farquhar).
!
! MDK, 25/10/26

USE conversions_mod, ONLY: zerodegc
USE c_rmol, ONLY: rmol
USE jules_vegetation_mod, ONLY: photo_collatz, photo_farquhar
USE pftparm, ONLY: alpha_elec, c3, deact_jmax, deact_vcmax, q10_leaf, tlow,    &
                   tupp

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

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

  ! (The soil-water factor of Vcmax and Jmax, l_som_vcmax_psi, is applied
  ! in calc_photo_parameters, after dark respiration.)

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_temp_responses

!#############################################################################
!#############################################################################

SUBROUTINE vcmax_psi_factor( ft, land_pts, veg_pts, veg_index, psi_rz,         &
                             f_vcmax_down )

! Lower photosynthetic capacity (Vcmax and Jmax) as the soil dries
! (l_som_vcmax_psi), by returning a factor f (0 to 1) that multiplies
! them. It uses the form of Zhou et al. (2013) on the root-zone (predawn)
! water potential psi, as in De Kauwe et al. (2015, Biogeosciences 12:
! 7503-7518, Eq. 4):
!
!   f = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi)))
!
! with psi and psi_f (psi_vcmax_f) in MPa and sf (sf_vcmax) in MPa-1. f is
! 1 in wet soil and about 0.5 at psi = psi_f.
!
! Two changes to De Kauwe et al. (2015). Without them f falls towards 0 in
! drought, which takes GPP to about 0, and capacity comes back in the
! timestep of the first rain.
!
!   1. A floor, fmin (psi_vcmax_fmin, 0 by default):
!
!        f_vcmax_down = fmin + (1 - fmin) f
!
!      so in the driest soil photosynthesis is limited by the stomata
!      closing as the plant's water supply fails.
!
!   2. A lag in recovery. As the soil dries, the factor follows f at once
!      (psi itself changes slowly). As the soil wets, it moves towards f
!      with an e-folding time of 5 days: 63 % of the way back after 5 days
!      and 90 % after about 11 days. The 5 days is a choice between two
!      known speeds of repair: fast reversible changes (Rubisco
!      activation, mesophyll conductance) take hours to days, and making
!      new protein takes a week or more. The lagged
!      factor is kept in sf_stom_state_mod and written to the dump
!      (fvcmax_psi_state), so a restart carries on from it. A run started
!      without it starts at f_vcmax_down.
!
! MDK, 12/11/26

USE pftparm, ONLY: psi_vcmax_f, sf_vcmax, psi_vcmax_fmin
USE timestep_mod, ONLY: timestep

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_rz(land_pts)
    ! Root-zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: f_vcmax_down(land_pts)
    ! Capacity factor (-); only veg points are set.

REAL(KIND=real_jlslsm), PARAMETER :: tau_recovery = 5.0 * 86400.0
    ! Recovery e-folding time (s).

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: psi_mpa
    ! Root-zone water potential (MPa).
REAL(KIND=real_jlslsm) :: f_zhou
    ! The Zhou et al. factor f (0-1), before the floor.
REAL(KIND=real_jlslsm) :: f_target
    ! This step's factor with the floor, fmin + (1 - fmin) f, before the lag
    ! in recovery.
REAL(KIND=real_jlslsm) :: w_recovery
    ! Weight of one step of the recovery, 1 - exp(-dt / tau_recovery).


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='VCMAX_PSI_FACTOR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL sf_stom_state_alloc( land_pts )
w_recovery = 1.0 - EXP( -timestep / tau_recovery )

DO m = 1,veg_pts
  l = veg_index(m)
  psi_mpa = MIN(psi_rz(l), 0.0) * 1.0e-6
  f_zhou = ( 1.0 + EXP( sf_vcmax(ft) * psi_vcmax_f(ft) * 1.0e-6 ) ) / ( 1.0 +  &
          EXP( sf_vcmax(ft) * ( psi_vcmax_f(ft) * 1.0e-6 - psi_mpa ) ) )
  f_zhou = MIN(MAX(f_zhou, 0.0_real_jlslsm), 1.0_real_jlslsm)
  f_target = psi_vcmax_fmin(ft) + (1.0 - psi_vcmax_fmin(ft)) * f_zhou
  IF ( f_vcmax_psi_state(l,ft) < 0.0 .OR. f_target <= f_vcmax_psi_state(l,ft)  &
                                                      ) THEN
    f_vcmax_psi_state(l,ft) = f_target
  ELSE
    f_vcmax_psi_state(l,ft) = f_vcmax_psi_state(l,ft) + w_recovery * (         &
                          f_target - f_vcmax_psi_state(l,ft) )
  END IF
  f_vcmax_down(l) = f_vcmax_psi_state(l,ft)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE vcmax_psi_factor

!#############################################################################
!#############################################################################

SUBROUTINE nsc_sink_factor( ft, land_pts, veg_pts, veg_index, psi_rz, beta )

! Lower photosynthetic capacity when drought stops growth and the plant's
! carbon store (the TRIFFID non-structural carbon store) has little room
! for the sugars that growth would have used (l_nsc_sink_feedback). It
! returns a factor beta (beta_min to 1) that multiplies the capacity.
! It has the form of QUINCY (Thum et al., 2019, GMD 12: 4781, eqs
! S43-S44), set against the room left in the store.
!
! Growth falls with the turgor factor f_turg (0 to 1), which goes from 0
! at the root-zone water potential nsc_psi_turg_zero to 1 at
! nsc_psi_turg_full. The carbon the store would have released to growth
! over t = nsc_sink_days days then stays in the store:
!
!   B    = W_pot (1 - f_turg) t / 360      (kg C m-2)
!   X    = B / R
!   beta = beta_min + (1 - beta_min) exp(-X^2)
!
! where W_pot is the store's release to growth at its last step with the
! turgor limit open (kg C m-2 per 360 days): the cold and season limits
! still apply, so outside the growing season, or in the cold, W_pot is 0
! and so is the feedback. R = C_max - store is the room left in the store
! (kg C m-2). X is the share of that room the held-back carbon would fill
! in t days, so the feedback is 63 % of its full size at X = 1. In wet soil
! (f_turg = 1) there is no feedback however full the store is. In a long
! drought the room shrinks and the feedback grows. Capacity falls by at
! most 1 - beta_min.
!
! beta is 1 where the store's maximum or content is not yet set, and the
! routine leaves beta unchanged without the turgor limit
! (nsc_psi_turg_zero not set) or before the store first exists.
!
! MDK, 12/11/26

USE triffid_nsc_mod, ONLY: nsc_pft, nsc_cmax_pft, nsc_wpot_pft
USE pftparm, ONLY: is_set, nsc_psi_turg_zero, nsc_psi_turg_full,              &
                   nsc_sink_days, nsc_sink_beta_min

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_rz(land_pts)
    ! Root-zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: beta(land_pts)
    ! Sink feedback factor (-); only veg points are set.

REAL(KIND=real_jlslsm), PARAMETER :: room_min = 1.0e-6
    ! Smallest room counted (kg C m-2), so a full store gives a large X.

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: store_room
    ! Room left in the store, R = C_max - store (kg C m-2).
REAL(KIND=real_jlslsm) :: c_held_back
    ! Carbon kept from growth over nsc_sink_days, B (kg C m-2).
REAL(KIND=real_jlslsm) :: room_fill_ratio
    ! Share of the room the held-back carbon would fill, X = B / R (-).
REAL(KIND=real_jlslsm) :: f_turg
    ! Turgor factor on growth (0-1).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='NSC_SINK_FACTOR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( ALLOCATED(nsc_pft) .AND. ALLOCATED(nsc_cmax_pft) .AND.                    &
     ALLOCATED(nsc_wpot_pft) .AND. is_set(nsc_psi_turg_zero(ft)) ) THEN
  DO m = 1,veg_pts
    l = veg_index(m)
    beta(l) = 1.0
    IF ( nsc_cmax_pft(l,ft) <= 0.0 .OR. nsc_pft(l,ft) < 0.0 ) CYCLE

    f_turg = MIN(MAX((1.0e-6 * psi_rz(l) - nsc_psi_turg_zero(ft))              &
                     / (nsc_psi_turg_full(ft) - nsc_psi_turg_zero(ft)),        &
                     0.0), 1.0)

    ! The carbon kept from growth over t days (the release is per 360
    ! days) against the room left in the store.
    c_held_back = nsc_wpot_pft(l,ft) * (1.0 - f_turg) * nsc_sink_days(ft) /    &
                  360.0
    store_room = MAX(nsc_cmax_pft(l,ft) - nsc_pft(l,ft), room_min)
    room_fill_ratio = c_held_back / store_room

    beta(l) = nsc_sink_beta_min(ft) + (1.0 - nsc_sink_beta_min(ft))            &
              * EXP( -room_fill_ratio**2 )
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE nsc_sink_factor

!#############################################################################
!#############################################################################

SUBROUTINE calc_photo_parameters( ft, land_pts, pft_photo_model, veg_pts,      &
                                  veg_index, denom, jmax_temp, jv25,           &
                                  nleaf, qtenf_term, vcmax_temp,               &
                                  jmax, rd_dark, vcmax, f_vcmax_down )

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

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: f_vcmax_down(land_pts)
    ! Soil-water capacity factor (l_som_vcmax_psi, vcmax_psi_factor):
    ! Vcmax and Jmax (Farquhar) are scaled by it after dark respiration is
    ! set, so it does not reach Rd, nor the stem and root maintenance
    ! respiration scaled from Rd.
    !
    ! MDK, 25/10/26

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

! The soil-water factor (l_som_vcmax_psi), after Rd.
!
! MDK, 25/10/26
IF ( PRESENT(f_vcmax_down) .AND. pft_photo_model == photo_farquhar ) THEN
  DO m = 1,veg_pts
    l        = veg_index(m)
    vcmax(l) = vcmax(l) * f_vcmax_down(l)
    jmax(l)  = jmax(l)  * f_vcmax_down(l)
  END DO
END IF

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
! Stop the stomatal optimisation choosing more transpiration than the
! soil can supply (l_som_supply_limit), by turning the supply into a cap
! on the stomatal conductance.
!
! stom_opt_mod costs
!
!   E_sample = dq gl pstar / (repsilon rmol T)      (mol H2O m-2 s-1),
!
! so E_sample <= E_supply is the same condition as
!
!   gl <= gl_supply = E_supply repsilon rmol T / (dq pstar).
!
! That is applied through gl_max, which every Ci search path already
! checks, so gs, A and E stay consistent: the optimiser picks the most
! profitable state the soil can supply. gl_max <= 0 means "no cap" in
! stom_opt_mod, so a zero supply is passed as TINY (no sample feasible:
! stomata close).
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE apply_supply_limit( land_pts, veg_pts, veg_index, e_supply, share,  &
                               dq, tstar, pstar, gl_max_in, gl_max_out,        &
                               r_bl, r_ca, g_o )

USE c_rmol, ONLY: rmol
USE planet_constants_mod, ONLY: repsilon
USE jules_vegetation_mod, ONLY: l_som_supply_limit, l_som_coupled_e
USE stom_opt_jls_mod, ONLY: g_of_geff
USE water_constants_mod, ONLY: Water_kg_per_mol

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

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
      ! The resistances and other-class conductance of stom_opt_mod, so that
      ! the cap is the g whose E through them is the supply: all three with
      ! l_som_coupled_e, otherwise the leaf boundary layer only.

REAL(KIND=real_jlslsm), PARAMETER :: m_h2o = Water_kg_per_mol
      ! Molar mass of water (kg mol-1).

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
    IF ( PRESENT(r_bl) ) THEN
      IF ( l_som_coupled_e ) THEN
        gl_sup = g_of_geff(gl_sup, r_bl(l), r_ca(l), g_o(l))
      ELSE
        ! The leaf boundary layer only (multilayer canopy).
        gl_sup = g_of_geff(gl_sup, r_bl(l), 0.0_real_jlslsm, 0.0_real_jlslsm)
      END IF
    END IF
    IF ( gl_max_in(l) > 0.0 ) gl_sup = MIN(gl_sup, gl_max_in(l))
    gl_max_out(l) = MAX(gl_sup, TINY(1.0_real_jlslsm))
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE apply_supply_limit

!#############################################################################
! Set up the inputs the stomatal optimisation (stom_opt_mod) of the sunlit
! or the shaded leaves in the two-leaf canopy needs to work out the
! transpiration of each stomatal conductance it tries: the humidity
! deficit the leaves see (dq_opt, kg kg-1), their boundary layer
! resistance (r_bl, s m-1), the canopy air resistance (r_ca, s m-1) and
! the conductance of the other paths into the canopy air (g_oth, m s-1).
!
! Without coupled transpiration (l_som_coupled_e), the deficit is the one
! from the leaves to the canopy air (dq_lt with l_leaf_temp, dqc
! otherwise) and the resistances and g_oth are 0.
!
! With it, the leaves lose water through the canopy air to the first
! atmospheric level. The deficit is from the leaf to the first level's
! air (qs(T_leaf) - q1 with l_leaf_temp, dq otherwise). r_bl is the
! leaves' boundary layer (l_leaf_temp only), r_ca the canopy air
! resistance, and g_oth the other leaves' stomatal conductance in series
! with their boundary layer, plus the soil's evaporation conductance g_s.
! The soil evaporates into the same canopy air at qs_soil = qs(tstar).
! Solving the canopy air humidity with these leaves, the other leaves
! (taken at these leaves' saturation humidity) and the soil gives
!
!   E = g' (D + r_ca g_s (qs_lt - qs_soil))
!       / (1 + r_ca (g' + g_other + g_s))
!
! where g' is these leaves' stomatal conductance in series with their
! boundary layer, D the deficit above, qs_lt the saturation humidity at
! these leaves' temperature and g_other the other leaves' conductance. So
! g_s is added to g_oth, and r_ca g_s (qs_lt - qs_soil) is added to the
! deficit. That term is 0 without l_leaf_temp, where the leaves are at
! tstar.
!
! MDK, 24/10/26
!#############################################################################
SUBROUTINE coupled_inputs( land_pts, veg_pts, veg_index, l_lt, dq_min, dqc,    &
                           dq, qs_lt, q1, dq_lt, gbw_this, gbw_other, r_ca_in, &
                           gl_other, g_soil, qs_soil, dq_opt, r_bl, r_ca,      &
                           g_oth )

USE jules_vegetation_mod, ONLY: l_som_coupled_e

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
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
! Give the coupled transpiration (l_som_coupled_e) of the two-leaf canopy
! its first guesses for the first pass of a timestep. Where there is a
! last timestep, they are its stomatal and boundary layer conductances of
! the sunlit and shaded leaves. Before that, the canopy conductance gc is
! shared between sunlit and shaded leaves by leaf area, and the boundary
! layer conductances are 0. (The saved canopy air resistance, slot 5 of
! cpl_state, is kept in the dump but not read.)
!
! MDK, 04/11/26
!#############################################################################
SUBROUTINE cpl_state_load( ft, land_pts, veg_pts, veg_index, lai, lai_sun,     &
                           lai_shd, gc, gl_sun, gl_shd, gbw_sun, gbw_shd )


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), lai_sun(land_pts),        &
                                      lai_shd(land_pts), gc(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: gl_sun(land_pts), gl_shd(land_pts),     &
                                       gbw_sun(land_pts), gbw_shd(land_pts)
INTEGER :: l, m
REAL(KIND=real_jlslsm) :: f_sun


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='CPL_STATE_LOAD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL sf_stom_state_alloc( land_pts )
gl_sun(:)  = 0.0
gl_shd(:)  = 0.0
gbw_sun(:) = 0.0
gbw_shd(:) = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( cpl_state(l,ft,1) >= 0.0 ) THEN
    gl_sun(l)  = cpl_state(l,ft,1)
    gl_shd(l)  = cpl_state(l,ft,2)
    gbw_sun(l) = cpl_state(l,ft,3)
    gbw_shd(l) = cpl_state(l,ft,4)
  ELSE
    f_sun = 0.5
    IF ( lai(l) > EPSILON(0.0) ) f_sun = lai_sun(l) / lai(l)
    gl_sun(l) = MAX(gc(l), 0.0) * f_sun
    gl_shd(l) = MAX(gc(l), 0.0) * (1.0 - f_sun)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cpl_state_load

!#############################################################################
! Keep this timestep's stomatal and boundary layer conductances of the
! sunlit and shaded leaves and the canopy air resistance (no lower than
! 0), so the coupled transpiration (l_som_coupled_e) of the two-leaf
! canopy can start from them next timestep (cpl_state_load).
!
! MDK, 24/10/26
!#############################################################################
SUBROUTINE cpl_state_save( ft, land_pts, veg_pts, veg_index, gl_sun, gl_shd,   &
                           gbw_sun, gbw_shd, ra_ca )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
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
