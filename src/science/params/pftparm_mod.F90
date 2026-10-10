! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds surface parameters for each Plant Functional Type (but
! not parameters that are only used by TRIFFID OR RED).


! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

MODULE pftparm

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Radiation and albedo parameters.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                        &
 orient(:)
                 ! Flag for leaf orientation: 1 for horizontal,
!                    0 for spherical.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 albsnc_max(:)                                                                 &
                 ! Snow-covered albedo for large LAI.
,albsnc_min(:)                                                                 &
                 ! Snow-covered albedo for zero LAI.
,albsnf_max(:)                                                                 &
                 ! Snow-free albedo for large LAI.
,albsnf_maxl(:)                                                                &
                 ! Min Snow-free albedo (max LAI) when scaled to obs
,albsnf_maxu(:)                                                                &
                 ! Max Snow-free albedo (max LAI) when scaled to obs
,alnir(:)                                                                      &
                 ! Leaf reflection coefficient for near infra-red.
,alnirl(:)                                                                     &
                 ! lower limit on alnir, when scaled to albedo obs
,alniru(:)                                                                     &
                 ! upper limit on alnir, when scaled to albedo obs
,alpar(:)                                                                      &
                 ! Leaf reflection coefficient for PAR.
,alparl(:)                                                                     &
                 ! lower limit on alpar, when scaled to albedo obs
,alparu(:)                                                                     &
                 ! upper limit on alpar, when scaled to albedo obs
,kext(:)                                                                       &
                 ! Light extinction coefficient - used to
!                    calculate weightings for soil and veg.
,kpar(:)                                                                       &
                 ! PAR Extinction coefficient
!                    (m2 leaf/m2 ground)
,lai_alb_lim(:)                                                                &
!                ! Lower limit on permitted LAI in albedo
,omega(:)                                                                      &
                 ! Leaf scattering coefficient for PAR.
,omegal(:)                                                                     &
                 ! lower limit on omega, when scaled to albedo obs
,omegau(:)                                                                     &
                 ! upper limit on omega, when scaled to albedo obs
,omnir(:)                                                                      &
                 ! Leaf scattering coefficient for near infra-red.
,omnirl(:)                                                                     &
                 ! lower limit on omnir, when scaled to albedo obs
,omniru(:)
                 ! upper limit on omnir, when scaled to albedo obs

!-----------------------------------------------------------------------------
! Parameters for phoyosynthesis and respiration.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                        &
 c3(:)           ! Flag for C3 types: 1 for C3 Plants,
!                    0 for C4 Plants.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 alpha(:)                                                                      &
                 ! Quantum efficiency of photosynthesis
!                    (mol CO2/mol PAR photons).
,can_struct_a(:)                                                               &
                 ! Pinty canopy structure factor corresponding to overhead
                 ! sun (dimensionless)
,dqcrit(:)                                                                     &
                 ! Critical humidity deficit (kg H2O/kg air), used with the
                 ! Jacobs closure.
,f0(:)                                                                         &
                 ! CI/CA for DQ = 0, used with the Jacobs closure.
,fd(:)                                                                         &
                 ! Dark respiration coefficient.
,g1_stomata(:)                                                                 &
                 ! Parameter g1 for the Medlyn et al. (2011) model of
                 ! stomatal conductance (kPa**0.5) - see Eqn.11 of
                 ! doi: 10.1111/j.1365-2486.2012.02790.x.
,kn(:)                                                                         &
                 ! Exponential for N profile in canopy, used with
                 ! can_rad_mod=4, 5 (decay is a function of layers).
,knl(:)                                                                        &
                 ! Decay coefficient for N profile in canopy, used with
                 ! can_rad_mod=6 and 7 (decay is a function of LAI; for the
                 ! two-leaf scheme the extkn of Wang & Leuning 1998).
,neff(:)                                                                       &
                ! Constant relating VCMAX and leaf N (mol/m2/s)
!                   from Schulze et al. 1994
!                   (AMAX = 0.4e-3 * NL  - assuming dry matter is
!                   40% carbon by mass)
!                   and Jacobs 1994:
!                   C3 : VCMAX = 2 * AMAX ;
!                   C4 : VCMAX = AMAX  ..
,nl0(:)                                                                        &
                 ! Top leaf nitrogen concentration
!                    (kg N/kg C).
,nr_nl(:)                                                                      &
                 ! Ratio of root nitrogen concentration to
!                    leaf nitrogen concentration.
,ns_nl(:)                                                                      &
                 ! Ratio of stem nitrogen concentration to
!                    leaf nitrogen concentration.
,r_grow(:)                                                                     &
                 ! Growth respiration fraction.
,tlow(:)                                                                       &
                 ! Lower temperature for photosynthesis (deg C).
,tupp(:)                                                                       &
                 ! Upper temperature for photosynthesis (deg C).
  !---------------------------------------------------------------------------
  ! Parameters for the Farquhar photosynthesis model.
  ! Jmax is the potential rate of electron transport.
  ! Vcmax is the maximum rate of carboxylation of Rubisco.
  !---------------------------------------------------------------------------
,act_jmax(:)                                                                   &
                 ! Activation energy for temperature response of Jmax
                 ! (J mol-1).
,act_vcmax(:)                                                                  &
                 ! Activation energy for temperature response of Vcmax
                 ! (J mol-1).
,alpha_elec(:)                                                                 &
                 ! Quantum yield of electron transport
                 ! (mol electrons/mol PAR photons).
,deact_jmax(:)                                                                 &
                 ! Deactivation energy for temperature response of Jmax
                 ! (J mol-1). This describes the rate of decrease
                 ! above the optimum temperature.
,deact_vcmax(:)                                                                &
                 ! Deactivation energy for temperature response of Vcmax
                 ! and Jmax (J mol-1). This describes the rate of decrease
                 ! above the optimum temperature.
,ds_jmax(:)                                                                    &
                 ! Entropy factor for temperature reponse of Jmax
                 ! (J mol-1 K-1).
,ds_vcmax(:)                                                                   &
                 ! Entropy factor for temperature reponse of Vcmax
                 ! (J mol-1 K-1).
,jv25_ratio(:)
                 ! Ratio of Jmax to Vcmax at 25 deg C
                 ! (mol electrons mol-1 CO2).

!-----------------------------------------------------------------------------
! Parameters for trait physiology
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
hw_sw(:)                                                                       &
                 ! Heart:Stemwood Ratio (kg N/kg N)
,lma(:)                                                                        &
                 ! Leaf mass by area (1/SLA), (kg leaf/ m2)
,nmass(:)                                                                      &
                 ! Leaf nitrogen dry weight (g N/g leaf)
,nr(:)                                                                         &
                 ! Root nitrogen concentration (kg N/kg C)
,nsw(:)                                                                        &
                 ! Stem nitrogen concentration (kg N/kg C)
,q10_leaf(:)                                                                   &
                 ! Factor for leaf respiration.
,rmass(:)                                                                      &
                 ! Root carbon dry weight (kg C/kg root) JBaguley
,vint(:)                                                                       &
                 ! Y intercept of the Narea to Vcmax relationship
                 ! from Kattge et al. (2009)
,vsl(:)
                 ! Slope of the Narea to Vcmax relationship
                 ! from Kattge et al. (2009)

!-----------------------------------------------------------------------------
! Allometric and other parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 a_wl(:)                                                                       &
                 ! Allometric coefficient relating the target
!                    woody biomass to the leaf area index
!                    (kg C/m2)
,a_ws(:)                                                                       &
                 ! Woody biomass as a multiple of live
!                    stem biomass.
,b_wl(:)                                                                       &
                 ! Allometric exponent relating the target
!                    woody biomass to the leaf area index.
,eta_sl(:)                                                                     &
                 ! Live stemwood coefficient (kg C m-1 (m2 leaf)-1)
,sigl(:)
                 ! Specific density of leaf carbon
!                    (kg C/m2 leaf).

!-----------------------------------------------------------------------------
! Phenology parameters.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 dgl_dm(:)                                                                     &
                 ! Rate of change of leaf turnover rate with
!                    moisture availability.
,dgl_dt(:)                                                                     &
                 ! Rate of change of leaf turnover rate with
!                    temperature (/K)
,fsmc_of(:)                                                                    &
                 ! Moisture availability below which leaves
!                    are dropped.
,g_leaf_0(:)                                                                   &
                 ! Minimum turnover rate for leaves (/360days).

,tleaf_of(:)
                 ! Temperature below which leaves are
!                    dropped (K)

!-----------------------------------------------------------------------------
! Parameters for hydrological, thermal and other "physical" characteristics.
!-----------------------------------------------------------------------------
INTEGER, ALLOCATABLE ::                                                        &
 fsmc_mod(:)
                   ! Flag for whether water stress is calculated from,
                   ! available water in layers weighted by root fraction (0),
                   ! available water in root zone (1),
                   ! or
                   ! soil to root resistivity and water potential (2).

LOGICAL, ALLOCATABLE ::                                                        &
 calc_rz_psi(:)
                   ! Flag for wheather the rootzone water potential is
                   ! calculated.

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 catch0(:)                                                                     &
                 ! Minimum canopy capacity (kg/m2).
,dcatch_dlai(:)                                                                &
                 ! Rate of change of canopy capacity with LAI.
,dust_veg_scj(:)                                                               &
                 ! Dust emission scaling factor for  each PFT
,dz0v_dh(:)                                                                    &
                 ! Rate of change of vegetation roughness
!                    length with height.
,emis_pft(:)                                                                   &
                 !  Surface emissivity
,fsmc_p0(:)                                                                    &
                 ! parameter in calculation of the
                 ! soil moisture at which the plant begins to experience
                 ! water stress
,glmin(:)                                                                      &
                 ! Minimum leaf conductance for H2O (m/s).
,gsoil_f(:)                                                                    &
                 ! Soil evaporation enhancement factor (no units).
,infil_f(:)                                                                    &
                 ! Infiltration enhancement factor.
,min_gl_pft(:)                                                                 & ! JBaguley
                 ! Minimum conductance to H2O (m/s) given to closed stomata
                 ! by the stomatal optimisation: a numerical floor, applied
                 ! as a canopy value in every canopy scheme (not scaled by
                 ! LAI). Physical night-time/cuticular loss is gcuticular
                 ! (l_som_cuticular_floor), per unit leaf area.
,min_rootc_pft(:)                                                              & ! JBaguley
                 ! Minimum root C mass per unit area (kg m-2) for each pft.
                 ! Used when calculating water stress using root resistivity
                 ! and water potential (fsmc_mod(i) == 2).
,psi_close(:)                                                                  &
                 ! soil matric potential (Pa) below which soil moisture
                 ! stress factor fsmc is zero. Should be negative.
,psi_open(:)                                                                   &
                 ! soil matric potential (Pa) above which soil moisture
                 ! stress factor fsmc is one. Should be negative.
,root_psi_crit(:)                                                              & ! JBaguley
                 ! Negative critical root water potential (Pa) above which
                 ! roots detach from soil.
,root_radi_pft(:)                                                              & ! JBaguley
                 ! fine root radius (m)
                 ! Used when calculating water stress using root resistivity
                 ! and water potential (fsmc_mod(i) == 2).
,rootc_density_pft(:)                                                          & ! JBaguley
                 ! root density, (kg m-3) specificaly root mas per unit root
                 ! volume.
                 ! Used when calculating water stress using root resistivity
                 ! and water potential (fsmc_mod(i) == 2).
,rootd_ft(:)                                                                   &
                 ! e-folding depth (m) of the root density.
,z0v(:)                                                                        &
                 ! Specified vegetation roughness length.
!                    used with l_spec_veg_z0 = .true.
,pft_big_leaf_corection_nitrogen_reduction_factor(:)                           &
                 ! Determins shape of nitrogen reduction when moving from the
                 ! top of the canopy to the bottom as a function of leaf area
                 ! index abovethe current height.
,pft_big_leaf_corection_light_atenuation_factor(:)
                 ! Determines the shape of the light attenuation as a
                 ! function of leaf area index above the current height.
                 ! Used to correct the big leaf model for the effect of
                 ! nitrogen reduction on light attenuation.



!-----------------------------------------------------------------------------
! Parameters for ozone damage
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 dfp_dcuo(:)                                                                   &
                 ! Plant type specific O3 sensitivity parameter
                 ! (nmol-1 m2 s).
,fl_o3_ct(:)
                 ! Critical flux of O3 to vegetation (nmol/m2/s).

!-----------------------------------------------------------------------------
! Parameters for BVOC emissions
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 aef(:)                                                                        &
                 ! Acetone Emission Factor (ugC/g/h)
,ci_st(:)                                                                      &
                 ! Internal CO2 partial pressure (Pa)
                 !   at standard conditions
,gpp_st(:)                                                                     &
                 ! Gross primary productivity (KgC/m2/s)
                 !   at standard conditions
,ief(:)                                                                        &
                 ! Isoprene Emission Factor (ugC/g/h)
                 ! See Pacifico et al., (2011) Atm. Chem. Phys.
,mef(:)                                                                        &
                 ! Methanol Emission Factor (ugC/g/h)
,tef(:)
                 ! (Mono-)Terpene Emission Factor (ugC/g/h)

!-----------------------------------------------------------------------------
! Parameters for INFERNO combustion
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 avg_ba(:)                                                                     &
                 ! Average PFT Burnt Area per fire
,ccleaf_min(:)                                                                 &
                 ! Leaf minimum combustion completeness (kg/kg)
,ccleaf_max(:)                                                                 &
                 ! Leaf maximum combustion completeness (kg/kg)
,ccwood_min(:)                                                                 &
                 ! Wood (or Stem) minimum combustion completeness (kg/kg)
,ccwood_max(:)                                                                 &
                 ! Wood (or Stem) maximum combustion completeness (kg/kg)
,fire_mort(:)
                 ! Fire mortality per PFT

!-----------------------------------------------------------------------------
! Parameters for INFERNO emissions
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 fef_bc(:)                                                                     &
                 ! Fire BC Emission Factor (g/kg)
,fef_ch4(:)                                                                    &
                 ! Fire CH4 Emission Factor (g/kg)
,fef_co(:)                                                                     &
                 ! Fire CO Emission Factor (g/kg)
,fef_co2(:)                                                                    &
                 ! Fire CO2 Emission Factor (g/kg)
                 ! See Thonicke et al., (2005,2010)
,fef_nox(:)                                                                    &
                 ! Fire NOx Emission Factor (g/kg)
,fef_oc(:)                                                                     &
                 ! Fire OC Emission Factor (g/kg)
,fef_so2(:)                                                                    &
                 ! Fire SO2 Emission Factor (g/kg)
,fef_c2h4(:)                                                                   &
                 ! Fire C2H4 Emission Factor (g/kg)
,fef_c2h6(:)                                                                   &
                 ! Fire C2H6 Emission Factor (g/kg)
,fef_c3h8(:)                                                                   &
                 ! Fire C2H8 Emission Factor (g/kg)
,fef_hcho(:)                                                                   &
                 ! Fire HCHO Emission Factor (g/kg)
,fef_mecho(:)                                                                  &
                 ! Fire MeCHO Emission Factor (g/kg)
,fef_nh3(:)                                                                    &
                 ! Fire NH3 Emission Factor (g/kg)
,fef_dms(:)
                 ! Fire DMS Emission Factor (g/kg)

!-----------------------------------------------------------------------------
! Parameters for SUGAR
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 sug_grec(:)                                                                   &
                 ! Turnover of structural carbon into NSC (KgC/m2/s)
,sug_g0(:)                                                                     &
                 ! Specific structural C growth rate (KgC/m2/s)
,sug_yg(:)
                 ! Growth yield fraction

!-----------------------------------------------------------------------------
! Parameters for stomatal optimisation model (som)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 leaf_crit(:)
                 ! Critical value of base som parameter. Used to determin
                 ! the range over which the stomatal optimisation model samples
                 ! the primary leaf parameter.
                 ! Physical meaning is ether:
                 !    -  Minimum possible leaf intercellular carbon (Pa).
                 !       for each set of ci sample points the smallest
                 !       value of ci is the maximum of leaf_crit and the
                 !       CO2 compensation point.
                 !         som_base_parm = 1
                 !    -  leaf critical water potential (Pa),
                 !         som_base_parm = 2
!-----------------------------------------------------------------------------
! Parameters for SOX
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 sox_a(:)                                                                      &
                 ! The shape parameter in the xylem vulnerability curve.
,sox_p50(:)                                                                    &
                 ! Xlem water potential at which xylem hydraulic
                 ! conductance is half its maximum value. (MPa)
,sox_rp_min(:)
                 ! Plant minimum hydraulic resistance. (m2 s MPa/mol)

INTEGER, ALLOCATABLE ::                                                        &
pft_conductance_model(:),                                                      &
                 ! Flag for the xylem conductance model used by the stomatal
                 !  optimisation model.
                 !      1: Cumulative Weibul distribution
                 !           k(psi) = kmax * exp((psi/b)^c)
                 !      2: SOX model
                 !           k(psi) = kmax / (1 + (psi/b)^c)
pft_xylem_impairment_model(:),                                                 &
                 ! Flag for the xylem impairment model used by the stomatal
                 !  optimisation model.
                 !       0: No xylem impairment
                 !       1: Reduce kmax
                 !       2: Whole trunk (historic minimum water potentials)
                 !       3: Embolism memory with recovery
ximpair_psi_driver(:),                                                         &
                 ! Water potential driving embolism in the memory impairment
                 !  model (pft_xylem_impairment_model = 3).
                 !       1: Leaf water potential
                 !       2: Mean of leaf and root zone water potentials
                 !       3: Root zone water potential
                 !       4: Stem water potential (outlet of the stem
                 !          segment; needs l_som_plant_segments)
ximpair_reset_mmdd(:),                                                         &
                 ! Month and day (month*100 + day, e.g. 401 = 1 April) on
                 !  which the memory impairment model's embolism is reset
                 !  (e.g. new earlywood in ring-porous species). 0 disables.
ximpair_growth_basis(:)
                 ! Carbon supply for the growth recovery term
                 !  (l_ximpair_rec_growth) of the memory impairment models:
                 !       1: ximpair_wood_alloc * canopy net photosynthesis
                 !       2: ximpair_wood_alloc * NPP (previous timestep)
                 !       3: TRIFFID gross wood production, from the change in
                 !          allometric wood carbon with canopy height plus
                 !          wood turnover (g_wood); needs l_triffid

REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 kmax_pft(:)                                                                   &
                 ! Maximum xylem conductance (mol m-2 s-1 Pa-1).
,P50(:)                                                                        &
                 ! Water potential at which 50% of the xylem conductance is
                 ! lost (Pa).
,P88(:)                                                                        &
                 ! Water potential at which 88% of the xylem conductance is
                 ! lost (Pa).
,conductance_b_pft(:)                                                              &
                 ! Sensetivity parameter, b, in the xylem conductance model
                 !  (Pa).
                 ! NOTE: This value is not directly input by the user, instead
                 !        it is calculated from the P50 and P88 values in
                 !        ptftparm_io_mod.F90.
,conductance_c_pft(:)                                                              &
                 ! Shape parameter, c, in the xylem conductance model.
                 ! NOTE: This value is not directly input by the user, instead
                 !        it is calculated from the P50 and P88 values in
                 !        ptftparm_io_mod.F90.
,seg_kfac(:,:)                                                                 &
                 ! Root / stem / leaf segments (l_som_plant_segments): segment
                 ! maximum conductance as a multiple of the whole-plant kmax,
                 ! sum(seg_frac) / seg_frac(s), so the three in series give
                 ! kmax when well watered. (npft, 3)
,conductance_b_seg(:,:)                                                        &
,conductance_c_seg(:,:)                                                        &
                 ! Cumulative Weibull b (Pa) and c of each segment, from the
                 ! segment P50/P88 (default: the PFT's P50/P88). (npft, 3)
,gcuticular(:)                                                                 &
                 ! Cuticular (minimum) leaf conductance to water vapour, per
                 ! unit leaf area (mmol H2O m-2 s-1), applied as a floor on
                 ! the canopy conductance when l_som_cuticular_floor.
,psi_nsl_onset(:)                                                              &
                 ! Leaf water potential (Pa) below which the nonstomatal
                 ! limitation starts (l_som_nsl). 0 (default) gives Dewar et
                 ! al. (2022) Eqn 3(b); the turgor loss point is an optional
                 ! variant.
,fsmc_q(:)                                                                     &
                 ! Curvature exponent q of the soil moisture stress factor
                 ! of each layer, fsmc = ((x - x_close)/(x_open - x_close))^q
                 ! (fsmc_layer); 1 (default) is the standard linear factor.
,psi_nsl0(:)                                                                   &
                 ! Leaf water potential (Pa) at which the nonstomatal
                 ! limitation reduces photosynthesis to zero (l_som_nsl);
                 ! psi_0 of Dewar et al. (2022). Must be < psi_nsl_onset.
,cmax_a(:)                                                                     &
                 ! CMax stomata (stomata_model = 6; Wolf et al. 2016,
                 ! Anderegg et al. 2018): curvature a of the carbon cost of
                 ! low leaf water potential, Theta = a/2 psi^2 + b |psi|
                 ! (umol CO2 m-2 leaf s-1 MPa-2). No default: must be set
                 ! (> 0) for stomata_model = 6.
,cmax_b(:)                                                                     &
                 ! CMax: linear term b of Theta (umol CO2 m-2 s-1 MPa-1).
                 ! Default 0: Sabot et al. (2022, JAMES) found it of low
                 ! influence and hard to constrain.
,cgain_varpi(:)                                                                &
                 ! CGain stomata (stomata_model = 7; Lu et al. 2020): carbon
                 ! cost varpi of losing the whole of the path conductance
                 ! (umol CO2 m-2 leaf s-1). No default: must be set (> 0)
                 ! for stomata_model = 7.
,sl_cica_well_watered(:)                                                       &
                 ! Supply-loss stomata (stomata_model = 9): ci/ca of the
                 ! carbon demand when water is not limiting,
                 ! ci = sl_cica_well_watered ca (Medlyn with D fixed at
                 ! 1 kPa: ci/ca = g1/(g1 + 1)). In (0, 1); default 0.8.
,psi_vcmax_f(:)                                                                &
                 ! Root-zone water potential (Pa) at which photosynthetic
                 ! capacity is about halved (l_som_vcmax_psi; psi_f of Zhou
                 ! et al. 2013).
,sf_vcmax(:)                                                                   &
                 ! Steepness of that down-regulation (MPa-1).
,or_z0soil_fac(:)                                                              &
                 ! Multiplier on the Or scheme's soil-surface roughness z0soil
                 ! beneath this PFT (l_soil_evap_or; CABLE: 0.01 min(1, LAI) +
                 ! 0.02 min(u*2/g, 1) m). In wet soil the Or resistance is
                 ! ~z0soil / D_vapour (~400 s m-1 at 1 cm), so this sets the
                 ! wet-soil evaporation rate. 1 = CABLE.
,soil_litter_depth(:)                                                          &
                 ! Litter layer depth over the soil (m) for the Or soil
                 ! evaporation scheme (l_soil_evap_or): a vapour diffusion
                 ! resistance depth / DvLitt in series with the soil
                 ! resistance, as CABLE's default-scheme litter (relitt; CABLE
                 ! depth = clitt * 0.003 m per t C ha-1). 0 = none.
,psi_vcmax_fmin(:)                                                             &
                 ! Floor of that down-regulation (-): f = fmin + (1 - fmin)
                 ! f_Zhou, so capacity never falls below fmin (0 = none, as
                 ! in De Kauwe et al. 2015; the floor is an adjustment).
,g1_tuzet(:)                                                                   &
                 ! DESICA (stomata_model = 8): slope of gs = g1 fw An / ca (-).
,sf_tuzet(:)                                                                   &
                 ! DESICA: sensitivity of the Tuzet closure (MPa-1).
,psi_f_tuzet(:)                                                                &
                 ! DESICA: reference leaf water potential of the Tuzet
                 ! closure, ~50% closure (Pa).
,cap_leaf(:)                                                                   &
,cap_stem(:)                                                                   &
                 ! DESICA: leaf and stem water capacitance per unit leaf area
                 ! (mol H2O m-2 leaf Pa-1).
,kcrit_fractional_loss(:)                                                      &
                 ! Critical fractional loss of xylem conductance.
,kcrit(:)                                                                      &
                 ! Critical xylem conductance (mmol m-2 s-1 MPa-1).
                 ! NOTE: This value is not directly input by the user, instead
                 !        it is calculated from kmax_pft and
                 !        kcrit_fractional_loss in ptftparm_io_mod.F90.
                 !         kcrit = kmax_pft * (1-kcrit_fractional_loss)
,psi_crit(:)                                                                  &
                 ! Critical water potential (Pa).
                 ! NOTE: This value is not directly input by the user, instead
                 !        it is calculated from the kcrit and the conductance
                 !        model in ptftparm_io_mod.F90.
,ximpair_tau_rec(:)                                                            &
                 ! Recovery timescale of embolism in the memory impairment
                 !  model (days). <= 0 disables recovery.
,ximpair_psi_refill(:)                                                         &
                 ! Water potential of the damage driver above which embolism
                 !  recovers in the memory impairment model (Pa).
,ximpair_wood_alloc(:)                                                          &
                 ! Fraction of canopy net assimilation that goes into new
                 !  conducting xylem, for recovery from xylem impairment with
                 !  l_ximpair_rec_growth.
,ximpair_leaf_sens(:)                                                          &
                 ! Sensitivity of the canopy to lasting xylem damage
                 !  (l_ximpair_leaf_loss): the phenological state is capped
                 !  at 1 - ximpair_leaf_sens * (1 - k_cap/kmax). 1 keeps leaf
                 !  area in proportion to the conducting capacity, 0 disables.
,ximpair_psi_growth(:)                                                         &
                 ! Root-zone (~predawn) water potential below which no new
                 !  conducting xylem grows, for the growth recovery term
                 !  (l_ximpair_rec_growth; Pa). Stem growth at Puechabon stops
                 !  near predawn -1.1 MPa (Lempereur et al. 2015, New Phytol
                 !  207: 579), so recovery falls in the growth windows (spring
                 !  and autumn), not in summer drought. Default -1e30: no gate.
,ximpair_rec_years(:)                                                         &
                 ! Years of typical growth for the growth recovery term
                 !  (l_ximpair_rec_growth)
                 !  to recover the loss of conductivity: the loss falls
                 !  exponentially with growth, e-folding over
                 !  ximpair_rec_years / 2 years of typical growth. <= 0 applies
                 !  the renewed fraction directly (fast recovery).
,ximpair_tau_stem(:)                                                          &
                 ! Per-segment memory (l_ximpair_seg_memory): e-folding time
                 !  (years; of typical growth with the growth clock, else of
                 !  calendar time) of the stem segment's loss of
                 !  conductivity, i.e. the residence time of conducting
                 !  sapwood (sapwood -> heartwood turnover; e.g. LPJ-GUESS
                 !  turnover_sap 0.05-0.1 yr-1, 0.075 for Q. ilex).
                 !  <= 0: no stem recovery.
,ximpair_tau_leaf(:)
                 ! As ximpair_tau_stem for the leaf (leaf and twig) segment,
                 !  renewed with the leaves: about the leaf lifespan.
                 !  <= 0: no leaf recovery.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PFTPARM'

CONTAINS

SUBROUTINE pftparm_alloc(npft)

USE missing_data_mod, ONLY: imdi, rmdi

!No USE statements other than Dr Hook
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: npft

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PFTPARM_ALLOC'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  ====pftparm module common====

! Radiation and albedo parameters.
ALLOCATE( albsnc_max(npft))
ALLOCATE( albsnc_min(npft))
ALLOCATE( albsnf_max(npft))
ALLOCATE( albsnf_maxl(npft))
ALLOCATE( albsnf_maxu(npft))
ALLOCATE( alnir(npft))
ALLOCATE( alnirl(npft))
ALLOCATE( alniru(npft))
ALLOCATE( alpar(npft))
ALLOCATE( alparl(npft))
ALLOCATE( alparu(npft))
ALLOCATE( kext(npft))
ALLOCATE( kpar(npft))
ALLOCATE( lai_alb_lim(npft))
ALLOCATE( omega(npft))
ALLOCATE( omegal(npft))
ALLOCATE( omegau(npft))
ALLOCATE( omnir(npft))
ALLOCATE( omnirl(npft))
ALLOCATE( omniru(npft))
ALLOCATE( orient(npft))

orient(:)       = imdi
albsnc_max(:)   = rmdi
albsnc_min(:)   = rmdi
albsnf_max(:)   = rmdi
albsnf_maxl(:)  = rmdi
albsnf_maxu(:)  = rmdi
alnir(:)        = rmdi
alnirl(:)       = rmdi
alniru(:)       = rmdi
alpar(:)        = rmdi
alparl(:)       = rmdi
alparu(:)       = rmdi
kext(:)         = rmdi
kpar(:)         = rmdi
lai_alb_lim(:)  = rmdi
omega(:)        = rmdi
omegal(:)       = rmdi
omegau(:)       = rmdi
omnir(:)        = rmdi
omnirl(:)       = rmdi
omniru(:)       = rmdi

! Photosynthesis and respiration parameters
ALLOCATE( act_jmax(npft))
ALLOCATE( act_vcmax(npft))
ALLOCATE( alpha(npft))
ALLOCATE( alpha_elec(npft))
ALLOCATE( c3(npft))
ALLOCATE( can_struct_a(npft))
ALLOCATE( deact_jmax(npft))
ALLOCATE( deact_vcmax(npft))
ALLOCATE( dqcrit(npft))
ALLOCATE( ds_jmax(npft))
ALLOCATE( ds_vcmax(npft))
ALLOCATE( f0(npft))
ALLOCATE( fd(npft))
ALLOCATE( g1_stomata(npft))
ALLOCATE( jv25_ratio(npft))
ALLOCATE( kn(npft))
ALLOCATE( knl(npft))
ALLOCATE( neff(npft))
ALLOCATE( nl0(npft))
ALLOCATE( nr_nl(npft))
ALLOCATE( ns_nl(npft))
ALLOCATE( r_grow(npft))
ALLOCATE( tlow(npft))
ALLOCATE( tupp(npft))

c3(:)           = imdi
act_jmax(:)     = rmdi
act_vcmax(:)    = rmdi
alpha(:)        = rmdi
alpha_elec(:)   = rmdi
can_struct_a(:) = rmdi
deact_jmax(:)   = rmdi
deact_vcmax(:)  = rmdi
dqcrit(:)       = rmdi
ds_jmax(:)      = rmdi
ds_vcmax(:)     = rmdi
f0(:)           = rmdi
fd(:)           = rmdi
g1_stomata(:)   = rmdi
jv25_ratio(:)   = rmdi
kn(:)           = rmdi
knl(:)          = rmdi
neff(:)         = rmdi
nl0(:)          = rmdi
nr_nl(:)        = rmdi
ns_nl(:)        = rmdi
r_grow(:)       = rmdi
tlow(:)         = rmdi
tupp(:)         = rmdi

! Traint physiology parameters
ALLOCATE( hw_sw(npft))
ALLOCATE( lma(npft))
ALLOCATE( nmass(npft))
ALLOCATE( nr(npft))
ALLOCATE( nsw(npft))
ALLOCATE( q10_leaf(npft))
ALLOCATE( rmass(npft))
ALLOCATE( vint(npft))
ALLOCATE( vsl(npft))

hw_sw(:)        = rmdi
lma(:)          = rmdi
nmass(:)        = rmdi
nr(:)           = rmdi
nsw(:)          = rmdi
q10_leaf(:)     = rmdi
rmass(:)        = 0.49 !JBaguley
vint(:)         = rmdi
vsl(:)          = rmdi

! Allometric parameters
ALLOCATE( a_wl(npft))
ALLOCATE( a_ws(npft))
ALLOCATE( b_wl(npft))
ALLOCATE( eta_sl(npft))
ALLOCATE( sigl(npft))

a_wl(:)         = rmdi
a_ws(:)         = rmdi
b_wl(:)         = rmdi
eta_sl(:)       = rmdi
sigl(:)         = rmdi

! Phenology parameters
ALLOCATE( dgl_dm(npft))
ALLOCATE( dgl_dt(npft))
ALLOCATE( fsmc_of(npft))
ALLOCATE( g_leaf_0(npft))
ALLOCATE( tleaf_of(npft))

dgl_dm(:)       = rmdi
dgl_dt(:)       = rmdi
fsmc_of(:)      = rmdi
g_leaf_0(:)     = rmdi
tleaf_of(:)     = rmdi

! Hydrological parameters
ALLOCATE( calc_rz_psi(npft)) ! JBaguley
ALLOCATE( catch0(npft))
ALLOCATE( dcatch_dlai(npft))
ALLOCATE( dust_veg_scj(npft))
ALLOCATE( dz0v_dh(npft))
ALLOCATE( emis_pft(npft))
ALLOCATE( fsmc_mod(npft))
ALLOCATE( fsmc_p0(npft))
ALLOCATE( glmin(npft))
ALLOCATE( gsoil_f(npft))
ALLOCATE( infil_f(npft))
ALLOCATE( min_rootc_pft(npft)) ! JBaguley
ALLOCATE( min_gl_pft(npft)) ! JBaguley
ALLOCATE( psi_close(npft))
ALLOCATE( psi_open(npft))
ALLOCATE( root_psi_crit(npft))  ! JBaguley
ALLOCATE( root_radi_pft(npft))  ! JBaguley
ALLOCATE( rootc_density_pft(npft))  ! JBaguley
ALLOCATE( rootd_ft(npft))
ALLOCATE( z0v(npft))
ALLOCATE( pft_big_leaf_corection_nitrogen_reduction_factor(npft))

fsmc_mod(:)     = imdi
catch0(:)       = rmdi
dcatch_dlai(:)  = rmdi
dust_veg_scj(:) = rmdi
dz0v_dh(:)      = rmdi
emis_pft(:)     = rmdi
fsmc_p0(:)      = rmdi
glmin(:)        = rmdi
gsoil_f(:)      = rmdi
infil_f(:)      = rmdi
psi_close(:)    = rmdi
psi_open(:)     = rmdi
rootd_ft(:)     = rmdi
z0v(:)          = rmdi
calc_rz_psi(:)       = .FALSE. ! JBaguley
min_gl_pft(:)        = 1.0e-9 ! numerical floor, see above
min_rootc_pft(:)     = 1.0 ! JBaguley M.Williams etal 2001
root_psi_crit(:)     =-0.1e6 ! JBaguley
root_radi_pft(:)     = 0.29e-3 ! Bonan et al. 2014 (was 0.0005, Williams 2001)
rootc_density_pft(:) = 0.31e3  ! Bonan et al. 2014 (was 0.5e3, Williams 2001)
pft_big_leaf_corection_nitrogen_reduction_factor(:) = 1.0

! Ozone damage parameters
ALLOCATE( dfp_dcuo(npft))
ALLOCATE( fl_o3_ct(npft))

dfp_dcuo(:) = rmdi
fl_o3_ct(:) = rmdi

! BVOC emission parameters
ALLOCATE( aef(npft))
ALLOCATE( ci_st(npft))
ALLOCATE( gpp_st(npft))
ALLOCATE( ief(npft))
ALLOCATE( mef(npft))
ALLOCATE( tef(npft))

aef(:)    = rmdi
ci_st(:)  = rmdi
gpp_st(:) = rmdi
ief(:)    = rmdi
mef(:)    = rmdi
tef(:)    = rmdi

! INFERNO combustion parameters
ALLOCATE( avg_ba(npft))
ALLOCATE( ccleaf_min(npft))
ALLOCATE( ccleaf_max(npft))
ALLOCATE( ccwood_min(npft))
ALLOCATE( ccwood_max(npft))
ALLOCATE( fire_mort(npft))

avg_ba(:)     = rmdi
ccleaf_min(:) = rmdi
ccleaf_max(:) = rmdi
ccwood_min(:) = rmdi
ccwood_max(:) = rmdi
fire_mort(:)  = rmdi

! INFERNO emission parameters
ALLOCATE( fef_bc(npft))
ALLOCATE( fef_ch4(npft))
ALLOCATE( fef_co(npft))
ALLOCATE( fef_co2(npft))
ALLOCATE( fef_nox(npft))
ALLOCATE( fef_oc(npft))
ALLOCATE( fef_so2(npft))
ALLOCATE( fef_c2h4(npft))
ALLOCATE( fef_c2h6(npft))
ALLOCATE( fef_c3h8(npft))
ALLOCATE( fef_mecho(npft))
ALLOCATE( fef_hcho(npft))
ALLOCATE( fef_nh3(npft))
ALLOCATE( fef_dms(npft))

fef_bc(:)   = rmdi
fef_ch4(:)  = rmdi
fef_co(:)   = rmdi
fef_co2(:)  = rmdi
fef_nox(:)  = rmdi
fef_oc(:)   = rmdi
fef_so2(:)  = rmdi
fef_c2h4(:) = rmdi
fef_c2h6(:) = rmdi
fef_c3h8(:) = rmdi
fef_mecho(:)= rmdi
fef_hcho(:) = rmdi
fef_nh3(:)  = rmdi
fef_dms(:)  = rmdi

! SUGAR parameters
ALLOCATE( sug_grec(npft))
ALLOCATE( sug_g0(npft))
ALLOCATE( sug_yg(npft))

sug_grec(:) = rmdi
sug_g0(:)   = rmdi
sug_yg(:)   = rmdi

! SOM parameters
ALLOCATE( leaf_crit(npft))
ALLOCATE( pft_conductance_model(npft))
ALLOCATE( pft_xylem_impairment_model(npft))
ALLOCATE( kcrit_fractional_loss(npft))
ALLOCATE( kcrit(npft))
ALLOCATE( psi_crit(npft))
ALLOCATE( kmax_pft(npft))
ALLOCATE( P50(npft))
ALLOCATE( P88(npft))
ALLOCATE( conductance_b_pft(npft))
ALLOCATE( conductance_c_pft(npft))
ALLOCATE( ximpair_psi_driver(npft))
ALLOCATE( seg_kfac(npft,3))
ALLOCATE( gcuticular(npft))
ALLOCATE( psi_nsl_onset(npft))
ALLOCATE( psi_nsl0(npft))
ALLOCATE( cmax_a(npft))
ALLOCATE( cmax_b(npft))
ALLOCATE( cgain_varpi(npft))
ALLOCATE( sl_cica_well_watered(npft))
ALLOCATE( fsmc_q(npft))
ALLOCATE( psi_vcmax_f(npft))
ALLOCATE( sf_vcmax(npft))
ALLOCATE( psi_vcmax_fmin(npft))
ALLOCATE( soil_litter_depth(npft))
ALLOCATE( or_z0soil_fac(npft))
ALLOCATE( g1_tuzet(npft))
ALLOCATE( sf_tuzet(npft))
ALLOCATE( psi_f_tuzet(npft))
ALLOCATE( cap_leaf(npft))
ALLOCATE( cap_stem(npft))
ALLOCATE( conductance_b_seg(npft,3))
ALLOCATE( conductance_c_seg(npft,3))
ALLOCATE( ximpair_reset_mmdd(npft))
ALLOCATE( ximpair_growth_basis(npft))
ALLOCATE( ximpair_tau_rec(npft))
ALLOCATE( ximpair_psi_refill(npft))
ALLOCATE( ximpair_wood_alloc(npft))
ALLOCATE( ximpair_leaf_sens(npft))
ALLOCATE( ximpair_psi_growth(npft))
ALLOCATE( ximpair_rec_years(npft))
ALLOCATE( ximpair_tau_stem(npft))
ALLOCATE( ximpair_tau_leaf(npft))

leaf_crit(:) = 0.0
pft_conductance_model(:) = 0
pft_xylem_impairment_model(:) = 0
kcrit_fractional_loss(:) = 0.95
! NOTE: kcrit is calculated in ptftparm_io_mod.F90 using
!        kmax_pft * kcrit_fractional_loss.
kcrit(:) = 0.0
! NOTE: psi_crit is calculated in ptftparm_io_mod.F90
psi_crit(:) = 0.0
kmax_pft(:) = 0.0
P50(:) = 0.0
P88(:) = 0.0
! NOTE: conductance_b_pft and conductance_c_pft are calculated in ptftparm_io_mod.F90
!        using the P50 and P88 values and the choice of conductance model.
conductance_b_pft(:) = 1.0
conductance_c_pft(:) = 1.0
ximpair_psi_driver(:) = 1
seg_kfac(:,:) = 1.0
gcuticular(:) = 3.0
psi_nsl_onset(:) = 0.0
psi_nsl0(:) = -3.0e6
cmax_a(:) = 0.0
cmax_b(:) = 0.0
cgain_varpi(:) = 0.0
sl_cica_well_watered(:) = 0.8
fsmc_q(:) = 1.0
psi_vcmax_f(:) = -2.0e6
sf_vcmax(:) = 2.0
psi_vcmax_fmin(:) = 0.0
soil_litter_depth(:) = 0.0
or_z0soil_fac(:) = 1.0
g1_tuzet(:) = 4.19
sf_tuzet(:) = 2.0
psi_f_tuzet(:) = -2.05e6
cap_leaf(:) = 83.3e-9
cap_stem(:) = 3000.0e-9
conductance_b_seg(:,:) = 1.0
conductance_c_seg(:,:) = 1.0
ximpair_reset_mmdd(:) = 0
ximpair_growth_basis(:) = 1
ximpair_tau_rec(:) = 0.0
ximpair_psi_refill(:) = -0.5e6
ximpair_wood_alloc(:) = 0.25
ximpair_leaf_sens(:) = 0.0
ximpair_psi_growth(:) = -1.0e30
ximpair_rec_years(:) = 0.0
ximpair_tau_stem(:) = 0.0
ximpair_tau_leaf(:) = 0.0

! SOX parameters
ALLOCATE( sox_a(npft))
ALLOCATE( sox_p50(npft))
ALLOCATE( sox_rp_min(npft))

sox_a(:)      = rmdi
sox_p50(:)    = rmdi
sox_rp_min(:) = rmdi

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pftparm_alloc


SUBROUTINE print_nlist_jules_pftparm()

USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('pftparm',                                                    &
    'Contents of namelist jules_pftparm')

#if !defined(UM_JULES)
WRITE(lineBuffer,*)' fsmc_mod = ',fsmc_mod
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_close = ',psi_close
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_open = ',psi_open
CALL jules_print('pftparm',lineBuffer)
#endif


WRITE(lineBuffer,*)' a_wl = ',a_wl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' a_ws = ',a_ws
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' act_jmax = ',act_jmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' act_vcmax = ',act_vcmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' aef = ',aef
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' albsnc_max = ',albsnc_max
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' albsnc_min = ',albsnc_min
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' albsnf_max = ',albsnf_max
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' alpha = ',alpha
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' alpha_elec = ',alpha_elec
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' alnir = ',alnir
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' alpar = ',alpar
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' avg_ba = ',avg_ba
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' b_wl = ',b_wl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' c3 = ',c3
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' can_struct_a = ',can_struct_a
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' catch0 = ',catch0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_max = ',ccleaf_max
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ccleaf_min = ',ccleaf_min
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ccwood_max = ',ccwood_max
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ccwood_min = ',ccwood_min
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ci_st = ',ci_st
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dcatch_dlai = ',dcatch_dlai
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' deact_vcmax = ',deact_jmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' deact_vcmax = ',deact_vcmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dfp_dcuo = ',dfp_dcuo
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dgl_dm = ',dgl_dm
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dgl_dt = ',dgl_dt
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dqcrit = ',dqcrit
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ds_jmax = ',ds_jmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ds_vcmax = ',ds_vcmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dust_veg_scj = ',dust_veg_scj
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' dz0v_dh = ',dz0v_dh
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' emis_pft = ',emis_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' eta_sl = ',eta_sl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' f0 = ',f0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fd = ',fd
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_bc = ',fef_bc
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_ch4 = ',fef_ch4
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_co = ',fef_co
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_co2 = ',fef_co2
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_nox = ',fef_nox
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_oc = ',fef_oc
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_so2 = ',fef_so2
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_c2h4 = ',fef_c2h4
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_c2h6 = ',fef_c2h6
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_c3h8 = ',fef_c3h8
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_hcho = ',fef_hcho
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_mecho = ',fef_mecho
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_nh3 = ',fef_nh3
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fef_dms = ',fef_dms
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fire_mort = ',fire_mort
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fl_o3_ct = ',fl_o3_ct
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fsmc_of = ',fsmc_of
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fsmc_p0 = ',fsmc_p0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sug_g0 = ',sug_g0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' g1_stomata = ',g1_stomata
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' g_leaf_0 = ',g_leaf_0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' glmin = ',glmin
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' gpp_st = ',gpp_st
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sug_grec = ',sug_grec
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' gsoil_f = ',gsoil_f
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' hw_sw = ',hw_sw
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ief = ',ief
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' infil_f = ',infil_f
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' jv25_ratio = ',jv25_ratio
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kext = ',kext
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kn = ',kn
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' knl = ',knl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kpar = ',kpar
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' lai_alb_lim = ',lai_alb_lim
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' lma = ',lma
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' mef = ',mef
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' neff = ',neff
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' nl0 = ',nl0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' nmass = ',nmass
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' nr_nl = ',nr_nl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ns_nl = ',ns_nl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' nsw = ',nsw
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' nr = ',nr
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' omega = ',omega
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' omnir = ',omnir
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' orient = ',orient
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' q10_leaf = ',q10_leaf
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' r_grow = ',r_grow
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' rootd_ft = ',rootd_ft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sigl = ',sigl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' tef = ',tef
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' tleaf_of = ',tleaf_of
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' tlow = ',tlow
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' tupp = ',tupp
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' vint = ',vint
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' vsl = ',vsl
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sug_yg = ',sug_yg
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' z0v = ',z0v
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sox_a = ',sox_a
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sox_p50 = ',sox_p50
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sox_rp_min = ',sox_rp_min
CALL jules_print('pftparm',lineBuffer)

! Stomatal optimisation / root water uptake parameters (JBaguley, gs_opt_dev).
! kmax_pft and kcrit are in model units (mol m-2 s-1 Pa-1); the namelist
! kmax_pft_io is in mmol m-2 s-1 MPa-1.
WRITE(lineBuffer,*)' calc_rz_psi = ',calc_rz_psi
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' min_gl_pft = ',min_gl_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' min_rootc_pft = ',min_rootc_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' root_psi_crit = ',root_psi_crit
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' root_radi_pft = ',root_radi_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' rootc_density_pft = ',rootc_density_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' rmass = ',rmass
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' leaf_crit = ',leaf_crit
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' pft_conductance_model = ',pft_conductance_model
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kcrit_fractional_loss = ',kcrit_fractional_loss
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kmax_pft = ',kmax_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' kcrit = ',kcrit
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' P50 = ',P50
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' P88 = ',P88
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' conductance_b_pft = ',conductance_b_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' conductance_c_pft = ',conductance_c_pft
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' gcuticular = ',gcuticular
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_nsl_onset = ',psi_nsl_onset
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_nsl0 = ',psi_nsl0
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' cmax_a = ',cmax_a
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' cmax_b = ',cmax_b
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' cgain_varpi = ',cgain_varpi
WRITE(lineBuffer,*)' sl_cica_well_watered = ',sl_cica_well_watered
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' fsmc_q = ',fsmc_q
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_vcmax_f = ',psi_vcmax_f
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sf_vcmax = ',sf_vcmax
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_vcmax_fmin = ',psi_vcmax_fmin
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' soil_litter_depth = ',soil_litter_depth
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' or_z0soil_fac = ',or_z0soil_fac
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' g1_tuzet = ',g1_tuzet
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' sf_tuzet = ',sf_tuzet
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' psi_f_tuzet = ',psi_f_tuzet
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' cap_leaf = ',cap_leaf
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' cap_stem = ',cap_stem
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' seg_kfac = ',seg_kfac
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' conductance_b_seg = ',conductance_b_seg
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' conductance_c_seg = ',conductance_c_seg
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' pft_xylem_impairment_model = ',pft_xylem_impairment_model
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_tau_rec = ',ximpair_tau_rec
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_psi_refill = ',ximpair_psi_refill
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_wood_alloc = ',ximpair_wood_alloc
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_leaf_sens = ',ximpair_leaf_sens
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_psi_growth = ',ximpair_psi_growth
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_rec_years = ',ximpair_rec_years
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_tau_stem = ',ximpair_tau_stem
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_tau_leaf = ',ximpair_tau_leaf
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_psi_driver = ',ximpair_psi_driver
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_reset_mmdd = ',ximpair_reset_mmdd
CALL jules_print('pftparm',lineBuffer)
WRITE(lineBuffer,*)' ximpair_growth_basis = ',ximpair_growth_basis
CALL jules_print('pftparm',lineBuffer)

CALL jules_print('pftparm',                                                    &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_pftparm


SUBROUTINE check_jules_pftparm(npft,nnpft)

USE jules_soil_biogeochem_mod, ONLY: l_layeredC, soil_bgc_model,               &
                                      soil_model_4pool

USE jules_vegetation_mod, ONLY: can_rad_mod, l_crop, l_trait_phys,             &
                                 l_use_pft_psi, l_bvoc_emis, l_inferno,        &
                                 l_o3_damage, l_trif_fire, photo_acclim_model, &
                                 photo_act_model, photo_act_pft,               &
                                 photo_farquhar, photo_johnson, photo_model,   &
                                 stomata_jacobs, stomata_medlyn,              &
                                 stomata_sox_analytical,                      &
                                 stomata_model, l_spec_veg_z0, l_sugar,        &
                                 l_scale_resp_pm, stomata_desica,              &
                                 som_psi_solver, psi_solver_lut,               &
                                 l_ximpair_seg_memory, l_som_plant_segments,   &
                                 l_ximpair_rec_growth, xylem_impairment_memory

USE jules_radiation_mod, ONLY: l_spec_albedo, l_albedo_obs, l_snow_albedo

USE missing_data_mod, ONLY: rmdi

USE ereport_mod,     ONLY: ereport
USE jules_print_mgr, ONLY: jules_print

IMPLICIT NONE

!Arguments
INTEGER, INTENT(IN) :: npft, nnpft

! Work variables
INTEGER :: ERROR  ! Error indicator

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_JULES_PFTPARM'

!-----------------------------------------------------------------------------
! Check that all required variables were present in the namelist.
! The namelist variables were initialised to rmdi.
! Some configurations don't need all parameters but in some cases these are
! still tested below.
!-----------------------------------------------------------------------------
ERROR = 0
#if !defined(UM_JULES)
! Trigger ignored in the UM for now use ifdef, but may be better regarding LFRic
! and other triggered off option to check the value if any > rmdi and then
! check if they should have a value in check_available_options or based on
! science options.
IF ( ANY( fsmc_mod(:) < 0 ) ) THEN  ! fsmc_mod was initialised to < 0
  ERROR = 1
  CALL jules_print(routinename, "No value for fsmc_mod")
END IF
#endif
IF ( ANY( orient(:) < 0 ) ) THEN  ! orient was initialised to < 0
  ERROR = 1
  CALL jules_print(routinename, "No value for orient")
END IF
IF ( ANY( ABS( kext(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for kext")
END IF
IF ( ANY( ABS( kpar(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for kpar")
END IF
IF ( ANY( ABS( lai_alb_lim(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for lai_alb_lim")
END IF
IF ( ANY( ABS( can_struct_a(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for can_struct_a")
END IF
IF ( ANY( ABS( gsoil_f(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for gsoil_f")
END IF
IF ( ANY( c3(:) < 0 ) ) THEN  ! c3 was initialised to < 0
  ERROR = 1
  CALL jules_print(routinename, "No value for c3")
END IF
IF ( ANY( ABS( alpha(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for alpha")
END IF
IF ( ANY( ABS( fd(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for fd")
END IF
IF ( ANY( ABS( nr_nl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for nr_nl")
END IF
IF ( ANY( ABS( ns_nl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for ns_nl")
END IF
IF ( ANY( ABS( r_grow(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for r_grow")
END IF
! Note that tlow and tupp are always required, for some PFTs at least.
! If using the Farquhar model for C3 plants, we still need tlow and
! tupp for C4 plants.
IF ( ANY( ABS( tlow(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for tlow")
END IF
IF ( ANY( ABS( tupp(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for tupp")
END IF

SELECT CASE ( photo_model )
CASE ( photo_farquhar, photo_johnson )
  !---------------------------------------------------------------------------
  ! First check parameters that are always required with this model.
  !---------------------------------------------------------------------------
  ! Note that these parameter values are not used for C4 plants, but
  ! here we're still checking that they have been provided.
  IF ( ANY( ABS( alpha_elec(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alpha_elec")
  END IF
  IF ( ANY( ABS( deact_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for deact_jmax")
  END IF
  IF ( ANY( ABS( deact_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for deact_vcmax")
  END IF
  !---------------------------------------------------------------------------
  ! Check parameters that depend on any chosen acclimation model.
  !---------------------------------------------------------------------------
  IF ( photo_acclim_model == 0 ) THEN
    ! No acclimation.
    IF ( ANY( ABS( ds_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for ds_jmax")
    END IF
    IF ( ANY( ABS( ds_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for ds_vcmax")
    END IF
    IF ( ANY( ABS( jv25_ratio(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for jv25_ratio")
    END IF
  END IF  !  photo_acclim_model == 0

  IF ( photo_acclim_model == 0 .OR.                                            &
      (photo_acclim_model /= 0 .AND. photo_act_model == photo_act_pft) ) THEN
    IF ( ANY( ABS( act_jmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for act_jmax")
    END IF
    IF ( ANY( ABS( act_vcmax(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for act_vcmax")
    END IF
  END IF  !  photo_acclim_model

END SELECT  !  photo_model

SELECT CASE ( stomata_model )
CASE ( stomata_jacobs )
  IF ( ANY( ABS( dqcrit(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for dqcrit")
  END IF
  IF ( ANY( ABS( f0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for f0")
  END IF
CASE ( stomata_medlyn )
  IF ( ANY( ABS( g1_stomata(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for g1_stomata")
  END IF
CASE ( stomata_sox_analytical )
  IF ( ANY( ABS( sox_p50(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sox_p50")
  END IF
  IF ( ANY( ABS( sox_rp_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sox_rp_min")
  END IF
  IF ( ANY( ABS( sox_a(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sox_a")
  END IF
END SELECT

IF ( .NOT. l_spec_albedo .AND. can_rad_mod == 1 ) THEN
  ! These don't need to be set
ELSE
  IF ( ANY( ABS( alnir(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alnir")
  END IF
  IF ( ANY( ABS( alpar(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alpar")
  END IF
  IF ( ANY( ABS( omega(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omega")
  END IF
  IF ( ANY( ABS( omnir(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omnir")
  END IF
END IF

IF ( l_albedo_obs ) THEN
  IF ( ANY( ABS( alnirl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alnirl")
  END IF
  IF ( ANY( ABS( alniru(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alniru")
  END IF
  IF ( ANY( ABS( alparl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alparl")
  END IF
  IF ( ANY( ABS( alparu(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for alparu")
  END IF
  IF ( ANY( ABS( omegal(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omegal")
  END IF
  IF ( ANY( ABS( omegau(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omegau")
  END IF
  IF ( ANY( ABS( omnirl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omnirl")
  END IF
  IF ( ANY( ABS( omniru(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for omniru")
  END IF
END IF

IF ( .NOT. l_spec_albedo ) THEN
  IF ( ANY( ABS( albsnf_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for albsnf_max")
  END IF
  IF ( l_albedo_obs ) THEN
    IF ( ANY( ABS( albsnf_maxl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for albsnf_maxl")
    END IF
    IF ( ANY( ABS( albsnf_maxu(:) - rmdi ) < EPSILON(1.0) ) ) THEN
      ERROR = 1
      CALL jules_print(routinename, "No value for albsnf_maxu")
    END IF
  END IF
END IF

IF ( .NOT. l_snow_albedo ) THEN
  IF ( ANY( ABS( albsnc_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for albsnc_max")
  END IF
  IF ( ANY( ABS( albsnc_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for albsnc_min")
  END IF
END IF

IF (l_trait_phys) THEN
  IF ( ANY( ABS( lma(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for lma")
  END IF
  IF ( ANY( ABS( nmass(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for nmass")
  END IF
  IF ( ANY( ABS( vsl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for vsl")
  END IF
  IF ( ANY( ABS( vint(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for vint")
  END IF
  IF ( ANY( ABS( nr(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for nr")
  END IF
  IF ( ANY( ABS( nsw(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for nsw")
  END IF
  IF ( ANY( ABS( hw_sw(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for hw_sw")
  END IF
ELSE
  IF ( ANY( ABS( neff(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for neff")
  END IF
  IF ( ANY( ABS( nl0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for nl0")
  END IF
  IF ( ANY( ABS( sigl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sigl")
  END IF
END IF !l_trait_phys

IF ( ANY( ABS( kn(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for kn")
END IF
IF ( ( can_rad_mod == 6 .OR. can_rad_mod == 7 ) .AND.                         &
     ANY( ABS( knl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for knl")
END IF
IF ( ANY( ABS( q10_leaf(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for q10_leaf")
END IF
IF ( ANY( ABS( a_wl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for a_wl")
END IF
IF ( ANY( ABS( a_ws(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for a_ws")
END IF
IF ( ANY( ABS( b_wl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for b_wl")
END IF
IF ( ANY( ABS( eta_sl(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for eta_sl")
END IF
IF ( ANY( ABS( g_leaf_0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for g_leaf_0")
END IF
IF ( ANY( ABS( dgl_dm(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for dgl_dm")
END IF
IF ( ANY( ABS( fsmc_of(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for fsmc_of")
END IF
IF ( ANY( ABS( dgl_dt(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for dgl_dt")
END IF
IF ( ANY( ABS( tleaf_of(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for tleaf_of")
END IF
IF ( ANY( ABS( catch0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for catch0")
END IF
IF ( ANY( ABS( dcatch_dlai(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for dcatch_dlai")
END IF
IF ( ANY( ABS( infil_f(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for infil_f")
END IF
IF ( ANY( ABS( glmin(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for glmin")
END IF
IF ( .NOT. l_spec_veg_z0) THEN
  IF ( ANY( ABS( dz0v_dh(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for dz0v_dh")
  END IF
ELSE
  IF ( ANY( ABS( z0v(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for z0v")
  END IF
END IF
IF ( ANY( ABS( rootd_ft(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for rootd_ft")
END IF

IF ( l_use_pft_psi ) THEN
  IF ( ANY( ABS( psi_close(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for psi_close")
  END IF
  IF ( ANY( ABS( psi_open(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for psi_open")
  END IF
ELSE
  IF ( ANY( ABS( fsmc_p0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fsmc_p0")
  END IF
END IF !l_use_pft_psi

IF ( l_bvoc_emis ) THEN
  IF ( ANY( ABS( ci_st(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ci_st")
  END IF
  IF ( ANY( ABS( gpp_st(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for gpp_st")
  END IF
  IF ( ANY( ABS( ief(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ief")
  END IF
  IF ( ANY( ABS( tef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for tef")
  END IF
  IF ( ANY( ABS( mef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for mef")
  END IF
  IF ( ANY( ABS( aef(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for aef")
  END IF
END IF

IF ( l_inferno ) THEN
  IF ( ANY( ABS( fef_co2(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_co2")
  END IF
  IF ( ANY( ABS( fef_co(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_co")
  END IF
  IF ( ANY( ABS( fef_ch4(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_ch4")
  END IF
  IF ( ANY( ABS( fef_nox(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_nox")
  END IF
  IF ( ANY( ABS( fef_so2(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_so2")
  END IF
  IF ( ANY( ABS( fef_c2h4(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_c2h4")
  END IF
  IF ( ANY( ABS( fef_c2h6(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_c2h6")
  END IF
  IF ( ANY( ABS( fef_c3h8(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_c3h8")
  END IF
  IF ( ANY( ABS( fef_hcho(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_hcho")
  END IF
  IF ( ANY( ABS( fef_mecho(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_mecho")
  END IF
  IF ( ANY( ABS( fef_nh3(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_nh3")
  END IF
  IF ( ANY( ABS( fef_dms(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fef_dms")
  END IF
  IF ( ANY( ABS( ccleaf_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ccleaf_min")
  END IF
  IF ( ANY( ABS( ccleaf_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ccleaf_max")
  END IF
  IF ( ANY( ABS( ccwood_min(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ccwood_min")
  END IF
  IF ( ANY( ABS( ccwood_max(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for ccwood_max")
  END IF
  IF ( ANY( ABS( avg_ba(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for avg_ba")
  END IF
END IF

IF ( l_trif_fire ) THEN
  IF ( ANY( ABS( fire_mort(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fire_mort")
  END IF
END IF

IF ( l_o3_damage ) THEN
  IF ( ANY( ABS( fl_o3_ct(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for fl_o3_ct")
  END IF
  IF ( ANY( ABS( dfp_dcuo(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for dfp_dcuo")
  END IF
END IF

IF ( l_sugar ) THEN
  IF ( ANY( ABS( sug_g0(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sug_g0")
  END IF
  IF ( ANY( ABS( sug_grec(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sug_grec")
  END IF
  IF ( ANY( ABS( sug_yg(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for sug_yg")
  END IF
END IF

IF ( ANY( ABS( emis_pft(:) - rmdi ) < EPSILON(1.0) ) ) THEN
  ERROR = 1
  CALL jules_print(routinename, "No value for emis_pft")
END IF

!-----------------------------------------------------------------------------
! Soil moisture availability from soil to root conductance (fsmc_mod = 2,
! JBaguley).
!-----------------------------------------------------------------------------
IF ( ANY( fsmc_mod(:) == 2 ) ) THEN
  IF ( ANY( ABS( min_rootc_pft(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for min_rootc_pft")
  ELSE IF ( ANY( min_rootc_pft(:) < 0.0 ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "Negative value for min_rootc_pft")
  END IF
  IF ( ANY( ABS( root_psi_crit(:) - rmdi ) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for root_psi_crit")
  ELSE IF ( ANY( root_psi_crit(:) > -EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "Positive or zero value for root_psi_crit")
  END IF
  IF ( ANY( root_radi_pft(:) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for root_radi_pft")
  END IF
  IF ( ANY( rootc_density_pft(:) < EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "No value for rootc_density_pft")
  END IF
  IF ( ANY( rmass(:) < EPSILON(1.0) ) .OR. ANY( rmass(:) > 1.0 ) ) THEN
    ERROR = 1
    CALL jules_print(routinename, "rmass must be larger than zero and at most one")
  END IF
END IF

IF ( ERROR /= 0 ) THEN
  CALL ereport(routinename, ERROR,                                             &
                 ": Variable(s) missing from namelist - see earlier " //       &
                 "message(s)")
END IF

!******************************************************************************
! Do we want this in the UM & LFRic???
!-----------------------------------------------------------------------------
! Check that glmin is >0.
! This ensures that wt_ext in subroutine soil_evap cannot become a NaN (which
! it would if gs=glmin and gsoil=0), or blow up, and might well be required
! elsewhere too.
!-----------------------------------------------------------------------------
ERROR = 0
IF ( ANY(glmin < 1.0e-10) ) THEN
  ERROR = -1
  CALL ereport(routinename, ERROR,                                             &
               "Increasing one or more values of glmin - very small " //       &
               "values can cause model to blow up or NaNs")
  WHERE ( glmin < 1.0e-10 )
    glmin = 1.0e-10
  END WHERE
END IF

IF ( l_crop ) THEN
  IF ( ANY( ABS( a_ws(nnpft+1: npft) - 1.0 ) > EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR, "crop tiles should have a_ws=1.0")
  END IF
END IF

IF ( l_use_pft_psi ) THEN
  IF ( ANY( psi_close(1: npft) > EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR, "psi_close should be negative")
  END IF
  IF ( ANY( psi_open(1: npft) > EPSILON(1.0) ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR, "psi_open should be negative")
  END IF
END IF

!-----------------------------------------------------------------------------
! fsmc_mod=1 should not be allowed with a layered 4-pool C model until this has
! been properly evaluated. (With fsmc_mod=1, subroutine root_frac does not
! return the exponential root profile that users might expect.)
!-----------------------------------------------------------------------------
IF ( l_layeredC .AND. ( soil_bgc_model == soil_model_4pool ) .AND.             &
     ANY( fsmc_mod(:) == 1 ) ) THEN
  ERROR = 1
  CALL ereport(routinename, ERROR,                                             &
               "fsmc_mod=1 is not allowed with l_layeredC and 4-pool C model")
END IF

!-----------------------------------------------------------------------------
! stomata_model = stomata_sox_analytical must be used with fsmc_mod = 1
! Cannot be run with l_scale_resp_pm
! Must be run with can_rad_mod = 1 (implementation for can_rad_mod = 6 ongoing)
!-----------------------------------------------------------------------------
IF ( stomata_model == stomata_sox_analytical ) THEN ! SOX
  IF ( l_scale_resp_pm ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'l_scale_resp_pm=T is incompatible with SOX (stomata_model=3)')
  END IF

  IF ( .NOT. ANY ( fsmc_mod(:) == 1 ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'SOX (stomata_model=3) must be used with fsmc_mod = 1')
  END IF

  IF ( .NOT. ( can_rad_mod == 1 ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'SOX (stomata_model=3) must be used with can_rad_mod = 1')
  END IF

  IF ( .NOT. ( photo_model == 3 ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'SOX (stomata_model=3) uses the SOX derivation of Collatz (photo_model=3)')
  END IF
END IF

!-----------------------------------------------------------------------------
! Xylem impairment is coded for the stomatal optimisation's flat and bounded
! Ci searches with the Taylor or Newton leaf-psi solvers: the lookup table
! (keyed by PFT) and DESICA use the intact PFT curve.
!-----------------------------------------------------------------------------
IF ( ANY( pft_xylem_impairment_model(:) /= 0 ) ) THEN
  IF ( som_psi_solver == psi_solver_lut ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'xylem impairment (pft_xylem_impairment_model /= 0) is not coded for '  // &
    'som_psi_solver = 3 (or l_som_fast)')
  END IF
  IF ( stomata_model == stomata_desica ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'xylem impairment (pft_xylem_impairment_model /= 0) is not coded for '  // &
    'DESICA (stomata_model=8)')
  END IF
  IF ( ANY( pft_xylem_impairment_model(1:npft) == 1 .OR.                      &
            pft_xylem_impairment_model(1:npft) == 4 ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'pft_xylem_impairment_model 1 (kmax) and 4 (kmax refit) are retired; '  // &
    'use 2 or 3 (their code is on branch xylem_impairment at 130a575)')
  END IF
END IF

!-----------------------------------------------------------------------------
! Per-segment memory: the memory model's caps on the stem and leaf segments.
! With the growth clock it needs the slow (running-mean) recovery,
! ximpair_rec_years > 0, which then only sets the averaging window.
!-----------------------------------------------------------------------------
IF ( l_ximpair_seg_memory ) THEN
  IF ( .NOT. l_som_plant_segments .OR.                                        &
       ANY( pft_xylem_impairment_model(1:npft) /= 0 .AND.                     &
            pft_xylem_impairment_model(1:npft) /= xylem_impairment_memory ) )  &
       THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'l_ximpair_seg_memory needs l_som_plant_segments and '                  // &
    'pft_xylem_impairment_model = 0 or 3 (memory)')
  END IF
  IF ( l_ximpair_rec_growth .AND.                                             &
       ANY( pft_xylem_impairment_model(1:npft) == xylem_impairment_memory   &
            .AND. ximpair_rec_years(1:npft) <= 0.0 ) ) THEN
    ERROR = 1
    CALL ereport(routinename, ERROR,                                           &
    'l_ximpair_seg_memory with l_ximpair_rec_growth needs '                 // &
    'ximpair_rec_years > 0 (the growth-clock averaging window)')
  END IF
END IF

END SUBROUTINE check_jules_pftparm

END MODULE pftparm
