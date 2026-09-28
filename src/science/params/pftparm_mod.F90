! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module holds surface parameters for each Plant Functional Type (but
! not parameters that are only used by TRIFFID).


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
                 ! can_rad_mod=6 (decay is a function of LAI).
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
,r_Cmass_frac(:)                                                               &
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
,min_glw_pft(:)                                                                 & ! JBaguley
                 ! Minimum leaf conductance to H2O (m/s)
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
,gcut(:)                                                                       &
                 ! Cuticular (minimum) leaf conductance to water vapour, per
                 ! unit leaf area (mmol H2O m-2 s-1), applied as a floor on
                 ! the canopy conductance when l_som_cuticular_floor.
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
,ximpair_leaf_weight(:)                                                        &
                 ! Weighting factor for the leaf conductance compaired to the
                 ! root conductance when calculating the new xyelem impairment.
,ximpair_new_kmax_weight(:)                                                    &
                 ! Weighting factor for the new impaired kmax compaired to the
                 ! current kmax when updating xyelem impairment.
,ximpair_threshold(:)                                                          &
                 ! Minimum change in kamx before the xylem impairment is
                 ! updated.
                 ! NOTE: This value is not directly input by the user, instead
                 !        it is calculated from kmax and a fractional user
                 !        input in ptftparm_io_mod.F90.
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
,ximpair_leaf_sens(:)
                 ! Sensitivity of the canopy to lasting xylem damage
                 !  (l_ximpair_leaf_loss): the phenological state is capped
                 !  at 1 - ximpair_leaf_sens * (1 - k_cap/kmax). 1 keeps leaf
                 !  area in proportion to the conducting capacity, 0 disables.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PFTPARM'

CONTAINS

SUBROUTINE pftparm_alloc(npft)

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

albsnc_max(:)   = 0.0
albsnc_min(:)   = 0.0
albsnf_max(:)   = 0.0
albsnf_maxl(:)  = 0.0
albsnf_maxu(:)  = 0.0
alnir(:)        = 0.0
alnirl(:)       = 0.0
alniru(:)       = 0.0
alpar(:)        = 0.0
alparl(:)       = 0.0
alparu(:)       = 0.0
kext(:)         = 0.0
kpar(:)         = 0.0
lai_alb_lim(:)  = 0.0
omega(:)        = 0.0
omegal(:)       = 0.0
omegau(:)       = 0.0
omnir(:)        = 0.0
omnirl(:)       = 0.0
omniru(:)       = 0.0
orient(:)       = 0.0

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

act_jmax(:)     = 0.0
act_vcmax(:)    = 0.0
alpha(:)        = 0.0
alpha_elec(:)   = 0.0
c3(:)           = 0.0
can_struct_a(:) = 0.0
deact_jmax(:)   = 0.0
deact_vcmax(:)  = 0.0
dqcrit(:)       = 0.0
ds_jmax(:)      = 0.0
ds_vcmax(:)     = 0.0
f0(:)           = 0.0
fd(:)           = 0.0
g1_stomata(:)   = 0.0
jv25_ratio(:)   = 0.0
kn(:)           = 0.0
knl(:)          = 0.0
neff(:)         = 0.0
nl0(:)          = 0.0
nr_nl(:)        = 0.0
ns_nl(:)        = 0.0
r_grow(:)       = 0.0
tlow(:)         = 0.0
tupp(:)         = 0.0

! Traint physiology parameters
ALLOCATE( hw_sw(npft))
ALLOCATE( lma(npft))
ALLOCATE( nmass(npft))
ALLOCATE( nr(npft))
ALLOCATE( nsw(npft))
ALLOCATE( q10_leaf(npft))
ALLOCATE( r_Cmass_frac(npft))
ALLOCATE( vint(npft))
ALLOCATE( vsl(npft))

hw_sw(:)        = 0.0
lma(:)          = 0.0
nmass(:)        = 0.0
nr(:)           = 0.0
nsw(:)          = 0.0
q10_leaf(:)     = 0.0
r_Cmass_frac(:) = 0.49 !JBaguley
vint(:)         = 0.0
vsl(:)          = 0.0

! Allometric parameters
ALLOCATE( a_wl(npft))
ALLOCATE( a_ws(npft))
ALLOCATE( b_wl(npft))
ALLOCATE( eta_sl(npft))
ALLOCATE( sigl(npft))

a_wl(:)         = 0.0
a_ws(:)         = 0.0
b_wl(:)         = 0.0
eta_sl(:)       = 0.0
sigl(:)         = 0.0

! Phenology parameters
ALLOCATE( dgl_dm(npft))
ALLOCATE( dgl_dt(npft))
ALLOCATE( fsmc_of(npft))
ALLOCATE( g_leaf_0(npft))
ALLOCATE( tleaf_of(npft))

dgl_dm(:)       = 0.0
dgl_dt(:)       = 0.0
fsmc_of(:)      = 0.0
g_leaf_0(:)     = 0.0
tleaf_of(:)     = 0.0

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
ALLOCATE( min_glw_pft(npft)) ! JBaguley
ALLOCATE( psi_close(npft))
ALLOCATE( psi_open(npft))
ALLOCATE( root_psi_crit(npft))  ! JBaguley
ALLOCATE( root_radi_pft(npft))  ! JBaguley
ALLOCATE( rootc_density_pft(npft))  ! JBaguley
ALLOCATE( rootd_ft(npft))
ALLOCATE( z0v(npft))
ALLOCATE( pft_big_leaf_corection_nitrogen_reduction_factor(npft))

calc_rz_psi(:)       = .FALSE. ! JBaguley
catch0(:)            = 0.0
dcatch_dlai(:)       = 0.0
dust_veg_scj(:)      = 0.0
dz0v_dh(:)           = 0.0
emis_pft(:)          = 0.0
fsmc_mod(:)          = 0.0
fsmc_p0(:)           = 0.0
glmin(:)             = 0.0
gsoil_f(:)           = 0.0
infil_f(:)           = 0.0
min_glw_pft(:)        = 0.0 ! JBaguley
min_rootc_pft(:)     = 1.0 ! JBaguley M.Williams etal 2001
psi_close(:)         = 0.0
psi_open(:)          = 0.0
root_psi_crit(:)     =-0.1e6 ! JBaguley
root_radi_pft(:)     = 0.0005 ! JBaguley M.Williams etal 2001
rootc_density_pft(:) = 0.5e3 ! JBaguley M.Williams etal 2001
rootd_ft(:)          = 0.0
z0v(:)               = 0.0
pft_big_leaf_corection_nitrogen_reduction_factor(:) = 1.0

! Ozone damage parameters
ALLOCATE( dfp_dcuo(npft))
ALLOCATE( fl_o3_ct(npft))

dfp_dcuo(:) = 0.0
fl_o3_ct(:) = 0.0

! BVOC emission parameters
ALLOCATE( aef(npft))
ALLOCATE( ci_st(npft))
ALLOCATE( gpp_st(npft))
ALLOCATE( ief(npft))
ALLOCATE( mef(npft))
ALLOCATE( tef(npft))

aef(:)    = 0.0
ci_st(:)  = 0.0
gpp_st(:) = 0.0
ief(:)    = 0.0
mef(:)    = 0.0
tef(:)    = 0.0

! INFERNO combustion parameters
ALLOCATE( avg_ba(npft))
ALLOCATE( ccleaf_min(npft))
ALLOCATE( ccleaf_max(npft))
ALLOCATE( ccwood_min(npft))
ALLOCATE( ccwood_max(npft))
ALLOCATE( fire_mort(npft))

avg_ba(:)     = 0.0
ccleaf_min(:) = 0.0
ccleaf_max(:) = 0.0
ccwood_min(:) = 0.0
ccwood_max(:) = 0.0
fire_mort(:)  = 0.0

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

fef_bc(:)   = 0.0
fef_ch4(:)  = 0.0
fef_co(:)   = 0.0
fef_co2(:)  = 0.0
fef_nox(:)  = 0.0
fef_oc(:)   = 0.0
fef_so2(:)  = 0.0
fef_c2h4(:) = 0.0
fef_c2h6(:) = 0.0
fef_c3h8(:) = 0.0
fef_mecho(:)= 0.0
fef_hcho(:) = 0.0
fef_nh3(:)  = 0.0
fef_dms(:)  = 0.0

! SUGAR parameters
ALLOCATE( sug_grec(npft))
ALLOCATE( sug_g0(npft))
ALLOCATE( sug_yg(npft))

sug_grec(:) = 0.0
sug_g0(:)   = 0.0
sug_yg(:)   = 0.0

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
ALLOCATE( ximpair_leaf_weight(npft))
ALLOCATE( ximpair_new_kmax_weight(npft))
ALLOCATE( ximpair_threshold(npft))
ALLOCATE( ximpair_psi_driver(npft))
ALLOCATE( seg_kfac(npft,3))
ALLOCATE( gcut(npft))
ALLOCATE( conductance_b_seg(npft,3))
ALLOCATE( conductance_c_seg(npft,3))
ALLOCATE( ximpair_reset_mmdd(npft))
ALLOCATE( ximpair_growth_basis(npft))
ALLOCATE( ximpair_tau_rec(npft))
ALLOCATE( ximpair_psi_refill(npft))
ALLOCATE( ximpair_wood_alloc(npft))
ALLOCATE( ximpair_leaf_sens(npft))

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
ximpair_leaf_weight(:) = 1.0
ximpair_new_kmax_weight(:) = 1.0
ximpair_threshold(:) = 0.0
ximpair_psi_driver(:) = 1
seg_kfac(:,:) = 1.0
gcut(:) = 3.0
conductance_b_seg(:,:) = 1.0
conductance_c_seg(:,:) = 1.0
ximpair_reset_mmdd(:) = 0
ximpair_growth_basis(:) = 1
ximpair_tau_rec(:) = 0.0
ximpair_psi_refill(:) = -0.5e6
ximpair_wood_alloc(:) = 0.25
ximpair_leaf_sens(:) = 0.0

! SOX parameters
ALLOCATE( sox_a(npft))
ALLOCATE( sox_p50(npft))
ALLOCATE( sox_rp_min(npft))

sox_a(:)      = 0.0
sox_p50(:)    = 0.0
sox_rp_min(:) = 0.0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE pftparm_alloc

END MODULE pftparm
