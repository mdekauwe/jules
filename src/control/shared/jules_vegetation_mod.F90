! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE jules_vegetation_mod

USE max_dimensions, ONLY: npft_max, nsurft_max
USE missing_data_mod, ONLY: rmdi, imdi

!-----------------------------------------------------------------------------
! Description:
!   Contains vegetation options and a namelist for setting them
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Module constants
!-----------------------------------------------------------------------------
INTEGER, PARAMETER ::                                                          &
  i_veg_vn_1b = 1,                                                             &
      ! Constant indicating fixed vegetation scheme
  i_veg_vn_2b = 2
      ! Constant indicating interactive vegetation scheme

! Parameters identifying alternative models of leaf photosynthesis.
! These should have unique values.
! In all cases photosynthesis of C4 plants uses the model of Collatz et
! al., 1992, Aust.J.Plant Physiol., 19, 519-538.
INTEGER, PARAMETER ::                                                          &
  photo_collatz = 1,                                                           &
    ! C3 plants use the model of Collatz et al., 1991, Agricultural and
    ! Forest Meteorology, 54, 107-136.
  photo_farquhar = 2,                                                          &
    ! C3 plants use the model of Farquhar et al. ,1980, Planta, 149: 78-90.
  photo_sox_collatz = 3,                                                       &
    ! C3 plants use the model of Collatz et al.as derived for use with the
    ! SOX stomata model.
  photo_johnson = 4
    ! C3 plants use the Farquhar et al. (1980) model with the electron
    ! transport of Johnson & Berry (2021), Photosynth. Res., 148: 101-136,
    ! through Cyt b6f. Uses the Farquhar parameters; see jb_photo_mod.

! Parameters identifying alternative models for the thermal response of
! photosynthetic capacity.
! These should have unique, non-zero values. A value of zero is used to
! indicate fixed parameters (effectively that no model is used).
INTEGER, PARAMETER ::                                                          &
  photo_adapt  = 1,                                                            &
    ! Thermal adaptation of photosynthesis only. Plant response to
    ! long-term mean "home" temperature varies geographically.
  photo_acclim = 2,                                                            &
    ! Thermal acclimation of photosynthesis only. Plant response to
    ! short-term mean "growth" temperature varies geographically and
    ! temporally.
  photo_adapt_acclim = 3
    ! Thermal adaptation and acclimation are both used.

INTEGER, PARAMETER ::                                                          &
  n_photo_coef = 3
    ! Number of coefficients in each of the thermal acclimation equations.

! Parameters identifying alternative ways for specifying the activation energies
! of Jmax and Vcmax when acclimation is switched on.  These must have unique
! values.
INTEGER, PARAMETER ::                                                          &
  photo_act_pft = 1,                                                           &
    ! Activation energies vary by PFT but not by land point or time. Values are
    ! taken from act_jmax_io and act_vcmax_io in the PFT parameter namelist.
  photo_act_gb = 2
    ! Activation energies vary by land point and/or time but not by PFT. Values
    ! are subject to acclimation as specified by photo_acclim_model.

! Parameters identifying alternative ways for using the ratio J25:V25 (the
! ratio of the maximum rate of carboxylation of Rubisco to the potential rate
! of electron transport, at 25degC).
! These should have unique values.
INTEGER, PARAMETER ::                                                          &
  jv_scale = 1,                                                                &
    ! J25 is found by scaling V25 by the given ratio J25/V25, that is, all
    ! the variation in the ratio comes from varying J25 (while V25 remains
    ! fixed).
  jv_ntotal = 2
    ! J25 and V25 are calculated assuming that the total amount of nitrogen
    ! allocated to photosynthesis remains constant, thus any change in J25
    ! requires a compensatory change in V25.

! Parameters identifying alternative stomatal conductance models.
! These should have unique values.
INTEGER, PARAMETER ::                                                          &
  stomata_jacobs = 1,                                                          &
    ! Use the original JULES model, including the Jacobs closure
    !   - see Eqn.9 of Best et al. (2011), doi:10.5194/gmd-4-677-2011.
  stomata_medlyn = 2,                                                          &
    ! Use the model of Medlyn et al. (2011) - see Eqn.11,
    !   doi: 10.1111/j.1365-2486.2010.02375.x.
  stomata_sox = 3,                                                             &
    ! Use the semi-analytical version of the SOX model (Eller et al 2020)
    ! doi: 10.1111/nph.16419 - Eqns. 4 & 5
  stomata_profit_max = 4,                                                      &
    ! Stomatal optimisation, profit maximisation (Sperry et al. 2017); see
    ! stom_opt_jls_mod.
  stomata_desica = 5,                                                          &
    ! DESICA: Tuzet et al. (2003) stomatal closure on leaf water potential,
    ! with leaf and stem water potentials (and stem storage) from the plant
    ! hydraulics of Xu et al. (2016), doi: 10.1111/nph.14009 (Notes S1).
    ! See desica_jls_mod.
  stomata_sox_profit = 6
    ! Stomatal optimisation with the SOX profit (Eller et al. 2018); see
    ! stom_opt_jls_mod. Not the semi-analytical SOX (stomata_sox).
!
! stomata_model is the one switch for the stomatal scheme. leaf_flux_mod
! (fsmc leaf path or stomatal optimisation) and som_profit_model (which
! profit) are derived from it in check_jules_vegetation. Setting
! leaf_flux_mod = 2 in the namelist (the old way to select the
! optimisation) still works, with a warning: stomata_model is then set
! from som_profit_model.

! Parameters identifying alternate models for determaning the net carbon
! uptake and stomatal conductance within plants.
! These should have unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  leaf_flux_fsmc = 1,                                                          &
    ! Use the fsmc factor to down regulate net photosynthesis based on
    ! the soil moisture stress factor fsmc.
  leaf_flux_stom_opt = 2
    ! Use the stomatal optimisation model to determin the optimal net
    ! photosynthesis and stomatal conductance.

! Parameters identifying alternate physical properties that can be uniformly
! itterated over when determaning the optimal stomatal conductance.
! These should be unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  som_base_parm_ci = 1,                                                        &
    ! Itterate uniformly over intercellular carbon.
  som_base_parm_psi = 2
    ! Iterate uniformly over leaf water potential.

! Parameters identifying different xylem conductance models.
! These should have unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  CW_conductance = 1,                                                          &
    ! Use the cumulative weibull function to calculate the xylem conductance.
    !   k(psi) = kmax * (1 - exp(-b * psi^c))
  SOX_conductance = 2
    ! Use the SOX conductance model when calculating xylem conductance.
    !   k(psi) = kmax / (1 + (psi / b)^c)

! Solvers for the leaf water potential given the transpiration rate
! (som_psi_solver). These should have unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  psi_solver_taylor = 1,                                                       &
    ! Use a zeroth order Taylor series expansion of the xylem conductence
    ! model to aproximate leaf water potential from transpiration rate.
  psi_solver_newton = 2,                                                       &
    ! Use the Newton Raphson method to aproximate leaf water potential
     ! from transpiration rate.
  psi_solver_lut = 3
    ! Invert a per-PFT lookup table of the supply function (the integral
    ! of the vulnerability curve) to get leaf water potential directly
    ! from transpiration rate (cumulative Weibull or SOX conductance).

! Ci searches of the stomatal optimisation (som_ci_search).
INTEGER, PARAMETER ::                                                          &
  som_ci_flat = 1,                                                             &
    ! Evaluate a flat grid of som_n_sample Ci values and take the best.
  som_ci_bounded = 2
    ! Root find for the upper edge ci_b of the feasible Ci range (gl, E and
    ! the loss of k all rise with Ci, so the feasible range is [ccp, ci_b]),
    ! which gives the CG/HC normalisation of the flat grid in the limit of a
    ! fine grid, then golden-section on [ccp, ci_b] (som_n_ci_golden_iter).
    ! See stom_opt_bounded_search. Fastest with som_psi_solver = 3, which
    ! gives the edge directly. SOX_profit_model always uses the flat grid.

! Parameters identifying different profit models for determaning the optimal
! stomatal conductance.
! These should have unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  profit_max_profit_model = 1,                                                 &
    ! Use the profit max profit model to determan the optimal
    !  stomatal conductance.
    !      Profit(psi) = CarbonGain(psi) - HydraulicCost(psi)
    !
    !      CarbonGain(psi) = Anet(psi) / max(Anet(psi))
    !
    !      HydraulicCost(psi) =   (max(k(psi)) - k(psi))
    !                           / (max(k(psi)) - kcrit))
  SOX_profit_model = 2
    ! Use the SOX profit model to determan the optimal stomatal conductance.
    !      Profit(psi) = CarbonGain(psi) * (1 - HydraulicCost(psi))
    !
    !      CarbonGain(psi) = Anet(psi)
    !
    !      HydraulicCost(psi) = 1 - (k(psi) / kmax)

!-----------------------------------------------------------------------------
! Items set in namelist
!-----------------------------------------------------------------------------
LOGICAL ::                                                                     &
  l_nrun_mid_trif = .FALSE.,                                                   &
      ! Switch for starting NRUN mid-way through a TRIFFID period
  l_trif_init_accum = .FALSE.,                                                 &
      ! Switch so that an NRUN will bit-compare with a CRUN when FALSE
  l_phenol = .FALSE.,                                                          &
      ! Switch for leaf phenology
  l_triffid = .FALSE.,                                                         &
      ! Switch for interactive veg model
  l_trif_eq = .FALSE.,                                                         &
      ! Switch for running TRIFFID in equilibrium mode
  l_veg_compete  = .FALSE.,                                                    &
      ! Switch for competing vegetation
      ! Setting l_triffid = .TRUE. and this as .FALSE. means that the carbon
      ! pools evolve but the PFT distribution does not change
  l_trait_phys  = .FALSE.,                                                     &
      ! .TRUE. for new trait-based PFTs (uses nmass & lma)
      ! .FALSE. for pre-new PFT configuration (nl0, sigl, and neff)
  l_ht_compete  = .FALSE.,                                                     &
      ! Switch for TRIFFID competition
      ! (T for height F for lotka)
      ! Must be true if npft > 5!
  l_landuse = .FALSE.,                                                         &
      ! Switch for landuse change that invokes wood product pools
  l_nitrogen = .FALSE.,                                                        &
      ! Switch for Nitrogen limiting NPP
  l_bvoc_emis = .FALSE.,                                                       &
      ! Switch to enable calculation of BVOC emissions
  l_o3_damage    = .FALSE.,                                                    &
      ! Switch for ozone damage
  l_prescsow = .FALSE.,                                                        &
      ! Only used if crop model is on ( ncpft > 0 )
      !   T => read in the sowing dates for each crop
      !   F => let the model determine sowing date
  l_croprotate = .FALSE.,                                                      &
      ! Only used if crop model is on ( ncpft > 0 )
      ! and l_prescsow is T
  l_recon = .TRUE.,                                                            &
      ! Used to switch on reconfiguration of veg fractions for TRIFFID
  l_trif_crop = .FALSE.,                                                       &
      ! switch to prevent crop and natural PFTs competing
  l_trif_biocrop = .FALSE.,                                                    &
      ! Switch to enable periodic harvesting of bioenergy crop PFTs
  l_inferno = .FALSE.,                                                         &
      ! Switch used to control whether the Interactive fire scheme is used
  l_trif_fire = .FALSE.,                                                       &
      ! Switch used to control whether interactive fire is used
      !   T => if l_inferno is also true, g_burn is calculated in INFERNO
!   and passed to TRIFFID to calculate emissions and vegetation
       !   dynamics
       !   T => if l_inferno is false, interactive fire is calculated via
!   ancillary if provided, and is 0 if not provided
       !   F => g_burn is calculated via ancillary if provided, and is 0 if
!   not provided
   l_use_pft_psi = .FALSE.,                                                    &
       ! Switch used to control what parameters are used in the calculation
       ! of the soil moisture stress factor
       !   T => use psi_close and psi_open
       !   F => use sm_wilt, sm_crit and fsmc_p0
   l_spec_veg_z0 = .FALSE.,                                                    &
       ! Switch used to specify the roughness length for vegetation rather
       ! than calculate it from the vegetation canopy height
   l_limit_canhc = .FALSE.,                                                    &
       ! Switch to limit the value of vegetation canopy heat capacity to
       ! the value specified in the subroutine CANCAP
  l_ag_expand = .FALSE.,                                                       &
       ! T => New crop and biocrop areas are automatically "planted" with
       !      appropriate PFTs.
       ! F => Expansion of new PFTs into new crop and biocrop areas depends on
       !      available NPP.
   l_sugar = .FALSE.,                                                          &
       ! Switch fo the non-structural carbohydrate model

! Switches for bug fixes.
    l_leaf_n_resp_fix = .FALSE.,                                               &
        ! Switch to use correct forms for canopy-average leaf nitrogen.
        ! This affects can_rad_mod = 1, 4 and 5, not 6 (which is correct).
    l_stem_resp_fix = .FALSE.,                                                 &
        ! Switch used to control whether LAI or LAI_BAL is used in stem
        ! resp calculation. Only affects non-crop PFTs with l_trait_phys=F.
    l_vegcan_soilfx = .FALSE.,                                                 &
        ! Switch to modify the canopy model to allow for conduction through
        ! the soil below vegetation.
    l_gleaf_fix = .FALSE.,                                                     &
        ! Used to fix a bug accumulating g_leaf_phen_ac in standalone JULES
        !
        ! This bug occurs in standalone JULES because veg2 is called on TRIFFID
        ! timesteps and veg1 is called on phenol timesteps
        !
        ! In the UM, veg2 (where the accumulation is working) is called on
        ! TRIFFID AND phenol timesteps if l_triffid = TRUE and veg1 is called
        ! on phenol timesteps if l_triffid = FALSE, hence this bug doesn't arise
        !
        ! This means we don't bother adding it to the UM namelist transfer
        ! between PEs
    l_scale_resp_pm = .FALSE.,                                                 &
        ! Switch for scaling whole plant maintenance respiraiton by the soil
        ! moisture stress factor. FALSE = Only scale leaf respiration;
        ! TRUE = scale whole plant respiration.
    l_vegdrag_surft(nsurft_max),                                               &
        ! Switch for using vegetation canopy drag scheme on each tile.
        ! Must be false for non-PFT tiles
    l_vegdrag_pft(npft_max),                                                   &
        ! Switch for using vegetation canopy drag scheme.
        ! This is what appears in the namelist, to ensure that true
        ! values can only given for pfts
    l_rsl_scalar = .FALSE.,                                                    &
        ! Switch for using roughness sublayer correction scheme in scalar
        ! variables. This is only valid when l_vegdrag_surft = .true.
    l_red = .FALSE.
        ! Switch for using the Robust Ecosystem Demography (RED).

DATA l_vegdrag_surft / nsurft_max * .FALSE. /
DATA l_vegdrag_pft / npft_max * .FALSE. /

INTEGER ::                                                                     &
  can_model = 4,                                                               &
      ! Switch for thermal vegetation
  can_rad_mod = 4,                                                             &
      ! Canopy radiation model
  ilayers = imdi,                                                              &
      ! Number of layers for canopy radiation model
  leaf_flux_mod = 1
      ! Leaf flux path (1: fsmc, 2: stomatal optimisation). Derived from
      ! stomata_model; as a namelist input it is deprecated (see above).
      ! JBaguley

! Stomatal optimisation model (som) integers
INTEGER ::                                                                     &
  som_base_parm = 1,                                                           &
      ! Switch used to determin the physical property that is uniformly
      ! itterated over when determaning the optimal stomatal conductance.
      ! JBaguley
  som_n_sample = 100,                                                          &
      ! Number of sample points used by the stomatal optimisation model.
      ! JBaguley
  som_n_ci_golden_iter = 16,                                                  &
      ! Maximum golden-section iterations of the bounded Ci search, on
      ! top of the initial two-point bracket setup.
  som_psi_solver = psi_solver_taylor,                                          &
      ! Leaf water potential from transpiration rate: 1 Taylor series,
      ! 2 Newton-Raphson, 3 lookup table of the supply function.
  som_ci_search = som_ci_flat,                                                 &
      ! Ci search of the stomatal optimisation: 1 flat grid, 2 bounded.
  som_psi_aprox_method = imdi,                                                 &
      ! Deprecated namelist input: use som_psi_solver (same values).
  som_profit_model = 1
      ! Profit used by the stomatal optimisation (1: profit max, 2: SOX).
      ! Derived from stomata_model (5 or 6); only read from the namelist
      ! with the deprecated leaf_flux_mod = 2. JBaguley

LOGICAL ::                                                                     &
  l_som_gain_gross = .FALSE.
      ! When .TRUE., the carbon gain in the profit-max search is normalised
      ! gross assimilation, (An + Rd)/max(An + Rd), instead of net An/max(An).
      ! Rd does not depend on gs, so it cannot be traded against water, but
      ! normalising net A stretches the gain when Rd is a large part of A
      ! (low light) and fails near the compensation point. GPP/NPP are
      ! unaffected (only the optimisation's gain changes).

LOGICAL ::                                                                     &
  l_som_cuticular_floor = .FALSE.
      ! When .TRUE., leaf water loss never falls below a cuticular floor:
      ! after the profit-max search (and for closed stomata, including at
      ! night) the canopy conductance is MAX(gs, gcut_io * LAI). The floor is
      ! an uncontrolled leak, so it is not part of the optimisation and adds
      ! no carbon; the extra water is taken from the soil through the
      ! normal evaporation path (still bounded by the soil-supply cap when
      ! l_som_supply_limit) and psi_leaf is re-solved for the total flux.

LOGICAL ::                                                                     &
  l_som_gravity = .FALSE.
      ! When .TRUE., the stomatal optimisation (stomata_model = 4 or 6) takes
      ! the gravitational drop rho_water g h to the canopy height h off the
      ! root zone water potential, so the plant path starts from
      ! psi_root_zone - rho_water g h (0.01 MPa per m), as DESICA's psi_h.

LOGICAL ::                                                                     &
  l_som_plant_segments = .FALSE.
      ! When .TRUE., the plant hydraulics are three segments in series
      ! (root, stem, leaf; as in GDAY gs_opt) instead of one: the segments
      ! share the whole-plant resistance by seg_frac_*_io (pft_params) and
      ! each has its own vulnerability curve (p50/p88_root/stem/leaf_io,
      ! default the PFT's p50_io/p88_io). Leaf psi comes from solving the
      ! segments downstream, and the hydraulic cost uses the whole-plant
      ! conductance k = -dE/dpsi_leaf (Sperry et al. 2017). With .FALSE.
      ! (default) the hydraulics are the single segment as before.

LOGICAL ::                                                                     &
  l_som_supply_limit = .FALSE.
      ! When .TRUE., the stomatal optimisation can only choose transpiration
      ! the soil can supply this timestep: smc (the available water that
      ! sf_evap already caps esoil at) per timestep. It enters the Ci search
      ! as a cap on gl (with gl_max), so on steps where the soil can't meet
      ! the demand gs and A are re-derived by the optimiser and stay
      ! consistent with the water actually used; other steps are unchanged.
LOGICAL ::                                                                     &
  l_som_root_supply = .FALSE.
      ! When .TRUE., the stomatal optimisation can only choose transpiration
      ! the roots can take up with the root held at root_psi_crit:
      !   E <= sum_layers soil_to_root_k * MAX(psi_soil - root_psi_crit, 0)
      !        / (rho_water g),
      ! so uptake falls as the soil (and rhizosphere conductance) dries and
      ! stops at root_psi_crit, which bounds the soil water potential the
      ! plant can produce. Uses the same cap on gl as l_som_supply_limit
      ! (with which it combines, taking the smaller supply).
LOGICAL ::                                                                     &
  l_leaf_temp = .FALSE.
      ! When .TRUE., the two-leaf canopy (can_rad_mod = 7) solves a leaf
      ! energy balance for the sunlit and the shaded leaf (Penman-Monteith;
      ! Leuning et al., 1995; Wang & Leuning, 1998) and runs photosynthesis
      ! and the stomatal optimisation of each at its own leaf temperature,
      ! instead of both at the surface temperature tstar. The tile energy
      ! balance is unchanged (it still solves one tstar from gc).
REAL(KIND=real_jlslsm) ::                                                      &
  leaf_width = 0.05
      ! Characteristic leaf width (m), for the leaf boundary layer
      ! conductance (l_leaf_temp).
INTEGER ::                                                                     &
  leaf_temp_iter = 3
      ! Number of passes of the two-leaf stomatal optimisation with
      ! l_leaf_temp (each pass is one stom_opt_mod call per leaf, then a
      ! leaf temperature update). iter (3) without l_leaf_temp.
REAL(KIND=real_jlslsm) ::                                                      &
  leaf_shelter = 1.0
      ! With l_leaf_temp, sheltering factor dividing the forced-convection
      ! leaf boundary-layer conductance (CABLE shelrb; 2 in CABLE).
INTEGER ::                                                                     &
  leaf_aero_model = 0
      ! With l_leaf_temp, the leaves' aerodynamic environment:
      ! 0: as the two-leaf model two_leaf_at_WTC (Wang & Leuning, 1998;
      !    Leuning et al., 1995): the leaves exchange directly with the
      !    level-1 air (t_c = tair, q_c = q1), at the level-1 wind speed.
      ! 1: CABLE's canopy aerodynamics (Raupach, 1994; Raupach et al., 1997,
      !    CSIRO SCAM): the canopy air coupled to level 1 by rt1 from the
      !    roughness-sublayer theory (6-9 s m-1 at FR-Pue for 2-3 m s-1 at
      !    12 m, neutral), with the wind at the canopy top, u* / (u*/u_h),
      !    declining as exp(-coexp L/2) into the canopy. u* and rt1 follow
      !    CABLE's Monin-Obukhov stability, iterated (4 times, as CABLE's
      !    niter) on the leaves' sensible and latent heat (the soil's are
      !    not known in the leaf solve). The top-leaf boundary-layer
      !    conductance is CABLE's gbvtop (Pohlhausen 0.7, viscosity of air;
      !    floor 0.05 mol m-2 s-1) and u* uses CABLE's 1 m s-1 minimum wind.
      !    Use with leaf_shelter = 2 for CABLE's shelrb; see also
      !    l_leaf_coexp_lai.
      ! 2: through a canopy air space coupled to level 1 by JULES's ra
      !    (neutral: physiol sets rib = 0), at the level-1 wind speed; for
      !    comparison. At FR-Pue (forcing at 12 m over a 5.5 m canopy) ra is
      !    30-45 s m-1 and the canopy air ~8 K warmer than the air at
      !    midday in summer, against an observed radiometric surface
      !    temperature 2-3 K above the air.
LOGICAL ::                                                                     &
  l_leaf_coexp_lai = .FALSE.
      ! With l_leaf_temp and leaf_aero_model = 1, apply the in-canopy wind
      ! extinction coefficient coexp to normalised height, as it is derived
      ! (Raupach, 1994; CSIRO SCAM eq. 3.14: u(z) = u_h exp(-coexp (1 -
      ! z/h))), i.e. coexp / LAI per unit cumulative leaf area. .FALSE.
      ! applies coexp per unit leaf area, as CABLE's gbhu does, which for
      ! LAI > 1 attenuates the wind LAI times too fast (at LAI 5 the
      ! forced-convection conductance of the canopy is 0.17 of the top
      ! leaf's times the leaf area, against 0.59 with .TRUE.).
LOGICAL ::                                                                     &
  l_leaf_temp_gc_eq = .TRUE.
      ! With l_leaf_temp, return to the surface energy balance the canopy
      ! conductance that gives the transpiration chosen by the stomatal
      ! optimisation (at the leaf temperatures) at the surface temperature,
      ! rather than the sum of the leaf stomatal conductances. Keeps the
      ! water used consistent with the hydraulics and the supply limits.
LOGICAL ::                                                                     &
  l_som_nsl = .FALSE.
      ! When .TRUE., a water-potential driven nonstomatal limitation (NSL)
      ! of photosynthesis in the profit-max optimisation (stomata_model = 4,
      ! som_ci_search = 2): gross photosynthesis is scaled by
      !   f = 1 - (psi_leaf - psi_nsl_onset) / (psi_nsl0 - psi_nsl_onset),
      !   clipped to [0, 1]
      ! (pft_params psi_nsl_onset_io, psi_nsl0_io). With psi_nsl_onset = 0
      ! (default) this is Dewar et al. (2022, New Phytol 233: 639) Eqn 3(b),
      ! A = (1 - psi_leaf / psi_0) A0, which acts at every psi_leaf; an onset
      ! at the turgor loss point is an optional variant (no limitation, and
      ! so no marginal cost, above it). f and psi_leaf are solved together for each Ci, and the
      ! limitation reduces both the optimiser's gain and the actual
      ! photosynthesis.
LOGICAL ::                                                                     &
  l_som_vcmax_psi = .FALSE.
      ! When .TRUE., photosynthetic capacity (Vcmax and Jmax) is down-
      ! regulated as the soil dries, by the root-zone water potential (the
      ! predawn proxy), with the Zhou et al. (2013, Agric. For. Meteorol.
      ! 182-183: 204) form f = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f -
      ! psi_rz))) (pft_params psi_vcmax_f_io, sf_vcmax_io; psi in MPa).
      ! Unlike l_som_nsl it does not act in wet soil and does not depend on
      ! the midday leaf psi. Profit max (stomata_model = 4) only.
LOGICAL ::                                                                     &
  l_som_fast = .FALSE.
      ! Deprecated namelist input: .TRUE. sets som_ci_search = 2 (bounded)
      ! and som_psi_solver = 3 (lookup table).
LOGICAL ::                                                                     &
  l_som_skip_search_wellwatered = .FALSE.
      ! When .TRUE., skip the full Ci search in
      ! profit_max_profit_model for any point where hydraulic cost is
      ! negligible even at the most water-demanding candidate Ci (closest
      ! to ca) - see som_hc_negligible_tol and the fast path at the top of
      ! stom_opt_mod's profit_max_profit_model branch.
      ! Off by default: the fast path places ci at
      ! ca - 0.001*(ca - ccp), much closer to ca than the top point of the
      ! flat grid ((ca - ccp)/som_n_sample below ca), so gl there is ~10x
      ! what the flat search itself would select at som_n_sample=100 -
      ! i.e. turning it on changes results, not just speed.

REAL(KIND=real_jlslsm) ::                                                      &
  som_hc_negligible_tol = 0.01,                                                &
      ! Threshold (as a fraction of the kmax-kcrit range) below which
      ! hydraulic cost at the most water-demanding candidate Ci is treated
      ! as negligible by l_som_skip_search_wellwatered.
  som_leaf_resist_frac = 0.5
      ! DESICA only (stomata_desica): fraction of whole-plant hydraulic
      ! resistance (1/kmax_pft) placed in the leaf segment, the remaining
      ! (1 - som_leaf_resist_frac) in the root-side segment. The two
      ! segments are in series, so each gets conductance kmax/frac and
      ! kmax/(1-frac) respectively, keeping the total root-to-leaf
      ! resistance equal to 1/kmax. Must be strictly between 0 and 1.
      ! (The multilayer stomatal optimisation used it for a shared canopy
      ! node until each layer was given its own full parallel path.)

REAL(KIND=real_jlslsm) ::                                                      &
  som_gl_max = 0.02
      ! Maximum leaf stomatal conductance for H2O (m/s, per unit leaf area)
      ! allowed by the stomatal optimisation; Ci samples whose gl exceeds it
      ! are treated as infeasible, like those failing the hydraulic tests.
      ! When transpiration is (nearly) free - VPD ~ 0 - profit(ci) keeps
      ! rising towards ci = ca, where gl = An/(ca - ci) diverges, so without
      ! a cap gl was set only by where the Ci search happened to stop
      ! (som_n_sample, the golden bracket, the 1e-2 Pa floor): canopy gc of
      ! 0.1-1.5 m/s at FR-Pue in low-VPD daylight. 0.02 m/s is ~0.8
      ! mol m-2 s-1, generous for most C3 leaves, so it only binds in that
      ! degenerate regime. Set <= 0 to disable. Every canopy scheme applies
      ! it per unit leaf area, so the canopy cap is som_gl_max * LAI (big
      ! leaf; two-leaf per class leaf area; multilayer per leaf).

REAL(KIND=real_jlslsm) ::                                                      &
  light_curvature_fvcb = 0.90
      ! Curvature (theta) of the light response of electron transport in the
      ! Farquhar model (Eq.4 of Medlyn et al. 2002, who used 0.9). Most
      ! models use 0.7 (e.g. von Caemmerer 2009, CLM5, FATES, MAESPA). With
      ! photo_johnson it only enters the conversion of Jmax to Vqmax (see
      ! jb_photo_mod), and should match the theta Jmax was derived with.

INTEGER ::                                                                     &
  ignition_method = 1,                                                         &
      ! Switch for the calculation method of INFERNO fire ignitions
      ! IGNITION_METHOD=1:Constant (1.67 per km2 per s)
      ! IGNITION_METHOD=2:Constant (Human - 1.5 per km2 per s)
      !                   Varying  (Lightning - see Pechony and Shindell,2009)
      ! IGNITION_METHOD=3:Vary Human and Lightning (Pechony and Shindell,2009)
  fsmc_shape = 0,                                                              &
      ! shape of the soil moisture stress function fsmc
      ! 0: piece-wise linear in vol. soil moisture.
      ! 1: piece-wise linear in soil potential.
  photo_acclim_model = imdi,                                                   &
      ! Chosen model for thermal response of photosynthetic capacity.
  photo_act_model = imdi,                                                      &
      ! Chosen model for activation energies of Jmax and Vcmax.
  photo_jv_model = imdi,                                                       &
      ! Chosen model for the variation of J25:V25 (the ratio of the maximum
      ! rate of carboxylation of Rubisco to the potential rate of electron
      ! transport, at 25degC).
  photo_model = imdi,                                                          &
      ! Chosen model of leaf photosynthesis.
  stomata_model = imdi
      ! Chosen model of stomatal conductance.

INTEGER, PARAMETER :: ignition_constant = 1
INTEGER, PARAMETER :: ignition_vary_natural = 2
INTEGER, PARAMETER :: ignition_vary_natural_human = 3

INTEGER ::                                                                     &
  phenol_period = imdi,                                                        &
      ! Update frequency for leaf phenology (days)
  triffid_period = imdi
      ! Update frequency for TRIFFID (days)

INTEGER :: errcode   ! error code to pass to ereport.

INTEGER :: n ! Loop counter

REAL(KIND=real_jlslsm) ::                                                      &
  frac_min  = rmdi,                                                            &
      ! Minimum areal fraction for PFTs.
  frac_seed = rmdi,                                                            &
      ! "Seed" fraction for PFTs.
  pow = rmdi,                                                                  &
      ! Power in sigmoidal function.
  cd_leaf = rmdi,                                                              &
      ! Leaf level drag coefficient
  c1_usuh = rmdi,                                                              &
      ! u*/U(h) at the top of dense canopy
  c2_usuh = rmdi,                                                              &
      ! u*/U(h) above bare soil
  c3_usuh = rmdi,                                                              &
      ! Used in the exponent of equation weighting dense and sparse
      ! vegetation to get u*/U(h) in neutral condition
  stanton_leaf = rmdi,                                                         &
      ! Leaf-level Stanton number
  !---------------------------------------------------------------------------
  ! Parameters used when the variation of J25:V25 (the ratio of the maximum
  ! rate of carboxylation of Rubisco to the potential rate of electron
  ! transport, at 25degC) is modelled assuming a constant allocation of
  ! nitrogen to photosynthetic components.
  ! Reference: Mercado et al., 2018, New Phytologist,
  ! https://doi.org/10.1111/nph.15100.
  !---------------------------------------------------------------------------
  n_alloc_jmax = rmdi,                                                         &
      ! Constant relating nitrogen allocation to Jmax
      ! (mol CO2 m-2 s-1 [kg m-2]-1).
      ! This is 5.3 in Eq.5 of Mercado et al. (2018).
  n_alloc_vcmax = rmdi,                                                        &
      ! Constant relating nitrogen allocation to Vcmax
      ! (mol CO2 m-2 s-1 [kg m-2]-1).
      ! This is 3.8 in Eq.5 of Mercado et al. (2018).
  !---------------------------------------------------------------------------
  ! Parameters used with temperature adaptation or acclimation of
  ! photosynthesis. Variables suffixed _coef contain coefficients for the
  ! adaptation/acclimation equations of the form,
  !
  ! Q = Q_coef(1) + Q_coef(2)*t_home + Q_coef(3)*t_growth
  !
  ! Jmax is the potential rate of electron transport.
  ! Vcmax is the maximum rate of carboxylation of Rubisco.
  !---------------------------------------------------------------------------
  dsj_coef(n_photo_coef) = rmdi,                                               &
      ! Rate of change with leaf temperature of the entropy factor for Jmax
      ! (J mol-1 K-1).
  dsv_coef(n_photo_coef) = rmdi,                                               &
      ! Rate of change with leaf temperature of the entropy factor for Vcmax
      ! (J mol-1 K-1).
  jv25_coef(n_photo_coef) = rmdi,                                              &
      ! The ratio Jmax:Vcmax at 25 degC (mol electrons mol-1 CO2).
  act_j_coef(n_photo_coef) = rmdi,                                             &
      ! Activation energy of Jmax (J mol-1).
  act_v_coef(n_photo_coef) = rmdi,                                             &
      ! Activation energy of Vcmax (J mol-1).
  n_day_photo_acclim = rmdi
      ! Time constant for exponential moving average of temperature used with
      ! thermal acclimation of photosynthesis (days).
      ! Given a step (down) function as input, the smoothed output has fallen
      ! to 1/e (37%) of the initial value after this number of days.

!-----------------------------------------------------------------------------
! Single namelist definition for UM and standalone
!-----------------------------------------------------------------------------
NAMELIST  / jules_vegetation/                                                  &
! UM only
    l_nrun_mid_trif, l_trif_init_accum,                                        &
! Shared
    l_phenol, l_triffid, l_trif_eq, l_veg_compete,                             &
    phenol_period, triffid_period, l_trait_phys, l_ht_compete,                 &
    l_bvoc_emis, l_o3_damage, can_model, can_rad_mod, ilayers, leaf_flux_mod,  &
    som_base_parm, som_n_sample, som_n_ci_golden_iter,                        &
    l_som_skip_search_wellwatered, som_hc_negligible_tol,                     &
    l_som_fast,                                                               &
    l_som_supply_limit, l_som_root_supply, l_som_nsl, l_som_plant_segments,   &
    l_som_vcmax_psi,                                                           &
    l_leaf_temp, leaf_width, leaf_temp_iter, l_leaf_temp_gc_eq,                &
    leaf_shelter, leaf_aero_model, l_leaf_coexp_lai,                           &
    l_som_gain_gross,                                                          &
    l_som_cuticular_floor, l_som_gravity,                                     &
    som_leaf_resist_frac, som_gl_max, light_curvature_fvcb,                   &
    som_psi_aprox_method, som_profit_model, som_psi_solver, som_ci_search,    &
    frac_min, frac_seed, pow, l_landuse, l_leaf_n_resp_fix, l_stem_resp_fix,   &
    l_nitrogen, l_vegcan_soilfx, l_trif_crop, l_trif_fire,                     &
    l_inferno, ignition_method, l_vegdrag_pft, l_rsl_scalar,                   &
    cd_leaf, c1_usuh, c2_usuh, c3_usuh, dsj_coef, dsv_coef, jv25_coef,         &
    act_j_coef, act_v_coef,                                                    &
    n_alloc_jmax, n_alloc_vcmax, n_day_photo_acclim,                           &
    stanton_leaf, photo_acclim_model, photo_act_model, photo_jv_model,         &
    photo_model, stomata_model, l_spec_veg_z0, l_limit_canhc,                  &
! Not used in the UM yet
    l_prescsow, l_recon,l_gleaf_fix,                                           &
    l_scale_resp_pm, l_use_pft_psi, fsmc_shape, l_croprotate,                  &
    l_trif_biocrop, l_ag_expand, l_sugar,                                      &
    ! RED variables
    l_red

!-----------------------------------------------------------------------------
! Items derived from namelist inputs
!-----------------------------------------------------------------------------
INTEGER :: i_veg_vn = imdi
    ! Switch to determine version of vegetation scheme
    ! Must be one of i_veg_vn_1b or i_veg_vn_2b

LOGICAL :: l_acclim = .FALSE.
    ! Indicates that both Farquhar photosynthesis and adaptation/acclimation
    ! are selected.  Used to make flow control in the UM simpler and safer.

LOGICAL :: l_crop = .FALSE.
    ! Indicates if the crop model is on
    ! This is derived from the value of ncpft

LOGICAL :: l_fapar_diag = .FALSE.
    ! Only true if fapar or apar is in one of the output profiles

LOGICAL :: l_fao_ref_evapotranspiration = .FALSE.
    ! Only true if fao_et0 is in one of the output profiles

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='JULES_VEGETATION_MOD'

CONTAINS


#if !defined(RIVERS_ONLY)
SUBROUTINE check_jules_vegetation()

USE ereport_mod, ONLY: ereport

USE jules_surface_types_mod, ONLY: npft, ncpft, nnpft

USE jules_surface_mod, ONLY: l_aggregate

USE jules_print_mgr, ONLY: jules_message, jules_print

!-----------------------------------------------------------------------------
! Description:
!   Checks JULES_VEGETATION namelist for consistency and calculates some
!   derived values.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------

IMPLICIT NONE


!-----------------------------------------------------------------------------
! The stomatal scheme: stomata_model is the switch; derive the internal
! leaf_flux_mod and som_profit_model from it. Accept the old namelist form
! (leaf_flux_mod = 2 with som_profit_model) with a warning.
!-----------------------------------------------------------------------------
SELECT CASE ( stomata_model )
CASE ( stomata_profit_max, stomata_sox_profit )
  leaf_flux_mod = leaf_flux_stom_opt
  IF ( stomata_model == stomata_profit_max ) THEN
    som_profit_model = profit_max_profit_model
  ELSE
    som_profit_model = SOX_profit_model
  END IF
CASE ( stomata_jacobs, stomata_medlyn, stomata_sox, stomata_desica )
  IF ( leaf_flux_mod == leaf_flux_stom_opt ) THEN
    IF ( stomata_model == stomata_sox .OR. stomata_model == stomata_desica ) THEN
      errcode = 101
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "leaf_flux_mod = 2 cannot be combined with " //             &
                   "stomata_model = 3 or 5; use stomata_model = 4 " //         &
                   "(profit max) or 6 (SOX profit)")
    END IF
    IF ( som_profit_model == SOX_profit_model ) THEN
      stomata_model = stomata_sox_profit
    ELSE
      stomata_model = stomata_profit_max
    END IF
    errcode = -101   ! warning
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "leaf_flux_mod is deprecated: use stomata_model = 4 " //      &
                 "(profit max) or 6 (SOX profit). Taking stomata_model " //    &
                 "from som_profit_model.")
  ELSE
    leaf_flux_mod = leaf_flux_fsmc
  END IF
CASE DEFAULT
  ! Reported below.
END SELECT

! Phenology or TRIFFID cannot be used with the aggregate surface scheme
IF ( l_aggregate .AND. (l_phenol .OR. l_triffid) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Phenology or TRIFFID cannot be used with the ' //              &
               'aggregated surface scheme (i.e. l_aggregate = true)')
END IF

! Check that phenol_period is specified if phenology is on
IF ( l_phenol .AND. phenol_period < 0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Phenology is on but phenol_period is not given')
END IF
! Same for triffid_period if TRIFFID is on
IF ( l_triffid .AND. triffid_period < 0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'TRIFFID is on but triffid_period is not given')
END IF

! Check options that depend on TRIFFID if it is not enabled
IF ( .NOT. l_triffid .AND. ANY( [ l_veg_compete, l_trif_eq, l_landuse,         &
   l_ht_compete, l_nitrogen, l_trif_crop, l_trif_fire, l_trif_biocrop,         &
   l_ag_expand ] ) ) THEN
  errcode = 101
  WRITE(jules_message,'(A,8(1x,L1))')                                          &
     'These should be false when l_triffid = F: l_veg_compete, ' //            &
     'l_trif_eq, l_landuse, l_ht_compete, l_nitrogen, l_trif_crop, ' //        &
     'l_trif_fire, l_trif_biocrop, l_ag_expand = ', l_veg_compete, l_trif_eq,  &
     l_landuse, l_ht_compete, l_nitrogen, l_trif_crop, l_trif_fire,            &
     l_trif_biocrop, l_ag_expand
  CALL ereport( "check_jules_vegetation", errcode, jules_message )
END IF

! Always make sure that a veg version is selected
! If TRIFFID is on, select interactive veg, otherwise select fixed veg
i_veg_vn = i_veg_vn_1b
IF ( l_triffid ) THEN
  i_veg_vn = i_veg_vn_2b
END IF

! Check can_model and can_rad_mod are suitable
IF ( can_model < 1 .OR. can_model > 4 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'can_model should be in range 1 to 4')
END IF

IF ( can_model == 4 .AND. l_aggregate ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'can_model=4 cannot be used with the aggregated ' //            &
               'surface scheme')
END IF

SELECT CASE ( can_rad_mod )
CASE ( 1, 4, 5, 6 )
  ! These are valid, so nothing to do.
CASE ( 7 )
  ! Two-leaf (sunlit/shaded big-leaf) canopy: stomatal optimisation only.
  IF ( leaf_flux_mod /= leaf_flux_stom_opt ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'can_rad_mod=7 requires stomata_model=4 or 6')
  END IF
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'can_rad_mod should be 1, 4, 5, 6 or 7')
END SELECT

! Check that the leaf_flux_mod is suitable. JBaguley
SELECT CASE( leaf_flux_mod)
CASE ( 1, 2)
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'leaf_flux_mod should be 1 or 2')
END SELECT

! Check that the som_base_parm is suitable. JBaguley
SELECT CASE( som_base_parm)
CASE ( 1, 2)
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_bas_parm should be 1 or 2')
END SELECT

! Check that the som_n_sample is suficently large. JBaguley
IF (som_n_sample < 10) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_n_sample should be grater than or equal to 10')
END IF

! Check that the bounded-search golden-section iteration count is sufficiently
! large (at least 1 iteration to refine past the initial two-point bracket).
IF (som_n_ci_golden_iter < 1) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_n_ci_golden_iter should be greater than or equal to 1')
END IF

! Check that som_hc_negligible_tol is a valid fraction.
IF (som_hc_negligible_tol < 0.0 .OR. som_hc_negligible_tol > 1.0) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_hc_negligible_tol should be between 0 and 1')
END IF

! Check that som_leaf_resist_frac is strictly between 0 and 1 (both
! segments need a finite conductance).
IF (som_leaf_resist_frac <= 0.0 .OR. som_leaf_resist_frac >= 1.0) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_leaf_resist_frac should be strictly between 0 and 1')
END IF

! Deprecated inputs: l_som_fast and som_psi_aprox_method.
IF ( som_psi_aprox_method /= imdi ) THEN
  som_psi_solver = som_psi_aprox_method
  errcode = -101   ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               "som_psi_aprox_method is deprecated: use som_psi_solver")
END IF
IF ( l_som_fast ) THEN
  som_ci_search  = som_ci_bounded
  som_psi_solver = psi_solver_lut
  errcode = -101   ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               "l_som_fast is deprecated: use som_ci_search = 2 and " //       &
               "som_psi_solver = 3 (set here)")
END IF

! Check that the som_psi_solver is suitable. JBaguley
SELECT CASE( som_psi_solver )
CASE ( psi_solver_taylor, psi_solver_newton, psi_solver_lut )
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
     'som_psi_solver should be Taylor series (1), Newton-Raphson (2) ' //      &
     'or lookup table (3)')
END SELECT

! Check that the som_ci_search is suitable.
SELECT CASE( som_ci_search )
CASE ( som_ci_flat, som_ci_bounded )
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
     'som_ci_search should be flat (1) or bounded (2)')
END SELECT

! Check that the som_profit_model is suitable. JBaguley
SELECT CASE( som_profit_model)
CASE ( 1, 2)
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_profit_model should be profit max (1) or SOX (2)')
END SELECT

! Check that the photosynthesis option is reasonable.
SELECT CASE ( photo_model )
CASE ( photo_collatz, photo_farquhar, photo_sox_collatz, photo_johnson )
  ! These are valid, nothing more to do.
CASE DEFAULT
  errcode = 101  !  a fatal error
  CALL ereport("check_jules_vegetation", errcode,                              &
               "Invalid value given for photo_model.")
END SELECT

!-----------------------------------------------------------------------------
! Check options for the Farquhar model. Johnson-Berry uses the same
! parameters and options.
!-----------------------------------------------------------------------------
IF (  photo_model == photo_farquhar .OR. photo_model == photo_johnson ) THEN

  IF ( light_curvature_fvcb <= 0.0 .OR. light_curvature_fvcb > 1.0 ) THEN
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "light_curvature_fvcb must be > 0 and <= 1")
  END IF

  ! The Farquhar model of photosynthesis has only been coded for certain
  ! values of can_rad_mod.
  SELECT CASE ( can_rad_mod )
  CASE ( 1, 5, 6, 7 )
    ! These are supported, nothing more to do.
  CASE DEFAULT
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "Farquhar model can only be used with can_rad_mod =1, 5, " // &
                 "6 or 7.")
  END SELECT

  ! Check that the option for acclimation of photosynthesis is reasonable.
  SELECT CASE ( photo_acclim_model )
  CASE ( 0 )
    ! This is valid, nothing more to do.
  CASE ( photo_adapt, photo_acclim, photo_adapt_acclim )
    ! Set a switch to indicate to the UM that extra input fields are needed.
    l_acclim = .TRUE.
  CASE DEFAULT
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "Invalid value given for photo_acclim_model.")
  END SELECT

  !---------------------------------------------------------------------------
  ! Check that the option for acclimation of activation energies is reasonable.
  !---------------------------------------------------------------------------
  SELECT CASE ( photo_act_model )
  CASE ( photo_act_pft )
    ! This is valid, nothing more to do.
  CASE ( photo_act_gb )
    IF ( photo_acclim_model == 0 ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "This value of photo_act_model is not permitted with " //   &
                   "photo_acclim_model == 0.")
    END IF
  CASE DEFAULT
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "Invalid value given for photo_act_model.")
  END SELECT

  !---------------------------------------------------------------------------
  ! Check that the option for variation of J25:V25 is reasonable.
  !---------------------------------------------------------------------------
  ! First we check that we recognise the given value, then we check it is
  ! allowed with the given configuration.
  SELECT CASE ( photo_jv_model )
  CASE ( jv_scale, jv_ntotal )

    ! These are valid. Check the value is allowed with this configuration.
    SELECT CASE ( photo_acclim_model )
    CASE ( 0 )
      ! Without adaptation or acclimation, all variation must come from J25.
      IF ( photo_jv_model /= jv_scale ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "This value of photo_jv_model cannot be used with " //    &
                     "this configuration.")
      END IF
    CASE DEFAULT
      ! With adaptation or acclimation, all values of photo_jv_model are
      ! allowed. ! Nothing more to do.
    END SELECT  !  photo_acclim_model

    !-------------------------------------------------------------------------
    ! Check that any further parameter values related to photo_jv_model
    ! have been provided.
    !-------------------------------------------------------------------------
    IF ( photo_jv_model == jv_ntotal ) THEN
      IF ( ABS( n_alloc_jmax - rmdi ) < EPSILON(rmdi) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "n_alloc_jmax needs to be specified.")
      END IF
      IF ( ABS( n_alloc_vcmax - rmdi ) < EPSILON(rmdi) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "n_alloc_vcmax needs to be specified.")
      END IF
    END IF  !  photo_jv_model

  CASE DEFAULT
    errcode = 101  !  a fatal error
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "Invalid value given for photo_jv_model.")

  END SELECT  !  photo_jv_model

  !---------------------------------------------------------------------------
  ! Check parameters that are used for both adaptation and acclimation.
  !---------------------------------------------------------------------------
  SELECT CASE ( photo_acclim_model )
  CASE ( photo_adapt, photo_acclim, photo_adapt_acclim )
    IF ( ANY(ABS( dsj_coef - rmdi ) < EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "dsj_coef needs to be specified.")
    END IF
    IF ( ANY(ABS( dsv_coef - rmdi ) < EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "dsv_coef needs to be specified.")
    END IF
    IF ( ANY(ABS( jv25_coef - rmdi ) < EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "jv25_coef needs to be specified.")
    END IF

    ! Check that the base coefficients are all strictly positive.
    IF ( ANY([dsj_coef(1), dsv_coef(1), jv25_coef(1)] < EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "Acclimation base coefficients must be positive values.")
    END IF

    ! These params are only needed for the gridbox activation energies model.
    IF ( photo_act_model == photo_act_gb ) THEN
      IF ( ANY(ABS( act_j_coef - rmdi ) < EPSILON(rmdi)) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "act_j_coef needs to be specified.")
      END IF
      IF ( ANY(ABS( act_v_coef - rmdi ) < EPSILON(rmdi)) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "act_v_coef needs to be specified.")
      END IF
      ! Check that the base coefficients are all strictly positive.
      IF ( ANY([act_j_coef(1), act_v_coef(1)] < EPSILON(rmdi)) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "Acclimation base coefficients must be positive values.")
      END IF
    END IF

  END SELECT  !  photo_acclim_model

  !---------------------------------------------------------------------------
  ! When running with adaptation-only or acclimation-only, check that the
  ! coefficients NOT in use are set to zero.
  ! ---------------------------------------------------------------------------
  SELECT CASE ( photo_acclim_model )
  CASE ( photo_adapt )

    IF ( ANY(ABS([dsj_coef(3), dsv_coef(3), jv25_coef(3)])                     &
             > EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "Acclimation coefficients must all be zero.")
    END IF

    IF ( photo_act_model == photo_act_gb ) THEN
      IF ( ANY(ABS([act_j_coef(3), act_v_coef(3)]) > EPSILON(rmdi)) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "Acclimation coefficients must all be zero.")
      END IF
    END IF

  CASE ( photo_acclim )

    IF ( ANY(ABS([dsj_coef(2), dsv_coef(2), jv25_coef(2)])                     &
             > EPSILON(rmdi)) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "Adaptation coefficients must all be zero.")
    END IF

    IF ( photo_act_model == photo_act_gb ) THEN
      IF ( ANY(ABS([act_j_coef(2), act_v_coef(2)]) > EPSILON(rmdi)) ) THEN
        errcode = 101  !  a fatal error
        CALL ereport("check_jules_vegetation", errcode,                        &
                     "Adaptation coefficients must all be zero.")
      END IF
    END IF

  CASE ( photo_adapt_acclim )
    ! No checks necessary.

  END SELECT  !  photo_acclim_model

  !---------------------------------------------------------------------------
  ! Check parameters that are used for acclimation (but not for adaptation).
  !---------------------------------------------------------------------------
  IF ( photo_acclim_model == photo_acclim .OR.                                 &
       photo_acclim_model == photo_adapt_acclim ) THEN

    IF ( ABS( n_day_photo_acclim - rmdi ) < EPSILON(rmdi) ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "n_day_photo_acclim needs to be specified.")
    END IF
    ! Check that the acclimation timescale is long enough to avoid problems
    ! when it appears in the denominator.  This timescale is typically tens of
    ! days, so set 1 day as a hard lower bound.
    IF ( n_day_photo_acclim < 1.0 ) THEN
      errcode = 101  !  a fatal error
      CALL ereport("check_jules_vegetation", errcode,                          &
                   "n_day_photo_acclim must not be less than 1.0 days.")
    END IF
  END IF  !  photo_acclim_model == photo_acclim

END IF  !  photo_model == photo_farquhar or photo_johnson

! Check that the stomatal conductance model is reasonable.
SELECT CASE ( stomata_model )
CASE ( stomata_jacobs, stomata_medlyn, stomata_sox, stomata_desica,           &
       stomata_profit_max, stomata_sox_profit )
  ! These are valid, so nothing to do.
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               "Invalid value for stomata_model" )
END SELECT

! DESICA sets ci from the Tuzet closure in leaf_limits, i.e. it uses the
! fsmc (leaf_flux_mod = 1) leaf path. Stress acts through psi_leaf only
! (fsmc is not applied).
IF ( stomata_model == stomata_desica ) THEN
  ! The within-step gs-psi_leaf solve (bisection on fw) is coded for the
  ! big leaf only.
  IF ( can_rad_mod /= 1 ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "stomata_model = 5 (DESICA) requires can_rad_mod = 1")
  END IF
END IF

IF ( l_triffid .AND. ( .NOT. l_phenol ) ) THEN
  errcode = -105 ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               "When triffid is on, l_phenol=T is recommended. You " //        &
               "have set l_phenol=F. The LAI will be set to the " //           &
               "balanced LAI.")
END IF

! Check triffid-crop options are sensible
IF ( l_trif_crop .AND. l_trif_eq ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'trif_crop and trif_eq are incompatible')
END IF

IF ( l_trif_biocrop .AND. ( .NOT. l_trif_crop) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'trif_biocrop requires trif_crop')
END IF

IF ( l_ag_expand .AND. ( .NOT. (l_trif_crop .AND. l_landuse) ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'ag_expand requires l_trif_crop and l_landuse')
END IF


IF (l_veg_compete .AND. ( .NOT. l_ht_compete ) .AND. ( nnpft /= 5 )) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_ht_compete=F requires 5 natural PFTs: ' //                   &
               'BT, NT, C3, C4, SH')
END IF

! Check crop options are sensible
l_crop = ncpft > 0

IF ( l_crop .AND. l_triffid ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Crop model and triffid are incompatible')
END IF

IF ( l_aggregate .AND. l_crop ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Crop model cannot be used with the aggregated surface ' //     &
               'scheme (i.e. l_aggregate = true)')
END IF

IF ( .NOT. l_crop ) THEN
  ! Note that in the UM, since ncpft_max = 0, l_crop is always .FALSE. and hence
  ! l_prescsow is .FALSE.
  l_prescsow = .FALSE.
  l_croprotate = .FALSE.
END IF

IF ( l_croprotate .AND. .NOT. l_prescsow) THEN
  ! Check that l_prescsow is T when l_croprotate is T
  ! Highlight that prescribed frac should also be provided
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_croprotate=T requires l_prescsow=T, ' //                     &
               'prescribed fractions and crop model to be active')
END IF

! Check a suitable ignition_method was given
IF ( ignition_method /= ignition_constant .AND.                                &
     ignition_method /= ignition_vary_natural .AND.                            &
     ignition_method /= ignition_vary_natural_human ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'ignition_method must be 1, 2 or 3')
END IF

IF ( fsmc_shape == 1 .AND. .NOT. l_use_pft_psi ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'fsmc_shape=1 requires l_use_pft_psi=T ')
END IF

IF ( l_som_plant_segments .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.     &
                                  ( som_psi_solver /= psi_solver_newton .AND.   &
                                    som_psi_solver /= psi_solver_lut ) .OR.     &
                                  ( can_rad_mod /= 1 .AND. can_rad_mod /= 7 ) ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_plant_segments requires stomata_model=4 or 6, ' //            &
               'som_psi_solver=2 or 3 and can_rad_mod=1 or 7')
END IF

IF ( l_som_cuticular_floor .AND.                                               &
     ( ( leaf_flux_mod /= leaf_flux_stom_opt .AND.                            &
         stomata_model /= stomata_desica ) .OR.                               &
       ( can_rad_mod /= 1 .AND. can_rad_mod /= 7 ) ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_cuticular_floor requires stomata_model=4, 5 or 6, ' //   &
               'and can_rad_mod=1 or 7')
END IF

IF ( l_som_gravity .AND. leaf_flux_mod /= leaf_flux_stom_opt ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_gravity requires stomata_model=4 or 6 (DESICA ' //       &
               'always includes gravity)')
END IF

IF ( l_som_supply_limit .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.       &
                                .NOT. l_use_pft_psi ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_supply_limit requires stomata_model=4 or 6 and ' //           &
               'l_use_pft_psi=T')
END IF

IF ( l_som_root_supply .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.        &
                               .NOT. l_use_pft_psi ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_root_supply requires stomata_model=4 or 6 and ' //       &
               'l_use_pft_psi=T')
END IF

IF ( l_leaf_temp .AND. ( can_rad_mod /= 7 .OR.                               &
                         leaf_flux_mod /= leaf_flux_stom_opt ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp is only coded for can_rad_mod=7 with ' //          &
               'stomata_model=4 or 6')
END IF

IF ( l_leaf_temp .AND. leaf_temp_iter < 1 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp needs leaf_temp_iter >= 1')
END IF

IF ( l_leaf_temp .AND. ( leaf_shelter <= 0.0 .OR. leaf_aero_model < 0 .OR.   &
                         leaf_aero_model > 2 ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp needs leaf_shelter > 0 and leaf_aero_model ' //    &
               '0, 1 or 2')
END IF

IF ( l_leaf_temp .AND. l_leaf_coexp_lai .AND. leaf_aero_model /= 1 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_coexp_lai is only used with leaf_aero_model = 1')
END IF

IF ( l_leaf_temp .AND. leaf_width <= 0.0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp needs leaf_width > 0')
END IF

IF ( l_som_vcmax_psi .AND. ( stomata_model /= stomata_profit_max .OR.         &
                             photo_model /= photo_farquhar ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_vcmax_psi requires stomata_model=4 and photo_model=2')
END IF

IF ( l_som_nsl .AND. ( stomata_model /= stomata_profit_max .OR.               &
                       som_ci_search /= som_ci_bounded ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_nsl requires stomata_model=4 and som_ci_search=2')
END IF

IF ( l_aggregate .AND. ANY(l_vegdrag_pft) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Vegetative drag scheme cannot be used with the ' //            &
               'aggregated surface scheme (i.e. l_aggregate = true)')
ELSE
  ! Copy the values for the given number of pfts across
  l_vegdrag_surft(1:npft) = l_vegdrag_pft(1:npft)
END IF

! Check crop options are sensible
IF ( l_crop .AND. l_red ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Crop model and RED are incompatible')
END IF

IF ( l_red .AND. .NOT. l_triffid ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'RED needs to be run with l_triffid=T')
END IF

END SUBROUTINE check_jules_vegetation
#endif

SUBROUTINE print_nlist_jules_vegetation()

USE jules_print_mgr, ONLY: jules_print
USE jules_surface_types_mod, ONLY: npft

IMPLICIT NONE

CHARACTER(LEN=50000) :: lineBuffer

CALL jules_print('jules_vegetation_mod',                                       &
                 'Contents of namelist jules_vegetation')

WRITE(lineBuffer,*)' l_nrun_mid_trif = ', l_nrun_mid_trif
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_init_accum = ', l_trif_init_accum
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_phenol = ',l_phenol
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_triffid = ',l_triffid
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_eq = ',l_trif_eq
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_veg_compete = ',l_veg_compete
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_ht_compete = ',l_ht_compete
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_crop = ',l_trif_crop
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_biocrop = ',l_trif_biocrop
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trif_fire = ',l_trif_fire
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_trait_phys = ',l_trait_phys
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_landuse = ',l_landuse
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_ag_expand = ',l_ag_expand
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_leaf_n_resp_fix = ',l_leaf_n_resp_fix
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_stem_resp_fix = ',l_stem_resp_fix
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_scale_resp_pm = ',l_scale_resp_pm
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_nitrogen = ',l_nitrogen
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_vegcan_soilfx = ',l_vegcan_soilfx
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' phenol_period = ',phenol_period
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' triffid_period = ',triffid_period
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_bvoc_emis = ',l_bvoc_emis
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_o3_damage = ',l_o3_damage
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_vegdrag_pft = ',l_vegdrag_pft(1:npft)
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_rsl_scalar = ',l_rsl_scalar
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' can_model = ',can_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' can_rad_mod = ',can_rad_mod
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' ilayers = ',ilayers
CALL jules_print('jules_vegetation_mod',lineBuffer)

!JBaguley
WRITE(lineBuffer,*) ' leaf_flux_mod = ', leaf_flux_mod
CALL jules_print('jules_vegetation_mod',lineBuffer)

!JBaguley
WRITE(lineBuffer,*) ' som_base_parm = ', som_base_parm
CALL jules_print('jules_vegetation_mod',lineBuffer)

!JBaguley
WRITE(lineBuffer,*) ' som_n_sample = ', som_n_sample
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' som_n_ci_golden_iter = ', som_n_ci_golden_iter
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_skip_search_wellwatered = ',                       &
                    l_som_skip_search_wellwatered
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' som_ci_search = ', som_ci_search
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_supply_limit = ', l_som_supply_limit
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_root_supply = ', l_som_root_supply
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_nsl = ', l_som_nsl
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_vcmax_psi = ', l_som_vcmax_psi
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_leaf_temp = ', l_leaf_temp
CALL jules_print('jules_vegetation_mod',lineBuffer)

IF ( l_leaf_temp ) THEN
  WRITE(lineBuffer,*) ' leaf_width = ', leaf_width
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_temp_iter = ', leaf_temp_iter
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' l_leaf_temp_gc_eq = ', l_leaf_temp_gc_eq
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_aero_model = ', leaf_aero_model
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_shelter = ', leaf_shelter
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' l_leaf_coexp_lai = ', l_leaf_coexp_lai
  CALL jules_print('jules_vegetation_mod',lineBuffer)
END IF

WRITE(lineBuffer,*) ' l_som_plant_segments = ', l_som_plant_segments
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_gain_gross = ', l_som_gain_gross
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_cuticular_floor = ', l_som_cuticular_floor
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_gravity = ', l_som_gravity
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' som_hc_negligible_tol = ', som_hc_negligible_tol
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' som_leaf_resist_frac = ', som_leaf_resist_frac
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' som_gl_max = ', som_gl_max
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' light_curvature_fvcb = ', light_curvature_fvcb
CALL jules_print('jules_vegetation_mod',lineBuffer)

!JBaguley
WRITE(lineBuffer,*) ' som_psi_solver = ', som_psi_solver
CALL jules_print('jules_vegetation_mod',lineBuffer)

!JBaguley
WRITE(lineBuffer,*) ' som_profit_model = ', som_profit_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' frac_min = ',frac_min
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' frac_seed = ',frac_seed
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' pow = ',pow
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' cd_leaf = ',cd_leaf
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c1_usuh = ',c1_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c2_usuh = ',c2_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' c3_usuh = ',c3_usuh
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' stanton_leaf = ',stanton_leaf
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' photo_model = ',photo_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' photo_acclim_model = ',photo_acclim_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' photo_act_model = ',photo_act_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' photo_jv_model = ',photo_jv_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' n_alloc_jmax = ',n_alloc_jmax
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' n_alloc_vcmax = ',n_alloc_vcmax
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' dsj_coef = ',dsj_coef
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' dsv_coef = ',dsv_coef
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' jv25_coef = ',jv25_coef
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' act_j_coef = ',act_j_coef
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' act_v_coef = ',act_v_coef
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' n_day_photo_acclim = ',n_day_photo_acclim
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' stomata_model = ',stomata_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_sugar = ',l_sugar
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*)' l_red = ',l_red
CALL jules_print('jules_vegetation_mod',lineBuffer)

CALL jules_print('jules_vegetation_mod',                                       &
    '- - - - - - end of namelist - - - - - -')

END SUBROUTINE print_nlist_jules_vegetation

#if defined(UM_JULES) && !defined(LFRIC)

SUBROUTINE read_nml_jules_vegetation (unitnumber)

! Description:
!  Read the JULES_VEGETATION namelist

USE setup_namelist,   ONLY: setup_nml_type
USE check_iostat_mod, ONLY: check_iostat
USE UM_parcore,       ONLY: mype


USE parkind1,         ONLY: jprb, jpim
USE yomhook,          ONLY: lhook, dr_hook

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb) :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_VEGETATION'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 21 ! +2 leaf_temp_iter/leaf_aero_model, was 16, +1 for som_n_ci_golden_iter,
                                 ! +2 for som_psi_solver/som_ci_search
INTEGER, PARAMETER :: n_real = 17 + (n_photo_coef * 5) ! +2 leaf_width/shelter, +4 for
                                  ! som_hc_negligible_tol/som_leaf_resist_frac/
                                  ! som_gl_max/light_curvature_fvcb
INTEGER, PARAMETER :: n_log = 42 + npft_max ! +1 l_som_vcmax_psi, +1 l_leaf_coexp_lai,
                                  ! +2 l_leaf_temp(_gc_eq), +1 for l_som_fast,
                                  ! +1 for
                                  ! l_som_gain_gross, +1 for
                                  ! l_som_cuticular_floor, +1 for
                                  ! l_som_gravity, +1 for
                                  ! l_som_skip_search_wellwatered, +1 for
                                  ! l_som_supply_limit, +1 for
                                  ! l_som_root_supply, +1 for l_som_nsl,
                                  ! +1 for
                                  ! l_som_plant_segments (trunk vn7.9: 29)

TYPE :: my_namelist
  SEQUENCE
  INTEGER :: phenol_period
  INTEGER :: triffid_period
  INTEGER :: can_model
  INTEGER :: can_rad_mod
  INTEGER :: leaf_temp_iter
  INTEGER :: leaf_aero_model
  INTEGER :: ilayers
  INTEGER :: leaf_flux_mod !JBaguley
  INTEGER :: som_base_parm !JBaguley
  INTEGER :: som_n_sample !JBaguley
  INTEGER :: som_n_ci_golden_iter
  INTEGER :: som_psi_aprox_method !JBaguley
  INTEGER :: som_psi_solver
  INTEGER :: som_ci_search
  INTEGER :: som_profit_model !JBaguley
  INTEGER :: ignition_method
  INTEGER :: photo_acclim_model
  INTEGER :: photo_act_model
  INTEGER :: photo_jv_model
  INTEGER :: photo_model
  INTEGER :: stomata_model
  REAL(KIND=real_jlslsm) :: frac_min
  REAL(KIND=real_jlslsm) :: frac_seed
  REAL(KIND=real_jlslsm) :: pow
  REAL(KIND=real_jlslsm) :: cd_leaf
  REAL(KIND=real_jlslsm) :: c1_usuh
  REAL(KIND=real_jlslsm) :: c2_usuh
  REAL(KIND=real_jlslsm) :: c3_usuh
  REAL(KIND=real_jlslsm) :: dsj_coef(n_photo_coef)
  REAL(KIND=real_jlslsm) :: dsv_coef(n_photo_coef)
  REAL(KIND=real_jlslsm) :: jv25_coef(n_photo_coef)
  REAL(KIND=real_jlslsm) :: act_j_coef(n_photo_coef)
  REAL(KIND=real_jlslsm) :: act_v_coef(n_photo_coef)
  REAL(KIND=real_jlslsm) :: n_alloc_jmax
  REAL(KIND=real_jlslsm) :: n_alloc_vcmax
  REAL(KIND=real_jlslsm) :: n_day_photo_acclim
  REAL(KIND=real_jlslsm) :: stanton_leaf
  REAL(KIND=real_jlslsm) :: som_hc_negligible_tol
  REAL(KIND=real_jlslsm) :: som_leaf_resist_frac
  REAL(KIND=real_jlslsm) :: som_gl_max
  REAL(KIND=real_jlslsm) :: light_curvature_fvcb
  REAL(KIND=real_jlslsm) :: leaf_width
  REAL(KIND=real_jlslsm) :: leaf_shelter
  LOGICAL :: l_som_skip_search_wellwatered
  LOGICAL :: l_som_fast
  LOGICAL :: l_som_supply_limit
  LOGICAL :: l_som_root_supply
  LOGICAL :: l_som_nsl
  LOGICAL :: l_som_vcmax_psi
  LOGICAL :: l_leaf_temp
  LOGICAL :: l_leaf_temp_gc_eq
  LOGICAL :: l_leaf_coexp_lai
  LOGICAL :: l_som_plant_segments
  LOGICAL :: l_som_gain_gross
  LOGICAL :: l_som_cuticular_floor
  LOGICAL :: l_som_gravity
  LOGICAL :: l_nrun_mid_trif
  LOGICAL :: l_trif_init_accum
  LOGICAL :: l_phenol
  LOGICAL :: l_triffid
  LOGICAL :: l_trif_eq
  LOGICAL :: l_veg_compete
  LOGICAL :: l_bvoc_emis
  LOGICAL :: l_o3_damage
  LOGICAL :: l_prescsow
  LOGICAL :: l_croprotate
  LOGICAL :: l_trait_phys
  LOGICAL :: l_ht_compete
  LOGICAL :: l_trif_crop
  LOGICAL :: l_trif_biocrop
  LOGICAL :: l_ag_expand
  LOGICAL :: l_trif_fire
  LOGICAL :: l_landuse
  LOGICAL :: l_nitrogen
  LOGICAL :: l_recon
  LOGICAL :: l_leaf_n_resp_fix
  LOGICAL :: l_stem_resp_fix
  LOGICAL :: l_scale_resp_pm
  LOGICAL :: l_vegcan_soilfx
  LOGICAL :: l_inferno
  LOGICAL :: l_vegdrag_pft(npft_max)
  LOGICAL :: l_rsl_scalar
  LOGICAL :: l_spec_veg_z0
  LOGICAL :: l_limit_canhc
  LOGICAL :: l_sugar
  LOGICAL :: l_red
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,               &
                    n_real_in = n_real, n_log_in = n_log)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_vegetation, IOSTAT = errorstatus,       &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_vegetation", iomessage)

  my_nml % phenol_period   = phenol_period
  my_nml % triffid_period  = triffid_period
  my_nml % can_model       = can_model
  my_nml % can_rad_mod     = can_rad_mod
  my_nml % ilayers         = ilayers
  my_nml % leaf_flux_mod   = leaf_flux_mod !JBaguley
  my_nml % som_base_parm   = som_base_parm !JBaguley
  my_nml % som_n_sample    = som_n_sample  !JBaguley
  my_nml % som_n_ci_golden_iter = som_n_ci_golden_iter
  my_nml % som_psi_aprox_method = som_psi_aprox_method !JBaguley
  my_nml % som_psi_solver = som_psi_solver
  my_nml % som_ci_search = som_ci_search
  my_nml % som_profit_model = som_profit_model !JBaguley
  my_nml % ignition_method = ignition_method
  my_nml % photo_acclim_model = photo_acclim_model
  my_nml % photo_act_model = photo_act_model
  my_nml % photo_jv_model     = photo_jv_model
  my_nml % photo_model     = photo_model
  my_nml % stomata_model   = stomata_model
  my_nml % frac_min        = frac_min
  my_nml % frac_seed       = frac_seed
  my_nml % pow             = pow
  my_nml % cd_leaf         = cd_leaf
  my_nml % c1_usuh         = c1_usuh
  my_nml % c2_usuh         = c2_usuh
  my_nml % c3_usuh         = c3_usuh
  my_nml % dsj_coef        = dsj_coef
  my_nml % dsv_coef        = dsv_coef
  my_nml % jv25_coef       = jv25_coef
  my_nml % act_j_coef      = act_j_coef
  my_nml % act_v_coef      = act_v_coef
  my_nml % n_alloc_jmax    = n_alloc_jmax
  my_nml % n_alloc_vcmax   = n_alloc_vcmax
  my_nml % n_day_photo_acclim = n_day_photo_acclim
  my_nml % stanton_leaf    = stanton_leaf
  my_nml % som_hc_negligible_tol = som_hc_negligible_tol
  my_nml % som_leaf_resist_frac = som_leaf_resist_frac
  my_nml % som_gl_max = som_gl_max
  my_nml % light_curvature_fvcb = light_curvature_fvcb
  my_nml % l_som_skip_search_wellwatered = l_som_skip_search_wellwatered
  my_nml % l_som_fast = l_som_fast
  my_nml % l_som_supply_limit = l_som_supply_limit
  my_nml % l_som_root_supply = l_som_root_supply
  my_nml % l_som_nsl = l_som_nsl
  my_nml % l_som_vcmax_psi = l_som_vcmax_psi
  my_nml % l_leaf_temp = l_leaf_temp
  my_nml % l_leaf_temp_gc_eq = l_leaf_temp_gc_eq
  my_nml % l_leaf_coexp_lai = l_leaf_coexp_lai
  my_nml % leaf_aero_model = leaf_aero_model
  my_nml % leaf_shelter = leaf_shelter
  my_nml % leaf_temp_iter = leaf_temp_iter
  my_nml % leaf_width = leaf_width
  my_nml % l_som_plant_segments = l_som_plant_segments
  my_nml % l_som_gain_gross = l_som_gain_gross
  my_nml % l_som_cuticular_floor = l_som_cuticular_floor
  my_nml % l_som_gravity = l_som_gravity
  my_nml % l_nrun_mid_trif = l_nrun_mid_trif
  my_nml % l_trif_init_accum   = l_trif_init_accum
  my_nml % l_phenol        = l_phenol
  my_nml % l_triffid       = l_triffid
  my_nml % l_trif_eq       = l_trif_eq
  my_nml % l_veg_compete   = l_veg_compete
  my_nml % l_bvoc_emis     = l_bvoc_emis
  my_nml % l_o3_damage     = l_o3_damage
  my_nml % l_prescsow      = l_prescsow
  my_nml % l_croprotate    = l_croprotate
  my_nml % l_trait_phys    = l_trait_phys
  my_nml % l_ht_compete    = l_ht_compete
  my_nml % l_trif_crop     = l_trif_crop
  my_nml % l_trif_biocrop  = l_trif_biocrop
  my_nml % l_ag_expand     = l_ag_expand
  my_nml % l_trif_fire     = l_trif_fire
  my_nml % l_landuse       = l_landuse
  my_nml % l_nitrogen      = l_nitrogen
  my_nml % l_recon         = l_recon
  my_nml % l_leaf_n_resp_fix = l_leaf_n_resp_fix
  my_nml % l_stem_resp_fix = l_stem_resp_fix
  my_nml % l_scale_resp_pm = l_scale_resp_pm
  my_nml % l_vegcan_soilfx = l_vegcan_soilfx
  my_nml % l_inferno       = l_inferno
  my_nml % l_vegdrag_pft   = l_vegdrag_pft
  my_nml % l_rsl_scalar    = l_rsl_scalar
  my_nml % l_spec_veg_z0   = l_spec_veg_z0
  my_nml % l_limit_canhc   = l_limit_canhc
  my_nml % l_sugar         = l_sugar
  my_nml % l_red           = l_red
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  phenol_period   = my_nml % phenol_period
  triffid_period  = my_nml % triffid_period
  can_model       = my_nml % can_model
  can_rad_mod     = my_nml % can_rad_mod
  ilayers         = my_nml % ilayers
  leaf_flux_mod   = my_nml % leaf_flux_mod !JBaguley
  som_base_parm   = my_nml % som_base_parm !JBaguley
  som_n_sample    = my_nml % som_n_sample  !JBaguley
  som_n_ci_golden_iter = my_nml % som_n_ci_golden_iter
  som_psi_aprox_method = my_nml % som_psi_aprox_method !JBaguley
  som_psi_solver = my_nml % som_psi_solver
  som_ci_search = my_nml % som_ci_search
  som_profit_model = my_nml % som_profit_model !JBaguley
  ignition_method = my_nml % ignition_method
  photo_acclim_model = my_nml % photo_acclim_model
  photo_act_model = my_nml % photo_act_model
  photo_jv_model     = my_nml % photo_jv_model
  photo_model     = my_nml % photo_model
  stomata_model   = my_nml % stomata_model
  frac_min        = my_nml % frac_min
  frac_seed       = my_nml % frac_seed
  pow             = my_nml % pow
  cd_leaf         = my_nml % cd_leaf
  c1_usuh         = my_nml % c1_usuh
  c2_usuh         = my_nml % c2_usuh
  c3_usuh         = my_nml % c3_usuh
  dsj_coef        = my_nml % dsj_coef
  dsv_coef        = my_nml % dsv_coef
  jv25_coef       = my_nml % jv25_coef
  act_j_coef      = my_nml % act_j_coef
  act_v_coef      = my_nml % act_v_coef
  n_alloc_jmax    = my_nml % n_alloc_jmax
  n_alloc_vcmax   = my_nml % n_alloc_vcmax
  n_day_photo_acclim = my_nml % n_day_photo_acclim
  stanton_leaf    = my_nml % stanton_leaf
  som_hc_negligible_tol = my_nml % som_hc_negligible_tol
  som_leaf_resist_frac = my_nml % som_leaf_resist_frac
  som_gl_max = my_nml % som_gl_max
  light_curvature_fvcb = my_nml % light_curvature_fvcb
  l_som_skip_search_wellwatered = my_nml % l_som_skip_search_wellwatered
  l_som_fast = my_nml % l_som_fast
  l_som_supply_limit = my_nml % l_som_supply_limit
  l_som_root_supply = my_nml % l_som_root_supply
  l_som_nsl = my_nml % l_som_nsl
  l_som_vcmax_psi = my_nml % l_som_vcmax_psi
  l_leaf_temp = my_nml % l_leaf_temp
  l_leaf_temp_gc_eq = my_nml % l_leaf_temp_gc_eq
  l_leaf_coexp_lai = my_nml % l_leaf_coexp_lai
  leaf_aero_model = my_nml % leaf_aero_model
  leaf_shelter = my_nml % leaf_shelter
  leaf_temp_iter = my_nml % leaf_temp_iter
  leaf_width = my_nml % leaf_width
  l_som_plant_segments = my_nml % l_som_plant_segments
  l_som_gain_gross = my_nml % l_som_gain_gross
  l_som_cuticular_floor = my_nml % l_som_cuticular_floor
  l_som_gravity = my_nml % l_som_gravity
  l_nrun_mid_trif = my_nml % l_nrun_mid_trif
  l_trif_init_accum = my_nml % l_trif_init_accum
  l_phenol        = my_nml % l_phenol
  l_triffid       = my_nml % l_triffid
  l_trif_eq       = my_nml % l_trif_eq
  l_veg_compete   = my_nml % l_veg_compete
  l_bvoc_emis     = my_nml % l_bvoc_emis
  l_o3_damage     = my_nml % l_o3_damage
  l_prescsow      = my_nml % l_prescsow
  l_croprotate    = my_nml % l_croprotate
  l_trait_phys    = my_nml % l_trait_phys
  l_ht_compete    = my_nml % l_ht_compete
  l_trif_crop     = my_nml % l_trif_crop
  l_trif_biocrop  = my_nml % l_trif_biocrop
  l_ag_expand     = my_nml % l_ag_expand
  l_trif_fire     = my_nml % l_trif_fire
  l_landuse       = my_nml % l_landuse
  l_nitrogen      = my_nml % l_nitrogen
  l_recon         = my_nml % l_recon
  l_leaf_n_resp_fix = my_nml % l_leaf_n_resp_fix
  l_stem_resp_fix = my_nml % l_stem_resp_fix
  l_scale_resp_pm = my_nml % l_scale_resp_pm
  l_vegcan_soilfx = my_nml % l_vegcan_soilfx
  l_inferno       = my_nml % l_inferno
  l_vegdrag_pft   = my_nml % l_vegdrag_pft
  l_rsl_scalar    = my_nml % l_rsl_scalar
  l_spec_veg_z0   = my_nml % l_spec_veg_z0
  l_limit_canhc   = my_nml % l_limit_canhc
  l_sugar         = my_nml % l_sugar
  l_red           = my_nml % l_red
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_vegetation
#endif

END MODULE jules_vegetation_mod
