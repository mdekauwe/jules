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
    ! Use the model of Medlyn et al. (2011, Global Change Biology) - see
    !   Eqn.11, doi: 10.1111/j.1365-2486.2010.02375.x.
  stomata_sox_analytical = 3,                                                  &
    ! Use the semi-analytical version of the SOX model (Eller et al 2020)
    ! doi: 10.1111/nph.16419 - Eqns. 4 & 5
  stomata_profit_max = 4,                                                      &
    ! Stomatal optimisation, profit maximisation (Sperry et al. 2017; De
    ! Kauwe et al. 2022, New Phytologist; Baguley et al. 2026,
    ! Biogeosciences); see stom_opt_jls_mod.
  stomata_sox_opt = 5,                                                         &
    ! SOX_opt: the SOX optimisation of Eller et al. (2018) solved as Sabot
    ! et al. (2022) do (their SOXopt), maximising
    !   A_n (k(psi_bar) - kcrit)/(k_ref - kcrit),
    ! k of the whole path at psi_bar, the mean of the root zone and leaf
    ! water potentials; see stom_opt_jls_mod. The halfway point is Eller et
    ! al.'s (2018, Eqs 2.7-2.8: f(psi_c,mid), psi_c,mid = (psi_c,pd +
    ! psi_c)/2, their approximation of the path-averaged conductance of
    ! Sperry et al.); Sabot et al. (2022, Eq. 16) write k at psi_leaf
    ! instead, which we do not follow. From Sabot et al. we take the kcrit
    ! subtraction (Eller: f/k_max); Eller's cost is also xylem-only, with
    ! the soil acting through the root psi, whereas k here includes the
    ! soil-to-root link. Not the semi-analytical SOX of Eller et al. (2020)
    ! (stomata_sox_analytical = 3).
  stomata_cmax = 6,                                                            &
    ! Carbon maximisation (CMax) of Wolf et al. (2016, PNAS 113: E7222) in
    ! the form of Anderegg et al. (2018): maximise net carbon gain
    !   A_n - Theta(psi_leaf),  Theta = a/2 psi_leaf^2 + b |psi_leaf|,
    ! a carbon cost of low leaf water potential in absolute units (no
    ! normalisation), on the profit-max hydraulics; see stom_opt_jls_mod.
  stomata_cgain = 7,                                                           &
    ! CGain of Lu et al. (2020) (Sabot et al. 2022, Eq. 15): maximise
    !   A_n - varpi (k_max - k(psi_leaf)) / k_max,
    ! the profit max's hydraulic cost in absolute carbon units (varpi,
    ! cgain_varpi_io) instead of the gain normalised by the instantaneous
    ! A_max; see stom_opt_jls_mod.
  stomata_desica = 8,                                                          &
    ! DESICA (De Kauwe et al. 2020): Tuzet et al. (2003) stomatal closure on
    ! leaf water potential, with leaf and stem water potentials (and stem
    ! storage) from the plant hydraulics of Xu et al. (2016),
    ! doi: 10.1111/nph.14009 (Notes S1). See desica_jls_mod.
  stomata_supply_loss = 9,                                                     &
    ! Carbon-coupled demand (ci = sl_cica_well_watered ca: Medlyn with D
    ! fixed, so VPD acts only through the hydraulics, not twice) regulated
    ! by the supply-loss rule of Sperry et al. (2016, Eqn 5) on the
    ! profit-max hydraulics; see stom_opt_supply_loss in stom_opt_jls_mod.
  stomata_cap = 10,                                                            &
    ! CAP of Dewar et al. (2018, New Phytol 217: 571; Sabot et al. 2022,
    ! Eqs 19-21): no stomatal cost; photosynthesis is limited
    ! nonstomatally by the leaf water potential, gross A x
    !   f = 1 - psi_leaf / psi_nsl0   (psi_nsl_onset = 0; clipped to [0, 1])
    ! (Dewar et al. scale Vcmax and J, i.e. gross A), and the stomata
    ! maximise this A_n on the profit-max hydraulics (k > kcrit still
    ! bounds it). Uses the l_som_nsl machinery, which it switches on;
    ! psi_nsl0_io is the parameter (psi_phi,lim of Sabot et al.). See
    ! stom_opt_bounded_search.
  stomata_pm2 = 11
    ! ProfitMax2, the "new model" of Wang et al. (2020, New Phytol 227:
    ! 311, Eqn 11): maximise
    !   A_n (1 - E / E_crit),  i.e. Theta = A_n E / E_crit,
    ! with E_crit the transpiration at which the soil-to-leaf conductance
    ! falls to kcrit (the edge of the supply function), on the same path
    ! and basis as the profit max. No cost parameter: the marginal water
    ! cost A_n / (E_crit - E) is > 0 at any E, without new cavitation (the
    ! profit max's cost is ~0 in wet soil), scales with the current A_n and
    ! rises as drought (or drought damage) shrinks E_crit. On the bounded
    ! Ci search; see stom_opt_bounded_search.
    !
    ! MDK, 22/10/26
!
! Renumbered 7 Oct 2026 (branch cmax): SOX_opt 6 -> 5, DESICA 5 -> 8, CMax 6,
! CGain 7, supply-loss 9. Namelists from before that use 5 for DESICA or 6
! for SOX_opt must be changed.
! CAP (10) added 8 Oct 2026 (branch stomata_fixes).
! ProfitMax2 (11) added 22 Oct 2026 (branch pm_water_cost). MDK, 22/10/26
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

! Paramiters identifying different xylem imparement models.
! These should have unique values. JBaguley
INTEGER, PARAMETER ::                                                          &
  xylem_impairment_none = 0,                                                   &
    ! No xylem imparement model is used.
  xylem_impairment_kmax = 1,                                                   &
    ! RETIRED (2026-10-10): reduced kmax with refitted b, c (on/off at
    ! UK-Ham). Number kept; the code is on xylem_impairment at 130a575.
  xylem_impairment_whole_trunk = 2,                                            &
    ! Conductance along the xylem limited by the historic minimum water
    ! potential profile.
  xylem_impairment_memory = 3,                                                 &
    ! Embolism memory with recovery: the vulnerability curve is capped at
    ! kmax * (1 - PLC_mem), where PLC_mem is the loss of conductivity at the
    ! most negative damage-driving water potential experienced, relaxing
    ! back when the plant rehydrates.
  xylem_impairment_kmax_refit = 4
    ! RETIRED (2026-10-10): a Weibull refitted to the memory model's capped
    ! curve (Mackay et al. 2015 / TREES), an approximation of model 3, which
    ! uses the capped curve exactly. Number kept; the code is on
    ! xylem_impairment at 130a575.

! Water potential driving embolism in the memory impairment model
! (ximpair_psi_driver).
INTEGER, PARAMETER ::                                                          &
  ximpair_driver_leaf = 1,                                                     &
    ! Leaf water potential.
  ximpair_driver_mean = 2,                                                     &
    ! Mean of the leaf and root zone water potentials.
  ximpair_driver_root = 3,                                                     &
    ! Root zone water potential (~ predawn; Mackay et al. 2015).
  ximpair_driver_stem = 4
    ! Stem water potential: the leaf-side end of the stem segment
    ! (l_som_plant_segments), i.e. the root zone potential less the drop
    ! through the root and stem segments at the chosen transpiration.

! Solvers for the leaf water potential given the transpiration rate
! (som_psi_solver). These should have unique values. JBaguley
! (1, a zeroth-order Taylor expansion of the conductance about the mid-path
! psi, was retired: it overstated k near E_crit. MDK, 25/10/26)
INTEGER, PARAMETER ::                                                          &
  psi_solver_newton = 2,                                                       &
    ! Use the Newton-Raphson method to approximate leaf water potential
     ! from transpiration rate.
  psi_solver_lut = 3
    ! Invert a per-PFT lookup table of the supply function (the integral
    ! of the vulnerability curve) to get leaf water potential directly
    ! from transpiration rate (cumulative Weibull or SOX conductance).

! Diagnostic leaf 13C discrimination (discrimination_model), computed from
! the ci of whichever stomatal model is in use. Diagnostic only: it does not
! feed back on any flux.
INTEGER, PARAMETER ::                                                          &
  discrimination_off = 0,                                                      &
    ! No 13C discrimination diagnostic.
  discrimination_simple = 1,                                                   &
    ! Farquhar et al. (1982) simple form, a + (b - a) ci/ca, b = 27 permil.
  discrimination_photoresp = 2
    ! As 1 with b = 28 permil and the photorespiratory term -f Gamma*/ca
    ! (Ubierna & Farquhar 2014), as in Lavergne et al. (2022) and
    ! Palmer et al. (2022) without the mesophyll term.

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
    ! gives the edge directly. sox_opt_profit_model uses it too (A_n (k - kcrit)
    ! is unimodal in Ci); SOX_opt has no xylem impairment (rejected at
    ! start-up). MDK, 24/10/26

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
  sox_opt_profit_model = 2,                                                    &
    ! SOX_opt (set from stomata_model = 5) to determine the optimal stomatal
    ! conductance.
  supply_loss_model = 3,                                                       &
    ! Not a profit: the supply-loss regulation of Sperry et al. (2016) (set
    ! from stomata_model = 9, see stomata_supply_loss).
    !      Profit(psi) = CarbonGain(psi) * (1 - HydraulicCost(psi))
    !
    !      CarbonGain(psi) = Anet(psi)
    !
    !      HydraulicCost(psi) = 1 - (k(psi) / kmax)
  cmax_profit_model = 4,                                                       &
    ! CMax (Wolf et al. 2016; set from stomata_model = 6), in absolute
    ! carbon units (mol CO2 m-2 s-1):
    !      Profit(psi) = Anet(psi) - Theta(psi_leaf)
    !      Theta(psi)  = cmax_a/2 psi^2 + cmax_b |psi|   (psi in MPa)
  cgain_profit_model = 5,                                                      &
    ! CGain (Lu et al. 2020; set from stomata_model = 7), mol CO2 m-2 s-1:
    !      Profit(psi) = Anet(psi) - varpi (k_max - k(psi_leaf)) / k_max
    ! with k the whole-path conductance -dE/dpsi_leaf (soil link included)
    ! and k_max the plant's unstressed conductance (psi = 0, without the
    ! soil link, so that drying soil adds to the cost rather than shrinking
    ! the reference; Lu et al. 2020: the plant's k_max).
  cap_profit_model = 6,                                                        &
    ! CAP (Dewar et al. 2018; set from stomata_model = 10):
    !      Profit(psi) = Anet(psi) with gross A x f(psi_leaf)
  pm2_profit_model = 7
    ! ProfitMax2 (Wang et al. 2020, Eqn 11; set from stomata_model = 11),
    ! mol CO2 m-2 s-1:
    !      Profit(E) = Anet(E) (1 - E / E_crit)
    !
    ! MDK, 22/10/26

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
   l_triffid_nsc = .FALSE.,                                                    &
       ! TRIFFID storage pool (Trugman et al. 2018, Ecol. Lett. 21: 1552):
       ! NPP enters a non-structural carbon store and TRIFFID grows from the
       ! carbon the store releases, W = W_max C_a / (k_c C_max + C_a)
       ! (Michaelis-Menten loading), not from the NPP of the step. The
       ! store buffers carbon deficits (it pays for negative NPP) and, with
       ! l_ximpair_rec_nsc, pays for xylem repair. Not with l_sugar or
       ! l_trif_eq. Without l_triffid it runs daily with the prescribed LAI
       ! and height (store and repair only). See triffid_nsc_mod and
       ! NOTES_nsc_storage.md.
       ! MDK, 17/10/26
   l_trif_hydmort = .FALSE.,                                                   &
       ! Hydraulic (drought) mortality in TRIFFID: the stand's lethal losses
       ! of conductance are spread logistically about hydmort_mu (width
       ! hydmort_width), so the fraction of the stand dead is S(p), with p
       ! the lasting loss of conductance of the memory impairment model
       ! (ximpair_lock). A ratchet: each TRIFFID step kills only the part of
       ! S the stand has not already lost, at the period's peak damage. The
       ! killed area goes to litter (lotka with l_veg_compete; otherwise it
       ! becomes bare soil). See triffid_hydmort_mod.
       ! MDK, 21/10/26

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
  som_n_ci_golden_iter = 16,                                                   &
      ! Maximum golden-section iterations of the bounded Ci search, on
      ! top of the initial two-point bracket setup.
  som_psi_solver = psi_solver_newton,                                          &
      ! Leaf water potential from transpiration rate: 2 Newton-Raphson,
      ! 3 lookup table of the supply function. (1, the Taylor series, was
      ! the default before 24/10/26 and is retired: it held k at the
      ! mid-path psi, so near failure it overstated k several-fold and could
      ! diverge. 1 stops the run.)
      !
      ! MDK, 25/10/26
  som_ci_search = som_ci_bounded,                                              &
      ! Ci search of the stomatal optimisation: 1 flat grid, 2 bounded.
      ! Default 2 (it was 1 before 25/10/26): the bounded search is faster,
      ! not limited to the grid, and the only one coded for hold_psi,
      ! the plant capacitance store, l_leaf_temp_ci, l_som_nsl and models 6, 7,
      ! 10 and 11. The flat grid is kept as a check where the profit may
      ! not have a single peak; utils/upgrade_nml_gce.py sets 1 in old
      ! namelists that relied on the default.
      !
      ! MDK, 25/10/26
  discrimination_model = discrimination_off,                                   &
      ! Leaf 13C discrimination diagnostic: 0 off, 1 simple,
      ! 2 with photorespiration.
  som_psi_aprox_method = imdi,                                                 &
      ! RETIRED input: setting it stops the run; use som_psi_solver (same
      ! values). Still read so namelists that leave it unset run.
      ! MDK, 24/10/26
  som_profit_model = imdi
      ! Profit used by the stomatal optimisation (1: profit max, 2: SOX).
      ! Derived from stomata_model; a namelist value is ignored, with a
      ! warning (it was read with the retired leaf_flux_mod = 2). Unset
      ! (imdi) in the namelist, and 1 for models without the optimisation,
      ! as before.
      !
      ! JBaguley; MDK, 25/10/26

LOGICAL ::                                                                     &
  l_som_gain_gross = .FALSE.
      ! When .TRUE., the carbon gain in the profit-max search is normalised
      ! gross assimilation, (An + Rd)/max(An + Rd), instead of net An/max(An).
      ! Rd does not depend on gs, so it cannot be traded against water, but
      ! normalising net A stretches the gain when Rd is a large part of A
      ! (low light) and fails near the compensation point. GPP/NPP are
      ! unaffected (only the optimisation's gain changes).

LOGICAL ::                                                                     &
  l_som_cuticular_loss = .FALSE.
      ! When .TRUE., leaf water loss never falls below a cuticular floor:
      ! the canopy conductance is MAX(gs, gcuticular_io * LAI). The floor is
      ! an uncontrolled leak, so it is not part of the optimisation and adds
      ! no carbon. Under the stomatal optimisation (profit max, with the
      ! capacitance store, which it turns on): lit points only, the full
      ! leak, its transpiration the tile's flux at that conductance, paid
      ! by the store and the soil (the soil-supply cap applies) through the
      ! store's balance, and cut where the leaf would end below the psi of
      ! the critical loss of conductance (kcrit_fractional_loss); the xylem
      ! impairment update is made after the leak, so it drives the loss of
      ! conductance. DESICA: its own leak, at night too, cut to the largest
      ! flux with k > kcrit.
      !
      ! MDK, 25/10/26

LOGICAL ::                                                                     &
  l_som_capacitance = .FALSE.
      ! Plant capacitance. When .TRUE., the stomatal optimisation (profit
      ! max) has a plant water store at the stem node of the root, stem
      ! and leaf segments, of capacitance C = cap_stem times the reference
      ! leaf area of the plant conductance, carried and dumped. The
      ! optimiser sees it: the store's flow out at the stem node, q_f, is
      ! held fixed while it searches, so the soil, root and stem carry
      ! E - q_f and the leaf E. A full store (morning) leaves the stem and
      ! leaf wetter at the same E, the hydraulic cost lower and the stomata
      ! more open; a refilling store (afternoon, night) the reverse. q_f
      ! is updated between the passes of the leaf-temperature (or humidity)
      ! iteration by the store's balance at the stem node over the step
      ! (backward Euler), committed at the final E, so the water is
      ! conserved; sf_evap takes the store's flow off the soil extraction.
      ! The store drives the leaf and stem psi, and so the embolism. Needs
      ! the plant segments with the Newton solver, the bounded Ci search,
      ! at least two passes, and the big-leaf, multilayer or two-leaf
      ! canopy. On with l_som_cuticular_loss whatever its value (the leak
      ! drains the store).
      !
      ! MDK, 25/10/26

LOGICAL ::                                                                     &
  l_som_capacitance_on = .FALSE.
      ! Derived (not in the namelist): the stomatal optimisation's plant
      ! capacitance is on (l_som_capacitance or l_som_cuticular_loss, with
      ! the stomatal optimisation).
      !
      ! MDK, 23/10/26

LOGICAL ::                                                                     &
  l_som_gravity = .FALSE.
      ! When .TRUE., the stomatal optimisation (stomata_model = 4 to 7) takes
      ! the gravitational drop rho_water g h to the canopy height h off the
      ! root zone water potential, so the plant path starts from
      ! psi_root_zone - rho_water g h (0.01 MPa per m), as DESICA's psi_h.

LOGICAL ::                                                                     &
  l_som_coupled_e = .TRUE.
      ! Not a namelist input since 25 Oct 2026: always on where it is
      ! supported, and check_jules_vegetation turns it off elsewhere (see
      ! the restrictions below). (On by default since 9 Oct 2026: a
      ! recalibration with it on (FR-Pue, profit max) gave the same
      ! parameters and fit as with it off; NOTE_coupled_e_handoff_2026-10-08.)
      ! It also carries the CO2 path through the leaf boundary layer.
      !
      ! MDK, 25/10/26
      ! When .TRUE., the stomatal optimisation (stomata_model = 4 to 7, 10, 11;
      ! can_rad_mod = 1 or 7) evaluates the transpiration of each trial
      ! stomatal conductance g with the resistances it passes through,
      !   E = D g' / (1 + r_ca (g' + g'_other)),   g' = g / (1 + g r_bl),
      ! instead of E = g D_c at a humidity deficit D_c held fixed from the
      ! last pass. D is from the leaf (two-leaf with l_leaf_temp: each
      ! class's leaf temperature; otherwise tstar) to the level-1 air q1;
      ! r_bl the leaf boundary layer (two-leaf with l_leaf_temp, from the
      ! leaf energy balance's gb; 0 otherwise); r_ca the canopy air to
      ! level 1 (JULES's ra; with l_leaf_temp the leaf_aero_model's, 0 for
      ! leaf_aero_model = 0); g'_other the other leaf class and the soil
      ! evaporation beneath the canopy (soil_evap's conductance, which
      ! shares the same path; with l_leaf_temp at qs(tstar), which adds
      ! r_ca g_soil (qs(T_leaf) - qs(tstar)) to D). The fixed-D_c
      ! form makes E rise linearly with g, which overstates the marginal
      ! water cost of opening by about (1 + r g) (cf. Sabot et al. 2022,
      ! Eqs 1 and 10; Wolf et al. 2016, Eq. 6). With xylem impairment the
      ! hydraulics of each trial E are those of the impaired curve.

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

INTEGER, PARAMETER ::                                                          &
  kplant_model_param  = 0,                                                     &
  kplant_model_ground = 1
      ! Values of kplant_model.
      ! MDK, 21/10/26
REAL(KIND=real_jlslsm), PARAMETER :: k_ground_ratio_max = 10.0
      ! Cap on LAI_ref / LAI (kplant_model = 1), also used as LAI -> 0.
      ! (Was in xylem_hydraulics_CW_jls_mod.)
      !
      ! MDK, 23/10/26
INTEGER ::                                                                     &
  kplant_model = kplant_model_param
      ! Whole-plant hydraulic conductance (needs l_som_plant_segments for 1):
      ! 0 (kplant_model_param, default): kmax_pft, per unit leaf area, at
      !   every LAI, so all three segments scale with leaf area.
      ! 1 (kplant_model_ground): the root and stem segments are held per
      !   unit ground area at their kmax_pft value for a reference LAI
      !   (lai_ref_k_io in pft_params; unset or <= 0: the TRIFFID balanced
      !   LAI from canopy height), while the leaf segment scales with LAI.
      !   Segment s in {root, stem} has kmax_pft * seg_kfac_s * LAI_ref / LAI
      !   (ratio capped at k_ground_ratio_max), so losing leaves raises the
      !   supply per remaining leaf and adding leaves lowers it (as in
      !   SurEau, where trunk, branch and root conductances are Ks x area /
      !   length). Xylem impairment stays a fraction of each segment's
      !   intact conductance, so a change in LAI neither repairs nor adds
      !   damage. kmax_pft keeps its meaning at LAI = LAI_ref.
      ! MDK, 21/10/26

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
      ! RETIRED (2026-10-05): the soil-to-root resistance is now always in
      ! series in the profit-max path (l_som_rhizo_series). Still read so
      ! old namelists run; .TRUE. gives a warning and is ignored. Was:
      ! When .TRUE., the stomatal optimisation can only choose transpiration
      ! the roots can take up with the root held at root_psi_crit:
      !   E <= sum_layers soil_to_root_k * MAX(psi_soil - root_psi_crit, 0)
      !        / (rho_water g),
      ! so uptake falls as the soil (and rhizosphere conductance) dries and
      ! stops at root_psi_crit, which bounds the soil water potential the
      ! plant can produce. Uses the same cap on gl as l_som_supply_limit
      ! (with which it combines, taking the smaller supply).
LOGICAL ::                                                                     &
  l_som_rhizo_series = .FALSE.
      ! Not a namelist switch: set .TRUE. in check_jules_vegetation for the
      ! profit max and DESICA. The plant hydraulic path then includes the
      ! soil-to-root (rhizosphere) conductance ksr in series with the plant,
      ! as in SPA, MAESPA (KTOT = 1/(TOTSOILRES + 1/PLANTK)) and Bonan et al.
      ! (2014): ksr is soil_to_root_k summed over the layers (parallel; soil
      ! part only, Bonan eq. A23), kept apart from the xylem kmax (the
      ! vulnerability curve, and any impairment, act on kmax only). Profit max:
      ! the root water potential falls with the flux,
      ! psi_root = psi_src - E/ksr, with ksr shared between the leaf paths by
      ! their share of the undamaged plant conductance (set_ksr_path); the
      ! marginal conductance used by the hydraulic cost and kcrit is that of the
      ! whole soil-to-leaf path (with l_som_plant_segments the soil link comes
      ! ahead of the root segment). DESICA: ksr is in series with the root-side
      ! conductance. Root (radial), stem and leaf resistance stay in kmax, so
      ! nothing is counted twice. psi_root_zone is the conductance-weighted soil
      ! psi (smc_ext), and the layer extraction weights are k_i (psi_i -
      ! psi_root) at the actual root potential (physiol). Applied for PFTs with
      ! fsmc_mod = 2 only (smc_ext computes soil_to_root_k only there). It
      ! replaces l_som_root_supply.
LOGICAL ::                                                                     &
  l_som_rhizo_mfp = .FALSE.
      ! RETIRED: the nonlinear (matric flux potential) soil-to-root link made
      ! no measurable difference at FR-Pue (recalibration 10 Oct 2026; drought
      ! tests 22 Oct 2026) and was removed. The code is kept on branch
      ! rhizo_mfp. Still read so old namelists run; .TRUE. gives a warning
      ! and is ignored, and the linear link psi_root = psi_src - E / ksr is
      ! used.
      !
      ! MDK, 23/10/26
! (l_som_rhizo_hr, hydraulic redistribution in the l_som_rhizo_mfp root
! balance, was removed: no namelist set it. MDK, 24/10/26)
LOGICAL ::                                                                     &
  l_root_mass_fixed = .FALSE.
      ! Root mass in the soil-to-root conductance (soil_to_root_k, Bonan et
      ! al. 2014 eq. A23; used by the l_som_root_supply cap and the root
      ! zone psi weights). .FALSE. (default): the vegetation's own root
      ! carbon, the one root respiration uses in sf_stom - lma*lai_bal with
      ! l_trait_phys, else sigl*lai_bal/rmass, rootc_cpft/rmass for crops -
      ! with lai_bal from the canopy height (so dynamic with TRIFFID, fixed
      ! with prescribed height) or from RED, floored at root_mass_min.
      ! .TRUE.: the fixed min_rootc_pft (kg m-2), the earlier behaviour.

REAL(KIND=real_jlslsm) ::                                                      &
  root_mass_min = 0.05
      ! Floor on the vegetation root dry mass (kg m-2; not used with
      ! l_root_mass_fixed), so a
      ! vanishing root carbon cannot drive the soil-to-root conductance to
      ! zero (or its LOG term to a division by zero).

LOGICAL ::                                                                     &
  l_leaf_temp = .FALSE.
      ! When .TRUE., the two-leaf canopy (can_rad_mod = 7) solves a leaf
      ! energy balance for the sunlit and the shaded leaf (Penman-Monteith;
      ! Leuning et al., 1995; Wang & Leuning, 1998) and runs photosynthesis
      ! and the stomatal optimisation of each at its own leaf temperature,
      ! instead of both at the surface temperature tstar. The tile energy
      ! balance is unchanged (it still solves one tstar from gc).
      ! The multilayer canopy (can_rad_mod = 5, 6) solves the same leaf
      ! energy balance as set by ml_leaf_temp_model; its passes are
      ! leaf_temp_iter, as two-leaf.
      ! MDK, 19/10/26
INTEGER, PARAMETER ::                                                          &
  ml_leaf_temp_aggregate = 1,                                                  &
  ml_leaf_temp_layers    = 2
      ! Values of ml_leaf_temp_model.
      ! MDK, 19/10/26
INTEGER ::                                                                     &
  ml_leaf_temp_model = ml_leaf_temp_layers
      ! With l_leaf_temp and the multilayer canopy (can_rad_mod = 5, 6), the
      ! leaf energy balance:
      ! 1 (ml_leaf_temp_aggregate): one sunlit and one shaded balance (the
      !   two-leaf one) from the sunlit and shaded leaf area, absorbed
      !   radiation and stomatal conductance summed over the layers; every
      !   layer's sunlit (shaded) leaves are at that sunlit (shaded)
      !   temperature.
      ! 2 (ml_leaf_temp_layers, default): as MAESPA, a balance for the sunlit
      !   and the shaded leaf of every layer, each from its own absorbed
      !   radiation, stomatal conductance and boundary layer (the wind
      !   falling with depth, see leaf_wind_ext), all exchanging with one
      !   canopy air.
      ! Not used otherwise.
      ! MDK, 19/10/26
REAL(KIND=real_jlslsm) ::                                                      &
  leaf_width = 0.05
      ! Characteristic leaf width (m), for the leaf boundary layer
      ! conductance (l_leaf_temp).
INTEGER, PARAMETER ::                                                          &
  leaf_temp_iter = 3
      ! Number of passes of the two-leaf stomatal optimisation with
      ! l_leaf_temp (each pass is one stom_opt_mod call per leaf, then a
      ! leaf temperature update). iter (3) without l_leaf_temp. No longer
      ! in the namelist: with 1 pass the leaf temperature was never updated
      ! and stayed at the first step's value for the whole run, and every
      ! run used 3.
      !
      ! MDK, 25/10/26
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
REAL(KIND=real_jlslsm) ::                                                      &
  leaf_wind_ext = 0.0
      ! With l_leaf_temp and leaf_aero_model = 0 or 2, the extinction
      ! coefficient of the wind in the canopy per unit cumulative leaf area
      ! above (dimensionless): u(L) = u_top exp(-leaf_wind_ext L), so the
      ! forced-convection boundary-layer conductance falls as
      ! exp(-leaf_wind_ext L / 2) (MAESPA's EXTWIND, default 1 there;
      ! Wang & Leuning, 1998), u_top being the level-1 wind. 0 (default):
      ! no attenuation, every leaf at the level-1 wind. Applies to the
      ! two-leaf (can_rad_mod 7, through the wind-weighted class leaf area)
      ! and the multilayer (5, 6, layer by layer) leaf energy balance.
      ! leaf_aero_model = 1 uses CABLE's own coexp instead.
      ! MDK, 19/10/26
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
  l_leaf_temp_ci = .FALSE.
      ! With l_leaf_temp (two-leaf canopy, can_rad_mod = 7) and the bounded
      ! Ci search: the optimiser sees each candidate's own leaf temperature.
      ! For each Ci the leaf temperature comes from the leaf energy balance
      ! linearised in transpiration about the last pass,
      !   T_leaf = t_c + (Rn - lambda E) / (cp gh),
      ! (the balance of leaf_temp_update_n at the canopy air of the last
      ! pass), and the photosynthesis parameters and the leaf-to-air
      ! humidity deficit follow T_leaf (log-linear and linear in T about
      ! the last pass). So opening the stomata cools the leaf in the gain
      ! and the transpiration, as in Sperry et al. (2017), where T_leaf is
      ! computed for each E on the supply curve. The first pass of each
      ! step uses the fixed leaf temperatures (no balance yet). Not with
      ! l_som_nsl. .FALSE.: the leaf temperature is fixed during the search
      ! and updated only between passes.
      !
      ! MDK, 24/10/26
REAL(KIND=real_jlslsm) ::                                                      &
  sl_kfail_frac = 0.0
      ! Supply-loss stomata (stomata_model = 9) only: the whole-path
      ! conductance at which the supply function fails (E_crit, P_crit), as
      ! a fraction of kmax. 0 (default): kcrit, as the profit max
      ! (kcrit_fractional_loss, 5 % of kmax at 0.95). Sperry et al. (2016)
      ! use 0.0005 (0.05 %, "physiological zero").
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
      ! so no marginal cost, above it). f and psi_leaf are solved together
      ! for each Ci, and the limitation reduces both the optimiser's gain and the actual
      ! photosynthesis.
LOGICAL ::                                                                     &
  l_som_vcmax_psi = .FALSE.
      ! When .TRUE., photosynthetic capacity (Vcmax and Jmax) is down-
      ! regulated as the soil dries, by the root-zone water potential (the
      ! predawn proxy), with the Zhou et al. (2013, Agric. For. Meteorol.
      ! 182-183: 204) form f = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f -
      ! psi_rz))) (pft_params psi_vcmax_f_io, sf_vcmax_io; psi in MPa), as
      ! applied in a land surface model by De Kauwe et al. (2015,
      ! Biogeosciences 12: 7503-7518, Eq. 4). Two adjustments to that
      ! formulation here: an optional floor (psi_vcmax_fmin_io) and a fixed
      ! 5-day recovery after rewetting (see vcmax_psi_factor in sf_stom).
      ! Unlike l_som_nsl it does not act in wet soil and does not depend on
      ! the midday leaf psi. Stomatal optimisation models 4-7 and 9-11 with
      ! photo_model = 2; C4 PFTs (Collatz) get no down-regulation, as it is
      ! applied in the Farquhar branch only. Dark respiration is not scaled
      ! (it is set from Vcmax before the factor), so neither is the stem and
      ! root maintenance respiration scaled from it.
      !
      ! MDK, 25/10/26
! How the stomata respond to xylem impairment (ximpair_stomatal_response).
! In both, the stomata decide on the intact (undamaged) path; they differ in
! what is held when the damaged path cannot do what the intact one would.
INTEGER, PARAMETER ::                                                          &
  ximpair_hold_gs = 0,                                                         &
    ! Hold gs: the hydraulic cost is the intact path's conductance at the
    ! (impaired) leaf psi, relative to the intact zero-flow conductance;
    ! E, psi and A on the impaired path at the chosen Ci. The damage lowers
    ! psi_leaf at about the same gs and E (J. Baguley's scheme).
  ximpair_hold_psi = 1
    ! Hold leaf psi (Sperry et al. 2016, Plant Cell Environ 39: 2155, Fig. 3
    ! and their C++ code): the optimisation runs on the intact path and sets
    ! the target psi_leaf; E is what the impaired path supplies at that
    ! psi_leaf, and gs and A follow. The damage lowers E, gs and A at the
    ! same psi_leaf. The target does not change with past damage.
    ! (A third option, stomata optimising on the impaired path, opened the
    ! stomata after damage and was removed; see the xylem_impairment branch
    ! at 5536eb8.)
INTEGER ::                                                                     &
  ximpair_stomatal_response = ximpair_hold_gs
      ! See ximpair_hold_*. ximpair_hold_psi needs som_ci_search = 2.
LOGICAL ::                                                                     &
      ! (l_ximpair_rec_growth, the memory model's time-based recovery from
      ! new xylem grown from carbon gain, was retired 23/10/26 in favour of
      ! the store repair, l_ximpair_rec_nsc; it is in the git history.
      ! MDK, 23/10/26)
  l_ximpair_leaf_loss = .FALSE.,                                               &
      ! Xylem impairment model 3 (memory) with TRIFFID phenology (l_phenol):
      ! cap the phenological state at what the damaged xylem can supply,
      ! phen <= 1 - ximpair_leaf_sens * (1 - k_cap/kmax), so the canopy sheds
      ! leaves after lasting hydraulic damage and regrows as it recovers -
      ! see phenol.
  l_ximpair_seg_memory = .FALSE.,                                              &
      ! Memory model (pft_xylem_impairment_model = 3) with segments
      ! (l_som_plant_segments): a separate cap for the stem and the leaf
      ! segment, each damaged by the water potential at its own leaf-side end
      ! (psi_stem, psi_leaf) on its own vulnerability curve
      ! (conductance_b/c_seg). The leaf segment recovers with leaf turnover
      ! (ximpair_tau_leaf), the stem by the store repair (l_ximpair_rec_nsc).
      ! ximpair_psi_driver is then not used. .FALSE.: one cap on both segments
      ! (as before). (MDK, 23/10/26: ximpair_tau_stem retired.)
  l_ximpair_rec_nsc = .FALSE.,                                                 &
      ! Memory model (pft_xylem_impairment_model = 3) with the TRIFFID
      ! storage pool (l_triffid_nsc): the stem's lost conductance recovers
      ! only by regrowing xylem, paid for from the store (Trugman et al.
      ! 2018, Ecol. Lett. 21: 1552). At each TRIFFID step a fraction u of
      ! the carbon released from the store rebuilds the lost conducting
      ! sapwood, u = X/(X + L) (1 - m_L L / W) (their Eqn 5, no xylem
      ! turnover), with X the live stemwood and L the leaf and fine-root
      ! carbon the damaged xylem supports; see triffid_nsc_mod. The only
      ! recovery of the stem (the time-based terms were retired 23/10/26);
      ! the leaf segment keeps ximpair_tau_leaf. Use with l_ximpair_leaf_loss
      ! (ximpair_leaf_sens = 1) for their leaf area in proportion to the
      ! functional xylem.
      ! MDK, 17/10/26
  l_ximpair_rec_leaf_turn = .FALSE.,                                           &
      ! Per-segment memory (l_ximpair_seg_memory) with the TRIFFID store
      ! (l_triffid_nsc): the leaf segment recovers by leaf renewal, not by
      ! the fixed e-folding time ximpair_tau_leaf. New leaves arrive intact,
      ! so the leaf segment's loss falls by the fraction of the canopy that
      ! is new: leaf production at the turnover rate g_leaf_0, made only when
      ! the store's growth limits allow (turgor, cold; a spring and an
      ! autumn flush at FR-Pue) and normalised by the running annual mean
      ! of the growth factor so a year's production still matches the
      ! turnover, plus, with phenology, the regrowth after shedding. No extra
      ! carbon: the store already pays for the leaves.
      ! MDK, 18/10/26
  l_phen_psi = .FALSE.,                                                        &
      ! Drought leaf shedding from water potential (Xu et al. 2016, New
      ! Phytol, Notes S1 Eqn S7), with phenology (l_phenol): a leaf biomass
      ! factor f_bleaf (0-1) falls by phen_dfbleaf a day once the predawn
      ! water potential (daily mean root-zone psi) has been below the turgor
      ! loss point phen_psi_tlp on phen_ndays consecutive days, rises by
      ! phen_dfbleaf a day once it has been above 0.5 phen_psi_tlp on
      ! phen_ndays days, and is held otherwise; the canopy is capped at
      ! f_bleaf (phen_max), shed leaves go to litter, and regrowth stops
      ! while the TRIFFID store is empty.
      ! MDK, 18/10/26
  l_nsc_wmax_norm = .FALSE.,                                                   &
      ! TRIFFID store with the turgor / cold limits: W_max (the running mean
      ! NPP) is divided by the running annual mean of the growth factor
      ! (nsc_ft_mean), so release is concentrated in the periods when growth
      ! is allowed (spring flush drawn from the store) while the annual
      ! release still matches NPP. .FALSE.: W_max = running mean NPP.
      ! MDK, 18/10/26
  l_nsc_xylem_rate_norm = .FALSE.,                                             &
      ! Store repair (l_ximpair_rec_nsc), diffuse-porous stem
      ! (nsc_xylem_type 1): the cambial cap on the repair,
      ! nsc_xylem_rate_max f_t X per year, is divided by the running annual
      ! mean of the growth factor <f_t> (nsc_ft_mean), as the leaf renewal
      ! (l_ximpair_rec_leaf_turn) and l_nsc_wmax_norm are, so that
      ! nsc_xylem_rate_max_io is the fraction of the sapwood the cambium can
      ! rebuild per year, made while growth is allowed. .FALSE.: the cap is
      ! scaled by f_t only, so its annual total is nsc_xylem_rate_max <f_t>
      ! (about a quarter of nsc_xylem_rate_max at FR-Pue, <f_t> ~ 0.27).
      ! The cambium's window (the share of the step after
      ! nsc_grow_start_doy, with the leaf flush) still applies on top, so
      ! the realised annual cap is nsc_xylem_rate_max times the share of
      ! the year's growth factor inside that window. Needs the turgor or
      ! cold limits (otherwise f_t = 1 and this does nothing). Ring-porous
      ! stems (nsc_xylem_type 2) are not affected: there
      ! nsc_xylem_rate_max is the earlywood ring's share of the sapwood,
      ! built once a year, not a rate. As for l_nsc_wmax_norm, <f_t> is
      ! floored at 0.05 (a hard-coded guard), and it starts at 1, so the
      ! cap is low in the first year of a run without spin-up. Switching
      ! this on from a dump without nsc_ft_mean needs nsc_ft_mean in the
      ! initial conditions.
      !
      ! MDK, 25/10/26
  l_nsc_xylem_ring_renewal = .FALSE.,                                          &
      ! Carbon-paid ring renewal for PFTs with proportional replacement
      ! (nsc_xylem_replace = 2; needs TRIFFID, not equilibrium TRIFFID, and
      ! l_ximpair_rec_nsc): the stem recovers by its new xylem ring, which
      ! is the wood TRIFFID grows from the store's release, instead of by a
      ! repair paid from the store. Each TRIFFID step the gross wood
      ! production G = dC_W / dt + g_wood C_W (the local wood increment
      ! plus the turnover replaced; the carbon to spreading is not wood)
      ! renews the share rho = nsc_xylem_cond_ratio G dt / X of the stem's
      ! conductance (X the conducting sapwood), and the whole stem's loss
      ! is diluted, L -> L (1 - rho): new xylem replaces old at random,
      ! intact and embolised alike (one pool; the two-pool state is
      ! folded in and kept at 0), as leaf renewal dilutes the leaf. No cap,
      ! no damage term and no background floor; no carbon is drawn or sent
      ! to litter (TRIFFID already grows this wood), so with the switch on
      ! the carbon budget is that of no repair and only the hydraulics (and
      ! what they feed) change. The store paces the renewal through its
      ! release W. Rules 0 and 1 keep the store-paid repair; ring-porous
      ! stems can not take rule 2. .FALSE.: rule 2 is repaired from the
      ! store, sized as for rule 1.
      !
      ! MDK, 25/10/26
  l_nsc_hold_drought = .FALSE.,                                                &
      ! TRIFFID store: the running mean of NPP behind W_max is held only
      ! while the recent water status imposes a stem loss above the
      ! background embolism by more than nsc_hold_plc,
      !   (1 - ximpair_fmin_acc) - nsc_plc_background > nsc_hold_plc,
      ! with ximpair_fmin_acc the worst intact stem fraction, a memory that
      ! decays with e-folding nsc_repair_window (needs l_ximpair_rec_nsc).
      ! This holds through dry summers with a loss above background +
      ! nsc_hold_plc and, after a large peak, until the memory decays
      ! (about nsc_repair_window ln(peak / threshold); into the next spring
      ! at FR-Pue), but it clears every wet season. Never held in spin-up
      ! (the impairment is frozen there). .FALSE.: held while the locked
      ! stem loss ximpair_lock is above nsc_hold_plc, which a slowly
      ! diluting legacy loss (rule 2, ring renewal) can keep for decades,
      ! freezing W_max: ring renewal needs this switch (a warning is given
      ! otherwise).
      !
      ! MDK, 25/10/26
  l_som_fast = .FALSE.
      ! RETIRED input: .TRUE. stops the run; use som_ci_search = 2 (bounded)
      ! and som_psi_solver = 3 (lookup table). Still read so namelists that
      ! set it .FALSE. run. MDK, 24/10/26
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
    som_base_parm, som_n_sample, som_n_ci_golden_iter,                         &
    l_som_skip_search_wellwatered, som_hc_negligible_tol,                      &
    l_som_fast,                                                                &
    l_ximpair_leaf_loss, l_ximpair_seg_memory,                                 &
    l_ximpair_rec_nsc, l_triffid_nsc, l_ximpair_rec_leaf_turn, l_phen_psi,     &
    l_nsc_wmax_norm, l_trif_hydmort,                                           &
    l_nsc_xylem_rate_norm, l_nsc_xylem_ring_renewal, l_nsc_hold_drought,       &
    ximpair_stomatal_response,                                                 &
    l_som_supply_limit, l_som_root_supply, l_som_nsl, l_som_plant_segments,    &
    sl_kfail_frac,                                                             &
    l_som_vcmax_psi, l_root_mass_fixed, root_mass_min,                         &
    l_leaf_temp, leaf_width, l_leaf_temp_gc_eq,                                &
    l_leaf_temp_ci,                                                            &
    leaf_shelter, leaf_aero_model, l_leaf_coexp_lai, leaf_wind_ext,            &
    ml_leaf_temp_model, kplant_model,                                          &
    l_som_gain_gross,                                                          &
    l_som_cuticular_loss,                                                      &
    l_som_capacitance,                                                         &
    l_som_gravity,                                                             &
    l_som_rhizo_mfp,                                                           &
    som_leaf_resist_frac, som_gl_max, light_curvature_fvcb,                    &
    som_psi_aprox_method, som_profit_model, som_psi_solver, som_ci_search,     &
    discrimination_model,                                                      &
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
! leaf_flux_mod and som_profit_model from it. The old namelist form
! (leaf_flux_mod = 2 with som_profit_model) stops (retired 24/10/26).
!-----------------------------------------------------------------------------
IF ( som_profit_model /= imdi ) THEN
  ! MDK, 25/10/26
  errcode = -101 ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_profit_model is set from stomata_model; the namelist ' //  &
               'value is ignored.')
END IF
som_profit_model = profit_max_profit_model
SELECT CASE ( stomata_model )
CASE ( stomata_profit_max, stomata_sox_opt, stomata_cmax, stomata_cgain,       &
       stomata_supply_loss, stomata_cap, stomata_pm2 )
  leaf_flux_mod = leaf_flux_stom_opt
  IF ( stomata_model == stomata_pm2 ) THEN
    ! MDK, 22/10/26
    som_profit_model = pm2_profit_model
  ELSE IF ( stomata_model == stomata_cap ) THEN
    ! CAP is the nonstomatal limitation without a stomatal cost.
    som_profit_model = cap_profit_model
    l_som_nsl = .TRUE.
  ELSE IF ( stomata_model == stomata_profit_max ) THEN
    som_profit_model = profit_max_profit_model
  ELSE IF ( stomata_model == stomata_cmax ) THEN
    som_profit_model = cmax_profit_model
  ELSE IF ( stomata_model == stomata_cgain ) THEN
    som_profit_model = cgain_profit_model
  ELSE IF ( stomata_model == stomata_supply_loss ) THEN
    som_profit_model = supply_loss_model
  ELSE
    som_profit_model = sox_opt_profit_model
  END IF
CASE ( stomata_jacobs, stomata_medlyn, stomata_sox_analytical,                 &
       stomata_desica )
  IF ( leaf_flux_mod == leaf_flux_stom_opt ) THEN
    ! The old form (leaf_flux_mod = 2 with som_profit_model) is retired.
    ! MDK, 24/10/26
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 "leaf_flux_mod = 2 is retired: set stomata_model " //         &
                 "instead (som_profit_model 1 -> 4 profit max, 2 -> 5 " //     &
                 "SOX_opt, 4 -> 6 CMax, 5 -> 7 CGain, 6 -> 10 CAP, " //        &
                 "7 -> 11 ProfitMax2)")
  ELSE
    leaf_flux_mod = leaf_flux_fsmc
  END IF
CASE DEFAULT
  ! Reported below.
END SELECT

! som_gl_max is in m/s (0.02 m/s is about 0.8 mol m-2 s-1). A value above
! 0.1 m/s is most likely in mol m-2 s-1, which in m/s turns the cap off.
! MDK, 25/10/26
IF ( leaf_flux_mod == leaf_flux_stom_opt .AND. som_gl_max > 0.1 ) THEN
  errcode = -101 ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_gl_max is in m/s (default 0.02, ~0.8 mol m-2 s-1); ' //    &
               'the value set looks like mol m-2 s-1.')
END IF

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
CASE ( 1, 5, 6 )
  ! These are valid, so nothing to do.
CASE ( 4 )
  ! can_rad_mod = 4 only has the fsmc leaf path (sf_stom calls leaf() for
  ! every layer): with a stomatal optimisation model it would silently run
  ! fsmc instead.
  ! MDK, 18/10/26
  IF ( leaf_flux_mod == leaf_flux_stom_opt ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'can_rad_mod=4 has no stomatal optimisation: '             // &
                 'stomata_model=4 to 7, 9 or 10 needs can_rad_mod=1, 5, '   // &
                 '6 or 7')
  END IF
CASE ( 7 )
  ! Two-leaf (sunlit/shaded big-leaf) canopy: stomatal optimisation only.
  IF ( leaf_flux_mod /= leaf_flux_stom_opt ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'can_rad_mod=7 requires stomata_model=4 to 7 or 9')
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
! Only 1 (Ci) is implemented: 2 (psi) only prints a message in stom_opt_mod
! and leaves every point closed.
! MDK, 22/10/26
SELECT CASE( som_base_parm)
CASE ( som_base_parm_ci )
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_base_parm should be 1 (2, psi, is not implemented)')
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

! Retired inputs: l_som_fast and som_psi_aprox_method (no longer aliases).
! MDK, 24/10/26
IF ( som_psi_aprox_method /= imdi ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               "som_psi_aprox_method is retired: use som_psi_solver " //       &
               "(same values)")
END IF
IF ( l_som_fast ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               "l_som_fast is retired: use som_ci_search = 2 and " //          &
               "som_psi_solver = 3")
END IF

! Check that the som_psi_solver is suitable. JBaguley
SELECT CASE( som_psi_solver )
CASE ( psi_solver_newton, psi_solver_lut )
  ! Valid values
CASE ( 1 )
  ! RETIRED: the Taylor series solver. MDK, 25/10/26
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
     'som_psi_solver = 1 (Taylor) is retired: use 2 (Newton) or 3 ' //         &
     '(lookup table); utils/upgrade_nml_gce.py converts it')
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
     'som_psi_solver should be Newton-Raphson (2) or lookup table (3)')
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

! Check that the discrimination_model is suitable.
SELECT CASE( discrimination_model )
CASE ( discrimination_off, discrimination_simple, discrimination_photoresp )
  ! Valid values
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
     'discrimination_model should be off (0), simple (1) or ' //               &
     'photorespiration (2)')
END SELECT

! Check that the som_profit_model is suitable. JBaguley
SELECT CASE( som_profit_model)
CASE ( profit_max_profit_model, sox_opt_profit_model, cmax_profit_model,       &
       cgain_profit_model, cap_profit_model, pm2_profit_model )
  ! Valid values
CASE ( supply_loss_model )
  ! Valid only as set from stomata_model = 9 (above).
  IF ( stomata_model /= stomata_supply_loss ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'som_profit_model = 3 is set by stomata_model = 9')
  END IF
CASE DEFAULT
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'som_profit_model should be profit max (1), SOX (2), ' //       &
               'CMax (4), CGain (5), CAP (6) or ProfitMax2 (7)')
END SELECT

! CMax and CGain: absolute-carbon profit on the bounded Ci search (no
! flat-grid or fast-path variant: the cost does not vanish where the
! profit max's does).
IF ( ( stomata_model == stomata_cmax .OR. stomata_model == stomata_cgain       &
       .OR. stomata_model == stomata_cap .OR. stomata_model == stomata_pm2 )   &
     .AND. ( som_base_parm /= som_base_parm_ci .OR.                            &
             som_ci_search /= som_ci_bounded ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'stomata_model = 6 (CMax), 7 (CGain), 10 (CAP) or 11 ' //       &
               '(ProfitMax2) requires som_base_parm = 1 and ' //               &
               'som_ci_search = 2')
END IF
! Coupled transpiration in the optimisation: big leaf and two-leaf only (the
! multilayer canopy and the supply-loss rule keep the fixed-deficit E).
! It is on by default, so elsewhere it is turned off rather than an error.
IF ( l_som_coupled_e .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.          &
                             stomata_model == stomata_supply_loss .OR.         &
                             ( can_rad_mod /= 1 .AND. can_rad_mod /= 7 ) ) )   &
  THEN
  l_som_coupled_e = .FALSE.
  ! Say so only for a stomatal optimisation run (it is on by default).
  IF ( leaf_flux_mod == leaf_flux_stom_opt )                                   &
    CALL jules_print('check_jules_vegetation',                                 &
                     'Coupled transpiration (l_som_coupled_e) is off: it ' //  &
                     'applies only to stomata_model = 4 to 7, 10 or 11 ' //    &
                     'with can_rad_mod 1 or 7')
END IF
! Supply-loss stomata: the profit-max hydraulics, on the Ci axis.
IF ( stomata_model == stomata_supply_loss .AND.                                &
     som_base_parm /= som_base_parm_ci ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'stomata_model = 9 (supply-loss) requires som_base_parm = 1')
END IF
! The loss rule needs the supply function and its slope k = -dE/dpsi_leaf
! consistently.
IF ( stomata_model == stomata_supply_loss .AND.                                &
     som_psi_solver /= psi_solver_newton .AND.                                 &
     som_psi_solver /= psi_solver_lut ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'stomata_model = 9 (supply-loss) requires som_psi_solver = ' // &
               '2 or 3 (3, the lookup table, recommended)')
END IF
IF ( sl_kfail_frac < 0.0 .OR. sl_kfail_frac >= 1.0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'sl_kfail_frac must be in [0, 1)')
END IF

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
CASE ( stomata_jacobs, stomata_medlyn, stomata_sox_analytical, stomata_desica, &
       stomata_profit_max, stomata_sox_opt, stomata_cmax, stomata_cgain,       &
       stomata_supply_loss, stomata_cap, stomata_pm2 )
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
                 "stomata_model = 8 (DESICA) requires can_rad_mod = 1")
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

! Segments with the multilayer canopy (can_rad_mod 5, 6): each leaf's
! parallel path is root, stem and leaf in series (MDK, 23/10/26).
IF ( l_som_plant_segments .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.     &
                                  ( som_psi_solver /= psi_solver_newton .AND.  &
                                    som_psi_solver /= psi_solver_lut ) .OR.    &
                                  can_rad_mod == 4 ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_plant_segments requires stomata_model=4 to 7 or 9, ' //  &
               'som_psi_solver=2 or 3 and can_rad_mod=1, 5, 6 or 7')
END IF
! kplant_model 1 (root and stem per unit ground area) is not set up for the
! multilayer canopy, whose leaves' paths are per unit leaf area.
!
! MDK, 23/10/26
IF ( kplant_model == kplant_model_ground .AND.                                 &
     ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'kplant_model=1 is not available with can_rad_mod=5 or 6')
END IF

! kplant_model 1 needs the segments (which imply the profit-max path and
! not can_rad_mod 4, checked above; DESICA has no segments).
! MDK, 21/10/26
IF ( ( kplant_model /= kplant_model_param .AND.                                &
       kplant_model /= kplant_model_ground ) .OR.                              &
     ( kplant_model == kplant_model_ground .AND.                               &
       .NOT. l_som_plant_segments ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'kplant_model must be 0, or 1 with l_som_plant_segments=T')
END IF

! The plant water store and its target on the transpiration.
!
! MDK, 23/10/26
l_som_capacitance_on = ( l_som_capacitance .OR. l_som_cuticular_loss ) .AND.   &
                       leaf_flux_mod == leaf_flux_stom_opt
IF ( l_som_capacitance .AND. leaf_flux_mod /= leaf_flux_stom_opt ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_capacitance requires the stomatal optimisation')
END IF
IF ( l_som_capacitance_on .AND. can_rad_mod == 4 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'Plant capacitance (l_som_capacitance or ' //                   &
               'l_som_cuticular_loss) is not available with can_rad_mod = 4')
END IF
! The plant capacitance store sits at the stem node of the plant segments,
! seen by the optimiser (its flow out held fixed in each pass of the
! search). It is coded for profit max with the root, stem and leaf segments
! and the Newton solver, the bounded Ci search and the big-leaf, multilayer
! and two-leaf canopies (with or without the cuticular loss).
!
! MDK, 25/10/26
IF ( l_som_capacitance_on ) THEN
  IF ( .NOT. l_som_plant_segments .OR.                                         &
       som_psi_solver /= psi_solver_newton ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'Plant capacitance (l_som_capacitance or ' //                 &
                 'l_som_cuticular_loss) requires l_som_plant_segments ' //    &
                 'and som_psi_solver = 2 (Newton)')
  END IF
  IF ( can_rad_mod /= 1 .AND. can_rad_mod /= 5 .AND. can_rad_mod /= 6 .AND.   &
       can_rad_mod /= 7 ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'Plant capacitance (l_som_capacitance or ' //                 &
                 'l_som_cuticular_loss) is only coded for can_rad_mod = ' //   &
                 '1, 5, 6 or 7')
  END IF
  IF ( stomata_model /= stomata_profit_max ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'Plant capacitance (l_som_capacitance or ' //                 &
                 'l_som_cuticular_loss) is only coded for profit max ' //     &
                 '(stomata_model = 4)')
  END IF
  IF ( som_ci_search /= som_ci_bounded ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'Plant capacitance (l_som_capacitance or ' //                 &
                 'l_som_cuticular_loss) requires som_ci_search = 2 ' //        &
                 '(bounded)')
  END IF
  IF ( l_som_skip_search_wellwatered ) THEN
    errcode = 101
    CALL ereport("check_jules_vegetation", errcode,                            &
                 'Plant capacitance (l_som_capacitance or ' //                 &
                 'l_som_cuticular_loss) is not coded with ' //                 &
                 'l_som_skip_search_wellwatered')
  END IF
END IF

! The stomatal optimisation's leak also runs with the multilayer canopy
! (can_rad_mod 5, 6); DESICA's with big leaf and two-leaf only.
!
! MDK, 23/10/26
IF ( l_som_cuticular_loss .AND.                                                &
     ( ( leaf_flux_mod /= leaf_flux_stom_opt .AND.                             &
         stomata_model /= stomata_desica ) .OR.                                &
       ( leaf_flux_mod == leaf_flux_stom_opt .AND. can_rad_mod == 4 ) .OR.     &
       ( stomata_model == stomata_desica .AND.                                 &
         can_rad_mod /= 1 .AND. can_rad_mod /= 7 ) ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_cuticular_loss requires stomata_model=4 to 9, ' //       &
               'and can_rad_mod=1, 5, 6 or 7 (1 for DESICA)')
END IF

IF ( l_som_gravity .AND. leaf_flux_mod /= leaf_flux_stom_opt ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_gravity requires stomata_model=4 to 7 or 9 (DESICA ' //  &
               'always includes gravity)')
END IF

IF ( l_som_supply_limit .AND. ( leaf_flux_mod /= leaf_flux_stom_opt .OR.       &
                                .NOT. l_use_pft_psi ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_supply_limit requires stomata_model=4 to 7 or 9 and ' // &
               'l_use_pft_psi=T')
END IF

IF ( l_som_root_supply ) THEN
  errcode = -101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_root_supply is retired and ignored: the soil-to-root ' //&
               'resistance is now always in series in the profit-max path')
  l_som_root_supply = .FALSE.
END IF

! The soil-to-root conductance in series (see l_som_rhizo_series) is part of
! the plant hydraulics of the profit max and DESICA.
l_som_rhizo_series = ( leaf_flux_mod == leaf_flux_stom_opt .OR.                &
                       stomata_model == stomata_desica )
! l_som_rhizo_mfp is retired (see its declaration).
! MDK, 23/10/26
IF ( l_som_rhizo_mfp ) THEN
  errcode = -101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_rhizo_mfp is retired and ignored: the soil-to-root ' //  &
               'link is linear (code kept on branch rhizo_mfp)')
  l_som_rhizo_mfp = .FALSE.
END IF

IF ( .NOT. l_root_mass_fixed .AND. root_mass_min <= 0.0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'root_mass_min must be > 0 (or set l_root_mass_fixed)')
END IF

! (can_rad_mod 5 and 6: the multilayer canopy, see ml_leaf_temp_model.
! MDK, 19/10/26)
IF ( l_leaf_temp .AND. ( can_rad_mod < 5 .OR. can_rad_mod > 7 .OR.             &
                         leaf_flux_mod /= leaf_flux_stom_opt ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp is only coded for can_rad_mod=5, 6 or 7 with ' //  &
               'stomata_model=4 to 7 or 9')
END IF

IF ( l_leaf_temp .AND. ( leaf_shelter <= 0.0 .OR. leaf_aero_model < 0 .OR.     &
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

! MDK, 19/10/26
IF ( l_leaf_temp .AND. ( can_rad_mod == 5 .OR. can_rad_mod == 6 ) .AND.        &
     ml_leaf_temp_model /= ml_leaf_temp_aggregate .AND.                        &
     ml_leaf_temp_model /= ml_leaf_temp_layers ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'ml_leaf_temp_model must be 1 (aggregate) or 2 (per layer)')
END IF

IF ( l_leaf_temp .AND. ( leaf_wind_ext < 0.0 .OR.                              &
     ( leaf_aero_model == 1 .AND. leaf_wind_ext > 0.0 ) ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp needs leaf_wind_ext >= 0, and 0 with ' //          &
               'leaf_aero_model = 1 (which uses CABLE''s coexp)')
END IF

IF ( l_leaf_temp .AND. leaf_width <= 0.0 ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp needs leaf_width > 0')
END IF

! MDK, 24/10/26
IF ( l_leaf_temp_ci .AND. ( .NOT. l_leaf_temp .OR. can_rad_mod /= 7 .OR.       &
     leaf_flux_mod /= leaf_flux_stom_opt .OR.                                  &
     som_ci_search /= som_ci_bounded .OR. l_som_nsl ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_leaf_temp_ci needs l_leaf_temp, the two-leaf canopy ' //     &
               '(can_rad_mod = 7), stomatal optimisation with the ' //         &
               'bounded Ci search (som_ci_search = 2), and not l_som_nsl')
END IF

IF ( ximpair_stomatal_response /= ximpair_hold_gs .AND.                        &
     ximpair_stomatal_response /= ximpair_hold_psi ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'ximpair_stomatal_response should be 0 (hold gs) or ' //        &
               '1 (hold leaf psi)')
END IF
IF ( ximpair_stomatal_response == ximpair_hold_psi .AND.                       &
     som_ci_search /= som_ci_bounded ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'ximpair_stomatal_response = 1 is coded for som_ci_search ' //  &
               '= 2 only')
END IF

IF ( l_som_vcmax_psi .AND. ( ( stomata_model /= stomata_profit_max .AND.       &
                               stomata_model /= stomata_sox_opt .AND.          &
                               stomata_model /= stomata_cmax .AND.             &
                               stomata_model /= stomata_cgain .AND.            &
                               stomata_model /= stomata_supply_loss .AND.      &
                               stomata_model /= stomata_cap .AND.              &
                               stomata_model /= stomata_pm2 ) .OR.             &
                             photo_model /= photo_farquhar ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_vcmax_psi requires stomata_model=4 to 7 or 9 to 11 ' //  &
               'and ' //                                                       &
               'photo_model=2')
END IF

IF ( l_som_nsl .AND. ( ( stomata_model /= stomata_profit_max .AND.             &
                         stomata_model /= stomata_cap ) .OR.                   &
                       som_ci_search /= som_ci_bounded ) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_nsl requires stomata_model=4 (or 10) and ' //            &
               'som_ci_search=2')
END IF

! Both water limits on photosynthesis: A = f_vc(psi_soil) f_nsl(psi_leaf)
! A0. Each is fitted on its own to the same kind of data, so together they
! count the drought limit about twice. MDK, 25/10/26
IF ( l_som_nsl .AND. l_som_vcmax_psi ) THEN
  errcode = -101 ! warning
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_som_nsl and l_som_vcmax_psi are both on: photosynthesis ' // &
               'is limited by soil and leaf water potential together; ' //     &
               'each was fitted alone.')
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

! TRIFFID storage pool and the xylem repair it pays for.
! MDK, 17/10/26
IF ( l_triffid_nsc .AND. l_triffid .AND. l_trif_eq ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_triffid_nsc is not coded for equilibrium TRIFFID ' //        &
               '(l_trif_eq=T)')
END IF
IF ( l_triffid_nsc .AND. l_sugar ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_triffid_nsc and l_sugar are two different NSC stores; ' //   &
               'use one')
END IF
! MDK, 18/10/26
IF ( l_ximpair_rec_leaf_turn .AND.                                             &
     .NOT. (l_triffid_nsc .AND. l_ximpair_seg_memory) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_ximpair_rec_leaf_turn needs the TRIFFID store and the '  //  &
               'per-segment memory (l_triffid_nsc, l_ximpair_seg_memory)')
END IF
! MDK, 21/10/26
IF ( l_trif_hydmort .AND. .NOT. l_triffid ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_trif_hydmort needs TRIFFID (l_triffid=T)')
END IF
IF ( l_trif_hydmort .AND. l_triffid .AND. l_trif_eq ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_trif_hydmort is not coded for equilibrium TRIFFID ' //       &
               '(l_trif_eq=T)')
END IF
IF ( l_nsc_wmax_norm .AND. .NOT. l_triffid_nsc ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_nsc_wmax_norm needs the TRIFFID store (l_triffid_nsc=T)')
END IF
! The stem recovery switches act on the store repair.
!
! MDK, 25/10/26
IF ( (l_nsc_xylem_rate_norm .OR. l_nsc_xylem_ring_renewal .OR.                &
      l_nsc_hold_drought) .AND. .NOT. l_ximpair_rec_nsc ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_nsc_xylem_rate_norm, l_nsc_xylem_ring_renewal and ' //       &
               'l_nsc_hold_drought need the store repair '              //     &
               '(l_ximpair_rec_nsc=T)')
END IF
IF ( l_nsc_xylem_ring_renewal .AND. .NOT. l_nsc_hold_drought ) THEN
  errcode = -101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_nsc_xylem_ring_renewal without l_nsc_hold_drought: the ' //  &
               'slowly diluting stem loss will hold W_max for years')
END IF
IF ( l_nsc_xylem_ring_renewal .AND. (.NOT. l_triffid .OR. l_trif_eq) ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_nsc_xylem_ring_renewal needs dynamic TRIFFID '        //     &
               '(l_triffid=T, l_trif_eq=F): the new ring is its wood growth')
END IF
IF ( l_phen_psi .AND. .NOT. l_phenol ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_phen_psi (drought leaf shedding) needs l_phenol=T')
END IF
IF ( l_ximpair_rec_nsc .AND. .NOT. l_triffid_nsc ) THEN
  errcode = 101
  CALL ereport("check_jules_vegetation", errcode,                              &
               'l_ximpair_rec_nsc needs the TRIFFID store (l_triffid_nsc=T)')
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

WRITE(lineBuffer,*) ' discrimination_model = ', discrimination_model
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_supply_limit = ', l_som_supply_limit
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_root_supply = ', l_som_root_supply
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_som_rhizo_series = ', l_som_rhizo_series
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_root_mass_fixed = ', l_root_mass_fixed
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' root_mass_min = ', root_mass_min
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_nsl = ', l_som_nsl
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' sl_kfail_frac = ', sl_kfail_frac
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
  WRITE(lineBuffer,*) ' l_leaf_temp_ci = ', l_leaf_temp_ci
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_aero_model = ', leaf_aero_model
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_shelter = ', leaf_shelter
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' l_leaf_coexp_lai = ', l_leaf_coexp_lai
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' leaf_wind_ext = ', leaf_wind_ext
  CALL jules_print('jules_vegetation_mod',lineBuffer)
  WRITE(lineBuffer,*) ' ml_leaf_temp_model = ', ml_leaf_temp_model
  CALL jules_print('jules_vegetation_mod',lineBuffer)
END IF

WRITE(lineBuffer,*) ' l_som_plant_segments = ', l_som_plant_segments
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' kplant_model = ', kplant_model   ! MDK, 21/10/26
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_gain_gross = ', l_som_gain_gross
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_cuticular_loss = ', l_som_cuticular_loss
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_capacitance = ', l_som_capacitance   ! MDK, 23/10/26
CALL jules_print('jules_vegetation_mod',lineBuffer)



WRITE(lineBuffer,*) ' l_ximpair_leaf_loss = ', l_ximpair_leaf_loss
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_ximpair_seg_memory = ', l_ximpair_seg_memory
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_ximpair_rec_nsc = ', l_ximpair_rec_nsc
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_triffid_nsc = ', l_triffid_nsc
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_trif_hydmort = ', l_trif_hydmort
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_ximpair_rec_leaf_turn = ', l_ximpair_rec_leaf_turn
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_phen_psi = ', l_phen_psi
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_nsc_wmax_norm = ', l_nsc_wmax_norm
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_nsc_xylem_rate_norm = ', l_nsc_xylem_rate_norm
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_nsc_xylem_ring_renewal = ',                            &
                    l_nsc_xylem_ring_renewal
CALL jules_print('jules_vegetation_mod',lineBuffer)
WRITE(lineBuffer,*) ' l_nsc_hold_drought = ', l_nsc_hold_drought
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' ximpair_stomatal_response = ', ximpair_stomatal_response
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_gravity = ', l_som_gravity
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_coupled_e = ', l_som_coupled_e
CALL jules_print('jules_vegetation_mod',lineBuffer)

WRITE(lineBuffer,*) ' l_som_rhizo_mfp = ', l_som_rhizo_mfp
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
INTEGER, PARAMETER :: n_int = 24
                                  ! -1 leaf_temp_iter (MDK, 25/10/26),
                                  ! -1 som_store_coupling removed (MDK, 25/10/26),
                                  ! +1 discrimination_model (MDK, 22/10/26), +1
                                  ! kplant_model (MDK, 21/10/26), +1
                                  ! ml_leaf_temp_model (MDK, 19/10/26), +1
                                  ! ximpair_stomatal_response, +2
                                  ! leaf_temp_iter/leaf_aero_model, was 16, +1
                                  ! for som_n_ci_golden_iter,
                                 ! +2 for som_psi_solver/som_ci_search
INTEGER, PARAMETER :: n_real = 20 + (n_photo_coef * 5)
                                  ! +1 leaf_wind_ext (MDK, 19/10/26), +1
                                  ! sl_kfail_frac, +1 root_mass_min, +2
                                  ! leaf_width/shelter, +4 for
                                  ! som_hc_negligible_tol/som_leaf_resist_frac/
                                  ! som_gl_max/light_curvature_fvcb
INTEGER, PARAMETER :: n_log = 57 + npft_max
                                  ! +3 l_nsc_xylem_rate_norm,
                                  ! l_nsc_xylem_ring_renewal,
                                  ! l_nsc_hold_drought (MDK, 25/10/26)
                                  ! -1 l_som_coupled_e (MDK, 25/10/26)
                                  ! -3 l_som_plant_store,
                                  ! l_som_cuticular_floor and
                                  ! l_som_leak_capacitance_only removed (MDK,
                                  ! 25/10/26),
                                  ! -1 l_som_rhizo_hr removed (MDK,
                                  ! 24/10/26),
                                  ! -1 l_ximpair_rec_growth retired (MDK,
                                  ! 23/10/26),
                                  ! +3 l_som_capacitance/plant_store
                                  ! (deprecated)/leak_capacitance_only (MDK
                                  ! 23/10/26), +1 l_trif_hydmort, +1
                                  ! l_som_cuticular_floor (deprecated name, MDK
                                  ! 21/10/26), +1 l_nsc_wmax_norm, +1
                                  ! l_phen_psi, +1 l_ximpair_rec_leaf_turn, +2
                                  ! l_ximpair_rec_nsc/l_triffid_nsc, +1
                                  ! l_ximpair_seg_memory, +2 l_som_rhizo_mfp/hr,
                                  ! +1 l_som_coupled_e, +1 l_root_mass_fixed, +1
                                  ! l_som_vcmax_psi, +1 l_leaf_coexp_lai,
                                  ! +2 l_leaf_temp(_gc_eq), +1 for l_som_fast,
                                  ! +1 for
                                  ! l_som_gain_gross, +1 for
                                  ! l_som_cuticular_loss, +1 for
                                  ! l_som_gravity, +1 for
                                  ! l_som_skip_search_wellwatered, +1 for
                                  ! l_som_supply_limit, +1 for
                                  ! l_som_root_supply, +1 for l_som_nsl,
                                  ! +1 for
                                  ! l_som_plant_segments, +2 for
                                  ! l_ximpair_rec_growth/
                                  ! l_ximpair_leaf_loss (trunk vn7.9: 29)

TYPE :: my_namelist
  SEQUENCE
  INTEGER :: phenol_period
  INTEGER :: triffid_period
  INTEGER :: can_model
  INTEGER :: can_rad_mod
  INTEGER :: leaf_aero_model
  INTEGER :: ml_leaf_temp_model
  INTEGER :: kplant_model
  INTEGER :: ilayers
  INTEGER :: leaf_flux_mod !JBaguley
  INTEGER :: som_base_parm !JBaguley
  INTEGER :: som_n_sample !JBaguley
  INTEGER :: som_n_ci_golden_iter
  INTEGER :: som_psi_aprox_method !JBaguley
  INTEGER :: som_psi_solver
  INTEGER :: som_ci_search
  INTEGER :: ximpair_stomatal_response
  INTEGER :: discrimination_model
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
  REAL(KIND=real_jlslsm) :: sl_kfail_frac
  REAL(KIND=real_jlslsm) :: light_curvature_fvcb
  REAL(KIND=real_jlslsm) :: leaf_width
  REAL(KIND=real_jlslsm) :: leaf_shelter
  REAL(KIND=real_jlslsm) :: leaf_wind_ext
  REAL(KIND=real_jlslsm) :: root_mass_min
  LOGICAL :: l_som_skip_search_wellwatered
  LOGICAL :: l_som_fast
  LOGICAL :: l_som_supply_limit
  LOGICAL :: l_som_root_supply
  LOGICAL :: l_root_mass_fixed
  LOGICAL :: l_som_nsl
  LOGICAL :: l_som_vcmax_psi
  LOGICAL :: l_leaf_temp
  LOGICAL :: l_leaf_temp_gc_eq
  LOGICAL :: l_leaf_temp_ci
  LOGICAL :: l_leaf_coexp_lai
  LOGICAL :: l_som_plant_segments
  LOGICAL :: l_som_gain_gross
  LOGICAL :: l_som_cuticular_loss
  LOGICAL :: l_som_capacitance   ! MDK, 23/10/26
  LOGICAL :: l_ximpair_leaf_loss
  LOGICAL :: l_ximpair_seg_memory
  LOGICAL :: l_ximpair_rec_nsc
  LOGICAL :: l_triffid_nsc
  LOGICAL :: l_trif_hydmort
  LOGICAL :: l_ximpair_rec_leaf_turn
  LOGICAL :: l_phen_psi
  LOGICAL :: l_nsc_wmax_norm
  LOGICAL :: l_nsc_xylem_rate_norm
  LOGICAL :: l_nsc_xylem_ring_renewal
  LOGICAL :: l_nsc_hold_drought
  LOGICAL :: l_som_gravity
  LOGICAL :: l_som_rhizo_mfp
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
  my_nml % ximpair_stomatal_response = ximpair_stomatal_response
  my_nml % discrimination_model = discrimination_model
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
  my_nml % sl_kfail_frac = sl_kfail_frac
  my_nml % light_curvature_fvcb = light_curvature_fvcb
  my_nml % l_som_skip_search_wellwatered = l_som_skip_search_wellwatered
  my_nml % l_som_fast = l_som_fast
  my_nml % l_som_supply_limit = l_som_supply_limit
  my_nml % l_som_root_supply = l_som_root_supply
  my_nml % l_root_mass_fixed = l_root_mass_fixed
  my_nml % root_mass_min = root_mass_min
  my_nml % l_som_nsl = l_som_nsl
  my_nml % l_som_vcmax_psi = l_som_vcmax_psi
  my_nml % l_leaf_temp = l_leaf_temp
  my_nml % l_leaf_temp_gc_eq = l_leaf_temp_gc_eq
  my_nml % l_leaf_temp_ci = l_leaf_temp_ci
  my_nml % l_leaf_coexp_lai = l_leaf_coexp_lai
  my_nml % leaf_aero_model = leaf_aero_model
  my_nml % ml_leaf_temp_model = ml_leaf_temp_model
  my_nml % kplant_model = kplant_model
  my_nml % leaf_shelter = leaf_shelter
  my_nml % leaf_width = leaf_width
  my_nml % leaf_wind_ext = leaf_wind_ext
  my_nml % l_som_plant_segments = l_som_plant_segments
  my_nml % l_som_gain_gross = l_som_gain_gross
  my_nml % l_som_cuticular_loss = l_som_cuticular_loss
  my_nml % l_som_capacitance = l_som_capacitance
  my_nml % l_ximpair_leaf_loss = l_ximpair_leaf_loss
  my_nml % l_ximpair_seg_memory = l_ximpair_seg_memory
  my_nml % l_ximpair_rec_nsc = l_ximpair_rec_nsc
  my_nml % l_triffid_nsc = l_triffid_nsc
  my_nml % l_trif_hydmort = l_trif_hydmort
  my_nml % l_ximpair_rec_leaf_turn = l_ximpair_rec_leaf_turn
  my_nml % l_phen_psi = l_phen_psi
  my_nml % l_nsc_wmax_norm = l_nsc_wmax_norm
  my_nml % l_nsc_xylem_rate_norm = l_nsc_xylem_rate_norm
  my_nml % l_nsc_xylem_ring_renewal = l_nsc_xylem_ring_renewal
  my_nml % l_nsc_hold_drought = l_nsc_hold_drought
  my_nml % l_som_gravity = l_som_gravity
  my_nml % l_som_rhizo_mfp = l_som_rhizo_mfp
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
  ximpair_stomatal_response = my_nml % ximpair_stomatal_response
  discrimination_model = my_nml % discrimination_model
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
  sl_kfail_frac = my_nml % sl_kfail_frac
  light_curvature_fvcb = my_nml % light_curvature_fvcb
  l_som_skip_search_wellwatered = my_nml % l_som_skip_search_wellwatered
  l_som_fast = my_nml % l_som_fast
  l_som_supply_limit = my_nml % l_som_supply_limit
  l_som_root_supply = my_nml % l_som_root_supply
  l_root_mass_fixed = my_nml % l_root_mass_fixed
  root_mass_min = my_nml % root_mass_min
  l_som_nsl = my_nml % l_som_nsl
  l_som_vcmax_psi = my_nml % l_som_vcmax_psi
  l_leaf_temp = my_nml % l_leaf_temp
  l_leaf_temp_gc_eq = my_nml % l_leaf_temp_gc_eq
  l_leaf_temp_ci = my_nml % l_leaf_temp_ci
  l_leaf_coexp_lai = my_nml % l_leaf_coexp_lai
  leaf_aero_model = my_nml % leaf_aero_model
  ml_leaf_temp_model = my_nml % ml_leaf_temp_model
  kplant_model = my_nml % kplant_model
  leaf_shelter = my_nml % leaf_shelter
  leaf_width = my_nml % leaf_width
  leaf_wind_ext = my_nml % leaf_wind_ext
  l_som_plant_segments = my_nml % l_som_plant_segments
  l_som_gain_gross = my_nml % l_som_gain_gross
  l_som_cuticular_loss = my_nml % l_som_cuticular_loss
  l_som_capacitance = my_nml % l_som_capacitance
  l_ximpair_leaf_loss = my_nml % l_ximpair_leaf_loss
  l_ximpair_seg_memory = my_nml % l_ximpair_seg_memory
  l_ximpair_rec_nsc = my_nml % l_ximpair_rec_nsc
  l_triffid_nsc = my_nml % l_triffid_nsc
  l_trif_hydmort = my_nml % l_trif_hydmort
  l_ximpair_rec_leaf_turn = my_nml % l_ximpair_rec_leaf_turn
  l_phen_psi = my_nml % l_phen_psi
  l_nsc_wmax_norm = my_nml % l_nsc_wmax_norm
  l_nsc_xylem_rate_norm = my_nml % l_nsc_xylem_rate_norm
  l_nsc_xylem_ring_renewal = my_nml % l_nsc_xylem_ring_renewal
  l_nsc_hold_drought = my_nml % l_nsc_hold_drought
  l_som_gravity = my_nml % l_som_gravity
  l_som_rhizo_mfp = my_nml % l_som_rhizo_mfp
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
