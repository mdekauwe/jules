! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE stom_opt_jls_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STOM_OPT_JLS_MOD'

PRIVATE stom_opt_mod_ci, stom_opt_profit_max_select, stom_opt_bounded_search, &
        stom_opt_supply_loss, g_eff_e
PUBLIC stom_opt_mod, g_of_geff, som_path_psi, stom_opt_at_e

! Diagnostic marginal water price at the optimum, lambda_eff = dA/dE
! (mol CO2 per mol H2O; = umol per mmol x 1e-3), from the bounded Ci
! search (som_ci_search = 2: profit max, CMax, CGain, CAP, ProfitMax2), by
! a finite difference about the chosen Ci on the search path; 0 for closed
! stomata and for the other searches. som_lambda_call: of the latest
! stom_opt_mod call, by land point; som_lambda_eff_pft: the canopy value
! per PFT, leaf-area weighted over the sunlit and shaded leaves by sf_stom
! (can_rad_mod = 7 only), for the output lambda_eff_pft.
!
! MDK, 22/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: som_lambda_call(:)
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: som_lambda_eff_pft(:,:)

CONTAINS

!-----------------------------------------------------------------------------
! Coupled transpiration (l_som_coupled_e): the effective conductance for
! water vapour from the leaf to the level-1 air, so that E = D g_eff, with D
! the deficit from the leaf to the level-1 air.
!
! Water vapour crosses the stomata and the leaf boundary layer in series,
! then the canopy air, which this leaf class shares with the other class:
!   g'    = g / (1 + g r_bl)
!   g_eff = g' / (1 + r_ca (g' + g_o))
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION g_eff_e( g, r_bl, r_ca, g_o ) RESULT( g_eff )
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  g,                                                                           &
                            ! Stomatal conductance for H2O (m/s, >= 0).
  r_bl,                                                                        &
                            ! Leaf boundary layer resistance (s/m).
  r_ca,                                                                        &
                            ! Canopy air to level-1 resistance (s/m).
  g_o
                            ! Effective conductance of the other leaf class
                            ! sharing the canopy air (m/s).
REAL(KIND=real_jlslsm) :: g_eff
                            ! Effective conductance (m/s).
REAL(KIND=real_jlslsm) :: gp
                            ! Stomata and boundary layer in series (m/s).
gp    = g / ( 1.0 + g * r_bl )
g_eff = gp / ( 1.0 + r_ca * ( gp + g_o ) )
END FUNCTION g_eff_e

!-----------------------------------------------------------------------------
! The inverse of g_eff_e: the stomatal conductance g whose effective
! conductance is x. Returns HUGE when no finite g reaches x, because the
! resistances in series cap g_eff below it.
!
! Used to turn a cap on transpiration into a cap on g: the hydraulic limit
! in the bounded search, and the soil supply limit in sf_stom.
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION g_of_geff( x, r_bl, r_ca, g_o ) RESULT( g )
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  x,                                                                           &
                            ! Effective conductance to reach (m/s).
  r_bl,                                                                        &
                            ! Leaf boundary layer resistance (s/m).
  r_ca,                                                                        &
                            ! Canopy air to level-1 resistance (s/m).
  g_o
                            ! Effective conductance of the other leaf class
                            ! sharing the canopy air (m/s).
REAL(KIND=real_jlslsm) :: g
                            ! Stomatal conductance for H2O (m/s).
REAL(KIND=real_jlslsm) :: gp, den
                            ! g' as in g_eff_e (m/s), and a denominator.
g   = HUGE(1.0_real_jlslsm)
den = 1.0 - r_ca * x
IF ( den <= 0.0 ) RETURN
gp  = x * ( 1.0 + r_ca * g_o ) / den
den = 1.0 - gp * r_bl
IF ( den <= 0.0 ) RETURN
g   = gp / den
END FUNCTION g_of_geff

! *********************************************************************
! Contains routines to calculate the optimal stomatal conductance and
! net carbon uptake.
!
! Profit maximisation in six steps (stomata_model = 4; Sperry et al. 2017).
! For each open leaf (or sunlit/shaded class, or canopy layer) the stomata
! choose the internal CO2, Ci, that maximises carbon gain minus hydraulic
! cost:
!
!   1. A(Ci)        net photosynthesis at Ci (Collatz or Farquhar)
!   2. gl(Ci)       = ratio A R T / (ca - Ci), the leaf conductance that
!                     supplies that CO2 (ratio = 1.6, water vs CO2)
!   3. E(Ci)        = gl x the leaf-to-air humidity deficit, the
!                     transpiration that follows
!   4. psi_leaf(E)  from the plant's supply function: the integral of its
!                     vulnerability curve from the soil to the leaf
!   5. k(E)         = -dE/dpsi_leaf, the plant conductance left at that E
!   6. profit       = A / A_max - (k_max - k) / (k_max - k_crit)
!
! States with k <= k_crit (hydraulic failure) or gl > gl_max are not
! allowed. Ci rather than gl is the search variable because A and gl are
! both explicit in Ci. The other stomatal models (CMax, CGain, SOX, CAP,
! ProfitMax2) keep steps 1-5 and change step 6 (see profit in
! stom_opt_bounded_search).
!
! MDK, 23/10/26
!
! References:
! - Sperry, J. S., Venturas, M. D., Anderegg, W. R. L., Mencuccini, M., Mackay,
!   D. S., Wang, Y., and Love, D. M. (2017) Predicting stomatal responses to
!   the environment from the optimization of photosynthetic gain and hydraulic
!   cost. Plant, Cell & Environment, 40: 816-830. doi: 10.1111/pce.12852.
! - De Kauwe, M.G., Sabot, M.E.B., Medlyn, B.E., Pitman, A.J., Meir, P.,
!   Cernusak, L.A., Gallagher, R.V., Ukkola, A.M., Rifai, S.W. and Choat, B.
!   (2022), Towards species-level forecasts of drought-induced tree mortality
!   risk. New Phytol, 235: 94-110. https://doi.org/10.1111/nph.18129
! - Baguley, J. C., Quaife, T., Rowland, L., Meir, P., Bell, M. C., Sabot, M.
!   E. B., and De Kauwe, M. G.: Limited modelled legacy effects of drought on
!   carbon and water fluxes in a temperate oak forest despite hydraulic
!   impairment. Biogeosciences, in review
! *********************************************************************

SUBROUTINE stom_opt_mod (                                                      &
! IN
        land_pts, som_base_parm, pft, open_pts, open_index,                    &
        pft_photo_model, veg_pts, veg_index,                                   &
        ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
        km, dq, qs, je, t_leaf, je_ratio, fapar_lf, kmax, kcrit,               &
        gl_max, ipar, l_multilayer,                                            &
        r_bl, r_ca, g_oth,                                                     &
        kmax_ref, conductance_b, conductance_c,                                &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment_in,             &
! IN OUT
        rd,                                                                    &
! OUT
        ci, al, el, flux_o3, fo3, gl, psi_leaf,                                &
        carbon_gain_out, hydraulic_cost_out, leaf_k,                           &
! OUT (optional)
        psi_stem,                                                              &
! IN (optional)
        lai_cost                                                               &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE jules_vegetation_mod, ONLY:                                                &
        som_base_parm_ci, som_base_parm_psi, som_n_sample,                     &
        profit_max_profit_model, sox_opt_profit_model, som_profit_model,       &
        cmax_profit_model, cgain_profit_model, supply_loss_model,              &
        cap_profit_model, som_ci_search, som_ci_bounded,                       &
        pm2_profit_model,                                                      &
        som_n_ci_golden_iter,                                                  &
        l_som_skip_search_wellwatered, som_hc_negligible_tol, l_som_nsl

USE pftparm, ONLY:                                                             &
        min_gl_pft

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_conductance_impaired_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: k_path_zero_flow, ksr_path,          &
                                       seg_kmax_fac   ! MDK, 21/10/26

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft,                       &
                   pft_xylem_impairment_model
USE jules_vegetation_mod, ONLY: xylem_impairment_none, l_som_plant_segments,  &
                                l_som_rhizo_series
USE pftparm, ONLY: conductance_b_seg, conductance_c_seg
USE model_time_mod, ONLY: is_spinup

LOGICAL, INTENT(IN) :: l_multilayer

LOGICAL, INTENT(IN) :: l_xylem_impairment_in
                            ! .TRUE. to apply the xylem impairment model
                            ! (pft_xylem_impairment_model) to the
                            ! vulnerability curve: leaf water potential,
                            ! xylem conductance and the kl > kcrit
                            ! feasibility test then use the impaired curve
                            ! (kmax, conductance_b, conductance_c and the
                            ! historic psi extremes), while the hydraulic
                            ! cost is measured against the unimpaired PFT
                            ! curve (kmax_ref, conductance_b_pft,
                            ! conductance_c_pft). JBaguley
                            ! .FALSE. uses kmax/conductance_b/conductance_c
                            ! directly, with no impairment - kmax_ref and
                            ! the psi extremes are then not used.

!-----------------------------------------------------------------------------
! IN integer variables.
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
!                           ! Number of land points to be processed
, som_base_parm                                                                &
!                           ! Switch used to determine the physical property
!                           ! that is uniformly iterated over when
!                           ! determining the optimal stomatal conductance.
!                           !     1 => leaf intercellular carbon (Ci)
!                           !     2 => leaf water potential (psi)
, pft                                                                          &
!                           ! Plant functional type
,open_pts                                                                      &
                            ! Number of land points with open stomata.
,open_index(land_pts)                                                          &
                            ! Index of land points with open stomata.
,pft_photo_model                                                               &
                            ! Indicates which photosynthesis model to use for
                            ! the current PFT.
,veg_pts                                                                       &
                            ! Number of vegetated points.
,veg_index(land_pts)
                            ! Index of vegetated points on the land grid.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_bl(land_pts),                                                              &
                            ! Leaf boundary layer resistance (s/m, on the
                            ! basis of gl), for l_som_coupled_e.
  r_ca(land_pts),                                                              &
                            ! Canopy air to level-1 resistance (s/m), for
                            ! l_som_coupled_e.
  g_oth(land_pts)
                            ! Effective conductance of the other leaf class
                            ! sharing the canopy air (m/s), for
                            ! l_som_coupled_e.

!-----------------------------------------------------------------------------
! IN real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 ca(land_pts)                                                                  &
                             ! Atmospheric CO2 pressure (Pa)
,psi_root_zone(land_pts)                                                       &
                             ! Root zone water potential (Pa)
,acr(land_pts)                                                                 &
                             ! Absorbed PAR (mol photons/m2/s).
,apar(land_pts)                                                                &
                             ! Absorbed PAR (W m-2).
,oi(land_pts)                                                                  &
                             ! Intercellular oxygen pressure (Pa).
,vcmax(land_pts)                                                               &
                             ! Maximum rate of carboxylation of Rubisco
                             !  (mol CO2/m2/s).
,kc(land_pts)                                                                  &
                             ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_pts)                                                                  &
                             ! Michaelis-Menten constant for O (Pa).
,ccp(land_pts)                                                                 &
                            ! Photorespiratory compensatory point (Pa).
,pstar(land_pts)                                                               &
                            ! Atmospheric pressure (Pa).
,km(land_pts)                                                                  &
                            ! A combination of Michaelis-Menten and other
                            ! terms. Only used with the Farquhar model.
,dq(land_pts)                                                                  &
                            ! Specific humidity deficit of air (kg H2O/kg air).
,qs(land_pts)                                                                  &
                            ! Saturated specific humidity
                            !      (kg H2O/kg air).
,je(land_pts)                                                                  &
                            ! Electron transport rate (mol m-2 s-1)
,t_leaf(land_pts)                                                              &
                            ! Leaf temperature (K)
,je_ratio(land_pts)                                                            &
                            ! Electron transport of the sunlit (or shaded)
                            ! leaves of this call over that of the canopy
                            ! (dimensionless). Scales J in multilayer calls;
                            ! 1 otherwise.
                            !
                            ! MDK, 21/10/26
,fapar_lf(land_pts)                                                            &
                            ! Fraction of the incident PAR absorbed by this
                            ! layer's sunlit (or shaded) leaves
                            ! (dimensionless). Used in multilayer calls.
                            !
                            ! MDK, 21/10/26
,kmax(land_pts)                                                                &
                            ! Maximum whole-path (soil-to-leaf) conductance
                            ! (mol m-2 s-1 Pa-1), on the basis of this call
                            ! (whole canopy, the sunlit or shaded leaves, or
                            ! a layer's leaves). The impaired value with
                            ! l_xylem_impairment.
                            !
                            ! MDK, 21/10/26
,kcrit(land_pts)                                                               &
                            ! Critical conductance (mol m-2 s-1 Pa-1, same
                            ! basis as kmax): states with k <= kcrit are
                            ! infeasible.
                            !
                            ! MDK, 21/10/26
,gl_max(land_pts)                                                              &
                            ! Maximum stomatal conductance for H2O (m/s, same
                            ! basis as gl: per leaf area for multilayer,
                            ! canopy for big-leaf). Samples with gl above it
                            ! are infeasible. <= 0 disables the cap.
,ipar(land_pts)                                                                &
                            ! Incident PAR (W m-2).
                            !
                            ! MDK, 21/10/26
,kmax_ref(land_pts)                                                            &
                            ! Unimpaired maximum conductance
                            ! (mol m-2 s-1 Pa-1), with the same scaling as
                            ! kmax (kmax is the impaired value when
                            ! l_xylem_impairment).
                            !
                            ! MDK, 21/10/26
,conductance_b(land_pts)                                                       &
                            ! Conductance parameter b for each land point
                            ! (Pa).
,conductance_c(land_pts)                                                       &
                            ! Conductance parameter c for each land point
                            ! (dimensionless).
                            !
                            ! MDK, 21/10/26
,psi_leaf_extreme(land_pts)                                                    &
                            ! Historic minimum leaf water potential (Pa).
,psi_root_extreme(land_pts)
                            ! Historic minimum root zone water potential (Pa).

!-----------------------------------------------------------------------------
! IN OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 rd(land_pts)
                            ! Dark respiration (mol CO2/m2/s).
                            ! Not modified: the ozone damage is not coded
                            ! here (flux_o3 and fo3 are returned as 0).
                            !
                            ! MDK, 21/10/26

!-----------------------------------------------------------------------------
! OUT real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ci(land_pts)                                                                  &
                            ! Internal CO2 pressure (Pa).
,al(land_pts)                                                                  &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,el(land_pts)                                                                  &
                            ! Transpiration rate (mol H2O/m2/s)
,flux_o3(land_pts)                                                             &
                            ! Flux of O3 to stomata (nmol O3/m2/s).
                            ! Always 0: ozone not coded.
                            !
                            ! MDK, 21/10/26
,fo3(land_pts)                                                                 &
                            ! Ozone exposure factor. Always 0: ozone not
                            ! coded.
                            !
                            ! MDK, 21/10/26
,gl(land_pts)                                                                  &
                            ! Leaf conductance for H2O (m/s).
,psi_leaf(land_pts)                                                            &
                            ! Leaf water potential (Pa)
,leaf_k(land_pts)                                                              &
                            ! Whole-path (soil-to-leaf) conductance
                            ! -dE/dpsi_leaf, in the units of kmax (Kplant).
,carbon_gain_out(land_pts)                                                     &
                            ! Carbon gain of the chosen state: A / max A
                            ! (dimensionless) for profit max, CMax, CGain and
                            ! CAP; A (mol CO2/m2/s) for the SOX flat search.
                            ! 0 where the stomata are closed.
                            !
                            ! MDK, 21/10/26
,hydraulic_cost_out(land_pts)
                            ! Hydraulic cost of the chosen state, as a
                            ! fraction (dimensionless): of the k range for
                            ! profit max, of max A for CMax, CGain and CAP,
                            ! and 1 - the factor on A for SOX. 0 where the
                            ! stomata are closed.
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem(land_pts)
                            ! Stem water potential (Pa): the leaf-side end of
                            ! the stem segment of the chosen sample, with
                            ! l_som_plant_segments and xylem impairment;
                            ! psi_root_zone otherwise.
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
 lai_cost(land_pts)
                            ! Leaf area of this call (m2 leaf m-2 ground):
                            ! the CMax / CGain cost is per unit leaf area.
                            ! Absent: kmax_ref / kmax_pft (the leaf area for
                            ! big leaf and multilayer, without kplant_model
                            ! = 1). The two-leaf sunlit and shaded leaves
                            ! pass theirs, as their kmax share is
                            ! N-weighted.
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm) ::                                                      &
 psi_stem_sample(0:som_n_sample, open_pts)
                            ! Stem water potential of each Ci sample.

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 optimal_index                                                                 &
                            ! Index of the optimal sample for a land point.
                            ! Used by the SOX flat search
                            ! (sox_opt_profit_model).
                            !
                            ! MDK, 21/10/26
,i,j,l,iseg,m                                                                  &
                            ! Iterators
,errcode
                            ! Error code to pass to ereport.

INTEGER ::                                                                     &
 optimal_index_flat(open_pts)
                            ! Index of the best sample for the flat search,
                            ! from the shared stom_opt_profit_max_select
                            ! helper.

LOGICAL ::                                                                     &
 l_good_sample(som_n_sample, open_pts)                                        &
                            ! Feasibility mask for the single-pass
                            ! sox_opt_profit_model search.
,l_good_sample_flat(som_n_sample, open_pts)
                            ! Feasibility mask for the flat search.

!-----------------------------------------------------------------------------
! Arrays containing results of sampling over the base parameter. Only need
! to do so for locations with open stomata.
! (profit_max_profit_model with som_ci_search = 2 uses
! stom_opt_bounded_search
! instead, which needs none of the sample arrays.)
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 ci_lo(land_pts)                                                               &
                            ! Lower bound of the Ci sampling range for the
                            ! pass about to be run.
,ci_hi(land_pts)
                            ! Upper bound of the Ci sampling range for the
                            ! pass about to be run.

! -- Flat-search sample arrays (sox_opt_profit_model, and
! -- profit_max_profit_model with som_ci_search = 1) --
REAL(KIND=real_jlslsm) ::                                                      &
 ci_sample(0:som_n_sample, open_pts)                                           &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:som_n_sample, open_pts)                                           &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:som_n_sample, open_pts)                                           &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:som_n_sample, open_pts)                                           &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
                            !
                            ! MDK, 21/10/26
,kl_hc_sample(0:som_n_sample, open_pts)                                        &
                            ! Conductance used for the hydraulic cost
                            ! (mol m-2 s-1 Pa-1) - see stom_opt_mod_ci.
                            !
                            ! MDK, 21/10/26
,kl_SOX(0:som_n_sample, open_pts)                                              &
                            ! Conductance for the SOX_opt model: the whole
                            ! path at the mean of the leaf and root zone
                            ! water potentials (mol m-2 s-1 Pa-1).
                            !
                            ! MDK, 21/10/26
,psi_sample(0:som_n_sample, open_pts)                                          &
                            ! Water potential at leaf (Pa)
,SOX_mean_psi(0:som_n_sample, open_pts)                                        &
                            ! Mean water potential between leaf and root zone
                            ! (Pa). Used in SOX_opt model.
,el_sample(0:som_n_sample, open_pts)                                           &
                            ! Transpiration rate (mol H2O/m2/s)
,profit(0:som_n_sample, open_pts)                                              &
                            ! Profit for each leaf state
,carbon_gain(0:som_n_sample, open_pts)                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(0:som_n_sample, open_pts)                                      &
                            ! Hydraulic cost for each leaf state
,max_kl(open_pts)                                                              &
                            ! Maximum conductance (mol m-2 s-1 Pa-1)
                            !
                            ! MDK, 21/10/26
,kl_hc_max(land_pts)                                                           &
                            ! Reference (maximum) xylem conductance for the
                            ! hydraulic cost with l_xylem_impairment: the
                            ! unimpaired conductance at the root zone water
                            ! potential. JBaguley
,b_ref(land_pts)                                                               &
,c_ref(land_pts)
                            ! Unimpaired PFT conductance_b/c.

! -- stom_opt_bounded_search outputs (scalar per open point) --
REAL(KIND=real_jlslsm) ::                                                      &
 ci_bnd(open_pts)                                                              &
,al_bnd(open_pts)                                                              &
,gl_bnd(open_pts)                                                              &
,kl_bnd(open_pts)                                                              &
,psi_bnd(open_pts)                                                             &
,el_bnd(open_pts)                                                              &
,carbon_gain_bnd(open_pts)                                                     &
,hydraulic_cost_bnd(open_pts)                                                  &
,psi_stem_bnd(open_pts)                                                        &
                            ! Stem segment leaf-side psi at the optimum (xylem
                            ! impairment with segments; else psi_root_zone).
,lambda_bnd(open_pts)
                            ! dA/dE at the optimum (see som_lambda_call).
                            !
                            ! MDK, 22/10/26

! -- Negligible-hydraulic-cost fast path (profit_max_profit_model only) --
! A single-point evaluation at the most water-demanding candidate Ci (near
! ca) - if even that point's hydraulic cost is negligible, the full Ci
! search is skipped for that point (see the fast-path block below).
REAL(KIND=real_jlslsm) ::                                                      &
 ci_e_fp(0:1, open_pts)                                                       &
,al_e_fp(0:1, open_pts)                                                       &
,gl_e_fp(0:1, open_pts)                                                       &
,kl_e_fp(0:1, open_pts)                                                       &
,kl_hc_e_fp(0:1, open_pts)                                                    &
,psi_e_fp(0:1, open_pts)                                                      &
,el_e_fp(0:1, open_pts)                                                       &
,hc_at_maxstress(open_pts)

LOGICAL :: l_fastpath(open_pts)

INTEGER ::                                                                     &
 open_pts_search                                                              &
                            ! Number of open points still needing the full
                            ! Ci search after the fast-path check below
                            ! (== open_pts when l_som_skip_search_wellwatered
                            ! is .FALSE.).
,open_index_search(open_pts)                                           &
                            ! Compressed subset of open_index holding only
                            ! the points not resolved by the fast path.
,open_pts_flat                                                                &
                            ! Number of points to run the flat search over
                            ! (all of open_index_search, or 0 with
                            ! som_ci_search = 2).
,open_index_flat(open_pts)
                            ! Compressed subset of open_index for the flat
                            ! search (same convention as open_index_search).

REAL(KIND=real_jlslsm) :: hc_ref(open_pts)
                            ! Reference conductance for the fast path's
                            ! hydraulic cost: kmax, or kl_hc_max with
                            ! l_xylem_impairment.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
REAL(KIND=real_jlslsm) :: kx_closed(land_pts)
                            ! Impaired xylem conductance at psi_root_zone,
                            ! for the closed-stomata leaf_k.
LOGICAL :: l_xylem_impairment
                            ! l_xylem_impairment_in, but only where an
                            ! impairment model is actually in use (not
                            ! pft_xylem_impairment_model = none, not in
                            ! spin-up), so that model 0 takes exactly the
                            ! unimpaired code path, including with
                            ! l_som_plant_segments (whose hydraulic cost
                            ! uses the segments' whole-plant conductance).

CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_MOD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_xylem_impairment = l_xylem_impairment_in .AND. .NOT. is_spinup .AND.     &
                     pft_xylem_impairment_model(pft) /= xylem_impairment_none

! ----------------------------------------------------------------------------
!  Parameter setup
! ----------------------------------------------------------------------------
ci(:)      = ca(:)

! Set to -rd to get net photosynthesis for points with closed stomata. This
! is overwritten for points with open stomata.
al(:)      = -rd(:)
el(:)      = 0.0
flux_o3(:) = 0.0
fo3(:)     = 0.0
gl(:)      = min_gl_pft(pft)
psi_leaf(:)= psi_root_zone(:)

! Closed stomata: the whole-path conductance at zero flow (xylem at
! psi_root_zone, in series with the soil-to-root ksr_path), the same
! definition as the optimiser's leaf_k for open points.
leaf_k(:)  = kmax
IF (l_xylem_impairment .AND. .NOT. l_som_plant_segments) THEN
  ! Impaired xylem at psi_root_zone (any impairment model), with the
  ! soil-to-root conductance in series.
  CALL leaf_conductance_impaired_jls( pft, land_pts, psi_root_zone, kmax_ref, &
                                      kmax, kcrit, conductance_b,              &
                                      conductance_c, psi_leaf_extreme,         &
                                      kx_closed )
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l),           &
                                 kx = kx_closed(l))
  END DO
ELSE IF (l_xylem_impairment) THEN
  ! Segments (memory model): the stem and leaf segments capped at
  ! kmax / kmax_ref of the intact ones, as leaf_psi_segments_jls.
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax_ref(l), psi_root_zone(l),       &
                  kcap_frac = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm)))
  END DO
ELSE
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))
  END DO
END IF

IF (PRESENT(psi_stem)) psi_stem(:) = psi_root_zone(:)

! Default outputs for points with closed stomata. Overwritten below for
! points with open stomata. Set before the early return, so that a call
! with no open points does not return the previous call's values.
carbon_gain_out(:) = 0.0
hydraulic_cost_out(:) = 0.0

! MDK, 22/10/26
!$OMP CRITICAL (som_lambda_call_alloc)
IF ( .NOT. ALLOCATED(som_lambda_call) ) ALLOCATE( som_lambda_call(land_pts) )
!$OMP END CRITICAL (som_lambda_call_alloc)
som_lambda_call(:) = 0.0

carbon_gain(:,:) = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:) = 0.0
kl_SOX(:,:) = 0.0

! If there are no land points with open stomata then no calculation is
! needed (leaf_k holds the closed-stomata value set above).
IF(0 == open_pts) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  return
END IF

! With xylem impairment the hydraulic cost is measured against the
! unimpaired PFT vulnerability curve, relative to its conductance at the
! current root zone water potential. JBaguley
b_ref(:) = conductance_b_pft(pft)
c_ref(:) = conductance_c_pft(pft)
kl_hc_max(:) = kmax_ref(:)
IF (l_xylem_impairment) THEN
  IF (l_som_plant_segments) THEN
    ! Segments (root, stem, leaf): the reference conductance is that of the
    ! intact plant when it is not transpiring. With no transpiration there
    ! is no drop in water potential along the plant, so every segment sits
    ! at the root zone water potential, and the segments add in series:
    !   1 / kl_hc_max = sum over segments of 1 / k_s(psi_root_zone),
    ! where k_s is the segment's maximum conductance reduced by its own
    ! vulnerability curve. This is the most the undamaged plant could
    ! conduct at today's soil water potential; the hydraulic cost of a
    ! stomatal opening is the loss of conductance from it.
    ! Vegetated points only (psi_root_zone means nothing elsewhere).
    !
    ! MDK, 21/10/26
    DO m = 1, veg_pts
      l = veg_index(m)
      kl_hc_max(l) = 0.0
      DO iseg = 1, 3
        ! The segment's maximum conductance is kmax_ref times its share,
        ! seg_kmax_fac: seg_kfac, and with kplant_model = 1 also the LAI
        ! scaling seg_kscale.
        !
        ! MDK, 21/10/26
        kl_hc_max(l) = kl_hc_max(l) + 1.0 / MAX(kmax_ref(l)                   &
                       * seg_kmax_fac(pft,l,iseg)                              &
                       * EXP( -(psi_root_zone(l)                               &
                                / conductance_b_seg(pft,iseg))                 &
                              **conductance_c_seg(pft,iseg) ),                 &
                       TINY(1.0_real_jlslsm))
      END DO
      kl_hc_max(l) = 1.0 / kl_hc_max(l)
    END DO
  ELSE
    CALL leaf_conductance_jls( pft, land_pts, psi_root_zone, kmax_ref, kcrit,  &
                               b_ref, c_ref, kl_hc_max )
  END IF

  ! The soil-to-root conductance in series (l_som_rhizo_series): the soil is
  ! not damaged, so the stomata see it, as in the unimpaired cost (where
  ! leaf_k is the whole soil-to-leaf path). Zero flow: the series sum.
  IF ( l_som_rhizo_series ) THEN
    DO m = 1, veg_pts
      l = veg_index(m)
      IF ( ksr_path(l) >= 0.0 ) THEN
        IF ( ksr_path(l) > TINY(1.0_real_jlslsm) .AND.                         &
             kl_hc_max(l) > TINY(1.0_real_jlslsm) ) THEN
          kl_hc_max(l) = 1.0 / ( 1.0 / kl_hc_max(l) + 1.0 / ksr_path(l) )
        ELSE
          kl_hc_max(l) = 0.0
        END IF
      END IF
    END DO
  END IF
END IF

! ----------------------------------------------------------------------------
!  Calculate leaf properties as a function of the selected parameter,
!    either ci or psi.
! ----------------------------------------------------------------------------
SELECT CASE ( som_base_parm )
  CASE (som_base_parm_ci)

    SELECT CASE (som_profit_model)

    CASE (profit_max_profit_model)

      !---------------------------------------------------------------------
      ! Negligible-hydraulic-cost fast path: HydraulicCost(ci) is
      ! monotonically non-decreasing in ci (higher ci => more open stomata
      ! => more transpiration => more negative leaf psi => lower xylem k),
      ! so its worst case over [ccp, ca] is at the most water-demanding
      ! candidate, ci closest to ca. If even *that* point's hydraulic cost
      ! is negligible, it is negligible everywhere in range, so
      ! profit(ci) = CarbonGain(ci) - (~0) is just CarbonGain, which keeps
      ! rising with ci - the search would converge on this same near-ca
      ! point anyway, just at O(som_n_sample) times the cost.
      ! Skip it and use this point directly. Points that fail the check
      ! (hydraulic cost not negligible, or even the near-ca point already
      ! infeasible under water stress) fall through to the normal search
      ! below, unaffected.
      !---------------------------------------------------------------------
      IF (l_som_skip_search_wellwatered) THEN
        ci_lo(:) = ca(:) - MAX((ca(:) - MAX(ccp(:), 0.0)) * 0.001_real_jlslsm, &
                               1.0e-3_real_jlslsm)
        ci_hi(:) = ca(:)

        CALL stom_opt_mod_ci(                                                  &
        ! IN
            land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,   &
            1, ci_lo, ci_hi,                                                   &
            rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,   &
            km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,     &
            l_multilayer, r_bl, r_ca, g_oth,                                  &
            kmax_ref, conductance_b, conductance_c,                           &
            psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,           &
        ! OUT
            ci_e_fp, al_e_fp, gl_e_fp, kl_e_fp, psi_e_fp, el_e_fp, kl_hc_e_fp  &
                )
      END IF

      open_pts_search = 0
      DO j = 1, open_pts
        l = veg_index(open_index(j))
        l_fastpath(j) = .FALSE.
        IF (l_xylem_impairment) THEN
          hc_ref(j) = kl_hc_max(l)
        ELSE
          hc_ref(j) = kmax(l)
        END IF

        ! (Not with l_som_nsl: the limitation makes the optimum interior even
        ! where the hydraulic cost is negligible.)
        IF (l_som_skip_search_wellwatered .AND. .NOT. l_som_nsl .AND.          &
            psi_e_fp(1,j) <= psi_root_zone(l) + TINY(psi_root_zone(l)) .AND.   &
            kl_e_fp(1,j) > kcrit(l) .AND.                                      &
            (gl_max(l) <= 0.0 .OR. gl_e_fp(1,j) <= gl_max(l))) THEN
          hc_at_maxstress(j) = (hc_ref(j) - kl_hc_e_fp(1,j)) /                 &
                                MAX(hc_ref(j) - kcrit(l), TINY(1.0_real_jlslsm))
          l_fastpath(j) = hc_at_maxstress(j) < som_hc_negligible_tol
        END IF

        IF (l_fastpath(j)) THEN
          ci(l) = MAX(0.0, ci_e_fp(1,j))
          al(l) = MAX(-rd(l), al_e_fp(1,j))
          gl(l) = MAX(0.0, gl_e_fp(1,j))

          ! gl = 0 means closed: al = -rd, as in set_closed.
          IF (gl(l) <= 0.0) al(l) = -rd(l)
          psi_leaf(l) = psi_e_fp(1,j)
          el(l) = el_e_fp(1,j)
          leaf_k(l) = kl_e_fp(1,j)

          ! The near-ca sample is the max-A one, so A / max A = 1.
          !
          ! MDK, 21/10/26
          carbon_gain_out(l) = 1.0_real_jlslsm
          hydraulic_cost_out(l) = hc_at_maxstress(j)
        ELSE
          open_pts_search = open_pts_search + 1

          ! NOTE: open_index_search must stay in the same veg_index-compressed
          ! convention as open_index (i.e. hold open_index(j), not l =
          ! veg_index(open_index(j))) - it is later passed as the open_index
          ! actual argument to stom_opt_mod_ci/stom_opt_bounded_search, both of
          ! which apply their own veg_index(open_index(j)) remap internally.
          ! Storing the already-remapped land point l here would cause that
          ! remap to be applied twice.
          open_index_search(open_pts_search) = open_index(j)
        END IF
      END DO

      IF (som_ci_search == som_ci_bounded) THEN
        !-------------------------------------------------------------------
        ! Root find for the feasible-range edge, then golden-section within
        ! it - see stom_opt_bounded_search. Never needs the flat search.
        !-------------------------------------------------------------------
        open_pts_flat = 0
        IF (open_pts_search > 0) THEN
          CALL bounded_search( open_pts_search, open_index_search, .TRUE. )
          CALL scatter_bnd( open_pts_search, open_index_search, .TRUE. )
        END IF
      ELSE
        ! Every point not resolved by the fast path goes to the flat search
        ! below.
        open_pts_flat = open_pts_search
        open_index_flat(1:open_pts_search) =                                   &
                                          open_index_search(1:open_pts_search)
      END IF

      !---------------------------------------------------------------------
      ! Flat grid: som_n_sample points over the full [ccp, ca] range, routed
      ! through the shared stom_opt_profit_max_select helper, run over
      ! open_index_flat.
      !---------------------------------------------------------------------
      IF (open_pts_flat > 0) THEN
        ci_lo(:) = MAX(ccp(:), 0.0)
        ci_hi(:) = ca(:)

        CALL stom_opt_mod_ci(                                                  &
        ! IN
            land_pts, pft, open_pts_flat, open_index_flat,                    &
            pft_photo_model, veg_index,                                       &
            som_n_sample, ci_lo, ci_hi,                                       &
            rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,  &
            km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,    &
            l_multilayer, r_bl, r_ca, g_oth,                                  &
            kmax_ref, conductance_b, conductance_c,                           &
            psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,           &
        ! OUT
            ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample,&
            kl_hc_sample,                                                     &
            psi_stem_sample = psi_stem_sample                                 &
                )

        CALL stom_opt_profit_max_select(                                      &
            land_pts, open_pts_flat, open_index_flat, veg_index,              &
            som_n_sample,                                                     &
            al_sample, kl_sample, psi_sample, gl_sample, psi_root_zone, kcrit,&
            gl_max, rd,                                                       &
            kl_hc_sample, kl_hc_max, l_xylem_impairment,                      &
        ! OUT
            l_good_sample_flat, carbon_gain, hydraulic_cost, profit,          &
            optimal_index_flat                                               &
                )

        DO j = 1, open_pts_flat
          l = veg_index(open_index_flat(j))
          ci(l) = MAX(0.0, ci_sample(optimal_index_flat(j),j))

          ! Closed state (index 0) keeps al = -rd, as for points that
          ! were never open; an open optimum is clipped at -rd (as closed).
          IF (optimal_index_flat(j) == 0) THEN
            al(l) = al_sample(0,j)
          ELSE
            al(l) = MAX(-rd(l), al_sample(optimal_index_flat(j),j))
          END IF
          gl(l) = MAX(0.0, gl_sample(optimal_index_flat(j),j))

          ! gl = 0 means closed: al = -rd, as in set_closed.
          IF (gl(l) <= 0.0) al(l) = -rd(l)
          psi_leaf(l) = psi_sample(optimal_index_flat(j),j)
          IF (PRESENT(psi_stem))                                              &
            psi_stem(l) = psi_stem_sample(optimal_index_flat(j),j)
          el(l) = el_sample(optimal_index_flat(j),j)
          leaf_k(l) = kl_sample(optimal_index_flat(j),j)
          carbon_gain_out(l) = carbon_gain(optimal_index_flat(j),j)
          hydraulic_cost_out(l) = hydraulic_cost(optimal_index_flat(j),j)
        END DO
      END IF

    CASE (cmax_profit_model, cgain_profit_model, cap_profit_model,            &
          pm2_profit_model)
      !-------------------------------------------------------------------
      ! The alternatives to profit max, all on the bounded Ci search for
      ! every open point (no fast path; see stom_opt_bounded_search). Each
      ! maximises net assimilation A_n (absolute carbon units) minus a cost:
      !   CMax (Wolf et al. 2016; stomata_model = 6): A_n - Theta(psi_leaf),
      !     a carbon cost of leaf water potential.
      !   CGain (Lu et al. 2020; stomata_model = 7):
      !     A_n - varpi (k_max - k) / k_max, a carbon cost of lost
      !     conductance. varpi (cgain_varpi, umol CO2 m-2 leaf s-1) is the
      !     cost of losing all of the path's conductance.
      !   CAP (Dewar et al. 2017, stomata_model = 10): A_n with the nonstomatal
      !     limitation (l_som_nsl) and no stomatal cost.
      !   ProfitMax2 (Wang et al. 2020, Eqn 11; stomata_model = 11):
      !     A_n (1 - E / E_crit), the cost A_n E / E_crit of the proximity
      !     to hydraulic failure, with no parameter.
      !
      ! References:
      ! -----------
      ! - Wolf, A, Anderegg, W.R.L. & Pacala, S.W. (2016) Optimal stomatal
      !   behavior with competition for water and risk of hydraulic impairment,
      !   Proc. Natl. Acad. Sci. U.S.A. 113 (46) E7222-E7230,
      !   https://doi.org/10.1073/pnas.1615144113
      ! - Lu, Y., Duursma, R.A., Farrior, C.E., Medlyn, B.E. and Feng, X.
      !   (2020), Optimal stomatal drought response shaped by competition for
      !   water and hydraulic risk can explain plant trait covariation. New
      !   Phytol, 225: 1206-1217. https://doi.org/10.1111/nph.16207
      ! - Dewar, R., Mauranen, A., Makela, A., Holtta, T., Medlyn, B., &
      !   Vesala, T. (2017). New insights into the covariation of stomatal,
      !   mesophyll and hydraulic conductances from optimization models
      !   incorporating nonstomatal limitations to photosynthesis. New
      !   Phytologist, 217, 571-585. https://doi.org/10.1111/nph.14848
      ! - Wang, Y., Sperry, J.S., Anderegg, W.R.L., Venturas, M.D. and
      !   Trugman, A.T. (2020), A theoretical and empirical assessment of
      !   stomatal optimization modeling. New Phytol, 227: 311-325.
      !   https://doi.org/10.1111/nph.16572
      !
      ! MDK, 22/10/26
      !-------------------------------------------------------------------
      ! psi_stem as the profit-max bounded path, so the per-segment stem
      ! damage (l_ximpair_seg_memory) follows the stem's own water
      ! potential, not the root zone's.
      !
      ! MDK, 21/10/26
      CALL bounded_search( open_pts, open_index, .TRUE. )
      CALL scatter_bnd( open_pts, open_index, .TRUE. )

    CASE (sox_opt_profit_model)
      !-------------------------------------------------------------------
      ! sox_opt_profit_model: with som_ci_search = 2 and no xylem
      ! impairment the bounded search (as the profit max; A_n (k - kcrit)
      ! is unimodal in Ci), otherwise a single-pass search at som_n_sample
      ! resolution over [ccp, ca]. (The bounded SOX_opt profit is not coded
      ! for the impaired curve: with impairment the flat grid is kept.)
      !-------------------------------------------------------------------
      IF (som_ci_search == som_ci_bounded .AND. .NOT. l_xylem_impairment) THEN
        CALL bounded_search( open_pts, open_index, .FALSE. )
        CALL scatter_bnd( open_pts, open_index, .FALSE. )
      ELSE
        carbon_gain(:,:) = 0.0
        hydraulic_cost(:,:) = 0.0
        profit(:,:) = 0.0
        kl_SOX(:,:) = 0.0
        l_good_sample(:,:) = .TRUE.
        ci_lo(:) = MAX(ccp(:), 0.0)
        ci_hi(:) = ca(:)

        CALL stom_opt_mod_ci(                                                  &
        ! IN
            land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,   &
            som_n_sample, ci_lo, ci_hi,                                        &
            rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,   &
            km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,     &
            l_multilayer, r_bl, r_ca, g_oth,                                   &
            kmax_ref, conductance_b, conductance_c,                            &
            psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,            &
        ! OUT
            ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample, &
            kl_hc_sample,                                                      &
            psi_stem_sample = psi_stem_sample                                  &
                )

        ! ------------------------------------------------------------------
        !  Make mask of those samples that are physically reasonable.
        ! ------------------------------------------------------------------
        DO j = 1, open_pts
          l = veg_index(open_index(j))

          ! Require psi <= psi_root_zone
          l_good_sample(:,j) = l_good_sample(:,j)                              &
              .and. psi_sample(1:,j) <= psi_root_zone(l)                       &
                                        + TINY(psi_root_zone(l))

          ! Require k > k_crit
          l_good_sample(:,j) = l_good_sample(:,j)                              &
                               .and. kl_sample(1:,j) > kcrit(l)

          ! Require gl <= gl_max (when the cap is enabled)
          IF (gl_max(l) > 0.0) l_good_sample(:,j) = l_good_sample(:,j)         &
                               .AND. gl_sample(1:,j) <= gl_max(l)
        END DO

        IF (l_xylem_impairment) THEN
          ! With xylem impairment the SOX hydraulic cost uses the unimpaired
          ! conductance at the leaf water potential, relative to the
          ! unimpaired maximum conductance kmax_ref. JBaguley
          kl_SOX(:,:) = kl_hc_sample(:,:)
          DO j = 1, open_pts
            l = veg_index(open_index(j))
            max_kl(j) = kmax_ref(l)
          END DO
        ELSE
          ! SOX_opt (Eller et al. 2018): the conductance is taken at psi_bar,
          ! the mean of the leaf and root zone water potentials. Here it is the
          ! whole path (root, stem and leaf segments and the soil-to-root link
          ! in series, all at psi_bar; k_path_zero_flow), the same path as
          ! psi_leaf and the k > kcrit test.
          !
          ! MDK, 21/10/26
          DO j = 1, open_pts
            l = veg_index(open_index(j))

            SOX_mean_psi(:,j) = 0.5 * (psi_sample(:,j) + psi_root_zone(l))
            DO i = 1, som_n_sample
              kl_SOX(i,j) = k_path_zero_flow(pft, l, kmax(l), SOX_mean_psi(i,j))
            END DO
          END DO

          ! Get the maximum xylem conductance for each land point.
          max_kl = MAXVAL(kl_SOX(1:,:), MASK = l_good_sample, DIM = 1)
        END IF

        ! CG = An
        carbon_gain(:,:) = al_sample(:,:)

        ! Loop over land points with open stomata
        DO j = 1, open_pts
          l = veg_index(open_index(j))

          ! HC = 1 - (k - kcrit) / (max_k - kcrit)
          !
          ! MDK, 21/10/26
          hydraulic_cost(:,j) = 1.0 - (kl_SOX(:,j)-kcrit(l))                   &
                                      /(max_kl(j)-kcrit(l))

          ! Profit = CG*(1-HC)
          profit(:,j) = carbon_gain(:,j)*(1.0-hydraulic_cost(:,j))
        END DO ! j open points

        ! ------------------------------------------------------------------
        !  Identify maximum profit and set leaf properties to those of the
        !   optimal state. See the NOTE above stom_opt_profit_max_select for
        !   what the zeroth (closed-stomata) sample represents.
        ! ------------------------------------------------------------------
        DO j = 1, open_pts
          l = veg_index(open_index(j))

          optimal_index = MAXLOC(profit(1:,j), DIM = 1,                        &
                                 MASK = l_good_sample(:,j))

          ci(l) = MAX(0.0, ci_sample(optimal_index,j))
          IF (optimal_index == 0) THEN
            al(l) = al_sample(0,j)
          ELSE
            al(l) = MAX(-rd(l), al_sample(optimal_index,j))
          END IF
          gl(l) = MAX(0.0, gl_sample(optimal_index,j))

          ! gl = 0 means closed: al = -rd, as in set_closed.
          IF (gl(l) <= 0.0) al(l) = -rd(l)
          psi_leaf(l) = psi_sample(optimal_index,j)

          ! psi_stem as the profit-max flat search, so the per-segment stem
          ! damage follows the stem's water potential, not the root zone's.
          !
          ! MDK, 21/10/26
          IF (PRESENT(psi_stem)) psi_stem(l) = psi_stem_sample(optimal_index,j)
          el(l) = el_sample(optimal_index,j)
          leaf_k(l) = kl_sample(optimal_index,j)

          carbon_gain_out(l) = carbon_gain(optimal_index,j)
          hydraulic_cost_out(l) = hydraulic_cost(optimal_index,j)
        END DO
      END IF  ! som_ci_search (SOX)

    CASE (supply_loss_model)
      !-------------------------------------------------------------------
      ! Supply-loss stomata (stomata_model = 9; Sperry and Love 2015): no
      ! optimisation. gl_max carries only the soil supply limit; see
      ! stom_opt_supply_loss. Its results come back in the bounded search's
      ! *_bnd arrays.
      !
      ! MDK, 21/10/26
      !-------------------------------------------------------------------
      CALL stom_opt_supply_loss(                                               &
      ! IN
          land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,     &
          rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,     &
          km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,           &
          gl_max, l_multilayer,                                                &
      ! OUT
          ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                             &
          el_bnd, carbon_gain_bnd, hydraulic_cost_bnd                          &
              )

      CALL scatter_bnd( open_pts, open_index, .FALSE. )

    CASE DEFAULT
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
                 'profit_model should be profit_max, SOX_opt, CMax, CGain, ' //&
                 'supply-loss, CAP or ProfitMax2')

    END SELECT ! som_profit_model

  ! (som_base_parm = 2, a search in psi, was never implemented; the check in
  ! jules_vegetation_mod rejects it, so only the Ci search is left.)
  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
               'som_base_parm should be 1 (2, psi, is not implemented)')
END SELECT

! Open points the optimiser closed (gl = 0) get the same minimum
! conductance as points closed by darkness (gl = min_gl_pft above), so a
! night-closed and a drought-closed canopy are treated alike. (They got 0.)
!
! MDK, 21/10/26
DO j = 1, open_pts
  l = veg_index(open_index(j))
  IF ( gl(l) <= 0.0 ) gl(l) = min_gl_pft(pft)
END DO

! leaf_k is Kplant: the whole-path (soil-to-leaf) conductance of the chosen
! state, also with xylem impairment.
! Ozone exposure (flux_o3, fo3) is not calculated (TODO: port from
! leaf_jls_mod).
!
! MDK, 21/10/26

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

!-----------------------------------------------------------------------------
! The bounded Ci search (stom_opt_bounded_search) on n points, idx in the
! open_index convention. l_lai: pass lai_cost (profit max, which does not
! use it, and CMax, CGain and CAP); SOX_opt's bounded search has always run
! without it. These were three copies of the same call.
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
SUBROUTINE bounded_search( n, idx, l_lai )
INTEGER, INTENT(IN) :: n, idx(:)
LOGICAL, INTENT(IN) :: l_lai

IF ( l_lai ) THEN
  CALL stom_opt_bounded_search(                                               &
  ! IN
      land_pts, pft, n, idx,                                                  &
      pft_photo_model, veg_index,                                             &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,        &
      km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,              &
      gl_max, l_multilayer, r_bl, r_ca, g_oth,                                &
      som_n_sample, som_n_ci_golden_iter,                                     &
      kmax_ref, conductance_b, conductance_c, psi_leaf_extreme,               &
      psi_root_extreme, l_xylem_impairment, kl_hc_max,                        &
  ! OUT
      ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                                &
      el_bnd, carbon_gain_bnd, hydraulic_cost_bnd,                            &
      psi_stem_g = psi_stem_bnd, lambda_g = lambda_bnd,                       &
  ! IN (optional)
      lai_cost = lai_cost                                                     &
          )
ELSE
  CALL stom_opt_bounded_search(                                               &
  ! IN
      land_pts, pft, n, idx,                                                  &
      pft_photo_model, veg_index,                                             &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,        &
      km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,              &
      gl_max, l_multilayer, r_bl, r_ca, g_oth,                                &
      som_n_sample, som_n_ci_golden_iter,                                     &
      kmax_ref, conductance_b, conductance_c, psi_leaf_extreme,               &
      psi_root_extreme, l_xylem_impairment, kl_hc_max,                        &
  ! OUT
      ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                                &
      el_bnd, carbon_gain_bnd, hydraulic_cost_bnd,                            &
      lambda_g = lambda_bnd                                                   &
          )
END IF

END SUBROUTINE bounded_search

!-----------------------------------------------------------------------------
! Copy the *_bnd results of n points (idx in the open_index convention) to
! the land-point outputs; l_stem: also psi_stem (not set by
! stom_opt_supply_loss, nor used by SOX_opt). These were four copies.
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
SUBROUTINE scatter_bnd( n, idx, l_stem )
INTEGER, INTENT(IN) :: n, idx(:)
LOGICAL, INTENT(IN) :: l_stem
INTEGER :: j, l

DO j = 1, n
  l = veg_index(idx(j))
  IF (l_stem .AND. PRESENT(psi_stem)) psi_stem(l) = psi_stem_bnd(j)
  ci(l) = ci_bnd(j)
  al(l) = al_bnd(j)
  gl(l) = gl_bnd(j)
  psi_leaf(l) = psi_bnd(j)
  el(l) = el_bnd(j)
  leaf_k(l) = kl_bnd(j)
  carbon_gain_out(l) = carbon_gain_bnd(j)
  hydraulic_cost_out(l) = hydraulic_cost_bnd(j)

  ! (Not set by the supply-loss rule. MDK, 22/10/26)
  IF ( som_profit_model /= supply_loss_model )                                &
    som_lambda_call(l) = lambda_bnd(j)
END DO

END SUBROUTINE scatter_bnd

END SUBROUTINE stom_opt_mod

!-----------------------------------------------------------------------------
! Leaf fluxes at a prescribed transpiration (plant capacitance).
!
! Returns the
!
!   ci, net assimilation al, conductance gl and transpiration el
!
! at which the leaf transpires e_target (mol H2O m-2 s-1, as el), from the
! same photosynthesis, gl(ci) and canopy coupling as stom_opt_mod_ci.
!
! E rises with ci, so:
!
!   1. the ci grid over [ccp, ca] is searched for the samples either side
!      of e_target;
!   2. a second grid inside that bracket refines it;
!   3. ci, al, gl and el are interpolated linearly in E.
!
! Samples above gl_max (when > 0, the soil-supply cap) are not used; above
! the largest usable E that sample is returned. The hydraulics are not
! applied: the capacitance store pays the difference from the optimiser's
! steady state.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_at_e(                                                      &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp,              &
        pstar, km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax,         &
        kcrit, gl_max, l_multilayer, r_bl, r_ca, g_oth,                        &
        kmax_ref, conductance_b, conductance_c,                                &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment, e_target,      &
! OUT
        ci, al, gl, el                                                         &
)

USE jules_vegetation_mod, ONLY: som_n_sample

IMPLICIT NONE

LOGICAL, INTENT(IN) :: l_multilayer
                            ! Called from the multilayer canopy.
LOGICAL, INTENT(IN) :: l_xylem_impairment
                            ! The impaired vulnerability curve is used (as
                            ! stom_opt_mod).
INTEGER, INTENT(IN) :: land_pts
                            ! Number of land points.
INTEGER, INTENT(IN) :: pft
                            ! Plant functional type.
INTEGER, INTENT(IN) :: open_pts
                            ! Number of points to solve.
INTEGER, INTENT(IN) :: open_index(land_pts)
                            ! The points to solve, as indices of veg_index.
INTEGER, INTENT(IN) :: pft_photo_model
                            ! Photosynthesis model of the PFT.
INTEGER, INTENT(IN) :: veg_index(land_pts)
                            ! Land index of each vegetated point.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),          &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts),   &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), qs(land_pts),    &
  je(land_pts), t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts),      &
  ipar(land_pts), kmax(land_pts), kcrit(land_pts), r_bl(land_pts),             &
  r_ca(land_pts), g_oth(land_pts), kmax_ref(land_pts),                         &
  conductance_b(land_pts), conductance_c(land_pts),                            &
  psi_leaf_extreme(land_pts), psi_root_extreme(land_pts)
                            ! Photosynthesis, coupling and hydraulic inputs,
                            ! exactly as passed to stom_opt_mod_ci.
REAL(KIND=real_jlslsm), INTENT(IN) :: gl_max(land_pts)
                            ! Soil-supply cap on gl (m s-1; <= 0: none).
REAL(KIND=real_jlslsm), INTENT(IN) :: e_target(land_pts)
                            ! Target transpiration (mol H2O m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: ci(land_pts)
                            ! Intercellular CO2 (Pa).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: al(land_pts)
                            ! Net assimilation (mol CO2 m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: gl(land_pts)
                            ! Conductance (m s-1).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: el(land_pts)
                            ! Transpiration (mol H2O m-2 s-1). Only the
                            ! points solved are set.
                            !
                            ! MDK, 23/10/26

INTEGER :: ipass, i, j, l
INTEGER :: ib(open_pts)
                            ! The first usable sample at or above the target
                            ! E, at each point.
INTEGER :: n_ok(open_pts)
                            ! The number of usable samples (within gl_max),
                            ! at each point.
REAL(KIND=real_jlslsm) :: ci_lo(land_pts), ci_hi(land_pts)
                            ! The ci range searched (Pa).
REAL(KIND=real_jlslsm) :: w
                            ! Interpolation weight between two samples (-).
REAL(KIND=real_jlslsm), DIMENSION(0:som_n_sample, open_pts) ::                &
  ci_s, al_s, gl_s, kl_s, psi_s, el_s, klhc_s
                            ! The ci samples and their A, gl, conductance,
                            ! psi_leaf and E (as stom_opt_mod_ci).
                            !
                            ! MDK, 23/10/26

IF ( open_pts == 0 ) RETURN

! Pass 1 searches the whole ci range [ccp, ca]; pass 2 the bracket around
! the target E that pass 1 found.
ci_lo(:) = MAX(ccp(:), 0.0)
ci_hi(:) = ca(:)
DO ipass = 1,2
  CALL stom_opt_mod_ci(                                                        &
      land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,         &
      som_n_sample, ci_lo, ci_hi,                                              &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,         &
      km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,           &
      l_multilayer, r_bl, r_ca, g_oth,                                         &
      kmax_ref, conductance_b, conductance_c,                                  &
      psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                  &
      ci_s, al_s, gl_s, kl_s, psi_s, el_s, klhc_s )
  DO j = 1,open_pts
    l = veg_index(open_index(j))

    ! Usable samples: 1..n_ok (ci rising), those within gl_max.
    n_ok(j) = som_n_sample
    IF ( gl_max(l) > 0.0 ) THEN
      DO i = 1,som_n_sample
        IF ( gl_s(i,j) > gl_max(l) ) THEN
          n_ok(j) = MAX(i - 1, 1)
          EXIT
        END IF
      END DO
    END IF

    ! First usable sample at or above the target.
    ib(j) = n_ok(j)
    DO i = 1,n_ok(j)
      IF ( el_s(i,j) >= e_target(l) ) THEN
        ib(j) = i
        EXIT
      END IF
    END DO

    ! After pass 1: narrow the ci range to the bracket around the target.
    IF ( ipass == 1 ) THEN
      i = MAX(ib(j), 2)
      ci_lo(l) = ci_s(i-1,j)
      ci_hi(l) = ci_s(i,j)
    END IF
  END DO
END DO

! Interpolate linearly in E between the two samples either side of the
! target. If the target is below the first sample, or above the largest
! usable E (the soil-supply cap or ci = ca), that sample is returned.
DO j = 1,open_pts
  l = veg_index(open_index(j))
  i = ib(j)
  IF ( i == 1 .OR. el_s(i,j) < e_target(l) ) THEN
    ci(l) = ci_s(i,j)
    al(l) = al_s(i,j)
    gl(l) = gl_s(i,j)
    el(l) = el_s(i,j)
  ELSE
    w = (e_target(l) - el_s(i-1,j))                                            &
        / MAX(el_s(i,j) - el_s(i-1,j), TINY(1.0_real_jlslsm))
    ci(l) = ci_s(i-1,j) + w * (ci_s(i,j) - ci_s(i-1,j))
    al(l) = al_s(i-1,j) + w * (al_s(i,j) - al_s(i-1,j))
    gl(l) = gl_s(i-1,j) + w * (gl_s(i,j) - gl_s(i-1,j))
    el(l) = el_s(i-1,j) + w * (el_s(i,j) - el_s(i-1,j))
  END IF
  gl(l) = MAX(gl(l), 0.0)
  al(l) = MAX(al(l), -rd(l))
END DO

END SUBROUTINE stom_opt_at_e

!-----------------------------------------------------------------------------
! Stomatal optimisation model that iterates over leaf ci
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_mod_ci(                                                    &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        n_sample, ci_lo, ci_hi,                                                &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp,              &
        pstar, km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax,         &
        kcrit, l_multilayer, r_bl, r_ca, g_oth,                               &
        kmax_ref, conductance_b, conductance_c,                               &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,               &
! OUT
        ci_sample, al_sample, gl_sample, kl_sample, psi_sample,el_sample,      &
        kl_hc_sample,                                                          &
! OUT (optional)
        psi_stem_sample                                                        &
)

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls, xylem_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_psi_impaired_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_segments_jls, ksr_path,      &
                                       som_psi_in_min, xylem_f,                &
                                       leaf_psi_segments_lut_jls

                                       ! MDK, 18/10/26
                                       ! MDK, 21/10/26

USE jules_vegetation_mod, ONLY:                                                &
        photo_collatz, photo_farquhar, photo_johnson, photo_model,             &
        som_psi_solver, psi_solver_lut,                                        &
        l_som_plant_segments, l_som_rhizo_series, l_som_coupled_e

USE jb_photo_mod, ONLY: jb_eta_scale

USE pftparm, ONLY:                                                             &
        c3, alpha, conductance_b_pft, conductance_c_pft

USE jules_surface_mod, ONLY: fwe_c3, fwe_c4
USE jules_surface_mod, ONLY: beta1, beta2, ratio

USE planet_constants_mod, ONLY: repsilon

USE c_rmol, ONLY: rmol

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

LOGICAL, INTENT(IN) :: l_multilayer
LOGICAL, INTENT(IN) :: l_xylem_impairment
                            ! See stom_opt_mod.

!-----------------------------------------------------------------------------
! IN integer variables
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
!                           ! Number of land points to be processed
, pft                                                                          &
                            ! Plant functional type
, open_pts                                                                     &
                            ! number of land points with open stomata
, open_index(land_pts)                                                         &
                            ! Index of each land point with open stomata
, pft_photo_model                                                              &
                            ! Indicates which photosynthesis model to use for
                            ! the current PFT.
, veg_index(land_pts)                                                          &
                            ! Index of vegetated points on the land grid.
, n_sample
                            ! Number of Ci sample points to use for this
                            ! call. Lets this routine be reused at different
                            ! resolutions/ranges in stom_opt_mod rather
                            ! than being fixed to the module-level
                            ! som_n_sample.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_bl(land_pts), r_ca(land_pts), g_oth(land_pts)
                            ! Coupled transpiration (see stom_opt_mod).

!-----------------------------------------------------------------------------
! IN real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 rd(land_pts)                                                                  &
                             ! Dark respiration (mol CO2/m2/s).
,ca(land_pts)                                                                  &
                             ! Atmospheric CO2 pressure (Pa)
,psi_root_zone(land_pts)                                                       &
                             ! Water potential in the root zone (Pa)
,ci_lo(land_pts)                                                               &
                             ! Lower bound of the Ci sampling range (Pa).
,ci_hi(land_pts)                                                               &
                             ! Upper bound of the Ci sampling range (Pa).
,acr(land_pts)                                                                 &
                             ! Absorbed PAR (mol photons/m2/s).
,apar(land_pts)                                                                &
                             ! Absorbed PAR (W m-2).
,oi(land_pts)                                                                  &
                             ! Intercellular oxygen pressure (Pa).
,vcmax(land_pts)                                                               &
                             ! Maximum rate of carboxylation of Rubisco
                             !  (mol CO2/m2/s).
,kc(land_pts)                                                                  &
                             ! Michaelis-Menten constant for CO2 (Pa).
,ko(land_pts)                                                                  &
                             ! Michaelis-Menten constant for O (Pa).
,ccp(land_pts)                                                                 &
                             ! Photorespiratory compensatory point (Pa).
,pstar(land_pts)                                                               &
                             ! Atmospheric pressure (Pa).
,km(land_pts)                                                                  &
                             ! A combination of Michaelis-Menten and other
                             ! terms. Only used with the Farquhar model.
,dq(land_pts)                                                                  &
                            ! Specific humidity deficit of air (kg H2O/kg air).
                            ! Leaf saturation vapour pressure (Pa).
,qs(land_pts)                                                                  &

                            ! Saturated specific humidity
                            !      (kg H2O/kg air).
,je(land_pts)                                                                  &
                             ! Electron transport rate (mol m-2 s-1)
,t_leaf(land_pts)                                                              &
                             ! Leaf temperature (K)
,je_ratio(land_pts)                                                            &
,fapar_lf(land_pts)                                                            &
,ipar (land_pts)                                                               &
,kmax(land_pts)                                                                &
,kcrit(land_pts)                                                               &
,kmax_ref(land_pts)                                                            &
,conductance_b(land_pts)                                                       &
,conductance_c(land_pts)                                                       &
,psi_leaf_extreme(land_pts)                                                    &
,psi_root_extreme(land_pts)
                            ! See stom_opt_mod.
!-----------------------------------------------------------------------------
! OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 ci_sample(0:n_sample, open_pts)                                             &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:n_sample, open_pts)                                             &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:n_sample, open_pts)                                             &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:n_sample, open_pts)                                             &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
                            !
                            ! MDK, 21/10/26
,psi_sample(0:n_sample, open_pts)                                            &
                            ! Water potential at leaf (Pa)
,el_sample(0:n_sample, open_pts)                                             &
                            ! Transpiration rate (mol H2O/m2/s)
,kl_hc_sample(0:n_sample, open_pts)
                            ! Conductance used for the hydraulic cost
                            ! (mol m-2 s-1 Pa-1): kl_sample, or with
                            ! l_xylem_impairment the conductance at
                            ! psi_sample on the unimpaired PFT curve.
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem_sample(0:n_sample, open_pts)
                            ! Stem water potential (Pa) of each sample: the
                            ! stem segment's leaf-side end with segments and
                            ! impairment, else psi_root_zone.

!-----------------------------------------------------------------------------
! Local integer variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 i,j,l                                                                         &
                            ! Internal iterators
,errcode
                            ! Error code to pass to ereport.

!-----------------------------------------------------------------------------
! Local real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 el_closed(1, open_pts)                                                      &
,psi_closed(1, open_pts)                                                     &
,kl_closed(1, open_pts)                                                      &
                            ! Zero-flow solve for the closed-stomata state
                            ! (entry 0): transpiration, leaf water
                            ! potential and xylem conductance.
,wcarb_sample(0:n_sample, open_pts)                                          &
                            ! Rubisco limited photosynthesis
!                           ! rate (mol CO2/m2/s).
,wlite_sample(0:n_sample, open_pts)                                          &
                            ! Electron transport limited photosynthesis
!                           ! rate (mol CO2/m2/s).
,wexpt_sample(0:n_sample, open_pts)                                          &
                            ! Export-limited gross photosynthetic rate
!                           ! (mol CO2/m2/s). Not used with Farquhar model.
,wl_sample(0:n_sample, open_pts)                                             &
                            ! Gross photosynthetic rate (mol CO2/m2/s).
,b1(0:n_sample, open_pts)                                                    &
,b2(0:n_sample, open_pts)                                                    &
,b3(0:n_sample, open_pts)                                                    &
                            ! Temporary variables for use in the quadratic
                            ! equation.
,beta1p2m4, beta2p2m4                                                          &
                            ! beta[12] ** 2 * 4.
,vpd(land_pts)                                                                 &
                            ! Vapour pressure deficit (Pa)
,vpd_factor                                                                    &
                            ! Factor used to convert specific humidity to
                            ! vapour pressure.
,glco2_sample(0:n_sample, open_pts)                                        &
                            ! Leaf conductance for CO2 (m/s)
,conv(land_pts)                                                                &
                            ! Factor for converting mol/m3 into Pa (J/m3)
,kmax_open(open_pts), kcrit_open(open_pts), b_open(open_pts), c_open(open_pts) &
                            ! The unimpaired PFT curve (kmax_ref, kcrit,
                            ! conductance_b/c_pft) gathered onto the
                            ! open-point index, for the conductance of the
                            ! hydraulic cost (kl_hc_sample).
                            !
                            ! MDK, 21/10/26
,kcap_frac_pts(land_pts)                                                       &
                            ! Segments with xylem impairment: impaired over
                            ! intact kmax (kmax / kmax_ref), the cap on the
                            ! stem and leaf segments.
                            !
                            ! MDK, 21/10/26
,psi_stem_work(0:n_sample, open_pts)
                            ! Stem water potential of each sample (Pa).

REAL(KIND=real_jlslsm) :: k_s, psi_in_s
                            ! l_som_rhizo_series: the path's soil-to-root
                            ! conductance and a sample's root water potential
                            ! (Pa), for the impaired hydraulic cost.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_MOD_CI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Setup parameters.
!-----------------------------------------------------------------------------

al_sample(:,:)      = 0.0
gl_sample(:,:)      = 0.0
kl_sample(:,:)      = 0.0
el_sample(:,:)      = 0.0
kl_hc_sample(:,:)   = 0.0

! Entry 0 of the sample arrays represents the closed-stomata state used as
! the fallback in stom_opt_mod (ci = ca, psi = psi_root_zone; see the note
! there). ci_sample and psi_sample are set here so the whole-column
! photosynthesis calculations below have a valid ci for it, but those
! calculations DO overwrite al/gl (and el) for entry 0 with the open-stomata
! values at ci = ca - so entry 0 is reset to the closed state again at the
! end of this routine, after all of them.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  ci_sample(0,j)  = ca(l)
  psi_sample(0,j) = psi_root_zone(l)

  DO i = 1, n_sample
    !NOTE: Array entry 1 is at ci_lo (hence i-1), and the last step is less
    !       than ci_hi to avoid infinite transpiration rates when ci_hi is
    !       the full-range bound ca(l) (see the callers in stom_opt_mod).
    ci_sample(i,j) = ci_lo(l) + (i-1) * (ci_hi(l) - ci_lo(l)) / n_sample
  END DO
END DO

!-----------------------------------------------------------------------------
! Calculate photosynthetic rates.
! Code based on that in leaf_limits_mod.F90 and leaf.F90 with loops over
! sample_ci added.
!-----------------------------------------------------------------------------
SELECT CASE ( pft_photo_model )
CASE ( photo_collatz )
  !---------------------------------------------------------------------------
  ! Use the Collatz models (for C3 or C4 plants).
  !---------------------------------------------------------------------------

  IF (c3(pft) == 1) THEN

!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,ci_sample,       &
!$OMP        ccp,kc,oi,ko,wlite_sample,pft,wexpt_sample,fwe_c3,alpha,acr,      &
!$OMP        n_sample,l_multilayer,apar,fapar_lf,ipar)
    DO j = 1,open_pts
      l = veg_index(open_index(j))

      ! The numbers in these equations are from Cox, HCTN 24,
      ! "Description ... Vegetation Model", equations 54 and 55.
      wcarb_sample(:,j) = vcmax(l) * (ci_sample(:,j) - ccp(l))                 &
                          / (ci_sample(:,j) + kc(l) * (1.0 + oi(l) / ko(l)))
      wlite_sample(:,j) = alpha(pft) * acr(l)                                  &
                          * (ci_sample(:,j) - ccp(l))                          &
                          / (ci_sample(:,j) + 2.0 * ccp(l))
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))

      ! Multilayer calls: rescale the light-limited rate from the canopy's
      ! absorbed PAR (apar) to that of this layer's sunlit or shaded leaves
      ! (fapar_lf * ipar).
      !
      ! MDK, 21/10/26
      IF (l_multilayer) THEN
        wlite_sample(:,j) = wlite_sample(:,j) / apar(l) * fapar_lf(l)          &
                            * ipar(l)
      END IF
      wexpt_sample(:,j) = fwe_c3 * vcmax(l)
    END DO
!$OMP END PARALLEL DO

  ELSE
    !  C4
!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,wlite_sample,    &
!$OMP        pft,wexpt_sample,pstar,alpha,fwe_c4,ci_sample,acr,n_sample,       &
!$OMP        l_multilayer,apar,fapar_lf,ipar)
    DO j = 1,open_pts
      l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l)
      wlite_sample(:,j) = alpha(pft) * acr(l)
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))

      ! Multilayer calls: rescale the light-limited rate from the canopy's
      ! absorbed PAR (apar) to that of this layer's sunlit or shaded leaves
      ! (fapar_lf * ipar).
      !
      ! MDK, 21/10/26
      IF (l_multilayer) THEN
        wlite_sample(:,j) = wlite_sample(:,j) / apar(l) * fapar_lf(l)          &
                            * ipar(l)
      END IF
      wexpt_sample(:,j) = fwe_c4 * vcmax(l) * ci_sample(:,j) / pstar(l)
    END DO
!$OMP END PARALLEL DO

  END IF  !  c3

  !---------------------------------------------------------------------------
  ! Calculate the co-limited rate of gross photosynthesis.
  !---------------------------------------------------------------------------
  beta1p2m4 = 4 * beta1 * beta1
  beta2p2m4 = 4 * beta2 * beta2

  b1(:,:) = beta1
  b2(:,:) = - (wcarb_sample(:,:) + wlite_sample(:,:))
  b3(:,:) = wcarb_sample(:,:) * wlite_sample(:,:)

  wl_sample(:,:) = -b2(:,:) / (2.0 * b1(:,:))                                  &
          - SQRT(b2(:,:) * b2(:,:) / beta1p2m4      - b3(:,:) / b1(:,:))

  b1(:,:) = beta2
  b2(:,:) = - (wl_sample(:,:) + wexpt_sample(:,:))
  b3(:,:) = wl_sample(:,:) * wexpt_sample(:,:)

  wl_sample(:,:) = -b2(:,:) / (2.0 * b1(:,:))                                  &
          - SQRT(b2(:,:) * b2(:,:) / beta2p2m4       - b3(:,:) / b1(:,:))

CASE ( photo_farquhar )
  !---------------------------------------------------------------------------
  ! Use the Farquhar model (for C3 plants only).
  !---------------------------------------------------------------------------
!$OMP PARALLEL DO IF(open_pts > 1)                                             &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(NONE)                                                            &
!$OMP PRIVATE(l,j,i)                                                           &
!$OMP SHARED(open_pts,veg_index,open_index,wcarb_sample,vcmax,wlite_sample,    &
!$OMP        ci_sample,ccp,km,je,n_sample,photo_model,l_multilayer,je_ratio)
  DO j = 1,open_pts
    l = veg_index(open_index(j))
    wcarb_sample(:,j) = vcmax(l) * ( ci_sample(:,j) - ccp(l) )                 &
                        / ( ci_sample(:,j) + km(l) )
    wlite_sample(:,j) = je(l) / 4.0 * ( ci_sample(:,j) - ccp(l) )              &
                        / ( ci_sample(:,j) + 2.0 * ccp(l) )
    ! Johnson-Berry: J depends on ci (see jb_photo_mod).
    IF ( photo_model == photo_johnson )                                        &
      wlite_sample(:,j) = wlite_sample(:,j)                                    &
                          * jb_eta_scale( ccp(l), ci_sample(:,j) )
    wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))

    ! Multilayer calls: the electron transport of this layer's sunlit or
    ! shaded leaves (je_ratio).
    !
    ! MDK, 21/10/26
    IF (l_multilayer) THEN
      wlite_sample(:,j) = wlite_sample(:,j) * je_ratio(l)
    END IF
  END DO
!$OMP END PARALLEL DO

  !---------------------------------------------------------------------------
  ! Calculate minimum of Rubisco and electron transport limited rates.
  !---------------------------------------------------------------------------
  wl_sample(:,:) = MIN( wcarb_sample(:,:), wlite_sample(:,:) )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'pft_photo_model should be photo_collatz or photo_farquhar')

END SELECT  !  pft_photo_model

! ----------------------------------------------------------------------------
! Calculate net photosynthetic rate (mol CO2 m-2 s-1).
!  al = wl - rd
! ----------------------------------------------------------------------------
DO j = 1,open_pts
  l = veg_index(open_index(j))

  al_sample(:,j) = wl_sample(:,j) - rd(l)
END DO

! ----------------------------------------------------------------------------
! Calculate stomatal conductance for H2O (m/s).
!  gl = ratio A R T_leaf / (ca - ci), with ratio = 1.6 the ratio of the
!  diffusivities of H2O and CO2 in air.
!
!  MDK, 21/10/26
! ----------------------------------------------------------------------------
DO j = 1,open_pts
  l = veg_index(open_index(j))

  !---------------------------------------------------------------------------
  ! Calculate the factor for converting mol/m3 into Pa (J/m3).
  !---------------------------------------------------------------------------
  conv(l) = rmol * t_leaf(l)

  !---------------------------------------------------------------------------
  ! Diagnose the leaf conductance
  !
  ! NOTE: (ca-ci) is designed to stay bounded away from 0 by the "n-1 step"
  ! construction of ci_sample above, but when a caller sets ci_hi(l) close
  ! to (or exactly at) ca(l) for a *narrow* range - as the fast-path single
  ! sample near ca does - this margin can shrink far more than the flat
  ! single-pass grid ever produced, overflowing this division to Infinity in
  ! practice.
  ! Floor the denominator so this can only ever produce a large-but-finite
  ! gl, which the psi/k feasibility mask downstream then correctly rejects
  ! as infeasible rather than an Infinity that can corrupt reductions and
  ! propagate into the canopy-level conductance passed out of this model.
  !---------------------------------------------------------------------------
  glco2_sample(:,j) = (al_sample(:,j) * conv(l))                             &
                     / MAX(ca(l) - ci_sample(:,j), 1.0e-2_real_jlslsm)
  gl_sample(:,j)    = ratio * glco2_sample(:,j)

END DO

! ----------------------------------------------------------------------------
! Calculate transpiration rate.
!   E = vpd * gl / pstar
! ----------------------------------------------------------------------------

! Calculate the factor for converting specific humidity to vapour pressure.
! Taken from leaf_limits_mod.F90, the 1.0e3 has been removed to give vpd
! in Pa.
vpd_factor = 1.0 / repsilon

DO j = 1,open_pts
  l = veg_index(open_index(j))

  ! Convert the specific humidity to vapour pressure. This formula is taken
  ! from leaf_limits_mod.F90.
  vpd(l) = dq(l) * pstar(l) * vpd_factor

  ! Calculate transpiration rate in m/s
  el_sample(:,j) = vpd(l) * gl_sample(:,j) / pstar(l)

  ! Convert from m/s to mol H2O/m2/s
  el_sample(:,j) = el_sample(:,j) * pstar(l) / (rmol * t_leaf(l))

  ! Coupled: E = D g_eff(g), dq the deficit to the level-1 air.
  IF ( l_som_coupled_e ) el_sample(:,j) = vpd(l)                               &
       * g_eff_e(MAX(gl_sample(:,j), 0.0_real_jlslsm), r_bl(l), r_ca(l),     &
                 g_oth(l)) / (rmol * t_leaf(l))

END DO

! Transpiration can't be negative
el_sample(:,:) = MAX(0.0, el_sample(:,:))

! Approximate the leaf water potential and conductance.
! NOTE: el_sample/psi_sample/kl_sample are dimensioned (0:n_sample,
!       open_pts) here (entry 0 is the closed-stomata state, handled above),
!       but leaf_psi_jls and everything it calls declare their dummy
!       arguments explicit-shape as (n_sample, open_pts). Passing the
!       whole arrays would silently misalign columns/rows via Fortran
!       sequence association (each open point's data would be read shifted
!       by one row, and bleed into the next open point's column). Passing
!       the (1:n_sample,:) sections instead makes the actual and dummy
!       argument shapes match exactly, so entry 0 is left untouched and
!       entries 1:n_sample line up correctly.
IF (l_xylem_impairment .AND. l_som_plant_segments) THEN
  ! Segments: the memory model's cap on the stem and leaf segments (the
  ! root is not damaged); the hydraulic cost uses the intact segments'
  ! whole-plant conductance at the same water potentials.
  kcap_frac_pts(:) = 1.0
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    kcap_frac_pts(l) = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm))
  END DO

  ! som_psi_solver = 3: the same capped segments from the supply tables.
  !
  ! MDK, 18/10/26
  IF ( som_psi_solver == psi_solver_lut ) THEN
    CALL leaf_psi_segments_lut_jls( pft, n_sample, land_pts, open_pts,         &
                              veg_index, open_index, el_sample(1:n_sample,:),  &
                              psi_root_zone, kmax_ref,                         &
                              psi_sample(1:n_sample,:),                        &
                              kl_sample(1:n_sample,:),                         &
                              kcap_frac = kcap_frac_pts,                       &
                              leaf_k_intact = kl_hc_sample(1:n_sample,:),      &
                              stem_psi = psi_stem_work(1:n_sample,:) )
  ELSE
  CALL leaf_psi_segments_jls( pft, n_sample, land_pts, open_pts, veg_index,   &
                              open_index, el_sample(1:n_sample,:),             &
                              psi_root_zone, kmax_ref, kcrit,                  &
                              psi_sample(1:n_sample,:),                        &
                              kl_sample(1:n_sample,:),                         &
                              kcap_frac = kcap_frac_pts,                       &
                              leaf_k_intact = kl_hc_sample(1:n_sample,:),      &
                              stem_psi = psi_stem_work(1:n_sample,:) )
  END IF
ELSE IF (l_xylem_impairment) THEN
  ! Leaf water potential and conductance on the impaired vulnerability
  ! curve. JBaguley
  CALL leaf_psi_impaired_jls( pft,                                             &
                              n_sample,                                        &
                              land_pts,                                        &
                              open_pts,                                        &
                              veg_index,                                       &
                              open_index,                                      &
                              el_sample(1:n_sample,:),                         &
                              psi_root_zone,                                   &
                              kmax_ref,                                        &
                              kmax,                                            &
                              kcrit,                                           &
                              conductance_b,                                   &
                              conductance_c,                                   &
                              psi_leaf_extreme,                                &
                              psi_root_extreme,                                &
                            ! INTENT OUT
                              psi_sample(1:n_sample,:),                        &
                              kl_sample(1:n_sample,:)                          &
    )

  ! The hydraulic cost uses the unimpaired conductance at the (impaired)
  ! leaf water potential. JBaguley
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    kmax_open(j)  = kmax_ref(l)
    kcrit_open(j) = kcrit(l)
  END DO
  b_open(:) = conductance_b_pft(pft)
  c_open(:) = conductance_c_pft(pft)

  CALL xylem_conductance_jls( pft, n_sample, open_pts,                         &
                              psi_sample(1:n_sample,:),                        &
                              kmax_open, kcrit_open, b_open, c_open,           &
                              kl_hc_sample(1:n_sample,:) )

  ! With the soil-to-root conductance in series (l_som_rhizo_series), the
  ! cost uses the whole path as the unimpaired cost does: the intact xylem
  ! at the same water potentials, in series with the soil (which is not
  ! damaged, so the stomata see it). As leaf_psi_jls, the marginal
  ! conductance of the path is k(psi_leaf) / (1 + k(psi_in) / K_s), with the
  ! root water potential_in = psi_src - E / K_s and k the intact curve.
  IF ( l_som_rhizo_series ) THEN
    DO j = 1, open_pts
      l = veg_index(open_index(j))
      k_s = ksr_path(l)
      IF ( k_s < 0.0 ) CYCLE
      DO i = 1, n_sample
        IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
          psi_in_s = MAX(psi_root_zone(l) - el_sample(i,j) / k_s,             &
                         som_psi_in_min)
          kl_hc_sample(i,j) = kl_hc_sample(i,j)                                &
                              / (1.0 + kmax_ref(l) * xylem_f(pft, psi_in_s)    &
                                       / k_s)
        ELSE IF ( el_sample(i,j) > 0.0 ) THEN
          kl_hc_sample(i,j) = 0.0
        END IF
      END DO
    END DO
  END IF
ELSE
  CALL leaf_psi_jls( pft,                                                      &
                     n_sample,                                                 &
                     land_pts,                                                 &
                     open_pts,                                                 &
                     veg_index,                                                &
                     open_index,                                               &
                     el_sample(1:n_sample,:),                                  &
                     psi_root_zone,                                            &
                     kmax,                                                     &
                     kcrit,                                                    &
                     conductance_b,                                            &
                     conductance_c,                                            &
                   ! INTENT OUT
                     psi_sample(1:n_sample,:),                                 &
                     kl_sample(1:n_sample,:)                                   &
    )

  kl_hc_sample(:,:) = kl_sample(:,:)
END IF

! Reset entry 0 to the closed-stomata state. The photosynthesis and
! conductance calculations above ran over the whole column, so they also
! filled entry 0 with the open state at ci = ca. Closed stomata have no
! conductance or transpiration, net assimilation -rd, and the leaf at the
! root zone water potential.
!
! The conductance of the closed state is that of the whole path when the
! plant is not transpiring, K(psi_root_zone). It is solved at E = 0 on
! the same path as the samples, so it includes the soil-to-root link and
! any impairment. It is not 0, which would report 100 % loss of
! conductance (PLC_pft) whenever the stomata shut, e.g. at twilight.
! Entry 0 is never a profit candidate, so this only changes the reported
! leaf_k.
!
! MDK, 21/10/26
el_closed(:,:) = 0.0

CALL som_path_psi( pft, 1, land_pts, open_pts, veg_index, open_index,          &
                   el_closed, psi_root_zone, kmax_ref, kmax, kcrit,            &
                   conductance_b, conductance_c, psi_leaf_extreme,             &
                   psi_root_extreme, l_xylem_impairment,                       &
                 ! INTENT OUT
                   psi_closed, kl_closed )

! Reset entry 0 to the closed-stomata state (the fallback when no Ci sample
! is worth opening for): net assimilation is dark respiration only (-rd),
! no conductance and no transpiration, the plant conductance is the zero-
! flow value from som_path_psi above, and with no flow the leaf sits at the
! root-zone water potential.
!
! Stem water potential samples (only when the caller asks for them): with
! xylem impairment on the segmented path, each Ci sample keeps the stem
! (root-stem node) potential from its own path solve; otherwise there is no
! separate stem node, so it is the root-zone potential. Entry 0 (closed, no
! flow) is the root-zone potential either way.
!
! MDK, 23/10/26
DO j = 1, open_pts
  l = veg_index(open_index(j))
  al_sample(0,j)  = -rd(l)
  gl_sample(0,j)  = 0.0
  el_sample(0,j)  = 0.0
  kl_sample(0,j)  = kl_closed(1,j)
  kl_hc_sample(0,j) = kl_sample(0,j)
  psi_sample(0,j) = psi_root_zone(l)
END DO

IF (PRESENT(psi_stem_sample)) THEN
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    IF (l_xylem_impairment .AND. l_som_plant_segments) THEN
      psi_stem_sample(1:,j) = psi_stem_work(1:,j)
    ELSE
      psi_stem_sample(1:,j) = psi_root_zone(l)
    END IF
    psi_stem_sample(0,j) = psi_root_zone(l)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod_ci

!-----------------------------------------------------------------------------
! Leaf water potential and whole-path conductance for a given transpiration
! e_leaf,
!
!   psi_leaf = psi_leaf(e_leaf)
!   k        = -dE / dpsi_leaf
!
! solved on the same path the optimiser uses:
!   - xylem impairment with segments: the root, stem and leaf segments, with
!     the stem and leaf capped by the memory model (kmax / kmax_ref); this
!     also gives psi_stem
!   - xylem impairment alone: the impaired whole-plant curve
!   - otherwise: the whole-plant vulnerability curve (b, c)
!
! Used for the closed-stomata state here and for the cuticular loss in
! sf_stom, so both solve on the same path as the Ci samples.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
SUBROUTINE som_path_psi( pft, n_e_leaf, land_pts, open_pts, veg_index,       &
                         open_index, e_leaf, psi_root_zone, kmax_ref, kmax,    &
                         kcrit, conductance_b, conductance_c,                  &
                         psi_leaf_extreme, psi_root_extreme,                   &
                         l_xylem_impairment, psi_leaf, leaf_k, psi_stem,      &
                         l_exact )

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_impairment_mod, ONLY: leaf_psi_impaired_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_segments_jls,                 &
                                       leaf_psi_segments_lut_jls
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut,               &
                                l_som_plant_segments

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pts),                                                  &
      ! Transpiration (mol H2O m-2 s-1, the units of kmax * psi).
  psi_root_zone(land_pts),                                                     &
      ! Source water potential (Pa).
  kmax_ref(land_pts), kmax(land_pts), kcrit(land_pts),                         &
      ! Intact and impaired maximum conductance, and the critical
      ! conductance, as passed to stom_opt_mod.
  conductance_b(land_pts), conductance_c(land_pts),                            &
      ! Each point's (impaired) curve.
  psi_leaf_extreme(land_pts), psi_root_extreme(land_pts)
      ! Impairment memory, as passed to stom_opt_mod.

LOGICAL, INTENT(IN) :: l_xylem_impairment
      ! Impairment active for this PFT (not in spin-up), as in stom_opt_mod.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  psi_leaf(n_e_leaf, open_pts), leaf_k(n_e_leaf, open_pts)
      ! Leaf water potential (Pa) and whole-path conductance.

REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  psi_stem(n_e_leaf, open_pts)
      ! Stem water potential at its leaf-side end (Pa); with segments and
      ! impairment only, otherwise left unset.

LOGICAL, INTENT(IN), OPTIONAL :: l_exact
      ! .TRUE.: with segments, solve each segment exactly (bisection, see
      ! leaf_psi_segments_jls), also with som_psi_solver = 3.

REAL(KIND=real_jlslsm) :: kcap_frac(land_pts)
INTEGER :: j, l
LOGICAL :: l_ex

l_ex = .FALSE.
IF ( PRESENT(l_exact) ) l_ex = l_exact

IF ( l_xylem_impairment .AND. l_som_plant_segments ) THEN
  kcap_frac(:) = 1.0
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    kcap_frac(l) = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm))
  END DO

  IF ( l_ex ) THEN
    IF ( PRESENT(psi_stem) ) THEN
      CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pts,           &
                                  veg_index, open_index, e_leaf,               &
                                  psi_root_zone, kmax_ref, kcrit, psi_leaf,    &
                                  leaf_k, kcap_frac = kcap_frac,               &
                                  stem_psi = psi_stem, l_exact = .TRUE. )
    ELSE
      CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pts,           &
                                  veg_index, open_index, e_leaf,               &
                                  psi_root_zone, kmax_ref, kcrit, psi_leaf,    &
                                  leaf_k, kcap_frac = kcap_frac,               &
                                  l_exact = .TRUE. )
    END IF
  ELSE IF ( som_psi_solver == psi_solver_lut ) THEN
    IF ( PRESENT(psi_stem) ) THEN
      CALL leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pts,       &
                                      veg_index, open_index, e_leaf,           &
                                      psi_root_zone, kmax_ref, psi_leaf,       &
                                      leaf_k, kcap_frac = kcap_frac,           &
                                      stem_psi = psi_stem )
    ELSE
      CALL leaf_psi_segments_lut_jls( pft, n_e_leaf, land_pts, open_pts,       &
                                      veg_index, open_index, e_leaf,           &
                                      psi_root_zone, kmax_ref, psi_leaf,       &
                                      leaf_k, kcap_frac = kcap_frac )
    END IF
  ELSE
    IF ( PRESENT(psi_stem) ) THEN
      CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pts,           &
                                  veg_index, open_index, e_leaf,               &
                                  psi_root_zone, kmax_ref, kcrit, psi_leaf,    &
                                  leaf_k, kcap_frac = kcap_frac,               &
                                  stem_psi = psi_stem )
    ELSE
      CALL leaf_psi_segments_jls( pft, n_e_leaf, land_pts, open_pts,           &
                                  veg_index, open_index, e_leaf,               &
                                  psi_root_zone, kmax_ref, kcrit, psi_leaf,    &
                                  leaf_k, kcap_frac = kcap_frac )
    END IF
  END IF
ELSE IF ( l_xylem_impairment ) THEN
  CALL leaf_psi_impaired_jls( pft, n_e_leaf, land_pts, open_pts, veg_index,   &
                              open_index, e_leaf, psi_root_zone, kmax_ref,     &
                              kmax, kcrit, conductance_b, conductance_c,       &
                              psi_leaf_extreme, psi_root_extreme,              &
                              psi_leaf, leaf_k )
ELSE
  CALL leaf_psi_jls( pft, n_e_leaf, land_pts, open_pts, veg_index,            &
                     open_index, e_leaf, psi_root_zone, kmax, kcrit,           &
                     conductance_b, conductance_c,                             &
                   ! INTENT OUT
                     psi_leaf, leaf_k )
END IF

END SUBROUTINE som_path_psi

!-----------------------------------------------------------------------------
! Given a Ci sample grid already evaluated by stom_opt_mod_ci (al_sample,
! kl_sample, psi_sample), builds the psi/k physical-feasibility mask,
! computes profit under the profit_max_profit_model formulation, and
! returns the index of the maximum-profit sample per open point.
!
! Returns optimal_index(j) = 0 when no sample for open point j is feasible;
! callers should then use sample index 0 (the closed-stomata state) rather
! than treat this as an error - see the note in stom_opt_mod.
!
! Used by the flat search's single pass to select the profit-maximising
! sample, keeping the mask/profit/MAXLOC logic out of stom_opt_mod itself.
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_profit_max_select(                                        &
! IN
        land_pts, open_pts, open_index, veg_index, n_sample,                  &
        al_sample, kl_sample, psi_sample, gl_sample, psi_root_zone, kcrit,    &
        gl_max, rd,                                                           &
        kl_hc_sample, kl_hc_max, l_xylem_impairment,                          &
! OUT
        l_good_sample, carbon_gain, hydraulic_cost, profit, optimal_index     &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: l_som_gain_gross

INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                    &
                            ! Number of land points.
, open_pts                                                                    &
                            ! Number of open points to select for.
, open_index(land_pts)                                                        &
                            ! Index of the open points on veg_index.
, veg_index(land_pts)                                                         &
                            ! Index of vegetated points on the land grid.
, n_sample
                            ! Number of Ci samples (entry 0 is the closed
                            ! state).
                            !
                            ! MDK, 21/10/26

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  al_sample(0:n_sample, open_pts)                                             &
                            ! Net photosynthesis of each sample
                            ! (mol CO2 m-2 s-1).
, kl_sample(0:n_sample, open_pts)                                             &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
, psi_sample(0:n_sample, open_pts)                                            &
                            ! Leaf water potential of each sample (Pa).
, gl_sample(0:n_sample, open_pts)                                             &
                            ! Stomatal conductance for H2O of each sample
                            ! (m/s).
, psi_root_zone(land_pts)                                                     &
                            ! Root zone water potential (Pa).
, kcrit(land_pts)                                                             &
                            ! Critical conductance (mol m-2 s-1 Pa-1).
, gl_max(land_pts)                                                            &
                            ! Cap on gl (m/s); <= 0 disables it.
                            !
                            ! MDK, 21/10/26
, kl_hc_sample(0:n_sample, open_pts)                                          &
                            ! Conductance used for the hydraulic cost - see
                            ! stom_opt_mod_ci.
, kl_hc_max(land_pts)                                                          &
                            ! Reference conductance for the hydraulic cost
                            ! with l_xylem_impairment - see stom_opt_mod.
, rd(land_pts)
                            ! Dark respiration (mol CO2 m-2 s-1); added back
                            ! to An for the gross-A gain (l_som_gain_gross).

LOGICAL, INTENT(IN) :: l_xylem_impairment
                            ! See stom_opt_mod.

LOGICAL, INTENT(OUT) :: l_good_sample(n_sample, open_pts)
                            ! .TRUE. for a feasible sample: psi_leaf <=
                            ! psi_root_zone, k > kcrit and gl <= gl_max.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  carbon_gain(0:n_sample, open_pts)                                           &
                            ! A / max A of each sample (dimensionless).
, hydraulic_cost(0:n_sample, open_pts)                                        &
                            ! Loss of conductance of each sample, as a
                            ! fraction of max_k - kcrit (dimensionless).
, profit(0:n_sample, open_pts)
                            ! carbon_gain - hydraulic_cost.

INTEGER, INTENT(OUT) :: optimal_index(open_pts)
                            ! Sample with the largest profit; 0 (closed) if
                            ! none is feasible.
                            !
                            ! MDK, 21/10/26

INTEGER :: j, l

REAL(KIND=real_jlslsm) :: max_al(open_pts), max_kl(open_pts)
REAL(KIND=real_jlslsm) :: gain_off(open_pts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_PROFIT_MAX_SELECT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------------
!  Make mask of those samples that are physically reasonable.
! ----------------------------------------------------------------------------
l_good_sample(:,:) = .TRUE.
DO j = 1, open_pts
  l = veg_index(open_index(j))

  ! Require psi <= psi_root_zone
  l_good_sample(:,j) = l_good_sample(:,j)                                    &
          .AND. psi_sample(1:,j) <= psi_root_zone(l) + TINY(psi_root_zone(l))

  ! Require k > k_crit
  l_good_sample(:,j) = l_good_sample(:,j)                                    &
                       .AND. kl_sample(1:,j) > kcrit(l)

  ! Require gl <= gl_max (when the cap is enabled)
  IF (gl_max(l) > 0.0) l_good_sample(:,j) = l_good_sample(:,j)               &
                       .AND. gl_sample(1:,j) <= gl_max(l)
END DO

! ----------------------------------------------------------------------------
!  Calculate the profit for each leaf state.
!    CG = An/max(An)   (An + Rd with l_som_gain_gross)
!    HC = (ki_max - k)/(ki_max - k_crit)
!    Profit = CG-HC
! ----------------------------------------------------------------------------
! gain offset: 0 (net An) or Rd (gross A), see l_som_gain_gross
DO j = 1, open_pts
  l = veg_index(open_index(j))
  gain_off(j) = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)
  max_al(j) = MAXVAL(al_sample(1:,j) + gain_off(j), MASK = l_good_sample(:,j))
END DO

IF (l_xylem_impairment) THEN
  ! Set the maximum conductance as equal to the current (unimpaired) root
  ! zone conductance. JBaguley
  DO j = 1, open_pts
    max_kl(j) = kl_hc_max(veg_index(open_index(j)))
  END DO
ELSE
  max_kl = MAXVAL(kl_hc_sample(1:,:), MASK = l_good_sample, DIM = 1)
END IF

carbon_gain(:,:)    = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:)         = 0.0

DO j = 1, open_pts
  l = veg_index(open_index(j))

  ! CG = An/max(An) would be a 0/0 (or x/0) whenever every feasible sample's
  ! net photosynthesis is <= 0 - e.g. right at dawn/dusk, when apar has only
  ! just ticked above the closed-stomata threshold. There is no carbon
  ! benefit to opening stomata in that state regardless of hydraulic cost,
  ! so skip the division and fall back to the closed-stomata state (mask
  ! every sample out; carbon_gain/hydraulic_cost/profit are left at the 0.0
  ! they were initialised to, and MAXLOC below returns 0 for this point).
  IF (max_al(j) > 0.0) THEN
    ! CG = An/max(An)
    carbon_gain(:,j) = (al_sample(:,j) + gain_off(j)) / max_al(j)

    ! HC = (ki_max - k)/(ki_max - k_crit). Without impairment max_kl(j) is
    ! the largest k of the feasible samples, all of which have k > kcrit,
    ! so the denominator is positive. With impairment max_kl(j) is
    ! kl_hc_max, not taken from a sample, and the MAX(..., TINY) is what
    ! stops a division by zero.
    !
    ! MDK, 21/10/26
    hydraulic_cost(:,j) = (max_kl(j) - kl_hc_sample(:,j))                      &
                          / MAX(max_kl(j) - kcrit(l), TINY(1.0_real_jlslsm))

    ! Profit = CG-HC
    profit(:,j) = carbon_gain(:,j) - hydraulic_cost(:,j)
  ELSE
    l_good_sample(:,j) = .FALSE.
  END IF

  ! ---------------------------------------------------------------------
  !  Identify maximum profit, excluding the closed stomata state.
  !  Returns 0 if all samples are not physically reasonable.
  ! ---------------------------------------------------------------------
  optimal_index(j) = MAXLOC(profit(1:,j), DIM = 1, MASK = l_good_sample(:,j))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_profit_max_select

!-----------------------------------------------------------------------------
! Bounded profit-max Ci search (som_ci_search = 2), one open point at a time.
!
! gl rises with Ci (A rises and ca - Ci falls), so E rises and k(psi_leaf)
! falls: the feasible samples of the flat grid form a single range
! [ccp, ci_b]. For each point:
!   1. ccp gives max k (E = 0); if it is infeasible every Ci is. A is
!      monotone in Ci, so if A <= 0 at the top of the range every feasible
!      A is too and the stomata are closed, as for the flat search.
!   2. ci_b by an Illinois (bracketed regula falsi) root find: with
!      psi_solver_lut on the cap that gl_max and k > kcrit put on gl, which
!      needs only photosynthesis (edge_by_gl_cap); otherwise on the
!      feasibility margin of full evaluations (edge_by_margin). Max A is
!      A(ci_b): the normalisation the flat grid gets from its feasible
!      samples in the limit of a fine grid.
!   3. golden-section maximisation of profit on [ccp, ci_b], where every Ci
!      is feasible (no cliff), to (ci_b - ccp)/opt_resol or n_iter
!      iterations.
! The best point seen (including ccp and ci_b) is returned. Each evaluation
! is scalar (eval_ci), so points stop as soon as they have converged.
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_bounded_search(                                            &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,       &
        km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,             &
        gl_max, l_multilayer, r_bl, r_ca, g_oth, n_top, n_iter,                &
        kmax_ref, conductance_b, conductance_c, psi_leaf_extreme,              &
        psi_root_extreme, l_xylem_impairment, kl_hc_max,                       &
! OUT
        ci_g, al_g, gl_g, kl_g, psi_g, el_g, carbon_gain_g, hydraulic_cost_g,  &
        psi_stem_g, lambda_g,                                                  &
! IN (optional)
        lai_cost                                                               &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE jules_vegetation_mod, ONLY: som_profit_model, cmax_profit_model,          &
                                cgain_profit_model, sox_opt_profit_model,      &
                                cap_profit_model, pm2_profit_model,            &
                                kplant_model, kplant_model_ground

                                ! MDK, 22/10/26
USE jules_vegetation_mod, ONLY: l_som_gain_gross, photo_collatz,               &
                                photo_farquhar, photo_johnson, photo_model,    &
                                CW_conductance,                                &
                                SOX_conductance,                               &
                                som_psi_solver, psi_solver_lut,                &
                                l_som_plant_segments, l_som_nsl,               &
                                ximpair_stomatal_response, ximpair_hold_gs,    &
                                ximpair_hold_psi,                              &
                                l_som_rhizo_series, l_som_coupled_e
USE xylem_hydraulics_CW_jls_mod, ONLY: k_path_zero_flow, kplant_scale
                                       ! MDK, 21/10/26
USE jb_photo_mod, ONLY: jb_eta_scale
USE pftparm, ONLY: c3, alpha, pft_conductance_model, psi_nsl_onset, psi_nsl0,  &
                   fsmc_mod, cmax_a, cmax_b, kmax_pft, cgain_varpi,            &
                   conductance_b_pft, conductance_c_pft
USE jules_surface_mod, ONLY: fwe_c3, fwe_c4, beta1, beta2, ratio
USE planet_constants_mod, ONLY: repsilon
USE c_rmol, ONLY: rmol
USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls, xylem_conductance_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: supply_lut_psi, supply_lut_e_crit,      &
                                       supply_lut_f, leaf_psi_segments_jls,    &
                                       ksr_path, som_psi_in_min, xylem_f,      &
                                       leaf_psi_segments_lut_jls

                                       ! MDK, 18/10/26
USE xylem_impairment_mod, ONLY: leaf_psi_impaired_jls
                                       ! MDK, 21/10/26

INTEGER, INTENT(IN) ::                                                         &
  land_pts, pft, open_pts, open_index(land_pts), pft_photo_model,             &
  veg_index(land_pts), n_top, n_iter

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),         &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts),  &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), je(land_pts),   &
  t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts), ipar(land_pts),   &
  kmax(land_pts), kcrit(land_pts), gl_max(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_bl(land_pts), r_ca(land_pts), g_oth(land_pts)
                            ! Coupled transpiration (see stom_opt_mod).

! Xylem impairment (as stom_opt_mod_ci): the leaf state on the impaired path
! (kmax, conductance_b/c and the extremes; segments: stem/leaf capped at
! kmax / kmax_ref), the hydraulic cost on the intact path (kmax_ref, the PFT
! curve) in series with the soil, relative to kl_hc_max.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  kmax_ref(land_pts), conductance_b(land_pts), conductance_c(land_pts),       &
  psi_leaf_extreme(land_pts), psi_root_extreme(land_pts), kl_hc_max(land_pts)
LOGICAL, INTENT(IN) :: l_xylem_impairment

LOGICAL, INTENT(IN) :: l_multilayer

! The chosen state of each open point:
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ci_g(open_pts),                                                              &
                            ! Internal CO2 pressure (Pa).
  al_g(open_pts),                                                              &
                            ! Net photosynthesis (mol CO2 m-2 s-1).
  gl_g(open_pts),                                                              &
                            ! Stomatal conductance for H2O (m/s).
  kl_g(open_pts),                                                              &
                            ! Whole-path conductance -dE/dpsi_leaf
                            ! (mol m-2 s-1 Pa-1).
  psi_g(open_pts),                                                             &
                            ! Leaf water potential (Pa).
  el_g(open_pts),                                                              &
                            ! Transpiration (mol H2O m-2 s-1).
  carbon_gain_g(open_pts),                                                     &
                            ! Carbon gain (dimensionless; see
                            ! carbon_gain_out in stom_opt_mod).
  hydraulic_cost_g(open_pts)
                            ! Hydraulic cost (dimensionless; see
                            ! hydraulic_cost_out in stom_opt_mod).
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: psi_stem_g(open_pts)
                            ! Stem segment leaf-side water potential at the
                            ! chosen state (xylem impairment with segments,
                            ! for the per-segment stem damage); psi_root_zone
                            ! otherwise and for closed points.
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: lai_cost(land_pts)
                            ! Leaf area of the call (see stom_opt_mod).
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: lambda_g(open_pts)
                            ! Marginal water price dA/dE at the chosen state
                            ! (mol CO2 per mol H2O; see som_lambda_call), on
                            ! the search path; 0 if closed.
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm) :: el_s1(1,1), psi_s1(1,1), kl_s1(1,1), st_s1(1,1)

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  edge_rtol = 1.0e-3,                                                          &
                            ! ci_b to edge_rtol * (ca - ci_b), so gl at the
                            ! edge (~ 1/(ca - Ci)) is good to ~0.1% ...
  opt_resol = 2000.0,                                                          &
                            ! Optimum Ci to (ci_b - ccp)/opt_resol.
  edge_margin = 1.0e-3,                                                        &
                            ! ... or once k is within this fraction of kmax
                            ! above kcrit, or gl within this fraction of
                            ! gl_max (gl ~ 1/(ca - Ci), so near ca the edge
                            ! needs a margin test, not just a Ci tolerance).
  golden_ratio = 0.6180339887498949_real_jlslsm
INTEGER, PARAMETER :: max_edge_iter = 20
INTEGER, PARAMETER :: n_nsl_iter = 12
                            ! Maximum Illinois steps for the nonstomatal
                            ! limitation factor f (l_som_nsl).
REAL(KIND=real_jlslsm), PARAMETER :: nsl_tol = 1.0e-3
                            ! Tolerance on f - f(psi_leaf(f)) (l_som_nsl).
INTEGER, PARAMETER :: n_ecrit_iter = 60
                            ! ProfitMax2: maximum doubling plus bisection
                            ! steps for E_crit.
REAL(KIND=real_jlslsm), PARAMETER :: ecrit_rtol = 1.0e-4
                            ! ProfitMax2: relative tolerance on E_crit.
                            !
                            ! MDK, 22/10/26
INTEGER, PARAMETER :: n_sperry_iter = 40
                            ! Hold leaf psi: maximum bisection steps.
REAL(KIND=real_jlslsm), PARAMETER :: sperry_psi_tol = 100.0
                            ! Hold leaf psi: tolerance on psi_leaf (Pa).

INTEGER :: i, j, l, side, idx1(land_pts)
LOGICAL :: l_intact_search
                            ! Hold leaf psi (ximpair_hold_psi): .TRUE.
                            ! while the search runs on the intact path.
REAL(KIND=real_jlslsm) :: psi_star, hc_star, cg_star, s_lo, s_hi, s_mid
                            ! Hold leaf psi: target psi_leaf and the intact-path
                            ! gain/cost, and the Ci bisection bracket.
REAL(KIND=real_jlslsm) :: h_ci, h_al, h_gl, h_el, h_psi, h_kl
LOGICAL :: l_h_ok
                            ! Hold leaf psi: the feasible evaluation on the
                            ! impaired path with the largest Ci, and whether
                            ! there is one.
                            !
                            ! MDK, 21/10/26
REAL(KIND=real_jlslsm) :: kcap_pts(land_pts)
                            ! Impairment with segments: stem/leaf cap,
                            ! kmax / kmax_ref.
REAL(KIND=real_jlslsm) :: b_curve(land_pts), c_curve(land_pts)
                            ! The PFT's intact vulnerability curve on
                            ! land_pts.
                            !
                            ! MDK, 21/10/26
LOGICAL ::                                                                     &
  l_lut,                                                                       &
                            ! Leaf psi straight from the supply table.
  ok_u,                                                                        &
                            ! The latest evaluation (eval_ci) is feasible.
  l_edge,                                                                      &
                            ! edge_by_gl_cap found the edge.
  l_cmax, l_cgain, l_sox, l_cap, l_pm2,                                        &
                            ! The profit model is CMax, CGain, SOX_opt, CAP
                            ! or ProfitMax2.
  l_abs
                            ! Any of those five: profit in absolute carbon
                            ! units, on net A.
                            !
                            ! MDK, 22/10/26

REAL(KIND=real_jlslsm) :: lam_a1, lam_e1
                            ! lambda_g: A and E at the lower Ci.
REAL(KIND=real_jlslsm) :: e_crit
                            ! ProfitMax2: transpiration at which the
                            ! conductance of the path the search runs on
                            ! falls to kcrit (mol H2O m-2 s-1, the basis of
                            ! el_u).
                            !
                            ! MDK, 22/10/26
REAL(KIND=real_jlslsm) :: lai_call, k_unstressed, k_sox_ref
                            ! k_sox_ref: SOX's whole-path k at zero flow
                            ! (psi_bar = psi_root_zone), which normalises
                            ! the SOX factor (a constant for the point, so
                            ! it does not move the optimum).
                            ! lai_call: leaf area of this call (m2 leaf m-2
                            ! ground): kmax_ref / kmax_pft, divided by
                            ! kplant_scale with kplant_model = 1, or
                            ! lai_cost when given. The CMax and CGain costs,
                            ! given per unit leaf area, are scaled by it to
                            ! the basis of A (whole canopy, the sunlit or
                            ! shaded leaves, or a layer's leaves).
                            ! k_unstressed: CGain's k_max, the plant's
                            ! conductance at psi = 0 (plant only, no soil
                            ! link).
                            !
                            ! MDK, 21/10/26

REAL(KIND=real_jlslsm) ::                                                      &
  g_off, max_al, max_kl, ci_lo, ci_top, tol1,                                  &
  ! Evaluation state (eval_ci); al0_u is A without the nonstomatal limitation;
  ! khc_u the conductance for the hydraulic cost (kl_u unless impaired)
  last_ci, al_u, gl_u, el_u, psi_u, kl_u, al0_u, khc_u,                        &
  ! Edge (last feasible) and the infeasible end of its bracket
  e_ci, e_al, e_gl, e_el, e_psi, e_kl, e_g, b_ci, b_g, g_u, c_ci, t_al,        &
  e_khc,                                                                       &
  ! Best seen
  best_f, best_ci, best_al, best_gl, best_el, best_psi, best_kl, best_khc,     &
  ! Golden-section bracket
  a, b, d_ci, fc, fd

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_BOUNDED_SEARCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Leaf psi from the supply-function table directly where it applies;
! otherwise one-point calls of leaf_psi_jls (open index list idx1).
l_lut = ( pft_conductance_model(pft) == CW_conductance .OR.                    &
          pft_conductance_model(pft) == SOX_conductance ) .AND.                &
        som_psi_solver == psi_solver_lut .AND. .NOT. l_som_plant_segments
! With the soil resistance in series the root water potential varies with E, so
! the direct table lookup (psi_root_zone as the root water potential) and the
! table's e_crit in edge_by_gl_cap do not apply: use leaf_psi_jls, which adds
! the soil link, and the margin-based edge.
IF ( l_som_rhizo_series .AND. fsmc_mod(pft) == 2 ) l_lut = .FALSE.

! The table is the intact PFT curve: not for the impaired path. (With
! som_psi_solver = 3 and the memory model the impaired path is still solved
! from the table, capped, in impaired_state; but the direct lookup and the
! intact e_crit of edge_by_gl_cap do not apply, so the edge is found by
! edge_by_margin.)
!
! MDK, 18/10/26
IF ( l_xylem_impairment ) l_lut = .FALSE.
idx1(:) = 1

! CMax (stomata_model = 6): profit A_n - Theta(psi_leaf) in absolute carbon
! units instead of CG - HC (see the profit function below).
l_cmax = som_profit_model == cmax_profit_model

! CGain (stomata_model = 7): A_n - varpi (k_max - k)/k_max.
l_cgain = som_profit_model == cgain_profit_model

! SOX_opt (stomata_model = 5 with som_ci_search = 2): Eller et al. (2018)
! solved as Sabot et al. (2022) do,
!   A_n (k(psi_bar) - kcrit)/(k_ref - kcrit),
! k of the whole path at psi_bar, the mean of the root zone and leaf water
! potentials. The halfway point is Eller et al.'s (2018, Eqs 2.7-2.8, for
! the path-averaged conductance); Sabot et al. (2022, Eq. 16) use k at
! psi_leaf, which we do not. The kcrit subtraction is Sabot et al.'s
! (Eller: k/k_max).
l_sox = som_profit_model == sox_opt_profit_model

! CAP (stomata_model = 10): A_n with the nonstomatal limitation (l_som_nsl,
! set with it), no stomatal cost.
l_cap = som_profit_model == cap_profit_model

! ProfitMax2 (stomata_model = 11): A_n (1 - E / E_crit) (Wang et al. 2020).
!
! MDK, 22/10/26
l_pm2 = som_profit_model == pm2_profit_model
l_abs = l_cmax .OR. l_cgain .OR. l_sox .OR. l_cap .OR. l_pm2
kcap_pts(:) = 1.0
b_curve(:) = conductance_b_pft(pft)
c_curve(:) = conductance_c_pft(pft)

DO j = 1, open_pts
  l = veg_index(open_index(j))
  IF ( PRESENT(psi_stem_g) ) psi_stem_g(j) = psi_root_zone(l)
  IF ( PRESENT(lambda_g) ) lambda_g(j) = 0.0
  idx1(1) = open_index(j)
  g_off = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)

  ! (CMax, CGain, SOX and CAP work on net A, as Wolf et al. 2016 and Lu et
  ! al. 2020.)
  !
  ! MDK, 21/10/26
  IF ( l_abs ) g_off = 0.0
  lai_call = 1.0

  ! The leaf area of the call from the intact plant (kmax_ref; = kmax
  ! without xylem impairment), so damage does not shrink the CMax or CGain
  ! cost, and without the kplant_model = 1 ground factor, which is not leaf
  ! area.
  !
  ! MDK, 21/10/26
  IF ( l_abs ) lai_call = kmax_ref(l)                                          &
                          / MAX(kmax_pft(pft), TINY(1.0_real_jlslsm))
  IF ( l_abs .AND. kplant_model == kplant_model_ground )                       &
    lai_call = lai_call / kplant_scale(l)
  ! The two-leaf sunlit and shaded leaves give their leaf area (their kmax
  ! share is N-weighted, not leaf area).
  !
  ! MDK, 21/10/26
  IF ( l_abs .AND. PRESENT(lai_cost) ) lai_call = lai_cost(l)

  k_unstressed = 1.0
  k_sox_ref = 1.0

  IF ( l_sox ) k_sox_ref = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))

  ! (The plant's k_max: without the soil-to-root link, which would shrink
  ! the reference, and so the cost of a given embolism, as the soil dries.)
  IF ( l_cgain ) k_unstressed = MAX(k_path_zero_flow(pft, l, kmax(l),          &
                                    0.0_real_jlslsm, l_plant_only = .TRUE.),   &
                                    TINY(1.0_real_jlslsm))
  ci_lo = MAX(ccp(l), 0.0)

  ! With the gl_max cap the edge is found to the 1e-2 Pa floor on ca - Ci
  ! used for gl; without it (gl unbounded near ca) stop at the top of the
  ! flat grid, as that search does.
  IF ( gl_max(l) > 0.0 ) THEN
    ci_top = MAX(ca(l) - 1.0e-2_real_jlslsm, ci_lo)
  ELSE
    ci_top = ca(l) - (ca(l) - ci_lo) / n_top
  END IF

  !---------------------------------------------------------------------------
  ! 1. Lower end (max k) and top of the range.
  !---------------------------------------------------------------------------
  IF ( l_xylem_impairment .AND. l_som_plant_segments )                         &
    kcap_pts(l) = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm))
  ! Hold leaf psi: search on the intact path (as no impairment, at kmax_ref);
  ! the impaired path is used after the search (section 4).
  l_intact_search = l_xylem_impairment .AND.                                  &
                    ximpair_stomatal_response == ximpair_hold_psi
  CALL eval_ci(ci_lo)
  max_kl = kl_u

  ! Impairment, hold gs: the cost is relative to the intact path at zero
  ! flow (as the flat search). Hold leaf psi: the intact path's (kl_u, the
  ! search being on it).
  IF ( l_xylem_impairment .AND. ximpair_stomatal_response == ximpair_hold_gs ) &
    max_kl = kl_hc_max(l)
  IF ( .NOT. ok_u ) THEN
    CALL set_closed()
    CYCLE
  END IF

  CALL store_best(-HUGE(1.0_real_jlslsm), ci_lo)
  e_ci = ci_lo; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
  e_kl = kl_u; e_khc = khc_u; e_g = margin()

  CALL eval_ci(ci_top)

  ! (al0_u = al_u without the nonstomatal limitation, which can make A at
  ! ci_top <= 0 although the optimum inside the range has A > 0.)
  IF ( al0_u + g_off <= 0.0 ) THEN
    CALL set_closed()
    CYCLE
  END IF
  t_al = al_u

  !---------------------------------------------------------------------------
  ! 2. Upper edge of the feasible range, ci_b (e_*), unless the whole range
  !    is feasible.
  !---------------------------------------------------------------------------
  IF ( ok_u ) THEN
    e_ci = ci_top; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
    e_kl = kl_u; e_khc = khc_u
  ELSE
    l_edge = .FALSE.

    ! (edge_by_gl_cap finds the edge from A without the hydraulics, which
    ! does not hold once A depends on psi_leaf.)
    IF ( l_lut .AND. .NOT. l_som_nsl ) CALL edge_by_gl_cap()
    IF ( .NOT. l_edge ) CALL edge_by_margin()
  END IF

  ! Normalisation; no carbon benefit from opening => closed, as for the flat
  ! search (stom_opt_profit_max_select).
  max_al = e_al + g_off

  ! With the nonstomatal limitation, normalise by the unlimited A at the edge
  ! so the gain-to-cost weighting is as without it (and the limitation
  ! lowers the gain fraction).
  IF ( l_som_nsl ) max_al = photo_al(e_ci) + g_off
  IF ( max_al <= 0.0 ) THEN
    CALL set_closed()
    CYCLE
  END IF

  ! ProfitMax2: E_crit of the point (needs the edge state e_*; overwrites
  ! the *_u evaluation state, which is not used again before eval_ci).
  !
  ! MDK, 22/10/26
  e_crit = 1.0
  IF ( l_pm2 ) CALL find_e_crit()

  ! Both ends of the feasible range are candidates.
  best_f = profit(best_al, best_khc, best_psi, best_el)
  IF ( profit(e_al, e_khc, e_psi, e_el) > best_f ) THEN
    best_f = profit(e_al, e_khc, e_psi, e_el); best_ci = e_ci; best_al = e_al
    best_gl = e_gl; best_el = e_el; best_psi = e_psi; best_kl = e_kl
    best_khc = e_khc
  END IF

  !---------------------------------------------------------------------------
  ! 3. Golden-section search for the maximum profit on [ccp, ci_b], stopping
  !    after n_iter iterations or once the bracket is below
  !    (ci_b - ccp)/opt_resol. (Brent's parabolic steps were tried here, but
  !    in low light profit is nearly flat over a wide Ci range and, at this
  !    precision, too noisy for them.)
  !---------------------------------------------------------------------------
  ! The bracket [a, b] runs from the lowest Ci (ccp) to the feasible edge
  ! e_ci. tol1 is the target bracket width.
  a = ci_lo
  b = e_ci
  tol1 = (b - a) / opt_resol

  IF ( tol1 > 0.0 ) THEN
    ! Two interior points at the golden-ratio positions (0.382 and 0.618 of
    ! the bracket), c_ci < d_ci, and their profits fc and fd.
    c_ci = b - golden_ratio * (b - a)
    d_ci = a + golden_ratio * (b - a)
    CALL eval_ci(c_ci)
    fc = profit_u()
    CALL eval_ci(d_ci)
    fd = profit_u()

    ! Each step drops the part of the bracket beyond the worse point, which
    ! cannot hold the maximum (profit is unimodal in Ci). The kept interior
    ! point lands on a golden-ratio position of the new bracket, so only one
    ! new profit evaluation is needed per step, and the bracket shrinks by
    ! a factor 0.618 each time.
    DO i = 1, n_iter
      IF ( b - a <= tol1 ) EXIT

      IF ( fc >= fd ) THEN
        ! Maximum lies in [a, d_ci]: d_ci becomes the new upper end, c_ci the
        ! new upper interior point; evaluate a new lower interior point.
        b = d_ci
        d_ci = c_ci; fd = fc
        c_ci = b - golden_ratio * (b - a)
        CALL eval_ci(c_ci)
        fc = profit_u()
      ELSE
        ! Maximum lies in [c_ci, b]: c_ci becomes the new lower end, d_ci the
        ! new lower interior point; evaluate a new upper interior point.
        a = c_ci
        c_ci = d_ci; fc = fd
        d_ci = a + golden_ratio * (b - a)
        CALL eval_ci(d_ci)
        fd = profit_u()
      END IF
    END DO
  END IF

  ! Diagnostic dA/dE at the optimum: central difference over
  ! +-1e-3 of the range about best_ci, clipped to [ci_lo, e_ci] (one-sided
  ! at an end). Only *_u is changed, not the chosen state.
  !
  ! MDK, 22/10/26
  IF ( PRESENT(lambda_g) .AND. e_ci > ci_lo ) THEN
    s_lo = MAX(best_ci - 1.0e-3_real_jlslsm * (e_ci - ci_lo), ci_lo)
    s_hi = MIN(best_ci + 1.0e-3_real_jlslsm * (e_ci - ci_lo), e_ci)
    CALL eval_ci(s_lo)
    lam_a1 = al_u; lam_e1 = el_u
    CALL eval_ci(s_hi)
    IF ( ABS(el_u - lam_e1) > TINY(1.0_real_jlslsm) )                          &
      lambda_g(j) = (al_u - lam_a1) / (el_u - lam_e1)
  END IF

  ci_g(j) = MAX(0.0, best_ci)
  al_g(j) = MAX(-rd(l), best_al)
  gl_g(j) = MAX(0.0, best_gl)

  ! gl = 0 means closed: al = -rd, as in set_closed.
  IF (gl_g(j) <= 0.0) al_g(j) = -rd(l)
  psi_g(j) = best_psi
  el_g(j) = best_el
  kl_g(j) = best_kl
  carbon_gain_g(j) = (best_al + g_off) / max_al

  ! (Guarded as the flat search: with impairment max_kl is kl_hc_max, not
  ! taken from a sample. MDK, 21/10/26)
  hydraulic_cost_g(j) = (max_kl - best_khc)                                    &
                        / MAX(max_kl - kcrit(l), TINY(1.0_real_jlslsm))
  ! CMax, CGain: the cost as a fraction of the maximum A on the range.
  ! ProfitMax2: Theta = A_n E / E_crit, as a fraction of max A.
  !
  ! MDK, 22/10/26
  IF ( l_cmax .OR. l_cgain .OR. l_pm2 )                                        &
    hydraulic_cost_g(j) = (best_al - profit(best_al, best_khc, best_psi,       &
                                            best_el)) / max_al
  ! SOX: 1 - the factor on A (as the flat-grid SOX).
  IF ( l_sox ) hydraulic_cost_g(j) = 1.0 - sox_factor(best_psi)

  ! CAP: the nonstomatal loss of A as a fraction of the maximum A.
  IF ( l_cap ) hydraulic_cost_g(j) = (photo_al(best_ci) - best_al) / max_al

  !---------------------------------------------------------------------------
  ! 4. Hold leaf psi: the regulated psi_leaf is that of the intact optimum;
  !    the impaired path delivers less at it. Find Ci in [ci_lo, best_ci]
  !    whose transpiration on the impaired path gives psi_leaf = psi_star
  !    (psi falls as Ci, gl and E rise), by bisection. If even ci_lo is
  !    below psi_star, take ci_lo. The gain and cost are those of the
  !    intact-path optimum.
  !---------------------------------------------------------------------------
  IF ( l_intact_search .AND. gl_g(j) > 0.0 ) THEN
    ! Keep the intact optimum: its leaf psi is the target, and its gain and
    ! cost are reported.
    psi_star = best_psi
    cg_star  = carbon_gain_g(j)
    hc_star  = hydraulic_cost_g(j)

    ! From here on eval_ci solves on the impaired path. Reset the record of
    ! the best feasible evaluation (keep_feasible keeps the largest feasible
    ! Ci seen).
    l_intact_search = .FALSE.
    l_h_ok = .FALSE.
    h_ci = -HUGE(1.0_real_jlslsm)

    ! Bracket [s_lo, s_hi] = [ci_lo, intact optimum Ci]. At the top end the
    ! impaired path gives a lower psi_leaf than the intact one at the same
    ! Ci, so psi_u < psi_star is expected there.
    s_lo = ci_lo
    s_hi = best_ci
    CALL eval_ci(s_hi)
    CALL keep_feasible()

    IF ( psi_u < psi_star ) THEN
      ! Check that the target lies inside the bracket: at ci_lo (stomata
      ! nearly shut) psi_leaf must be above psi_star. If not, even ci_lo
      ! cannot hold psi_star and ci_lo is taken (below).
      CALL eval_ci(s_lo)
      CALL keep_feasible()

      IF ( psi_u > psi_star ) THEN
        ! Bisection on Ci for psi_leaf = psi_star: psi_leaf falls as Ci (and
        ! gl and E) rise, so move the end on the same side of psi_star as
        ! the midpoint. Stop after n_sperry_iter halvings or within
        ! sperry_psi_tol (Pa) of the target.
        DO i = 1, n_sperry_iter
          s_mid = 0.5 * (s_lo + s_hi)
          CALL eval_ci(s_mid)
          CALL keep_feasible()
          IF ( psi_u > psi_star ) THEN
            s_lo = s_mid
          ELSE
            s_hi = s_mid
          END IF

          IF ( ABS(psi_u - psi_star) < sperry_psi_tol ) EXIT
        END DO
      END IF
    END IF

    ! On the impaired path the last evaluation can be infeasible (k <=
    ! kcrit, or gl > gl_max). Take instead the feasible evaluation with the
    ! largest Ci (feasibility only falls as Ci rises), else ci_lo if that is
    ! feasible on the impaired path (section 1 checked it on the intact
    ! one), else close the point. (It returned the last state as open.)
    !
    ! MDK, 21/10/26
    IF ( .NOT. ok_u .AND. .NOT. l_h_ok ) THEN
      CALL eval_ci(ci_lo)
      CALL keep_feasible()
    END IF

    IF ( .NOT. ok_u .AND. l_h_ok ) THEN
      last_ci = h_ci; al_u = h_al; gl_u = h_gl; el_u = h_el
      psi_u = h_psi; kl_u = h_kl
    END IF

    IF ( ok_u .OR. l_h_ok ) THEN
      ci_g(j) = MAX(0.0, last_ci)
      al_g(j) = MAX(-rd(l), al_u)
      gl_g(j) = MAX(0.0, gl_u)
      IF (gl_g(j) <= 0.0) al_g(j) = -rd(l)
      psi_g(j) = psi_u
      el_g(j) = el_u
      kl_g(j) = kl_u
      carbon_gain_g(j) = cg_star
      hydraulic_cost_g(j) = hc_star
    ELSE
      CALL set_closed()
    END IF
  END IF

  ! Stem leaf-side psi at the chosen transpiration on the (impaired) segment
  ! chain, for the per-segment stem damage (diagnostic only: the state is
  ! already chosen).
  IF ( PRESENT(psi_stem_g) .AND. l_xylem_impairment .AND.                      &
       l_som_plant_segments .AND. gl_g(j) > 0.0 ) THEN
    el_s1(1,1) = el_g(j)

    ! MDK, 18/10/26
    IF ( som_psi_solver == psi_solver_lut ) THEN
      CALL leaf_psi_segments_lut_jls( pft, 1, land_pts, 1, veg_index, idx1,    &
                                      el_s1, psi_root_zone, kmax_ref, psi_s1,  &
                                      kl_s1, kcap_frac = kcap_pts,             &
                                      stem_psi = st_s1 )
    ELSE
    CALL leaf_psi_segments_jls( pft, 1, land_pts, 1, veg_index, idx1, el_s1,   &
                                psi_root_zone, kmax_ref, kcrit, psi_s1, kl_s1, &
                                kcap_frac = kcap_pts, stem_psi = st_s1 )
    END IF
    psi_stem_g(j) = st_s1(1,1)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

  !---------------------------------------------------------------------------
  ! Leaf state at internal CO2 ci for land point l: the same photosynthesis,
  ! conductance, transpiration and leaf-psi calculations as stom_opt_mod_ci
  ! for one sample, and its feasibility (ok_u) as stom_opt_profit_max_select.
  !---------------------------------------------------------------------------
  SUBROUTINE eval_ci(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: ag0, f_lo, f_hi, f_new, g_lo, g_hi, g_new
  INTEGER :: it, side

  last_ci = ci
  al0_u = photo_al(ci)
  CALL hydraulic_state(ci, al0_u)

  !-------------------------------------------------------------------------
  ! Nonstomatal limitation (l_som_nsl)
  !
  ! Gross photosynthesis is scaled by a factor of leaf water potential,
  !
  !   A = f(psi_leaf) A_unlimited
  !
  ! where psi_leaf comes from the transpiration that the scaled A gives.
  ! So f has to satisfy
  !
  !   f = f(psi_leaf(f))
  !
  ! For fixed Ci, f(psi_leaf(f)) does not increase with f, so
  !
  !   g(f) = f - f(psi_leaf(f))
  !
  ! increases with f and has a single root on [0, 1]:
  !   - g(1) > 0 here (the unlimited state, just evaluated);
  !   - at f = 0 there is no transpiration, so psi_leaf = psi_root_zone and
  !     g(0) = -f(psi_root_zone), which needs no hydraulics.
  !
  ! The root is found by the Illinois method (regula falsi) to |g| < nsl_tol,
  ! usually in 3-5 steps. Dark respiration Rd is not scaled. The factor f is
  ! nsl_factor (Dewar et al. 2022, Eqn 3(b), with an onset).
  !
  ! Reference:
  ! - Dewar, R., Holtta, T. and Salmon, Y. (2022), Exploring optimal
  !   stomatal control under alternative hypotheses for the regulation of
  !   plant sources and sinks. New Phytol, 233: 639-654.
  !   https://doi.org/10.1111/nph.17795
  !
  ! MDK, 23/10/26
  !-------------------------------------------------------------------------
  IF ( l_som_nsl ) THEN
    ! Nothing to do if the unlimited state is not limited (f = 1).
    IF ( nsl_factor(psi_u) < 1.0 ) THEN
      ! Unlimited gross photosynthesis, and the bracket for f:
      !   f_hi = 1: g_hi = 1 - f(psi_leaf at full A) > 0
      !   f_lo = 0: no transpiration, psi_leaf = psi_root_zone,
      !             g_lo = -f(psi_root_zone) <= 0
      ag0  = al0_u + rd(l)
      f_hi = 1.0
      g_hi = 1.0 - nsl_factor(psi_u)
      f_lo = 0.0
      g_lo = -nsl_factor(psi_root_zone(l))

      IF ( g_lo >= 0.0 ) THEN
        ! The soil alone stops photosynthesis (psi_root <= psi_nsl0).
        CALL hydraulic_state(ci, -rd(l))
      ELSE
        ! side records which end was replaced last (-1 low, +1 high), for
        ! the Illinois step.
        side = 0
        f_new = f_lo
        DO it = 1, n_nsl_iter
          ! Regula falsi: the root of the straight line through the two
          ! bracket ends, then the hydraulics and g at that f.
          f_new = f_lo - g_lo * (f_hi - f_lo) / (g_hi - g_lo)
          CALL hydraulic_state(ci, f_new * ag0 - rd(l))
          g_new = f_new - nsl_factor(psi_u)
          IF ( ABS(g_new) < nsl_tol ) EXIT

          ! Replace the end with the same sign as g_new, so the root stays
          ! bracketed. Illinois: if the same end is replaced twice in a row,
          ! halve g at the other end, so that end moves too and the bracket
          ! does not stall on one side.
          IF ( g_new < 0.0 ) THEN
            f_lo = f_new; g_lo = g_new
            IF ( side == -1 ) g_hi = 0.5 * g_hi
            side = -1
          ELSE
            f_hi = f_new; g_hi = g_new
            IF ( side == 1 ) g_lo = 0.5 * g_lo
            side = 1
          END IF
        END DO

        ! (The state left by hydraulic_state is that of f_new.)
      END IF
    END IF
  END IF

  ! A Ci sample is feasible when the leaf is no wetter than the root zone,
  ! the path conductance is above its critical value kcrit, and the leaf
  ! conductance is within gl_max (when gl_max is set).
  ok_u = psi_u <= psi_root_zone(l) + TINY(psi_root_zone(l))                    &
         .AND. kl_u > kcrit(l)                                                 &
         .AND. (gl_max(l) <= 0.0 .OR. gl_u <= gl_max(l))
  END SUBROUTINE eval_ci

  ! Conductance, transpiration, leaf psi and xylem k at internal CO2 ci for
  ! net photosynthesis al. Sets al_u, gl_u, el_u, psi_u, kl_u and khc_u, the
  ! conductance for the hydraulic cost (kl_u, or with impairment and
  ! ximpair_hold_gs the intact path's, khc1).
  !
  ! MDK, 21/10/26
  SUBROUTINE hydraulic_state(ci, al)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci, al
  REAL(KIND=real_jlslsm) :: vpd

  al_u = al
  gl_u = ratio * (al_u * rmol * t_leaf(l)) / MAX(ca(l) - ci, 1.0e-2_real_jlslsm)
  vpd = dq(l) * pstar(l) / repsilon
  el_u = MAX(0.0, vpd * gl_u / pstar(l) * pstar(l) / (rmol * t_leaf(l)))
  IF ( l_som_coupled_e ) el_u = vpd * g_eff_e(MAX(gl_u, 0.0_real_jlslsm),     &
                                    r_bl(l), r_ca(l), g_oth(l))                &
                                / (rmol * t_leaf(l))
  CALL psi_k_at_e()
  END SUBROUTINE hydraulic_state

  ! Leaf psi, conductance and the hydraulic-cost conductance (psi_u, kl_u,
  ! khc_u) at transpiration el_u, on the path the search runs on (split out
  ! of hydraulic_state for ProfitMax2's E_crit).
  !
  ! MDK, 22/10/26
  SUBROUTINE psi_k_at_e()
  REAL(KIND=real_jlslsm) :: el1(1,1), psi1(1,1), kl1(1,1), khc1(1,1)

  IF ( l_lut ) THEN
    psi_u = supply_lut_psi(pft, psi_root_zone(l),                              &
                           el_u / MAX(kmax(l), TINY(1.0_real_jlslsm)))
    kl_u = kmax(l) * supply_lut_f(pft, psi_u)
  ELSE IF ( l_xylem_impairment .AND. .NOT. l_intact_search ) THEN
    el1(1,1) = el_u
    CALL impaired_state( el1, psi1, kl1, khc1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  ELSE IF ( l_intact_search ) THEN
    ! Hold-psi search: the intact path (kmax_ref, PFT curve).
    el1(1,1) = el_u
    CALL leaf_psi_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax_ref, kcrit, b_curve, c_curve,       &
                       psi1, kl1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  ELSE
    el1(1,1) = el_u
    CALL leaf_psi_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax, kcrit, b_curve, c_curve,           &
                       psi1, kl1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  END IF
  khc_u = kl_u
  IF ( l_xylem_impairment .AND. .NOT. l_intact_search .AND.                  &
       ximpair_stomatal_response == ximpair_hold_gs ) khc_u = khc1(1,1)
  END SUBROUTINE psi_k_at_e

  !---------------------------------------------------------------------------
  ! ProfitMax2: find E_crit, the transpiration at which the conductance of
  ! the path the search runs on (psi_k_at_e) falls to kcrit,
  !
  !   kl(E_crit) = kcrit
  !
  ! Method:
  !   1. start from the edge transpiration e_el (feasible, kl > kcrit) and
  !      double E until kl <= kcrit, which brackets E_crit (the first upper
  !      guess is at least kmax x 0.1 MPa, in case e_el is ~0);
  !   2. bisect the bracket to a relative tolerance ecrit_rtol.
  !
  ! The upper edge of the Ci range is set by kcrit only when gl_max and
  ! ci_top do not bind, so e_el is a lower bound on E_crit, not E_crit
  ! itself. E_crit >= e_el > 0 is used. With no transpiration at the edge
  ! (no demand), a floor keeps E / E_crit finite (the cost is then ~0).
  !
  ! MDK, 23/10/26
  !---------------------------------------------------------------------------
  SUBROUTINE find_e_crit()
  REAL(KIND=real_jlslsm) :: e_lo, e_hi
  INTEGER :: it

  e_lo = MAX(e_el, 0.0_real_jlslsm)

  ! (kmax x 0.1 MPa: a scale for E when there is no flow at the edge.)
  e_hi = MAX(2.0 * e_lo, 1.0e5_real_jlslsm * kmax(l), TINY(1.0_real_jlslsm))

  ! 1. Bracket: double e_hi until the path conductance at that transpiration
  !    falls to kcrit. e_lo is always a transpiration with kl > kcrit.
  DO it = 1, n_ecrit_iter
    el_u = e_hi
    CALL psi_k_at_e()
    IF ( .NOT. ( kl_u > kcrit(l) ) ) EXIT
    e_lo = e_hi
    e_hi = 2.0 * e_hi
  END DO

  ! 2. Bisection: halve [e_lo, e_hi], keeping kl > kcrit at e_lo and
  !    kl <= kcrit at e_hi, until the bracket is within ecrit_rtol of e_hi.
  DO it = 1, n_ecrit_iter
    IF ( e_hi - e_lo <= ecrit_rtol * e_hi ) EXIT
    el_u = 0.5 * (e_lo + e_hi)
    CALL psi_k_at_e()
    IF ( kl_u > kcrit(l) ) THEN
      e_lo = el_u
    ELSE
      e_hi = el_u
    END IF
  END DO

  ! E_crit is the feasible end of the bracket, never below the edge
  ! transpiration, and above zero so that E / E_crit stays finite.
  e_crit = MAX(e_lo, e_el, TINY(1.0_real_jlslsm))
  END SUBROUTINE find_e_crit

  !---------------------------------------------------------------------------
  ! Xylem impairment, one transpiration el1 at point l (as stom_opt_mod_ci):
  ! leaf psi and conductance on the impaired path, and the conductance for
  ! the hydraulic cost, the intact path (kmax_ref, PFT curve) at the same
  ! water potentials in series with the soil (k_int(psi_leaf) /
  ! (1 + k_int(psi_in) / K_s), psi_in = psi_src - E / K_s; segments: the
  ! intact chain, leaf_k_intact).
  !---------------------------------------------------------------------------
  SUBROUTINE impaired_state( el1, psi1, kl1, khc1 )
  REAL(KIND=real_jlslsm), INTENT(IN)  :: el1(1,1)
  REAL(KIND=real_jlslsm), INTENT(OUT) :: psi1(1,1), kl1(1,1), khc1(1,1)
  REAL(KIND=real_jlslsm) :: kmx1(1), kcr1(1), b1(1), c1(1), psi_in1, k_s1

  IF ( l_som_plant_segments .AND. som_psi_solver == psi_solver_lut ) THEN
    ! The capped segments from the supply tables.
    !
    ! MDK, 18/10/26
    CALL leaf_psi_segments_lut_jls( pft, 1, land_pts, 1, veg_index, idx1, el1, &
                                    psi_root_zone, kmax_ref, psi1, kl1,        &
                                    kcap_frac = kcap_pts, leaf_k_intact = khc1 )
  ELSE IF ( l_som_plant_segments ) THEN
    CALL leaf_psi_segments_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,     &
                                psi_root_zone, kmax_ref, kcrit, psi1, kl1,     &
                                kcap_frac = kcap_pts, leaf_k_intact = khc1 )
  ELSE
    CALL leaf_psi_impaired_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,     &
                                psi_root_zone, kmax_ref, kmax, kcrit,          &
                                conductance_b, conductance_c,                  &
                                psi_leaf_extreme, psi_root_extreme,            &
                                psi1, kl1 )
    kmx1(1) = kmax_ref(l)
    kcr1(1) = kcrit(l)
    b1(1)   = conductance_b_pft(pft)
    c1(1)   = conductance_c_pft(pft)
    CALL xylem_conductance_jls( pft, 1, 1, psi1, kmx1, kcr1, b1, c1, khc1 )
    IF ( l_som_rhizo_series ) THEN
      k_s1 = ksr_path(l)
      IF ( k_s1 > TINY(1.0_real_jlslsm) ) THEN
        psi_in1 = MAX(psi_root_zone(l) - el1(1,1) / k_s1, som_psi_in_min)
        khc1(1,1) = khc1(1,1)                                                  &
                    / (1.0 + kmax_ref(l) * xylem_f(pft, psi_in1) / k_s1)
      ELSE IF ( k_s1 >= 0.0 .AND. el1(1,1) > 0.0 ) THEN
        khc1(1,1) = 0.0
      END IF
    END IF
  END IF
  END SUBROUTINE impaired_state

  ! Nonstomatal limitation factor at leaf psi (Pa): Dewar et al. (2022)
  ! Eqn 3(b) generalised to an onset, f = 1 - (psi - psi_on)/(psi0 - psi_on)
  ! clipped to [0, 1] (see l_som_nsl).
  !
  ! Reference:
  ! - Dewar, R., Holtta, T. and Salmon, Y. (2022), Exploring optimal
  !   stomatal control under alternative hypotheses for the regulation of
  !   plant sources and sinks. New Phytol, 233: 639-654.
  !   https://doi.org/10.1111/nph.17795
  !
  ! MDK, 21/10/26
  REAL(KIND=real_jlslsm) FUNCTION nsl_factor(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  nsl_factor = 1.0 - (psi - psi_nsl_onset(pft))                               &
                     / (psi_nsl0(pft) - psi_nsl_onset(pft))
  nsl_factor = MIN(MAX(nsl_factor, 0.0_real_jlslsm), 1.0_real_jlslsm)
  END FUNCTION nsl_factor

  ! Net photosynthesis at internal CO2 ci (as stom_opt_mod_ci).
  REAL(KIND=real_jlslsm) FUNCTION photo_al(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: wcarb, wlite, wexpt, wl, b1, b2, b3
  INTEGER :: errcode

  SELECT CASE ( pft_photo_model )
  CASE ( photo_collatz )
    IF (c3(pft) == 1) THEN
      wcarb = vcmax(l) * (ci - ccp(l)) / (ci + kc(l) * (1.0 + oi(l) / ko(l)))
      wlite = alpha(pft) * acr(l) * (ci - ccp(l)) / (ci + 2.0 * ccp(l))
      wlite = MAX(wlite, TINY(1.0e0))
      IF (l_multilayer) wlite = wlite / apar(l) * fapar_lf(l) * ipar(l)
      wexpt = fwe_c3 * vcmax(l)
    ELSE
      wcarb = vcmax(l)
      wlite = MAX(alpha(pft) * acr(l), TINY(1.0e0))
      IF (l_multilayer) wlite = wlite / apar(l) * fapar_lf(l) * ipar(l)
      wexpt = fwe_c4 * vcmax(l) * ci / pstar(l)
    END IF
    b1 = beta1
    b2 = -(wcarb + wlite)
    b3 = wcarb * wlite
    wl = -b2 / (2.0 * b1) - SQRT(b2 * b2 / (4 * beta1 * beta1) - b3 / b1)
    b1 = beta2
    b2 = -(wl + wexpt)
    b3 = wl * wexpt
    wl = -b2 / (2.0 * b1) - SQRT(b2 * b2 / (4 * beta2 * beta2) - b3 / b1)

  CASE ( photo_farquhar )
    wcarb = vcmax(l) * (ci - ccp(l)) / (ci + km(l))
    wlite = je(l) / 4.0 * (ci - ccp(l)) / (ci + 2.0 * ccp(l))
    IF (photo_model == photo_johnson) wlite = wlite * jb_eta_scale(ccp(l), ci)
    wlite = MAX(wlite, TINY(1.0e0))
    IF (l_multilayer) wlite = wlite * je_ratio(l)
    wl = MIN(wcarb, wlite)

  CASE DEFAULT
    wl = 0.0
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
                 'pft_photo_model should be photo_collatz or photo_farquhar')
  END SELECT

  photo_al = wl - rd(l)
  END FUNCTION photo_al

  !---------------------------------------------------------------------------
  ! Edge of the feasible range from the cap it puts on gl. With the supply
  ! table, k > kcrit <=> E < kmax * supply_lut_e_crit, and E = gl * vpd/(R T),
  ! so feasibility is gl <= g_cap = MIN(gl_max, E_crit * R T / vpd), i.e.
  !   f(Ci) = g_cap * (ca - Ci) - ratio * R T * A(Ci) >= 0,
  ! which is decreasing and close to linear in Ci and needs only A, not the
  ! hydraulics. Illinois root find for f = 0 between the feasible e_* end and
  ! ci_top, then one full evaluation at the edge. Leaves l_edge = .FALSE. (for
  ! edge_by_margin) if the table and the full evaluation disagree.
  !---------------------------------------------------------------------------
  SUBROUTINE edge_by_gl_cap()
  REAL(KIND=real_jlslsm) :: conv_e, g_cap, rk, fa, fb, fc_e, ea, eb, ec
  INTEGER :: it, sd

  conv_e = dq(l) * pstar(l) / repsilon / (rmol * t_leaf(l))
  g_cap = HUGE(1.0_real_jlslsm)
  IF ( gl_max(l) > 0.0 ) g_cap = gl_max(l)
  IF ( conv_e > 0.0 ) THEN
    ! (Coupled: the g whose transpiration is E_crit.)
    IF ( l_som_coupled_e ) THEN
      g_cap = MIN(g_cap, g_of_geff(kmax(l)                                     &
       * supply_lut_e_crit(pft, psi_root_zone(l),                              &
                           kcrit(l) / MAX(kmax(l), TINY(1.0_real_jlslsm)))     &
       / conv_e, r_bl(l), r_ca(l), g_oth(l)))
    ELSE
      g_cap = MIN(g_cap, kmax(l)                                               &
       * supply_lut_e_crit(pft, psi_root_zone(l),                              &
                           kcrit(l) / MAX(kmax(l), TINY(1.0_real_jlslsm)))     &
       / conv_e)
    END IF
  END IF

  IF ( g_cap >= HUGE(1.0_real_jlslsm) ) RETURN

  ! The leaf conductance a Ci implies is gl = ratio A R T_leaf / (ca - Ci),
  ! so gl <= g_cap is the same as
  !
  !   f(Ci) = g_cap (ca - Ci) - rk A(Ci) >= 0,   rk = ratio R T_leaf
  !
  ! which needs only the photosynthesis A(Ci), no hydraulics. (ca - Ci) is
  ! floored at 1e-2 Pa, as in gl.
  rk = ratio * rmol * t_leaf(l)

  ! Bracket: the current edge e_ci (feasible, f >= 0) and the top of the
  ! range ci_top (f < 0). Give up if they do not bracket a root.
  ea = e_ci
  fa = g_cap * MAX(ca(l) - ea, 1.0e-2_real_jlslsm) - rk * e_al
  eb = ci_top
  fb = g_cap * MAX(ca(l) - eb, 1.0e-2_real_jlslsm) - rk * t_al
  IF ( fa < 0.0 .OR. fb >= 0.0 ) RETURN

  ! Illinois root find for f = 0 (as in eval_ci): a regula falsi step,
  ! bisection if it falls outside the bracket, and the end with the same
  ! sign replaced; if the same end is replaced twice in a row, f at the
  ! other end is halved so the bracket keeps shrinking from both sides.
  ! sd records the end replaced last. Stop when the bracket is within
  ! edge_rtol of the distance to ca.
  sd = 0
  DO it = 1, max_edge_iter
    IF ( eb - ea <= edge_rtol * (ca(l) - ea) ) EXIT
    ec = (ea * fb - eb * fa) / (fb - fa)
    IF ( .NOT. ( ec > ea .AND. ec < eb ) ) ec = 0.5 * (ea + eb)
    fc_e = g_cap * MAX(ca(l) - ec, 1.0e-2_real_jlslsm) - rk * photo_al(ec)
    IF ( fc_e >= 0.0 ) THEN
      ea = ec; fa = fc_e
      IF ( sd == 1 ) fb = 0.5 * fb
      sd = 1
    ELSE
      eb = ec; fb = fc_e
      IF ( sd == -1 ) fa = 0.5 * fa
      sd = -1
    END IF
  END DO

  ! Full state at the edge (ea, the feasible end): a full evaluation with
  ! the hydraulics. If interpolation in the supply table leaves it just
  ! outside, step back towards lower Ci (up to 3 times). The edge is only
  ! moved if the new point lies above the current e_ci.
  DO it = 1, 3
    IF ( ea <= e_ci ) RETURN
    CALL eval_ci(ea)
    IF ( ok_u ) THEN
      e_ci = ea; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
      e_kl = kl_u; e_khc = khc_u
      l_edge = .TRUE.
      RETURN
    END IF
    ea = ea - edge_rtol * (ca(l) - ea)
  END DO
  END SUBROUTINE edge_by_gl_cap

  !---------------------------------------------------------------------------
  ! Edge of the feasible range by an Illinois root find on the feasibility
  ! margin of full evaluations (any hydraulics).
  !
  ! The margin (see margin()) is the smaller of the scaled distances to the
  ! two limits,
  !
  !   (kl - kcrit) / kmax   and   (gl_max - gl) / gl_max  (if gl_max set)
  !
  ! positive when the Ci is feasible, negative when not. The edge is where
  ! it crosses zero, between the feasible e_ci and the infeasible ci_top.
  !---------------------------------------------------------------------------
  SUBROUTINE edge_by_margin()
  INTEGER :: it

  ! edge_by_gl_cap may have left a different evaluation in *_u.
  IF ( l_lut ) CALL eval_ci(ci_top)

  ! Upper end of the bracket: ci_top (infeasible), and its margin.
  b_ci = ci_top
  b_g = margin()

  ! Illinois root find (as in edge_by_gl_cap). Each step is a full
  ! evaluation with the hydraulics. Stop when the bracket is within
  ! edge_rtol of the distance to ca, or once the feasible end is within
  ! edge_margin of a limit (close enough to the edge).
  side = 0
  DO it = 1, max_edge_iter
    IF ( b_ci - e_ci <= edge_rtol * (ca(l) - e_ci) .OR. e_g < edge_margin )    &
      EXIT

    ! Regula falsi step, or bisection if it falls outside the bracket.
    c_ci = (e_ci * b_g - b_ci * e_g) / (b_g - e_g)
    IF ( .NOT. ( c_ci > e_ci .AND. c_ci < b_ci ) ) c_ci = 0.5 * (e_ci + b_ci)
    CALL eval_ci(c_ci)
    g_u = margin()

    ! A feasible point becomes the new edge (with its full state); an
    ! infeasible one the new upper end. If the same end moves twice in a
    ! row, halve the margin at the other end (Illinois).
    IF ( ok_u ) THEN
      e_ci = c_ci; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
      e_kl = kl_u; e_khc = khc_u; e_g = g_u
      IF ( side == 1 ) b_g = 0.5 * b_g
      side = 1
    ELSE
      b_ci = c_ci; b_g = g_u
      IF ( side == -1 ) e_g = 0.5 * e_g
      side = -1
    END IF
  END DO
  END SUBROUTINE edge_by_margin

  ! Feasibility margin of the latest evaluation: > 0 where feasible,
  ! decreasing in Ci, for the Illinois edge search.
  REAL(KIND=real_jlslsm) FUNCTION margin()
  margin = (kl_u - kcrit(l)) / kmax(l)
  IF ( gl_max(l) > 0.0 ) margin = MIN(margin, (gl_max(l) - gl_u) / gl_max(l))
  IF ( ok_u ) THEN
    margin = MAX(margin, TINY(1.0_real_jlslsm))
  ELSE
    margin = MIN(margin, -TINY(1.0_real_jlslsm))
  END IF
  END FUNCTION margin

  ! Profit of a state with net photosynthesis al_in, hydraulic-cost
  ! conductance kl_in, leaf water potential psi_in and transpiration el_in.
  ! The search maximises it over Ci.
  !
  !   profit max   CG - HC = (A_n + g_off) / A_max
  !                          - (k_max - k) / (k_max - k_crit)
  !                (Sperry et al. 2017; as stom_opt_profit_max_select,
  !                dimensionless)
  !
  !   CMax         A_n - Theta(psi_leaf)
  !                (Wolf et al. 2016, Eq. 8)
  !
  !   CGain        A_n - varpi (k_max - k) / k_max
  !                (Lu et al. 2020)
  !
  !   SOX          A_n f_SOX(psi_leaf)
  !                (Eller et al. 2018)
  !
  !   CAP          A_n, which carries the nonstomatal limitation
  !
  !   ProfitMax2   A_n (1 - E / E_crit)
  !                (Wang et al. 2020, Eqn 11)
  !
  ! The last five are in mol CO2 m-2 s-1, on the basis of the call (whole
  ! canopy, the sunlit or shaded leaves, or a layer's leaves).
  !
  ! MDK, 23/10/26
  REAL(KIND=real_jlslsm) FUNCTION profit(al_in, kl_in, psi_in, el_in)
  REAL(KIND=real_jlslsm), INTENT(IN) :: al_in, kl_in, psi_in, el_in
  IF ( l_pm2 ) THEN
    ! ProfitMax2 (Wang et al. 2020): no cost parameter; the gain falls to
    ! zero as E reaches E_crit (find_e_crit).
    profit = al_in * (1.0 - el_in / e_crit)
  ELSE IF ( l_cmax ) THEN
    ! CMax (Wolf et al. 2016): a carbon cost Theta that rises as psi_leaf
    ! falls, in the same units as A.
    profit = al_in - theta_cmax(psi_in)
  ELSE IF ( l_cap ) THEN
    ! CAP: carbon gain only (al_in carries the nonstomatal limitation).
    profit = al_in
  ELSE IF ( l_sox ) THEN
    ! SOX_opt (Eller et al. 2018): A scaled by the hydraulic factor at
    ! psi_leaf.
    profit = al_in * sox_factor(psi_in)
  ELSE IF ( l_cgain ) THEN
    ! Lu et al. (2020) / Sabot et al. (2022, Eq. 15), varpi per unit leaf
    ! area (umol m-2 s-1) scaled to this call.
    profit = al_in - lai_call * 1.0e-6_real_jlslsm * cgain_varpi(pft)          &
                     * MAX(k_unstressed - kl_in, 0.0_real_jlslsm) / k_unstressed
  ELSE
    ! Profit max (Sperry et al. 2017): relative carbon gain minus relative
    ! hydraulic cost, both 0-1. g_off = Rd with l_som_gain_gross (gross
    ! gain), else 0; max_al and max_kl are the maxima over the Ci range.
    profit = (al_in + g_off) / max_al - (max_kl - kl_in) / (max_kl - kcrit(l))
  END IF
  END FUNCTION profit

  ! SOX factor (k(psi_bar) - kcrit)/(k_ref - kcrit), psi_bar the mean of
  ! the root zone and leaf water potentials, k of the whole path.
  REAL(KIND=real_jlslsm) FUNCTION sox_factor(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  sox_factor = ( k_path_zero_flow(pft, l, kmax(l),                            &
                                   0.5 * (psi + psi_root_zone(l)))             &
                 - kcrit(l) ) / MAX(k_sox_ref - kcrit(l), TINY(1.0_real_jlslsm))
  sox_factor = MAX(sox_factor, 0.0_real_jlslsm)
  END FUNCTION sox_factor

  ! CMax carbon cost of leaf water potential psi (Pa), Anderegg et al.
  ! (2018) / Sabot et al. (2022, Eq. 12): Theta = a/2 psi^2 + b |psi| with
  ! psi in MPa and a, b per unit leaf area in umol m-2 s-1 MPa-2 and MPa-1;
  ! returned in mol CO2 m-2 s-1 on the basis of this call (x lai_call).
  REAL(KIND=real_jlslsm) FUNCTION theta_cmax(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  REAL(KIND=real_jlslsm) :: p_mpa
  p_mpa = ABS(MIN(psi, 0.0_real_jlslsm)) * 1.0e-6_real_jlslsm
  theta_cmax = lai_call * 1.0e-6_real_jlslsm                                   &
               * ( 0.5_real_jlslsm * cmax_a(pft) * p_mpa * p_mpa               &
                   + cmax_b(pft) * p_mpa )
  END FUNCTION theta_cmax

  ! Profit of the latest evaluation (-HUGE if infeasible), kept if best.
  REAL(KIND=real_jlslsm) FUNCTION profit_u()
  IF ( ok_u ) THEN
    profit_u = profit(al_u, khc_u, psi_u, el_u)
  ELSE
    profit_u = -HUGE(1.0_real_jlslsm)
  END IF

  IF ( profit_u > best_f ) CALL store_best(profit_u, last_ci)
  END FUNCTION profit_u

  ! Keep the latest evaluation as the best so far, with profit f at Ci ci:
  ! copy the state the last eval_ci call left (al_u, gl_u, el_u, psi_u,
  ! kl_u, khc_u: net photosynthesis, leaf conductance, transpiration, leaf
  ! psi and the path conductances) into the best_* variables.
  !
  ! MDK, 23/10/26
  SUBROUTINE store_best(f, ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: f, ci
  best_f = f; best_ci = ci; best_al = al_u; best_gl = gl_u; best_el = el_u
  best_psi = psi_u; best_kl = kl_u; best_khc = khc_u
  END SUBROUTINE store_best

  !---------------------------------------------------------------------------
  ! Hold leaf psi (section 4): keep the latest evaluation if it is feasible
  ! and has the largest Ci so far.
  !
  ! MDK, 21/10/26
  !---------------------------------------------------------------------------
  SUBROUTINE keep_feasible()
  IF ( ok_u .AND. last_ci >= h_ci ) THEN
    l_h_ok = .TRUE.
    h_ci  = last_ci
    h_al  = al_u
    h_gl  = gl_u
    h_el  = el_u
    h_psi = psi_u
    h_kl  = kl_u
  END IF
  END SUBROUTINE keep_feasible

  SUBROUTINE set_closed()

  ! Closed stomata: zero flow, so the conductance is K(psi_root_zone) on the
  ! same curve as hydraulic_state, not 0 (0 made PLC_pft report 100 %; see
  ! stom_opt_mod_ci).
  REAL(KIND=real_jlslsm) :: el0(1,1), psi0(1,1), kl0(1,1)
  REAL(KIND=real_jlslsm) :: khc0(1,1)
  IF ( l_lut ) THEN
    kl0(1,1) = kmax(l) * supply_lut_f(pft, psi_root_zone(l))
  ELSE IF ( l_xylem_impairment ) THEN
    el0(1,1) = 0.0
    CALL impaired_state( el0, psi0, kl0, khc0 )
  ELSE
    el0(1,1) = 0.0
    CALL leaf_psi_jls( pft, 1, land_pts, 1, veg_index, idx1, el0,              &
                       psi_root_zone, kmax, kcrit, b_curve, c_curve,           &
                       psi0, kl0 )
  END IF
  ci_g(j) = ca(l)
  al_g(j) = -rd(l)
  gl_g(j) = 0.0
  kl_g(j) = kl0(1,1)
  psi_g(j) = psi_root_zone(l)
  el_g(j) = 0.0
  carbon_gain_g(j) = 0.0
  hydraulic_cost_g(j) = 0.0
  END SUBROUTINE set_closed

END SUBROUTINE stom_opt_bounded_search

!-----------------------------------------------------------------------------
! Supply-loss stomata (stomata_model = 9), one open point at a time: the
! regulation of Sperry et al. (2016, Eqn 5, with its saturation rule) on the
! profit-max hydraulics, with a carbon-coupled demand in place of their fixed
! G_max.
!
!   1. Demand. A fixed ci/ca (Medlyn et al. 2011 with the water price fixed,
!      so VPD acts only through the hydraulics):
!
!        ci_d = sl_cica_well_watered ca
!        A_d  = A(ci_d)
!        gl_d = ratio A_d R T / (ca - ci_d)
!        E'   = gl_d D                     (unregulated transpiration)
!
!      A_d <= 0 => closed.
!
!   2. Loss rule. On the supply function of the path (leaf_psi_jls: single
!      curve or segments, soil-to-root link, gravity in psi_root_zone),
!
!        P' = psi(E')
!        dP = (P0 - P') k(P') / k(P0),   k = -dE / dpsi_leaf
!
!      with P0 = psi_root_zone (zero flow). dP rises with E' to a maximum
!      and then falls; past the maximum it is held there (Sperry et al.
!      2016). A demand beyond E_crit (k(P') <= k_fail) is past the maximum.
!
!   3. Regulated state. The regulated E, gl and ci solve
!
!        psi(E) = P0 - dP
!        gl     = E / D  <= gl_d
!        A(ci) R T ratio = gl (ca - ci),   ci in [ccp, ci_d]
!
! k_fail is kcrit, or sl_kfail_frac * kmax. No som_gl_max: gl_max carries
! only the soil supply limit (l_som_supply_limit), applied to the regulated
! gl. Every search here is bracketed (bisection, Illinois or golden section
! after a coarse scan), so nothing can diverge near E_crit; the supply
! function itself comes from leaf_psi_jls (table, or Newton with
! som_psi_solver = 2).
!
! References:
! - Sperry, J.S. and Love, D.M. (2015), What plant hydraulics can tell us
!   about responses to climate-change droughts. New Phytol, 207: 14-27.
!   https://doi.org/10.1111/nph.13354
! - Sperry, J.S., Wang, Y., Wolfe, B.T., Mackay, D.S., Anderegg, W.R.L.,
!   McDowell, N.G. and Pockman, W.T. (2016), Pragmatic hydraulic theory
!   predicts stomatal responses to climatic water deficits. New Phytol,
!   212: 577-589. https://doi.org/10.1111/nph.14059
!
! MDK, 21/10/26
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_supply_loss(                                               &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,       &
        km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,             &
        gl_max, l_multilayer,                                                  &
! OUT
        ci_g, al_g, gl_g, kl_g, psi_g, el_g, carbon_gain_g, hydraulic_cost_g   &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE jules_vegetation_mod, ONLY: photo_collatz, photo_farquhar, photo_johnson,  &
                                photo_model, CW_conductance, SOX_conductance,  &
                                som_psi_solver, psi_solver_lut,                &
                                l_som_plant_segments, l_som_rhizo_series,      &
                                sl_kfail_frac
USE jb_photo_mod, ONLY: jb_eta_scale
USE pftparm, ONLY: c3, alpha, pft_conductance_model, fsmc_mod,                 &
                   sl_cica_well_watered, conductance_b_pft, conductance_c_pft
USE jules_surface_mod, ONLY: fwe_c3, fwe_c4, beta1, beta2, ratio
USE planet_constants_mod, ONLY: repsilon
USE c_rmol, ONLY: rmol
USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: supply_lut_psi, supply_lut_f

INTEGER, INTENT(IN) ::                                                         &
  land_pts, pft, open_pts, open_index(land_pts), pft_photo_model,              &
  veg_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),          &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts),   &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), je(land_pts),    &
  t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts), ipar(land_pts),    &
  kmax(land_pts), kcrit(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) :: gl_max(land_pts)
                            ! Cap on gl (m s-1) from sf_stom. With
                            ! stomata_model = 9 only the soil supply limit
                            ! (l_som_supply_limit; apply_supply_limit), as
                            ! som_gl_max is not used. <= 0: no cap.

LOGICAL, INTENT(IN) :: l_multilayer

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ci_g(open_pts), al_g(open_pts), gl_g(open_pts), kl_g(open_pts),              &
  psi_g(open_pts), el_g(open_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) :: carbon_gain_g(open_pts)
                            ! A / A_d, the fraction of the demand's net
                            ! photosynthesis kept.

REAL(KIND=real_jlslsm), INTENT(OUT) :: hydraulic_cost_g(open_pts)
                            ! 1 - k(P') / k(P0), the loss of conductance
                            ! the unregulated demand would cause.

REAL(KIND=real_jlslsm), PARAMETER :: de_slope = 1.0e-2
                            ! Relative step in E' for the slope of dP(E').

REAL(KIND=real_jlslsm), PARAMETER :: e_rtol = 1.0e-4
                            ! Regulated E to e_rtol * E' (or dP to e_rtol).

REAL(KIND=real_jlslsm), PARAMETER :: emax_rtol = 1.0e-3
                            ! E at the maximum of dP to emax_rtol * E'.

REAL(KIND=real_jlslsm), PARAMETER :: ci_rtol = 1.0e-4
                            ! ci to ci_rtol * (ci_d - ccp).

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  golden_ratio = 0.6180339887498949_real_jlslsm
                            ! (sqrt(5) - 1) / 2, for the golden-section search.

INTEGER, PARAMETER :: max_iter = 40
                            ! Maximum iterations of each root find and search.

INTEGER, PARAMETER :: n_scan = 8
                            ! Samples of dP(E) on [0, E_hi] before the golden
                            ! section, so that it brackets the largest one
                            ! (dP need not be unimodal with segments and the
                            ! soil link).

INTEGER :: j, l, it, side, idx1(land_pts)
INTEGER :: i_scan, i_best

LOGICAL :: l_lut, l_sat

REAL(KIND=real_jlslsm) :: b_curve(land_pts), c_curve(land_pts)
                            ! The PFT's vulnerability curve on land_pts, for
                            ! leaf_psi_jls (no xylem impairment here).

REAL(KIND=real_jlslsm) ::                                                      &
  rt, vpd, conv_e, k_fail, p0, k0, ci_d, al_d, gl_d, e_d, psi_d, k_d,          &
  dp, dp_d, psi_t, psi_u, k_u, e_reg, psi_reg, k_reg, gl_reg,                  &
  ! Illinois brackets
  xa, xb, xc, fa, fb, fc,                                                      &
  ! Golden section on dP(E), its upper end and the scan
  ga, gb, gc_e, gd_e, fgc, fgd, e_max, dp_max, e_hi, e_lo, f_scan,             &
  e_tol

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_SUPPLY_LOSS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Supply function from the table directly where it applies (as
! stom_opt_bounded_search); otherwise one-point calls of leaf_psi_jls.
l_lut = ( pft_conductance_model(pft) == CW_conductance .OR.                    &
          pft_conductance_model(pft) == SOX_conductance ) .AND.                &
        som_psi_solver == psi_solver_lut .AND. .NOT. l_som_plant_segments
IF ( l_som_rhizo_series .AND. fsmc_mod(pft) == 2 ) l_lut = .FALSE.
idx1(:) = 1
b_curve(:) = conductance_b_pft(pft)
c_curve(:) = conductance_c_pft(pft)

DO j = 1, open_pts
  l = veg_index(open_index(j))
  idx1(1) = open_index(j)
  rt = rmol * t_leaf(l)
  vpd = dq(l) * pstar(l) / repsilon

  ! E (mol m-2 s-1) = conv_e * gl (m s-1), as stom_opt_mod_ci.
  conv_e = vpd / rt
  p0 = psi_root_zone(l)
  k_fail = kcrit(l)
  IF ( sl_kfail_frac > 0.0 ) k_fail = sl_kfail_frac * kmax(l)

  !---------------------------------------------------------------------------
  ! 1. Demand.
  !---------------------------------------------------------------------------
  ci_d = sl_cica_well_watered(pft) * ca(l)
  al_d = photo_al(ci_d)
  IF ( al_d <= 0.0 .OR. ci_d <= MAX(ccp(l), 0.0) ) THEN
    CALL set_closed()
    CYCLE
  END IF
  gl_d = ratio * al_d * rt / MAX(ca(l) - ci_d, 1.0e-2_real_jlslsm)
  e_d = MAX(0.0, conv_e * gl_d)

  !---------------------------------------------------------------------------
  ! 2. Zero flow (P0, k0) and the demand on the supply function (P', k').
  !---------------------------------------------------------------------------
  CALL supply(0.0_real_jlslsm, psi_u, k0)
  IF ( k0 <= k_fail ) THEN
    ! The path has failed at zero flow: no supply.
    CALL set_closed()
    CYCLE
  END IF

  IF ( e_d <= TINY(1.0_real_jlslsm) ) THEN
    ! (Next to) no evaporative demand: nothing to regulate.
    CALL set_state(ci_d, al_d, gl_d, e_d, p0, k0)
    carbon_gain_g(j) = 1.0
    hydraulic_cost_g(j) = 0.0
    CYCLE
  END IF

  CALL supply(e_d, psi_d, k_d)
  hydraulic_cost_g(j) = MIN(MAX(1.0 - k_d / k0, 0.0), 1.0)

  ! dP at E' and whether E' is past the maximum of dP(E'): beyond E_crit,
  ! or dP falling there.
  l_sat = k_d <= k_fail
  dp_d = 0.0
  IF ( .NOT. l_sat ) THEN
    dp_d = dp_of(psi_d, k_d)
    CALL supply(e_d * (1.0 + de_slope), psi_u, k_u)
    l_sat = dp_of(psi_u, k_u) <= dp_d
  END IF

  IF ( l_sat ) THEN
    !-------------------------------------------------------------------------
    ! Maximum of dP(E) on [0, E_hi],
    !
    !   E_hi = MIN(E', E_crit)
    !
    ! dP starts at zero and rises at first,
    !
    !   dP(0) = 0,   d(dP)/dE = (1 + dP' dk/dE) / k0 = 1 / k0  at E = 0
    !
    ! then falls.
    !  a. If E' is beyond E_crit (k(E') <= k_fail), find E_crit by
    !     bisection on k(E), so the search and its tolerance scale with the
    !     usable range, not E'.
    !  b. n_scan samples bracket the largest dP (in case of more than one
    !     maximum), then golden section on that bracket to emax_rtol * E_hi.
    !     States beyond E_crit score -1.
    !-------------------------------------------------------------------------
    e_hi = e_d
    IF ( k_d <= k_fail ) THEN
      e_lo = 0.0
      DO it = 1, max_iter
        IF ( e_hi - e_lo <= 1.0e-2_real_jlslsm * emax_rtol * e_hi ) EXIT
        CALL supply(0.5 * (e_lo + e_hi), psi_u, k_u)
        IF ( k_u > k_fail ) THEN
          e_lo = 0.5 * (e_lo + e_hi)
        ELSE
          e_hi = 0.5 * (e_lo + e_hi)
        END IF
      END DO
      e_hi = e_lo
    END IF

    IF ( e_hi <= 0.0 ) THEN
      CALL set_closed()
      CYCLE
    END IF
    i_best = n_scan
    dp_max = -HUGE(1.0_real_jlslsm)
    DO i_scan = 1, n_scan
      f_scan = dp_scored(e_hi * i_scan / n_scan)
      IF ( f_scan > dp_max ) THEN
        dp_max = f_scan; i_best = i_scan
      END IF
    END DO
    e_max = e_hi * i_best / n_scan
    ga = e_hi * (i_best - 1) / n_scan
    gb = e_hi * MIN(i_best + 1, n_scan) / n_scan
    gc_e = gb - golden_ratio * (gb - ga)
    gd_e = ga + golden_ratio * (gb - ga)
    fgc = dp_scored(gc_e)
    fgd = dp_scored(gd_e)
    DO it = 1, max_iter
      IF ( gb - ga <= emax_rtol * e_hi ) EXIT
      IF ( fgc >= fgd ) THEN
        gb = gd_e
        gd_e = gc_e; fgd = fgc
        gc_e = gb - golden_ratio * (gb - ga)
        fgc = dp_scored(gc_e)
      ELSE
        ga = gc_e
        gc_e = gd_e; fgc = fgd
        gd_e = ga + golden_ratio * (gb - ga)
        fgd = dp_scored(gd_e)
      END IF
    END DO

    IF ( fgc > dp_max ) THEN
      e_max = gc_e; dp_max = fgc
    END IF

    IF ( fgd > dp_max ) THEN
      e_max = gd_e; dp_max = fgd
    END IF

    IF ( dp_max <= 0.0 ) THEN
      ! No point of the supply function below E' is usable.
      CALL set_closed()
      CYCLE
    END IF
    dp = dp_max
    xb = e_max
  ELSE
    dp = dp_d
    xb = e_d
  END IF

  !---------------------------------------------------------------------------
  ! 3. Regulated E: psi(E) = P0 - dP, by Illinois on [0, xb], where
  !    psi(0) - psi_t = dP >= 0 and psi(xb) - psi_t <= 0 (dP <= dP'(xb)).
  !---------------------------------------------------------------------------
  psi_t = p0 - dp
  xa = 0.0
  fa = dp
  e_tol = e_rtol * xb
  CALL supply(xb, psi_u, k_u)
  fb = psi_u - psi_t
  e_reg = xb; psi_reg = psi_u; k_reg = k_u
  IF ( fb < 0.0 .AND. dp > 0.0 ) THEN
    side = 0
    DO it = 1, max_iter
      xc = (xa * fb - xb * fa) / (fb - fa)
      IF ( .NOT. ( xc > xa .AND. xc < xb ) ) xc = 0.5 * (xa + xb)
      CALL supply(xc, psi_u, k_u)
      fc = psi_u - psi_t
      e_reg = xc; psi_reg = psi_u; k_reg = k_u
      IF ( ABS(fc) <= e_rtol * dp .OR. xb - xa <= e_tol ) EXIT
      IF ( fc > 0.0 ) THEN
        xa = xc; fa = fc
        IF ( side == 1 ) fb = 0.5 * fb
        side = 1
      ELSE
        xb = xc; fb = fc
        IF ( side == -1 ) fa = 0.5 * fa
        side = -1
      END IF
    END DO
  END IF

  gl_reg = MIN(gl_d, e_reg / MAX(conv_e, TINY(1.0_real_jlslsm)))

  ! Soil supply limit (l_som_supply_limit, via gl_max): cap gl and take E,
  ! psi_leaf and k at the capped transpiration. A cap at ~0 (no water to
  ! supply) closes, as the profit max's infeasible samples do.
  IF ( gl_max(l) > 0.0 .AND. gl_reg > gl_max(l) ) THEN
    IF ( gl_max(l) <= 10.0 * TINY(1.0_real_jlslsm) ) THEN
      CALL set_closed()
      CYCLE
    END IF
    gl_reg = gl_max(l)
    e_reg = conv_e * gl_reg
    CALL supply(e_reg, psi_reg, k_reg)
  END IF

  IF ( gl_reg <= 0.0 ) THEN
    CALL set_closed()
    CYCLE
  END IF

  IF ( gl_reg >= gl_d ) THEN
    CALL set_state(ci_d, al_d, gl_d, e_reg, psi_reg, k_reg)
    carbon_gain_g(j) = 1.0
    CYCLE
  END IF

  !---------------------------------------------------------------------------
  ! 4. ci for the regulated gl: h(ci) = gl (ca - ci) - ratio R T A(ci) = 0
  !    on [ccp, ci_d], h(ccp) = ratio R T rd > 0 and h(ci_d) =
  !    (gl - gl_d)(ca - ci_d) < 0, decreasing. Illinois.
  !---------------------------------------------------------------------------
  xa = MAX(ccp(l), 0.0)
  fa = gl_reg * MAX(ca(l) - xa, 1.0e-2_real_jlslsm) - ratio * rt * photo_al(xa)
  xb = ci_d
  fb = (gl_reg - gl_d) * MAX(ca(l) - ci_d, 1.0e-2_real_jlslsm)
  xc = xb
  IF ( fa > 0.0 .AND. fb < 0.0 ) THEN
    side = 0
    DO it = 1, max_iter
      IF ( xb - xa <= ci_rtol * (ci_d - MAX(ccp(l), 0.0)) ) EXIT
      xc = (xa * fb - xb * fa) / (fb - fa)
      IF ( .NOT. ( xc > xa .AND. xc < xb ) ) xc = 0.5 * (xa + xb)
      fc = gl_reg * MAX(ca(l) - xc, 1.0e-2_real_jlslsm)                        &
           - ratio * rt * photo_al(xc)
      IF ( fc > 0.0 ) THEN
        xa = xc; fa = fc
        IF ( side == 1 ) fb = 0.5 * fb
        side = 1
      ELSE
        xb = xc; fb = fc
        IF ( side == -1 ) fa = 0.5 * fa
        side = -1
      END IF
    END DO
    xc = 0.5 * (xa + xb)
  END IF

  ! A consistent with gl and ci (gl = ratio A R T/(ca - ci)).
  CALL set_state(xc,                                                           &
                 gl_reg * MAX(ca(l) - xc, 1.0e-2_real_jlslsm) / (ratio * rt),  &
                 gl_reg, e_reg, psi_reg, k_reg)
  carbon_gain_g(j) = al_g(j) / al_d
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

  ! Leaf psi and whole-path k = -dE/dpsi_leaf at transpiration e (the units
  ! of kmax times Pa), as stom_opt_bounded_search's hydraulic_state.
  SUBROUTINE supply(e, psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: e
  REAL(KIND=real_jlslsm), INTENT(OUT) :: psi, k
  REAL(KIND=real_jlslsm) :: el1(1,1), psi1(1,1), kl1(1,1)
  IF ( l_lut ) THEN
    psi = supply_lut_psi(pft, psi_root_zone(l),                                &
                         e / MAX(kmax(l), TINY(1.0_real_jlslsm)))
    k = kmax(l) * supply_lut_f(pft, psi)
  ELSE
    el1(1,1) = e
    CALL leaf_psi_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax, kcrit, b_curve, c_curve,           &
                       psi1, kl1 )
    psi = psi1(1,1)
    k = kl1(1,1)
  END IF
  END SUBROUTINE supply

  ! Regulated drop dP = (P0 - psi) k/k0 (Sperry et al. 2016, Eqn 5).
  REAL(KIND=real_jlslsm) FUNCTION dp_of(psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi, k
  dp_of = MAX(p0 - psi, 0.0) * k / k0
  END FUNCTION dp_of

  ! dP at e, or -1 beyond E_crit (k <= k_fail).
  REAL(KIND=real_jlslsm) FUNCTION dp_scored(e)
  REAL(KIND=real_jlslsm), INTENT(IN) :: e
  REAL(KIND=real_jlslsm) :: psi, k
  CALL supply(e, psi, k)
  IF ( k > k_fail ) THEN
    dp_scored = dp_of(psi, k)
  ELSE
    dp_scored = -1.0
  END IF
  END FUNCTION dp_scored

  SUBROUTINE set_state(ci, al, gl, e, psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci, al, gl, e, psi, k
  ci_g(j) = MAX(0.0, ci)
  al_g(j) = MAX(-rd(l), al)
  gl_g(j) = MAX(0.0, gl)
  el_g(j) = e
  psi_g(j) = psi
  kl_g(j) = k
  END SUBROUTINE set_state

  SUBROUTINE set_closed()

  ! Closed stomata: zero flow, so the conductance is k at P0 (see
  ! stom_opt_bounded_search).
  REAL(KIND=real_jlslsm) :: psi0, k0c
  CALL supply(0.0_real_jlslsm, psi0, k0c)
  ci_g(j) = ca(l)
  al_g(j) = -rd(l)
  gl_g(j) = 0.0
  kl_g(j) = k0c
  psi_g(j) = psi_root_zone(l)
  el_g(j) = 0.0
  carbon_gain_g(j) = 0.0
  hydraulic_cost_g(j) = 0.0
  END SUBROUTINE set_closed

  ! Net photosynthesis at internal CO2 ci (as stom_opt_bounded_search).
  REAL(KIND=real_jlslsm) FUNCTION photo_al(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: wcarb, wlite, wexpt, wl, b1, b2, b3
  INTEGER :: errcode

  SELECT CASE ( pft_photo_model )
  CASE ( photo_collatz )
    IF (c3(pft) == 1) THEN
      wcarb = vcmax(l) * (ci - ccp(l)) / (ci + kc(l) * (1.0 + oi(l) / ko(l)))
      wlite = alpha(pft) * acr(l) * (ci - ccp(l)) / (ci + 2.0 * ccp(l))
      wlite = MAX(wlite, TINY(1.0e0))
      IF (l_multilayer) wlite = wlite / apar(l) * fapar_lf(l) * ipar(l)
      wexpt = fwe_c3 * vcmax(l)
    ELSE
      wcarb = vcmax(l)
      wlite = MAX(alpha(pft) * acr(l), TINY(1.0e0))
      IF (l_multilayer) wlite = wlite / apar(l) * fapar_lf(l) * ipar(l)
      wexpt = fwe_c4 * vcmax(l) * ci / pstar(l)
    END IF
    b1 = beta1
    b2 = -(wcarb + wlite)
    b3 = wcarb * wlite
    wl = -b2 / (2.0 * b1) - SQRT(b2 * b2 / (4 * beta1 * beta1) - b3 / b1)
    b1 = beta2
    b2 = -(wl + wexpt)
    b3 = wl * wexpt
    wl = -b2 / (2.0 * b1) - SQRT(b2 * b2 / (4 * beta2 * beta2) - b3 / b1)

  CASE ( photo_farquhar )
    wcarb = vcmax(l) * (ci - ccp(l)) / (ci + km(l))
    wlite = je(l) / 4.0 * (ci - ccp(l)) / (ci + 2.0 * ccp(l))
    IF (photo_model == photo_johnson) wlite = wlite * jb_eta_scale(ccp(l), ci)
    wlite = MAX(wlite, TINY(1.0e0))
    IF (l_multilayer) wlite = wlite * je_ratio(l)
    wl = MIN(wcarb, wlite)

  CASE DEFAULT
    wl = 0.0
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
                 'pft_photo_model should be photo_collatz or photo_farquhar')
  END SELECT

  photo_al = wl - rd(l)
  END FUNCTION photo_al

END SUBROUTINE stom_opt_supply_loss

END MODULE stom_opt_jls_mod
