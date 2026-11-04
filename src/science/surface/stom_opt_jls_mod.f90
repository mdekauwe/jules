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

PRIVATE stom_opt_mod_ci, stom_opt_profit_max_select, stom_opt_bounded_search,  &
        stom_opt_supply_loss, g_eff_e, som_feasible
PRIVATE som_path_psi
PUBLIC stom_opt_mod, g_of_geff, e_of_g,                                        &
       set_leaf_temp_ci, clear_leaf_temp_ci

! Diagnostic marginal water price at the optimum, lambda_eff = dA/dE
! (mol CO2 per mol H2O; = umol per mmol x 1e-3), from the bounded Ci
! search (som_ci_search = 2: profit max, CMax, CGain, CAP, ProfitMax2), by
! a finite difference about the chosen Ci on the search path. It is 0 for
! closed stomata, for the other searches, and for open points taken by the
! well-watered fast path (l_som_skip_search_wellwatered). stom_opt_mod
! returns it per land point (optional argument lambda_eff);
! som_lambda_eff_pft: the canopy value per PFT, leaf-area weighted over the
! sunlit and shaded leaves by sf_stom (can_rad_mod = 7 only), for the
! output lambda_eff_pft.
!
! MDK, 24/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, PUBLIC :: som_lambda_eff_pft(:,:)

! Leaf temperature of each candidate Ci (l_leaf_temp_ci), for the next call
! of the bounded search; set by sf_stom with set_leaf_temp_ci before each
! leaf class's call and switched off after it (clear_leaf_temp_ci), as
! set_ksr_path. At point l, with E the class's transpiration (mol m-2 s-1,
! the basis of the search's E) and T0 the leaf temperature the call is
! given (t_leaf):
!
!   T_leaf(E) = ltc_a - ltc_b E                      (leaf energy balance)
!   dq(T)     = dq(T0) + ltc_s (T - T0)              (deficit, kg kg-1)
!   p(T)      = p(T0) exp(ltc_dl(p) (T - T0))        (photosynthesis
!               parameters: vcmax, je, kc, ko, ccp, km, rd)
!
! ltc_b = 0: the point keeps T0.
!
! MDK, 24/10/26
INTEGER, PARAMETER, PRIVATE :: n_ltc_par = 7
                            ! vcmax, je, kc, ko, ccp, km, rd.
LOGICAL, SAVE, PRIVATE :: l_ltc = .FALSE.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PRIVATE ::                          &
  ltc_a(:), ltc_b(:), ltc_s(:), ltc_dl(:,:)

! CO2 through the leaf boundary layer (l_som_coupled_e, or where r_bl > 0):
! gl_co2_bl gives this as the stomatal conductance of a Ci that the boundary
! layer alone cannot supply (the Ci is then infeasible).
!
! MDK, 04/11/26
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: gl_unreach = 1.0e3
                            ! (m s-1)
REAL(KIND=real_jlslsm), PARAMETER, PRIVATE :: r_bl_co2_fac = 1.37
                            ! Boundary-layer resistance to CO2 over that to
                            ! water vapour (the ratio of their
                            ! diffusivities to the 2/3 power).

CONTAINS

!-----------------------------------------------------------------------------
! Return the stomatal conductance to water vapour, gl (m s-1), that lets net
! photosynthesis a (mol CO2 m-2 s-1) run at internal CO2 ci, when the CO2
! crosses the leaf boundary layer and then the stomata, one after the other.
! The Ci searches use it to turn each candidate Ci into a stomatal
! conductance with coupled transpiration (l_som_coupled_e) or a leaf
! boundary layer resistance (r_bl > 0), where water vapour takes the same
! route out of the leaf.
!
! The CO2 supply through the two resistances in series is
!
!   A = (ca - ci) / (R T (1.6 / gl + 1.37 r_bl))
!
! where ca and ci are the CO2 partial pressures of the air and inside the
! leaf (Pa), R T is the gas constant times leaf temperature (J mol-1),
! r_bl is the boundary-layer resistance to water vapour (s m-1), 1.6 is the
! ratio of the diffusivities of water vapour and CO2 in air, and 1.37 turns
! r_bl into the boundary-layer resistance to CO2. Solving for gl,
!
!   1.6 / gl = (ca - ci) / (R T A) - 1.37 r_bl
!
! (Leuning et al. 1995, as in CABLE). The argument dc is ca - ci, which the
! caller keeps above a small floor, and rt is R T.
!
! If the boundary layer alone lets through less CO2 than a needs, no
! stomatal opening can supply it, and gl_unreach is returned, which the
! searches treat as infeasible. When a <= 0, or there is no boundary layer,
! the formula without the boundary layer is used,
!
!   gl = 1.6 a R T / (ca - ci)
!
! which gives gl <= 0 for a <= 0.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION gl_co2_bl( a, dc, rt, r_bl ) RESULT( gl )

USE jules_surface_mod, ONLY: ratio

IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) :: a, dc, rt, r_bl
REAL(KIND=real_jlslsm) :: gl, r_s

IF ( a <= 0.0 .OR. r_bl <= 0.0 ) THEN
  gl = ratio * a * rt / dc
ELSE
  r_s = dc / ( a * rt ) - r_bl_co2_fac * r_bl
  IF ( r_s > ratio / gl_unreach ) THEN
    gl = ratio / r_s
  ELSE
    gl = gl_unreach
  END IF
END IF

END FUNCTION gl_co2_bl

!-----------------------------------------------------------------------------
! Turn on a leaf temperature for each candidate Ci (l_leaf_temp_ci) for the
! next call of the bounded search, and store the coefficients sf_stom has
! worked out for each vegetated point:
!
!   T_leaf(E) = a - b E
!   dq(T)     = dq(T0) + s (T - T0)
!   p(T)      = p(T0) exp(dl (T - T0))
!
! The first is the leaf energy balance: the leaf temperature (K) falls as
! the leaf transpires E (mol H2O m-2 s-1). The second is the humidity
! deficit (kg kg-1) at leaf temperature T, with T0 the leaf temperature the
! call is given. The third scales each photosynthesis parameter p, with dl
! (K-1) in the order vcmax, je, kc, ko, ccp, km, rd.
!
! A negative b is set to 0, which keeps the point at T0. clear_leaf_temp_ci
! turns this off again after the call.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE set_leaf_temp_ci( land_pts, veg_pts, veg_index, a, b, s, dl )

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: a(land_pts), b(land_pts), s(land_pts),  &
                                      dl(land_pts,n_ltc_par)

INTEGER :: l, m

IF ( .NOT. ALLOCATED(ltc_a) ) ALLOCATE( ltc_a(land_pts), ltc_b(land_pts),     &
                                        ltc_s(land_pts),                       &
                                        ltc_dl(land_pts,n_ltc_par) )
ltc_a(:)    = 0.0
ltc_b(:)    = 0.0
ltc_s(:)    = 0.0
ltc_dl(:,:) = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  ltc_a(l)    = a(l)
  ltc_b(l)    = MAX(b(l), 0.0)
  ltc_s(l)    = s(l)
  ltc_dl(l,:) = dl(l,:)
END DO
l_ltc = .TRUE.

END SUBROUTINE set_leaf_temp_ci

!-----------------------------------------------------------------------------
! Turn off the leaf temperature of each candidate Ci (see set_leaf_temp_ci),
! so later calls of the bounded search use the leaf temperature they are
! given.
!-----------------------------------------------------------------------------
SUBROUTINE clear_leaf_temp_ci()
IMPLICIT NONE
l_ltc = .FALSE.
END SUBROUTINE clear_leaf_temp_ci

!-----------------------------------------------------------------------------
! Return the effective conductance for water vapour, g_eff (m s-1), from
! inside the leaf to the air at the first model level, so that the
! transpiration is E = D g_eff, with D the humidity deficit from the leaf
! to that air. Coupled transpiration (l_som_coupled_e) and a leaf boundary
! layer resistance (r_bl > 0) use it to work out the transpiration of each
! stomatal conductance the search tries.
!
! Water vapour leaves through the stomata (conductance g) and then the leaf
! boundary layer (resistance r_bl), one after the other, with combined
! conductance g'. It then crosses the canopy air to the first model level
! (resistance r_ca). The sunlit and shaded leaves share the canopy air, so
! the other group's conductance into it, g_o, also loads that resistance:
!
!   g'    = g / (1 + g r_bl)
!   g_eff = g' / (1 + r_ca (g' + g_o))
!
! MDK, 04/11/26
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
! Return the stomatal conductance g (m s-1) whose effective conductance from
! the leaf to the first-level air is x (m s-1), the inverse of g_eff_e. It
! turns a cap on transpiration into a cap on the stomatal conductance. Undoing
! the two steps of g_eff_e in turn,
!
!   g' = x (1 + r_ca g_o) / (1 - r_ca x)
!   g  = g' / (1 - g' r_bl)
!
! If either denominator is 0 or less, the resistances in series keep the
! effective conductance below x however wide the stomata open, and HUGE is
! returned.
!
! Used to turn a limit on transpiration into a limit on g: the hydraulic
! limit in the bounded search, and the soil supply limit in sf_stom.
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
! Stomatal optimisation: for each group of open leaves (the whole canopy,
! the sunlit or the shaded leaves, or the leaves of one canopy layer),
! choose how far the stomata open by weighing the carbon they let in
! against the harm the water they let out does to the plant's water
! transport.
!
! The search runs over the internal CO2 partial pressure, Ci (Pa), because
! photosynthesis and stomatal conductance both follow directly from it. For
! each candidate Ci:
!
!   1. A(Ci): net photosynthesis (mol CO2 m-2 s-1), from the Collatz or
!      Farquhar model.
!   2. gl(Ci) = 1.6 A R T / (ca - Ci): the stomatal conductance to water
!      vapour (m s-1) that lets that CO2 in, with ca the CO2 partial
!      pressure of the air (Pa), R T the gas constant times leaf
!      temperature, and 1.6 the ratio of the diffusivities of water vapour
!      and CO2 in air.
!   3. E(Ci): the transpiration that follows (mol H2O m-2 s-1), gl times
!      the humidity deficit from the leaf to the air.
!   4. psi_leaf(E): the leaf water potential (Pa, negative) needed to pull
!      E from the soil, from the plant's supply function (the integral of
!      its vulnerability curve from the soil to the leaf).
!   5. k(E) = -dE/dpsi_leaf: the conductance of the soil-to-leaf path left
!      at that E (mol m-2 s-1 Pa-1).
!   6. Profit max (stomata_model = 4; Sperry et al. 2017) scores the Ci as
!
!        profit = A / A_max - (k_max - k) / (k_max - k_crit)
!
!      the carbon gain as a fraction of the largest A on offer, minus the
!      fraction of the usable conductance lost. k_max is the conductance
!      of the path with no flow and k_crit the conductance at which the
!      path is taken to fail.
!
! A Ci is ruled out if k <= k_crit (hydraulic failure) or gl is above the
! cap gl_max (when set). The other stomatal models (CMax, CGain, SOX_opt,
! CAP and ProfitMax2) use steps 1 to 5 and replace the score in step 6
! (see the function profit in stom_opt_bounded_search). Supply-loss
! stomata (stomata_model = 9) use steps 1 to 5 and a fixed rule in place of
! a score (stom_opt_supply_loss).
!
! stom_opt_mod first sets every point to the closed-stomata state, then
! runs the chosen model on the points with open stomata. For profit max,
! points where even the widest opening costs almost no conductance skip
! the search (fast_path).
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
!
! MDK, 04/11/26
! *********************************************************************

SUBROUTINE stom_opt_mod (                                                      &
! IN
        land_pts, pft, open_pts, open_index,                                   &
        pft_photo_model, veg_pts, veg_index,                                   &
        ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
        km, dq, je, t_leaf, je_ratio, fapar_lf, kmax, kcrit,                   &
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
        psi_stem, lambda_eff,                                                  &
! IN (optional)
        lai_cost                                                               &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE jules_vegetation_mod, ONLY:                                                &
        som_n_sample,                                                          &
        profit_max_profit_model, sox_opt_profit_model, som_profit_model,       &
        cmax_profit_model, cgain_profit_model, supply_loss_model,              &
        cap_profit_model, som_ci_search, som_ci_bounded,                       &
        pm2_profit_model,                                                      &
        som_n_ci_golden_iter,                                                  &
        l_som_skip_search_wellwatered, som_hc_negligible_tol, l_som_nsl,       &
        xylem_impairment_none, l_som_plant_segments

USE pftparm, ONLY:                                                             &
        min_gl_pft, conductance_b_pft, conductance_c_pft,                      &
        pft_xylem_impairment_model

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_conductance_impaired_jls
USE plant_hydraulic_path_mod, ONLY: k_path_zero_flow, ksr_path,                &
                                    l_store_q_path, l_shared_stem_path,        &
                                    l_shared_whole_path, grav_node_drop
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
                            ! conductance_c_pft).
                            ! .FALSE. uses kmax/conductance_b/conductance_c
                            ! directly, with no impairment - kmax_ref and
                            ! the psi extremes are then not used.

!-----------------------------------------------------------------------------
! IN integer variables.
!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
!                           ! Number of land points to be processed
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
                            ! basis of gl), for l_som_coupled_e and the
                            ! multilayer canopy with l_leaf_temp. Where it
                            ! is > 0 the CO2 and water vapour cross it in
                            ! series with the stomata.
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
                            ! Specific humidity deficit from the leaf to the
                            ! air (kg H2O/kg air): to the level-1 air with
                            ! l_som_coupled_e or r_bl > 0.
,je(land_pts)                                                                  &
                            ! Electron transport rate (mol m-2 s-1)
,t_leaf(land_pts)                                                              &
                            ! Leaf temperature (K)
,je_ratio(land_pts)                                                            &
                            ! Electron transport of the sunlit (or shaded)
                            ! leaves of this call over that of the canopy
                            ! (dimensionless). Scales J in multilayer calls;
                            ! 1 otherwise.
,fapar_lf(land_pts)                                                            &
                            ! Fraction of the incident PAR absorbed by this
                            ! layer's sunlit (or shaded) leaves
                            ! (dimensionless). Used in multilayer calls.
,kmax(land_pts)                                                                &
                            ! Maximum whole-path (soil-to-leaf) conductance
                            ! (mol m-2 s-1 Pa-1), on the basis of this call
                            ! (whole canopy, the sunlit or shaded leaves, or
                            ! per unit leaf area of a layer's leaves). The
                            ! impaired value with l_xylem_impairment.
,kcrit(land_pts)                                                               &
                            ! Critical conductance (mol m-2 s-1 Pa-1, same
                            ! basis as kmax): states with k <= kcrit are
                            ! infeasible.
,gl_max(land_pts)                                                              &
                            ! Maximum stomatal conductance for H2O (m/s, same
                            ! basis as gl: per leaf area for multilayer, the
                            ! class total per ground area for the two-leaf
                            ! sunlit or shaded leaves, the canopy per ground
                            ! area for big leaf). Samples with gl above it
                            ! are infeasible. <= 0 disables the cap.
,ipar(land_pts)                                                                &
                            ! Incident PAR (W m-2).
,kmax_ref(land_pts)                                                            &
                            ! Unimpaired maximum conductance
                            ! (mol m-2 s-1 Pa-1), with the same scaling as
                            ! kmax (kmax is the impaired value when
                            ! l_xylem_impairment).
,conductance_b(land_pts)                                                       &
                            ! Conductance parameter b for each land point
                            ! (Pa).
,conductance_c(land_pts)                                                       &
                            ! Conductance parameter c for each land point
                            ! (dimensionless).
,psi_leaf_extreme(land_pts)                                                    &
                            ! Historic minimum leaf water potential (Pa).
,psi_root_extreme(land_pts)
                            ! Historic minimum root zone water potential (Pa).

!-----------------------------------------------------------------------------
! IN OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 rd(land_pts)
                            ! Dark respiration (mol CO2/m2/s). At open
                            ! points of the bounded search with
                            ! l_leaf_temp_ci it is returned at the chosen
                            ! candidate's leaf temperature, the Rd inside its
                            ! net photosynthesis al, so al + rd is its gross
                            ! photosynthesis; otherwise not modified (the
                            ! ozone damage is not coded here: flux_o3 and fo3
                            ! are returned as 0).

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
,fo3(land_pts)                                                                 &
                            ! Ozone exposure factor. Always 0: ozone not
                            ! coded.
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
,hydraulic_cost_out(land_pts)
                            ! Hydraulic cost of the chosen state, as a
                            ! fraction (dimensionless): of the k range for
                            ! profit max, of max A for CMax, CGain, CAP and
                            ! ProfitMax2, and 1 - the factor on A for SOX.
                            ! 0 where the stomata are closed.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem(land_pts)
                            ! Stem water potential (Pa): the leaf-side end of
                            ! the stem segment of the chosen sample, with
                            ! l_som_plant_segments and xylem impairment;
                            ! psi_root_zone otherwise.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 lambda_eff(land_pts)
                            ! Marginal water price dA/dE at the optimum (mol
                            ! CO2 per mol H2O; see the module header), from
                            ! the bounded Ci search; 0 otherwise. Under hold
                            ! leaf psi it is taken at the intact-path optimum,
                            ! before hold_leaf_psi moves Ci, so it is that
                            ! optimum's dA/dE.
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
 lai_cost(land_pts)
                            ! Leaf area of this call (m2 leaf m-2 ground):
                            ! the CMax / CGain cost is per unit leaf area.
                            ! Absent: kmax_ref / kmax_pft (the leaf area for
                            ! big leaf, without kplant_model = 1). The
                            ! two-leaf sunlit and shaded leaves pass theirs,
                            ! as their kmax share is N-weighted, and the
                            ! multilayer layers pass 1 (per unit leaf area).

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 optimal_index                                                                 &
                            ! Index of the optimal sample for a land point.
                            ! Used by the SOX flat search
                            ! (sox_opt_profit_model).
,i,j,l,m                                                                       &
                            ! Iterators
,errcode
                            ! Error code to pass to ereport.

INTEGER, ALLOCATABLE ::                                                        &
 optimal_index_flat(:)
                            ! Index of the best sample for the flat search,
                            ! from the shared stom_opt_profit_max_select
                            ! helper.

LOGICAL, ALLOCATABLE ::                                                        &
 l_good_sample(:,:)                                                            &
                            ! Feasibility mask for the single-pass
                            ! sox_opt_profit_model search.
,l_good_sample_flat(:,:)
                            ! Feasibility mask for the flat search.

!-----------------------------------------------------------------------------
! Arrays holding the Ci samples. Only needed for points with open stomata.
! (som_ci_search = 2 uses stom_opt_bounded_search instead, which needs none
! of the sample arrays.)
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
REAL(KIND=real_jlslsm), ALLOCATABLE ::                                         &
 ci_sample(:,:)                                           &
                            ! Internal CO2 pressure (Pa).
,al_sample(:,:)                                           &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(:,:)                                           &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(:,:)                                           &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
,kl_hc_sample(:,:)                                        &
                            ! Conductance used for the hydraulic cost
                            ! (mol m-2 s-1 Pa-1) - see stom_opt_mod_ci.
,kl_SOX(:,:)                                              &
                            ! Conductance for the SOX_opt model: the whole
                            ! path at the mean of the leaf and root zone
                            ! water potentials (mol m-2 s-1 Pa-1).
,psi_sample(:,:)                                          &
                            ! Water potential at leaf (Pa)
,SOX_mean_psi(:,:)                                        &
                            ! Mean water potential between leaf and root zone
                            ! (Pa). Used in SOX_opt model.
,el_sample(:,:)                                           &
                            ! Transpiration rate (mol H2O/m2/s)
,profit(:,:)                                              &
                            ! Profit for each leaf state
,carbon_gain(:,:)                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(:,:)                                      &
                            ! Hydraulic cost for each leaf state
,psi_stem_sample(:,:)
                            ! Stem water potential of each Ci sample (Pa).
! (The sample arrays are allocated only by the flat searches, to keep them
! off the stack.)
REAL(KIND=real_jlslsm) ::                                                      &
 k_top(open_pts)                                                               &
                            ! SOX flat search: the largest whole-path
                            ! conductance at psi_bar over the feasible
                            ! samples (mol m-2 s-1 Pa-1).
,kl_hc_max(land_pts)                                                           &
                            ! Reference conductance for the hydraulic cost:
                            ! the intact path with no flow at the root zone
                            ! water potential, with l_xylem_impairment, the
                            ! store's flow held fixed (l_store_q_path) or the
                            ! shared root and stem (l_shared_stem_path);
                            ! kmax_ref otherwise.
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
,lambda_bnd(open_pts)                                                          &
                            ! dA/dE at the optimum (see lambda_eff).
,rd_bnd(open_pts)
                            ! Dark respiration in al at the optimum (see
                            ! rd).

! -- Negligible-hydraulic-cost fast path (profit_max_profit_model only) --
! A single-point evaluation at the most water-demanding candidate Ci (near
! ca) - if even that point's hydraulic cost is negligible, the full Ci
! search is skipped for that point (see the fast-path block below).
REAL(KIND=real_jlslsm) ::                                                      &
 ci_e_fp(0:1, open_pts)                                                        &
,al_e_fp(0:1, open_pts)                                                        &
,gl_e_fp(0:1, open_pts)                                                        &
,kl_e_fp(0:1, open_pts)                                                        &
,kl_hc_e_fp(0:1, open_pts)                                                     &
,psi_e_fp(0:1, open_pts)                                                       &
,el_e_fp(0:1, open_pts)                                                        &
,psi_stem_e_fp(0:1, open_pts)                                                  &
                            ! Stem water potential of the fast-path sample
                            ! (Pa), for psi_stem at fast-path points.
,hc_at_maxstress(open_pts)

LOGICAL :: l_fastpath(open_pts)

INTEGER ::                                                                     &
 open_pts_search                                                               &
                            ! Number of open points still needing the full
                            ! Ci search after the fast-path check below
                            ! (== open_pts when l_som_skip_search_wellwatered
                            ! is .FALSE.).
,open_index_search(open_pts)                                                   &
                            ! Compressed subset of open_index holding only
                            ! the points not resolved by the fast path.
,open_pts_flat                                                                 &
                            ! Number of points to run the flat search over
                            ! (all of open_index_search, or 0 with
                            ! som_ci_search = 2).
,open_index_flat(open_pts)
                            ! Compressed subset of open_index for the flat
                            ! search (same convention as open_index_search).

REAL(KIND=real_jlslsm) :: k0_fp(open_pts)
                            ! Reference conductance for the fast path's
                            ! hydraulic cost, as the search's: the zero-
                            ! flow conductance of the path at the root zone
                            ! water potential (with the soil link), or
                            ! kl_hc_max with l_xylem_impairment, so the
                            ! cost does not jump between fast-path and
                            ! searched points.

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

l_xylem_impairment = l_xylem_impairment_in .AND. .NOT. is_spinup .AND.         &
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
! The leaf at zero flow: the root zone's water potential, less the gravity
! drop at the root-stem node where it is set (grav_node_drop).
psi_leaf(:)= psi_root_zone(:)
DO m = 1, veg_pts
  l = veg_index(m)
  psi_leaf(l) = psi_root_zone(l) - grav_node_drop(l)
END DO

! Closed stomata: the whole-path conductance at zero flow (xylem at
! psi_root_zone, in series with the soil-to-root ksr_path), the same
! definition as the optimiser's leaf_k for open points.
leaf_k(:)  = kmax
IF (l_xylem_impairment .AND. .NOT. l_som_plant_segments) THEN
  ! Impaired xylem at psi_root_zone (any impairment model), with the
  ! soil-to-root conductance in series.
  CALL leaf_conductance_impaired_jls( pft, land_pts, psi_root_zone, kmax_ref,  &
                                      kmax, kcrit, conductance_b,              &
                                      conductance_c, psi_leaf_extreme,         &
                                      kx_closed )
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l),            &
                                 kx = kx_closed(l))
  END DO
ELSE IF (l_xylem_impairment) THEN
  ! Segments (memory model): the stem and leaf segments capped at
  ! kmax / kmax_ref of the intact ones, as psi_leaf_segments_jls.
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax_ref(l), psi_root_zone(l),        &
                  kcap_frac = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm)))
  END DO
ELSE
  DO m = 1, veg_pts
    l = veg_index(m)
    leaf_k(l) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))
  END DO
END IF

IF (PRESENT(psi_stem)) psi_stem(:) = psi_leaf(:)

! Default outputs for points with closed stomata. Overwritten below for
! points with open stomata. Set before the early return, so that a call
! with no open points does not return the previous call's values.
carbon_gain_out(:) = 0.0
hydraulic_cost_out(:) = 0.0
IF ( PRESENT(lambda_eff) ) lambda_eff(:) = 0.0

! If there are no land points with open stomata then no calculation is
! needed (leaf_k holds the closed-stomata value set above).
IF(0 == open_pts) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! With xylem impairment the hydraulic cost is measured against the
! unimpaired PFT vulnerability curve, relative to its conductance at the
! current root zone water potential.
b_ref(:) = conductance_b_pft(pft)
c_ref(:) = conductance_c_pft(pft)
kl_hc_max(:) = kmax_ref(:)
IF (l_xylem_impairment) THEN
  IF (l_som_plant_segments) THEN
    ! Segments (root, stem, leaf): the reference conductance is that of the
    ! intact plant when it is not transpiring. With no transpiration there
    ! is no drop in water potential along the plant, so every segment sits
    ! at the root zone water potential, and the segments add in series:
    !
    !     1 / kl_hc_max = sum over segments of 1 / k_s(psi_root_zone)
    !
    ! where k_s is the segment's maximum conductance (kmax_ref times its
    ! share, seg_kmax_fac) reduced by its own vulnerability curve. This is
    ! the most the undamaged plant could conduct at today's soil water
    ! potential; the hydraulic cost of a stomatal opening is the loss of
    ! conductance from it. k_path_zero_flow does this sum, and adds the
    ! soil-to-root conductance in series (ksr_path): the soil is
    ! not damaged, so the stomata see it. Vegetated points only
    ! (psi_root_zone means nothing elsewhere).
    DO m = 1, veg_pts
      l = veg_index(m)
      kl_hc_max(l) = k_path_zero_flow(pft, l, kmax_ref(l), psi_root_zone(l))
    END DO
  ELSE
    CALL leaf_conductance_jls( pft, land_pts, psi_root_zone, kmax_ref, kcrit,  &
                               b_ref, c_ref, kl_hc_max )

    ! The soil-to-root conductance in series: the soil is not damaged, so
    ! the stomata see it, as in the unimpaired cost (where leaf_k is the
    ! whole soil-to-leaf path). Zero flow: the series sum.
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

! The store's flow held fixed (l_store_q_path, the capacitance store): at
! E near 0 the shifted chain carries -q_f, back towards the soil, a state
! the store's balance never allows. So the cost is measured against the
! unshifted path at zero flow (as with impairment), not the conductance
! at the low end of the Ci range. Without impairment the search's path
! (kmax); with impairment kl_hc_max is already that of the intact path.
IF ( l_store_q_path .AND. .NOT. l_xylem_impairment ) THEN
  DO m = 1, veg_pts
    l = veg_index(m)
    kl_hc_max(l) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))
  END DO
END IF

! The shared root and stem (l_shared_stem_path, l_som_shared_stem): the
! cost of either leaf class is the whole plant's (see psi_k_at_e), so the
! reference is the intact whole plant (shared root and stem, the sunlit
! and shaded leaves together, kmax_ref) at zero flow, in series with the
! soil link, with or without impairment.
IF ( l_shared_stem_path ) THEN
  l_shared_whole_path = .TRUE.
  DO m = 1, veg_pts
    l = veg_index(m)
    kl_hc_max(l) = k_path_zero_flow(pft, l, kmax_ref(l), psi_root_zone(l))
  END DO
  l_shared_whole_path = .FALSE.
END IF

! ----------------------------------------------------------------------------
!  Choose Ci by the profit model.
! ----------------------------------------------------------------------------

SELECT CASE (som_profit_model)

CASE (profit_max_profit_model)

  CALL fast_path()

  IF (som_ci_search == som_ci_bounded) THEN
    !-------------------------------------------------------------------
    ! Root find for the feasible-range edge, then golden-section within
    ! it - see stom_opt_bounded_search. Never needs the flat search.
    !-------------------------------------------------------------------
    open_pts_flat = 0
    IF (open_pts_search > 0) THEN
      CALL bounded_search( open_pts_search, open_index_search )
      CALL scatter_bnd( open_pts_search, open_index_search, .TRUE. )
    END IF
  ELSE
    ! Every point not resolved by the fast path goes to the flat search
    ! below.
    open_pts_flat = open_pts_search
    open_index_flat(1:open_pts_search) =                                       &
                                      open_index_search(1:open_pts_search)
  END IF

  IF (open_pts_flat > 0) CALL profit_max_flat()

CASE (cmax_profit_model, cgain_profit_model, cap_profit_model,                 &
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
  !-------------------------------------------------------------------
  ! psi_stem as the profit-max bounded path, so the per-segment stem
  ! damage (l_ximpair_seg_memory) follows the stem's own water
  ! potential, not the root zone's.
  CALL bounded_search( open_pts, open_index )
  CALL scatter_bnd( open_pts, open_index, .TRUE. )

CASE (sox_opt_profit_model)
  !-------------------------------------------------------------------
  ! sox_opt_profit_model: with som_ci_search = 2 the bounded search (as
  ! the profit max), otherwise a
  ! single-pass search at som_n_sample resolution over [ccp, ca]. SOX_opt
  ! cannot run with xylem impairment (rejected at start-up, see
  ! pftparm_io), so neither search has an impaired branch.
  !-------------------------------------------------------------------
  IF (som_ci_search == som_ci_bounded) THEN
    CALL bounded_search( open_pts, open_index )
    CALL scatter_bnd( open_pts, open_index, .FALSE. )
  ELSE
    CALL sox_flat_search()
  END IF  ! som_ci_search (SOX)

CASE (supply_loss_model)
  !-------------------------------------------------------------------
  ! Supply-loss stomata (stomata_model = 9; Sperry and Love 2015): no
  ! optimisation. gl_max carries only the soil supply limit; see
  ! stom_opt_supply_loss. Its results come back in the bounded search's
  ! *_bnd arrays.
  !-------------------------------------------------------------------
  CALL stom_opt_supply_loss(                                                   &
  ! IN
      land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,         &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,         &
      km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,               &
      gl_max, l_multilayer,                                                    &
  ! OUT
      ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                                 &
      el_bnd, carbon_gain_bnd, hydraulic_cost_bnd                              &
          )

  CALL scatter_bnd( open_pts, open_index, .FALSE. )

CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
             'profit_model should be profit_max, SOX_opt, CMax, CGain, ' //    &
             'supply-loss, CAP or ProfitMax2')

END SELECT ! som_profit_model

! Open points the optimiser closed (gl = 0) get the same minimum
! conductance as points closed by darkness (gl = min_gl_pft above), so a
! night-closed and a drought-closed canopy are treated alike.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  IF ( gl(l) <= 0.0 ) gl(l) = min_gl_pft(pft)
END DO

! leaf_k is Kplant: the whole-path (soil-to-leaf) conductance of the chosen
! state, also with xylem impairment.
! Ozone exposure (flux_o3, fo3) is not calculated (TODO: port from
! leaf_jls_mod).

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

!-----------------------------------------------------------------------------
! Run the bounded Ci search (stom_opt_bounded_search) on n open points and
! put its results in the *_bnd arrays. idx holds, like open_index, the
! points' positions in the list of vegetated points (veg_index). Every
! model that uses the bounded search shares this one call: only CMax and
! CGain read lai_cost, and the stem water potential (psi_stem_bnd) is set
! only with xylem impairment on the segmented path.
!-----------------------------------------------------------------------------
SUBROUTINE bounded_search( n, idx )
INTEGER, INTENT(IN) :: n, idx(:)

CALL stom_opt_bounded_search(                                                  &
! IN
    land_pts, pft, n, idx,                                                     &
    pft_photo_model, veg_index,                                                &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
    km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,                 &
    gl_max, l_multilayer, r_bl, r_ca, g_oth,                                   &
    som_n_sample, som_n_ci_golden_iter,                                        &
    kmax_ref, conductance_b, conductance_c, psi_leaf_extreme,                  &
    psi_root_extreme, l_xylem_impairment, kl_hc_max,                           &
! OUT
    ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                                   &
    el_bnd, carbon_gain_bnd, hydraulic_cost_bnd,                               &
    psi_stem_g = psi_stem_bnd, lambda_g = lambda_bnd, rd_g = rd_bnd,           &
! IN (optional)
    lai_cost = lai_cost                                                        &
        )

END SUBROUTINE bounded_search

!-----------------------------------------------------------------------------
! Copy the bounded search's results (*_bnd) for n open points (idx as in
! bounded_search) to the outputs on the land-point grid. l_stem: also copy
! the stem water potential, which stom_opt_supply_loss does not set and
! SOX_opt does not use. The marginal water price (lambda_eff) is copied
! for every model except supply-loss stomata, and so is the dark respiration
! of the chosen state (rd).
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

  ! (Not set by the supply-loss rule.)
  IF ( PRESENT(lambda_eff) .AND. som_profit_model /= supply_loss_model )       &
    lambda_eff(l) = lambda_bnd(j)
  IF ( som_profit_model /= supply_loss_model ) rd(l) = rd_bnd(j)
END DO

END SUBROUTINE scatter_bnd

!-----------------------------------------------------------------------------
! Profit max: skip the search at points where water costs the plant almost
! nothing.
!
! The hydraulic cost only grows as Ci rises (wider stomata, more
! transpiration, lower leaf water potential, less conductance), so it is
! largest at the top of the Ci range, just below ca. With
! l_som_skip_search_wellwatered, one Ci just below ca is evaluated at each
! open point. If that state is feasible and its cost,
!
!   (k_ref - k) / (k_ref - k_crit)
!
! is below som_hc_negligible_tol, the profit is in effect the carbon gain
! alone, which rises with Ci, so the search would end at this same Ci and
! the point takes this state directly. k_ref is the conductance of the
! path with no flow (of the intact path with xylem impairment). The
! shortcut is not used with the nonstomatal limitation (l_som_nsl), which
! can put the best Ci lower. The other open points are listed in
! open_index_search for the full search.
!-----------------------------------------------------------------------------
SUBROUTINE fast_path()

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
  ci_lo(:) = ca(:) - MAX((ca(:) - MAX(ccp(:), 0.0)) * 0.001_real_jlslsm,       &
                         1.0e-3_real_jlslsm)
  ci_hi(:) = ca(:)

  CALL stom_opt_mod_ci(                                                        &
  ! IN
      land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,         &
      1, ci_lo, ci_hi,                                                         &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,         &
      km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,               &
      l_multilayer, r_bl, r_ca, g_oth,                                         &
      kmax_ref, conductance_b, conductance_c,                                  &
      psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                  &
  ! OUT
      ci_e_fp, al_e_fp, gl_e_fp, kl_e_fp, psi_e_fp, el_e_fp, kl_hc_e_fp,       &
      psi_stem_sample = psi_stem_e_fp                                          &
          )
END IF

open_pts_search = 0
DO j = 1, open_pts
  l = veg_index(open_index(j))
  l_fastpath(j) = .FALSE.
  IF (l_xylem_impairment) THEN
    k0_fp(j) = kl_hc_max(l)
  ELSE
    k0_fp(j) = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))
  END IF

  ! (Not with l_som_nsl: the limitation makes the optimum interior even
  ! where the hydraulic cost is negligible.)
  IF (l_som_skip_search_wellwatered .AND. .NOT. l_som_nsl .AND.                &
      psi_e_fp(1,j) <= psi_root_zone(l) + TINY(psi_root_zone(l)) .AND.         &
      kl_e_fp(1,j) > kcrit(l) .AND.                                            &
      (gl_max(l) <= 0.0 .OR. gl_e_fp(1,j) <= gl_max(l))) THEN
    hc_at_maxstress(j) = (k0_fp(j) - kl_hc_e_fp(1,j)) /                       &
                          MAX(k0_fp(j) - kcrit(l), TINY(1.0_real_jlslsm))
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

    ! The stem water potential of the sample, as the searches set it, so
    ! the per-segment stem damage follows the stem, not the root zone.
    IF (PRESENT(psi_stem)) psi_stem(l) = psi_stem_e_fp(1,j)

    ! The near-ca sample is the max-A one, so A / max A = 1.
    carbon_gain_out(l) = 1.0_real_jlslsm
    hydraulic_cost_out(l) = hc_at_maxstress(j)
  ELSE
    open_pts_search = open_pts_search + 1

    ! open_index_search holds positions in veg_index, as open_index does.
    open_index_search(open_pts_search) = open_index(j)
  END IF
END DO

END SUBROUTINE fast_path

!-----------------------------------------------------------------------------
! Choose the stomatal opening by profit max on an even grid of som_n_sample Ci
! values from ccp (the CO2 compensation point without dark respiration, Pa) to
! ca, for the points in open_index_flat: evaluate every Ci (stom_opt_mod_ci),
! pick the feasible one with the largest profit (stom_opt_profit_max_select)
! and copy it to the outputs.
!-----------------------------------------------------------------------------
SUBROUTINE profit_max_flat()

CALL alloc_flat()

!---------------------------------------------------------------------
! Flat grid: som_n_sample points over the full [ccp, ca] range, routed
! through the shared stom_opt_profit_max_select helper, run over
! open_index_flat.
!---------------------------------------------------------------------
ci_lo(:) = MAX(ccp(:), 0.0)
ci_hi(:) = ca(:)

CALL stom_opt_mod_ci(                                                          &
! IN
    land_pts, pft, open_pts_flat, open_index_flat,                             &
    pft_photo_model, veg_index,                                                &
    som_n_sample, ci_lo, ci_hi,                                                &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
    km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,                 &
    l_multilayer, r_bl, r_ca, g_oth,                                           &
    kmax_ref, conductance_b, conductance_c,                                    &
    psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                    &
! OUT
    ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample,         &
    kl_hc_sample,                                                              &
    psi_stem_sample = psi_stem_sample                                          &
        )

CALL stom_opt_profit_max_select(                                               &
    land_pts, open_pts_flat, open_index_flat, veg_index,                       &
    som_n_sample,                                                              &
    al_sample, kl_sample, psi_sample, gl_sample, psi_root_zone, kcrit,         &
    gl_max, rd,                                                                &
    kl_hc_sample, kl_hc_max, l_xylem_impairment,                               &
! OUT
    l_good_sample_flat, carbon_gain, hydraulic_cost, profit,                   &
    optimal_index_flat                                                         &
        )

DO j = 1, open_pts_flat
  CALL take_sample( j, veg_index(open_index_flat(j)), optimal_index_flat(j) )
END DO

END SUBROUTINE profit_max_flat

!-----------------------------------------------------------------------------
! Allocate the arrays of the grid (flat) searches, once per call of
! stom_opt_mod. The bounded search does not need them.
!-----------------------------------------------------------------------------
SUBROUTINE alloc_flat()
IF ( ALLOCATED(ci_sample) ) RETURN
ALLOCATE( ci_sample(0:som_n_sample, open_pts),                                 &
          al_sample(0:som_n_sample, open_pts),                                 &
          gl_sample(0:som_n_sample, open_pts),                                 &
          kl_sample(0:som_n_sample, open_pts),                                 &
          kl_hc_sample(0:som_n_sample, open_pts),                              &
          kl_SOX(0:som_n_sample, open_pts),                                    &
          psi_sample(0:som_n_sample, open_pts),                                &
          SOX_mean_psi(0:som_n_sample, open_pts),                              &
          el_sample(0:som_n_sample, open_pts),                                 &
          profit(0:som_n_sample, open_pts),                                    &
          carbon_gain(0:som_n_sample, open_pts),                               &
          hydraulic_cost(0:som_n_sample, open_pts),                            &
          psi_stem_sample(0:som_n_sample, open_pts),                           &
          l_good_sample(som_n_sample, open_pts),                               &
          l_good_sample_flat(som_n_sample, open_pts),                          &
          optimal_index_flat(open_pts) )
END SUBROUTINE alloc_flat

!-----------------------------------------------------------------------------
! Copy Ci sample k of open point j of a grid search to the outputs at land
! point l. Sample 0 is the closed state. A closed result (k = 0, or a
! stomatal conductance of 0) has net photosynthesis -rd, the dark
! respiration, and a carbon gain and hydraulic cost of 0. An open result
! has net photosynthesis no lower than -rd.
!-----------------------------------------------------------------------------
SUBROUTINE take_sample( j, l, k )
INTEGER, INTENT(IN) :: j, l, k

ci(l) = MAX(0.0, ci_sample(k,j))

! Closed state (index 0) keeps al = -rd, as for points that were never
! open; an open optimum is clipped at -rd (as closed).
IF (k == 0) THEN
  al(l) = al_sample(0,j)
ELSE
  al(l) = MAX(-rd(l), al_sample(k,j))
END IF

gl(l) = MAX(0.0, gl_sample(k,j))
psi_leaf(l) = psi_sample(k,j)
! psi_stem, so the per-segment stem damage follows the stem's water
! potential, not the root zone's.
IF (PRESENT(psi_stem)) psi_stem(l) = psi_stem_sample(k,j)
el(l) = el_sample(k,j)
leaf_k(l) = kl_sample(k,j)
carbon_gain_out(l) = carbon_gain(k,j)
hydraulic_cost_out(l) = hydraulic_cost(k,j)

! gl = 0 means closed: al = -rd, as in set_closed, no gain or cost.
IF (k == 0 .OR. gl(l) <= 0.0) THEN
  IF (gl(l) <= 0.0) al(l) = -rd(l)
  carbon_gain_out(l) = 0.0
  hydraulic_cost_out(l) = 0.0
END IF

END SUBROUTINE take_sample

!-----------------------------------------------------------------------------
! Choose the stomatal opening by SOX_opt (Eller et al. 2018) on an even grid of
! Ci from ccp to ca, when the bounded search is off. Each feasible sample is
! scored
!
!   profit = A (k(psi_bar) - kcrit) / (k_top - kcrit)
!
! where A is net photosynthesis, psi_bar the mean of the leaf and root-zone
! water potentials, k(psi_bar) the conductance of the whole soil-to-leaf
! path with every part at psi_bar, and k_top the largest such k over the
! feasible samples. The sample with the largest profit is taken. SOX_opt
! cannot be run with xylem impairment (it is rejected at start-up).
!-----------------------------------------------------------------------------
SUBROUTINE sox_flat_search()

CALL alloc_flat()

carbon_gain(:,:) = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:) = 0.0
kl_SOX(:,:) = 0.0
ci_lo(:) = MAX(ccp(:), 0.0)
ci_hi(:) = ca(:)

CALL stom_opt_mod_ci(                                                          &
! IN
    land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,           &
    som_n_sample, ci_lo, ci_hi,                                                &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
    km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,                 &
    l_multilayer, r_bl, r_ca, g_oth,                                           &
    kmax_ref, conductance_b, conductance_c,                                    &
    psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                    &
! OUT
    ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample,         &
    kl_hc_sample,                                                              &
    psi_stem_sample = psi_stem_sample                                          &
        )

! ------------------------------------------------------------------
!  Make mask of those samples that are physically reasonable.
! ------------------------------------------------------------------
DO j = 1, open_pts
  l = veg_index(open_index(j))
  l_good_sample(:,j) = som_feasible( som_n_sample, psi_sample(1:,j),           &
                                     kl_sample(1:,j), gl_sample(1:,j),         &
                                     psi_root_zone(l), kcrit(l), gl_max(l),    &
                                     .TRUE. )
END DO

! SOX_opt (Eller et al. 2018): the conductance is taken at psi_bar,
! the mean of the leaf and root zone water potentials. Here it is the
! whole path (root, stem and leaf segments and the soil-to-root link
! in series, all at psi_bar; k_path_zero_flow), the same path as
! psi_leaf and the k > kcrit test.
DO j = 1, open_pts
  l = veg_index(open_index(j))

  SOX_mean_psi(:,j) = 0.5 * (psi_sample(:,j) + psi_root_zone(l))
  DO i = 1, som_n_sample
    kl_SOX(i,j) = k_path_zero_flow(pft, l, kmax(l), SOX_mean_psi(i,j))
  END DO
END DO

! Get the maximum xylem conductance for each land point.
k_top = MAXVAL(kl_SOX(1:,:), MASK = l_good_sample, DIM = 1)

! CG = An
carbon_gain(:,:) = al_sample(:,:)

! Loop over land points with open stomata
DO j = 1, open_pts
  l = veg_index(open_index(j))

  ! HC = 1 - (k - kcrit) / (max_k - kcrit)
  hydraulic_cost(:,j) = 1.0 - (kl_SOX(:,j)-kcrit(l))                           &
                              /(k_top(j)-kcrit(l))

  ! Profit = CG*(1-HC)
  profit(:,j) = carbon_gain(:,j)*(1.0-hydraulic_cost(:,j))
END DO ! j open points

! ------------------------------------------------------------------
!  Identify maximum profit and set leaf properties to those of the
!   optimal state. For what the zeroth sample is, see the closed state,
!   entry 0, in stom_opt_mod_ci.
! ------------------------------------------------------------------
DO j = 1, open_pts
  l = veg_index(open_index(j))

  optimal_index = MAXLOC(profit(1:,j), DIM = 1,                                &
                         MASK = l_good_sample(:,j))
  CALL take_sample( j, l, optimal_index )
END DO

END SUBROUTINE sox_flat_search

END SUBROUTINE stom_opt_mod

!-----------------------------------------------------------------------------
! Return the transpiration E (mol H2O m-2 s-1) through a conductance to water
! vapour g (m s-1), at specific-humidity deficit dq (kg kg-1), surface pressure
! P = pstar (Pa) and temperature T = t (K), so a stomatal conductance can be
! turned into the water the leaf loses. It is g times the molar density of air,
! P / (R T), times the deficit as a mole fraction, vpd / P:
!
!   E   = g (P / (R T)) (vpd / P) = g vpd / (R T)
!   vpd = dq P / repsilon
!
! with vpd the vapour pressure deficit (Pa), R the gas constant and
! repsilon the ratio of the molar masses of water and dry air. The Ci
! searches and sf_stom all use this one function.
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION e_of_g( g, dq, pstar, t ) RESULT( e )

USE planet_constants_mod, ONLY: repsilon
USE c_rmol, ONLY: rmol

IMPLICIT NONE

REAL(KIND=real_jlslsm), INTENT(IN) :: g
                            ! Conductance (m s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: dq
                            ! Specific-humidity deficit (kg kg-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: pstar
                            ! Surface pressure (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: t
                            ! Temperature (K).
REAL(KIND=real_jlslsm) :: e
                            ! Transpiration (mol H2O m-2 s-1).

e = dq * pstar / repsilon * g / (rmol * t)

END FUNCTION e_of_g

!-----------------------------------------------------------------------------
! Return the net photosynthesis A = W - rd (mol CO2 m-2 s-1) of a leaf at
! internal CO2 partial pressure ci (Pa): the gross rate W minus the dark
! respiration rd. Every Ci search in this module uses this one function to
! find the carbon gain of a candidate Ci. It follows leaf_limits_mod.F90 and
! leaf.F90.
!
! Collatz (C3 or C4): three limiting rates, set by Rubisco (wcarb), by
! light (wlite) and by the export of products (wexpt), are blended in two
! steps, each taking the smaller root of a quadratic:
!
!   beta1 W1^2 - (wcarb + wlite) W1 + wcarb wlite = 0
!   beta2 W^2  - (W1 + wexpt) W     + W1 wexpt    = 0
!
! (Cox, HCTN 24, equations 54 and 55). beta1 and beta2 set how sharply one
! limit gives way to the next.
!
! Farquhar (C3): W is the smaller of the Rubisco and light limited rates.
! With Johnson-Berry photosynthesis (photo_model) the light-limited rate is
! scaled by jb_eta_scale(ccp, ci).
!
! In a layer of the multilayer canopy (l_multilayer) the light-limited rate
! is rescaled to that layer's sunlit or shaded leaves: by
! fapar_lf ipar / apar (Collatz) or by je_ratio (Farquhar).
!
! MDK, 23/10/26
!-----------------------------------------------------------------------------
FUNCTION leaf_net_photo( pft, pft_photo_model, l_multilayer, ci, vcmax, kc,    &
                         ko, oi, ccp, km, je, acr, apar, fapar_lf, ipar,       &
                         je_ratio, pstar, rd ) RESULT( al )

USE jules_vegetation_mod, ONLY: photo_collatz, photo_farquhar,                 &
                                photo_johnson, photo_model
USE jb_photo_mod, ONLY: jb_eta_scale
USE pftparm, ONLY: c3, alpha
USE jules_surface_mod, ONLY: fwe_c3, fwe_c4, beta1, beta2
USE ereport_mod, ONLY: ereport

IMPLICIT NONE

INTEGER, INTENT(IN) :: pft, pft_photo_model
                            ! Plant functional type and its photosynthesis
                            ! model (photo_collatz or photo_farquhar).
LOGICAL, INTENT(IN) :: l_multilayer
                            ! .TRUE.: a layer of the multilayer canopy.
REAL(KIND=real_jlslsm), INTENT(IN) :: ci
                            ! Internal CO2 partial pressure (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: vcmax, kc, ko, oi, ccp, km, je, acr,     &
                                      apar, fapar_lf, ipar, je_ratio, pstar, rd
                            ! The point's photosynthesis inputs, as the
                            ! arguments of stom_opt_mod of the same names.
REAL(KIND=real_jlslsm) :: al
                            ! Net photosynthesis (mol CO2 m-2 s-1).

REAL(KIND=real_jlslsm) :: wcarb, wlite, wexpt, wl, b1, b2, b3
                            ! Rubisco, light and export limited rates, the
                            ! gross rate (mol CO2 m-2 s-1) and the
                            ! quadratic's coefficients.
INTEGER :: errcode

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEAF_NET_PHOTO'

SELECT CASE ( pft_photo_model )
CASE ( photo_collatz )
  IF (c3(pft) == 1) THEN
    wcarb = vcmax * (ci - ccp) / (ci + kc * (1.0 + oi / ko))
    wlite = alpha(pft) * acr * (ci - ccp) / (ci + 2.0 * ccp)
    wlite = MAX(wlite, TINY(1.0e0))
    IF (l_multilayer) wlite = wlite / apar * fapar_lf * ipar
    wexpt = fwe_c3 * vcmax
  ELSE
    wcarb = vcmax
    wlite = MAX(alpha(pft) * acr, TINY(1.0e0))
    IF (l_multilayer) wlite = wlite / apar * fapar_lf * ipar
    wexpt = fwe_c4 * vcmax * ci / pstar
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
  wcarb = vcmax * (ci - ccp) / (ci + km)
  wlite = je / 4.0 * (ci - ccp) / (ci + 2.0 * ccp)
  IF (photo_model == photo_johnson) wlite = wlite * jb_eta_scale(ccp, ci)
  wlite = MAX(wlite, TINY(1.0e0))
  IF (l_multilayer) wlite = wlite * je_ratio
  wl = MIN(wcarb, wlite)

CASE DEFAULT
  wl = 0.0
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'pft_photo_model should be photo_collatz or photo_farquhar')
END SELECT

al = wl - rd

END FUNCTION leaf_net_photo

!-----------------------------------------------------------------------------
! Work out the leaf states on an even grid of Ci, which the grid searches and
! the fast path then choose from. At each open point n_sample values of Ci are
! spaced evenly from ci_lo up to, but not including, ci_hi, so the top sample
! stays below ca when ci_hi = ca (where gl would be infinite). For each sample
! it returns net photosynthesis A, the stomatal conductance gl that supplies
! the CO2, the transpiration E, and from the plant's water path (som_path_psi)
! the leaf water potential, the path conductance k and the conductance the
! hydraulic cost uses (the intact path's with xylem impairment). With
! l_som_coupled_e, or where there is a leaf boundary layer resistance
! (r_bl > 0), the CO2 and the water vapour also cross the leaf boundary
! layer (gl_co2_bl, g_eff_e).
!
! Entry 0 of each output is the closed-stomata state: Ci = ca, A = -rd, no
! conductance or transpiration, the leaf at the root-zone water potential,
! and k that of the path with no flow.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_mod_ci(                                                    &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        n_sample, ci_lo, ci_hi,                                                &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp,              &
        pstar, km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax,             &
        kcrit, l_multilayer, r_bl, r_ca, g_oth,                                &
        kmax_ref, conductance_b, conductance_c,                                &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                &
! OUT
        ci_sample, al_sample, gl_sample, kl_sample, psi_sample,el_sample,      &
        kl_hc_sample,                                                          &
! OUT (optional)
        psi_stem_sample                                                        &
)

USE jules_vegetation_mod, ONLY: l_som_plant_segments, l_som_coupled_e

USE jules_surface_mod, ONLY: ratio

USE c_rmol, ONLY: rmol

USE plant_hydraulic_path_mod, ONLY: grav_node_drop

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
 ci_sample(0:n_sample, open_pts)                                               &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:n_sample, open_pts)                                               &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:n_sample, open_pts)                                               &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:n_sample, open_pts)                                               &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
,psi_sample(0:n_sample, open_pts)                                              &
                            ! Water potential at leaf (Pa)
,el_sample(0:n_sample, open_pts)                                               &
                            ! Transpiration rate (mol H2O/m2/s)
,kl_hc_sample(0:n_sample, open_pts)
                            ! Conductance used for the hydraulic cost
                            ! (mol m-2 s-1 Pa-1): kl_sample, or with
                            ! l_xylem_impairment the conductance at
                            ! psi_sample on the unimpaired PFT curve.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem_sample(0:n_sample, open_pts)
                            ! Stem water potential (Pa) of each sample: the
                            ! stem segment's leaf-side end with segments and
                            ! impairment, else psi_root_zone.

!-----------------------------------------------------------------------------
! Local integer variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 i,j,l
                            ! Internal iterators

!-----------------------------------------------------------------------------
! Local real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 el_closed(1, open_pts)                                                        &
,psi_closed(1, open_pts)                                                       &
,kl_closed(1, open_pts)                                                        &
                            ! Zero-flow solve for the closed-stomata state
                            ! (entry 0): transpiration, leaf water
                            ! potential and xylem conductance.
,conv                                                                          &
                            ! Factor for converting mol/m3 into Pa (J/m3), of
                            ! a point.
,psi_stem_work(0:n_sample, open_pts)
                            ! Stem water potential of each sample (Pa).

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

! Entry 0 is the closed state (see the reset at the end of this routine).
! Its ci and psi are set here so the photosynthesis below has a valid ci.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  ci_sample(0,j)  = ca(l)
  psi_sample(0,j) = psi_root_zone(l) - grav_node_drop(l)

  DO i = 1, n_sample
    !NOTE: Array entry 1 is at ci_lo (hence i-1), and the last step is less
    !       than ci_hi to avoid infinite transpiration rates when ci_hi is
    !       the full-range bound ca(l) (see the callers in stom_opt_mod).
    ci_sample(i,j) = ci_lo(l) + (i-1) * (ci_hi(l) - ci_lo(l)) / n_sample
  END DO
END DO

!-----------------------------------------------------------------------------
! Net photosynthesis of every sample (leaf_net_photo).
!-----------------------------------------------------------------------------
DO j = 1,open_pts
  l = veg_index(open_index(j))
  DO i = 0, n_sample
    al_sample(i,j) = leaf_net_photo( pft, pft_photo_model, l_multilayer,       &
                                     ci_sample(i,j), vcmax(l), kc(l), ko(l),   &
                                     oi(l), ccp(l), km(l), je(l), acr(l),      &
                                     apar(l), fapar_lf(l), ipar(l),            &
                                     je_ratio(l), pstar(l), rd(l) )
  END DO
END DO

! ----------------------------------------------------------------------------
! Calculate stomatal conductance for H2O (m/s),
!
!   gl = ratio A R T_leaf / (ca - ci),
!
! with ratio = 1.6 the ratio of the diffusivities of H2O and CO2 in air.
! ----------------------------------------------------------------------------
DO j = 1,open_pts
  l = veg_index(open_index(j))

  !---------------------------------------------------------------------------
  ! Calculate the factor for converting mol/m3 into Pa (J/m3).
  !---------------------------------------------------------------------------
  conv = rmol * t_leaf(l)

  !---------------------------------------------------------------------------
  ! Diagnose the leaf conductance. ca - ci is floored at 1e-2 Pa, so a
  ! sample next to ca gives a large, finite gl, which the feasibility test
  ! then rejects.
  !---------------------------------------------------------------------------
  gl_sample(:,j) = ratio * ( (al_sample(:,j) * conv)                           &
                             / MAX(ca(l) - ci_sample(:,j), 1.0e-2_real_jlslsm) )

  ! Coupled, or with a leaf boundary layer: the CO2 crosses the boundary
  ! layer too (gl_co2_bl).
  IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
    gl_sample(:,j) = gl_co2_bl(al_sample(:,j),                                 &
                               MAX(ca(l) - ci_sample(:,j),                     &
                                   1.0e-2_real_jlslsm), conv, r_bl(l))
  END IF

END DO

! ----------------------------------------------------------------------------
! Calculate transpiration rate (e_of_g): E = gl vpd / (R T).
! ----------------------------------------------------------------------------
DO j = 1,open_pts
  l = veg_index(open_index(j))
  el_sample(:,j) = e_of_g(gl_sample(:,j), dq(l), pstar(l), t_leaf(l))

  ! Coupled, or through a leaf boundary layer: E = D g_eff(g), dq the
  ! deficit to the level-1 air.
  IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
    el_sample(:,j) = e_of_g(g_eff_e(MAX(gl_sample(:,j), 0.0_real_jlslsm),     &
                                    r_bl(l), r_ca(l), g_oth(l)),               &
                            dq(l), pstar(l), t_leaf(l))
  END IF
END DO

! Transpiration can't be negative
el_sample(:,:) = MAX(0.0, el_sample(:,:))

! Approximate the leaf water potential and conductance. The solvers take
! (n_sample, open_pts) arrays, so the sections 1:n_sample are passed.
! One call on the path the optimiser uses (som_path_psi): segments, the
! impaired curve or the PFT curve. kl_hc_sample is the conductance the
! hydraulic cost uses (the intact path under impairment).
CALL som_path_psi( pft, n_sample, land_pts, open_pts, veg_index, open_index,   &
                   el_sample(1:n_sample,:), psi_root_zone, kmax_ref, kmax,     &
                   kcrit, conductance_b, conductance_c, psi_leaf_extreme,      &
                   psi_root_extreme, l_xylem_impairment,                       &
                   psi_sample(1:n_sample,:), kl_sample(1:n_sample,:),          &
                   psi_stem = psi_stem_work(1:n_sample,:),                     &
                   leaf_k_intact = kl_hc_sample(1:n_sample,:) )

! Reset entry 0 to the closed state, which the calculations above filled
! with the open state at ci = ca: no conductance or transpiration, net
! assimilation -rd, and the leaf at the root zone water potential. Its
! conductance is that of the whole path at E = 0, K(psi_root_zone), on the
! same path as the samples (with the soil-to-root link and any
! impairment), so a shut canopy does not report 100 % loss of conductance
! (PLC_pft). Entry 0 is never a profit candidate.
el_closed(:,:) = 0.0

CALL som_path_psi( pft, 1, land_pts, open_pts, veg_index, open_index,          &
                   el_closed, psi_root_zone, kmax_ref, kmax, kcrit,            &
                   conductance_b, conductance_c, psi_leaf_extreme,             &
                   psi_root_extreme, l_xylem_impairment,                       &
                 ! INTENT OUT
                   psi_closed, kl_closed )

! Stem water potential samples (only when the caller asks for them): with
! xylem impairment on the segmented path, each Ci sample keeps the stem
! (root-stem node) potential from its own path solve; otherwise it is the
! stem at zero flow, the root-zone potential less the gravity drop at the
! root-stem node (grav_node_drop, 0 without it). Entry 0 (closed, no flow)
! is the stem at zero flow either way.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  al_sample(0,j)  = -rd(l)
  gl_sample(0,j)  = 0.0
  el_sample(0,j)  = 0.0
  kl_sample(0,j)  = kl_closed(1,j)
  kl_hc_sample(0,j) = kl_sample(0,j)
  psi_sample(0,j) = psi_root_zone(l) - grav_node_drop(l)
END DO

IF (PRESENT(psi_stem_sample)) THEN
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    IF (l_xylem_impairment .AND. l_som_plant_segments) THEN
      psi_stem_sample(1:,j) = psi_stem_work(1:,j)
    ELSE
      psi_stem_sample(1:,j) = psi_root_zone(l) - grav_node_drop(l)
    END IF
    psi_stem_sample(0,j) = psi_root_zone(l) - grav_node_drop(l)
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod_ci

!-----------------------------------------------------------------------------
! Work out the leaf water potential psi_leaf (Pa) and the conductance of the
! soil-to-leaf path, k = -dE/dpsi_leaf (mol m-2 s-1 Pa-1), which give the
! hydraulic cost of a candidate Ci, for given transpiration rates e_leaf
! (mol H2O m-2 s-1), on the same path the optimiser uses:
!
!   - xylem impairment with plant segments: root, stem and leaf in series,
!     with the stem and leaf capped at kmax / kmax_ref of their intact
!     conductance (the damage remembered from earlier droughts). This also
!     gives the stem water potential, psi_stem.
!   - xylem impairment alone: the impaired whole-plant vulnerability curve.
!   - otherwise: the whole-plant vulnerability curve (conductance_b and
!     conductance_c).
!
! leaf_k_intact, when asked for, is the conductance the hydraulic cost
! uses: that of the intact (undamaged) path at the same water potentials,
! with the soil-to-root link where the path has one. It equals leaf_k
! without impairment.
!
! The Ci samples, the closed state and the one-Ci evaluations of the
! bounded search all call this routine, so they all solve the same path.
! (The cuticular loss in sf_stom calls the segment solver directly.)
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
SUBROUTINE som_path_psi( pft, n_e_leaf, land_pts, open_pts, veg_index,         &
                         open_index, e_leaf, psi_root_zone, kmax_ref, kmax,    &
                         kcrit, conductance_b, conductance_c,                  &
                         psi_leaf_extreme, psi_root_extreme,                   &
                         l_xylem_impairment, psi_leaf, leaf_k, psi_stem,       &
                         leaf_k_intact )

USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls, xylem_conductance_jls
USE xylem_impairment_mod, ONLY: psi_leaf_impaired_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: xylem_f
USE plant_hydraulic_path_mod, ONLY: psi_leaf_segments_jls,                     &
                                    psi_leaf_segments_lut_jls, ksr_path,       &
                                    som_psi_in_min
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut,                &
                                l_som_plant_segments
USE pftparm, ONLY: conductance_b_pft, conductance_c_pft
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

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


REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
  leaf_k_intact(n_e_leaf, open_pts)
      ! The conductance the hydraulic cost uses: that of the intact
      ! (undamaged) path at the same water potentials. Equal to leaf_k
      ! without xylem impairment.

REAL(KIND=real_jlslsm) :: kcap_frac(land_pts)
REAL(KIND=real_jlslsm) ::                                                      &
  kmax_open(open_pts), kcrit_open(open_pts), b_open(open_pts),                 &
  c_open(open_pts)
      ! Intact curve for each open point (impairment without segments).
REAL(KIND=real_jlslsm) :: k_s, psi_in_s
      ! Soil-to-root conductance of the path, and the root water potential
      ! of one sample.
INTEGER :: i, j, l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='SOM_PATH_PSI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Absent optional arguments (psi_stem, leaf_k_intact) are passed on as
! absent, so each solver is called once.
IF ( l_xylem_impairment .AND. l_som_plant_segments ) THEN
  ! Segments: the memory model's cap on the stem and leaf segments (the
  ! root is not damaged); leaf_k_intact is the intact segments' whole-plant
  ! conductance at the same water potentials. (Only the open points'
  ! kcap_frac are read, so the rest of the grid is not set: one-sample
  ! callers pay for one point, not the whole grid.)
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    kcap_frac(l) = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm))
  END DO

  IF ( som_psi_solver == psi_solver_lut ) THEN
    ! som_psi_solver = 3: the same capped segments from the supply tables.
    CALL psi_leaf_segments_lut_jls( pft, n_e_leaf, land_pts, open_pts,         &
                                    veg_index, open_index, e_leaf,             &
                                    psi_root_zone, kmax_ref, psi_leaf,         &
                                    leaf_k, kcap_frac = kcap_frac,             &
                                    leaf_k_intact = leaf_k_intact,             &
                                    stem_psi = psi_stem )
  ELSE
    CALL psi_leaf_segments_jls( pft, n_e_leaf, land_pts, open_pts,             &
                                veg_index, open_index, e_leaf,                 &
                                psi_root_zone, kmax_ref, kcrit, psi_leaf,      &
                                leaf_k, kcap_frac = kcap_frac,                 &
                                leaf_k_intact = leaf_k_intact,                 &
                                stem_psi = psi_stem )
  END IF

ELSE IF ( l_xylem_impairment ) THEN
  ! Leaf water potential and conductance on the impaired vulnerability
  ! curve.
  CALL psi_leaf_impaired_jls( pft, n_e_leaf, land_pts, open_pts, veg_index,    &
                              open_index, e_leaf, psi_root_zone, kmax_ref,     &
                              kmax, kcrit, conductance_b, conductance_c,       &
                              psi_leaf_extreme, psi_root_extreme,              &
                              psi_leaf, leaf_k )

  IF ( PRESENT(leaf_k_intact) ) THEN
    ! The hydraulic cost uses the unimpaired conductance at the (impaired)
    ! leaf water potential.
    DO j = 1, open_pts
      l = veg_index(open_index(j))
      kmax_open(j)  = kmax_ref(l)
      kcrit_open(j) = kcrit(l)
    END DO
    b_open(:) = conductance_b_pft(pft)
    c_open(:) = conductance_c_pft(pft)

    CALL xylem_conductance_jls( pft, n_e_leaf, open_pts, psi_leaf,             &
                                kmax_open, kcrit_open, b_open, c_open,         &
                                leaf_k_intact )

    ! With the soil-to-root conductance in series, the cost uses the whole
    ! path as the unimpaired cost does: the intact xylem at the same water
    ! potentials, in series with the soil (which is not damaged, so the
    ! stomata see it). As psi_leaf_jls, the marginal conductance of the
    ! path is
    !
    !   k(psi_leaf) / (1 + k(psi_in) / K_s),   psi_in = psi_src - E / K_s
    !
    ! with k the intact curve.
    DO j = 1, open_pts
      l = veg_index(open_index(j))
      k_s = ksr_path(l)
      IF ( k_s < 0.0 ) CYCLE
      DO i = 1, n_e_leaf
        IF ( k_s > TINY(1.0_real_jlslsm) ) THEN
          psi_in_s = MAX(psi_root_zone(l) - e_leaf(i,j) / k_s,                 &
                         som_psi_in_min)
          leaf_k_intact(i,j) = leaf_k_intact(i,j)                              &
                               / (1.0 + kmax_ref(l) * xylem_f(pft, psi_in_s)   &
                                        / k_s)
        ELSE IF ( e_leaf(i,j) > 0.0 ) THEN
          leaf_k_intact(i,j) = 0.0
        END IF
      END DO
    END DO
  END IF

ELSE
  CALL psi_leaf_jls( pft, n_e_leaf, land_pts, open_pts, veg_index,             &
                     open_index, e_leaf, psi_root_zone, kmax, kcrit,           &
                     conductance_b, conductance_c,                             &
                   ! INTENT OUT
                     psi_leaf, leaf_k )

  IF ( PRESENT(leaf_k_intact) ) leaf_k_intact(:,:) = leaf_k(:,:)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE som_path_psi

!-----------------------------------------------------------------------------
! Choose the stomatal opening by profit max from a grid of Ci samples already
! evaluated by stom_opt_mod_ci: mark the feasible samples, score each one, and
! return the best sample of each open point.
!
! A sample is feasible if the leaf is no wetter than the root zone, the
! path conductance is above kcrit, and gl is within gl_max (when set). The
! score is
!
!   profit = (A + g_off) / A_max - (k_max - k) / (k_max - kcrit)
!
! with A the net photosynthesis, g_off = rd for the gross gain
! (l_som_gain_gross) and 0 otherwise, A_max the largest A + g_off over the
! samples the hydraulics allow (leaving out the gl cap), k the conductance
! the cost uses, and k_max the largest such k over the feasible samples
! (with xylem impairment, kl_hc_max: the intact path with no flow). Both
! terms run from 0 to 1, so no exchange rate between carbon and water is
! needed.
!
! optimal_index(j) is 0 when no sample is feasible or A_max <= 0. The
! caller then takes sample 0, the closed state (see stom_opt_mod_ci).
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_profit_max_select(                                         &
! IN
        land_pts, open_pts, open_index, veg_index, n_sample,                   &
        al_sample, kl_sample, psi_sample, gl_sample, psi_root_zone, kcrit,     &
        gl_max, rd,                                                            &
        kl_hc_sample, kl_hc_max, l_xylem_impairment,                           &
! OUT
        l_good_sample, carbon_gain, hydraulic_cost, profit, optimal_index      &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: l_som_gain_gross

INTEGER, INTENT(IN) ::                                                         &
  land_pts                                                                     &
                            ! Number of land points.
, open_pts                                                                     &
                            ! Number of open points to select for.
, open_index(land_pts)                                                         &
                            ! Index of the open points on veg_index.
, veg_index(land_pts)                                                          &
                            ! Index of vegetated points on the land grid.
, n_sample
                            ! Number of Ci samples (entry 0 is the closed
                            ! state).

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  al_sample(0:n_sample, open_pts)                                              &
                            ! Net photosynthesis of each sample
                            ! (mol CO2 m-2 s-1).
, kl_sample(0:n_sample, open_pts)                                              &
                            ! Whole-path conductance -dE/dpsi_leaf of each
                            ! sample (mol m-2 s-1 Pa-1).
, psi_sample(0:n_sample, open_pts)                                             &
                            ! Leaf water potential of each sample (Pa).
, gl_sample(0:n_sample, open_pts)                                              &
                            ! Stomatal conductance for H2O of each sample
                            ! (m/s).
, psi_root_zone(land_pts)                                                      &
                            ! Root zone water potential (Pa).
, kcrit(land_pts)                                                              &
                            ! Critical conductance (mol m-2 s-1 Pa-1).
, gl_max(land_pts)                                                             &
                            ! Cap on gl (m/s); <= 0 disables it.
, kl_hc_sample(0:n_sample, open_pts)                                           &
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

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  carbon_gain(0:n_sample, open_pts)                                            &
                            ! A / max A of each sample (dimensionless).
, hydraulic_cost(0:n_sample, open_pts)                                         &
                            ! Loss of conductance of each sample, as a
                            ! fraction of max_k - kcrit (dimensionless).
, profit(0:n_sample, open_pts)
                            ! carbon_gain - hydraulic_cost.

INTEGER, INTENT(OUT) :: optimal_index(open_pts)
                            ! Sample with the largest profit; 0 (closed) if
                            ! none is feasible.

INTEGER :: j, l

REAL(KIND=real_jlslsm) :: max_al(open_pts), k0_ref(open_pts)
    ! Gain scale (largest A + gain_off over the hydraulically feasible
    ! samples) and the reference conductance of the hydraulic cost.
REAL(KIND=real_jlslsm) :: gain_off(open_pts)
    ! 0 (net A) or rd (gross A, l_som_gain_gross).
LOGICAL :: l_hyd_sample(n_sample)
    ! Hydraulically feasible samples (no gl cap), for max_al.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_PROFIT_MAX_SELECT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------------
!  Make mask of those samples that are physically reasonable.
! ----------------------------------------------------------------------------
DO j = 1, open_pts
  l = veg_index(open_index(j))
  l_good_sample(:,j) = som_feasible( n_sample, psi_sample(1:,j),               &
                                     kl_sample(1:,j), gl_sample(1:,j),         &
                                     psi_root_zone(l), kcrit(l), gl_max(l),    &
                                     .TRUE. )
END DO

! ----------------------------------------------------------------------------
!  Calculate the profit for each leaf state.
!
!    CG = An/max(An)   (An + Rd with l_som_gain_gross)
!    HC = (k0_ref - k)/(k0_ref - k_crit)

!    Profit = CG-HC
!
! Why this form (Sperry et al. 2017): HC is the fraction of the usable
! plant conductance that opening to this Ci would lose, 0 with no flow and
! 1 at hydraulic failure (k = k_crit). The gain CG is also 0-1, so
! profit = gain - cost needs no exchange rate between carbon and water:
! the optimum is where the marginal gain equals the marginal loss of
! conductance.
! ----------------------------------------------------------------------------
! gain offset: 0 (net An) or Rd (gross A), see l_som_gain_gross
DO j = 1, open_pts
  l = veg_index(open_index(j))
  gain_off(j) = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)

  ! The gain normalisation over the hydraulically feasible samples, without
  ! the gl cap (as Sperry et al. 2017; see scale_and_ends).
  l_hyd_sample(:) = som_feasible( n_sample, psi_sample(1:,j),                  &
                                  kl_sample(1:,j), gl_sample(1:,j),            &
                                  psi_root_zone(l), kcrit(l), 0.0_real_jlslsm, &
                                  .TRUE. )
  max_al(j) = MAXVAL(al_sample(1:,j) + gain_off(j), MASK = l_hyd_sample(:))
END DO

! (The store coupled to the optimiser needs the bounded search, so its flow
! is never held fixed in this flat one.)
IF (l_xylem_impairment) THEN
  ! Set the maximum conductance as equal to the current (unimpaired) root
  ! zone conductance.
  DO j = 1, open_pts
    k0_ref(j) = kl_hc_max(veg_index(open_index(j)))
  END DO
ELSE
  k0_ref = MAXVAL(kl_hc_sample(1:,:), MASK = l_good_sample, DIM = 1)
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

    ! HC = (k0_ref - k)/(k0_ref - k_crit). Without impairment k0_ref(j) is
    ! the largest k of the feasible samples, all of which have k > kcrit,
    ! so the denominator is positive. With impairment k0_ref(j) is
    ! kl_hc_max, not taken from a sample, and the MAX(..., TINY) is what
    ! stops a division by zero.
    hydraulic_cost(:,j) = (k0_ref(j) - kl_hc_sample(:,j))                      &
                          / MAX(k0_ref(j) - kcrit(l), TINY(1.0_real_jlslsm))

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
! Choose the stomatal opening of each open point, one point at a time, by
! searching for the Ci with the largest score inside the feasible Ci range (the
! bounded Ci search, som_ci_search = 2). Profit max, SOX_opt, CMax, CGain, CAP
! and ProfitMax2 use it. The score each model maximises is chosen by l_cmax,
! l_cgain, l_sox, l_cap and l_pm2 (see the function profit).
!
! As Ci rises, A rises and ca - Ci falls, so gl and E rise and the path
! conductance k falls. The feasible Ci therefore form one range, from ccp
! up to an upper edge ci_b. For each point:
!
!   1. Evaluate both ends of the full range (range_ends). ccp has the
!      largest k (no transpiration); if it is infeasible, every Ci is, and
!      the stomata close. They also close if A <= 0 at the top of the
!      range, since A only rises with Ci.
!   2. Find the edge ci_b (find_edge) by a bracketed search that keeps one
!      feasible and one infeasible end and places each new guess where the
!      straight line between them crosses the limit (the Illinois form of
!      regula falsi), on how far full evaluations are from the limits
!      (edge_by_margin).
!   3. Set the gain scale A_max (scale_and_ends): A at the edge set by the
!      hydraulics, found beyond the gl cap when the cap binds first
!      (hydraulic_edge), as in Sperry et al. (2017).
!   4. Find the Ci of largest profit between ccp and ci_b (golden_max) by a
!      golden-section search, which shrinks the bracket by the same
!      fraction, 0.618, at each step, until it is narrower than
!      (ci_b - ccp) / opt_resol or after n_iter steps. Every Ci in this
!      range is feasible, so the score has no sudden drop. The search
!      assumes one peak, which holds when A is concave in Ci and the cost
!      convex in E (Weibull c >= 1); where that may fail (l_multi_peak) a
!      coarse scan picks the peak first.
!   5. With xylem impairment and ximpair_hold_psi, steps 1 to 4 run on the
!      intact path, and Ci is then lowered until the impaired path gives
!      the same leaf water potential (hold_leaf_psi).
!
! The best Ci seen at any step, both ends included, is returned, with the
! marginal water price dA/dE there (lambda_diag). Each evaluation is for one
! point (eval_ci), so each point stops as soon as it has converged. With
! l_leaf_temp_ci the dark respiration in the chosen A, at the candidate's
! own leaf temperature, is returned too (rd_g), so the gross
! photosynthesis A + Rd is that of the same leaf temperature.
!
! MDK, 04/11/26
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
        psi_stem_g, lambda_g, rd_g,                                            &
! IN (optional)
        lai_cost                                                               &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: som_profit_model, cmax_profit_model,           &
                                cgain_profit_model, sox_opt_profit_model,      &
                                cap_profit_model, pm2_profit_model,            &
                                l_som_gain_gross,                              &
                                l_som_plant_segments, l_som_nsl,               &
                                ximpair_stomatal_response, ximpair_hold_gs,    &
                                ximpair_hold_psi, l_som_coupled_e,             &
                                photo_johnson
USE plant_hydraulic_path_mod, ONLY: l_kscale
USE plant_hydraulic_path_mod, ONLY: grav_node_drop
USE plant_hydraulic_path_mod, ONLY: k_path_zero_flow, kplant_scale,            &
                                    l_store_q_path, l_shared_stem_path,        &
                                    l_shared_whole_path, l_shared_gain_path,   &
                                    shared_kcrit_plant, shared_amax_other,     &
                                    shared_amax_self
USE pftparm, ONLY: psi_nsl_onset, psi_nsl0,                                    &
                   cmax_a, cmax_b, kmax_pft, cgain_varpi,                      &
                   conductance_b_pft, conductance_c_pft,                       &
                   conductance_c_seg
USE jules_surface_mod, ONLY: ratio
USE c_rmol, ONLY: rmol
USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls

INTEGER, INTENT(IN) ::                                                         &
  land_pts, pft, open_pts, open_index(land_pts), pft_photo_model,              &
  veg_index(land_pts), n_top, n_iter
                            ! n_top: grid size that sets ci_top without a gl
                            ! cap (som_n_sample); n_iter: largest number of
                            ! golden-section steps.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),          &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts),   &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), je(land_pts),    &
  t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts), ipar(land_pts),    &
  kmax(land_pts), kcrit(land_pts), gl_max(land_pts)
                            ! As the arguments of stom_opt_mod of the same
                            ! names (units and bases there).

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_bl(land_pts), r_ca(land_pts), g_oth(land_pts)
                            ! Coupled transpiration (see stom_opt_mod).

! Xylem impairment (as stom_opt_mod_ci): the leaf state on the impaired path
! (kmax, conductance_b/c and the extremes; segments: stem/leaf capped at
! kmax / kmax_ref), the hydraulic cost on the intact path (kmax_ref, the PFT
! curve) in series with the soil, relative to kl_hc_max.
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  kmax_ref(land_pts), conductance_b(land_pts), conductance_c(land_pts),        &
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
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: psi_stem_g(open_pts)
                            ! Stem segment leaf-side water potential at the
                            ! chosen state (xylem impairment with segments,
                            ! for the per-segment stem damage); psi_root_zone
                            ! otherwise and for closed points.
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: lai_cost(land_pts)
                            ! Leaf area of the call (see stom_opt_mod).
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: lambda_g(open_pts)
                            ! Marginal water price dA/dE at the chosen state
                            ! (mol CO2 per mol H2O; see lambda_eff in
                            ! stom_opt_mod), on the search path; 0 if closed.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: rd_g(open_pts)
                            ! Dark respiration inside al_g (mol CO2 m-2 s-1):
                            ! with l_leaf_temp_ci, rd at the chosen
                            ! candidate's leaf temperature; rd where the
                            ! point is closed or al_g is floored at -rd.
REAL(KIND=real_jlslsm) :: el_s1(1,1), psi_s1(1,1), kl_s1(1,1), st_s1(1,1)
                            ! stem_psi_of_choice: transpiration, leaf psi,
                            ! path conductance and stem psi of the chosen
                            ! state, as one-sample arrays.

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
INTEGER, PARAMETER :: n_hold_iter = 40
                            ! Hold leaf psi: maximum bisection steps.
REAL(KIND=real_jlslsm), PARAMETER :: hold_psi_tol = 100.0
                            ! Hold leaf psi: tolerance on psi_leaf (Pa).

INTEGER :: i, j, l, side, idx1(land_pts)
                            ! Loop counters, the open point and its land
                            ! point, the end of the edge bracket replaced
                            ! last (edge_by_margin), and the one-point
                            ! open index for the path solvers.
LOGICAL :: l_intact_search
                            ! Hold leaf psi (ximpair_hold_psi): .TRUE.
                            ! while the search runs on the intact path.
LOGICAL :: l_no_glcap
                            ! .TRUE. while the hydraulic edge is found for
                            ! the gain normalisation: the gl cap (gl_max)
                            ! is left out of the feasibility test.
REAL(KIND=real_jlslsm) :: hy_ci, hy_al
                            ! Ci and A at the hydraulic edge (the gain
                            ! normalisation, as Sperry et al. 2017).
LOGICAL :: l_shut
                            ! The point was closed by a step of the point
                            ! loop.
REAL(KIND=real_jlslsm) :: psi_star, hc_star, cg_star, s_lo, s_hi, s_mid
                            ! Hold leaf psi: target psi_leaf and the intact-path
                            ! gain/cost, and the Ci bisection bracket.
REAL(KIND=real_jlslsm) :: lam_ci1, lam_ci2
                            ! lambda_g: the lower and upper Ci of the
                            ! difference.
REAL(KIND=real_jlslsm) :: h_ci, h_al, h_gl, h_el, h_psi, h_kl
LOGICAL :: l_h_ok
                            ! Hold leaf psi: the feasible evaluation on the
                            ! impaired path with the largest Ci, and whether
                            ! there is one.
REAL(KIND=real_jlslsm) :: b_curve(land_pts), c_curve(land_pts)
                            ! The PFT's intact vulnerability curve on
                            ! land_pts, for the hold-psi search on the
                            ! intact path. (The other evaluations use the
                            ! point's own curve, conductance_b/c, as the
                            ! flat search does.)
LOGICAL ::                                                                     &
  ok_u,                                                                        &
                            ! The latest evaluation (eval_ci) is feasible.
  l_cmax, l_cgain, l_sox, l_cap, l_pm2,                                        &
                            ! The profit model is CMax, CGain, SOX_opt, CAP
                            ! or ProfitMax2.
  l_abs
                            ! Any of those five: profit in absolute carbon
                            ! units, on net A.

REAL(KIND=real_jlslsm) :: lam_a1, lam_e1
                            ! lambda_g: A and E at the lower Ci.
REAL(KIND=real_jlslsm) :: e_crit
                            ! ProfitMax2: transpiration at which the
                            ! conductance of the path the search runs on
                            ! falls to kcrit (mol H2O m-2 s-1, the basis of
                            ! el_u).
REAL(KIND=real_jlslsm) :: lai_call
                            ! Leaf area of this call (m2 leaf m-2 ground):
                            ! kmax_ref / kmax_pft, divided by kplant_scale
                            ! with kplant_model = 1, or lai_cost when given.
                            ! The CMax and CGain costs, given per unit leaf
                            ! area, are scaled by it to the basis of A
                            ! (whole canopy, the sunlit or shaded leaves, or
                            ! per unit leaf area of a layer's leaves).
REAL(KIND=real_jlslsm) :: k0_cgain
                            ! CGain's k_max, the plant's conductance at
                            ! psi = 0 (plant only, no soil link).
REAL(KIND=real_jlslsm) :: k0_sox
                            ! SOX's whole-path k at zero flow (psi_bar =
                            ! psi_root_zone), which normalises the SOX
                            ! factor (a constant for the point, so it does
                            ! not move the optimum).
REAL(KIND=real_jlslsm) :: kcrit_hc
                            ! kcrit of the hydraulic cost: kcrit, or with
                            ! the shared root and stem the whole plant's
                            ! (shared_kcrit_plant).

REAL(KIND=real_jlslsm) ::                                                      &
  ! Gain offset (0, or rd with l_som_gain_gross), gain scale A_max, reference
  ! conductance of the hydraulic cost (zero flow), the Ci range and the
  ! golden-section stop width
  g_off, max_al, k0_ref, ci_lo, ci_top, tol1,                                  &
  ! Evaluation state (eval_ci); al0_u is A without the nonstomatal limitation;
  ! khc_u the conductance for the hydraulic cost (kl_u unless impaired)
  last_ci, al_u, gl_u, el_u, psi_u, kl_u, al0_u, khc_u,                        &
  ! Rd inside al_u with l_leaf_temp_ci (state_at_leaf_temp)
  rd_u,                                                                        &
  ! Edge (last feasible) and the infeasible end of its bracket, with their
  ! feasibility margins (mrg_e, mrg_b; mrg_u that of the latest evaluation)
  e_ci, e_al, e_gl, e_el, e_psi, e_kl, mrg_e, b_ci, mrg_b, mrg_u,              &
  e_khc,                                                                       &
  ! edge_by_margin: the Ci tried
  ci_try,                                                                      &
  ! Best seen
  best_f, best_ci, best_al, best_gl, best_el, best_psi, best_kl, best_khc,     &
  ! Golden-section bracket, its interior points and their profits
  a, b, c_ci, d_ci, fc, fd
INTEGER :: i_best
REAL(KIND=real_jlslsm) :: f_scan, f_best_scan, d_scan
INTEGER, PARAMETER :: n_scan = 5
                            ! Coarse scan before the golden-section search
                            ! where profit may have several peaks (see
                            ! golden_max).
LOGICAL :: l_ltc_pt
                            ! l_leaf_temp_ci at this point (see eval_ci).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_BOUNDED_SEARCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Leaf psi from one-point calls of the path solvers (open index list idx1).
! With the soil-to-root link in series the root water potential falls with
! E, so the supply table is used only inside those solvers.
idx1(:) = 1

! CMax (stomata_model = 6): profit A_n - Theta(psi_leaf) in absolute carbon
! units instead of CG - HC (see the profit function below).
l_cmax = som_profit_model == cmax_profit_model

! CGain (stomata_model = 7): A_n - varpi (k_max - k)/k_max.
l_cgain = som_profit_model == cgain_profit_model

! SOX_opt (stomata_model = 5 with som_ci_search = 2): Eller et al. (2018)
! solved as Sabot et al. (2022) do,
!
!   A_n (k(psi_bar) - kcrit)/(k_ref - kcrit),
!
! k of the whole path at psi_bar, the mean of the root zone and leaf water
! potentials. The halfway point is Eller et al.'s (2018, Eqs 2.7-2.8, for
! the path-averaged conductance); Sabot et al. (2022, Eq. 16) use k at
! psi_leaf instead. The kcrit subtraction is Sabot et al.'s (Eller:
! k/k_max).
l_sox = som_profit_model == sox_opt_profit_model

! CAP (stomata_model = 10): A_n with the nonstomatal limitation (l_som_nsl,
! set with it), no stomatal cost.
l_cap = som_profit_model == cap_profit_model

! ProfitMax2 (stomata_model = 11): A_n (1 - E / E_crit) (Wang et al. 2020).
l_pm2 = som_profit_model == pm2_profit_model
l_abs = l_cmax .OR. l_cgain .OR. l_sox .OR. l_cap .OR. l_pm2
l_no_glcap = .FALSE.
b_curve(:) = conductance_b_pft(pft)
c_curve(:) = conductance_c_pft(pft)

DO j = 1, open_pts
  l = veg_index(open_index(j))

  ! The method, step by step. Each step is a contained routine below, which
  ! works on this loop's j and l and the host's variables.
  CALL point_setup()
  CALL range_ends( l_shut )
  IF ( l_shut ) CYCLE
  CALL find_edge()
  CALL scale_and_ends( l_shut )
  IF ( l_shut ) CYCLE
  CALL golden_max()
  CALL lambda_diag()
  CALL chosen_state()
  IF ( l_intact_search .AND. gl_g(j) > 0.0 ) CALL hold_leaf_psi()
  CALL stem_psi_of_choice()
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

  !---------------------------------------------------------------------------
  ! Set up one point: the leaf area of the call (lai_call, for the CMax and
  ! CGain costs, which are given per unit leaf area), the reference
  ! conductances of SOX_opt (the path with no flow) and CGain (the plant
  ! with no flow at zero water potential), the kcrit of the hydraulic cost,
  ! and the Ci range [ci_lo, ci_top]. ci_lo is ccp. With a gl cap, ci_top is
  ! 0.01 Pa below ca; without one, it is the top sample of the grid search,
  ! ca - (ca - ccp) / n_top.
  !---------------------------------------------------------------------------
  SUBROUTINE point_setup()
  IF ( PRESENT(psi_stem_g) )                                                   &
    psi_stem_g(j) = psi_root_zone(l) - grav_node_drop(l)
  IF ( PRESENT(lambda_g) ) lambda_g(j) = 0.0
  idx1(1) = open_index(j)
  g_off = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)

  ! (CMax, CGain, SOX and CAP work on net A, as Wolf et al. 2016 and Lu et
  ! al. 2020.)
  IF ( l_abs ) g_off = 0.0
  lai_call = 1.0

  ! The leaf area of the call from the intact plant (kmax_ref; = kmax
  ! without xylem impairment), so damage does not shrink the CMax or CGain
  ! cost, and without the kplant_model = 1 ground factor or the sapwood-area
  ! factor (l_kscale), which are not leaf area.
  IF ( l_abs ) lai_call = kmax_ref(l)                                          &
                          / MAX(kmax_pft(pft), TINY(1.0_real_jlslsm))
  IF ( l_abs .AND. l_kscale )                                                  &
    lai_call = lai_call / kplant_scale(l)

  ! The two-leaf sunlit and shaded leaves give their leaf area (their kmax
  ! share is N-weighted, not leaf area).
  IF ( l_abs .AND. PRESENT(lai_cost) ) lai_call = lai_cost(l)

  k0_cgain = 1.0
  k0_sox = 1.0
  kcrit_hc = kcrit(l)
  IF ( l_shared_stem_path ) kcrit_hc = shared_kcrit_plant(l)

  IF ( l_sox ) k0_sox = k_path_zero_flow(pft, l, kmax(l), psi_root_zone(l))

  ! (The plant's k_max: without the soil-to-root link, which would shrink
  ! the reference, and so the cost of a given embolism, as the soil dries.)
  IF ( l_cgain ) k0_cgain = MAX(k_path_zero_flow(pft, l, kmax(l),              &
                                0.0_real_jlslsm, l_plant_only = .TRUE.),       &
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
  END SUBROUTINE point_setup

  !---------------------------------------------------------------------------
  ! Step 1: evaluate the bottom (ccp) and top (ci_top) of the Ci range. The
  ! bottom gives the largest path conductance, which is the k_max of the
  ! hydraulic cost unless the cost uses the intact path with no flow
  ! (kl_hc_max: impairment with ximpair_hold_gs, the store's flow held
  ! fixed, or the shared root and stem). l_shut is .TRUE. when the point is
  ! closed: the bottom is infeasible, or A is 0 or less at the top. (That
  ! test uses A without the nonstomatal limitation, which could otherwise
  ! hide an open optimum.)
  !---------------------------------------------------------------------------
  SUBROUTINE range_ends( l_shut )
  LOGICAL, INTENT(OUT) :: l_shut

  l_shut = .FALSE.

  !---------------------------------------------------------------------------
  ! 1. Lower end (max k) and top of the range.
  !---------------------------------------------------------------------------
  ! Hold leaf psi: search on the intact path (as no impairment, at kmax_ref);
  ! the impaired path is used after the search (section 4).
  l_intact_search = l_xylem_impairment .AND.                                   &
                    ximpair_stomatal_response == ximpair_hold_psi
  CALL eval_ci(ci_lo)
  k0_ref = kl_u

  ! Impairment, hold gs: the cost is relative to the intact path at zero
  ! flow (as the flat search). Hold leaf psi: the intact path's (kl_u, the
  ! search being on it).
  IF ( l_xylem_impairment .AND. ximpair_stomatal_response == ximpair_hold_gs ) &
    k0_ref = kl_hc_max(l)

  ! The store's flow held fixed: the unshifted path at zero flow (see
  ! kl_hc_max in stom_opt_mod). The shared root and stem: the same.
  IF ( l_store_q_path .OR. l_shared_stem_path ) k0_ref = kl_hc_max(l)

  IF ( .NOT. ok_u ) THEN
    CALL set_closed()
    l_shut = .TRUE.
    RETURN
  END IF

  CALL store_best(-HUGE(1.0_real_jlslsm), ci_lo)
  e_ci = ci_lo
  e_al = al_u
  e_gl = gl_u
  e_el = el_u
  e_psi = psi_u
  e_kl = kl_u
  e_khc = khc_u
  mrg_e = margin()

  CALL eval_ci(ci_top)

  ! (al0_u = al_u without the nonstomatal limitation, which can make A at
  ! ci_top <= 0 although the optimum inside the range has A > 0.)
  IF ( al0_u + g_off <= 0.0 ) THEN
    CALL set_closed()
    l_shut = .TRUE.
    RETURN
  END IF
  END SUBROUTINE range_ends

  !---------------------------------------------------------------------------
  ! Step 2: find the upper edge of the feasible Ci range, which bounds the
  ! search for the best Ci, and hold it in e_*. If the top of the range is
  ! feasible, the edge is the top. Otherwise it is found from full
  ! evaluations (edge_by_margin).
  !---------------------------------------------------------------------------
  SUBROUTINE find_edge()
  !---------------------------------------------------------------------------
  ! 2. Upper edge of the feasible range, ci_b (e_*), unless the whole range
  !    is feasible.
  !---------------------------------------------------------------------------
  IF ( ok_u ) THEN
    e_ci = ci_top
    e_al = al_u
    e_gl = gl_u
    e_el = el_u
    e_psi = psi_u
    e_kl = kl_u
    e_khc = khc_u
  ELSE
    CALL edge_by_margin()
  END IF
  END SUBROUTINE find_edge

  !---------------------------------------------------------------------------
  ! Step 3: work out the gain scale A_max (max_al), which the scores use to
  ! scale A, and try both ends of the feasible range as candidates for the best
  ! Ci.
  !
  ! A_max is A at the hydraulic edge (k = kcrit, or ci_top), the largest A
  ! along the supply curve up to failure (Sperry et al. 2017). If the gl
  ! cap set the edge, the hydraulic edge is found beyond it
  ! (hydraulic_edge). With the nonstomatal limitation, A_max is A without
  ! the limitation. With the shared root and stem and the whole-plant gain
  ! (l_shared_gain_path), the other group of leaves' A_max is added.
  ! ProfitMax2 also needs E_crit (find_e_crit). l_shut is .TRUE. when
  ! A_max <= 0: opening gains no carbon, so the point is closed.
  !---------------------------------------------------------------------------
  SUBROUTINE scale_and_ends( l_shut )
  LOGICAL, INTENT(OUT) :: l_shut

  l_shut = .FALSE.

  ! Normalisation; no carbon benefit from opening => closed, as for the flat
  ! search (stom_opt_profit_max_select).
  !
  ! The gain is normalised by A at the hydraulic edge (k = kcrit, or ci_top),
  ! as Sperry et al. (2017): the largest A along the supply curve up to
  ! failure. The gl cap (som_gl_max, the soil-supply cap) only cuts off the
  ! search range. When the cap set the edge (its margin is the smaller), the
  ! hydraulic edge is found beyond it: A at the capped edge would inflate the
  ! gain fraction and move the optimum towards opening wherever a cap binds
  ! before the hydraulics.
  hy_ci = e_ci
  hy_al = e_al
  IF ( gl_max(l) > 0.0 .AND. e_ci < ci_top ) THEN
    IF ( (gl_max(l) - e_gl) / gl_max(l)                                        &
         < (e_kl - kcrit(l)) / kmax(l) ) CALL hydraulic_edge()
  END IF
  max_al = hy_al + g_off

  ! With the nonstomatal limitation, normalise by the unlimited A at the edge
  ! so the gain-to-cost weighting is as without it (and the limitation
  ! lowers the gain fraction).
  IF ( l_som_nsl ) max_al = photo_al(hy_ci) + g_off
  IF ( max_al <= 0.0 ) THEN
    CALL set_closed()
    l_shut = .TRUE.
    RETURN
  END IF

  ! The shared root and stem: this class's A at its edge for the other
  ! class's search, and with the whole-plant gain (l_shared_gain_path)
  ! the normalisation by both classes' A at their edges,
  !
  !   gain = (A + g_off) / (max_al + max_al_o),
  !
  ! max_al_o the other class's, held fixed.
  IF ( l_shared_stem_path ) shared_amax_self(l) = max_al
  IF ( l_shared_gain_path ) max_al = max_al + shared_amax_other(l)

  ! ProfitMax2: E_crit of the point (needs the edge state e_*; overwrites
  ! the *_u evaluation state, which is not used again before eval_ci).
  e_crit = 1.0
  IF ( l_pm2 ) CALL find_e_crit()

  ! Both ends of the feasible range are candidates.
  best_f = profit(best_al, best_khc, best_psi, best_el)
  IF ( profit(e_al, e_khc, e_psi, e_el) > best_f ) THEN
    best_f = profit(e_al, e_khc, e_psi, e_el)
    best_ci = e_ci
    best_al = e_al
    best_gl = e_gl
    best_el = e_el
    best_psi = e_psi
    best_kl = e_kl
    best_khc = e_khc
  END IF
  END SUBROUTINE scale_and_ends

  !---------------------------------------------------------------------------
  ! Find the edge of the Ci range set by the hydraulics alone (hy_ci, hy_al),
  ! beyond a gl cap that binds first, for the gain scale only. The edge search
  ! runs again from the capped edge to ci_top with the cap left out of the test
  ! (l_no_glcap). The capped edge e_* is put back afterwards, so the Ci range
  ! of the search does not change.
  !---------------------------------------------------------------------------
  SUBROUTINE hydraulic_edge()
  REAL(KIND=real_jlslsm) :: s_ci, s_al, s_gl, s_el, s_psi, s_kl, s_khc, s_mrg
                            ! The capped edge e_* and its margin, kept while
                            ! the hydraulic edge is found.

  s_ci  = e_ci
  s_al  = e_al
  s_gl  = e_gl
  s_el  = e_el
  s_psi = e_psi
  s_kl  = e_kl
  s_khc = e_khc
  s_mrg   = mrg_e
  l_no_glcap = .TRUE.

  CALL eval_ci(ci_top)
  IF ( ok_u ) THEN
    ! Hydraulically feasible up to ci_top.
    hy_ci = ci_top
    hy_al = al_u
  ELSE
    ! From the capped edge (hydraulically feasible) to ci_top.
    mrg_e = MAX((e_kl - kcrit(l)) / kmax(l), TINY(1.0_real_jlslsm))
    CALL edge_by_margin()
    hy_ci = e_ci
    hy_al = e_al
  END IF

  l_no_glcap = .FALSE.
  e_ci  = s_ci
  e_al  = s_al
  e_gl  = s_gl
  e_el  = s_el
  e_psi = s_psi
  e_kl  = s_kl
  e_khc = s_khc
  mrg_e   = s_mrg
  END SUBROUTINE hydraulic_edge

  !---------------------------------------------------------------------------
  ! Step 4: find the Ci of largest profit between ccp and the edge by
  ! golden-section search. Two candidate points sit at 0.382 and 0.618 of the
  ! bracket. Each step drops the part beyond the worse one, and the kept point
  ! falls on a golden-ratio position of the new bracket, so each step costs one
  ! evaluation. Where profit may have more than one peak (l_multi_peak), a
  ! coarse scan of n_scan points first narrows the bracket to the two
  ! neighbours of the best one.
  !---------------------------------------------------------------------------
  SUBROUTINE golden_max()
  !---------------------------------------------------------------------------
  ! 3. Golden-section search for the maximum profit on [ccp, ci_b], stopping
  !    after n_iter iterations or once the bracket is below
  !    (ci_b - ccp)/opt_resol. (Not Brent's parabolic steps: in low light
  !    profit is nearly flat over a wide Ci range and, at this precision, too
  !    noisy for them.)
  !---------------------------------------------------------------------------
  ! The bracket [a, b] runs from the lowest Ci (ccp) to the feasible edge
  ! e_ci. tol1 is the target bracket width.
  a = ci_lo
  b = e_ci
  tol1 = (b - a) / opt_resol

  ! Where profit need not have a single peak in Ci (a Weibull c < 1 on the
  ! path, Johnson-Berry photosynthesis, or the nonstomatal limitation), a
  ! coarse scan first: n_scan evenly spaced interior points, and the
  ! golden-section bracket is narrowed to the two neighbours of the best.
  ! The search then finds the peak the scan picked (store_best keeps the
  ! best Ci evaluated either way).
  IF ( tol1 > 0.0 .AND. l_multi_peak() ) THEN
    i_best = 0
    f_best_scan = -HUGE(1.0_real_jlslsm)
    DO i = 1, n_scan
      CALL eval_ci(a + REAL(i) * (b - a) / REAL(n_scan + 1))
      f_scan = profit_u()
      IF ( f_scan > f_best_scan ) THEN
        f_best_scan = f_scan
        i_best = i
      END IF
    END DO
    IF ( i_best > 0 ) THEN
      d_scan = (b - a) / REAL(n_scan + 1)
      b = a + REAL(i_best + 1) * d_scan
      a = a + REAL(i_best - 1) * d_scan
    END IF
  END IF

  IF ( tol1 > 0.0 ) THEN
    ! Two interior points at the golden-ratio positions (0.382 and 0.618 of
    ! the bracket), c_ci < d_ci, and their profits fc and fd.
    c_ci = b - golden_ratio * (b - a)
    d_ci = a + golden_ratio * (b - a)
    CALL eval_ci(c_ci)
    fc = profit_u()
    CALL eval_ci(d_ci)
    fd = profit_u()

    ! Golden section assumes one maximum on [a, b]. Monotone A and k do not
    ! give that by themselves: it holds when A is concave in Ci and the cost
    ! convex in E (Weibull c >= 1), and l_multi_peak marks the cases where
    ! it may not, which the coarse scan above handles. If it still fails
    ! (e.g. a very flat profit in low light), the search returns the best Ci
    ! it has evaluated (store_best keeps it).
    !
    ! Each step drops the part of the bracket beyond the worse point, which
    ! cannot hold the maximum when there is one peak. The kept interior
    ! point lands on a golden-ratio position of the new bracket, so only one
    ! new profit evaluation is needed per step, and the bracket shrinks by
    ! a factor 0.618 each time.
    DO i = 1, n_iter
      IF ( b - a <= tol1 ) EXIT

      IF ( fc >= fd ) THEN
        ! Maximum lies in [a, d_ci]: d_ci becomes the new upper end, c_ci the
        ! new upper interior point; evaluate a new lower interior point.
        b = d_ci
        d_ci = c_ci
        fd = fc
        c_ci = b - golden_ratio * (b - a)
        CALL eval_ci(c_ci)
        fc = profit_u()
      ELSE
        ! Maximum lies in [c_ci, b]: c_ci becomes the new lower end, d_ci the
        ! new lower interior point; evaluate a new upper interior point.
        a = c_ci
        c_ci = d_ci
        fc = fd
        d_ci = a + golden_ratio * (b - a)
        CALL eval_ci(d_ci)
        fd = profit_u()
      END IF
    END DO
  END IF
  END SUBROUTINE golden_max

  ! Say whether profit may have more than one peak in Ci at this point, so that
  ! golden_max scans the range coarsely first. It is .TRUE. with the
  ! nonstomatal limitation, with Johnson-Berry photosynthesis, or with a
  ! vulnerability curve of shape parameter c < 1 (the whole plant's, or any
  ! segment's with plant segments).
  LOGICAL FUNCTION l_multi_peak()
  l_multi_peak = l_som_nsl .OR. pft_photo_model == photo_johnson .OR.          &
                 conductance_c(l) < 1.0
  IF ( l_som_plant_segments )                                                  &
    l_multi_peak = l_multi_peak .OR. MINVAL(conductance_c_seg(pft,:)) < 1.0
  END FUNCTION l_multi_peak

  !---------------------------------------------------------------------------
  ! Work out, for output, the marginal water price at the chosen Ci,
  ! lambda_g = dA/dE (mol CO2 per mol H2O), from A and E at two Ci a small
  ! step (0.1 % of the feasible range) either side of it, kept inside the
  ! range. Diagnostic only: the chosen state is not changed.
  !---------------------------------------------------------------------------
  SUBROUTINE lambda_diag()
  ! Diagnostic dA/dE at the optimum: central difference over
  ! +-1e-3 of the range about best_ci, clipped to [ci_lo, e_ci] (one-sided
  ! at an end). Only *_u is changed, not the chosen state.
  IF ( PRESENT(lambda_g) .AND. e_ci > ci_lo ) THEN
    lam_ci1 = MAX(best_ci - 1.0e-3_real_jlslsm * (e_ci - ci_lo), ci_lo)
    lam_ci2 = MIN(best_ci + 1.0e-3_real_jlslsm * (e_ci - ci_lo), e_ci)
    CALL eval_ci(lam_ci1)
    lam_a1 = al_u
    lam_e1 = el_u
    CALL eval_ci(lam_ci2)
    IF ( ABS(el_u - lam_e1) > TINY(1.0_real_jlslsm) )                          &
      lambda_g(j) = (al_u - lam_a1) / (el_u - lam_e1)
  END IF
  END SUBROUTINE lambda_diag

  !---------------------------------------------------------------------------
  ! Copy the best state found to the outputs, with its carbon gain and
  ! hydraulic cost as fractions. For CMax, CGain and ProfitMax2 the cost is
  ! the carbon taken off A, over A_max; for SOX_opt it is 1 minus the factor
  ! on A; for CAP it is the loss of A to the nonstomatal limitation, over
  ! A_max. A stomatal conductance of 0 counts as closed (A = -rd).
  !---------------------------------------------------------------------------
  SUBROUTINE chosen_state()
  ci_g(j) = MAX(0.0, best_ci)
  al_g(j) = MAX(-rd(l), best_al)
  gl_g(j) = MAX(0.0, best_gl)

  ! gl = 0 means closed: al = -rd, as in set_closed.
  IF (gl_g(j) <= 0.0) al_g(j) = -rd(l)

  ! The Rd inside al_g: the candidate's own with l_leaf_temp_ci.
  IF ( PRESENT(rd_g) ) THEN
    rd_g(j) = rd(l)
    IF ( gl_g(j) > 0.0 .AND. best_al >= -rd(l) ) rd_g(j) = rd_at_ci(best_ci)
  END IF
  psi_g(j) = best_psi
  el_g(j) = best_el
  kl_g(j) = best_kl
  carbon_gain_g(j) = (best_al + g_off) / max_al

  ! (Guarded as the flat search: with impairment k0_ref is kl_hc_max, not
  ! taken from a sample.)
  hydraulic_cost_g(j) = (k0_ref - best_khc)                                    &
                        / MAX(k0_ref - kcrit_hc, TINY(1.0_real_jlslsm))

  ! CMax, CGain: the cost as a fraction of the maximum A on the range.
  ! ProfitMax2: Theta = A_n E / E_crit, as a fraction of max A.
  IF ( l_cmax .OR. l_cgain .OR. l_pm2 )                                        &
    hydraulic_cost_g(j) = (best_al - profit(best_al, best_khc, best_psi,       &
                                            best_el)) / max_al

  ! SOX: 1 - the factor on A (as the flat-grid SOX).
  IF ( l_sox ) hydraulic_cost_g(j) = 1.0 - sox_factor(best_psi)

  ! CAP: the nonstomatal loss of A as a fraction of the maximum A.
  IF ( l_cap ) hydraulic_cost_g(j) = (photo_al(best_ci) - best_al) / max_al
  END SUBROUTINE chosen_state

  !---------------------------------------------------------------------------
  ! Step 5, with xylem impairment and ximpair_hold_psi: narrow the stomata until
  ! the leaf water potential on the damaged path is back at psi_star, that of
  ! the optimum found on the intact path. The damaged path carries less water
  ! at psi_star. Ci is lowered from that optimum towards ccp until the impaired
  ! path gives psi_leaf = psi_star, by halving the bracket each step
  ! (bisection; psi_leaf falls as Ci rises). If even ccp leaves the leaf below
  ! psi_star, ccp is taken. If the result is infeasible on the impaired path,
  ! the feasible guess with the largest Ci is used, and with none the point is
  ! closed. The reported gain and cost are those of the intact-path optimum.
  !---------------------------------------------------------------------------
  SUBROUTINE hold_leaf_psi()
  !---------------------------------------------------------------------------
  ! 4. Hold leaf psi: the regulated psi_leaf is that of the intact optimum;
  !    the impaired path delivers less at it. Find Ci in [ci_lo, best_ci]
  !    whose transpiration on the impaired path gives psi_leaf = psi_star
  !    (psi falls as Ci, gl and E rise), by bisection. If even ci_lo is
  !    below psi_star, take ci_lo. The gain and cost are those of the
  !    intact-path optimum.
  !---------------------------------------------------------------------------
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
      ! the midpoint. Stop after n_hold_iter halvings or within
      ! hold_psi_tol (Pa) of the target.
      DO i = 1, n_hold_iter
        s_mid = 0.5 * (s_lo + s_hi)
        CALL eval_ci(s_mid)
        CALL keep_feasible()
        IF ( psi_u > psi_star ) THEN
          s_lo = s_mid
        ELSE
          s_hi = s_mid
        END IF

        IF ( ABS(psi_u - psi_star) < hold_psi_tol ) EXIT
      END DO
    END IF
  END IF

  ! On the impaired path the last evaluation can be infeasible (k <=
  ! kcrit, or gl > gl_max). Take instead the feasible evaluation with the
  ! largest Ci (feasibility only falls as Ci rises), else ci_lo if that is
  ! feasible on the impaired path (section 1 checked it on the intact
  ! one), else close the point.
  IF ( .NOT. ok_u .AND. .NOT. l_h_ok ) THEN
    CALL eval_ci(ci_lo)
    CALL keep_feasible()
  END IF

  IF ( .NOT. ok_u .AND. l_h_ok ) THEN
    last_ci = h_ci
    al_u = h_al
    gl_u = h_gl
    el_u = h_el
    psi_u = h_psi
    kl_u = h_kl
  END IF

  IF ( ok_u .OR. l_h_ok ) THEN
    ci_g(j) = MAX(0.0, last_ci)
    al_g(j) = MAX(-rd(l), al_u)
    gl_g(j) = MAX(0.0, gl_u)
    IF (gl_g(j) <= 0.0) al_g(j) = -rd(l)
    IF ( PRESENT(rd_g) ) THEN
      rd_g(j) = rd(l)
      IF ( gl_g(j) > 0.0 .AND. al_u >= -rd(l) ) rd_g(j) = rd_at_ci(last_ci)
    END IF
    psi_g(j) = psi_u
    el_g(j) = el_u
    kl_g(j) = kl_u
    carbon_gain_g(j) = cg_star
    hydraulic_cost_g(j) = hc_star
  ELSE
    CALL set_closed()
  END IF
  END SUBROUTINE hold_leaf_psi

  !---------------------------------------------------------------------------
  ! Work out the stem water potential at the chosen transpiration (psi_stem_g),
  ! on the impaired segmented path, for the record of stem damage. Only with
  ! xylem impairment, plant segments and open stomata; the state itself is
  ! already chosen.
  !---------------------------------------------------------------------------
  SUBROUTINE stem_psi_of_choice()
  ! Stem leaf-side psi at the chosen transpiration on the (impaired) segment
  ! chain, for the per-segment stem damage (diagnostic only: the state is
  ! already chosen).
  IF ( PRESENT(psi_stem_g) .AND. l_xylem_impairment .AND.                      &
       l_som_plant_segments .AND. gl_g(j) > 0.0 ) THEN
    el_s1(1,1) = el_g(j)
    ! The same path solve as the search (som_path_psi).
    CALL som_path_psi( pft, 1, land_pts, 1, veg_index, idx1, el_s1,            &
                       psi_root_zone, kmax_ref, kmax, kcrit, conductance_b,    &
                       conductance_c, psi_leaf_extreme, psi_root_extreme,      &
                       .TRUE., psi_s1, kl_s1, psi_stem = st_s1 )
    psi_stem_g(j) = st_s1(1,1)
  END IF
  END SUBROUTINE stem_psi_of_choice

  !---------------------------------------------------------------------------
  ! Work out the leaf state at one candidate Ci for point l, for the search to
  ! score: A, gl, E, the leaf water potential and the path conductances, worked
  ! out as stom_opt_mod_ci does for one sample, and whether the Ci is feasible
  ! (ok_u). With l_leaf_temp_ci the candidate has its own leaf temperature
  ! (state_at_leaf_temp). With the nonstomatal limitation (l_som_nsl), A is
  ! scaled by a factor of the leaf water potential, and that factor has to
  ! agree with the water potential the scaled A gives; it is found by a
  ! bracketed straight-line search (Illinois), set out in the body.
  !---------------------------------------------------------------------------
  SUBROUTINE eval_ci(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: ag0, f_lo, f_hi, f_new, g_lo, g_hi, g_new
  INTEGER :: it, side_f
                            ! Iteration count, and the end replaced last
                            ! (named so as not to hide the host's side).

  last_ci = ci

  ! l_leaf_temp_ci: A, gl and E at the candidate's own leaf temperature.
  l_ltc_pt = .FALSE.
  IF ( l_ltc ) l_ltc_pt = ltc_b(l) > 0.0
  IF ( l_ltc_pt ) THEN
    CALL state_at_leaf_temp(ci)
    CALL psi_k_at_e()
  ELSE
    al0_u = photo_al(ci)
    CALL hydraulic_state(ci, al0_u)
  END IF

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
  !-------------------------------------------------------------------------
  IF ( l_som_nsl ) THEN
    ! Nothing to do if the unlimited state is not limited (f = 1).
    IF ( nsl_factor(psi_u) < 1.0 ) THEN
      ! Unlimited gross photosynthesis, and the bracket for f:
      !
      !   f_hi = 1: g_hi = 1 - f(psi_leaf at full A) > 0
      !   f_lo = 0: no transpiration, psi_leaf = psi_root_zone,
      !             g_lo = -f(psi_root_zone) <= 0
      !
      ! With the store's flow held fixed (l_store_q_path) the leaf at no
      ! transpiration is not at psi_root_zone: the store keeps it wetter
      ! while it gives, so g_lo is taken at the leaf psi of the path at
      ! E = 0 (the closed state, A = -Rd).
      ag0  = al0_u + rd(l)
      f_hi = 1.0
      g_hi = 1.0 - nsl_factor(psi_u)
      f_lo = 0.0
      g_lo = -nsl_factor(psi_root_zone(l) - grav_node_drop(l))
      ! The same with the shared root and stem, where the other class's
      ! flow keeps the leaf below psi_root_zone.
      IF ( l_store_q_path .OR. l_shared_stem_path ) THEN
        CALL hydraulic_state(ci, -rd(l))
        g_lo = -nsl_factor(psi_u)
      END IF

      IF ( g_lo >= 0.0 ) THEN
        ! The leaf at no transpiration (the soil, or the store) alone stops
        ! photosynthesis (psi <= psi_nsl0).
        CALL hydraulic_state(ci, -rd(l))
      ELSE
        ! The Illinois method (a modified regula falsi), used for every root
        ! find in this module: it keeps a bracket, so it cannot diverge or
        ! leave the feasible range; it needs no derivative (the hydraulics
        ! give none cheaply); and it converges superlinearly where plain
        ! regula falsi would stall with one end fixed. side_f records which
        ! end was replaced last (-1 low, +1 high), for the Illinois step.
        side_f = 0
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
            f_lo = f_new
            g_lo = g_new
            IF ( side_f == -1 ) g_hi = 0.5 * g_hi
            side_f = -1
          ELSE
            f_hi = f_new
            g_hi = g_new
            IF ( side_f == 1 ) g_lo = 0.5 * g_lo
            side_f = 1
          END IF
        END DO

        ! (The state left by hydraulic_state is that of f_new.)
      END IF
    END IF
  END IF

  ! A Ci sample is feasible when the leaf is no wetter than the root zone,
  ! the path conductance is above its critical value kcrit, and the leaf
  ! conductance is within gl_max (when gl_max is set). With the store's flow
  ! held fixed (l_store_q_path) the test is not made: the store is a second
  ! source, wetter than the soil while it gives. (Where E >= q_f the chain
  ! carries flow to the leaf, so the leaf cannot be wetter than the soil
  ! anyway; a strict test there would reject samples near E = q_f on
  ! rounding and shut points.)
  ! (Not the gl cap while l_no_glcap: the hydraulic edge.)
  ok_u = ( psi_u <= psi_root_zone(l) + TINY(psi_root_zone(l))                  &
           .OR. l_store_q_path )                                               &
         .AND. kl_u > kcrit(l)                                                 &
         .AND. (gl_max(l) <= 0.0 .OR. l_no_glcap .OR. gl_u <= gl_max(l))      &
         .AND. gl_u < gl_unreach
  END SUBROUTINE eval_ci

  ! Work out how much water a candidate Ci costs: for net photosynthesis al at
  ! internal CO2 ci, the stomatal conductance that supplies it,
  ! gl = 1.6 al R T / (ca - ci) (through the leaf boundary layer too with
  ! l_som_coupled_e or r_bl > 0), the transpiration E it gives, and then the
  ! leaf water potential and conductances at that E (psi_k_at_e). E is
  ! floored at 0 only without l_som_coupled_e and with r_bl = 0.
  SUBROUTINE hydraulic_state(ci, al)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci, al

  al_u = al
  gl_u = ratio * (al_u * rmol * t_leaf(l)) / MAX(ca(l) - ci, 1.0e-2_real_jlslsm)
  ! (Coupled, or with a leaf boundary layer: the CO2 crosses the boundary
  ! layer too.)
  IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
    gl_u = gl_co2_bl(al_u, MAX(ca(l) - ci, 1.0e-2_real_jlslsm),               &
                     rmol * t_leaf(l), r_bl(l))
  END IF
  el_u = MAX(0.0, e_of_g(gl_u, dq(l), pstar(l), t_leaf(l)))
  IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
    el_u = e_of_g(g_eff_e(MAX(gl_u, 0.0_real_jlslsm), r_bl(l), r_ca(l),        &
                          g_oth(l)), dq(l), pstar(l), t_leaf(l))
  END IF
  CALL psi_k_at_e()
  END SUBROUTINE hydraulic_state

  !-------------------------------------------------------------------------
  ! Work out A, gl and E at internal CO2 ci with the leaf at its own
  ! temperature for this candidate (l_leaf_temp_ci), as a candidate that
  ! transpires more has a cooler leaf. The leaf temperature T has to satisfy
  ! both the leaf energy balance and the transpiration at T:
  !
  !   T     = ltc_a - ltc_b E            (K; more transpiration cools)
  !   E     = c_e dq(T)                  (mol H2O m-2 s-1)
  !   dq(T) = dq0 + ltc_s (T - T0)       (kg kg-1)
  !
  ! where c_e is the transpiration per unit humidity deficit through gl
  ! (through g_eff with l_som_coupled_e or r_bl > 0), and dq0 and T0 are
  ! the deficit and leaf temperature the call is given. For a given gl this
  ! gives
  !
  !   T = (ltc_a - ltc_b c_e (dq0 - ltc_s T0)) / (1 + ltc_b c_e ltc_s)
  !
  ! A, and so gl, depend on T through the photosynthesis parameters, each
  ! scaled from T0 by exp(ltc_dl (T - T0)). The loop evaluates A a fixed
  ! n_ltc_iter times, at T0 and then at n_ltc_iter - 1 updated temperatures,
  ! with no tolerance test, so that the number of updates, and so the
  ! profit, does not jump from one Ci to the next; the temperature worked
  ! out after the last evaluation is not used. T stays within dt_ltc of T0.
  ! E is floored at 0, so a leaf cooled below the dew point does not count
  ! as taking up water.
  !-------------------------------------------------------------------------
  SUBROUTINE state_at_leaf_temp(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: tl, tl_new, f(n_ltc_par), c_e
  INTEGER :: it
  INTEGER, PARAMETER :: n_ltc_iter = 4
  REAL(KIND=real_jlslsm), PARAMETER :: dt_ltc = 10.0
                            ! Iterations and the largest change from T0
                            ! (K).

  tl = t_leaf(l)
  DO it = 1, n_ltc_iter
    f(:) = EXP(ltc_dl(l,:) * (tl - t_leaf(l)))
    al_u = leaf_net_photo( pft, pft_photo_model, l_multilayer, ci,             &
                           vcmax(l) * f(1), kc(l) * f(3), ko(l) * f(4), oi(l), &
                           ccp(l) * f(5), km(l) * f(6), je(l) * f(2), acr(l),  &
                           apar(l), fapar_lf(l), ipar(l), je_ratio(l),         &
                           pstar(l), rd(l) * f(7) )
    gl_u = ratio * (al_u * rmol * tl) / MAX(ca(l) - ci, 1.0e-2_real_jlslsm)
    IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
      gl_u = gl_co2_bl(al_u, MAX(ca(l) - ci, 1.0e-2_real_jlslsm), rmol * tl,  &
                       r_bl(l))
    END IF
    IF ( l_som_coupled_e .OR. r_bl(l) > 0.0 ) THEN
      c_e = e_of_g(g_eff_e(MAX(gl_u, 0.0_real_jlslsm), r_bl(l), r_ca(l),       &
                           g_oth(l)), 1.0_real_jlslsm, pstar(l), tl)
    ELSE
      c_e = e_of_g(MAX(gl_u, 0.0_real_jlslsm), 1.0_real_jlslsm, pstar(l), tl)
    END IF
    tl_new = ( ltc_a(l) - ltc_b(l) * c_e * (dq(l) - ltc_s(l) * t_leaf(l)) )    &
             / ( 1.0 + ltc_b(l) * c_e * ltc_s(l) )
    tl_new = MIN(MAX(tl_new, t_leaf(l) - dt_ltc), t_leaf(l) + dt_ltc)
    IF ( it == n_ltc_iter ) EXIT
    tl = tl_new
  END DO

  ! The state at the last T at which A and gl were evaluated, and the Rd
  ! inside that A.
  al0_u = al_u
  rd_u  = rd(l) * f(7)

  ! E floored at 0 with or without coupled E, as in hydraulic_state: a leaf
  ! cooled below the air's dew point (dawn, dusk, net longwave loss) would
  ! give E < 0, a leaf psi above psi_root_zone, and so an infeasible Ci
  ! where E = 0 keeps it open.
  el_u = MAX(0.0, c_e * (dq(l) + ltc_s(l) * (tl - t_leaf(l))))
  END SUBROUTINE state_at_leaf_temp

  ! Work out the leaf water potential psi_u, path conductance kl_u and
  ! conductance for the hydraulic cost khc_u at transpiration el_u, on the path
  ! the search runs on: the impaired path, the intact path (during a
  ! hold-psi search) or the point's own vulnerability curve. khc_u is
  ! kl_u, except with impairment and ximpair_hold_gs, where it is the intact
  ! path's at the same water potentials. With the shared root and stem it is
  ! the whole plant's conductance at the whole canopy's flow. Kept apart from
  ! hydraulic_state so that find_e_crit can call it with E alone.
  SUBROUTINE psi_k_at_e()
  REAL(KIND=real_jlslsm) :: el1(1,1), psi1(1,1), kl1(1,1), khc1(1,1)

  IF ( l_xylem_impairment .AND. .NOT. l_intact_search ) THEN
    el1(1,1) = el_u
    CALL impaired_state( el1, psi1, kl1, khc1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  ELSE IF ( l_intact_search ) THEN
    ! Hold-psi search: the intact path (kmax_ref, PFT curve).
    el1(1,1) = el_u
    CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax_ref, kcrit, b_curve, c_curve,       &
                       psi1, kl1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  ELSE
    ! The point's own curve (conductance_b/c), as the flat search
    ! (stom_opt_mod_ci); the same as the PFT's unless a model moves b or c.
    el1(1,1) = el_u
    CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax, kcrit, conductance_b,              &
                       conductance_c, psi1, kl1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  END IF

  khc_u = kl_u
  IF ( l_xylem_impairment .AND. .NOT. l_intact_search .AND.                    &
       ximpair_stomatal_response == ximpair_hold_gs ) khc_u = khc1(1,1)

  ! The shared root and stem (l_shared_stem_path): the cost conductance is
  ! the whole plant's at the whole canopy's flow E + E_o, on the same
  ! basis as above (the search's path; the intact segments at the impaired
  ! water potentials with impairment and ximpair_hold_gs). The leaf psi
  ! and the conductance of the feasibility test stay the class's own.
  ! With a cost common to both classes, of their total flow, and their
  ! own gains, the passes (sunlit, then shaded, each with the other's flow
  ! held) seek the joint optimum, where the marginal gain of E is the same
  ! in both classes. (Not after a hold-psi search, section 4, where the
  ! cost is not used.)
  IF ( l_shared_stem_path .AND. .NOT. ( l_xylem_impairment .AND.               &
       .NOT. l_intact_search .AND.                                             &
       ximpair_stomatal_response == ximpair_hold_psi ) ) THEN
    el1(1,1) = el_u
    l_shared_whole_path = .TRUE.
    IF ( l_xylem_impairment .AND. .NOT. l_intact_search ) THEN
      CALL impaired_state( el1, psi1, kl1, khc1 )
      khc_u = kl1(1,1)
      IF ( ximpair_stomatal_response == ximpair_hold_gs ) khc_u = khc1(1,1)
    ELSE IF ( l_intact_search ) THEN
      CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,            &
                         psi_root_zone, kmax_ref, kcrit, b_curve, c_curve,     &
                         psi1, kl1 )
      khc_u = kl1(1,1)
    ELSE
      CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,            &
                         psi_root_zone, kmax, kcrit, conductance_b,            &
                         conductance_c, psi1, kl1 )
      khc_u = kl1(1,1)
    END IF
    l_shared_whole_path = .FALSE.
  END IF
  END SUBROUTINE psi_k_at_e

  !---------------------------------------------------------------------------
  ! Find E_crit (mol H2O m-2 s-1) for ProfitMax2: the transpiration at which
  ! the conductance of the path the search runs on (psi_k_at_e) falls to kcrit,
  !
  !   kl(E_crit) = kcrit
  !
  !   1. Start from the transpiration at the edge of the Ci range, e_el,
  !      where kl > kcrit, and double E until kl <= kcrit, so that E_crit
  !      lies between the last two values. The first upper value is at
  !      least kmax times 0.1 MPa, in case e_el is near 0.
  !   2. Halve that bracket (bisection) until its width is within
  !      ecrit_rtol of its upper end.
  !
  ! e_el is a lower limit on E_crit, since the edge is set by kcrit only
  ! when gl_max and ci_top do not bind first. E_crit is never below e_el,
  ! and is kept above 0 so that E / E_crit stays finite when there is no
  ! transpiration at the edge (the cost is then near 0).
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
  ! Work out, with xylem impairment, the leaf water potential psi1 and
  ! conductance kl1 on the impaired path at one transpiration el1 at point l,
  ! and the conductance for the hydraulic cost, khc1. khc1 is that of the
  ! intact path (kmax_ref and the PFT's curve) at the same water potentials,
  ! with the soil-to-root link in series:
  !
  !   khc1   = k_int(psi_leaf) / (1 + k_int(psi_in) / K_s)
  !   psi_in = psi_src - E / K_s
  !
  ! where k_int is the intact curve, K_s the soil-to-root conductance,
  ! psi_src the root-zone water potential and psi_in the root's. With plant
  ! segments the intact segments give khc1. The same path solve as
  ! stom_opt_mod_ci (som_path_psi), with the solver som_psi_solver sets.
  !---------------------------------------------------------------------------
  SUBROUTINE impaired_state( el1, psi1, kl1, khc1 )
  REAL(KIND=real_jlslsm), INTENT(IN)  :: el1(1,1)
  REAL(KIND=real_jlslsm), INTENT(OUT) :: psi1(1,1), kl1(1,1), khc1(1,1)

  ! The same path solve as stom_opt_mod_ci (som_path_psi). (The .TRUE. is
  ! l_xylem_impairment.)
  CALL som_path_psi( pft, 1, land_pts, 1, veg_index, idx1, el1,                &

                     psi_root_zone, kmax_ref, kmax, kcrit, conductance_b,      &
                     conductance_c, psi_leaf_extreme, psi_root_extreme,        &
                     .TRUE., psi1, kl1, leaf_k_intact = khc1 )
  END SUBROUTINE impaired_state

  ! Return the factor f that lowers gross photosynthesis as the leaf dries, at
  ! leaf water potential psi (Pa), for the nonstomatal limitation (l_som_nsl).
  ! It is Dewar et al. (2022) Eqn 3(b) with an onset,
  !
  !   f = 1 - (psi - psi_on) / (psi0 - psi_on)
  !
  ! clipped to [0, 1]. psi_on (psi_nsl_onset) is the water potential at
  ! which the limitation starts and psi0 (psi_nsl0) the one at which
  ! photosynthesis stops.
  !
  ! Reference:
  ! - Dewar, R., Holtta, T. and Salmon, Y. (2022), Exploring optimal
  !   stomatal control under alternative hypotheses for the regulation of
  !   plant sources and sinks. New Phytol, 233: 639-654.
  !   https://doi.org/10.1111/nph.17795
  REAL(KIND=real_jlslsm) FUNCTION nsl_factor(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  ! (The denominator is safe: pftparm_io checks psi_nsl0 < psi_nsl_onset
  ! when l_som_nsl is on.)
  nsl_factor = 1.0 - (psi - psi_nsl_onset(pft))                                &
                     / (psi_nsl0(pft) - psi_nsl_onset(pft))
  nsl_factor = MIN(MAX(nsl_factor, 0.0_real_jlslsm), 1.0_real_jlslsm)
  END FUNCTION nsl_factor

  ! Return the net photosynthesis (mol CO2 m-2 s-1) at internal CO2 ci (Pa) for
  ! point l, at the leaf temperature the call is given (leaf_net_photo).
  REAL(KIND=real_jlslsm) FUNCTION photo_al(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  photo_al = leaf_net_photo( pft, pft_photo_model, l_multilayer, ci,           &
                             vcmax(l), kc(l), ko(l), oi(l), ccp(l), km(l),     &
                             je(l), acr(l), apar(l), fapar_lf(l), ipar(l),     &
                             je_ratio(l), pstar(l), rd(l) )
  END FUNCTION photo_al

  !---------------------------------------------------------------------------
  ! Step 2 from full evaluations (any hydraulics): find where the
  ! feasibility margin (the function margin) crosses 0, between the
  ! feasible edge so far, e_ci, and the infeasible top of the range, ci_top.
  ! The margin is the smaller of the scaled distances to the two limits,
  !
  !   (kl - kcrit) / kmax   and   (gl_max - gl) / gl_max  (when gl_max set)
  !
  ! above 0 for a feasible Ci and below 0 for an infeasible one. The search
  ! keeps one end of each kind and places each guess where the straight
  ! line between the ends crosses 0 (regula falsi). If the same end moves
  ! twice in a row, the margin at the other end is halved so both ends keep
  ! moving (the Illinois form). It stops when the bracket is within
  ! edge_rtol of the distance to ca, or the feasible end is within
  ! edge_margin of a limit.
  !---------------------------------------------------------------------------
  SUBROUTINE edge_by_margin()
  INTEGER :: it

  ! Upper end of the bracket: ci_top (infeasible), and its margin.
  b_ci = ci_top
  mrg_b = margin()

  ! Illinois root find: a regula falsi step, bisection if it falls outside
  ! the bracket, and the end with the same sign replaced; if the same end
  ! is replaced twice in a row, the margin at the other end is halved so the
  ! bracket keeps shrinking from both sides. Each step is a full
  ! evaluation with the hydraulics. Stop when the bracket is within
  ! edge_rtol of the distance to ca, or once the feasible end is within
  ! edge_margin of a limit (close enough to the edge).
  side = 0
  DO it = 1, max_edge_iter
    IF ( b_ci - e_ci <= edge_rtol * (ca(l) - e_ci) .OR.                        &
         mrg_e < edge_margin ) EXIT

    ! Regula falsi step, or bisection if it falls outside the bracket.
    ci_try = (e_ci * mrg_b - b_ci * mrg_e) / (mrg_b - mrg_e)
    IF ( .NOT. ( ci_try > e_ci .AND. ci_try < b_ci ) )                         &
      ci_try = 0.5 * (e_ci + b_ci)
    CALL eval_ci(ci_try)
    mrg_u = margin()

    ! A feasible point becomes the new edge (with its full state); an
    ! infeasible one the new upper end. If the same end moves twice in a
    ! row, halve the margin at the other end (Illinois).
    IF ( ok_u ) THEN
      e_ci = ci_try
      e_al = al_u
      e_gl = gl_u
      e_el = el_u
      e_psi = psi_u
      e_kl = kl_u
      e_khc = khc_u
      mrg_e = mrg_u
      IF ( side == 1 ) mrg_b = 0.5 * mrg_b
      side = 1
    ELSE
      b_ci = ci_try
      mrg_b = mrg_u
      IF ( side == -1 ) mrg_e = 0.5 * mrg_e
      side = -1
    END IF
  END DO
  END SUBROUTINE edge_by_margin

  ! Return the feasibility margin of the latest evaluation, which
  ! edge_by_margin uses to find the edge of the feasible Ci range: the smaller
  ! of (kl - kcrit) / kmax and, when the gl cap applies,
  ! (gl_max - gl) / gl_max. It falls as Ci rises. It is forced above 0 for a
  ! feasible Ci and below 0 for an infeasible one, so its sign always agrees
  ! with ok_u.
  REAL(KIND=real_jlslsm) FUNCTION margin()
  margin = (kl_u - kcrit(l)) / kmax(l)
  IF ( gl_max(l) > 0.0 .AND. .NOT. l_no_glcap )                                &
    margin = MIN(margin, (gl_max(l) - gl_u) / gl_max(l))
  IF ( ok_u ) THEN
    margin = MAX(margin, TINY(1.0_real_jlslsm))
  ELSE
    margin = MIN(margin, -TINY(1.0_real_jlslsm))
  END IF
  END FUNCTION margin

  ! Return the score the search maximises over Ci, for a state with net
  ! photosynthesis al_in (A), conductance for the hydraulic cost kl_in (k),
  ! leaf water potential psi_in and transpiration el_in (E):
  !
  !   profit max   (A + g_off) / A_max - (k_max - k) / (k_max - k_crit)
  !                (Sperry et al. 2017; dimensionless)
  !   CMax         A - Theta(psi_leaf)            (Wolf et al. 2016, Eq. 8)
  !   CGain        A - varpi (k_max - k) / k_max  (Lu et al. 2020)
  !   SOX_opt      A f_SOX(psi_leaf)              (Eller et al. 2018)
  !   CAP          A, which carries the nonstomatal limitation
  !   ProfitMax2   A (1 - E / E_crit)             (Wang et al. 2020, Eqn 11)
  !
  ! For profit max, g_off is rd for the gross gain (l_som_gain_gross) and 0
  ! otherwise, A_max is the gain scale, k_max the reference conductance and
  ! k_crit the kcrit of the cost (the whole plant's with the shared root and
  ! stem). Theta is the CMax carbon cost (theta_cmax). For CGain, varpi
  ! (cgain_varpi) is the carbon cost of losing all the plant's conductance
  ! and k_max the plant's conductance at zero water potential. f_SOX is the
  ! SOX_opt factor (sox_factor), and E_crit the transpiration at which the
  ! path conductance falls to kcrit (find_e_crit). The last five scores are
  ! in mol CO2 m-2 s-1, on the basis of the call (whole canopy, the sunlit
  ! or shaded leaves, or a layer's leaves).
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
                     * MAX(k0_cgain - kl_in, 0.0_real_jlslsm) / k0_cgain
  ELSE
    ! Profit max (Sperry et al. 2017): relative carbon gain minus relative
    ! hydraulic cost, both 0-1. g_off = Rd with l_som_gain_gross (gross
    ! gain), else 0; max_al and k0_ref are the maxima over the Ci range.
    ! (k0_ref > kcrit at any feasible point; the guard is as for the
    ! diagnostic cost in stom_opt_bounded_search.) kcrit_hc: kcrit, or the
    ! whole plant's with the shared root and stem.
    profit = (al_in + g_off) / max_al - (k0_ref - kl_in)                       &
             / MAX(k0_ref - kcrit_hc, TINY(1.0_real_jlslsm))
  END IF
  END FUNCTION profit

  ! Return the SOX_opt factor on A at leaf water potential psi (Pa), which
  ! lowers A as the path loses conductance,
  !
  !   f_SOX = (k(psi_bar) - kcrit) / (k_ref - kcrit),   no lower than 0
  !
  ! where psi_bar is the mean of psi and the root-zone water potential,
  ! k(psi_bar) the conductance of the whole soil-to-leaf path with every
  ! part at psi_bar, and k_ref that conductance with no flow (psi_bar at the
  ! root-zone water potential).
  REAL(KIND=real_jlslsm) FUNCTION sox_factor(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  sox_factor = ( k_path_zero_flow(pft, l, kmax(l),                             &
                                   0.5 * (psi + psi_root_zone(l)))             &
                 - kcrit(l) ) / MAX(k0_sox - kcrit(l), TINY(1.0_real_jlslsm))
  sox_factor = MAX(sox_factor, 0.0_real_jlslsm)
  END FUNCTION sox_factor

  ! Return the CMax carbon cost of leaf water potential psi (Pa), from Anderegg
  ! et al. (2018) and Sabot et al. (2022, Eq. 12),
  !
  !   Theta = a/2 psi^2 + b |psi|
  !
  ! with psi in MPa (taken as 0 when above 0) and a and b (cmax_a, cmax_b)
  ! per unit leaf area, in umol m-2 s-1 MPa-2 and umol m-2 s-1 MPa-1.
  ! Returned in mol CO2 m-2 s-1 on the basis of this call (times lai_call).
  REAL(KIND=real_jlslsm) FUNCTION theta_cmax(psi)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi
  REAL(KIND=real_jlslsm) :: p_mpa
  p_mpa = ABS(MIN(psi, 0.0_real_jlslsm)) * 1.0e-6_real_jlslsm
  theta_cmax = lai_call * 1.0e-6_real_jlslsm                                   &
               * ( 0.5_real_jlslsm * cmax_a(pft) * p_mpa * p_mpa               &
                   + cmax_b(pft) * p_mpa )
  END FUNCTION theta_cmax

  ! Return the profit of the latest evaluation, or -HUGE if it is infeasible.
  ! If it is the best so far, the state is kept (store_best).
  REAL(KIND=real_jlslsm) FUNCTION profit_u()
  IF ( ok_u ) THEN
    profit_u = profit(al_u, khc_u, psi_u, el_u)
  ELSE
    profit_u = -HUGE(1.0_real_jlslsm)
  END IF

  IF ( profit_u > best_f ) CALL store_best(profit_u, last_ci)
  END FUNCTION profit_u

  ! Keep the latest evaluation as the best so far, with profit f at Ci ci:
  ! copy the state the last eval_ci call left (net photosynthesis, stomatal
  ! conductance, transpiration, leaf water potential and the two path
  ! conductances) into the best_* variables.
  SUBROUTINE store_best(f, ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: f, ci
  best_f = f
  best_ci = ci
  best_al = al_u
  best_gl = gl_u
  best_el = el_u
  best_psi = psi_u
  best_kl = kl_u
  best_khc = khc_u
  END SUBROUTINE store_best

  ! The dark respiration inside A at candidate ci: with l_leaf_temp_ci, Rd at
  ! the candidate's own leaf temperature (state_at_leaf_temp, the same
  ! evaluation the search made at ci); rd otherwise. The evaluation state
  ! (al_u, gl_u, el_u, al0_u) is kept.
  REAL(KIND=real_jlslsm) FUNCTION rd_at_ci(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  REAL(KIND=real_jlslsm) :: s_al, s_gl, s_el, s_al0
  rd_at_ci = rd(l)
  l_ltc_pt = .FALSE.
  IF ( l_ltc ) l_ltc_pt = ltc_b(l) > 0.0
  IF ( .NOT. l_ltc_pt ) RETURN
  s_al  = al_u
  s_gl  = gl_u
  s_el  = el_u
  s_al0 = al0_u
  CALL state_at_leaf_temp(ci)
  rd_at_ci = rd_u
  al_u  = s_al
  gl_u  = s_gl
  el_u  = s_el
  al0_u = s_al0
  END FUNCTION rd_at_ci

  !---------------------------------------------------------------------------
  ! Step 5 (hold_leaf_psi): keep the latest evaluation if it is feasible
  ! and has the largest Ci so far, as the fallback when the last guess is
  ! infeasible.
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

  !---------------------------------------------------------------------------
  ! Close the stomata at this point: Ci = ca, A = -rd, no conductance or
  ! transpiration, the leaf at the root-zone water potential, and a carbon
  ! gain and hydraulic cost of 0.
  !---------------------------------------------------------------------------
  SUBROUTINE set_closed()

  ! Closed stomata: zero flow, so the conductance is K(psi_root_zone) on the
  ! same curve as hydraulic_state, not 0 (0 would make PLC_pft report 100 %;
  ! see stom_opt_mod_ci).
  REAL(KIND=real_jlslsm) :: el0(1,1), psi0(1,1), kl0(1,1)
  REAL(KIND=real_jlslsm) :: khc0(1,1)

  IF ( l_xylem_impairment ) THEN
    el0(1,1) = 0.0
    CALL impaired_state( el0, psi0, kl0, khc0 )
  ELSE
    ! (The point's own curve, as eval_ci.)
    el0(1,1) = 0.0
    CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el0,              &
                       psi_root_zone, kmax, kcrit, conductance_b,              &
                       conductance_c, psi0, kl0 )
  END IF

  ci_g(j) = ca(l)
  al_g(j) = -rd(l)
  IF ( PRESENT(rd_g) ) rd_g(j) = rd(l)
  gl_g(j) = 0.0
  kl_g(j) = kl0(1,1)
  psi_g(j) = psi_root_zone(l) - grav_node_drop(l)
  el_g(j) = 0.0
  carbon_gain_g(j) = 0.0
  hydraulic_cost_g(j) = 0.0
  END SUBROUTINE set_closed

END SUBROUTINE stom_opt_bounded_search

!-----------------------------------------------------------------------------
! Set the stomatal opening of supply-loss stomata (stomata_model = 9), one open
! point at a time. The stomata follow the rule of Sperry et al. (2016, Eqn 5),
! and there is no optimisation. The leaf water potential may fall below the
! root zone's by the drop the unregulated demand would cause, scaled by the
! fraction of the path's conductance that demand would leave. The rule runs on
! the same plant hydraulics as profit max, and the demand comes from
! photosynthesis where Sperry et al. use a fixed maximum conductance.
!
!   1. Demand. A fixed ratio of internal to air CO2 (as Medlyn et al. 2011
!      with a fixed water price, so the humidity deficit acts only through
!      the hydraulics):
!
!        ci_d = sl_cica_well_watered ca
!        A_d  = A(ci_d)
!        gl_d = 1.6 A_d R T / (ca - ci_d)
!        E'   = gl_d D
!
!      E' is the transpiration with no regulation (mol H2O m-2 s-1) and D
!      the humidity deficit in mol m-3 (vpd / (R T), see e_of_g). A_d <= 0
!      closes the stomata.
!
!   2. Allowed drop in water potential. On the supply function of the path
!      (leaf water potential psi(E) and conductance k = -dE/dpsi_leaf),
!
!        P' = psi(E')
!        dP = (P0 - P') k(P') / k(P0)
!
!      with P0 = psi_root_zone, the leaf water potential with no flow. dP
!      first rises with E' and then falls. Past its maximum, or when E' is
!      beyond the flow at which k falls to k_fail, dP is held at that
!      maximum (Sperry et al. 2016), found by a coarse scan and then a
!      golden-section search.
!
!   3. Regulated transpiration E, the root of psi(E) = P0 - dP, by a
!      bracketed straight-line search (Illinois). Then gl = E / D, no more
!      than gl_d and no more than gl_max.
!
!   4. Ci for that gl, the root of
!
!        gl (ca - Ci) - 1.6 R T A(Ci) = 0,   Ci between ccp and ci_d
!
!      by the same kind of search.
!
! k_fail is kcrit, or sl_kfail_frac kmax when that is set. som_gl_max is
! not used: gl_max carries only the soil supply limit (l_som_supply_limit).
! Every search here keeps a bracket, so none can diverge near the failure
! point.
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
! MDK, 04/11/26
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
USE jules_vegetation_mod, ONLY: sl_kfail_frac
USE pftparm, ONLY: sl_cica_well_watered, conductance_b_pft, conductance_c_pft
USE jules_surface_mod, ONLY: ratio
USE c_rmol, ONLY: rmol
USE xylem_hydraulics_jls_mod, ONLY: psi_leaf_jls
USE plant_hydraulic_path_mod, ONLY: grav_node_drop

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

LOGICAL :: l_sat

REAL(KIND=real_jlslsm) :: b_curve(land_pts), c_curve(land_pts)
                            ! The PFT's vulnerability curve on land_pts, for
                            ! psi_leaf_jls (no xylem impairment here).

REAL(KIND=real_jlslsm) ::                                                      &
  rt, conv_e, k_fail, p0, k0, ci_d, al_d, gl_d, e_d, psi_d, k_d,               &
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

! Supply function from one-point calls of psi_leaf_jls (as
! stom_opt_bounded_search).
idx1(:) = 1
b_curve(:) = conductance_b_pft(pft)
c_curve(:) = conductance_c_pft(pft)

DO j = 1, open_pts
  l = veg_index(open_index(j))
  idx1(1) = open_index(j)
  rt = rmol * t_leaf(l)

  ! E (mol m-2 s-1) = conv_e * gl (m s-1), as stom_opt_mod_ci.
  conv_e = e_of_g(1.0_real_jlslsm, dq(l), pstar(l), t_leaf(l))
  p0 = psi_root_zone(l) - grav_node_drop(l)
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
        dp_max = f_scan
        i_best = i_scan
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
        gd_e = gc_e
        fgd = fgc
        gc_e = gb - golden_ratio * (gb - ga)
        fgc = dp_scored(gc_e)
      ELSE
        ga = gc_e
        gc_e = gd_e
        fgc = fgd
        gd_e = ga + golden_ratio * (gb - ga)
        fgd = dp_scored(gd_e)
      END IF
    END DO

    IF ( fgc > dp_max ) THEN
      e_max = gc_e
      dp_max = fgc
    END IF

    IF ( fgd > dp_max ) THEN
      e_max = gd_e
      dp_max = fgd
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
  ! 3. Regulated E, the root of
  !
  !      psi(E) = P0 - dP,
  !
  !    by Illinois on [0, xb], where psi(0) - psi_t = dP >= 0 and
  !    psi(xb) - psi_t <= 0 (dP <= dP'(xb)). (Illinois as in eval_ci.)
  !---------------------------------------------------------------------------
  psi_t = p0 - dp
  xa = 0.0
  fa = dp
  e_tol = e_rtol * xb
  CALL supply(xb, psi_u, k_u)
  fb = psi_u - psi_t
  e_reg = xb
  psi_reg = psi_u
  k_reg = k_u
  IF ( fb < 0.0 .AND. dp > 0.0 ) THEN
    side = 0
    DO it = 1, max_iter
      xc = (xa * fb - xb * fa) / (fb - fa)
      IF ( .NOT. ( xc > xa .AND. xc < xb ) ) xc = 0.5 * (xa + xb)
      CALL supply(xc, psi_u, k_u)
      fc = psi_u - psi_t
      e_reg = xc
      psi_reg = psi_u
      k_reg = k_u
      IF ( ABS(fc) <= e_rtol * dp .OR. xb - xa <= e_tol ) EXIT
      IF ( fc > 0.0 ) THEN
        xa = xc
        fa = fc
        IF ( side == 1 ) fb = 0.5 * fb
        side = 1
      ELSE
        xb = xc
        fb = fc
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
  ! 4. ci for the regulated gl, the root of
  !
  !      h(ci) = gl (ca - ci) - ratio R T A(ci) = 0
  !
  !    on [ccp, ci_d], h(ccp) = ratio R T rd > 0 and h(ci_d) =
  !    (gl - gl_d)(ca - ci_d) < 0, decreasing. Illinois (as in eval_ci).
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
        xa = xc
        fa = fc
        IF ( side == 1 ) fb = 0.5 * fb
        side = 1
      ELSE
        xb = xc
        fb = fc
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

  ! Return the supply function at transpiration e (mol H2O m-2 s-1): the leaf
  ! water potential psi (Pa) and the path conductance k = -dE/dpsi_leaf, from
  ! psi_leaf_jls on the PFT's vulnerability curve.
  SUBROUTINE supply(e, psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: e
  REAL(KIND=real_jlslsm), INTENT(OUT) :: psi, k
  REAL(KIND=real_jlslsm) :: el1(1,1), psi1(1,1), kl1(1,1)
  el1(1,1) = e
  CALL psi_leaf_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,                &
                     psi_root_zone, kmax, kcrit, b_curve, c_curve,             &
                     psi1, kl1 )
  psi = psi1(1,1)
  k = kl1(1,1)
  END SUBROUTINE supply



  ! Return the allowed drop in leaf water potential (Pa) for a state with leaf
  ! water potential psi and path conductance k (Sperry et al. 2016, Eqn 5):
  !
  !   dP = (P0 - psi) k / k0
  !
  ! with P0 the root-zone water potential and k0 the conductance with no
  ! flow.
  REAL(KIND=real_jlslsm) FUNCTION dp_of(psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: psi, k
  dp_of = MAX(p0 - psi, 0.0) * k / k0
  END FUNCTION dp_of

  ! Return the allowed drop dP at transpiration e, or -1 beyond the failure
  ! point (k <= k_fail), so the search for the largest dP never picks such a
  ! state.
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

  ! Copy a regulated state to the outputs of open point j, with Ci and gl
  ! no lower than 0 and net photosynthesis no lower than -rd.
  SUBROUTINE set_state(ci, al, gl, e, psi, k)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci, al, gl, e, psi, k
  ci_g(j) = MAX(0.0, ci)
  al_g(j) = MAX(-rd(l), al)
  gl_g(j) = MAX(0.0, gl)
  el_g(j) = e
  psi_g(j) = psi
  kl_g(j) = k
  END SUBROUTINE set_state

  ! Close the stomata at this point: Ci = ca, A = -rd, no conductance or
  ! transpiration, the leaf at the root-zone water potential, and a carbon
  ! gain and hydraulic cost of 0.
  SUBROUTINE set_closed()

  ! Closed stomata: zero flow, so the conductance is k at P0 (see
  ! stom_opt_bounded_search).
  REAL(KIND=real_jlslsm) :: psi0, k0c

  CALL supply(0.0_real_jlslsm, psi0, k0c)
  ci_g(j) = ca(l)
  al_g(j) = -rd(l)
  gl_g(j) = 0.0
  kl_g(j) = k0c
  psi_g(j) = psi_root_zone(l) - grav_node_drop(l)
  el_g(j) = 0.0
  carbon_gain_g(j) = 0.0
  hydraulic_cost_g(j) = 0.0
  END SUBROUTINE set_closed

  ! Return the net photosynthesis (mol CO2 m-2 s-1) at internal CO2 ci (Pa) for
  ! point l (leaf_net_photo).
  REAL(KIND=real_jlslsm) FUNCTION photo_al(ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci
  photo_al = leaf_net_photo( pft, pft_photo_model, l_multilayer, ci,           &
                             vcmax(l), kc(l), ko(l), oi(l), ccp(l), km(l),     &
                             je(l), acr(l), apar(l), fapar_lf(l), ipar(l),     &
                             je_ratio(l), pstar(l), rd(l) )
  END FUNCTION photo_al

END SUBROUTINE stom_opt_supply_loss

!-----------------------------------------------------------------------------
! Mark which Ci samples of one point are feasible, for the grid searches and
! stom_opt_profit_max_select to choose from: the path conductance is above
! kcrit, the leaf is no wetter than the root zone (when l_psi_test), gl is
! within gl_max (when gl_max > 0), and the leaf boundary layer can supply the
! CO2 (gl below gl_unreach, see gl_co2_bl).
!
! MDK, 25/10/26
!-----------------------------------------------------------------------------
PURE FUNCTION som_feasible( n, psi_s, kl_s, gl_s, psi_rz, kcr, gl_max,         &
                            l_psi_test ) RESULT( l_ok )
INTEGER, INTENT(IN) :: n
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_s(n), kl_s(n), gl_s(n), psi_rz,      &
                                      kcr, gl_max
LOGICAL, INTENT(IN) :: l_psi_test
LOGICAL :: l_ok(n)

l_ok(:) = kl_s(:) > kcr
IF ( l_psi_test ) l_ok(:) = l_ok(:) .AND. psi_s(:) <= psi_rz + TINY(psi_rz)
IF ( gl_max > 0.0 ) l_ok(:) = l_ok(:) .AND. gl_s(:) <= gl_max
! A Ci the leaf boundary layer cannot supply (gl_co2_bl).
l_ok(:) = l_ok(:) .AND. gl_s(:) < gl_unreach

END FUNCTION som_feasible

END MODULE stom_opt_jls_mod
