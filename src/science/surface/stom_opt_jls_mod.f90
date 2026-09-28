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

PRIVATE stom_opt_mod_ci, stom_opt_profit_max_select, stom_opt_golden_search
PUBLIC stom_opt_mod

CONTAINS

! *********************************************************************
! Contains routines to calculate the optimal stomatal conductance and
! net carbon uptake.
!
! References:
! TODO: Add reference to profit max paper.
! TODO: ADD reference to JULES paper.
! *********************************************************************

SUBROUTINE stom_opt_mod (                                                      &
! IN
        land_pts, som_base_parm, pft, open_pts, open_index,                    &
        pft_photo_model, veg_index,                                            &
        ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,           &
        km, dq, qs, je, t_leaf, je_ratio, fapar_lf, kmax, kcrit,               &
        gl_max, ipar, l_multilayer,                                            &
        kmax_ref, conductance_b, conductance_c,                                &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment_in,             &
! IN OUT
        rd,                                                                    &
! OUT
        ci, al, el, flux_o3, fo3, gl, psi_leaf,                                &
        carbon_gain_out, hydraulic_cost_out, leaf_k                            &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE jules_vegetation_mod, ONLY:                                                &
        som_base_parm_ci, som_base_parm_psi, som_n_sample,                     &
        profit_max_profit_model, SOX_profit_model, som_profit_model,          &
        som_ci_search_method, ci_search_flat,                                 &
        ci_search_golden, som_n_ci_prescan, som_n_ci_golden_iter,             &
        l_som_skip_search_wellwatered, som_hc_negligible_tol

USE pftparm, ONLY:                                                             &
        min_glw_pft, kcrit_fractional_loss

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls, leaf_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_conductance_impaired_jls

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft,                        &
                   pft_xylem_impairment_model
USE jules_vegetation_mod, ONLY: xylem_impairment_none
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
,veg_index(land_pts)
                            ! Index of vegetated points on the land grid.

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
,kmax(land_pts)                                                                &
,kcrit(land_pts)                                                               &
,gl_max(land_pts)                                                              &
                            ! Maximum stomatal conductance for H2O (m/s, same
                            ! basis as gl: per leaf area for multilayer,
                            ! canopy for big-leaf). Samples with gl above it
                            ! are infeasible. <= 0 disables the cap.
,ipar(land_pts)                                                                &
,kmax_ref(land_pts)                                                            &
                            ! Unimpaired maximum xylem conductance (m/s),
                            ! with the same scaling as kmax (kmax is the
                            ! impaired value when l_xylem_impairment).
,conductance_b(land_pts)                                                       &
                            ! Conductance parameter b for each land point
                            ! (Pa).
,conductance_c(land_pts)                                                       &
                            ! Conductance parameter c for each land point.
,psi_leaf_extreme(land_pts)                                                    &
                            ! Historic minimum leaf water potential (Pa).
,psi_root_extreme(land_pts)
                            ! Historic minimum root zone water potential (Pa).

!-----------------------------------------------------------------------------
! IN OUT real variables
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
 rd(land_pts)              ! Dark respiration (mol CO2/m2/s).
                             ! This is modified only if l_o3_damage=.TRUE..

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
,fo3(land_pts)                                                                 &
                            ! Ozone exposure factor.
,gl(land_pts)                                                                  &
                            ! Leaf conductance for H2O (m/s).
,psi_leaf(land_pts)                                                            &
                            ! Leaf water potential (Pa)
,leaf_k(land_pts)
                            ! Xylem conductance at leaf water potential (m/s)

! TEMPORARY: output variables for testing
REAL(KIND=real_jlslsm) ::                                                      &
 carbon_gain_out(land_pts)                                                     &
                            ! Carbon gain for each leaf state
,hydraulic_cost_out(land_pts)
                            ! Hydraulic cost for each leaf state

!-----------------------------------------------------------------------------
! Local scalar variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 optimal_index                                                                 &
                            ! Holds index of the optimal stomatal conductance
                            !  for each land point. Used by SOX_profit_model
                            !  and ci_search_flat's single-pass search.
,i,j,l                                                                         &
                            ! Iterators
,errcode
                            ! Error code to pass to ereport.

INTEGER ::                                                                     &
 optimal_index_flat(open_pts)
                            ! Index of the best sample for ci_search_flat,
                            ! from the shared stom_opt_profit_max_select
                            ! helper (replaces the bare "optimal_index"
                            ! scalar loop previously used inline here).

LOGICAL ::                                                                     &
 l_good_sample(som_n_sample, open_pts)                                        &
                            ! Feasibility mask for the single-pass
                            ! SOX_profit_model search.
,l_good_sample_flat(som_n_sample, open_pts)
                            ! Feasibility mask for ci_search_flat.

!-----------------------------------------------------------------------------
! Arrays containing results of sampling over the base parameter. Only need
! to do so for locations with open stomata.
!
! NOTE on som_ci_search_method (profit_max_profit_model only -
! SOX_profit_model is unaffected and always uses the flat single-pass
! search at som_n_sample resolution, using the "single-pass" arrays below):
!   ci_search_flat: unchanged flat grid, but now routed through the shared
!     stom_opt_profit_max_select helper (same arrays, same som_n_sample
!     resolution, no behaviour change from before this rework).
!   ci_search_golden: see stom_opt_golden_search - works on scalar
!     (open_pts)-shaped arrays, not a big sample grid, so needs no large
!     array declarations here.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 ci_lo(land_pts)                                                               &
                            ! Lower bound of the Ci sampling range for the
                            ! pass about to be run.
,ci_hi(land_pts)
                            ! Upper bound of the Ci sampling range for the
                            ! pass about to be run.

! -- Single-pass sample arrays (SOX_profit_model, always; ci_search_flat) --
REAL(KIND=real_jlslsm) ::                                                      &
 ci_sample(0:som_n_sample, open_pts)                                           &
                            ! Internal CO2 pressure (Pa).
,al_sample(0:som_n_sample, open_pts)                                           &
                            ! Net Leaf photosynthesis (mol CO2/m2/s).
,gl_sample(0:som_n_sample, open_pts)                                           &
                            ! Leaf conductance for H2O (m/s).
,kl_sample(0:som_n_sample, open_pts)                                           &
                            ! Xylem conductance at leaf (m/s).
,kl_hc_sample(0:som_n_sample, open_pts)                                        &
                            ! Xylem conductance at leaf used for the
                            ! hydraulic cost (m/s) - see stom_opt_mod_ci.
,kl_SOX(0:som_n_sample, open_pts)                                              &
                            ! Xylem conductance for SOX profit model.
                            !  Calculated from the average of the leaf and
                            !  root zone water potentials. (m/s)
,psi_sample(0:som_n_sample, open_pts)                                          &
                            ! Water potential at leaf (Pa)
,SOX_mean_psi(0:som_n_sample, open_pts)                                        &
                            ! Mean water potential between leaf and root zone
                            ! (Pa). Used in SOX profit model.
,el_sample(0:som_n_sample, open_pts)                                           &
                            ! Transpiration rate (mol H2O/m2/s)
,profit(0:som_n_sample, open_pts)                                              &
                            ! Profit for each leaf state
,carbon_gain(0:som_n_sample, open_pts)                                         &
                            ! Carbon gain for each leaf state
,hydraulic_cost(0:som_n_sample, open_pts)                                      &
                            ! Hydraulic cost for each leaf state
,max_al(open_pts)                                                              &
                            ! Maximum net leaf photosynthesis (mol CO2/m2/s)
,max_kl(open_pts)                                                              &
                            ! Maximum xylem conductance (m/s)
,kmax_open(open_pts)                                                          &
                            ! kmax gathered onto the open-point index, for
                            ! the SOX xylem_conductance_jls call.
,kcrit_open(open_pts)                                                          &
                            ! kcrit gathered onto the open-point index, for
                            ! the SOX xylem_conductance_jls call.
,b_open(open_pts)                                                              &
,c_open(open_pts)                                                              &
                            ! conductance_b/c gathered onto the open-point
                            ! index, for the SOX xylem_conductance_jls call.
,kl_hc_max(land_pts)                                                           &
                            ! Reference (maximum) xylem conductance for the
                            ! hydraulic cost with l_xylem_impairment: the
                            ! unimpaired conductance at the root zone water
                            ! potential. JBaguley
,b_ref(land_pts)                                                               &
,c_ref(land_pts)
                            ! Unimpaired PFT conductance_b/c.

! -- ci_search_golden outputs (scalar per open point) --
REAL(KIND=real_jlslsm) ::                                                      &
 ci_golden(open_pts)                                                           &
,al_golden(open_pts)                                                          &
,gl_golden(open_pts)                                                          &
,kl_golden(open_pts)                                                          &
,psi_golden(open_pts)                                                         &
,el_golden(open_pts)                                                          &
,carbon_gain_golden(open_pts)                                                 &
,hydraulic_cost_golden(open_pts)                                              &
,profit_golden(open_pts)

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
                            ! Number of points to run the flat search over:
                            ! all of open_index_search for ci_search_flat,
                            ! only golden's fallback points for
                            ! ci_search_golden.
,open_index_flat(open_pts)
                            ! Compressed subset of open_index for the flat
                            ! search (same convention as open_index_search).

REAL(KIND=real_jlslsm) :: hc_ref(open_pts)
                            ! Reference conductance for the fast path's
                            ! hydraulic cost: kmax, or kl_hc_max with
                            ! l_xylem_impairment.

LOGICAL :: l_golden_fallback(open_pts)
                            ! Points stom_opt_golden_search could not
                            ! resolve - see that routine.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
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
gl(:)      = min_glw_pft(pft)
psi_leaf(:)= psi_root_zone(:)
leaf_k(:)  = kmax

carbon_gain(:,:) = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:) = 0.0
kl_SOX(:,:) = 0.0

! If there are no land points with open stomata then no calculation is
! needed, other than (with xylem impairment) the leaf xylem conductance.
IF(0 == open_pts) THEN
  IF (l_xylem_impairment) THEN
    ! Calculate the xylem conductance at the leaf water potential.
    CALL leaf_conductance_impaired_jls( pft, land_pts, psi_leaf, kmax_ref,     &
                                        kmax, kcrit, conductance_b,            &
                                        conductance_c, psi_leaf_extreme,       &
                                        leaf_k )
  END IF
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
  CALL leaf_conductance_jls( pft, land_pts, psi_root_zone, kmax_ref, kcrit,    &
                             b_ref, c_ref, kl_hc_max )
END IF

! ----------------------------------------------------------------------------
!  Default outputs for points with closed stomata. Overwritten below for
!  points with open stomata.
! ----------------------------------------------------------------------------
carbon_gain_out(:) = 0.0
hydraulic_cost_out(:) = 0.0

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
      ! point anyway, just at O(som_n_sample) or O(golden) times the cost.
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
            l_multilayer,                                                     &
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
        IF (l_som_skip_search_wellwatered .AND.                                &
            psi_e_fp(1,j) <= psi_root_zone(l) + TINY(psi_root_zone(l)) .AND.   &
            kl_e_fp(1,j) > kcrit(l) .AND.                                      &
            (gl_max(l) <= 0.0 .OR. gl_e_fp(1,j) <= gl_max(l))) THEN
          hc_at_maxstress(j) = (hc_ref(j) - kl_hc_e_fp(1,j)) /                 &
                                MAX(hc_ref(j) - kcrit(l), TINY(1.0_real_jlslsm))
          l_fastpath(j) = hc_at_maxstress(j) < som_hc_negligible_tol
        END IF

        IF (l_fastpath(j)) THEN
          ci(l) = MAX(0.0, ci_e_fp(1,j))
          al(l) = MAX(0.0, al_e_fp(1,j))
          gl(l) = MAX(0.0, gl_e_fp(1,j))
          psi_leaf(l) = psi_e_fp(1,j)
          el(l) = el_e_fp(1,j)
          leaf_k(l) = kl_e_fp(1,j)
          carbon_gain_out(l) = 1.0_real_jlslsm
          hydraulic_cost_out(l) = hc_at_maxstress(j)
        ELSE
          open_pts_search = open_pts_search + 1
          ! NOTE: open_index_search must stay in the same veg_index-compressed
          ! convention as open_index (i.e. hold open_index(j), not l =
          ! veg_index(open_index(j))) - it is later passed as the open_index
          ! actual argument to stom_opt_mod_ci/stom_opt_golden_search, both of
          ! which apply their own veg_index(open_index(j)) remap internally.
          ! Storing the already-remapped land point l here would cause that
          ! remap to be applied twice.
          open_index_search(open_pts_search) = open_index(j)
        END IF
      END DO

      SELECT CASE (som_ci_search_method)

      CASE (ci_search_flat)
        ! Every point not resolved by the fast path goes to the flat search
        ! below.
        open_pts_flat = open_pts_search
        open_index_flat(1:open_pts_search) =                                   &
                                          open_index_search(1:open_pts_search)

      CASE (ci_search_golden)
        !-------------------------------------------------------------------
        ! Golden-section search - see stom_opt_golden_search for the
        ! algorithm and its caveats. Run only over the points the fast path
        ! above didn't already resolve. Points where the golden prescan
        ! could not resolve the feasible Ci region (l_golden_fallback, see
        ! stom_opt_golden_search) are handed on to the flat search below
        ! instead of taking golden's result.
        !-------------------------------------------------------------------
        open_pts_flat = 0
        IF (open_pts_search > 0) THEN
          CALL stom_opt_golden_search(                                        &
          ! IN
              land_pts, pft, open_pts_search, open_index_search,              &
              pft_photo_model, veg_index,                                     &
              rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,&
              km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,  &
              gl_max,                                                         &
              l_multilayer, som_n_ci_prescan, som_n_ci_golden_iter,           &
              kmax_ref, conductance_b, conductance_c,                         &
              psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,         &
              kl_hc_max,                                                      &
          ! OUT
              ci_golden, al_golden, gl_golden, kl_golden, psi_golden,         &
              el_golden, carbon_gain_golden, hydraulic_cost_golden,           &
              l_golden_fallback                                               &
                  )

          DO j = 1, open_pts_search
            IF (l_golden_fallback(j)) THEN
              open_pts_flat = open_pts_flat + 1
              open_index_flat(open_pts_flat) = open_index_search(j)
            ELSE
              l = veg_index(open_index_search(j))
              ci(l) = ci_golden(j)
              al(l) = al_golden(j)
              gl(l) = gl_golden(j)
              psi_leaf(l) = psi_golden(j)
              el(l) = el_golden(j)
              leaf_k(l) = kl_golden(j)
              carbon_gain_out(l) = carbon_gain_golden(j)
              hydraulic_cost_out(l) = hydraulic_cost_golden(j)
            END IF
          END DO
        END IF

      CASE DEFAULT
        errcode = 101  !  a hard error
        CALL ereport(RoutineName, errcode,                                     &
             'som_ci_search_method should be flat (1) or golden (2)')

      END SELECT ! som_ci_search_method

      !---------------------------------------------------------------------
      ! Flat grid: som_n_sample points over the full [ccp, ca] range, routed
      ! through the shared stom_opt_profit_max_select helper, run over
      ! open_index_flat - every remaining point for ci_search_flat, or only
      ! golden's fallback points for ci_search_golden.
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
            l_multilayer,                                                     &
            kmax_ref, conductance_b, conductance_c,                           &
            psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,           &
        ! OUT
            ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample,&
            kl_hc_sample                                                      &
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
          ! were never open; an open optimum is clipped at 0 as before.
          IF (optimal_index_flat(j) == 0) THEN
            al(l) = al_sample(0,j)
          ELSE
            al(l) = MAX(0.0, al_sample(optimal_index_flat(j),j))
          END IF
          gl(l) = MAX(0.0, gl_sample(optimal_index_flat(j),j))
          psi_leaf(l) = psi_sample(optimal_index_flat(j),j)
          el(l) = el_sample(optimal_index_flat(j),j)
          leaf_k(l) = kl_sample(optimal_index_flat(j),j)
          carbon_gain_out(l) = carbon_gain(optimal_index_flat(j),j)
          hydraulic_cost_out(l) = hydraulic_cost(optimal_index_flat(j),j)
        END DO
      END IF

    CASE (SOX_profit_model)
      !-------------------------------------------------------------------
      ! SOX_profit_model: unaffected by som_ci_search_method, always a
      ! single-pass search at som_n_sample resolution over [ccp, ca].
      !-------------------------------------------------------------------
      l_good_sample(:,:) = .TRUE.
      ci_lo(:) = MAX(ccp(:), 0.0)
      ci_hi(:) = ca(:)

      CALL stom_opt_mod_ci(                                                    &
      ! IN
          land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,     &
          som_n_sample, ci_lo, ci_hi,                                          &
          rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,     &
          km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,       &
          l_multilayer,                                                       &
          kmax_ref, conductance_b, conductance_c,                             &
          psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,             &
      ! OUT
          ci_sample, al_sample, gl_sample, kl_sample, psi_sample, el_sample,   &
          kl_hc_sample                                                         &
              )

      ! ------------------------------------------------------------------
      !  Make mask of those samples that are physically reasonable.
      ! ------------------------------------------------------------------
      DO j = 1, open_pts
        l = veg_index(open_index(j))

        ! Require psi <= psi_root_zone
        l_good_sample(:,j) = l_good_sample(:,j)                                &
            .and. psi_sample(1:,j) <= psi_root_zone(l) + TINY(psi_root_zone(l))

        ! Require k > k_crit
        l_good_sample(:,j) = l_good_sample(:,j)                                &
                             .and. kl_sample(1:,j) > kcrit(l)

        ! Require gl <= gl_max (when the cap is enabled)
        IF (gl_max(l) > 0.0) l_good_sample(:,j) = l_good_sample(:,j)           &
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
        ! Sox uses the xylem conductance for the average of the
        !  leaf and root zone water potentials.
        DO j = 1, open_pts
          l = veg_index(open_index(j))

          SOX_mean_psi(:,j) = 0.5 * (psi_sample(:,j) + psi_root_zone(l))
          kmax_open(j) = kmax(l)
          kcrit_open(j) = kcrit(l)
          b_open(j) = conductance_b(l)
          c_open(j) = conductance_c(l)
        END DO

        ! See the note above the leaf_psi_jls call in stom_opt_mod_ci: pass
        ! the (1:som_n_sample,:) sections so the actual/dummy shapes match
        ! exactly and avoid the same sequence-association misalignment. Entry
        ! 0 of kl_SOX is left at the 0.0 it was initialised to above.
        CALL xylem_conductance_jls(pft, som_n_sample, open_pts,                &
                                   SOX_mean_psi(1:som_n_sample,:),             &
                                   kmax_open, kcrit_open, b_open, c_open,      &
                                   kl_SOX(1:som_n_sample,:))

        ! Get the maximum xylem conductance for each land point.
        max_kl = MAXVAL(kl_SOX(1:,:), MASK = l_good_sample, DIM = 1)
      END IF

      ! CG = An
      carbon_gain(:,:) = al_sample(:,:)

      ! Loop over land points with open stomata
      DO j = 1, open_pts
          l = veg_index(open_index(j))
          ! HC = 1 - k/ki_max
          hydraulic_cost(:,j) = 1 - (kl_SOX(:,j)-kcrit(l))                     &
                                    /(max_kl(j)-kcrit(l))

          ! Profit = CG*(1-HC)
          profit(:,j) = carbon_gain(:,j)*(1-hydraulic_cost(:,j))
      END DO ! i Land points

      ! ------------------------------------------------------------------
      !  Identify maximum profit and set leaf properties to those of the
      !   optimal state. See the NOTE above stom_opt_profit_max_select for
      !   what the zeroth (closed-stomata) sample represents.
      ! ------------------------------------------------------------------
      DO j = 1, open_pts
        l = veg_index(open_index(j))

        optimal_index = MAXLOC(profit(1:,j), DIM = 1, MASK = l_good_sample(:,j))

        ci(l) = MAX(0.0, ci_sample(optimal_index,j))
        IF (optimal_index == 0) THEN
          al(l) = al_sample(0,j)
        ELSE
          al(l) = MAX(0.0, al_sample(optimal_index,j))
        END IF
        gl(l) = MAX(0.0, gl_sample(optimal_index,j))
        psi_leaf(l) = psi_sample(optimal_index,j)
        el(l) = el_sample(optimal_index,j)
        leaf_k(l) = kl_sample(optimal_index,j)

        carbon_gain_out(l) = carbon_gain(optimal_index,j)
        hydraulic_cost_out(l) = hydraulic_cost(optimal_index,j)
      END DO

    CASE DEFAULT
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
                 'profit_model should be profit_max or SOX')

    END SELECT ! som_profit_model

  CASE (som_base_parm_psi)
    WRITE(*,*) "psi not implimented yet"

  CASE DEFAULT
    errcode = 101  !  a hard error
    CALL ereport(RoutineName, errcode,                                         &
               'som_base_parm should be 1 or 2')
END SELECT

IF (l_xylem_impairment) THEN
  ! Calculate the xylem conductance at the leaf water potential, for all
  ! points (open and closed), on the impaired vulnerability curve.
  ! JBaguley
  CALL leaf_conductance_impaired_jls( pft, land_pts, psi_leaf, kmax_ref,       &
                                      kmax, kcrit, conductance_b,              &
                                      conductance_c, psi_leaf_extreme,         &
                                      leaf_k )
END IF

!TODO: Redetermine if stomata are open or closed?

!TODO: Calculate ozone exposure here. Create subroutine using code in
!TODO:  leaf_jls_mod.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod

!-----------------------------------------------------------------------------
! Stomatal optimisation model that iterates over leaf ci
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_mod_ci(                                                    &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        n_sample, ci_lo, ci_hi,                                                &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp,              &
        pstar, km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax,         &
        kcrit, l_multilayer,                                                  &
        kmax_ref, conductance_b, conductance_c,                               &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,               &
! OUT
        ci_sample, al_sample, gl_sample, kl_sample, psi_sample,el_sample,      &
        kl_hc_sample                                                           &
)

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls, xylem_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_psi_impaired_jls

USE jules_vegetation_mod, ONLY:                                                &
        photo_collatz, photo_farquhar, CW_conductance,                         &
        SOX_conductance, som_psi_aprox_method

USE pftparm, ONLY:                                                             &
        leaf_crit, c3, alpha, pft_conductance_model,                           &
        conductance_b_pft, conductance_c_pft

USE jules_surface_mod, ONLY: fwe_c3, fwe_c4
USE jules_surface_mod, ONLY: beta1, beta2, ratio, ratio_o3

USE planet_constants_mod, ONLY: repsilon, one_minus_epsilon

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
                            ! resolutions/ranges (flat/golden Ci searches
                            ! in stom_opt_mod) rather than being fixed to
                            ! the module-level som_n_sample.

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
                            ! Xylem conductance at leaf (m/s).
,psi_sample(0:n_sample, open_pts)                                            &
                            ! Water potential at leaf (Pa)
,el_sample(0:n_sample, open_pts)                                             &
                            ! Transpiration rate (mol H2O/m2/s)
,kl_hc_sample(0:n_sample, open_pts)
                            ! Xylem conductance at leaf used for the
                            ! hydraulic cost (m/s): kl_sample, or with
                            ! l_xylem_impairment the conductance at
                            ! psi_sample on the unimpaired PFT curve.

!-----------------------------------------------------------------------------
! Local integer variables.
!-----------------------------------------------------------------------------
INTEGER ::                                                                     &
 i,j,k,l                                                                       &
                            ! Internal iterators
,errcode
                            ! Error code to pass to ereport.

!-----------------------------------------------------------------------------
! Local real variables.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) ::                                                      &
 wcarb_sample(0:n_sample, open_pts)                                          &
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
,kl_sample_new, psi_sample_new                                                 &
                            ! Temporary variables for use when estimating
                            ! leaf water potential from transpiration rate.
,vpd(land_pts)                                                                 &
                            ! Vapour pressure deficit (Pa)
,vpd_factor                                                                    &
                            ! Factor used to convert specific humidity to
                            ! vapour pressure.
,glco2_sample(0:n_sample, open_pts)                                        &
                            ! Leaf conductance for CO2 (m/s)
,conv(land_pts)                                                                &
                            ! Factor for converting mol/m3 into Pa (J/m3)
,conductance_conversion_limit                                                  &
                            ! Limit for the change in conductance when
                            ! estimating leaf water potential from
                            ! transpiration rate.
,kmax_open(open_pts), kcrit_open(open_pts), b_open(open_pts), c_open(open_pts)
                            ! Unimpaired PFT curve gathered onto the
                            ! open-point index, for kl_hc_sample.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_MOD_CI'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Setup parameters.
! NOTE: The zeroth entry to the sample array for each land point corresponds
!        to closed stomata and is only used when determining the optimal leaf
!        state in stom_opt_mod above. This is why all loops are over the
!        range 1:n_sample, NOT 0:n_sample.
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
!$OMP        n_sample)
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
      ! scaling for multi-layer calls, 29 Apr, MGDK
      IF (l_multilayer) THEN
          wlite_sample(:,j) = wlite_sample(:,j) / apar(l) * fapar_lf(l) * ipar(l)
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
!$OMP        pft,wexpt_sample,pstar,alpha,fwe_c4,ci_sample,acr,n_sample)
    DO j = 1,open_pts
      l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l)
      wlite_sample(:,j) = alpha(pft) * acr(l)
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))
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
!$OMP        ci_sample,ccp,km,je,n_sample)
  DO j = 1,open_pts
    l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l) * ( ci_sample(:,j) - ccp(l) )               &
                          / ( ci_sample(:,j) + km(l) )
      wlite_sample(:,j) = je(l) / 4.0 * ( ci_sample(:,j) - ccp(l) )            &
                          / ( ci_sample(:,j) + 2.0 * ccp(l) )
      wlite_sample(:,j) = MAX(wlite_sample(:,j), TINY(1.0e0))

      ! scaling 29 Apr, MGDK
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
! Calculate net photosynthetic rate, mol CO2 /m2/s.
!  al = wl - rd
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
DO j = 1,open_pts
  l = veg_index(open_index(j))

  al_sample(:,j) = wl_sample(:,j) - rd(l)
END DO

! ----------------------------------------------------------------------------
! Calculate stomatal conductance in units of m/s.
!  gsw = (1/1.6)*(A/(ca-ci))
! ----------------------------------------------------------------------------
!TODO: add open mp to loop
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
  ! to (or exactly at) ca(l) for a *narrow* range - as ci_search_golden's
  ! converging bracket legitimately can - this margin can shrink far more
  ! than the flat single-pass grid ever produced, overflowing this division
  ! to Infinity in practice.
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
!TODO: add open mp to loop

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
IF (l_xylem_impairment) THEN
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

! Reset entry 0 to the closed-stomata state. The whole-column (:,j)
! photosynthesis/conductance calculations above also evaluated entry 0 at
! ci = ca, leaving al_sample(0,j) at the net assimilation rate for fully open
! stomata and gl_sample(0,j) = al * conv / MAX(ca - ci, 1e-2) - i.e. near the
! maximum photosynthesis and a huge conductance. So whenever stom_opt_mod
! fell back to "closed" (no feasible sample, e.g. under severe drought),
! it actually reported maximal carbon uptake and gc of O(1-10) m/s, exactly
! when the hydraulics said the stomata should shut. Closed stomata: no
! conductance or transpiration, net assimilation = -rd, leaf at the root
! zone water potential.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  ci_sample(0,j)  = ca(l)
  al_sample(0,j)  = -rd(l)
  gl_sample(0,j)  = 0.0
  el_sample(0,j)  = 0.0
  kl_sample(0,j)  = 0.0
  kl_hc_sample(0,j) = 0.0
  psi_sample(0,j) = psi_root_zone(l)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_mod_ci

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
! Used by ci_search_flat's single pass to select the profit-maximising
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
, open_pts                                                                    &
, open_index(land_pts)                                                        &
, veg_index(land_pts)                                                         &
, n_sample

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  al_sample(0:n_sample, open_pts)                                             &
, kl_sample(0:n_sample, open_pts)                                             &
, psi_sample(0:n_sample, open_pts)                                            &
, gl_sample(0:n_sample, open_pts)                                             &
, psi_root_zone(land_pts)                                                     &
, kcrit(land_pts)                                                             &
, gl_max(land_pts)                                                            &
, kl_hc_sample(0:n_sample, open_pts)                                          &
                            ! Conductance used for the hydraulic cost - see
                            ! stom_opt_mod_ci.
, kl_hc_max(land_pts)                                                          &
                            ! Reference conductance for the hydraulic cost
                            ! with l_xylem_impairment - see stom_opt_mod.
, rd(land_pts)
                            ! Dark respiration; added back to An for the
                            ! gross-A gain (l_som_gain_gross).

LOGICAL, INTENT(IN) :: l_xylem_impairment

LOGICAL, INTENT(OUT) :: l_good_sample(n_sample, open_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  carbon_gain(0:n_sample, open_pts)                                           &
, hydraulic_cost(0:n_sample, open_pts)                                        &
, profit(0:n_sample, open_pts)

INTEGER, INTENT(OUT) :: optimal_index(open_pts)

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
!    CG = An/max(An)
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

    ! HC = (ki_max - k)/(ki_max - k_crit). max_kl(j) > kcrit(l) is
    ! guaranteed here because l_good_sample already required
    ! kl_sample > kcrit(l) for at least the sample max_kl was taken from,
    ! so this denominator can be small under severe stress but never zero.
    hydraulic_cost(:,j) =   (max_kl(j) - kl_hc_sample(:,j))                  &
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
! Golden-section search for the profit-maximising Ci, over [ccp, ca], for
! profit_max_profit_model. profit(ci) = CG(ci) - HC(ci) is expected to rise
! then fall (possibly via a hard cliff to infeasible where kl drops below
! kcrit, rather than a smooth interior decline) as ci increases - i.e.
! "unimodal" in the weak sense golden-section needs. This assumption has no
! safety net of its own (unlike a grid search, which samples the whole
! domain and so cannot miss an entire alternative region): a genuinely
! bimodal or non-monotonic profit(ci) would make this converge confidently
! to the wrong answer with no way to detect it. Mitigated here by:
!   (a) tracking the best profit actually seen across every point sampled
!       (prescan and golden-section alike), not just trusting the final
!       bracket, so a stray better sample earlier in the search still wins;
!   (b) the som_n_ci_prescan coarse prescan itself, which - being a full
!       flat pass, just at low resolution - would need to be actively
!       misled by aliasing to miss a materially better distant region.
! Uses n_prescan + 2 + n_iter total Ci evaluations (each evaluation = one
! stom_opt_mod_ci call at n_sample=1, batched across all open points).
!-----------------------------------------------------------------------------
SUBROUTINE stom_opt_golden_search(                                            &
! IN
        land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,       &
        rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,       &
        km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,         &
        gl_max, l_multilayer, n_prescan, n_iter,                              &
        kmax_ref, conductance_b, conductance_c,                               &
        psi_leaf_extreme, psi_root_extreme, l_xylem_impairment, kl_hc_max,    &
! OUT
        ci_g, al_g, gl_g, kl_g, psi_g, el_g, carbon_gain_g, hydraulic_cost_g, &
        l_fallback                                                            &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: l_som_gain_gross

INTEGER, INTENT(IN) ::                                                         &
  land_pts, pft, open_pts, open_index(land_pts), pft_photo_model,             &
  veg_index(land_pts), n_prescan, n_iter

REAL(KIND=real_jlslsm), INTENT(IN) ::                                         &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),        &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts), &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), qs(land_pts),  &
  je(land_pts), t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts),    &
  ipar(land_pts), kmax(land_pts), kcrit(land_pts), gl_max(land_pts),        &
  kmax_ref(land_pts), conductance_b(land_pts), conductance_c(land_pts),      &
  psi_leaf_extreme(land_pts), psi_root_extreme(land_pts), kl_hc_max(land_pts)
                            ! See stom_opt_mod.

LOGICAL, INTENT(IN) :: l_multilayer, l_xylem_impairment

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                        &
  ci_g(open_pts), al_g(open_pts), gl_g(open_pts), kl_g(open_pts),            &
  psi_g(open_pts), el_g(open_pts), carbon_gain_g(open_pts),                  &
  hydraulic_cost_g(open_pts)

LOGICAL, INTENT(OUT) :: l_fallback(open_pts)
                            ! .TRUE. where the prescan found fewer than
                            ! min_good_prescan feasible samples with
                            ! positive net photosynthesis: the feasible Ci
                            ! region is then too narrow for the prescan to
                            ! resolve, max_al_ref is unreliable (or <= 0,
                            ! which would force every golden evaluation to
                            ! -HUGE and close the stomata), and the caller
                            ! should use the flat search for this point
                            ! instead of the ci_g/al_g/... values returned
                            ! here.

! Local
INTEGER, PARAMETER :: min_good_prescan = 2
REAL(KIND=real_jlslsm), PARAMETER :: golden_ratio = 0.6180339887498949_real_jlslsm

INTEGER :: i, j, l, errcode

REAL(KIND=real_jlslsm) :: ci_lo(land_pts), ci_hi(land_pts)

! -- Prescan (flat, coarse, only used for max_al_ref/max_kl_ref) --
REAL(KIND=real_jlslsm) ::                                                     &
  ci_p(0:n_prescan, open_pts), al_p(0:n_prescan, open_pts),                  &
  gl_p(0:n_prescan, open_pts), kl_p(0:n_prescan, open_pts),                  &
  psi_p(0:n_prescan, open_pts), el_p(0:n_prescan, open_pts),                &
  kl_hc_p(0:n_prescan, open_pts)
LOGICAL :: good_p(n_prescan, open_pts)
REAL(KIND=real_jlslsm) :: max_al_ref(open_pts), max_kl_ref(open_pts)
REAL(KIND=real_jlslsm) :: g_off(open_pts)
                            ! Gain offset: 0 (net An) or Rd (gross A).

! -- Single-point evaluation scratch (n_sample=1 each call) --
REAL(KIND=real_jlslsm) ::                                                     &
  ci_e(0:1, open_pts), al_e(0:1, open_pts), gl_e(0:1, open_pts),             &
  kl_e(0:1, open_pts), psi_e(0:1, open_pts), el_e(0:1, open_pts),          &
  kl_hc_e(0:1, open_pts)

! -- Golden-section bracket state, per open point --
REAL(KIND=real_jlslsm) ::                                                     &
  a_pt(open_pts), b_pt(open_pts), c_pt(open_pts), d_pt(open_pts),            &
  fc(open_pts), fd(open_pts), f_new(open_pts)
REAL(KIND=real_jlslsm) ::                                                     &
  al_c(open_pts), gl_c(open_pts), kl_c(open_pts), psi_c(open_pts),           &
  el_c(open_pts)
REAL(KIND=real_jlslsm) ::                                                     &
  al_d(open_pts), gl_d(open_pts), kl_d(open_pts), psi_d(open_pts),           &
  el_d(open_pts)
REAL(KIND=real_jlslsm) ::                                                     &
  best_f(open_pts), best_ci(open_pts), best_al(open_pts), best_gl(open_pts), &
  best_kl(open_pts), best_psi(open_pts), best_el(open_pts),                 &
  best_kl_hc(open_pts)
                            ! Conductance used for the hydraulic cost (see
                            ! stom_opt_mod_ci) at the best point seen.
LOGICAL :: update_c(open_pts), feasible(open_pts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='STOM_OPT_GOLDEN_SEARCH'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
! Stage 0: coarse prescan over the full range, purely to get max_al_ref/
! max_kl_ref for the CG/HC normalisation (golden-section itself doesn't
! naturally produce these the way a batch grid search does).
!-----------------------------------------------------------------------------
ci_lo(:) = MAX(ccp(:), 0.0)
ci_hi(:) = ca(:)

CALL stom_opt_mod_ci(                                                          &
    land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,          &
    n_prescan, ci_lo, ci_hi,                                                  &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,          &
    km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,            &
    l_multilayer,                                                            &
    kmax_ref, conductance_b, conductance_c,                                  &
    psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                  &
    ci_p, al_p, gl_p, kl_p, psi_p, el_p, kl_hc_p                              &
        )

good_p(:,:) = .TRUE.
DO j = 1, open_pts
  l = veg_index(open_index(j))
  good_p(:,j) = good_p(:,j)                                                  &
          .AND. psi_p(1:,j) <= psi_root_zone(l) + TINY(psi_root_zone(l))     &
          .AND. kl_p(1:,j) > kcrit(l)
  IF (gl_max(l) > 0.0) good_p(:,j) = good_p(:,j) .AND. gl_p(1:,j) <= gl_max(l)
END DO

! gain offset: 0 (net An) or Rd (gross A), see l_som_gain_gross
DO j = 1, open_pts
  l = veg_index(open_index(j))
  g_off(j) = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)
  max_al_ref(j) = MAXVAL(al_p(1:,j) + g_off(j), MASK = good_p(:,j))
END DO
IF (l_xylem_impairment) THEN
  ! See the matching note in stom_opt_profit_max_select. JBaguley
  DO j = 1, open_pts
    max_kl_ref(j) = kl_hc_max(veg_index(open_index(j)))
  END DO
ELSE
  max_kl_ref = MAXVAL(kl_hc_p(1:,:), MASK = good_p, DIM = 1)
END IF

! Hand a point back to the flat search when the feasible Ci region is too
! narrow for the prescan to resolve (see the l_fallback declaration). Under
! severe water stress the feasible region is a sliver just above ccp; if it
! is narrower than ~1/n_prescan of [ccp, ca] the only feasible prescan
! sample is at ccp itself, where net photosynthesis is <= 0, and golden
! would otherwise close the stomata where a finer grid finds an open
! optimum. Points whose ccp sample is itself infeasible (kl <= kcrit even
! at zero transpiration) are not handed back - nothing in [ccp, ca] can be
! feasible there, so the flat search would also close the stomata.
DO j = 1, open_pts
  l_fallback(j) = good_p(1,j) .AND.                                          &
                  COUNT(good_p(:,j) .AND. al_p(1:,j) > 0.0) < min_good_prescan
END DO

! Seed "best seen" from the prescan itself, in case golden-section's own
! bracket never revisits the true optimum (see the module-level note above
! this subroutine) - this also correctly seeds the closed-stomata fallback
! (best_f = -HUGE) for any point with no feasible prescan sample or
! max_al_ref <= 0.
best_f(:) = -HUGE(1.0_real_jlslsm)
DO j = 1, open_pts
  l = veg_index(open_index(j))
  IF (max_al_ref(j) > 0.0) THEN
    DO i = 1, n_prescan
      IF (good_p(i,j)) THEN
        f_new(j) = (al_p(i,j) + g_off(j))/max_al_ref(j)                      &
                 - (max_kl_ref(j)-kl_hc_p(i,j))/(max_kl_ref(j)-kcrit(l))
        IF (f_new(j) > best_f(j)) THEN
          best_f(j) = f_new(j)
          best_ci(j) = ci_p(i,j); best_al(j) = al_p(i,j)
          best_gl(j) = gl_p(i,j); best_kl(j) = kl_p(i,j)
          best_psi(j) = psi_p(i,j); best_el(j) = el_p(i,j)
          best_kl_hc(j) = kl_hc_p(i,j)
        END IF
      END IF
    END DO
  END IF
END DO

!-----------------------------------------------------------------------------
! Stage 1: initialise the golden-section bracket [a,b] = [ccp, ca] and its
! two interior points c < d, and evaluate both (2 separate n_sample=1
! calls, since the ci_sample stepping formula can't place two arbitrary
! points c/d in a single n_sample=2 call - see the note on stom_opt_mod_ci).
!-----------------------------------------------------------------------------
DO j = 1, open_pts
  l = veg_index(open_index(j))
  a_pt(j) = MAX(ccp(l), 0.0)
  b_pt(j) = ca(l)
  c_pt(j) = b_pt(j) - golden_ratio * (b_pt(j) - a_pt(j))
  d_pt(j) = a_pt(j) + golden_ratio * (b_pt(j) - a_pt(j))
  ci_lo(l) = c_pt(j)
  ci_hi(l) = ca(l)
END DO

CALL stom_opt_mod_ci(                                                          &
    land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,          &
    1, ci_lo, ci_hi,                                                          &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,          &
    km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,            &
    l_multilayer,                                                            &
    kmax_ref, conductance_b, conductance_c,                                  &
    psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                  &
    ci_e, al_e, gl_e, kl_e, psi_e, el_e, kl_hc_e                              &
        )

DO j = 1, open_pts
  l = veg_index(open_index(j))
  al_c(j) = al_e(1,j); gl_c(j) = gl_e(1,j); kl_c(j) = kl_e(1,j)
  psi_c(j) = psi_e(1,j); el_c(j) = el_e(1,j)
  feasible(j) = (psi_c(j) <= psi_root_zone(l) + TINY(psi_root_zone(l)))       &
               .AND. (kl_c(j) > kcrit(l))                                     &
               .AND. (gl_max(l) <= 0.0 .OR. gl_c(j) <= gl_max(l))
  IF (feasible(j) .AND. max_al_ref(j) > 0.0) THEN
    fc(j) = (al_c(j) + g_off(j))/max_al_ref(j)                               &
           - (max_kl_ref(j)-kl_hc_e(1,j))/(max_kl_ref(j)-kcrit(l))
  ELSE
    fc(j) = -HUGE(1.0_real_jlslsm)
  END IF
  IF (fc(j) > best_f(j)) THEN
    best_f(j) = fc(j); best_ci(j) = c_pt(j); best_al(j) = al_c(j)
    best_gl(j) = gl_c(j); best_kl(j) = kl_c(j); best_psi(j) = psi_c(j)
    best_el(j) = el_c(j); best_kl_hc(j) = kl_hc_e(1,j)
  END IF
  ci_lo(l) = d_pt(j)
END DO

CALL stom_opt_mod_ci(                                                          &
    land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,          &
    1, ci_lo, ci_hi,                                                          &
    rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,          &
    km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,            &
    l_multilayer,                                                            &
    kmax_ref, conductance_b, conductance_c,                                  &
    psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                  &
    ci_e, al_e, gl_e, kl_e, psi_e, el_e, kl_hc_e                              &
        )

DO j = 1, open_pts
  l = veg_index(open_index(j))
  al_d(j) = al_e(1,j); gl_d(j) = gl_e(1,j); kl_d(j) = kl_e(1,j)
  psi_d(j) = psi_e(1,j); el_d(j) = el_e(1,j)
  feasible(j) = (psi_d(j) <= psi_root_zone(l) + TINY(psi_root_zone(l)))       &
               .AND. (kl_d(j) > kcrit(l))                                     &
               .AND. (gl_max(l) <= 0.0 .OR. gl_d(j) <= gl_max(l))
  IF (feasible(j) .AND. max_al_ref(j) > 0.0) THEN
    fd(j) = (al_d(j) + g_off(j))/max_al_ref(j)                               &
           - (max_kl_ref(j)-kl_hc_e(1,j))/(max_kl_ref(j)-kcrit(l))
  ELSE
    fd(j) = -HUGE(1.0_real_jlslsm)
  END IF
  IF (fd(j) > best_f(j)) THEN
    best_f(j) = fd(j); best_ci(j) = d_pt(j); best_al(j) = al_d(j)
    best_gl(j) = gl_d(j); best_kl(j) = kl_d(j); best_psi(j) = psi_d(j)
    best_el(j) = el_d(j); best_kl_hc(j) = kl_hc_e(1,j)
  END IF
END DO

!-----------------------------------------------------------------------------
! Stage 2: golden-section iterations. Each iteration evaluates exactly one
! new point per open point (whichever of c/d that point needs refreshed,
! decided independently per point by its own fc vs fd comparison), batched
! into a single stom_opt_mod_ci call across all open points.
!-----------------------------------------------------------------------------
DO i = 1, n_iter

  DO j = 1, open_pts
    l = veg_index(open_index(j))
    ! NOTE: ">=" not ">". Under real water stress the feasible Ci region
    ! (kl > kcrit) is a narrow low-Ci sliver with everything above it
    ! infeasible (profit = -HUGE), so both c and d can easily land in the
    ! infeasible zone early on - a genuine tie (fc == fd == -HUGE). ">"
    ! would then take the ELSE branch below, which shrinks the bracket
    ! towards *higher* Ci - moving further from the feasible region, since
    ! feasibility here always sits on the low-Ci side (lower Ci means less
    ! transpiration demand, hence higher kl, for any given psi_root_zone).
    ! ">=" makes a tie shrink towards *lower* Ci instead, which is the
    ! correct direction whenever both probes are infeasible, and is a
    ! no-op (arbitrary but harmless) tie-break in the well-watered case
    ! where a true numerical tie between two feasible, similar-profit
    ! points can also occur (see the low-stress "plateau" discussion).
    update_c(j) = fc(j) >= fd(j)
    IF (update_c(j)) THEN
      b_pt(j) = d_pt(j)
      d_pt(j) = c_pt(j); fd(j) = fc(j)
      al_d(j) = al_c(j); gl_d(j) = gl_c(j); kl_d(j) = kl_c(j)
      psi_d(j) = psi_c(j); el_d(j) = el_c(j)
      c_pt(j) = b_pt(j) - golden_ratio * (b_pt(j) - a_pt(j))
      ci_lo(l) = c_pt(j)
    ELSE
      a_pt(j) = c_pt(j)
      c_pt(j) = d_pt(j); fc(j) = fd(j)
      al_c(j) = al_d(j); gl_c(j) = gl_d(j); kl_c(j) = kl_d(j)
      psi_c(j) = psi_d(j); el_c(j) = el_d(j)
      d_pt(j) = a_pt(j) + golden_ratio * (b_pt(j) - a_pt(j))
      ci_lo(l) = d_pt(j)
    END IF
    ci_hi(l) = ca(l)
  END DO

  CALL stom_opt_mod_ci(                                                        &
      land_pts, pft, open_pts, open_index, pft_photo_model, veg_index,        &
      1, ci_lo, ci_hi,                                                        &
      rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,        &
      km, dq, qs, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,          &
      l_multilayer,                                                          &
      kmax_ref, conductance_b, conductance_c,                                &
      psi_leaf_extreme, psi_root_extreme, l_xylem_impairment,                &
      ci_e, al_e, gl_e, kl_e, psi_e, el_e, kl_hc_e                            &
          )

  DO j = 1, open_pts
    l = veg_index(open_index(j))
    feasible(j) = (psi_e(1,j) <= psi_root_zone(l) + TINY(psi_root_zone(l)))   &
                 .AND. (kl_e(1,j) > kcrit(l))                                 &
                 .AND. (gl_max(l) <= 0.0 .OR. gl_e(1,j) <= gl_max(l))
    IF (feasible(j) .AND. max_al_ref(j) > 0.0) THEN
      f_new(j) = (al_e(1,j) + g_off(j))/max_al_ref(j)                        &
               - (max_kl_ref(j)-kl_hc_e(1,j))/(max_kl_ref(j)-kcrit(l))
    ELSE
      f_new(j) = -HUGE(1.0_real_jlslsm)
    END IF

    IF (update_c(j)) THEN
      fc(j) = f_new(j); al_c(j) = al_e(1,j); gl_c(j) = gl_e(1,j)
      kl_c(j) = kl_e(1,j); psi_c(j) = psi_e(1,j); el_c(j) = el_e(1,j)
    ELSE
      fd(j) = f_new(j); al_d(j) = al_e(1,j); gl_d(j) = gl_e(1,j)
      kl_d(j) = kl_e(1,j); psi_d(j) = psi_e(1,j); el_d(j) = el_e(1,j)
    END IF

    IF (f_new(j) > best_f(j)) THEN
      best_f(j) = f_new(j); best_ci(j) = ci_e(1,j); best_al(j) = al_e(1,j)
      best_gl(j) = gl_e(1,j); best_kl(j) = kl_e(1,j)
      best_psi(j) = psi_e(1,j); best_el(j) = el_e(1,j)
      best_kl_hc(j) = kl_hc_e(1,j)
    END IF
  END DO

END DO

!-----------------------------------------------------------------------------
! Final selection: the best point seen anywhere in the search (prescan,
! initial bracket, or golden-section iterations), or the closed-stomata
! state if nothing was ever feasible - matching the same fallback
! philosophy as stom_opt_mod's other search methods (see the note above
! stom_opt_profit_max_select).
!-----------------------------------------------------------------------------
DO j = 1, open_pts
  l = veg_index(open_index(j))
  IF (best_f(j) > -HUGE(1.0_real_jlslsm)) THEN
    ci_g(j) = MAX(0.0, best_ci(j))
    al_g(j) = MAX(0.0, best_al(j))
    gl_g(j) = MAX(0.0, best_gl(j))
    psi_g(j) = best_psi(j)
    el_g(j) = best_el(j)
    kl_g(j) = best_kl(j)
    carbon_gain_g(j) = (best_al(j) + g_off(j)) / max_al_ref(j)
    hydraulic_cost_g(j) = (max_kl_ref(j) - best_kl_hc(j))                    &
                         / (max_kl_ref(j) - kcrit(l))
  ELSE
    ci_g(j) = ca(l)
    al_g(j) = -rd(l)
    gl_g(j) = 0.0
    kl_g(j) = 0.0
    psi_g(j) = psi_root_zone(l)
    el_g(j) = 0.0
    carbon_gain_g(j) = 0.0
    hydraulic_cost_g(j) = 0.0
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stom_opt_golden_search

END MODULE stom_opt_jls_mod
