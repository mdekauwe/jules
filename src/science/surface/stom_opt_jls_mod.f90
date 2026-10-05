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

PRIVATE stom_opt_mod_ci, stom_opt_profit_max_select, stom_opt_bounded_search
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
        carbon_gain_out, hydraulic_cost_out, leaf_k,                           &
! OUT (optional)
        psi_stem                                                               &
)

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE jules_vegetation_mod, ONLY:                                                &
        som_base_parm_ci, som_base_parm_psi, som_n_sample,                     &
        profit_max_profit_model, SOX_profit_model, som_profit_model,          &
        som_ci_search, som_ci_bounded, som_n_ci_golden_iter,                  &
        l_som_skip_search_wellwatered, som_hc_negligible_tol, l_som_nsl

USE pftparm, ONLY:                                                             &
        min_gl_pft, kcrit_fractional_loss

USE xylem_hydraulics_jls_mod, ONLY: xylem_conductance_jls, leaf_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_conductance_impaired_jls

USE pftparm, ONLY: conductance_b_pft, conductance_c_pft,                        &
                   pft_xylem_impairment_model
USE jules_vegetation_mod, ONLY: xylem_impairment_none, l_som_plant_segments
USE pftparm, ONLY: seg_kfac, conductance_b_seg, conductance_c_seg
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
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem(land_pts)
                            ! Stem water potential (Pa): the outlet of the
                            ! stem segment of the chosen sample, with
                            ! l_som_plant_segments and xylem impairment (flat
                            ! search); psi_root_zone otherwise.
REAL(KIND=real_jlslsm) ::                                                      &
 psi_stem_sample(0:som_n_sample, open_pts)
                            ! Stem water potential of each Ci sample.

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
                            !  for each land point. Used by SOX_profit_model.
,i,j,l,iseg                                                                    &
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
                            ! SOX_profit_model search.
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

! -- Flat-search sample arrays (SOX_profit_model, and profit_max_profit_model
! -- with som_ci_search = 1) --
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

! -- stom_opt_bounded_search outputs (scalar per open point) --
REAL(KIND=real_jlslsm) ::                                                      &
 ci_bnd(open_pts)                                                              &
,al_bnd(open_pts)                                                              &
,gl_bnd(open_pts)                                                              &
,kl_bnd(open_pts)                                                              &
,psi_bnd(open_pts)                                                             &
,el_bnd(open_pts)                                                              &
,carbon_gain_bnd(open_pts)                                                     &
,hydraulic_cost_bnd(open_pts)

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
leaf_k(:)  = kmax
IF (PRESENT(psi_stem)) psi_stem(:) = psi_root_zone(:)

carbon_gain(:,:) = 0.0
hydraulic_cost(:,:) = 0.0
profit(:,:) = 0.0
kl_SOX(:,:) = 0.0

! If there are no land points with open stomata then no calculation is
! needed, other than (with xylem impairment) the leaf xylem conductance.
IF(0 == open_pts) THEN
  IF (l_xylem_impairment .AND. .NOT. l_som_plant_segments) THEN
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
  IF (l_som_plant_segments) THEN
    ! Segments in series: the intact whole-plant conductance with no flow,
    ! 1 / sum_s 1 / k_s(psi_root_zone).
    DO l = 1, land_pts
      kl_hc_max(l) = 0.0
      DO iseg = 1, 3
        kl_hc_max(l) = kl_hc_max(l) + 1.0 / MAX(kmax_ref(l) * seg_kfac(pft,iseg) &
                       * EXP( -(psi_root_zone(l) / conductance_b_seg(pft,iseg)) &
                              **conductance_c_seg(pft,iseg) ), TINY(1.0_real_jlslsm))
      END DO
      kl_hc_max(l) = 1.0 / kl_hc_max(l)
    END DO
  ELSE
    CALL leaf_conductance_jls( pft, land_pts, psi_root_zone, kmax_ref, kcrit,  &
                               b_ref, c_ref, kl_hc_max )
  END IF
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
          CALL stom_opt_bounded_search(                                       &
          ! IN
              land_pts, pft, open_pts_search, open_index_search,              &
              pft_photo_model, veg_index,                                     &
              rd, ca, psi_root_zone, acr, apar, oi, vcmax, kc, ko, ccp, pstar,&
              km, dq, je, t_leaf, je_ratio, fapar_lf, ipar, kmax, kcrit,      &
              gl_max, l_multilayer, som_n_sample, som_n_ci_golden_iter,       &
          ! OUT
              ci_bnd, al_bnd, gl_bnd, kl_bnd, psi_bnd,                        &
              el_bnd, carbon_gain_bnd, hydraulic_cost_bnd                     &
                  )

          DO j = 1, open_pts_search
            l = veg_index(open_index_search(j))
            ci(l) = ci_bnd(j)
            al(l) = al_bnd(j)
            gl(l) = gl_bnd(j)
            psi_leaf(l) = psi_bnd(j)
            el(l) = el_bnd(j)
            leaf_k(l) = kl_bnd(j)
            carbon_gain_out(l) = carbon_gain_bnd(j)
            hydraulic_cost_out(l) = hydraulic_cost_bnd(j)
          END DO
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
            l_multilayer,                                                     &
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

    CASE (SOX_profit_model)
      !-------------------------------------------------------------------
      ! SOX_profit_model: unaffected by som_ci_search, always a
      ! single-pass search at som_n_sample resolution over [ccp, ca].
      !-------------------------------------------------------------------
      carbon_gain(:,:) = 0.0
      hydraulic_cost(:,:) = 0.0
      profit(:,:) = 0.0
      kl_SOX(:,:) = 0.0
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
          al(l) = MAX(-rd(l), al_sample(optimal_index,j))
        END IF
        gl(l) = MAX(0.0, gl_sample(optimal_index,j))
        ! gl = 0 means closed: al = -rd, as in set_closed.
        IF (gl(l) <= 0.0) al(l) = -rd(l)
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

IF (l_xylem_impairment .AND. .NOT. l_som_plant_segments) THEN
  ! Calculate the xylem conductance at the leaf water potential, for all
  ! points (open and closed), on the impaired vulnerability curve. (With
  ! segments leaf_k is the chosen sample's whole-plant conductance.)
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
        kl_hc_sample,                                                          &
! OUT (optional)
        psi_stem_sample                                                        &
)

USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls, xylem_conductance_jls
USE xylem_impairment_mod, ONLY: leaf_psi_impaired_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_segments_jls

USE jules_vegetation_mod, ONLY:                                                &
        photo_collatz, photo_farquhar, photo_johnson, photo_model,             &
        CW_conductance, SOX_conductance, som_psi_solver,                       &
        l_som_plant_segments

USE jb_photo_mod, ONLY: jb_eta_scale

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
                            ! resolutions/ranges in stom_opt_mod rather
                            ! than being fixed to the module-level
                            ! som_n_sample.

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
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL ::                               &
 psi_stem_sample(0:n_sample, open_pts)
                            ! Stem water potential (Pa) of each sample: the
                            ! stem segment outlet with segments and
                            ! impairment, else psi_root_zone.

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
,kmax_open(open_pts), kcrit_open(open_pts), b_open(open_pts), c_open(open_pts) &
,kcap_frac_pts(land_pts)                                                       &
,psi_stem_work(0:n_sample, open_pts)
                            ! Stem water potential of each sample.
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
!$OMP        ci_sample,ccp,km,je,n_sample,photo_model)
  DO j = 1,open_pts
    l = veg_index(open_index(j))
      wcarb_sample(:,j) = vcmax(l) * ( ci_sample(:,j) - ccp(l) )               &
                          / ( ci_sample(:,j) + km(l) )
      wlite_sample(:,j) = je(l) / 4.0 * ( ci_sample(:,j) - ccp(l) )            &
                          / ( ci_sample(:,j) + 2.0 * ccp(l) )
      ! Johnson-Berry: J depends on ci (see jb_photo_mod).
      IF ( photo_model == photo_johnson )                                      &
        wlite_sample(:,j) = wlite_sample(:,j)                                  &
                            * jb_eta_scale( ccp(l), ci_sample(:,j) )
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
IF (l_xylem_impairment .AND. l_som_plant_segments) THEN
  ! Segments: the memory model's cap on the stem and leaf segments (the
  ! root is not damaged); the hydraulic cost uses the intact segments'
  ! whole-plant conductance at the same water potentials.
  kcap_frac_pts(:) = 1.0
  DO j = 1, open_pts
    l = veg_index(open_index(j))
    kcap_frac_pts(l) = kmax(l) / MAX(kmax_ref(l), TINY(1.0_real_jlslsm))
  END DO
  CALL leaf_psi_segments_jls( pft, n_sample, land_pts, open_pts, veg_index,   &
                              open_index, el_sample(1:n_sample,:),             &
                              psi_root_zone, kmax_ref, kcrit,                  &
                              psi_sample(1:n_sample,:),                        &
                              kl_sample(1:n_sample,:),                         &
                              kcap_frac = kcap_frac_pts,                       &
                              leaf_k_intact = kl_hc_sample(1:n_sample,:),      &
                              stem_psi = psi_stem_work(1:n_sample,:) )
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

  ! Closed-stomata conductance (entry 0, below): the vulnerability curve at
  ! zero flow, K(psi_root_zone). With impairment, stom_opt_mod recomputes
  ! leaf_k for every point at the end, so only this path needs it.
  el_closed(:,:) = 0.0
  CALL leaf_psi_jls( pft, 1, land_pts, open_pts, veg_index, open_index,        &
                     el_closed, psi_root_zone, kmax, kcrit,                    &
                     conductance_b, conductance_c,                             &
                   ! INTENT OUT
                     psi_closed, kl_closed )

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
  ! K(psi_root_zone), not 0: 0 made PLC_pft report 100 % whenever the
  ! stomata shut (e.g. at twilight, when net A <= 0 at every ci). Entry 0 is
  ! never a profit candidate, so this only changes the reported leaf_k.
  IF (l_xylem_impairment) THEN
    kl_sample(0,j)  = 0.0
  ELSE
    kl_sample(0,j)  = kl_closed(1,j)
  END IF
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
        gl_max, l_multilayer, n_top, n_iter,                                   &
! OUT
        ci_g, al_g, gl_g, kl_g, psi_g, el_g, carbon_gain_g, hydraulic_cost_g   &
)

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE jules_vegetation_mod, ONLY: l_som_gain_gross, photo_collatz,               &
                                photo_farquhar, photo_johnson, photo_model,    &
                                CW_conductance,                                &
                                SOX_conductance,                               &
                                som_psi_solver, psi_solver_lut,           &
                                l_som_plant_segments, l_som_nsl
USE jb_photo_mod, ONLY: jb_eta_scale
USE pftparm, ONLY: c3, alpha, pft_conductance_model, conductance_b_pft,        &
                   conductance_c_pft, psi_nsl_onset, psi_nsl0
USE jules_surface_mod, ONLY: fwe_c3, fwe_c4, beta1, beta2, ratio
USE planet_constants_mod, ONLY: repsilon
USE c_rmol, ONLY: rmol
USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_hydraulics_CW_jls_mod, ONLY: supply_lut_psi, supply_lut_e_crit,       &
                                       supply_lut_f

INTEGER, INTENT(IN) ::                                                         &
  land_pts, pft, open_pts, open_index(land_pts), pft_photo_model,             &
  veg_index(land_pts), n_top, n_iter

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  rd(land_pts), ca(land_pts), psi_root_zone(land_pts), acr(land_pts),         &
  apar(land_pts), oi(land_pts), vcmax(land_pts), kc(land_pts), ko(land_pts),  &
  ccp(land_pts), pstar(land_pts), km(land_pts), dq(land_pts), je(land_pts),   &
  t_leaf(land_pts), je_ratio(land_pts), fapar_lf(land_pts), ipar(land_pts),   &
  kmax(land_pts), kcrit(land_pts), gl_max(land_pts)

LOGICAL, INTENT(IN) :: l_multilayer

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ci_g(open_pts), al_g(open_pts), gl_g(open_pts), kl_g(open_pts),            &
  psi_g(open_pts), el_g(open_pts), carbon_gain_g(open_pts),                  &
  hydraulic_cost_g(open_pts)

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

INTEGER :: i, j, l, side, idx1(land_pts)
REAL(KIND=real_jlslsm) :: b_curve(land_pts), c_curve(land_pts)
                            ! The PFT's intact vulnerability curve on land_pts
                            ! (l_som_fast is not coded for xylem impairment).
LOGICAL :: l_lut, ok_u, l_edge

REAL(KIND=real_jlslsm) ::                                                      &
  g_off, max_al, max_kl, ci_lo, ci_top, tol1,                                  &
  ! Evaluation state (eval_ci); al0_u is A without the nonstomatal limitation
  last_ci, al_u, gl_u, el_u, psi_u, kl_u, al0_u,                               &
  ! Edge (last feasible) and the infeasible end of its bracket
  e_ci, e_al, e_gl, e_el, e_psi, e_kl, e_g, b_ci, b_g, g_u, c_ci, t_al,        &
  ! Best seen
  best_f, best_ci, best_al, best_gl, best_el, best_psi, best_kl,               &
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
idx1(:) = 1
b_curve(:) = conductance_b_pft(pft)
c_curve(:) = conductance_c_pft(pft)

DO j = 1, open_pts
  l = veg_index(open_index(j))
  idx1(1) = open_index(j)
  g_off = MERGE(rd(l), 0.0_real_jlslsm, l_som_gain_gross)
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
  CALL eval_ci(ci_lo)
  max_kl = kl_u
  IF ( .NOT. ok_u ) THEN
    CALL set_closed()
    CYCLE
  END IF
  CALL store_best(-HUGE(1.0_real_jlslsm), ci_lo)
  e_ci = ci_lo; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
  e_kl = kl_u; e_g = margin()

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
    e_kl = kl_u
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

  ! Both ends of the feasible range are candidates.
  best_f = profit(best_al, best_kl)
  IF ( profit(e_al, e_kl) > best_f ) THEN
    best_f = profit(e_al, e_kl); best_ci = e_ci; best_al = e_al
    best_gl = e_gl; best_el = e_el; best_psi = e_psi; best_kl = e_kl
  END IF

  !---------------------------------------------------------------------------
  ! 3. Golden-section search for the maximum profit on [ccp, ci_b], stopping
  !    after n_iter iterations or once the bracket is below
  !    (ci_b - ccp)/opt_resol. (Brent's parabolic steps were tried here, but
  !    in low light profit is nearly flat over a wide Ci range and, at this
  !    precision, too noisy for them.)
  !---------------------------------------------------------------------------
  a = ci_lo
  b = e_ci
  tol1 = (b - a) / opt_resol
  IF ( tol1 > 0.0 ) THEN
    c_ci = b - golden_ratio * (b - a)
    d_ci = a + golden_ratio * (b - a)
    CALL eval_ci(c_ci)
    fc = profit_u()
    CALL eval_ci(d_ci)
    fd = profit_u()
    DO i = 1, n_iter
      IF ( b - a <= tol1 ) EXIT
      IF ( fc >= fd ) THEN
        b = d_ci
        d_ci = c_ci; fd = fc
        c_ci = b - golden_ratio * (b - a)
        CALL eval_ci(c_ci)
        fc = profit_u()
      ELSE
        a = c_ci
        c_ci = d_ci; fc = fd
        d_ci = a + golden_ratio * (b - a)
        CALL eval_ci(d_ci)
        fd = profit_u()
      END IF
    END DO
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
  hydraulic_cost_g(j) = (max_kl - best_kl) / (max_kl - kcrit(l))
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
  ! Nonstomatal limitation (l_som_nsl): gross photosynthesis scaled by
  ! f(psi_leaf), with psi_leaf from the transpiration that f * A gives. For
  ! fixed Ci, f(psi_leaf(f)) is non-increasing in f, so
  ! g(f) = f - f(psi_leaf(f)) is increasing with a single root on [0, 1]:
  ! g(1) > 0 here (the unlimited state, just evaluated), and at f = 0 there
  ! is no transpiration, psi_leaf = psi_root_zone, so g(0) = -f(psi_root)
  ! needs no hydraulics. Illinois (regula falsi) to |g| < nsl_tol, usually
  ! in 3-5 steps. Rd is not scaled.
  !-------------------------------------------------------------------------
  IF ( l_som_nsl ) THEN
    IF ( nsl_factor(psi_u) < 1.0 ) THEN
      ag0  = al0_u + rd(l)
      f_hi = 1.0
      g_hi = 1.0 - nsl_factor(psi_u)
      f_lo = 0.0
      g_lo = -nsl_factor(psi_root_zone(l))
      IF ( g_lo >= 0.0 ) THEN
        ! The soil alone stops photosynthesis (psi_root <= psi_nsl0).
        CALL hydraulic_state(ci, -rd(l))
      ELSE
        side = 0
        f_new = f_lo
        DO it = 1, n_nsl_iter
          f_new = f_lo - g_lo * (f_hi - f_lo) / (g_hi - g_lo)
          CALL hydraulic_state(ci, f_new * ag0 - rd(l))
          g_new = f_new - nsl_factor(psi_u)
          IF ( ABS(g_new) < nsl_tol ) EXIT
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

  ok_u = psi_u <= psi_root_zone(l) + TINY(psi_root_zone(l))                    &
         .AND. kl_u > kcrit(l)                                                 &
         .AND. (gl_max(l) <= 0.0 .OR. gl_u <= gl_max(l))
  END SUBROUTINE eval_ci

  ! Conductance, transpiration, leaf psi and xylem k at internal CO2 ci for
  ! net photosynthesis al (sets al_u, gl_u, el_u, psi_u, kl_u).
  SUBROUTINE hydraulic_state(ci, al)
  REAL(KIND=real_jlslsm), INTENT(IN) :: ci, al
  REAL(KIND=real_jlslsm) :: vpd
  REAL(KIND=real_jlslsm) :: el1(1,1), psi1(1,1), kl1(1,1)

  al_u = al
  gl_u = ratio * (al_u * rmol * t_leaf(l)) / MAX(ca(l) - ci, 1.0e-2_real_jlslsm)
  vpd = dq(l) * pstar(l) / repsilon
  el_u = MAX(0.0, vpd * gl_u / pstar(l) * pstar(l) / (rmol * t_leaf(l)))

  IF ( l_lut ) THEN
    psi_u = supply_lut_psi(pft, psi_root_zone(l),                              &
                           el_u / MAX(kmax(l), TINY(1.0_real_jlslsm)))
    kl_u = kmax(l) * supply_lut_f(pft, psi_u)
  ELSE
    el1(1,1) = el_u
    CALL leaf_psi_jls( pft, 1, land_pts, 1, veg_index, idx1, el1,              &
                       psi_root_zone, kmax, kcrit, b_curve, c_curve,           &
                       psi1, kl1 )
    psi_u = psi1(1,1)
    kl_u = kl1(1,1)
  END IF
  END SUBROUTINE hydraulic_state

  ! Nonstomatal limitation factor at leaf psi (Pa): Dewar et al. (2022)
  ! Eqn 3(b) generalised to an onset, f = 1 - (psi - psi_on)/(psi0 - psi_on)
  ! clipped to [0, 1] (see l_som_nsl).
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
  IF ( conv_e > 0.0 ) g_cap = MIN(g_cap, kmax(l)                              &
       * supply_lut_e_crit(pft, psi_root_zone(l),                              &
                           kcrit(l) / MAX(kmax(l), TINY(1.0_real_jlslsm))) / conv_e)
  IF ( g_cap >= HUGE(1.0_real_jlslsm) ) RETURN

  rk = ratio * rmol * t_leaf(l)
  ea = e_ci
  fa = g_cap * MAX(ca(l) - ea, 1.0e-2_real_jlslsm) - rk * e_al
  eb = ci_top
  fb = g_cap * MAX(ca(l) - eb, 1.0e-2_real_jlslsm) - rk * t_al
  IF ( fa < 0.0 .OR. fb >= 0.0 ) RETURN

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

  ! Full state at the edge; step back a little if interpolation in the
  ! table leaves it just outside.
  DO it = 1, 3
    IF ( ea <= e_ci ) RETURN
    CALL eval_ci(ea)
    IF ( ok_u ) THEN
      e_ci = ea; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
      e_kl = kl_u
      l_edge = .TRUE.
      RETURN
    END IF
    ea = ea - edge_rtol * (ca(l) - ea)
  END DO
  END SUBROUTINE edge_by_gl_cap

  !---------------------------------------------------------------------------
  ! Edge of the feasible range by an Illinois root find on the feasibility
  ! margin of full evaluations (any hydraulics).
  !---------------------------------------------------------------------------
  SUBROUTINE edge_by_margin()
  INTEGER :: it

  ! edge_by_gl_cap may have left a different evaluation in *_u.
  IF ( l_lut ) CALL eval_ci(ci_top)
  b_ci = ci_top
  b_g = margin()
  side = 0
  DO it = 1, max_edge_iter
    IF ( b_ci - e_ci <= edge_rtol * (ca(l) - e_ci) .OR. e_g < edge_margin )    &
      EXIT
    c_ci = (e_ci * b_g - b_ci * e_g) / (b_g - e_g)
    IF ( .NOT. ( c_ci > e_ci .AND. c_ci < b_ci ) ) c_ci = 0.5 * (e_ci + b_ci)
    CALL eval_ci(c_ci)
    g_u = margin()
    IF ( ok_u ) THEN
      e_ci = c_ci; e_al = al_u; e_gl = gl_u; e_el = el_u; e_psi = psi_u
      e_kl = kl_u; e_g = g_u
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

  ! Profit = CG - HC, as stom_opt_profit_max_select.
  REAL(KIND=real_jlslsm) FUNCTION profit(al_in, kl_in)
  REAL(KIND=real_jlslsm), INTENT(IN) :: al_in, kl_in
  profit = (al_in + g_off) / max_al - (max_kl - kl_in) / (max_kl - kcrit(l))
  END FUNCTION profit

  ! Profit of the latest evaluation (-HUGE if infeasible), kept if best.
  REAL(KIND=real_jlslsm) FUNCTION profit_u()
  IF ( ok_u ) THEN
    profit_u = profit(al_u, kl_u)
  ELSE
    profit_u = -HUGE(1.0_real_jlslsm)
  END IF
  IF ( profit_u > best_f ) CALL store_best(profit_u, last_ci)
  END FUNCTION profit_u

  SUBROUTINE store_best(f, ci)
  REAL(KIND=real_jlslsm), INTENT(IN) :: f, ci
  best_f = f; best_ci = ci; best_al = al_u; best_gl = gl_u; best_el = el_u
  best_psi = psi_u; best_kl = kl_u
  END SUBROUTINE store_best

  SUBROUTINE set_closed()
  ! Closed stomata: zero flow, so the conductance is K(psi_root_zone) on the
  ! same curve as hydraulic_state, not 0 (0 made PLC_pft report 100 %; see
  ! stom_opt_mod_ci).
  REAL(KIND=real_jlslsm) :: el0(1,1), psi0(1,1), kl0(1,1)
  IF ( l_lut ) THEN
    kl0(1,1) = kmax(l) * supply_lut_f(pft, psi_root_zone(l))
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

END MODULE stom_opt_jls_mod
