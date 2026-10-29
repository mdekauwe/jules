! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE triffid_nsc_mod
!-----------------------------------------------------------------------------
! TRIFFID storage pool (l_triffid_nsc), after Trugman et al. (2018, Ecol.
! Lett. 21: 1552-1560, Eqns 1a-c and 5, Table S2).
!
! NPP goes into a non-structural carbon store, C (kg C m-2 of the PFT), and
! TRIFFID grows from the carbon the store releases (phloem loading), not
! from the NPP of the step:
!   dC/dt = NPP - W,
!   W     = W_max C_a / (k_c C_max + C_a)          (Michaelis-Menten),
! with
!   C_max = nsc_flush (leaf + root carbon at the balanced LAI)
!           (Trugman et al.: enough for two flushes of leaves and fine roots),
!   C_a   = C - nsc_reserve_frac C_max, the accessible store (they take a
!           store below 10 % of its initial value as death by carbon
!           starvation),
!   k_c   = nsc_kc C_max (their k_c = 0.02 of the initial store),
!   W_max = nsc_wmax_fac * <NPP>, the running mean of NPP (e-folding time
!           nsc_tau_load years): the sink the plant is adapted to. Trugman
!           et al. set W_max by tree size and hold it through recovery; the
!           running mean is the JULES stand-in, as TRIFFID has no tree size,
!           and it is not held.
! The store is held between the reserve and C_max: carbon above C_max goes to
! growth (overflow), so a well-supplied plant grows from its NPP as without
! the store. When the accessible store is used up the carbon deficit passes
! to TRIFFID as negative growth (the vegetation shrinks), which is the JULES
! form of their starvation.
!
! Stem recovery (l_ximpair_rec_nsc, with the memory impairment model and
! dynamic TRIFFID): ring renewal. The stem's lost conductance recovers as
! the leaf's does, by new tissue at the rate the plant pays for it, and
! nothing is drawn from the store for it. The new ring is the wood TRIFFID
! grows in the step from the store's release W (after the leaves and fine
! roots), its gross wood production, net of the wood TRIFFID lost earlier
! in the ring year (the net ring count, below)
!   G = dC_W / dt + g_wood C_W     (vegcarb: the local wood increment plus
!                                   the replacement of wood turnover; the
!                                   carbon to spreading is not wood),
! so the carbon budget is that of no recovery. The ring's share of the
! stem conductance,
!   rho = nsc_xylem_cond_ratio G dt / X,
! with X the conducting sapwood (nsc_sapwood_frac woodc, or with
! nsc_la_sa > 0 the pipe model on the balanced leaf area as in FATES,
! X = LAI_bal / nsc_la_sa h nsc_wood_density cmass) and
! nsc_xylem_cond_ratio the new xylem's conductance per unit carbon
! relative to the mean sapwood, dilutes the whole stem's loss,
! L -> L (1 - rho), as leaf renewal dilutes the leaf's (random
! replacement of old and new xylem, intact and embolised alike, one pool;
! no cap, damage term or background floor). The store paces it through W:
! NPP, the turgor, cold and season gates, the flush and the accessible
! store (triffid_nsc_ring_renewal). With l_ximpair_leaf_loss the canopy
! is shed in proportion to the stem's functional xylem (Trugman et al.
! Eqn 2).
! (The store-paid repair, Trugman et al. Eqn 5, with its replacement
! rules 0 and 1, ring-porous earlywood rebuild, cambial cap, background
! loss and growth gate, nsc_xylem_replace, nsc_xylem_type,
! nsc_xylem_rate_max, nsc_plc_background, nsc_ring_mmdd and
! ximpair_psi_growth, was removed on 26/10/26; it is in the git history
! at 6133017.)
!
! Net ring count: while the release is shut (the season, cold and turgor
! limits) W = 0 but TRIFFID still charges the turnover, so the plant
! shrinks along its allometry, almost all of it in wood, and regrows when
! the release opens again: G < 0 in the shut steps and G > 0 in the
! regrowth. A deficit D (nsc_ring_deficit, in dumps) carries the
! shortfall: a step with G < 0 adds -G dt to D and renews nothing; a step
! with G > 0 repays D first, and only the rest is new ring. D is set to 0
! at the step that holds nsc_grow_end_doy (the ring year ends with the
! growing season), or, without a season end, fades with an e-folding time
! of a year, so a lasting loss of wood does not hold back the rings of
! later years. A ring year's new ring is then its net gross wood
! production, or 0. (Until 29/10/26 each step's G was clipped at 0, so the
! regrowth counted as new rings: at FR-Pue 209 against a net 142 g C m-2 a
! year, and renewal even of a stem whose wood was shrinking under severe
! drought. A build of 48fb366 reproduces it.)
!
! Sapwood state (l_nsc_sapwood_state, with ring renewal): the
! conducting stem sapwood X_s (nsc_sapwood, kg C m-2, in dumps) is a
! prognostic pool instead of the pipe model on the current leaf area,
!   dX_s/dt = G_s - k_h X_s,   G_s = nsc_stem_frac G,
! with k_h = nsc_heartwood_rate the share of the sapwood that becomes
! heartwood each year (about 1 / the number of rings in the sapwood) and
! nsc_stem_frac the stem's share of TRIFFID's gross wood production (the
! rest is branches and coarse roots). It starts from the pipe model at the
! first step of the main run, the ring's share is rho = w_cond G_s dt / X_s
! (X_s before the step), and in equilibrium rho = w_cond k_h: the
! recovery rate is set by a measurable sapwood turnover, not by the
! allometry, and the sapwood no longer follows the leaf area at once. As in
! MEDFATE (De Caceres et al. 2023), the sapwood is a pool with growth and
! heartwood turnover. No carbon moves: X_s is a share of TRIFFID's wood.
!
! Growth threshold (nsc_grow_thresh set): the release W for growth draws
! only on the store above nsc_grow_thresh C_max (MEDFATE grows sapwood only
! from starch above half its capacity), so after a deep drawdown the store
! refills before growth, and so ring renewal, resumes. The reserve still
! pays the carbon deficits and the leaf flush and regrowth still draw down
! to it (leaves before wood). The normalised release keeps W = W_max for a
! full store.
!
! Ring renewal replaces damaged xylem first (l_nsc_ring_additive): the
! ring's share is taken off the loss, L -> MAX(L - rho, 0), as in MEDFATE,
! instead of diluting it, L -> L (1 - rho). It assumes the sapwood retired
! to heartwood is the embolised part (old inner rings, vessels sealed by
! tyloses); the multiplicative default assumes it has the mean loss. The
! update is in update_xylem_impairment_memory.
!
! (The W_max hold, which froze the running mean of NPP while the stem
! took new drought damage, with nsc_hold_plc, nsc_repair_window and the
! memory ximpair_fmin_acc, was removed on 29/10/26: at FR-Pue it changed
! GPP, growth and damage by < 0.3 %. It is in the git history at 891f580.)
!
! Turgor-limited sink (nsc_psi_turg_zero set): growth needs turgor (Lockhart),
! so the release is limited by the plant's water status, not only by the
! store. With f_t the mean over the store period of
!   (psi_rz - nsc_psi_turg_zero) / (nsc_psi_turg_full - nsc_psi_turg_zero)
! (0 to 1; psi_rz the root-zone water potential, ~predawn; nsc_turg_acc),
!   W = W_max f_t [C_a / (k_c C_max + C_a)] / [C_a,max / (k_c C_max + C_a,max)]
! with C_a,max = (1 - nsc_reserve_frac) C_max: the full store releases W_max
! whatever k_c, so a softer release (larger nsc_kc) only shapes how loading
! falls as the store is drawn down. In drought the sink, not the supply,
! limits growth (Lempereur et al. 2015, Puechabon: stem growth stops at a
! predawn psi near -1.1 MPa while assimilation goes on), so the store fills
! in drought instead of paying for the intact sink. Carbon above C_max still
! goes to growth unless respired (nsc_overflow_resp, below). Growth
! respiration is not charged at release: the JULES NPP the store receives is
! already net of it.
!
! Cold limit (nsc_t_grow_zero set): cambial growth stops in the cold (wood
! formation needs daily mean temperatures above ~5-8 degC, Rossi et al.
! 2008; Lempereur et al. 2015 at Puechabon), so at each time step the turgor
! factor is multiplied by
!   (T_d - nsc_t_grow_zero) / (nsc_t_grow_full - nsc_t_grow_zero)   (0 to 1),
! with T_d the running mean of the air temperature over about a day (degC;
! e-folding 1 day, nsc_t_run, in dumps): a daily mean, so a sunny winter
! midday does not switch growth on. Either limit alone also turns on the
! normalised release.
!
! Growing season (nsc_grow_start_doy and / or nsc_grow_end_doy set): the
! growth factor is 0 before the start day and from the end day to the year
! end (days of the year). Trees can be dormant in winter even when it is
! mild, so the store refills from winter NPP; e.g. FR-Pue Q. ilex day 129
! to day 296 (Lempereur et al. 2015: stem growth starts in early-mid May,
! autumn growth ends in late October). Like the other limits it turns on the
! normalised release.
!
! Leaf flush (nsc_flush_start_doy set): leaves flush before stem growth
! starts and the flush draws on the store (Q. ilex at Puechabon: budburst
! about day 105, leaves open about day 127, Limousin et al. 2012; full size
! 6-7 weeks after budburst, Misson et al. 2011), so
! the store's release starts at budburst, not at nsc_grow_start_doy, and
! the year's new leaf carbon, g_leaf_0 x leaf carbon, is drawn from the
! accessible store evenly over nsc_flush_days days (part of W, like the
! regrowth after drought shedding). The draw only moves when carbon leaves
! the store: over a year the release still balances NPP.
!
! Leaf renewal (l_ximpair_rec_leaf_turn, with the per-segment memory): new
! leaves arrive intact, so the leaf segment's loss falls by the fraction of
! the canopy that is new, at each time step
!   dR = g_leaf_0 f_t / <f_t> dt,
! with f_t the growth factor above (1 without the limits) and <f_t> its
! running mean (e-folding 360 days, nsc_ft_mean, in dumps): leaves are made
! only when growth is allowed, and a year's production matches the turnover
! g_leaf_0. With phenology, regrowth after shedding renews the canopy too
! (phenol). This replaces the leaf segment's e-folding time ximpair_tau_leaf.
!
! Drought leaf shedding (l_phen_psi; Xu et al. 2016, New Phytol, Notes S1
! Eqn S7; phen_psi_update): the daily mean root-zone psi (predawn proxy)
! counts consecutive days below the turgor loss point phen_psi_tlp
! (phen_nlow) and above half of it (phen_nhigh); once either reaches
! phen_ndays the leaf biomass factor phen_fbleaf falls / rises by
! phen_dfbleaf a day (bounded by phen_fbleaf_min and 1), and phenol caps
! the canopy at it. Regrowth in phenol is paid from the store: its leaf
! carbon (dphen x leaf carbon at the balanced LAI), limited to the
! accessible store, is queued in nsc_regrow_acc and drawn at the next store
! step, before the store's own release, as part of W (TRIFFID holds leaf
! carbon at the balanced LAI, so the shed leaves were already charged as
! litter against W; over a shed-and-regrow cycle the store, not NPP, carries
! the cost, and the drawn carbon goes to TRIFFID growth).
!
! Release in the growth windows (l_nsc_wmax_norm, with the turgor or cold
! limit): W_max is divided by <f_t>, the running annual mean of the growth
! factor (nsc_ft_mean), so the store releases faster when growth is allowed
! (spring flush drawn from the store) and the annual release still matches
! the mean NPP.
!
! Capacity: nsc_cmax > 0 sets C_max directly (kg C m-2 of the PFT), e.g.
! from measured NSC concentrations times biomass, in place of nsc_flush
! (leaf + root carbon).
!
! Overflow respiration (nsc_overflow_resp > 0): of the carbon above C_max,
! that share is respired (nsc_resp_pft) instead of going to growth, so a
! sink-limited plant with a full store respires its surplus. The cap is
! applied after the store's input (NPP). nsc_resp_pft is added to the plant
! respiration and NPP outputs (extract_var, lagged by one store step), is
! part of resp_p_actual through nsc_dc, and leaves the carbon-conservation
! check of veg2.
!
! Starvation (nsc_g_starv > 0, with l_veg_compete): an area disturbance,
! nsc_starv_pft, added to g_area in lotka, at the full rate when the
! accessible store is used up and 0 for C_a >= nsc_kc C_max. When the PFT
! loses area its store goes with it, to litter (triffid_nsc_area); when it
! gains area the store is spread over the new area.
!
! Without TRIFFID (l_triffid off; prescribed LAI and height) the store runs
! once a day from sf_stom (triffid_nsc_daily): the day's NPP is summed
! (nsc_npp_day), and at the last time step of the day the same store step
! runs for that PFT with a one-day step. Leaf, root and wood carbon come
! from the prescribed height through the TRIFFID allometry
! (calc_c_comps_triffid at the balanced LAI). The release W leaves the store
! (there is no growth to feed), there is no stem recovery (it needs TRIFFID
! wood growth), and nsc_starv is a diagnostic only (the PFT fractions are
! fixed).
!
! MDK, 29/10/26
!
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: nsc_xylem_prod_pft, nsc_xylem_renew_pft   ! MDK, 29/10/26
PUBLIC :: nsc_ring_deficit   ! MDK, 29/10/26
PUBLIC :: nsc_sapwood   ! MDK, 29/10/26
PUBLIC :: nsc_pft, nsc_npp_mean, nsc_npp_wt, nsc_load_pft,                     &
          nsc_dc_pft, nsc_starv_pft, nsc_npp_day, nsc_turg_acc, nsc_resp_pft,  &
          nsc_t_run, nsc_ft_mean, phen_psi_acc, phen_nlow, phen_nhigh,         &
          phen_fbleaf, nsc_regrow_acc, phen_psi_update,                        &
          triffid_nsc_alloc, triffid_nsc_store, triffid_nsc_area,              &
          triffid_nsc_daily, triffid_nsc_turgor, triffid_nsc_ring_renewal

! Store (kg C m-2 of the PFT); < 0 = not yet set (filled to C_max at the
! first TRIFFID step). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_pft(:,:)
! Running mean of NPP (kg C m-2 per 360 days) and its weight, for W_max
! (bias-corrected mean = nsc_npp_mean / nsc_npp_wt). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_wt(:,:)
! Diagnostics of the last TRIFFID step (kg C m-2 per 360 days): the carbon
! released from the store (W, which TRIFFID grows from) and the net flux
! into the store (NPP - W).
!
! MDK, 26/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_load_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_dc_pft(:,:)
! Overflow respiration of the last store step (kg C m-2 per 360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_resp_pft(:,:)
! Ring renewal (stem recovery), last TRIFFID step: the new ring, TRIFFID's
! gross wood production G (kg C m-2 of the PFT per 360 days; the wood
! TRIFFID grows, not extra; net of nsc_ring_deficit), and the share of the
! stem conductance it renews, w_cond G / X (per 360 days). Diagnostics.
!
! MDK, 29/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_prod_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_renew_pft(:,:)
! Net ring count: the wood TRIFFID lost earlier in the ring year, which its
! regrowth repays before new rings count (kg C m-2 of the PFT). In dumps
! with ring renewal.
!
! MDK, 29/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_ring_deficit(:,:)
! Sapwood state (l_nsc_sapwood_state): the conducting stem sapwood X_s (kg C
! m-2 of the PFT; < 0 = not yet set: the pipe model at the first step). In
! dumps.
!
! MDK, 29/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_sapwood(:,:)
! Carbon-starvation mortality of the last TRIFFID step (per 360 days), used
! by lotka.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_starv_pft(:,:)
! Without TRIFFID: NPP summed over the day so far (kg C m-2 of the PFT). In
! dumps when l_triffid is off.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_day(:,:)
! Turgor-limited sink: the turgor factor summed over the store period so far,
! as a fraction of the period (< -1e29 = not yet accumulating: f_t = 1 at the
! next store step). In dumps when the sink is on.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_turg_acc(:,:)
! Cold limit: running mean of the air temperature, e-folding 1 day (K;
! < 0 = not yet set: the air temperature at the first step). In dumps when
! the cold limit is on.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_t_run(:,:)
! Leaf renewal: running mean of the growth factor, e-folding 360 days (< 0 =
! not yet set: 1). In dumps with l_ximpair_rec_leaf_turn.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_ft_mean(:,:)
! Drought leaf shedding (l_phen_psi), in dumps: root-zone psi summed over
! the day so far (Pa s), consecutive dry and wet days, the leaf biomass
! factor (0-1), and the regrowth carbon queued for the store (kg C m-2).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_psi_acc(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_nlow(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_nhigh(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_fbleaf(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_regrow_acc(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='TRIFFID_NSC_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Allocate the store and its diagnostics (n_land_pts, npft), if not done.
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(nsc_pft)) THEN
  ALLOCATE(nsc_pft(n_land_pts, npft))
  nsc_pft(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_npp_mean)) THEN
  ALLOCATE(nsc_npp_mean(n_land_pts, npft))
  nsc_npp_mean(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_npp_wt)) THEN
  ALLOCATE(nsc_npp_wt(n_land_pts, npft))
  nsc_npp_wt(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_load_pft)) THEN
  ALLOCATE(nsc_load_pft(n_land_pts, npft))
  nsc_load_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_dc_pft)) THEN
  ALLOCATE(nsc_dc_pft(n_land_pts, npft))
  nsc_dc_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_resp_pft)) THEN
  ALLOCATE(nsc_resp_pft(n_land_pts, npft))
  nsc_resp_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_xylem_prod_pft)) THEN   ! MDK, 25/10/26
  ALLOCATE(nsc_xylem_prod_pft(n_land_pts, npft))
  nsc_xylem_prod_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_xylem_renew_pft)) THEN   ! MDK, 25/10/26
  ALLOCATE(nsc_xylem_renew_pft(n_land_pts, npft))
  nsc_xylem_renew_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_ring_deficit)) THEN   ! MDK, 29/10/26
  ALLOCATE(nsc_ring_deficit(n_land_pts, npft))
  nsc_ring_deficit(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_sapwood)) THEN   ! MDK, 29/10/26
  ALLOCATE(nsc_sapwood(n_land_pts, npft))
  nsc_sapwood(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_starv_pft)) THEN
  ALLOCATE(nsc_starv_pft(n_land_pts, npft))
  nsc_starv_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_npp_day)) THEN
  ALLOCATE(nsc_npp_day(n_land_pts, npft))
  nsc_npp_day(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_turg_acc)) THEN
  ALLOCATE(nsc_turg_acc(n_land_pts, npft))
  nsc_turg_acc(:,:) = -1.0e30
END IF
IF (.NOT. ALLOCATED(nsc_t_run)) THEN
  ALLOCATE(nsc_t_run(n_land_pts, npft))
  nsc_t_run(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_ft_mean)) THEN
  ALLOCATE(nsc_ft_mean(n_land_pts, npft))
  nsc_ft_mean(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(phen_psi_acc)) THEN
  ALLOCATE(phen_psi_acc(n_land_pts, npft), phen_nlow(n_land_pts, npft),        &
           phen_nhigh(n_land_pts, npft), phen_fbleaf(n_land_pts, npft),        &
           nsc_regrow_acc(n_land_pts, npft))
  phen_psi_acc(:,:)   = 0.0
  phen_nlow(:,:)      = 0.0
  phen_nhigh(:,:)     = 0.0
  phen_fbleaf(:,:)    = 1.0
  nsc_regrow_acc(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_alloc

! ---------------------------------------------------------------------
! One TRIFFID step of the store. On return npp holds W, the carbon released
! from the store (kg C m-2 per 360 days), which TRIFFID grows from. The
! ring renewal (stem recovery) is triffid_nsc_ring_renewal, after vegcarb.
!
! MDK, 29/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_store( land_pts, trif_pts, trif_index, r_gamma,         &
                              leafc, rootc, npp, pft_only )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE pftparm, ONLY: is_set, nsc_in_flush_window, nsc_limits_on   ! MDK, 24/10/26
USE jules_vegetation_mod, ONLY: l_nsc_wmax_norm
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE pftparm, ONLY: g_leaf_0,                                                   &
                   nsc_flush, nsc_kc, nsc_reserve_frac, nsc_tau_load,          &
                   nsc_wmax_fac, nsc_g_starv,                                  &
                   nsc_psi_turg_zero, nsc_t_grow_zero, nsc_overflow_resp,      &
                   nsc_cmax, nsc_flush_start_doy, nsc_flush_days
USE pftparm, ONLY: nsc_grow_thresh   ! MDK, 29/10/26
USE model_time_mod, ONLY: current_time

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  trif_pts,                                                                    &
    ! Number of points on which TRIFFID may operate.
  trif_index(land_pts)
    ! Indices of land points on which TRIFFID may operate.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_gamma,                                                                     &
    ! Inverse TRIFFID timestep (/360 days).
  leafc(land_pts,npft),                                                        &
    ! Leaf carbon at the balanced LAI (kg C m-2).
  rootc(land_pts,npft)
    ! Root carbon (kg C m-2).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  npp(land_pts,npft)
    ! IN NPP, OUT the carbon released from the store, W
    ! (kg C m-2 per 360 days).

INTEGER, INTENT(IN), OPTIONAL :: pft_only
    ! Only this PFT (the daily step without TRIFFID); all PFTs if absent.

! Local variables
INTEGER :: n, t, l, n1, n2
REAL(KIND=real_jlslsm) ::                                                      &
  dt,                                                                          &
    ! TRIFFID timestep (360-day years).
  a_mean,                                                                      &
    ! Weight of this step in the running mean of NPP.
  c_max, c_res, c_acc,                                                         &
    ! Store capacity, reserve and accessible store (kg C m-2).
  w_max, w, draw,                                                              &
    ! Maximum and actual loading, net draw on the store (per 360 days).
  c_old,                                                                       &
    ! Store at the start of the step (kg C m-2).
  f_turg,                                                                      &
    ! Turgor-limited sink: turgor factor (0-1).
  c_thr, c_grow, c_grow_max,                                                   &
    ! The store the release for growth draws on: its floor (the reserve, or
    ! the growth threshold nsc_grow_thresh), the store above it, and that
    ! when full (kg C m-2).
    !
    ! MDK, 29/10/26
  w_regrow
    ! Phenology regrowth drawn from the store (per 360 days).
    !
    ! MDK, 26/10/26
LOGICAL :: l_turg, l_flush
    ! Growth limits on, leaf flush on.
    !
    ! MDK, 29/10/26
INTEGER :: doy_now, k, d, nd
    ! Day of year now; day counter, day of year and number of days of
    ! this store step.
    !
    ! MDK, 26/10/26
REAL(KIND=real_jlslsm) :: f_flush, w_flush
    ! Leaf flush: share of the year's flush in this step, flush draw (per
    ! 360 days).
    !
    ! MDK, 26/10/26
INTEGER :: n_days
    ! Days in this year.
    !
    ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_STORE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL triffid_nsc_alloc( land_pts )

dt = 1.0 / r_gamma
n1 = 1
n2 = nnpft
IF (PRESENT(pft_only)) THEN
  n1 = pft_only
  n2 = pft_only
END IF

DO n = n1,n2

  a_mean = MIN(dt / MAX(nsc_tau_load(n), TINY(1.0_real_jlslsm)), 1.0)
  l_turg = nsc_limits_on(n)
  ! Leaf flush: days of this store step (the nd days up to today) in the
  ! flush.
  !
  ! MDK, 26/10/26
  l_flush = is_set(nsc_flush_start_doy(n))
  f_flush = 0.0
  IF (l_flush) THEN
    n_days  = days_in_year(current_time%year, l_360, l_leap)
    doy_now = day_of_year(current_time%year, current_time%month,               &
                          current_time%day, l_360, l_leap)
    nd = MAX(NINT(dt * 360.0), 1)
    DO k = 0,nd-1
      d = MODULO(doy_now - k - 1, n_days) + 1
      ! The flush window may cross the year end.
      IF (nsc_in_flush_window(n, d, n_days))                                   &
        f_flush = f_flush + 1.0 / REAL(NINT(nsc_flush_days(n)))
    END DO
  END IF

  DO t = 1,trif_pts
    l = trif_index(t)

    c_max = nsc_flush(n) * (leafc(l,n) + rootc(l,n))
    IF (nsc_cmax(n) > 0.0) c_max = nsc_cmax(n)
    c_res = nsc_reserve_frac(n) * c_max
    IF (nsc_pft(l,n) < 0.0) nsc_pft(l,n) = c_max
    c_old = nsc_pft(l,n)

    ! Maximum loading from the running mean of NPP (bias-corrected).
    !
    ! MDK, 29/10/26
    nsc_npp_mean(l,n) = nsc_npp_mean(l,n) + a_mean * (npp(l,n)                 &
                        - nsc_npp_mean(l,n))
    nsc_npp_wt(l,n)   = nsc_npp_wt(l,n) + a_mean * (1.0 - nsc_npp_wt(l,n))
    w_max = 0.0
    IF (nsc_npp_wt(l,n) > 0.0)                                                 &
      w_max = nsc_wmax_fac(n) * MAX(nsc_npp_mean(l,n) / nsc_npp_wt(l,n), 0.0)
    ! Release concentrated in the growth windows.
    IF (l_nsc_wmax_norm .AND. l_turg .AND. nsc_ft_mean(l,n) > 0.0)             &
      w_max = w_max / MAX(nsc_ft_mean(l,n), 0.05)

    ! Leaf regrowth queued by phenol (l_phen_psi), paid first from the
    ! accessible store (phenol limited it to that), as part of the release.
    w_regrow = 0.0
    IF (nsc_regrow_acc(l,n) > 0.0) THEN
      w_regrow = MIN(nsc_regrow_acc(l,n), MAX(nsc_pft(l,n) - c_res, 0.0))
      nsc_pft(l,n) = nsc_pft(l,n) - w_regrow
      w_regrow = w_regrow / dt
      nsc_regrow_acc(l,n) = 0.0
    END IF

    ! Leaf flush: this step's share of the year's new leaves, from the
    ! accessible store, as part of the release.
    ! Phenology regrowth paid above in the flush window counts towards
    ! it, so a deciduous canopy rebuilt by phenol is not paid twice: the
    ! larger of the two is charged.
    !
    ! MDK, 20/10/26
    w_flush = 0.0
    IF (f_flush > 0.0) THEN
      w_flush = MIN(MAX(g_leaf_0(n) * leafc(l,n) * f_flush - w_regrow * dt,    &
                        0.0),                                                  &
                    MAX(nsc_pft(l,n) - c_res, 0.0))
      nsc_pft(l,n) = nsc_pft(l,n) - w_flush
      w_flush = w_flush / dt
    END IF

    ! Michaelis-Menten loading on the accessible store.
    c_acc = MAX(nsc_pft(l,n) - c_res, 0.0)
    ! The release for growth draws only on the store above the growth
    ! threshold, if set (the reserve still pays the deficits below).
    !
    ! MDK, 29/10/26
    c_thr = c_res
    IF (is_set(nsc_grow_thresh(n)))                                            &
      c_thr = MAX(c_res, nsc_grow_thresh(n) * c_max)
    c_grow = MAX(nsc_pft(l,n) - c_thr, 0.0)
    w = 0.0
    f_turg = 1.0
    IF (l_turg) THEN
      ! Turgor-limited sink, normalised so the full store releases W_max.
      IF (is_set(nsc_turg_acc(l,n)))                                           &
        f_turg = MIN(MAX(nsc_turg_acc(l,n), 0.0), 1.0)
      nsc_turg_acc(l,n) = 0.0
      c_grow_max = c_max - c_thr
      IF (c_grow > 0.0 .AND. c_grow_max > 0.0)                                 &
        w = w_max * f_turg * c_grow / (nsc_kc(n) * c_max + c_grow)             &
            * (nsc_kc(n) * c_max + c_grow_max) / c_grow_max
    ELSE
      IF (c_grow > 0.0) w = w_max * c_grow / (nsc_kc(n) * c_max + c_grow)
    END IF

    ! The store can not be drawn below the reserve: a deficit beyond the
    ! accessible store passes to TRIFFID (w = npp < 0 when it is empty).
    draw = w - npp(l,n)
    IF (draw > 0.0) draw = MIN(draw, c_acc / dt)
    w = npp(l,n) + draw
    nsc_pft(l,n) = nsc_pft(l,n) - draw * dt

    ! The regrowth paid above is part of the release.
    w = w + w_regrow + w_flush

    ! Overflow above the capacity goes to growth, less the share respired.
    nsc_resp_pft(l,n) = 0.0
    IF (nsc_pft(l,n) > c_max) THEN
      nsc_resp_pft(l,n) = nsc_overflow_resp(n) * (nsc_pft(l,n) - c_max) / dt
      w = w + (1.0 - nsc_overflow_resp(n)) * (nsc_pft(l,n) - c_max) / dt
      nsc_pft(l,n) = c_max
    END IF

    ! Carbon starvation, ramping up as the accessible store runs out.
    nsc_starv_pft(l,n) = 0.0
    IF (nsc_g_starv(n) > 0.0 .AND. c_max > 0.0) THEN
      c_acc = MAX(nsc_pft(l,n) - c_res, 0.0)
      nsc_starv_pft(l,n) = nsc_g_starv(n)                                      &
                           * MAX(1.0 - c_acc / (nsc_kc(n) * c_max), 0.0)
    END IF

    npp(l,n)            = w
    nsc_load_pft(l,n)   = w
    nsc_dc_pft(l,n)     = (nsc_pft(l,n) - c_old) * r_gamma

  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_store

! ---------------------------------------------------------------------
! The store after TRIFFID's change of PFT area (lotka). Area lost takes its
! store with it, to litter (lit_store); area gained shares the store, which
! is spread over the new area.
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_area( land_pts, trif_pts, trif_index, r_gamma,          &
                             frac, dfrac, lit_store )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  trif_pts,                                                                    &
    ! Number of points on which TRIFFID may operate.
  trif_index(land_pts)
    ! Indices of land points on which TRIFFID may operate.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_gamma,                                                                     &
    ! Inverse TRIFFID timestep (/360 days).
  frac(land_pts,ntype),                                                        &
    ! Areal fractions after the TRIFFID step.
  dfrac(land_pts,npft)
    ! Change in the areal fractions over the step.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  lit_store(land_pts,npft)
    ! Store carbon of the area lost, to litter (kg C m-2 of land per 360
    ! days).

INTEGER :: n, t, l
REAL(KIND=real_jlslsm) :: frac_old
    ! Areal fraction before the TRIFFID step.
    !
    ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_AREA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lit_store(:,:) = 0.0
IF (.NOT. ALLOCATED(nsc_pft)) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO n = 1,nnpft
  DO t = 1,trif_pts
    l = trif_index(t)
    IF (nsc_pft(l,n) <= 0.0) CYCLE
    frac_old = frac(l,n) - dfrac(l,n)
    IF (dfrac(l,n) < 0.0) THEN
      lit_store(l,n) = - dfrac(l,n) * nsc_pft(l,n) * r_gamma
    ELSE IF (dfrac(l,n) > 0.0 .AND. frac(l,n) > 0.0) THEN
      nsc_pft(l,n) = nsc_pft(l,n) * frac_old / frac(l,n)
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_area

! ---------------------------------------------------------------------
! Turgor-limited sink: called from sf_stom every time step for PFT ft, adds
! this time step's turgor factor, from the root-zone water potential, to
! nsc_turg_acc as a fraction of the store period (the TRIFFID period, or a
! day without TRIFFID).
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_turgor( land_pts, veg_pts, veg_index, ft, timestep,     &
                               psi_rz, tair )

USE jules_vegetation_mod, ONLY: l_triffid, triffid_period,                     &
                                l_ximpair_rec_leaf_turn, l_nsc_wmax_norm
USE pftparm, ONLY: is_set, nsc_limits_on, nsc_season_on   ! MDK, 24/10/26
USE pftparm, ONLY: nsc_psi_turg_full, nsc_psi_turg_zero, nsc_t_grow_zero,      &
                   nsc_t_grow_full, g_leaf_0, nsc_grow_start_doy,              &
                   nsc_grow_end_doy, nsc_flush_start_doy
USE model_time_mod, ONLY: current_time, is_spinup   ! is_spinup: MDK, 25/10/26
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                   &
                                          ximpair_seg_state_alloc
USE water_constants_mod, ONLY: tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  veg_pts,                                                                     &
    ! Number of points of this PFT.
  veg_index(land_pts),                                                         &
    ! Land points of this PFT.
  ft
    ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model time step (s).
  psi_rz(land_pts),                                                            &
    ! Root-zone water potential (Pa).
  tair(land_pts)
    ! Air temperature (K).

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: period, f, a_t, a_yr, d_renew, f_season, doy,        &
                          doy_start
    ! Store period (s); growth factor (0-1); weights of this time step in
    ! the daily and annual running means; leaf renewal this time step;
    ! growing-season factor (0 or 1); day of year now and of the start of
    ! the release.
    !
    ! MDK, 24/10/26
LOGICAL :: l_lim, l_season
    ! Growth limits on; growing season set.
    !
    ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_TURGOR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

l_season = nsc_season_on(ft)
l_lim = nsc_limits_on(ft)
IF (.NOT. l_lim .AND. .NOT. l_ximpair_rec_leaf_turn) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL triffid_nsc_alloc( land_pts )
IF (l_ximpair_rec_leaf_turn) CALL ximpair_seg_state_alloc( land_pts )

period = 86400.0
IF (l_triffid) period = 86400.0 * REAL(triffid_period)
a_t = 1.0 - EXP(-timestep / 86400.0)
a_yr = 1.0 - EXP(-timestep / (360.0 * 86400.0))
! Growing season: 0 outside the start and end days of the year.
!
! MDK, 20/10/26
f_season = 1.0
IF (l_season) THEN
  doy = REAL(day_of_year(current_time%year, current_time%month,                &
                         current_time%day, l_360, l_leap))
  ! With the leaf flush the release starts at budburst.
  !
  ! MDK, 20/10/26
  doy_start = nsc_grow_start_doy(ft)
  IF (is_set(nsc_flush_start_doy(ft))) doy_start = nsc_flush_start_doy(ft)
  IF (is_set(doy_start) .AND. doy < doy_start) f_season = 0.0
  IF (is_set(nsc_grow_end_doy(ft)) .AND. doy >= nsc_grow_end_doy(ft))          &
    f_season = 0.0
END IF

DO j = 1,veg_pts
  l = veg_index(j)
  f = 1.0
  IF (is_set(nsc_psi_turg_zero(ft)))                                           &
    f = MIN(MAX((1.0e-6 * psi_rz(l) - nsc_psi_turg_zero(ft))                   &
                / (nsc_psi_turg_full(ft) - nsc_psi_turg_zero(ft)), 0.0), 1.0)
  IF (is_set(nsc_t_grow_zero(ft))) THEN
    ! Daily running mean of the air temperature.
    IF (nsc_t_run(l,ft) < 0.0) nsc_t_run(l,ft) = tair(l)
    nsc_t_run(l,ft) = nsc_t_run(l,ft) + a_t * (tair(l) - nsc_t_run(l,ft))
    f = f * MIN(MAX((nsc_t_run(l,ft) - tm - nsc_t_grow_zero(ft))               &
                    / (nsc_t_grow_full(ft) - nsc_t_grow_zero(ft)), 0.0), 1.0)
  END IF
  f = f * f_season
  IF (l_lim) THEN
    IF (.NOT. is_set(nsc_turg_acc(l,ft))) nsc_turg_acc(l,ft) = 0.0
    nsc_turg_acc(l,ft) = nsc_turg_acc(l,ft) + f * timestep / period
  END IF
  ! Running annual mean of the growth factor (leaf renewal, l_nsc_wmax_norm).
  IF (l_ximpair_rec_leaf_turn .OR. l_nsc_wmax_norm) THEN
    IF (nsc_ft_mean(l,ft) < 0.0) nsc_ft_mean(l,ft) = 1.0
    nsc_ft_mean(l,ft) = nsc_ft_mean(l,ft) + a_yr * (f - nsc_ft_mean(l,ft))
  END IF
  ! Leaf renewal: production at the turnover rate, made while growing. Not
  ! in spin-up, where the stem damage is frozen (the update returns), so the
  ! leaf cap does not heal on its own either. MDK, 25/10/26
  IF (l_ximpair_rec_leaf_turn .AND. .NOT. is_spinup) THEN
    d_renew = MIN(g_leaf_0(ft) * f / MAX(nsc_ft_mean(l,ft), 0.05)              &
                  * timestep / (360.0 * 86400.0), 1.0)
    ximpair_kcap_leaf(l,ft) = 1.0 - (1.0 - ximpair_kcap_leaf(l,ft))            &
                              * (1.0 - d_renew)
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_turgor

! ---------------------------------------------------------------------
! Drought leaf shedding (l_phen_psi): called from sf_stom every time step for
! PFT ft. Sums the root-zone psi over the day; at the last time step of the
! day updates the consecutive dry / wet day counts and the leaf biomass
! factor phen_fbleaf (Xu et al. 2016, Eqn S7).
! ---------------------------------------------------------------------
SUBROUTINE phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,        &
                            psi_rz )

USE model_time_mod, ONLY: current_time
USE pftparm, ONLY: is_set   ! MDK, 24/10/26
USE pftparm, ONLY: phen_psi_tlp, phen_dfbleaf, phen_ndays, phen_fbleaf_min

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  veg_pts,                                                                     &
    ! Number of points of this PFT.
  veg_index(land_pts),                                                         &
    ! Land points of this PFT.
  ft
    ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model time step (s).
  psi_rz(land_pts)
    ! Root-zone water potential (Pa).

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: psi_d
    ! Daily mean root-zone water potential (MPa).
    !
    ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PHEN_PSI_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. is_set(phen_psi_tlp(ft))) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL triffid_nsc_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)
  phen_psi_acc(l,ft) = phen_psi_acc(l,ft) + psi_rz(l) * timestep
END DO

! The day's update at its last time step.
IF (REAL(current_time%time) + timestep < 86400.0 - 0.5) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO j = 1,veg_pts
  l = veg_index(j)
  ! Daily mean root-zone psi (MPa), the predawn proxy (a part day after a
  ! mid-day start counts as wetter: runs here start at the end of a day).
  psi_d = 1.0e-6 * phen_psi_acc(l,ft) / 86400.0
  phen_psi_acc(l,ft) = 0.0
  IF (psi_d < phen_psi_tlp(ft)) THEN
    phen_nlow(l,ft) = phen_nlow(l,ft) + 1.0
  ELSE
    phen_nlow(l,ft) = 0.0
  END IF
  IF (psi_d > 0.5 * phen_psi_tlp(ft)) THEN
    phen_nhigh(l,ft) = phen_nhigh(l,ft) + 1.0
  ELSE
    phen_nhigh(l,ft) = 0.0
  END IF
  IF (phen_nlow(l,ft) >= phen_ndays(ft)) THEN
    phen_fbleaf(l,ft) = phen_fbleaf(l,ft) - phen_dfbleaf(ft)
  ELSE IF (phen_nhigh(l,ft) >= phen_ndays(ft)) THEN
    phen_fbleaf(l,ft) = phen_fbleaf(l,ft) + phen_dfbleaf(ft)
  END IF
  phen_fbleaf(l,ft) = MIN(MAX(phen_fbleaf(l,ft), phen_fbleaf_min(ft)), 1.0)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE phen_psi_update

! ---------------------------------------------------------------------
! The store without TRIFFID (l_triffid off), called from sf_stom every time
! step for PFT ft after NPP: sums the day's NPP, and at the last time step of
! the day runs the store step for ft with a one-day step, with leaf and
! root carbon from the prescribed height (TRIFFID allometry).
!
! MDK, 26/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_daily( land_pts, veg_pts, veg_index, ft, timestep,      &
                              ht, npp )

USE jules_surface_types_mod, ONLY: npft
USE model_time_mod, ONLY: current_time
USE calc_c_comps_triffid_mod, ONLY: calc_c_comps_triffid

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  veg_pts,                                                                     &
    ! Number of points of this PFT.
  veg_index(land_pts),                                                         &
    ! Land points of this PFT.
  ft
    ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model time step (s).
  ht(land_pts),                                                                &
    ! Canopy height (m).
  npp(land_pts)
    ! NPP of this time step (kg C m-2 s-1).

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: c_veg, lai_b, woodc
    ! Vegetation carbon (kg C m-2), balanced LAI and wood carbon (kg C m-2)
    ! from the allometry.
    !
    ! MDK, 26/10/26
REAL(KIND=real_jlslsm), ALLOCATABLE :: leafc(:,:), rootc(:,:), npp_d(:,:)
    ! Inputs and outputs of triffid_nsc_store for this PFT (see there).
    !
    ! MDK, 26/10/26

REAL(KIND=real_jlslsm), PARAMETER :: sec_per_day = 86400.0
REAL(KIND=real_jlslsm), PARAMETER :: days_per_trif_year = 360.0

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_DAILY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL triffid_nsc_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)
  nsc_npp_day(l,ft) = nsc_npp_day(l,ft) + npp(l) * timestep
END DO

! Store step at the last time step of the day.
IF (REAL(current_time%time) + timestep < sec_per_day - 0.5) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

ALLOCATE(leafc(land_pts,npft), rootc(land_pts,npft), npp_d(land_pts,npft))
leafc(:,:) = 0.0
rootc(:,:) = 0.0
npp_d(:,:) = 0.0
DO j = 1,veg_pts
  l = veg_index(j)
  CALL calc_c_comps_triffid(ft, ht(l), lai_b, leafc(l,ft), rootc(l,ft),        &
                            woodc, c_veg)
  ! The day's NPP (kg C m-2) as a rate per 360 days (TRIFFID units).
  npp_d(l,ft) = nsc_npp_day(l,ft) * days_per_trif_year
  nsc_npp_day(l,ft) = 0.0
END DO

CALL triffid_nsc_store( land_pts, veg_pts, veg_index, days_per_trif_year,      &
                        leafc, rootc, npp_d, pft_only = ft )

DEALLOCATE(leafc, rootc, npp_d)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_daily

! ---------------------------------------------------------------------
! Stem recovery by ring renewal (l_ximpair_rec_nsc): called by triffid
! after vegcarb. The new ring is TRIFFID's gross wood production of
! the step, G = dwood / dt + g_wood C_W,old (wood_litc, the turnover
! replaced), per 360 days, never negative. Its share of the stem
! conductance, nsc_xylem_cond_ratio G dt / X (X the conducting sapwood
! before the step's growth, as the store used), is queued in
! ximpair_repair_pend; the next impairment update dilutes the whole stem's
! loss by it. Not in spin-up (the impairment is frozen there). No carbon
! moves here. G is net of the deficit nsc_ring_deficit: wood TRIFFID lost
! earlier in the ring year (G < 0 while the release was shut) is regrown
! before any new ring counts (see the module header).
!
! MDK, 29/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_ring_renewal( land_pts, trif_pts, trif_index,           &
                                     r_gamma, ht, lai_bal, woodc, dwood,       &
                                     wood_litc )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE jules_vegetation_mod, ONLY: l_ximpair_rec_nsc, xylem_impairment_memory
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_state   ! MDK, 29/10/26
USE pftparm, ONLY: pft_xylem_impairment_model, nsc_xylem_cond_ratio
USE pftparm, ONLY: is_set, nsc_grow_end_doy, nsc_heartwood_rate,               &
                   nsc_stem_frac   ! MDK, 29/10/26
USE calc_c_comps_triffid_mod, ONLY: sapwood_c
USE model_time_mod, ONLY: is_spinup
USE model_time_mod, ONLY: current_time   ! MDK, 29/10/26
USE datetime_utils_mod, ONLY: day_of_year, days_in_year   ! MDK, 29/10/26
USE time_info_mod, ONLY: l_360, l_leap   ! MDK, 29/10/26
USE xylem_impairment_memory_mod, ONLY: ximpair_repair_pend,                    &
                                       ximpair_memory_alloc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Total number of land points.
  trif_pts,                                                                    &
    ! Number of points on which TRIFFID may operate.
  trif_index(land_pts)
    ! Indices of land points on which TRIFFID may operate.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  r_gamma,                                                                     &
    ! Inverse TRIFFID timestep (/360 days).
  ht(land_pts,npft),                                                           &
    ! Canopy height before the step's growth (m).
  lai_bal(land_pts,npft),                                                      &
    ! Balanced LAI before the step's growth.
  woodc(land_pts,npft),                                                        &
    ! Wood carbon after the step's growth (kg C m-2).
  dwood(land_pts,npft),                                                        &
    ! Wood increment of the step (kg C m-2).
  wood_litc(land_pts,npft)
    ! Wood turnover of the step, g_wood C_W (kg C m-2 per 360 days).

INTEGER :: n, t, l
REAL(KIND=real_jlslsm) :: g_wood_prod, x_s
    ! Gross wood production (kg C m-2 per 360 days); conducting sapwood
    ! (kg C m-2).
    !
    ! MDK, 25/10/26
REAL(KIND=real_jlslsm) :: dt, g_step, f_decay
    ! Net ring count: TRIFFID timestep (360-day years), gross wood
    ! production of the step (kg C m-2, < 0 when wood is lost) and the
    ! deficit's decay over the step without a season end.
    !
    ! MDK, 29/10/26
REAL(KIND=real_jlslsm) :: g_ring
    ! The new conducting xylem of the step: G, or with the sapwood state
    ! the stem's share of it, G_s (kg C m-2 per 360 days).
    !
    ! MDK, 29/10/26
INTEGER :: doy_now, n_days, nd, k, d, doy_end
    ! Net ring count: day of year now, days in this year, days of this
    ! step, day counter, a day of the step and the season's end day.
    !
    ! MDK, 29/10/26
LOGICAL :: l_ring_end
    ! Net ring count: this step holds the end of the growing season.
    !
    ! MDK, 29/10/26

REAL(KIND=real_jlslsm), PARAMETER :: x_min = 1.0e-6
    ! Minimum conducting xylem (kg C m-2).
    !
    ! MDK, 25/10/26
REAL(KIND=real_jlslsm), PARAMETER :: tau_deficit = 1.0
    ! Net ring count without a season end: e-folding time of the deficit
    ! (360-day years), one ring year.
    !
    ! MDK, 29/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_RING_RENEWAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_ximpair_rec_nsc .AND. .NOT. is_spinup) THEN
  CALL triffid_nsc_alloc( land_pts )
  CALL ximpair_memory_alloc( land_pts )
  ! Net ring count: the days of this step (the nd days up to today, as the
  ! leaf flush in triffid_nsc_store counts them).
  !
  ! MDK, 29/10/26
  dt      = 1.0 / r_gamma
  n_days  = days_in_year(current_time%year, l_360, l_leap)
  doy_now = day_of_year(current_time%year, current_time%month,                 &
                        current_time%day, l_360, l_leap)
  nd      = MAX(NINT(dt * 360.0), 1)
  f_decay = EXP(-dt / tau_deficit)
  DO n = 1,nnpft
    IF (pft_xylem_impairment_model(n) /= xylem_impairment_memory) CYCLE
    ! Net ring count: the ring year ends at the step that holds the end of
    ! the growing season.
    !
    ! MDK, 29/10/26
    l_ring_end = .FALSE.
    IF (is_set(nsc_grow_end_doy(n))) THEN
      doy_end = MIN(NINT(nsc_grow_end_doy(n)), n_days)
      DO k = 0,nd-1
        d = MODULO(doy_now - k - 1, n_days) + 1
        IF (d == doy_end) l_ring_end = .TRUE.
      END DO
    END IF
    DO t = 1,trif_pts
      l = trif_index(t)
      nsc_xylem_prod_pft(l,n)  = 0.0
      nsc_xylem_renew_pft(l,n) = 0.0
      x_s = sapwood_c(n, ht(l,n), lai_bal(l,n), woodc(l,n) - dwood(l,n))
      ! Sapwood state: the pipe model sets it at the first step of the
      ! main run; then the state is the sapwood the ring renews. (It can
      ! not fall to x_min: that would take centuries without growth.)
      !
      ! MDK, 29/10/26
      IF (l_nsc_sapwood_state) THEN
        IF (nsc_sapwood(l,n) < 0.0) nsc_sapwood(l,n) = x_s
        x_s = nsc_sapwood(l,n)
      END IF
      IF (x_s <= x_min) CYCLE
      ! A new ring year starts with no deficit; without a season end the
      ! deficit fades over a year.
      !
      ! MDK, 29/10/26
      IF (is_set(nsc_grow_end_doy(n))) THEN
        IF (l_ring_end) nsc_ring_deficit(l,n) = 0.0
      ELSE
        nsc_ring_deficit(l,n) = nsc_ring_deficit(l,n) * f_decay
      END IF
      ! Wood lost in the step is owed to the regrowth (no new ring); the
      ! regrowth repays the deficit before the rest counts as new ring.
      !
      ! MDK, 29/10/26
      g_step = dwood(l,n) + wood_litc(l,n) * dt
      IF (g_step < 0.0) THEN
        nsc_ring_deficit(l,n) = nsc_ring_deficit(l,n) - g_step
        g_wood_prod = 0.0
      ELSE
        g_wood_prod = MAX(g_step - nsc_ring_deficit(l,n), 0.0) * r_gamma
        nsc_ring_deficit(l,n) = MAX(nsc_ring_deficit(l,n) - g_step, 0.0)
      END IF
      nsc_xylem_prod_pft(l,n)  = g_wood_prod
      ! With the sapwood state the stem's share of the new wood renews the
      ! stem sapwood, which loses its oldest part to heartwood.
      !
      ! MDK, 29/10/26
      g_ring = g_wood_prod
      IF (l_nsc_sapwood_state) g_ring = nsc_stem_frac(n) * g_wood_prod
      nsc_xylem_renew_pft(l,n) = nsc_xylem_cond_ratio(n) * g_ring / x_s
      ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n)                      &
                                 + nsc_xylem_renew_pft(l,n) / r_gamma
      IF (l_nsc_sapwood_state)                                                 &
        nsc_sapwood(l,n) = nsc_sapwood(l,n)                                    &
                           * EXP(-nsc_heartwood_rate(n) / r_gamma)             &
                           + g_ring / r_gamma
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_ring_renewal


END MODULE triffid_nsc_mod
