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
! NPP goes into a non-structural carbon store C (kg C m-2 of the PFT), and
! TRIFFID grows from the carbon the store releases (phloem loading), W, not
! from the NPP of the step:
!
!   dC/dt = NPP - W,
!   W     = W_max C_a / (k_c C_max + C_a)          (Michaelis-Menten),
!
! with
!   C_max = nsc_flush (leaf + root carbon at the balanced LAI), enough for
!           two flushes of leaves and fine roots, or nsc_cmax if set,
!   C_a   = C - nsc_reserve_frac C_max, the accessible store,
!   k_c   = nsc_kc C_max,
!   W_max = nsc_wmax_fac <NPP>, the running mean of NPP (e-folding
!           nsc_tau_load years): the sink the plant is adapted to (Trugman
!           et al. set it by tree size, which TRIFFID does not have).
!
! Carbon above C_max goes to growth, less the share nsc_overflow_resp
! respired. When the accessible store is used up the deficit passes to
! TRIFFID as negative growth (the vegetation shrinks), the JULES form of
! carbon starvation; with nsc_g_starv > 0 (and l_veg_compete) it also
! drives an area loss in lotka.
!
! Growth limits (triffid_nsc_turgor). The release is scaled by a growth
! factor f (0-1), the product of a turgor limit (nsc_psi_turg_zero set:
! growth needs turgor, so in drought the sink, not the supply, limits
! growth),
!
!   f_psi = (psi_rz - psi_0) / (psi_1 - psi_0),
!
! with psi_rz the root-zone (~predawn) water potential, psi_0 =
! nsc_psi_turg_zero and psi_1 = nsc_psi_turg_full; a cold limit
! (nsc_t_grow_zero set),
!
!   f_T = (T_d - T_0) / (T_1 - T_0),
!
! with T_d the daily mean air temperature, T_0 = nsc_t_grow_zero and T_1 =
! nsc_t_grow_full; and the growing season (0 outside nsc_grow_start_doy to
! nsc_grow_end_doy). With any limit on, the release is normalised so a
! full store still releases W_max,
!
!   W = W_max f [C_a / (k_c C_max + C_a)] / [C_m / (k_c C_max + C_m)],
!
! with f the mean over the store period and C_m the accessible store when
! full. With l_nsc_wmax_norm, W_max is also divided by the running annual
! mean of f, so the store releases faster while growth is allowed.
!
! Leaves. The leaf flush (nsc_flush_start_doy set) draws the year's new
! leaf carbon from the store over nsc_flush_days, before stem growth
! starts. Drought leaf shedding (l_phen_psi, Xu et al. 2016;
! phen_psi_update) sheds and regrows the canopy on the root-zone psi, and
! the store pays for the regrowth. Leaf renewal (l_ximpair_rec_leaf_turn)
! renews the leaf segment's conductance by the canopy's new fraction.
!
! Stem recovery by ring renewal (l_ximpair_rec_nsc, with the memory
! impairment model and dynamic TRIFFID; triffid_nsc_ring_renewal). The new
! ring is TRIFFID's gross wood production from the release,
!
!   G = dC_W/dt + g_wood C_W,
!
! and its share of the stem conductance comes off the stem's loss L,
!
!   rho = nsc_xylem_cond_ratio G dt / X,    L -> MAX(L - rho, 0),
!
! with X the conducting sapwood (nsc_sapwood_frac of the wood, or the pipe
! model on the balanced leaf area with nsc_la_sa > 0, as in FATES). No
! carbon is drawn from the store for it: the store paces the renewal
! through W. The ring replaces damaged xylem first (l_nsc_ring_additive,
! the default, as in MEDFATE); with it off the ring dilutes the loss,
! L -> L (1 - rho).
!
! Net wood growth. While the release is shut TRIFFID still charges the
! turnover, so the wood shrinks (G < 0) and regrows when the release opens.
! A deficit D carries the wood lost, which the regrowth repays before any
! new ring counts. At the end of the ring year (the step that holds
! nsc_grow_end_doy, or the last day of the year) its net production S
! settles the count,
!
!   D <- MAX(D + MIN(S, 0), 0),    S <- 0,
!
! so the count over the years is the sum of MAX(S, 0).
!
! Sapwood pool (l_nsc_sapwood_pool, on by default with ring renewal). The
! conducting sapwood X_s of TRIFFID's wood is a prognostic pool, as in
! MEDFATE (De Caceres et al. 2023), instead of the pipe model:
!
!   dX_s/dt = G - (k_h + g_t) X_s,
!
! with k_h = nsc_heartwood_rate the share of the sapwood that becomes
! heartwood a year and g_t the wood turnover TRIFFID applies (its wood
! litter over its wood, g_wood), which takes the same share of the
! sapwood as of the wood. The pool starts, where the dump has none, at
! its steady state for the current wood C_W (G = g_t C_W in balance),
!
!   X_s = C_W g_t / (k_h + g_t),
!
! so it does not drift through the run. The ring renews about
! k_h + g_t of the sapwood a year, so a measurable sapwood turnover, not
! the allometry, sets the recovery rate. A ring year with a net loss of
! wood also takes |S| from X_s.
!
! Growth threshold (nsc_grow_thresh set): growth draws only on the store
! above nsc_grow_thresh C_max, so after a deep drawdown the store refills
! before growth and ring renewal resume; the reserve still pays the
! deficits and the leaves.
!
! Without TRIFFID (l_triffid off) the store runs once a day
! (triffid_nsc_daily), with leaf, root and wood carbon from the prescribed
! height; the release leaves the store and there is no stem recovery.
!
! MDK, 31/10/26
!
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: nsc_xylem_prod_pft, nsc_xylem_renew_pft
PUBLIC :: nsc_wood_lost
PUBLIC :: nsc_net_wood_growth
PUBLIC :: nsc_sapwood
PUBLIC :: nsc_pft, nsc_npp_mean, nsc_npp_wt, nsc_load_pft,                     &
          nsc_dc_pft, nsc_starv_pft, nsc_npp_day, nsc_turg_acc, nsc_resp_pft,  &
          nsc_t_run, nsc_ft_mean, phen_psi_acc, phen_nlow, phen_nhigh,         &
          phen_fbleaf, nsc_regrow_acc, phen_psi_update,                        &
          triffid_nsc_alloc, triffid_nsc_store, triffid_nsc_area,              &
          triffid_nsc_daily, triffid_nsc_turgor, triffid_nsc_ring_renewal
PUBLIC :: nsc_cmax_pft, nsc_wpot_pft

! Store (kg C m-2 of the PFT); < 0 = not yet set (filled to C_max at the
! first TRIFFID step). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_pft(:,:)
! The store's capacity C_max at the last store step (kg C m-2 of the PFT),
! for the sink feedback on photosynthesis (l_nsc_sink_feedback); < 0 = not
! yet set (no feedback until the first store step, also after a restart).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_cmax_pft(:,:)
! The release for growth at the last store step with the growth limits
! open (kg C m-2 of the PFT per 360 days), for the sink feedback: what the
! store would give growth if the turgor limit allowed.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_wpot_pft(:,:)
! Running mean of NPP (kg C m-2 per 360 days) and its weight, for W_max
! (bias-corrected mean = nsc_npp_mean / nsc_npp_wt). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_wt(:,:)
! Diagnostics of the last TRIFFID step (kg C m-2 per 360 days): the carbon
! released from the store (W, which TRIFFID grows from) and the net flux
! into the store (NPP - W).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_load_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_dc_pft(:,:)
! Overflow respiration of the last store step (kg C m-2 per 360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_resp_pft(:,:)
! Ring renewal (stem recovery), last TRIFFID step: the new ring, TRIFFID's
! gross wood production G (kg C m-2 of the PFT per 360 days; the wood
! TRIFFID grows, not extra; net of nsc_wood_lost), and the share of the
! stem conductance it renews, w_cond G / X (per 360 days). Diagnostics.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_prod_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_renew_pft(:,:)
! Ring renewal: the wood TRIFFID lost and has not yet regrown, plus the
! over-count carried from the last ring year, which the regrowth repays
! before new rings count (kg C m-2 of the PFT). In dumps with ring renewal.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_wood_lost(:,:)
! Ring renewal: the ring year's net wood growth so far (kg C m-2 of the
! PFT; < 0 when the year has lost wood), which settles the count at the end
! of the ring year. In dumps with ring renewal.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_net_wood_growth(:,:)
! Sapwood pool (l_nsc_sapwood_pool): the conducting sapwood X_s (kg C m-2
! of the PFT; < 0 = not yet set: the pipe model at the first step). In
! dumps.
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
!
! MDK, 30/10/26
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
IF (.NOT. ALLOCATED(nsc_xylem_prod_pft)) THEN
  ALLOCATE(nsc_xylem_prod_pft(n_land_pts, npft))
  nsc_xylem_prod_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_xylem_renew_pft)) THEN
  ALLOCATE(nsc_xylem_renew_pft(n_land_pts, npft))
  nsc_xylem_renew_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_wood_lost)) THEN
  ALLOCATE(nsc_wood_lost(n_land_pts, npft))
  nsc_wood_lost(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_net_wood_growth)) THEN
  ALLOCATE(nsc_net_wood_growth(n_land_pts, npft))
  nsc_net_wood_growth(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_sapwood)) THEN
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
IF (.NOT. ALLOCATED(nsc_cmax_pft)) THEN
  ALLOCATE(nsc_cmax_pft(n_land_pts, npft))
  nsc_cmax_pft(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_wpot_pft)) THEN
  ALLOCATE(nsc_wpot_pft(n_land_pts, npft))
  nsc_wpot_pft(:,:) = 0.0
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
! MDK, 30/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_store( land_pts, trif_pts, trif_index, r_gamma,         &
                              leafc, rootc, npp, pft_only )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE pftparm, ONLY: is_set, nsc_in_flush_window, nsc_limits_on
USE jules_vegetation_mod, ONLY: l_nsc_wmax_norm
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE pftparm, ONLY: g_leaf_0,                                                   &
                   nsc_flush, nsc_kc, nsc_reserve_frac, nsc_tau_load,          &
                   nsc_wmax_fac, nsc_g_starv,                                  &
                   nsc_psi_turg_zero, nsc_t_grow_zero, nsc_overflow_resp,      &
                   nsc_cmax, nsc_flush_start_doy, nsc_flush_days
USE pftparm, ONLY: nsc_grow_thresh
USE model_time_mod, ONLY: current_time

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points on which TRIFFID may operate.
  trif_pts,                                                                    &
  ! Indices of land points on which TRIFFID may operate.
  trif_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Inverse TRIFFID timestep (/360 days).
  r_gamma,                                                                     &
  ! Leaf carbon at the balanced LAI (kg C m-2).
  leafc(land_pts,npft),                                                        &
  ! Root carbon (kg C m-2).
  rootc(land_pts,npft)

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  ! IN NPP, OUT the carbon released from the store, W
  ! (kg C m-2 per 360 days).
  npp(land_pts,npft)

! Only this PFT (the daily step without TRIFFID); all PFTs if absent.
INTEGER, INTENT(IN), OPTIONAL :: pft_only

! Local variables
INTEGER :: n, t, l, n1, n2
REAL(KIND=real_jlslsm) ::                                                      &
  ! TRIFFID timestep (360-day years).
  dt,                                                                          &
  ! Weight of this step in the running mean of NPP.
  a_mean,                                                                      &
  ! Store capacity, reserve and accessible store (kg C m-2).
  c_max, c_res, c_acc,                                                         &
  ! Maximum and actual loading, net draw on the store (per 360 days).
  w_max, w, draw,                                                              &
  ! Store at the start of the step (kg C m-2).
  c_old,                                                                       &
  ! Turgor-limited sink: turgor factor (0-1).
  f_turg,                                                                      &
  ! The store the release for growth draws on: its floor (the reserve, or
  ! the growth threshold nsc_grow_thresh), the store above it, and that
  ! when full (kg C m-2).
  c_thr, c_grow, c_grow_max,                                                   &
  ! Phenology regrowth drawn from the store (per 360 days).
  w_regrow
! Growth limits on, leaf flush on.
LOGICAL :: l_turg, l_flush
! Day of year now; day counter, day of year and number of days of
! this store step.
INTEGER :: doy_now, k, d, nd
! Leaf flush: share of the year's flush in this step, flush draw (per
! 360 days).
REAL(KIND=real_jlslsm) :: f_flush, w_flush
! Days in this year.
INTEGER :: n_days

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

    ! Capacity: nsc_flush times the leaf and root carbon, or nsc_cmax set
    ! directly (e.g. measured NSC concentrations times biomass). The store
    ! starts full.
    c_max = nsc_flush(n) * (leafc(l,n) + rootc(l,n))
    IF (nsc_cmax(n) > 0.0) c_max = nsc_cmax(n)
    nsc_cmax_pft(l,n) = c_max
    c_res = nsc_reserve_frac(n) * c_max
    IF (nsc_pft(l,n) < 0.0) nsc_pft(l,n) = c_max
    c_old = nsc_pft(l,n)

    ! Maximum loading from the running mean of NPP (bias-corrected).
    nsc_npp_mean(l,n) = nsc_npp_mean(l,n) + a_mean * (npp(l,n)                 &
                        - nsc_npp_mean(l,n))
    nsc_npp_wt(l,n)   = nsc_npp_wt(l,n) + a_mean * (1.0 - nsc_npp_wt(l,n))
    w_max = 0.0
    IF (nsc_npp_wt(l,n) > 0.0)                                                 &
      w_max = nsc_wmax_fac(n) * MAX(nsc_npp_mean(l,n) / nsc_npp_wt(l,n), 0.0)

    ! l_nsc_wmax_norm: W_max over the running annual mean of the growth
    ! factor, so the release is concentrated in the growth windows and the
    ! annual release still matches the mean NPP.
    IF (l_nsc_wmax_norm .AND. l_turg .AND. nsc_ft_mean(l,n) > 0.0)             &
      w_max = w_max / MAX(nsc_ft_mean(l,n), 0.05)

    ! Leaf regrowth after drought shedding, queued by phenol (l_phen_psi),
    ! paid first from the accessible store (phenol limited it to that), as
    ! part of the release: over a shed-and-regrow cycle the store, not NPP,
    ! carries the cost.
    w_regrow = 0.0
    IF (nsc_regrow_acc(l,n) > 0.0) THEN
      w_regrow = MIN(nsc_regrow_acc(l,n), MAX(nsc_pft(l,n) - c_res, 0.0))
      nsc_pft(l,n) = nsc_pft(l,n) - w_regrow
      w_regrow = w_regrow / dt
      nsc_regrow_acc(l,n) = 0.0
    END IF

    ! Leaf flush: this step's share of the year's new leaves, g_leaf_0
    ! leafc, from the accessible store, as part of the release. The draw
    ! only moves carbon in time: over a year the release still balances
    ! NPP. Phenology regrowth paid above in the flush window counts towards
    ! it, so a deciduous canopy rebuilt by phenol is not paid twice: the
    ! larger of the two is charged.
    w_flush = 0.0
    IF (f_flush > 0.0) THEN
      w_flush = MIN(MAX(g_leaf_0(n) * leafc(l,n) * f_flush - w_regrow * dt,    &
                        0.0),                                                  &
                    MAX(nsc_pft(l,n) - c_res, 0.0))
      nsc_pft(l,n) = nsc_pft(l,n) - w_flush
      w_flush = w_flush / dt
    END IF

    ! Michaelis-Menten loading on the accessible store. The release for
    ! growth draws only on the store above the growth threshold, if set,
    ! as MEDFATE grows sapwood only from starch above half its capacity;
    ! the reserve still pays the deficits below.
    c_acc = MAX(nsc_pft(l,n) - c_res, 0.0)
    c_thr = c_res
    IF (is_set(nsc_grow_thresh(n)))                                            &
      c_thr = MAX(c_res, nsc_grow_thresh(n) * c_max)
    c_grow = MAX(nsc_pft(l,n) - c_thr, 0.0)

    ! With growth limits the release is scaled by the step's mean growth
    ! factor and normalised so a full store still releases W_max, whatever
    ! k_c: a softer release (larger nsc_kc) only shapes how loading falls
    ! as the store is drawn down.
    w = 0.0
    f_turg = 1.0
    IF (l_turg) THEN
      IF (is_set(nsc_turg_acc(l,n)))                                           &
        f_turg = MIN(MAX(nsc_turg_acc(l,n), 0.0), 1.0)
      nsc_turg_acc(l,n) = 0.0
      c_grow_max = c_max - c_thr
      IF (c_grow > 0.0 .AND. c_grow_max > 0.0)                                 &
        w = w_max * f_turg * c_grow / (nsc_kc(n) * c_max + c_grow)             &
            * (nsc_kc(n) * c_max + c_grow_max) / c_grow_max
      ! The release with the limits open, for the sink feedback.
      nsc_wpot_pft(l,n) = 0.0
      IF (f_turg > 0.0) nsc_wpot_pft(l,n) = w / f_turg
      IF (f_turg <= 0.0 .AND. c_grow > 0.0 .AND. c_grow_max > 0.0)             &
        nsc_wpot_pft(l,n) = w_max * c_grow / (nsc_kc(n) * c_max + c_grow)      &
                            * (nsc_kc(n) * c_max + c_grow_max) / c_grow_max
    ELSE
      IF (c_grow > 0.0) w = w_max * c_grow / (nsc_kc(n) * c_max + c_grow)
      nsc_wpot_pft(l,n) = w
    END IF

    ! Growth threshold: the step's release for growth can not take the
    ! store below the threshold (the deficits still draw to the reserve).
    IF (is_set(nsc_grow_thresh(n)))                                            &
      w = MIN(w, MAX(npp(l,n) + c_grow / dt, 0.0))

    ! The store can not be drawn below the reserve: a deficit beyond the
    ! accessible store passes to TRIFFID (w = npp < 0 when it is empty).
    draw = w - npp(l,n)
    IF (draw > 0.0) draw = MIN(draw, c_acc / dt)
    w = npp(l,n) + draw
    nsc_pft(l,n) = nsc_pft(l,n) - draw * dt

    ! The regrowth paid above is part of the release.
    w = w + w_regrow + w_flush

    ! Overflow above the capacity goes to growth, less the share
    ! nsc_overflow_resp respired: a sink-limited plant with a full store
    ! respires its surplus. nsc_resp_pft is added to the plant respiration
    ! and NPP outputs (extract_var) and leaves veg2's carbon check.
    nsc_resp_pft(l,n) = 0.0
    IF (nsc_pft(l,n) > c_max) THEN
      nsc_resp_pft(l,n) = nsc_overflow_resp(n) * (nsc_pft(l,n) - c_max) / dt
      w = w + (1.0 - nsc_overflow_resp(n)) * (nsc_pft(l,n) - c_max) / dt
      nsc_pft(l,n) = c_max
    END IF

    ! Carbon starvation (nsc_g_starv, with l_veg_compete): an area
    ! disturbance for lotka, at the full rate when the accessible store is
    ! used up and 0 for C_a >= nsc_kc C_max.
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
! is spread over the new area. The wood lost (nsc_wood_lost) and the ring
! year's net wood growth (nsc_net_wood_growth) are spread over area gained
! too: the new cover has lost no wood and grown none this ring year.
!
! MDK, 30/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_area( land_pts, trif_pts, trif_index, r_gamma,          &
                             frac, dfrac, lit_store )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points on which TRIFFID may operate.
  trif_pts,                                                                    &
  ! Indices of land points on which TRIFFID may operate.
  trif_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Inverse TRIFFID timestep (/360 days).
  r_gamma,                                                                     &
  ! Areal fractions after the TRIFFID step.
  frac(land_pts,ntype),                                                        &
  ! Change in the areal fractions over the step.
  dfrac(land_pts,npft)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Store carbon of the area lost, to litter (kg C m-2 of land per 360
  ! days).
  lit_store(land_pts,npft)

INTEGER :: n, t, l
! Areal fraction before the TRIFFID step.
REAL(KIND=real_jlslsm) :: frac_old

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
    frac_old = frac(l,n) - dfrac(l,n)

    ! Ring renewal: the new cover has lost no wood and grown none.
    IF (dfrac(l,n) > 0.0 .AND. frac(l,n) > 0.0 .AND.                           &
        ALLOCATED(nsc_wood_lost)) THEN
      nsc_wood_lost(l,n)       = nsc_wood_lost(l,n) * frac_old / frac(l,n)
      nsc_net_wood_growth(l,n) = nsc_net_wood_growth(l,n) * frac_old           &
                                 / frac(l,n)
    END IF

    ! The store: lost with the area lost, spread over the area gained.
    IF (nsc_pft(l,n) <= 0.0) CYCLE
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
! Growth limits on the store's release, called from sf_stom every time
! step for PFT ft. The growth factor of the time step,
!
!   f = f_psi f_T f_season,
!
! is added to nsc_turg_acc as a fraction of the store period (the TRIFFID
! period, or a day without TRIFFID), so the store step uses its mean.
!
! f_psi: growth needs turgor, so in drought the sink, not the supply,
! limits growth (Lempereur et al. 2015: stem growth stops at a predawn psi
! near -1.1 MPa while assimilation goes on). f_T: wood formation needs
! daily mean temperatures above ~5-8 degC (Rossi et al. 2008), so it uses
! a daily running mean of the air temperature and a sunny winter midday
! does not switch growth on. f_season: trees can be dormant in winter even
! when it is mild, so the store refills from winter NPP.
!
! Leaf renewal (l_ximpair_rec_leaf_turn) also runs here: new leaves arrive
! intact, so the leaf segment's loss falls by the new fraction of the
! canopy,
!
!   dR = g_leaf_0 f / <f> dt,
!
! with <f> the running annual mean of f (nsc_ft_mean): leaves are made only
! while growth is allowed, and a year's production matches g_leaf_0.
!
! MDK, 30/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_turgor( land_pts, veg_pts, veg_index, ft, timestep,     &
                               psi_rz, tair )

USE jules_vegetation_mod, ONLY: l_triffid, triffid_period,                     &
                                l_ximpair_rec_leaf_turn, l_nsc_wmax_norm
USE pftparm, ONLY: is_set, nsc_limits_on, nsc_season_on
USE pftparm, ONLY: nsc_psi_turg_full, nsc_psi_turg_zero, nsc_t_grow_zero,      &
                   nsc_t_grow_full, g_leaf_0, nsc_grow_start_doy,              &
                   nsc_grow_end_doy, nsc_flush_start_doy
USE model_time_mod, ONLY: current_time, is_spinup
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                   &
                                          ximpair_seg_state_alloc
USE water_constants_mod, ONLY: tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points of this PFT.
  veg_pts,                                                                     &
  ! Land points of this PFT.
  veg_index(land_pts),                                                         &
  ! Plant functional type.
  ft

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Model time step (s).
  timestep,                                                                    &
  ! Root-zone water potential (Pa).
  psi_rz(land_pts),                                                            &
  ! Air temperature (K).
  tair(land_pts)

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: period, f, a_t, a_yr, d_renew, f_season, doy,        &
                          doy_start
    ! Store period (s); growth factor (0-1); weights of this time step in
    ! the daily and annual running means; leaf renewal this time step;
    ! growing-season factor (0 or 1); day of year now and of the start of
    ! the release.
! Growth limits on; growing season set.
LOGICAL :: l_lim, l_season

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
f_season = 1.0
IF (l_season) THEN
  doy = REAL(day_of_year(current_time%year, current_time%month,                &
                         current_time%day, l_360, l_leap))
  ! With the leaf flush the release starts at budburst.
  doy_start = nsc_grow_start_doy(ft)
  IF (is_set(nsc_flush_start_doy(ft))) doy_start = nsc_flush_start_doy(ft)
  IF (is_set(doy_start) .AND. doy < doy_start) f_season = 0.0
  IF (is_set(nsc_grow_end_doy(ft)) .AND. doy >= nsc_grow_end_doy(ft))          &
    f_season = 0.0
END IF

DO j = 1,veg_pts
  l = veg_index(j)

  ! Turgor: the root-zone psi (MPa) between nsc_psi_turg_zero and
  ! nsc_psi_turg_full.
  f = 1.0
  IF (is_set(nsc_psi_turg_zero(ft)))                                           &
    f = MIN(MAX((1.0e-6 * psi_rz(l) - nsc_psi_turg_zero(ft))                   &
                / (nsc_psi_turg_full(ft) - nsc_psi_turg_zero(ft)), 0.0), 1.0)

  ! Cold: the daily running mean of the air temperature (degC) between
  ! nsc_t_grow_zero and nsc_t_grow_full.
  IF (is_set(nsc_t_grow_zero(ft))) THEN
    IF (nsc_t_run(l,ft) < 0.0) nsc_t_run(l,ft) = tair(l)
    nsc_t_run(l,ft) = nsc_t_run(l,ft) + a_t * (tair(l) - nsc_t_run(l,ft))
    f = f * MIN(MAX((nsc_t_run(l,ft) - tm - nsc_t_grow_zero(ft))               &
                    / (nsc_t_grow_full(ft) - nsc_t_grow_zero(ft)), 0.0), 1.0)
  END IF
  f = f * f_season

  ! The store step uses the mean over its period.
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
  ! leaf cap does not heal on its own either.
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
! Drought leaf shedding (l_phen_psi; Xu et al. 2016, New Phytol, Notes S1
! Eqn S7), called from sf_stom every time step for PFT ft. The daily mean
! root-zone psi (a predawn proxy) counts consecutive days below the turgor
! loss point phen_psi_tlp and above half of it; once either count reaches
! phen_ndays the leaf biomass factor phen_fbleaf falls or rises by
! phen_dfbleaf a day (between phen_fbleaf_min and 1). phenol caps the
! canopy at it and queues the regrowth's leaf carbon for the store
! (nsc_regrow_acc), so the store, not NPP, pays for the regrowth.
!
! MDK, 30/10/26
! ---------------------------------------------------------------------
SUBROUTINE phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,        &
                            psi_rz )

USE model_time_mod, ONLY: current_time
USE pftparm, ONLY: is_set
USE pftparm, ONLY: phen_psi_tlp, phen_dfbleaf, phen_ndays, phen_fbleaf_min

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points of this PFT.
  veg_pts,                                                                     &
  ! Land points of this PFT.
  veg_index(land_pts),                                                         &
  ! Plant functional type.
  ft

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Model time step (s).
  timestep,                                                                    &
  ! Root-zone water potential (Pa).
  psi_rz(land_pts)

INTEGER :: j, l
! Daily mean root-zone water potential (MPa).
REAL(KIND=real_jlslsm) :: psi_d

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
! root carbon from the prescribed height (TRIFFID allometry). The release
! leaves the store (there is no growth to feed), there is no stem recovery
! (it needs TRIFFID wood growth), and nsc_starv is a diagnostic only (the
! PFT fractions are fixed).
!
! MDK, 30/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_daily( land_pts, veg_pts, veg_index, ft, timestep,      &
                              ht, npp )

USE jules_surface_types_mod, ONLY: npft
USE model_time_mod, ONLY: current_time
USE calc_c_comps_triffid_mod, ONLY: calc_c_comps_triffid

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points of this PFT.
  veg_pts,                                                                     &
  ! Land points of this PFT.
  veg_index(land_pts),                                                         &
  ! Plant functional type.
  ft

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Model time step (s).
  timestep,                                                                    &
  ! Canopy height (m).
  ht(land_pts),                                                                &
  ! NPP of this time step (kg C m-2 s-1).
  npp(land_pts)

INTEGER :: j, l
! Vegetation carbon (kg C m-2), balanced LAI and wood carbon (kg C m-2)
! from the allometry.
REAL(KIND=real_jlslsm) :: c_veg, lai_b, woodc
! Inputs and outputs of triffid_nsc_store for this PFT (see there).
REAL(KIND=real_jlslsm), ALLOCATABLE :: leafc(:,:), rootc(:,:), npp_d(:,:)

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
! Stem recovery by ring renewal (l_ximpair_rec_nsc), called by triffid
! after vegcarb. The new ring is TRIFFID's gross wood production of the
! step (per 360 days),
!
!   G = dwood / dt + g_wood C_W,old,
!
! the wood increment plus the turnover replaced (wood_litc), counted net
! (below). The ring's share of the stem conductance,
!
!   rho = nsc_xylem_cond_ratio G dt / X,
!
! with X the conducting sapwood before the step's growth, is queued in
! ximpair_repair_pend for the next impairment update, which takes it off
! the stem's loss (l_nsc_ring_additive, the default) or dilutes the loss
! by it. With the sapwood pool X is the pool, which gains G and loses
! heartwood and its share of TRIFFID's wood turnover.
! Not in spin-up, where the impairment is frozen. No carbon moves here.
!
! MDK, 31/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_ring_renewal( land_pts, trif_pts, trif_index,           &
                                     r_gamma, ht, lai_bal, woodc, dwood,       &
                                     wood_litc )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE jules_vegetation_mod, ONLY: l_ximpair_rec_nsc, xylem_impairment_memory
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_pool
USE pftparm, ONLY: pft_xylem_impairment_model, nsc_xylem_cond_ratio
USE pftparm, ONLY: is_set, nsc_grow_end_doy, nsc_heartwood_rate
USE calc_c_comps_triffid_mod, ONLY: sapwood_c
USE model_time_mod, ONLY: is_spinup
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE xylem_impairment_memory_mod, ONLY: ximpair_repair_pend,                    &
                                       ximpair_memory_alloc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of points on which TRIFFID may operate.
  trif_pts,                                                                    &
  ! Indices of land points on which TRIFFID may operate.
  trif_index(land_pts)

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Inverse TRIFFID timestep (/360 days).
  r_gamma,                                                                     &
  ! Canopy height before the step's growth (m).
  ht(land_pts,npft),                                                           &
  ! Balanced LAI before the step's growth.
  lai_bal(land_pts,npft),                                                      &
  ! Wood carbon after the step's growth (kg C m-2).
  woodc(land_pts,npft),                                                        &
  ! Wood increment of the step (kg C m-2).
  dwood(land_pts,npft),                                                        &
  ! Wood turnover of the step, g_wood C_W (kg C m-2 per 360 days).
  wood_litc(land_pts,npft)

INTEGER :: n, t, l
! Gross wood production (kg C m-2 per 360 days); conducting sapwood
! (kg C m-2).
REAL(KIND=real_jlslsm) :: g_wood_prod, x_s
! TRIFFID timestep (360-day years) and the step's gross wood
! production (kg C m-2, < 0 when wood is lost).
REAL(KIND=real_jlslsm) :: dt, g_step
! Sapwood pool: the wood before the step (kg C m-2) and TRIFFID's wood
! turnover rate (per 360 days).
REAL(KIND=real_jlslsm) :: wood_old, g_turn
! Day of year now, days in this year, days of this step, day counter,
! a day of the step and the ring year's end day.
INTEGER :: doy_now, n_days, nd, k, d, doy_end
! This step ends the ring year (it holds the end of the growing season
! or, without one, the last day of the year).
LOGICAL :: l_ring_end

! Minimum conducting xylem (kg C m-2): a smaller sapwood renews
! nothing.
REAL(KIND=real_jlslsm), PARAMETER :: x_min = 1.0e-6

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_RING_RENEWAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_ximpair_rec_nsc .AND. .NOT. is_spinup) THEN

  CALL triffid_nsc_alloc( land_pts )
  CALL ximpair_memory_alloc( land_pts )

  ! The days of this step: the nd days up to today, as the leaf flush in
  ! triffid_nsc_store counts them.
  dt      = 1.0 / r_gamma
  n_days  = days_in_year(current_time%year, l_360, l_leap)
  doy_now = day_of_year(current_time%year, current_time%month,                 &
                        current_time%day, l_360, l_leap)
  nd      = MAX(NINT(dt * 360.0), 1)

  DO n = 1,nnpft
    IF (pft_xylem_impairment_model(n) /= xylem_impairment_memory) CYCLE

    ! The ring year ends at the step that holds the end of the growing
    ! season (nsc_grow_end_doy) or, without one, the last day of the year.
    doy_end = n_days
    IF (is_set(nsc_grow_end_doy(n)))                                           &
      doy_end = MIN(NINT(nsc_grow_end_doy(n)), n_days)
    l_ring_end = .FALSE.
    DO k = 0,nd-1
      d = MODULO(doy_now - k - 1, n_days) + 1
      IF (d == doy_end) l_ring_end = .TRUE.
    END DO

    DO t = 1,trif_pts
      l = trif_index(t)
      nsc_xylem_prod_pft(l,n)  = 0.0
      nsc_xylem_renew_pft(l,n) = 0.0

      ! The conducting sapwood X from the allometry, or with the sapwood
      ! pool the prognostic X_s, which starts (where it is not yet set) at
      ! its steady state for the wood before the step,
      !
      !   X_s = C_W g_t / (k_h + g_t),   g_t = wood_litc / C_W.
      x_s = sapwood_c(n, ht(l,n), lai_bal(l,n), woodc(l,n) - dwood(l,n))
      IF (l_nsc_sapwood_pool) THEN
        wood_old = woodc(l,n) - dwood(l,n)
        g_turn   = 0.0
        IF (wood_old > 0.0) g_turn = wood_litc(l,n) / wood_old
        IF (nsc_sapwood(l,n) < 0.0)                                            &
          nsc_sapwood(l,n) = MAX(wood_old, 0.0) * g_turn                       &
                             / (nsc_heartwood_rate(n) + g_turn)
        x_s = nsc_sapwood(l,n)
      END IF

      ! Net wood growth. While the store's release is shut TRIFFID still
      ! charges the turnover, so the wood shrinks (G < 0) and regrows when
      ! the release opens. A step with G < 0 adds its loss to the deficit
      ! D and renews nothing; a step with G > 0 repays D first, and only
      ! the rest is new ring. S sums the ring year's G dt.
      g_step = dwood(l,n) + wood_litc(l,n) * dt
      IF (g_step < 0.0) THEN
        nsc_wood_lost(l,n) = nsc_wood_lost(l,n) - g_step
        g_wood_prod = 0.0
      ELSE
        g_wood_prod = MAX(g_step - nsc_wood_lost(l,n), 0.0) * r_gamma
        nsc_wood_lost(l,n) = MAX(nsc_wood_lost(l,n) - g_step, 0.0)
      END IF
      nsc_net_wood_growth(l,n)       = nsc_net_wood_growth(l,n) + g_step
      nsc_xylem_prod_pft(l,n) = g_wood_prod

      ! The ring's share of the stem conductance, queued for the next
      ! impairment update. A sapwood at or below x_min renews nothing, but
      ! its count and state go on.
      IF (x_s > x_min) THEN
        nsc_xylem_renew_pft(l,n) = nsc_xylem_cond_ratio(n) * g_wood_prod      &
                                   / x_s
        ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n)                    &
                                   + nsc_xylem_renew_pft(l,n) / r_gamma
      END IF

      ! Sapwood pool: the oldest sapwood becomes heartwood, TRIFFID's
      ! wood turnover takes its share and the new ring is added,
      !
      !   dX_s/dt = G - (k_h + g_t) X_s,   k_h = nsc_heartwood_rate,
      !
      ! so in equilibrium X_s = G / (k_h + g_t).
      IF (l_nsc_sapwood_pool)                                                  &
        nsc_sapwood(l,n) = nsc_sapwood(l,n)                                    &
                           * EXP(-(nsc_heartwood_rate(n) + g_turn) / r_gamma)  &
                           + g_wood_prod / r_gamma

      ! End of the ring year: its net production S settles the count,
      !
      !   D <- MAX(D + MIN(S, 0), 0),   S <- 0,
      !
      ! so D keeps only the over-count (growth counted before a loss that
      ! was not regrown) for next year's growth to repay, and a lasting
      ! loss of wood does not hold back later rings. With the sapwood
      ! pool a net loss of wood also leaves the sapwood, |S| (all rings
      ! alike).
      IF (l_ring_end) THEN
        nsc_wood_lost(l,n) = MAX(nsc_wood_lost(l,n)                            &
                                    + MIN(nsc_net_wood_growth(l,n), 0.0), 0.0)
        IF (l_nsc_sapwood_pool)                                                &
          nsc_sapwood(l,n) = MAX(nsc_sapwood(l,n)                              &
                                 + MIN(nsc_net_wood_growth(l,n), 0.0), 0.0)
        nsc_net_wood_growth(l,n) = 0.0
      END IF

    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_ring_renewal


END MODULE triffid_nsc_mod
