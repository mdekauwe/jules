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
!           running mean is the JULES stand-in, as TRIFFID has no tree size.
! The store is held between the reserve and C_max: carbon above C_max goes to
! growth (overflow), so a well-supplied plant grows from its NPP as without
! the store. When the accessible store is used up the carbon deficit passes
! to TRIFFID as negative growth (the vegetation shrinks), which is the JULES
! form of their starvation.
!
! Xylem repair (l_ximpair_rec_nsc, with the memory impairment model): the
! stem's lost conductance recovers only by regrowing xylem paid for from the
! store. With the loss of conductance PLC = ximpair_lock, the functional
! xylem is f_x = 1 - PLC of X = nsc_sapwood_frac woodc, the conducting
! sapwood (Trugman et al. Table S2: sapwood fraction SA = 0.6 of the stem
! biomass; 1 / a_ws would be the JULES live stemwood, eta_sl h LAI_bal,
! which at FR-Pue is about a sixth of that), or, with nsc_la_sa > 0, the
! pipe model on the balanced leaf area as in FATES, X = LAI_bal /
! nsc_la_sa h nsc_wood_density cmass (from measured leaf:sapwood area and
! wood density, independent of the TRIFFID wood), and the leaves and fine roots
! the plant keeps are L = f_x (leaf + root) (Trugman et al. Eqn 2, leaf area
! in proportion to functional xylem; l_ximpair_leaf_loss sheds the canopy to
! match). Their optimal share of W rebuilding xylem, with no xylem turnover
! (their m_X = 0), is
!   u = X / (X + L) (1 - m_L L / W)                (Eqn 5)
! with m_L L the leaf and fine-root turnover (g_leaf_0, g_root) of what is
! kept. Here m_L L is the turnover TRIFFID charges against the growth W
! - repair: the leaves kept, (1 - PLC) g_leaf_0 leafc, all the roots,
! g_root rootc (TRIFFID does not shed them), and the wood, g_wood woodc
! (the repair carbon is litter, not wood, so wood turnover is not paid by
! the repair; their m_X = 0). Only the loss the current water status does
! not re-impose is repaired, PLC_r = PLC - (1 - f_min), with f_min the
! lowest intact stem conductance fraction seen (ximpair_fmin_acc), whose
! memory decays back to intact with e-folding time nsc_repair_window days
! (0 = since the last store step). Loss up to nsc_plc_background is
! background embolism the store does not pay for: the stem is rebuilt to
! that loss, not to 0. The repair, min(u W dt, PLC_r X), restores
! conductance in proportion,
! dPLC = -repair / X, at the next xylem impairment update. TRIFFID's
! allometry fixes the wood carbon, so the repair carbon leaves as litter
! (the embolised conduits it replaces are retired). With
! ximpair_psi_growth set, the repair is scaled by the fraction of time since
! the last TRIFFID step that the root zone was wetter than it (no stem growth
! in drought; ximpair_grow_acc).
!
! Wood anatomy (nsc_xylem_type). 1, diffuse-porous (the default, as above):
! repair whenever growth is allowed, on all the conducting xylem X, at most
! nsc_xylem_rate_max f_t X per year (f_t the growth factor of the turgor
! and cold limits, 1 without them; the cambium builds a ring a year, so a
! few % of the sapwood a year; 0 = no cap). 2, ring-porous: water moves
! mostly in the current earlywood ring, X_r = nsc_xylem_rate_max X, built
! each spring from the store. Once a year, at the first store step on or
! after nsc_ring_mmdd, the repairable loss is rebuilt in full, at a cost of
! PLC_r X_r drawn from the store (as part of W, to litter), or in part if
! the accessible store is short; there is no repair the rest of the year.
!
! W_max is held while the stem is damaged (ximpair_lock > nsc_hold_plc), as
! Trugman et al. hold it through recovery: the damaged plant keeps loading
! for the intact sink, which is what runs the store down.
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
! MDK, 20/10/26
!
! Leaf flush (nsc_flush_start_doy set): leaves flush before stem growth
! starts and the flush draws on the store (Q. ilex at Puechabon: budburst
! about day 105, leaves open about day 127, Limousin et al. 2012; full size
! 6-7 weeks after budburst, Misson et al. 2011), so
! the store's release starts at budburst, not at nsc_grow_start_doy, and
! the year's new leaf carbon, g_leaf_0 x leaf carbon, is drawn from the
! accessible store evenly over nsc_flush_days days (part of W, like the
! regrowth after drought shedding). The draw only moves when carbon leaves
! the store: over a year the release still balances NPP. The growing
! season then gates only the cambium: xylem repair and its cap are scaled
! by the share of the store step's days inside it.
! MDK, 20/10/26
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
! applied after all inputs to the store (NPP and the optional extra input
! c_in). nsc_resp_pft is added to the plant respiration and NPP outputs
! (extract_var, lagged by one store step), is part of resp_p_actual through
! nsc_dc, and leaves the carbon-conservation check of veg2.
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
! (there is no growth to feed), the repair rebuilds the stem cap as with
! TRIFFID, and nsc_starv is a diagnostic only (the PFT fractions are fixed).
!
! MDK, 17/10/26
!
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: nsc_pft, nsc_npp_mean, nsc_npp_wt, nsc_load_pft, nsc_repair_pft,   &
          nsc_dc_pft, nsc_starv_pft, nsc_npp_day, nsc_turg_acc, nsc_resp_pft, &
          nsc_t_run, nsc_ft_mean, phen_psi_acc, phen_nlow, phen_nhigh,        &
          phen_fbleaf, nsc_regrow_acc, phen_psi_update,                       &
          triffid_nsc_alloc, triffid_nsc_store, triffid_nsc_area,             &
          triffid_nsc_daily, triffid_nsc_turgor

! Store (kg C m-2 of the PFT); < 0 = not yet set (filled to C_max at the
! first TRIFFID step). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_pft(:,:)
! Running mean of NPP (kg C m-2 per 360 days) and its weight, for W_max
! (bias-corrected mean = nsc_npp_mean / nsc_npp_wt). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_wt(:,:)
! Diagnostics of the last TRIFFID step (kg C m-2 per 360 days): the carbon
! released from the store (W, which TRIFFID grows from), the part of it
! that rebuilt xylem, and the net flux into the store (NPP - W).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_load_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_repair_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_dc_pft(:,:)
! Overflow respiration of the last store step (kg C m-2 per 360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_resp_pft(:,:)
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

INTEGER, INTENT(IN) :: n_land_pts

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
IF (.NOT. ALLOCATED(nsc_repair_pft)) THEN
  ALLOCATE(nsc_repair_pft(n_land_pts, npft))
  nsc_repair_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_dc_pft)) THEN
  ALLOCATE(nsc_dc_pft(n_land_pts, npft))
  nsc_dc_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_resp_pft)) THEN
  ALLOCATE(nsc_resp_pft(n_land_pts, npft))
  nsc_resp_pft(:,:) = 0.0
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
  ALLOCATE(phen_psi_acc(n_land_pts, npft), phen_nlow(n_land_pts, npft),      &
           phen_nhigh(n_land_pts, npft), phen_fbleaf(n_land_pts, npft),      &
           nsc_regrow_acc(n_land_pts, npft))
  phen_psi_acc(:,:)   = 0.0
  phen_nlow(:,:)      = 0.0
  phen_nhigh(:,:)     = 0.0
  phen_fbleaf(:,:)    = 1.0
  nsc_regrow_acc(:,:) = 0.0
END IF

END SUBROUTINE triffid_nsc_alloc

! ---------------------------------------------------------------------
! One TRIFFID step of the store. On return npp holds W, the carbon released
! from the store (kg C m-2 per 360 days), and repair the part of W that
! rebuilds xylem (TRIFFID grows from W - repair; the repair goes to
! litter). The repaired fraction of conductance is queued in
! ximpair_repair_pend for the next xylem impairment update.
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_store( land_pts, trif_pts, trif_index, r_gamma,       &
                              woodc, leafc, rootc, ht, lai_bal, npp, repair,  &
                              pft_only, c_in )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE jules_vegetation_mod, ONLY: l_ximpair_rec_nsc, xylem_impairment_memory,  &
                                l_triffid, l_nsc_wmax_norm
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE calc_c_comps_triffid_mod, ONLY: sapwood_c   ! MDK, 20/10/26
USE pftparm, ONLY: g_leaf_0, pft_xylem_impairment_model,                     &
                   nsc_flush, nsc_kc, nsc_reserve_frac, nsc_tau_load,        &
                   nsc_wmax_fac, nsc_sapwood_frac, nsc_hold_plc,             &
                   nsc_g_starv, ximpair_psi_growth, nsc_repair_window,       &
                   nsc_plc_background, nsc_la_sa, nsc_wood_density,          &
                   nsc_psi_turg_zero, nsc_t_grow_zero, nsc_overflow_resp,    &
                   nsc_xylem_type, nsc_xylem_rate_max, nsc_ring_mmdd,        &
                   nsc_cmax, nsc_grow_start_doy, nsc_grow_end_doy,           &
                   nsc_flush_start_doy, nsc_flush_days
USE jules_surface_mod, ONLY: cmass
USE model_time_mod, ONLY: current_time
USE trif, ONLY: g_root, g_wood
USE xylem_impairment_memory_mod, ONLY: ximpair_lock, ximpair_repair_pend,    &
                                       ximpair_grow_acc, ximpair_fmin_acc,   &
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
  woodc(land_pts,npft),                                                        &
    ! Wood carbon (kg C m-2).
  leafc(land_pts,npft),                                                        &
    ! Leaf carbon at the balanced LAI (kg C m-2).
  rootc(land_pts,npft),                                                        &
    ! Root carbon (kg C m-2).
  ht(land_pts,npft),                                                           &
    ! Canopy height (m).
  lai_bal(land_pts,npft)
    ! Balanced LAI (for the pipe-model xylem, nsc_la_sa > 0).

REAL(KIND=real_jlslsm), INTENT(IN OUT) ::                                      &
  npp(land_pts,npft)
    ! IN NPP, OUT the carbon released from the store, W
    ! (kg C m-2 per 360 days).

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  repair(land_pts,npft)
    ! Carbon rebuilding xylem (kg C m-2 per 360 days), part of W.

INTEGER, INTENT(IN), OPTIONAL :: pft_only
    ! Only this PFT (the daily step without TRIFFID); all PFTs if absent.

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: c_in(land_pts,npft)
    ! Extra carbon into the store besides NPP (kg C m-2 per 360 days), taken
    ! before the cap, so it can overflow.

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
  plc, x_s, l_tot, turn_l, u,                                                 &
  f_grow, plc_rep,                                                             &
  f_turg, c_acc_max, w_regrow
    ! Repair: loss of conductance, conducting xylem (kg C m-2), leaf and fine
    ! root of the intact plant (kg C m-2), turnover TRIFFID charges against
    ! growth (per 360 days), share of W, growth gate, repairable loss.
    ! Turgor-limited sink: turgor factor, accessible store when full (kg C
    ! m-2).
LOGICAL :: l_rep, l_mem, l_turg, l_ring_now, l_flush
INTEGER :: doy_now, doy_ring, k, d, nd
REAL(KIND=real_jlslsm) :: f_flush, f_stem, w_flush
    ! Leaf flush: share of the year's flush in this step, share of the
    ! step's days in the growing season (cambium), flush draw (per 360 days).
REAL(KIND=real_jlslsm) :: x_r
    ! Ring-porous: this store step holds the earlywood rebuild; day of year
    ! now and of nsc_ring_mmdd; earlywood ring (kg C m-2).
INTEGER :: n_days

REAL(KIND=real_jlslsm), PARAMETER :: x_min = 1.0e-6
    ! Minimum conducting xylem for repair (kg C m-2).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_NSC_STORE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL triffid_nsc_alloc( land_pts )

dt = 1.0 / r_gamma
repair(:,:) = 0.0
n1 = 1
n2 = nnpft
IF (PRESENT(pft_only)) THEN
  n1 = pft_only
  n2 = pft_only
END IF

DO n = n1,n2

  l_mem = pft_xylem_impairment_model(n) == xylem_impairment_memory
  l_rep = l_ximpair_rec_nsc .AND. l_mem
  IF (l_rep) CALL ximpair_memory_alloc( land_pts )

  a_mean = MIN(dt / MAX(nsc_tau_load(n), TINY(1.0_real_jlslsm)), 1.0)
  l_turg = nsc_psi_turg_zero(n) > -1.0e29 .OR. nsc_t_grow_zero(n) > -1.0e29 &
           .OR. nsc_grow_start_doy(n) > -1.0e29                          &
           .OR. nsc_grow_end_doy(n) > -1.0e29                            &
           .OR. nsc_flush_start_doy(n) > -1.0e29
  ! Leaf flush: days of this store step (the nd days up to today) in the
  ! flush and in the growing season.
  ! MDK, 20/10/26
  l_flush = nsc_flush_start_doy(n) > -1.0e29
  f_flush = 0.0
  f_stem  = 1.0
  IF (l_flush) THEN
    n_days  = days_in_year(current_time%year, l_360, l_leap)
    doy_now = day_of_year(current_time%year, current_time%month,             &
                          current_time%day, l_360, l_leap)
    nd = MAX(NINT(dt * 360.0), 1)
    f_stem = 0.0
    DO k = 0,nd-1
      d = MODULO(doy_now - k - 1, n_days) + 1
      ! The flush window may cross the year end.
      IF (MODULO(d - NINT(nsc_flush_start_doy(n)), n_days)                    &
          < NINT(nsc_flush_days(n)))                                          &
        f_flush = f_flush + 1.0 / REAL(NINT(nsc_flush_days(n)))
      ! Only the start of the stem window: its end is already in the
      ! growth factor (the release season).
      IF (nsc_grow_start_doy(n) < -1.0e29 .OR.                                &
          REAL(d) >= nsc_grow_start_doy(n))                                   &
        f_stem = f_stem + 1.0 / REAL(nd)
    END DO
  END IF
  ! Ring-porous: is the earlywood date within this store step?
  l_ring_now = .FALSE.
  IF (nsc_xylem_type(n) == 2) THEN
    n_days   = days_in_year(current_time%year, l_360, l_leap)
    doy_now  = day_of_year(current_time%year, current_time%month,            &
                           current_time%day, l_360, l_leap)
    IF (l_360) THEN
      doy_ring = day_of_year(current_time%year, nsc_ring_mmdd(n) / 100,      &
                             MIN(MOD(nsc_ring_mmdd(n), 100), 30), l_360,     &
                             l_leap)
    ELSE
      doy_ring = day_of_year(current_time%year, nsc_ring_mmdd(n) / 100,      &
                             MOD(nsc_ring_mmdd(n), 100), l_360, l_leap)
    END IF
    doy_ring = MIN(doy_ring, n_days)
    ! Within this store step, wrapping over the new year.
    l_ring_now = MODULO(doy_now - doy_ring, n_days) < MAX(NINT(dt * 360.0), 1)
  END IF

  DO t = 1,trif_pts
    l = trif_index(t)

    c_max = nsc_flush(n) * (leafc(l,n) + rootc(l,n))
    IF (nsc_cmax(n) > 0.0) c_max = nsc_cmax(n)
    c_res = nsc_reserve_frac(n) * c_max
    IF (nsc_pft(l,n) < 0.0) nsc_pft(l,n) = c_max
    c_old = nsc_pft(l,n)

    ! Stem loss of conductance (0 before the first impairment update).
    plc = 0.0
    IF (l_mem .AND. ALLOCATED(ximpair_lock))                                   &
      plc = MIN(MAX(ximpair_lock(l,n), 0.0), 1.0)

    ! Maximum loading from the running mean of NPP (bias-corrected), held
    ! while the stem is damaged.
    IF (plc <= nsc_hold_plc(n)) THEN
      nsc_npp_mean(l,n) = nsc_npp_mean(l,n) + a_mean * (npp(l,n)            &
                          - nsc_npp_mean(l,n))
      nsc_npp_wt(l,n)   = nsc_npp_wt(l,n) + a_mean * (1.0 - nsc_npp_wt(l,n))
    END IF
    w_max = 0.0
    IF (nsc_npp_wt(l,n) > 0.0)                                                 &
      w_max = nsc_wmax_fac(n) * MAX(nsc_npp_mean(l,n) / nsc_npp_wt(l,n), 0.0)
    ! Release concentrated in the growth windows.
    IF (l_nsc_wmax_norm .AND. l_turg .AND. nsc_ft_mean(l,n) > 0.0)            &
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
    ! MDK, 20/10/26
    ! Phenology regrowth paid above in the flush window counts towards
    ! it, so a deciduous canopy rebuilt by phenol is not paid twice: the
    ! larger of the two is charged.
    w_flush = 0.0
    IF (f_flush > 0.0) THEN
      w_flush = MIN(MAX(g_leaf_0(n) * leafc(l,n) * f_flush - w_regrow * dt,  &
                        0.0),                                                 &
                    MAX(nsc_pft(l,n) - c_res, 0.0))
      nsc_pft(l,n) = nsc_pft(l,n) - w_flush
      w_flush = w_flush / dt
    END IF

    ! Michaelis-Menten loading on the accessible store.
    c_acc = MAX(nsc_pft(l,n) - c_res, 0.0)
    w = 0.0
    f_turg = 1.0
    IF (l_turg) THEN
      ! Turgor-limited sink, normalised so the full store releases W_max.
      IF (nsc_turg_acc(l,n) > -1.0e29)                                         &
        f_turg = MIN(MAX(nsc_turg_acc(l,n), 0.0), 1.0)
      nsc_turg_acc(l,n) = 0.0
      c_acc_max = c_max - c_res
      IF (c_acc > 0.0 .AND. c_acc_max > 0.0)                                   &
        w = w_max * f_turg * c_acc / (nsc_kc(n) * c_max + c_acc)               &
            * (nsc_kc(n) * c_max + c_acc_max) / c_acc_max
    ELSE
      IF (c_acc > 0.0) w = w_max * c_acc / (nsc_kc(n) * c_max + c_acc)
    END IF

    ! The store can not be drawn below the reserve: a deficit beyond the
    ! accessible store passes to TRIFFID (w = npp < 0 when it is empty).
    draw = w - npp(l,n)
    IF (draw > 0.0) draw = MIN(draw, c_acc / dt)
    w = npp(l,n) + draw
    nsc_pft(l,n) = nsc_pft(l,n) - draw * dt

    ! The regrowth paid above is part of the release.
    w = w + w_regrow + w_flush

    ! Extra input, before the cap.
    IF (PRESENT(c_in)) nsc_pft(l,n) = nsc_pft(l,n) + MAX(c_in(l,n), 0.0) * dt

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

    ! Xylem repair from the released carbon, only while the root zone is
    ! wetter than ximpair_psi_growth (if set).
    f_grow  = 1.0
    plc_rep = 0.0
    IF (l_rep) THEN
      IF (ximpair_psi_growth(n) > -1.0e29)                                     &
        f_grow = MIN(MAX(ximpair_grow_acc(l,n), 0.0), 1.0)
      ! With the leaf flush, only in the growing season (cambium).
      f_grow = f_grow * f_stem
      ximpair_grow_acc(l,n) = 0.0
      ! Only the loss the water status of the recent past does not
      ! re-impose. The memory of the worst stem conductance decays back to
      ! intact with e-folding time nsc_repair_window (days; 0 = reset).
      plc_rep = MAX(plc - (1.0 - MIN(MAX(ximpair_fmin_acc(l,n), 0.0), 1.0)),   &
                    0.0)
      ! Loss up to nsc_plc_background is background embolism, not repaid.
      plc_rep = MAX(MIN(plc_rep, plc - nsc_plc_background(n)), 0.0)
      IF (nsc_repair_window(n) > 0.0) THEN
        ximpair_fmin_acc(l,n) = 1.0 - (1.0 - ximpair_fmin_acc(l,n))            &
                                * EXP(-dt * 360.0 / nsc_repair_window(n))
      ELSE
        ximpair_fmin_acc(l,n) = 1.0
      END IF
    END IF
    IF (l_rep .AND. (w > 0.0 .OR. nsc_xylem_type(n) == 2)) THEN
      ! Pipe model on the balanced leaf area, or a fraction of the wood.
      x_s = sapwood_c(n, ht(l,n), lai_bal(l,n), woodc(l,n))   ! MDK, 20/10/26
      IF (nsc_xylem_type(n) == 2) THEN
        ! Ring-porous: the earlywood rebuild, once a year, from the store
        ! (not gated by the turgor or growth limits, nor capped: the new
        ! ring is the conducting system, built at the spring date).
        x_r = nsc_xylem_rate_max(n) * x_s
        IF (l_ring_now .AND. plc_rep > 0.0 .AND. x_r > x_min) THEN
          repair(l,n) = MIN(plc_rep * x_r,                                     &
                            MAX(nsc_pft(l,n) - c_res, 0.0)) / dt
          nsc_pft(l,n) = nsc_pft(l,n) - repair(l,n) * dt
          w = w + repair(l,n)
          ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n)                  &
                                     + repair(l,n) * dt / x_r
        END IF
      ELSE IF (plc_rep > 0.0 .AND. x_s > x_min) THEN
        ! X and L both scale with f_x, so X / (X + L) is that of the intact
        ! plant; the turnover is what TRIFFID charges against the growth
        ! W - repair (leaves kept, all roots, wood).
        l_tot  = leafc(l,n) + rootc(l,n)
        IF (l_triffid) THEN
          turn_l = (1.0 - plc) * g_leaf_0(n) * leafc(l,n)                      &
                   + g_root(n) * rootc(l,n) + g_wood(n) * woodc(l,n)
        ELSE
          ! Without TRIFFID (no TRIFFID turnover to keep up; the trif
          ! parameters are not read): their leaf and fine-root turnover of
          ! what is kept, at g_leaf_0, and m_X = 0.
          turn_l = (1.0 - plc) * g_leaf_0(n) * l_tot
        END IF
        ! The leaf flush's carbon is for leaves: not shared with xylem.
        u = x_s / (x_s + l_tot) * MAX(1.0 - turn_l / MAX(w - w_flush,        &
                                                         TINY(1.0)), 0.0)
        repair(l,n) = MIN(f_grow * u * MAX(w - w_flush, 0.0),                 &
                          plc_rep * x_s / dt)
        ! Cambial production limit.
        IF (nsc_xylem_rate_max(n) > 0.0)                                       &
          repair(l,n) = MIN(repair(l,n), nsc_xylem_rate_max(n) * f_turg      &
                                         * f_stem * x_s)
        ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n)                    &
                                   + repair(l,n) * dt / x_s
      END IF
    END IF

    npp(l,n)            = w
    nsc_load_pft(l,n)   = w
    nsc_repair_pft(l,n) = repair(l,n)
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
SUBROUTINE triffid_nsc_area( land_pts, trif_pts, trif_index, r_gamma,        &
                             frac, dfrac, lit_store )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype

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

lit_store(:,:) = 0.0
IF (.NOT. ALLOCATED(nsc_pft)) RETURN

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

END SUBROUTINE triffid_nsc_area

! ---------------------------------------------------------------------
! Turgor-limited sink: called from sf_stom every time step for PFT ft, adds
! this time step's turgor factor, from the root-zone water potential, to
! nsc_turg_acc as a fraction of the store period (the TRIFFID period, or a
! day without TRIFFID).
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_turgor( land_pts, veg_pts, veg_index, ft, timestep,   &
                               psi_rz, tair )

USE jules_vegetation_mod, ONLY: l_triffid, triffid_period,                    &
                                l_ximpair_rec_leaf_turn, l_nsc_wmax_norm
USE pftparm, ONLY: nsc_psi_turg_full, nsc_psi_turg_zero, nsc_t_grow_zero,   &
                   nsc_t_grow_full, g_leaf_0, nsc_grow_start_doy,             &
                   nsc_grow_end_doy, nsc_flush_start_doy
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                  &
                                          ximpair_seg_state_alloc
USE water_constants_mod, ONLY: tm

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
REAL(KIND=real_jlslsm) :: period, f, a_t, a_yr, d_renew, f_season, doy,     &
                          doy_start
LOGICAL :: l_lim, l_season

l_season = nsc_grow_start_doy(ft) > -1.0e29 .OR. nsc_grow_end_doy(ft) > -1.0e29 &
           .OR. nsc_flush_start_doy(ft) > -1.0e29
l_lim = nsc_psi_turg_zero(ft) > -1.0e29 .OR. nsc_t_grow_zero(ft) > -1.0e29   &
        .OR. l_season
IF (.NOT. l_lim .AND. .NOT. l_ximpair_rec_leaf_turn) RETURN
CALL triffid_nsc_alloc( land_pts )
IF (l_ximpair_rec_leaf_turn) CALL ximpair_seg_state_alloc( land_pts )

period = 86400.0
IF (l_triffid) period = 86400.0 * REAL(triffid_period)
a_t = 1.0 - EXP(-timestep / 86400.0)
a_yr = 1.0 - EXP(-timestep / (360.0 * 86400.0))
! Growing season: 0 outside the start and end days of the year.
! MDK, 20/10/26
f_season = 1.0
IF (l_season) THEN
  doy = REAL(day_of_year(current_time%year, current_time%month,              &
                         current_time%day, l_360, l_leap))
  ! With the leaf flush the release starts at budburst.
  ! MDK, 20/10/26
  doy_start = nsc_grow_start_doy(ft)
  IF (nsc_flush_start_doy(ft) > -1.0e29) doy_start = nsc_flush_start_doy(ft)
  IF (doy_start > -1.0e29 .AND. doy < doy_start) f_season = 0.0
  IF (nsc_grow_end_doy(ft) > -1.0e29 .AND. doy >= nsc_grow_end_doy(ft))      &
    f_season = 0.0
END IF

DO j = 1,veg_pts
  l = veg_index(j)
  f = 1.0
  IF (nsc_psi_turg_zero(ft) > -1.0e29)                                         &
    f = MIN(MAX((1.0e-6 * psi_rz(l) - nsc_psi_turg_zero(ft))                   &
                / (nsc_psi_turg_full(ft) - nsc_psi_turg_zero(ft)), 0.0), 1.0)
  IF (nsc_t_grow_zero(ft) > -1.0e29) THEN
    ! Daily running mean of the air temperature.
    IF (nsc_t_run(l,ft) < 0.0) nsc_t_run(l,ft) = tair(l)
    nsc_t_run(l,ft) = nsc_t_run(l,ft) + a_t * (tair(l) - nsc_t_run(l,ft))
    f = f * MIN(MAX((nsc_t_run(l,ft) - tm - nsc_t_grow_zero(ft))               &
                    / (nsc_t_grow_full(ft) - nsc_t_grow_zero(ft)), 0.0), 1.0)
  END IF
  f = f * f_season
  IF (l_lim) THEN
    IF (nsc_turg_acc(l,ft) < -1.0e29) nsc_turg_acc(l,ft) = 0.0
    nsc_turg_acc(l,ft) = nsc_turg_acc(l,ft) + f * timestep / period
  END IF
  ! Running annual mean of the growth factor (leaf renewal, l_nsc_wmax_norm).
  IF (l_ximpair_rec_leaf_turn .OR. l_nsc_wmax_norm) THEN
    IF (nsc_ft_mean(l,ft) < 0.0) nsc_ft_mean(l,ft) = 1.0
    nsc_ft_mean(l,ft) = nsc_ft_mean(l,ft) + a_yr * (f - nsc_ft_mean(l,ft))
  END IF
  ! Leaf renewal: production at the turnover rate, made while growing.
  IF (l_ximpair_rec_leaf_turn) THEN
    d_renew = MIN(g_leaf_0(ft) * f / MAX(nsc_ft_mean(l,ft), 0.05)             &
                  * timestep / (360.0 * 86400.0), 1.0)
    ximpair_kcap_leaf(l,ft) = 1.0 - (1.0 - ximpair_kcap_leaf(l,ft))          &
                              * (1.0 - d_renew)
  END IF
END DO

END SUBROUTINE triffid_nsc_turgor

! ---------------------------------------------------------------------
! Drought leaf shedding (l_phen_psi): called from sf_stom every time step for
! PFT ft. Sums the root-zone psi over the day; at the last time step of the
! day updates the consecutive dry / wet day counts and the leaf biomass
! factor phen_fbleaf (Xu et al. 2016, Eqn S7).
! ---------------------------------------------------------------------
SUBROUTINE phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,      &
                            psi_rz )

USE model_time_mod, ONLY: current_time
USE pftparm, ONLY: phen_psi_tlp, phen_dfbleaf, phen_ndays, phen_fbleaf_min

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

IF (phen_psi_tlp(ft) < -1.0e29) RETURN
CALL triffid_nsc_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)
  phen_psi_acc(l,ft) = phen_psi_acc(l,ft) + psi_rz(l) * timestep
END DO

! The day's update at its last time step.
IF (REAL(current_time%time) + timestep < 86400.0 - 0.5) RETURN

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

END SUBROUTINE phen_psi_update

! ---------------------------------------------------------------------
! The store without TRIFFID (l_triffid off), called from sf_stom every time
! step for PFT ft after NPP: sums the day's NPP, and at the last time step of
! the day runs the store step for ft with a one-day step, with leaf, root
! and wood carbon from the prescribed height (TRIFFID allometry).
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_daily( land_pts, veg_pts, veg_index, ft, timestep,    &
                              ht, npp )

USE jules_surface_types_mod, ONLY: npft
USE model_time_mod, ONLY: current_time
USE calc_c_comps_triffid_mod, ONLY: calc_c_comps_triffid

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
REAL(KIND=real_jlslsm) :: c_veg, lai_b
REAL(KIND=real_jlslsm), ALLOCATABLE :: woodc(:,:), leafc(:,:), rootc(:,:),   &
                                       npp_d(:,:), repair(:,:), ht_p(:,:),    &
                                       laib_p(:,:)

REAL(KIND=real_jlslsm), PARAMETER :: sec_per_day = 86400.0
REAL(KIND=real_jlslsm), PARAMETER :: days_per_trif_year = 360.0

CALL triffid_nsc_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)
  nsc_npp_day(l,ft) = nsc_npp_day(l,ft) + npp(l) * timestep
END DO

! Store step at the last time step of the day.
IF (REAL(current_time%time) + timestep < sec_per_day - 0.5) RETURN

ALLOCATE(woodc(land_pts,npft), leafc(land_pts,npft), rootc(land_pts,npft), &
         npp_d(land_pts,npft), repair(land_pts,npft), ht_p(land_pts,npft),   &
         laib_p(land_pts,npft))
woodc(:,:) = 0.0
ht_p(:,:) = 0.0
laib_p(:,:) = 0.0
leafc(:,:) = 0.0
rootc(:,:) = 0.0
npp_d(:,:) = 0.0
DO j = 1,veg_pts
  l = veg_index(j)
  CALL calc_c_comps_triffid(ft, ht(l), lai_b, leafc(l,ft), rootc(l,ft),    &
                            woodc(l,ft), c_veg)
  ht_p(l,ft)   = ht(l)
  laib_p(l,ft) = lai_b
  ! The day's NPP (kg C m-2) as a rate per 360 days (TRIFFID units).
  npp_d(l,ft) = nsc_npp_day(l,ft) * days_per_trif_year
  nsc_npp_day(l,ft) = 0.0
END DO

CALL triffid_nsc_store( land_pts, veg_pts, veg_index, days_per_trif_year,   &
                        woodc, leafc, rootc, ht_p, laib_p, npp_d, repair,    &
                        pft_only = ft )

DEALLOCATE(woodc, leafc, rootc, npp_d, repair, ht_p, laib_p)

END SUBROUTINE triffid_nsc_daily

END MODULE triffid_nsc_mod
