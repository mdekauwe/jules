! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE triffid_nsc_mod
!-----------------------------------------------------------------------------
! Non-structural carbon store for TRIFFID (l_triffid_nsc), after Trugman et
! al. (2018, Ecol. Lett. 21: 1552-1560, Eqns 1a-c and 5, Table S2), and the
! stem recovery, sapwood and heartwood that depend on it.
!
! 1. The store
!
! NPP goes into a store of non-structural carbon (sugars and starch), C
! (kg C m-2 of the PFT's area), and TRIFFID grows from the carbon the store
! releases, W, instead of from the NPP of the step:
!
!   dC/dt = NPP - W,
!   W     = W_max C_a / (k_c C_max + C_a),
!
! with NPP and W in kg C m-2 per 360 days. The release saturates as the
! store fills (a Michaelis-Menten form), and
!   C_max is the store's capacity: nsc_flush times the leaf plus root
!         carbon at the balanced LAI (enough for nsc_flush flushes of
!         leaves and fine roots), or nsc_cmax if that is set;
!   C_a   is the accessible store, C - nsc_reserve_frac C_max. The reserve
!         below it is never released for growth;
!   k_c   is nsc_kc (no units): the release is half of W_max when
!         C_a = k_c C_max;
!   W_max is nsc_wmax_fac times the running mean of NPP (e-folding time
!         nsc_tau_load years), the growth the plant is adapted to.
!         Trugman et al. set it by tree size, which TRIFFID does not have.
!
! The store starts full. Carbon above C_max goes to growth, less the share
! nsc_overflow_resp, which is respired. When NPP is negative and the
! accessible store is used up, the deficit passes to TRIFFID as negative
! growth and the vegetation shrinks: this is carbon starvation. With
! nsc_g_starv > 0 and competition on (l_veg_compete), starvation also
! removes area in lotka, at the rate nsc_g_starv (per 360 days) when the
! accessible store is empty, falling linearly to 0 at C_a = k_c C_max.
!
! 2. Growth limits (triffid_nsc_turgor)
!
! Growth can stop while photosynthesis goes on. The release is scaled by a
! growth factor f (0-1), the mean over the store step of
!
!   f = f_psi f_T f_season.
!
! f_psi is the turgor limit, used when nsc_psi_turg_zero is set. Cells need
! turgor to expand, so in drought growth stops before photosynthesis does:
!
!   f_psi = (psi_rz - psi_0) / (psi_1 - psi_0),   kept within 0-1,
!
! with psi_rz the root-zone water potential (MPa, close to predawn), psi_0
! = nsc_psi_turg_zero, where growth stops, and psi_1 = nsc_psi_turg_full,
! where it is unlimited. f_T is the cold limit, used when nsc_t_grow_zero
! is set:
!
!   f_T = (T_d - T_0) / (T_1 - T_0),   kept within 0-1,
!
! with T_d a running mean of the air temperature over about a day (degC),
! T_0 = nsc_t_grow_zero and T_1 = nsc_t_grow_full. f_season is 1 from the
! start of the growing season to nsc_grow_end_doy and 0 outside it. The
! season starts at nsc_flush_start_doy when a leaf flush is set (then
! nsc_grow_start_doy is not used here), and at nsc_grow_start_doy
! otherwise. With any limit on, the release is scaled so that a
! full store with f = 1 still releases W_max:
!
!   W = W_max f [C_a / (k_c C_max + C_a)] / [C_m / (k_c C_max + C_m)],
!
! with C_m the accessible store when the store is full. With
! l_nsc_wmax_norm, W_max is also divided by the running annual mean of f,
! so the store releases faster while growth is allowed and the yearly
! release still matches the mean NPP.
!
! 3. Growth threshold (nsc_grow_thresh set)
!
! The release for growth draws only on the store above nsc_grow_thresh
! C_max, which takes the place of the reserve in C_a. After a deep
! drawdown the store then refills before growth and ring renewal start
! again. The store down to the reserve still pays the deficits and the
! leaves.
!
! 4. Leaves
!
! With a leaf flush (nsc_flush_start_doy set), the year's new leaf carbon
! is drawn from the store over nsc_flush_days, before stem growth starts.
! Drought leaf shedding (l_phen_psi, Xu et al. 2016; phen_psi_update) sheds
! the canopy when the root zone stays dry and regrows it when it is wet
! again, and the store pays for the regrowth. With leaf renewal
! (l_ximpair_rec_leaf_turn), new leaves arrive with intact xylem, so the
! leaf's loss of conductance falls by the share of the canopy that is new.
!
! 5. Stem recovery by ring renewal (l_ximpair_rec_nsc)
!
! Used with the memory impairment model and dynamic TRIFFID
! (triffid_nsc_ring_renewal). Drought damage to the stem xylem is lasting,
! and the stem recovers as new wood with intact xylem is laid down. The new
! ring is TRIFFID's gross wood production from the release,
!
!   G = dC_W/dt + g_wood C_W   (kg C m-2 per 360 days),
!
! the change in the wood carbon C_W plus the wood turnover it replaces
! (g_wood, per 360 days). Its share of the stem conductance comes off the
! stem's loss of conductance L (0-1),
!
!   rho = nsc_xylem_cond_ratio G dt / X,    L -> MAX(L - rho, 0),
!
! with dt the TRIFFID step (360-day years) and X the conducting sapwood
! (kg C m-2). X is the sapwood pool (section 6) or, without it, TRIFFID's
! live stemwood eta_sl h LAI_bal (section 8), the pipe model, in which the
! sapwood area is proportional to the leaf area it supplies.
! No carbon is drawn from the store for the ring: the store
! sets its pace through W. The ring replaces damaged xylem first
! (l_nsc_ring_additive); with that off the ring dilutes the
! loss, L -> L (1 - rho).
!
! While the release is shut TRIFFID still removes the wood turnover, so the
! wood shrinks (G < 0) and regrows when the release opens again. A deficit
! D (kg C m-2) records the wood lost, and the regrowth repays D before any
! new ring counts. S is the ring year's net wood production so far. The
! ring year ends at the TRIFFID step that holds nsc_grow_end_doy or,
! without it, the last day of the year, and then
!
!   D <- MAX(D + MIN(S, 0), 0),    S <- 0,
!
! so over the years the wood counted as new rings adds up to the sum of
! MAX(S, 0).
!
! 6. Sapwood pool (l_nsc_sapwood_pool, on by default)
!
! The conducting sapwood X_s (kg C m-2) inside TRIFFID's wood is its own
! pool, as in MEDFATE (De Caceres et al. 2023):
!
!   dX_s/dt = G - (k_h + g_t) X_s,
!
! with k_h = nsc_heartwood_rate the share of the sapwood that becomes
! heartwood (wood that no longer conducts) each year, and g_t the wood
! turnover rate TRIFFID applies (its wood litter over its wood), which
! takes the same share of the sapwood as of the wood. Where the dump has
! no pool, the pool starts at its steady state for the current wood C_W,
!
!   X_s = C_W g_t / (k_h + g_t),
!
! so it does not drift through the run. The new ring renews about
! k_h + g_t of the sapwood each year, so the sapwood turnover, which can be
! measured, sets the recovery rate. A ring year with a net loss of wood
! (S < 0) also takes |S| from X_s.
!
! 7. Ring cohorts (l_nsc_sapwood_rings, with the sapwood pool)
!
! The sapwood pool is followed ring by ring. Each ring's share of the
! sapwood and its loss of conductance, and the growth year of the youngest
! ring, are kept in dumps (ximpair_ring_frac, ximpair_ring_loss,
! ximpair_ring_year). There are at most nsc_n_ring rings, and the oldest
! two merge when a new one starts. At each TRIFFID step TRIFFID's wood
! turnover takes from all rings alike, the new wood G dt joins the youngest
! ring with intact xylem, and heartwood forms from the oldest rings first.
! At the end of the ring year a new ring starts. A drought raises each
! conducting ring's loss to at least the loss it causes. The stem's loss
! is the mean of the rings' losses, weighted by their share of the sapwood
! and by their conductivity, which falls with ring age a (in winters
! passed, so a ring stays new until its calendar year ends):
!
!   c(a) = c_old + (1 - c_old) EXP(-a / e),
!
! with c_old = nsc_ring_c_old and e = nsc_ring_cond_efold (years). With
! these unset, all rings conduct alike. In a ring-porous stem
! (nsc_ring_func_doy set), the youngest ring's early wood conducts only
! from that day of the year; before it the flow goes through the older
! rings, with last year's early wood at c_old after the winter. The
! youngest ring holds only the wood grown so far, so the early wood is not
! paid for from the store ahead of growth. With the rings, recovery comes
! from the rings themselves: new rings dilute the stem's loss, and a
! drought's damage leaves with the rings it hit when they become
! heartwood. Nothing is queued, and l_nsc_ring_additive is ignored. The
! rings start, at the first update of the main run or from a dump without
! them, in the steady state of these dynamics, each with the stem's
! current loss: a ring grown m years ago holds a share in proportion to
! EXP(-g_t m), as many rings as hold G / (k_h + g_t) of new wood, and the
! youngest ring holds the part of this year's ring grown so far.
!
! 8. Pipe-model target (l_nsc_pipe_target, with the sapwood pool)
!
! The conducting sapwood has a target size from the pipe model, on the
! leaf area it supplied, last ring year's peak LAI L_pk (the balanced LAI
! until the first ring year has ended), and on the canopy height h (m).
! With the pipe heartwood rule or the Trugman response (section 10), L_pk
! is the peak of the leaf area the climate allows, the balanced LAI less
! drought leaf shedding (l_phen_psi), so leaves shed for damage do not
! lower the target:
!
!   X_t = eta_sl h L_pk   (kg C m-2),
!
! TRIFFID's live stemwood (live_stemwood) on the peak LAI. eta_sl is the
! sapwood carbon per unit leaf area and height,
!
!   eta_sl = F rho cmass / LA:SA,
!
! with LA:SA the leaf area per unit sapwood area, F the form factor of the
! wood (its volume over basal area x height), rho the wood density and
! cmass the carbon fraction of dry mass, so it is set from data (with a_ws
! scaled to keep a_ws eta_sl, which TRIFFID's growth depends on). X_t is
! at most the wood C_W, so a_ws should be at least 1. This cap guards
! against inputs that do not fit together: a warning names the point and
! the eta_sl that would fit the first time it binds, and nsc_xt_cap_count
! (an output) counts the TRIFFID steps at which it bound.
!
! nsc_heartwood_model sets how the sapwood turns to heartwood. With the
! pipe rule (1, default) the pipe model is kept, as in LPJ (Sitch et al.
! 2003): the sapwood above the target becomes heartwood at once, oldest
! rings first, and below it heartwood forms at k_h. Over a TRIFFID step,
! with X' the sapwood after TRIFFID's turnover,
!
!   H = MAX(X' (1 - EXP(-k_h dt)), X' - X_t).
!
! TRIFFID's allocation is fixed, so a sapwood below its target is not
! grown faster. With the fixed rate (0) heartwood forms at k_h whatever
! the target. The heartwood rate over k_h is an output
! (nsc_heartwood_factor). The pool starts at the target. The sapwood is
! at most the wood, and the heartwood is the rest,
!
!   H = C_W - X_s   (nsc_heartwood, an output).
!
! TRIFFID's wood carbon is not changed.
!
! 9. Conductance from the sapwood area (l_nsc_k_sapwood, with the
!    pipe-model target; nsc_k_area)
!
! The plant's maximum hydraulic conductance, kmax, is scaled by the
! conducting sapwood area against its target,
!
!   f_A = MIN(1, SUM_i X_i c_i / (X_t c_ref)),
!
! with X_i the sapwood carbon in ring i, c_i its conductivity by age and
! season, and c_ref the mean conductivity of a stand in equilibrium on the
! same day (the ring structure the rings start with, above, the youngest
! grown as far as the season has gone), so a stand at its target has
! f_A = 1 on any date.
! Without the rings f_A = MIN(1, X_s / X_t). f_A counts the sapwood's
! capacity only: drought damage in the rings lowers the conductance
! through the impairment memory. The ratio before the MIN is an output
! (nsc_k_area_raw).
!
! With the Trugman response (section 10) the stem conducts through its
! sapwood area, whatever the leaf area. The leaf area the sapwood can
! supply (nsc_lai_sapwood),
!
!   L_A = X_A / (eta_sl h),
!   X_A = X_s SUM_i f_i c_i / c_ref,
!
! sets the stem's conductance per unit ground area, kmax L_A, so its
! conductance per unit leaf area is kmax L_A / LAI and rises when leaves
! are shed. At the target and a full canopy L_A = LAI.
!
! 10. Response to damage (nsc_damage_response)
!
! With the Trugman response (1; needs the rings, the pipe-model target,
! l_phenol, fixed cover, l_nsc_k_sapwood and kplant_model = 1), a tree
! with lasting stem damage matches its leaves to its working xylem and
! rebuilds the xylem first, as in the optimal allocation of Trugman et
! al. (2018):
!   - The canopy loses the share of its leaves that matches the share
!     of the sapwood's conductance lost: the phenological LAI is at most
!     the balanced LAI times
!       SUM_i f_i c_i (1 - L_i) / SUM_i f_i c_i,
!     so an undamaged stem never caps the canopy, whatever its sapwood
!     (the heartwood rule sets that). Trugman et al. tie leaves to the
!     working xylem, and their allocation keeps the two in ratio; TRIFFID's
!     is fixed, so the cap is taken against the tree's undamaged state.
!     phenol sheds the leaves above the cap as leaf litter, at the rate of
!     phenological leaf drop, and the canopy regrows by the usual
!     phenology as the xylem recovers. The cap (LAI) is an output
!     (nsc_lai_sap_cap).
!   - The stem conducts through its whole sapwood (section 9), so the
!     supply per remaining leaf rises when leaves are shed.
!   - Repairing damaged xylem has first call on the store's release for
!     growth (as u, the share of the release to xylem, in Trugman et al.
!     eqn 5), so it obeys the growth limits as growth does: at each store
!     step up to that release is claimed, at most the damaged sapwood
!     that repair can restore, SUM_i X_i MAX(L_i - L_p, 0), with L_p the
!     highest loss the water potential has imposed since the last
!     TRIFFID step (damage at or below it would come back at once while
!     the drought lasts). The claim restores damaged xylem in the rings
!     in place, in proportion to their repairable damaged carbon. No ring
!     changes size, so a ring conducts at most as well as it did
!     undamaged. The claim is taken before TRIFFID splits its carbon
!     between growth and spreading (lambda_g), so a step's claim comes
!     out of growth by lambda_g and out of the spreading carbon (litter
!     with fixed cover) by the rest, and the balanced LAI then settles
!     lower. The claim goes to litter as the damaged wood it replaces, so
!     the sapwood and TRIFFID's wood carbon do not change. Outputs: the
!     claim (nsc_repair_claim) and the repair made (nsc_xylem_repair),
!     which is less when the step's heartwood or new wood has already
!     removed damaged wood; the rest of the claim still goes to litter.
! The target is set on the leaf area the climate allows (section 8), so
! the leaves shed for damage do not shrink the sapwood as well.
!
! 11. Without TRIFFID (l_triffid off)
!
! The store runs once a day (triffid_nsc_daily), with the leaf, root and
! wood carbon from the prescribed canopy height. The release leaves the
! store, as there is no growth to feed, and there is no stem recovery.
!
! MDK, 04/11/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: nsc_xylem_prod_pft, nsc_xylem_renew_pft
PUBLIC :: nsc_wood_lost
PUBLIC :: nsc_net_wood_growth
PUBLIC :: nsc_sapwood
PUBLIC :: nsc_heartwood, nsc_sapwood_target, nsc_lai_peak, nsc_lai_peak_prev
PUBLIC :: nsc_k_area_frac, nsc_k_area
PUBLIC :: nsc_k_area_raw, nsc_lai_sapwood, nsc_lai_sap_cap
PUBLIC :: nsc_xylem_renew_cw, nsc_heartwood_factor, nsc_xt_cap_count
PUBLIC :: nsc_xylem_repair_pft, nsc_repair_claim
PUBLIC :: nsc_f_a_min

! Smallest sapwood-area factor on kmax (and on the stem with
! the Trugman response), so the plant keeps a conducting path.
REAL(KIND=real_jlslsm), PARAMETER :: nsc_f_a_min = 0.01
PUBLIC :: nsc_pft, nsc_npp_mean, nsc_npp_wt, nsc_load_pft,                     &
          nsc_dc_pft, nsc_starv_pft, nsc_npp_day, nsc_turg_acc, nsc_resp_pft,  &
          nsc_t_run, nsc_ft_mean, phen_psi_acc, phen_nlow, phen_nhigh,         &
          phen_fbleaf, nsc_regrow_acc, phen_psi_update,                        &
          triffid_nsc_alloc, triffid_nsc_store, triffid_nsc_area,              &
          triffid_nsc_daily, triffid_nsc_turgor, triffid_nsc_ring_renewal
PUBLIC :: nsc_cmax_pft, nsc_wpot_pft, nsc_tseas_acc

! Store (kg C m-2 of the PFT); < 0 = not yet set (filled to C_max at the
! first TRIFFID step). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_pft(:,:)
! The store's capacity C_max at the last store step (kg C m-2 of the PFT),
! for the sink feedback on photosynthesis (l_nsc_sink_feedback); < 0 = not
! yet set (no feedback until the first store step). In dumps with the sink
! feedback.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_cmax_pft(:,:)
! The release for growth at the last store step with the turgor limit
! open (kg C m-2 of the PFT per 360 days), for the sink feedback: what the
! store would give growth if the turgor limit allowed, with the cold and
! season limits still applied (0 outside the growing season). In dumps
! with the sink feedback.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_wpot_pft(:,:)
! The growth factor without its turgor part (the cold and season factors)
! summed over the store period so far, as a fraction of the period, for
! nsc_wpot_pft (< -1e29 = not yet accumulating: 1 at the next store step).
! In dumps with the sink feedback and the growth limits.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_tseas_acc(:,:)
! Running mean of NPP (kg C m-2 per 360 days) and its weight, for W_max
! (bias-corrected mean = nsc_npp_mean / nsc_npp_wt). In dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_mean(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_npp_wt(:,:)
! Diagnostics of the last TRIFFID step (kg C m-2 per 360 days): the carbon
! released from the store (W, which TRIFFID grows from) and the net flux
! into the store (NPP - W).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_load_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_dc_pft(:,:)
! Overflow respiration of the last store step (kg C m-2 per 360 days). In
! dumps.
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
! Pipe-model target (l_nsc_pipe_target): the heartwood H = C_W - X_s of
! the last TRIFFID step (kg C m-2 of the PFT; < 0 = not yet set), a
! diagnostic; the target sapwood X_t of the last TRIFFID step (kg C m-2;
! 0 = not yet set), in dumps; and whether a point has been warned that
! its target was capped at the wood.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_heartwood(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_sapwood_target(:,:)
LOGICAL, ALLOCATABLE, SAVE :: nsc_xt_warned(:,:)
! Pipe-model target: the ring year's peak LAI so far and last ring year's
! peak LAI (< 0 = not yet set), in dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_lai_peak(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_lai_peak_prev(:,:)
! Conductance from the sapwood area (l_nsc_k_sapwood): the factor f_A on
! the plant's kmax of the last step (1 = full area). A diagnostic, set by
! sf_stom.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_k_area_frac(:,:)
! The sapwood-area ratio before the MIN of f_A of the last step (can
! exceed 1; 1 until set). A diagnostic, set by sf_stom.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_k_area_raw(:,:)
! Canopy limited by the sapwood (Trugman response): the largest LAI
! the functional sapwood can supply at the last phenology step (< 0 = no
! cap). A diagnostic, set by phenol.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_lai_sap_cap(:,:)
! Ring cohorts: the new ring's share of the stem's conductance weight,
! the share of the stem's loss it dilutes, of the last TRIFFID step (per
! 360 days). A diagnostic.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_renew_cw(:,:)
! Xylem repair (Trugman response): damaged sapwood replaced at the last
! TRIFFID step from the release claimed for it (kg C m-2 of the PFT per
! 360 days). A diagnostic.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xylem_repair_pft(:,:)
! Trugman response: the release claimed for xylem repair at the last
! store step (kg C m-2 of the PFT per 360 days), taken from the release
! for growth before TRIFFID grows; triffid adds it to the vegetation's
! input and, as the damaged wood it replaces, to the litter. Output
! nsc_repair_claim.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_repair_claim(:,:)
! Pipe-model target: the heartwood factor g of the last TRIFFID step, and
! the number of TRIFFID steps since the start (or restart) at which the
! target was capped at the wood. Diagnostics.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_heartwood_factor(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: nsc_xt_cap_count(:,:)
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
! Allocates the store's state and diagnostics for all land points and
! PFTs, if not already allocated, at their values for a fresh start.
! Several start below 0, which marks them as not yet set; each is filled
! at its first use.
!
! MDK, 01/11/26
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
IF (.NOT. ALLOCATED(nsc_heartwood)) THEN
  ALLOCATE(nsc_heartwood(n_land_pts, npft))
  nsc_heartwood(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_sapwood_target)) THEN
  ALLOCATE(nsc_sapwood_target(n_land_pts, npft))
  nsc_sapwood_target(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_lai_peak)) THEN
  ALLOCATE(nsc_lai_peak(n_land_pts, npft))
  ALLOCATE(nsc_lai_peak_prev(n_land_pts, npft))
  nsc_lai_peak(:,:)      = 0.0
  nsc_lai_peak_prev(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_xt_warned)) THEN
  ALLOCATE(nsc_xt_warned(n_land_pts, npft))
  nsc_xt_warned(:,:) = .FALSE.
END IF
IF (.NOT. ALLOCATED(nsc_k_area_frac)) THEN
  ALLOCATE(nsc_k_area_frac(n_land_pts, npft))
  nsc_k_area_frac(:,:) = 1.0
END IF
IF (.NOT. ALLOCATED(nsc_k_area_raw)) THEN
  ALLOCATE(nsc_k_area_raw(n_land_pts, npft))
  nsc_k_area_raw(:,:) = 1.0
END IF
IF (.NOT. ALLOCATED(nsc_lai_sap_cap)) THEN
  ALLOCATE(nsc_lai_sap_cap(n_land_pts, npft))
  nsc_lai_sap_cap(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(nsc_xylem_renew_cw)) THEN
  ALLOCATE(nsc_xylem_renew_cw(n_land_pts, npft))
  nsc_xylem_renew_cw(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_xylem_repair_pft)) THEN
  ALLOCATE(nsc_xylem_repair_pft(n_land_pts, npft))
  nsc_xylem_repair_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_repair_claim)) THEN
  ALLOCATE(nsc_repair_claim(n_land_pts, npft))
  nsc_repair_claim(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(nsc_heartwood_factor)) THEN
  ALLOCATE(nsc_heartwood_factor(n_land_pts, npft))
  ALLOCATE(nsc_xt_cap_count(n_land_pts, npft))
  nsc_heartwood_factor(:,:) = 1.0
  nsc_xt_cap_count(:,:)     = 0.0
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
IF (.NOT. ALLOCATED(nsc_tseas_acc)) THEN
  ALLOCATE(nsc_tseas_acc(n_land_pts, npft))
  nsc_tseas_acc(:,:) = -1.0e30
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
! Pass the PFT's carbon through its non-structural carbon store, so that
! growth uses stored carbon and does not follow each step's NPP directly.
! It runs one store step for each natural PFT (or only pft_only), at each
! TRIFFID step or, without TRIFFID, once a day. On entry npp is the NPP of
! the step; on return it is W, the carbon released from the store, which
! TRIFFID grows from (both kg C m-2 per 360 days). In order:
!
!   1. Sets the capacity C_max (the store starts full) and updates the
!      running mean of NPP that sets W_max.
!   2. Pays, from the accessible store, the leaf regrowth queued after
!      drought shedding and then this step's share of the leaf flush.
!      Regrowth in the flush window counts towards the flush, so the
!      leaves are not paid for twice.
!   3. Works out the release for growth on the store above its floor (the
!      reserve or the growth threshold), scaled by the mean growth factor
!      when growth limits are on. The release with only the turgor limit
!      open (the cold and season limits applied) is kept in nsc_wpot_pft
!      for the sink feedback on photosynthesis.
!   4. Takes W - NPP from the store, never below the reserve: a deficit
!      beyond the accessible store passes to TRIFFID as negative growth.
!      With a growth threshold, the release for growth is first capped so
!      it cannot take the store below the threshold.
!   5. Sends carbon above C_max to growth, less the share
!      nsc_overflow_resp, which is respired, and sets the starvation
!      disturbance for lotka.
!   6. With the Trugman response, gives repairing damaged xylem first
!      call on the release (nsc_repair_claim); TRIFFID grows from the rest.
!
! MDK, 04/11/26
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
USE model_time_mod, ONLY: is_spinup
USE jules_vegetation_mod, ONLY: nsc_damage_response, nsc_damage_trugman
USE xylem_impairment_memory_mod, ONLY: ximpair_ring_damaged

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
  ! Turgor-limited sink: turgor factor (0-1), and the cold and season
  ! factor (0-1) of the release with the turgor limit open.
  f_turg, f_ts,                                                                &
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
      ! The release with the turgor limit open, for the sink feedback: the
      ! period's mean cold and season factor f_ts in place of the mean
      ! growth factor, so it is 0 outside the growing season.
      f_ts = 1.0
      IF (is_set(nsc_tseas_acc(l,n)))                                          &
        f_ts = MIN(MAX(nsc_tseas_acc(l,n), 0.0), 1.0)
      nsc_tseas_acc(l,n) = 0.0
      nsc_wpot_pft(l,n) = 0.0
      IF (c_grow > 0.0 .AND. c_grow_max > 0.0)                                 &
        nsc_wpot_pft(l,n) = w_max * f_ts * c_grow                              &
                            / (nsc_kc(n) * c_max + c_grow)                     &
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

    ! Trugman response: repairing damaged xylem has first call on the
    ! release for growth, so it obeys the growth limits as growth does
    ! (the leaf regrowth, the flush and the overflow below are not used).
    ! The claim is at most that release and the damaged sapwood repair can
    ! restore, X_s SUM_i f_i MAX(L_i - L_p, 0) (ximpair_ring_damaged).
    ! Not in spin-up or before the rings and the sapwood are set.
    nsc_repair_claim(l,n) = 0.0
    IF (nsc_damage_response == nsc_damage_trugman .AND. .NOT. is_spinup      &
        .AND. .NOT. PRESENT(pft_only) .AND. w > 0.0) THEN
      IF (ALLOCATED(nsc_sapwood)) THEN
        IF (nsc_sapwood(l,n) > 0.0) THEN
          nsc_repair_claim(l,n) = MIN(w, nsc_sapwood(l,n)                     &
                                         * ximpair_ring_damaged(l, n) / dt)
          w = w - nsc_repair_claim(l,n)
        END IF
      END IF
    END IF

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
! Adjusts the store after TRIFFID changes a PFT's area (lotka). The store
! is per unit area of the PFT. Area lost takes its store with it, to
! litter (lit_store, kg C m-2 of land per 360 days). Area gained shares
! the existing store, so the store per unit area is multiplied by
! frac_old / frac, the PFT's fraction before the step over its fraction
! after. The wood lost and not yet regrown (nsc_wood_lost) and the ring
! year's net wood growth (nsc_net_wood_growth) are diluted the same way
! when area is gained, as the new cover has lost no wood and grown none
! this ring year.
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
! Slow the store's release of carbon for growth when the plant is too dry
! or too cold to grow, or out of its growing season. Called from sf_stom
! every time step for PFT ft, it works out the step's growth factor
!
!   f = f_psi f_T f_season   (0-1; see the module header)
!
! and adds f times the step's share of the store period (the TRIFFID
! period, or one day without TRIFFID) to nsc_turg_acc, so the store step
! uses the mean of f over its period.
!
! f_psi: growth needs turgor, and in drought it stops before
! photosynthesis does (Lempereur et al. 2015: stem growth stops at a
! predawn water potential near -1.1 MPa while assimilation goes on).
! f_T: wood formation needs daily mean temperatures above about 5-8 degC
! (Rossi et al. 2008), so f_T uses a running mean of the air temperature
! with an e-folding time of one day, and a sunny winter midday does not
! start growth. f_season: trees can stay dormant through a mild winter, so
! the store refills from winter NPP. With a leaf flush the season starts
! on the flush start day.
!
! The routine also keeps <f>, the running annual mean of f (e-folding 360
! days), used by leaf renewal and l_nsc_wmax_norm. With leaf renewal
! (l_ximpair_rec_leaf_turn), new leaves arrive with intact xylem, so the
! leaf's remaining loss of conductance is multiplied each step by
! (1 - dR), with dR the share of the canopy that is new,
!
!   dR = g_leaf_0 f / <f> dt,
!
! g_leaf_0 the leaf turnover rate (per 360 days) and dt the time step
! (360-day years). Leaves are made only while growth is allowed, and a
! year's production still matches g_leaf_0. Leaf renewal does not run in
! spin-up, where the stem damage is held fixed.
!
! MDK, 01/11/26
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
                          doy_start, f_ts
    ! Store period (s); growth factor (0-1); weights of this time step in
    ! the daily and annual running means; leaf renewal this time step;
    ! growing-season factor (0 or 1); day of year now and of the start of
    ! the release; the growth factor without its turgor part.
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
  f_ts = 1.0
  IF (is_set(nsc_t_grow_zero(ft))) THEN
    IF (nsc_t_run(l,ft) < 0.0) nsc_t_run(l,ft) = tair(l)
    nsc_t_run(l,ft) = nsc_t_run(l,ft) + a_t * (tair(l) - nsc_t_run(l,ft))
    f_ts = MIN(MAX((nsc_t_run(l,ft) - tm - nsc_t_grow_zero(ft))                &
                   / (nsc_t_grow_full(ft) - nsc_t_grow_zero(ft)), 0.0), 1.0)
    f = f * f_ts
  END IF
  f = f * f_season
  f_ts = f_ts * f_season

  ! The store step uses the mean over its period, of f and of the factor
  ! without turgor (for the sink feedback's release with turgor open).
  IF (l_lim) THEN
    IF (.NOT. is_set(nsc_turg_acc(l,ft))) nsc_turg_acc(l,ft) = 0.0
    nsc_turg_acc(l,ft) = nsc_turg_acc(l,ft) + f * timestep / period
    IF (.NOT. is_set(nsc_tseas_acc(l,ft))) nsc_tseas_acc(l,ft) = 0.0
    nsc_tseas_acc(l,ft) = nsc_tseas_acc(l,ft) + f_ts * timestep / period
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
! Shed leaves when the soil stays dry (l_phen_psi; Xu et al. 2016, New
! Phytol, Notes S1 Eqn S7), so a droughted plant loses water through less
! leaf area. It is called from sf_stom every time step for PFT ft, and does
! nothing unless phen_psi_tlp is set.
!
! Each step adds the root-zone water potential to a daily sum. At the last
! time step of the day the daily mean (MPa, a stand-in for predawn) is
! compared with the turgor loss point phen_psi_tlp (MPa, negative). A day
! below phen_psi_tlp is dry, and a day above 0.5 phen_psi_tlp is wet; each
! count of days in a row goes back to 0 when the run breaks. After
! phen_ndays dry days in a row, the leaf biomass factor phen_fbleaf (0-1,
! the share of the full canopy allowed) falls by phen_dfbleaf a day; after
! phen_ndays wet days in a row it rises by phen_dfbleaf a day. It stays
! between phen_fbleaf_min and 1. phenol caps the canopy at this factor and
! queues the leaf carbon of the regrowth (nsc_regrow_acc), which the store
! pays for.
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
! Run the non-structural carbon store when TRIFFID is off (l_triffid off).
! It is called from sf_stom every time step for PFT ft after NPP, adds the
! step's NPP to a daily sum, and at the last time step of the day runs
! the store step (triffid_nsc_store) for ft with a one-day step. The leaf
! and root carbon come from the prescribed canopy height through
! TRIFFID's allometry. The released
! carbon leaves the store, as there is no growth to feed. Stem recovery
! needs TRIFFID's wood growth, so it does not run, and the starvation term
! is only a diagnostic, as the PFT fractions are fixed.
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
! Return the factor f_A (0-1) that scales the plant's kmax at land point
! l for PFT n with how much conducting sapwood it has (l_nsc_k_sapwood),
!
!   f_A = MIN(1, X_s SUM_i f_i c_i / (X_t c_ref)),
!
! with X_s the sapwood and X_t its pipe-model target (both kg C m-2), f_i
! each ring's share of the sapwood, c_i its conductivity by age and season
! (ximpair_ring_cweight gives the sum), and c_ref the mean conductivity of
! a stand in equilibrium on the same day (ximpair_ring_cref). Without the
! rings f_A = MIN(1, X_s / X_t). f_A is 1 until the sapwood and its target
! are set (in spin-up and before the first ring step), and at least
! nsc_f_a_min, so the plant keeps a conducting path. f_raw, if present, is
! the ratio before the MIN and the floor.
!
! MDK, 01/11/26
! ---------------------------------------------------------------------
FUNCTION nsc_k_area( l, n, f_raw ) RESULT( f_a )

USE jules_vegetation_mod, ONLY: l_nsc_sapwood_rings
USE xylem_impairment_memory_mod, ONLY: ximpair_ring_frac, ximpair_ring_year,  &
                                       ximpair_ring_cweight, ximpair_ring_cref

! Land point and PFT.
INTEGER, INTENT(IN) :: l, n
! The ratio before the MIN and the floor.
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: f_raw
! Area factor on kmax (0-1).
REAL(KIND=real_jlslsm) :: f_a

! Sapwood and its target (kg C m-2), the rings' equilibrium mean
! conductivity and the ratio.
REAL(KIND=real_jlslsm) :: x_s, x_t, c_ref, r

f_a = 1.0
IF (PRESENT(f_raw)) f_raw = 1.0
IF (.NOT. ALLOCATED(nsc_sapwood) .OR. .NOT. ALLOCATED(nsc_sapwood_target))     &
  RETURN
x_s = nsc_sapwood(l,n)
x_t = nsc_sapwood_target(l,n)
IF (x_s < 0.0 .OR. x_t <= 0.0) RETURN

r = x_s / x_t
IF (l_nsc_sapwood_rings .AND. ALLOCATED(ximpair_ring_frac)) THEN
  IF (ximpair_ring_frac(l,n,1) >= 0.0 .AND.                                   &
      ximpair_ring_year(l,n) >= 0.0) THEN
    c_ref = ximpair_ring_cref(l, n)
    r = 1.0
    IF (c_ref > 0.0) THEN
      r   = x_s * ximpair_ring_cweight(l, n) / (x_t * c_ref)
      f_a = MIN(1.0, r)
    END IF
    f_a = MAX(f_a, nsc_f_a_min)
    IF (PRESENT(f_raw)) f_raw = r
    RETURN
  END IF
END IF
f_a = MAX(MIN(1.0, x_s / x_t), nsc_f_a_min)
IF (PRESENT(f_raw)) f_raw = r

END FUNCTION nsc_k_area

! ---------------------------------------------------------------------
! Return the leaf area index the sapwood of PFT n at land point l can
! supply under the pipe model (the Trugman response to damage),
!
!   L_A = X_A / (eta_sl h),
!
! the inverse of the pipe-model target (live_stemwood), with h the canopy
! height (m) and X_A = X_s SUM_i f_i c_i / c_ref (kg C m-2) the sapwood
! weighted by how well its rings conduct, as in f_A (nsc_k_area) before
! the MIN; without the rings X_A = X_s. Drought damage is not taken off.
! A stand at its target with a full canopy has L_A = the LAI the target
! was set on. Returns -1 when the sapwood, its target or the pipe model is
! not set (in spin-up and before the first ring step).
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
FUNCTION nsc_lai_sapwood( l, n, ht ) RESULT( lai_s )

USE jules_vegetation_mod, ONLY: l_nsc_sapwood_rings
USE pftparm, ONLY: eta_sl
USE xylem_impairment_memory_mod, ONLY: ximpair_ring_frac, ximpair_ring_year,  &
                                       ximpair_ring_cweight, ximpair_ring_cref

! Land point and PFT.
INTEGER, INTENT(IN) :: l, n
! Canopy height (m).
REAL(KIND=real_jlslsm), INTENT(IN) :: ht
! Leaf area index the sapwood can supply (m2 m-2; -1 = not set).
REAL(KIND=real_jlslsm) :: lai_s

! The sapwood that counts (kg C m-2) and the rings' equilibrium mean
! conductivity.
REAL(KIND=real_jlslsm) :: x_a, c_ref

lai_s = -1.0
IF (.NOT. ALLOCATED(nsc_sapwood) .OR. .NOT. ALLOCATED(nsc_sapwood_target))     &
  RETURN
IF (nsc_sapwood(l,n) < 0.0 .OR. nsc_sapwood_target(l,n) <= 0.0) RETURN
IF (eta_sl(n) <= 0.0 .OR. ht <= 0.0) RETURN

x_a = nsc_sapwood(l,n)
IF (l_nsc_sapwood_rings .AND. ALLOCATED(ximpair_ring_frac)) THEN
  IF (ximpair_ring_frac(l,n,1) < 0.0 .OR. ximpair_ring_year(l,n) < 0.0)      &
    RETURN
  c_ref = ximpair_ring_cref(l, n)
  IF (c_ref <= 0.0) RETURN
  x_a = x_a * ximpair_ring_cweight(l, n) / c_ref
END IF

lai_s = x_a / (eta_sl(n) * ht)

END FUNCTION nsc_lai_sapwood

! ---------------------------------------------------------------------
! Let the stem recover from drought damage as new sapwood rings take over
! from damaged ones (l_ximpair_rec_nsc). It is called by triffid after
! vegcarb at each TRIFFID step. It does not run in spin-up, where
! the impairment is held fixed, and it moves no carbon. For each natural
! PFT that uses the memory impairment model:
!
!   1. The new ring is TRIFFID's gross wood production of the step,
!
!        G = dwood / dt + g_wood C_W,old   (kg C m-2 per 360 days),
!
!      the wood increment plus the wood turnover it replaced (wood_litc),
!      with dt the TRIFFID step (360-day years) and C_W,old the wood
!      before the step. G counts only after the deficit D of wood lost
!      and not yet regrown is repaid (see the module header).
!   2. Without the ring cohorts, the ring's share of the stem conductance,
!
!        rho = nsc_xylem_cond_ratio G dt / X,
!
!      with X the conducting sapwood before the step's growth, is queued
!      in ximpair_repair_pend. The next impairment update takes it off
!      the stem's loss (l_nsc_ring_additive) or dilutes the loss by it.
!   3. With the sapwood pool, X is the pool X_s, which gains G and loses
!      heartwood and its share of TRIFFID's wood turnover.
!   4. With the ring cohorts nothing is queued: the new wood joins the
!      youngest ring, the turnover takes from all rings alike and
!      heartwood takes the oldest rings first.
!   5. With the pipe-model target and the pipe heartwood rule, the
!      sapwood above the target (on last ring year's peak leaf area)
!      becomes heartwood, and the heartwood is the rest of the wood. The
!      target is at most the wood; the steps at which this cap binds are
!      counted.
!   6. With the Trugman response, the release the store step claimed
!      for repair replaces damaged sapwood in the rings
!      (ximpair_ring_repair).
!   7. At the end of the ring year, the year's net wood growth S settles
!      the deficit D, a net loss of wood also comes off the sapwood, and
!      with the pipe-model target this ring year's peak LAI becomes last
!      year's.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_nsc_ring_renewal( land_pts, trif_pts, trif_index,           &
                                     r_gamma, ht, lai_bal, lai, woodc, dwood,  &
                                     wood_litc )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE jules_vegetation_mod, ONLY: l_ximpair_rec_nsc, xylem_impairment_memory
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_pool
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_rings, l_nsc_pipe_target
USE jules_vegetation_mod, ONLY: nsc_heartwood_model, nsc_heartwood_pipe,      &
                                nsc_damage_response, nsc_damage_trugman,       &
                                l_phen_psi
USE pftparm, ONLY: pft_xylem_impairment_model, nsc_xylem_cond_ratio
USE pftparm, ONLY: eta_sl, phen_psi_tlp
USE pftparm, ONLY: is_set, nsc_grow_end_doy, nsc_heartwood_rate
USE calc_c_comps_triffid_mod, ONLY: live_stemwood
USE model_time_mod, ONLY: is_spinup
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap
USE xylem_impairment_memory_mod, ONLY: ximpair_repair_pend,                    &
                                       ximpair_memory_alloc
USE xylem_impairment_memory_mod, ONLY: ximpair_lock, ximpair_area_w,           &
                                       ximpair_ring_frac, ximpair_ring_init,   &
                                       ximpair_ring_grow, ximpair_ring_repair, &
                                       ximpair_ring_next_year,                 &
                                       ximpair_ring_peak
USE logging_mod, ONLY: log_warn
USE string_utils_mod, ONLY: to_string

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
  ! Current (phenological) LAI.
  lai(land_pts,npft),                                                          &
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
! The calendar year in which a ring started at this step grows.
INTEGER :: g_year
! This step ends the ring year (it holds the end of the growing season
! or, without one, the last day of the year).
LOGICAL :: l_ring_end
! The leaf area the sapwood serves and that the target is set on, the
! pipe-model target sapwood (kg C m-2), the heartwood rate (per 360
! days), the sapwood after TRIFFID's turnover and the heartwood of the
! step (kg C m-2).
REAL(KIND=real_jlslsm) :: lai_ref, lai_pot, x_t, rate_h, x_turn, h_step
! The target is set on the leaf area the climate allows (pipe heartwood
! rule or Trugman response).
LOGICAL :: l_pot
! The new ring's share of the conductance weight.
REAL(KIND=real_jlslsm) :: w_cw
! Xylem repair: damaged sapwood replaced at this step (kg C m-2).
REAL(KIND=real_jlslsm) :: r_rep

! Share of the sapwood the pipe heartwood rule always leaves after a
! step, so a target far below the sapwood cannot empty it.
REAL(KIND=real_jlslsm), PARAMETER :: h_keep = 1.0e-3

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

  l_pot = nsc_heartwood_model == nsc_heartwood_pipe .OR.                      &
          nsc_damage_response == nsc_damage_trugman

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
    ! The calendar year in which a ring started at this step grows.
    g_year = ximpair_ring_next_year(n)

    DO t = 1,trif_pts
      l = trif_index(t)
      nsc_xylem_prod_pft(l,n)  = 0.0
      nsc_xylem_renew_pft(l,n) = 0.0
      nsc_xylem_renew_cw(l,n)  = 0.0
      wood_old = woodc(l,n) - dwood(l,n)

      ! The leaf area the sapwood serves: the balanced LAI, or with the
      ! pipe-model target last ring year's peak LAI (the balanced LAI until
      ! the first ring year has ended). This ring year's peak so far: of
      ! the LAI, or with the pipe heartwood rule or the Trugman response of
      ! the leaf area the climate allows, the balanced LAI less drought
      ! leaf shedding (l_phen_psi), so leaves shed for damage do not lower
      ! the target.
      lai_ref = lai_bal(l,n)
      IF (l_nsc_pipe_target) THEN
        lai_pot = lai(l,n)
        IF (l_pot) THEN
          lai_pot = lai_bal(l,n)
          IF (l_phen_psi .AND. is_set(phen_psi_tlp(n)))                       &
            lai_pot = lai_pot * MIN(MAX(phen_fbleaf(l,n), 0.01), 1.0)
        END IF
        nsc_lai_peak(l,n) = MAX(nsc_lai_peak(l,n), lai_pot)
        IF (nsc_lai_peak_prev(l,n) >= 0.0) lai_ref = nsc_lai_peak_prev(l,n)
      END IF

      ! The conducting sapwood X from the allometry, or with the sapwood
      ! pool the prognostic X_s, which starts (where it is not yet set) at
      ! its steady state for the wood before the step,
      !
      !   X_s = C_W g_t / (k_h + g_t),   g_t = wood_litc / C_W,
      !
      ! or with the pipe-model target at the target. TRIFFID's live
      ! stemwood on lai_ref is the pipe-model target X_t, at most the wood.
      x_s = live_stemwood(n, ht(l,n), lai_ref)
      x_t = x_s
      IF (l_nsc_pipe_target .AND. x_t > wood_old) THEN
        ! The cap is a guard: eta_sl (from LA:SA, form factor and wood
        ! density) and the height give more sapwood than the wood, so a_ws
        ! is below 1. The first time, say by how much and what eta_sl
        ! would fit; count every step.
        IF (.NOT. nsc_xt_warned(l,n))                                          &
          CALL log_warn(RoutineName,                                           &
                        "Pipe-model target capped at the wood carbon " //      &
                        "at land point " //                                    &
                        TRIM(to_string(l)) // ", PFT " //                      &
                        TRIM(to_string(n)) // ": X_t / C_W = " //              &
                        TRIM(to_string(x_t / MAX(wood_old, 1.0e-12))) //       &
                        "; eta_sl_io at most " //                              &
                        TRIM(to_string(eta_sl(n) * MAX(wood_old, 0.0)          &
                                       / x_t)) // " would fit, with " //       &
                        "a_ws_io raised by the same factor.")
        nsc_xt_warned(l,n) = .TRUE.
        nsc_xt_cap_count(l,n) = nsc_xt_cap_count(l,n) + 1.0
        x_t = MAX(wood_old, 0.0)
      END IF
      IF (l_nsc_sapwood_pool) THEN
        g_turn   = 0.0
        IF (wood_old > 0.0) g_turn = wood_litc(l,n) / wood_old
        IF (nsc_sapwood(l,n) < 0.0) THEN
          IF (l_nsc_pipe_target) THEN
            nsc_sapwood(l,n) = x_t
          ELSE
            nsc_sapwood(l,n) = MAX(wood_old, 0.0) * g_turn                     &
                               / (nsc_heartwood_rate(n) + g_turn)
          END IF
        END IF
        x_s = nsc_sapwood(l,n)
      END IF

      ! Heartwood. With the pipe heartwood rule (and the pipe-model target)
      ! the sapwood above the target becomes heartwood at once, and below
      ! it heartwood forms at k_h. Over the step, with X' the sapwood after
      ! TRIFFID's turnover,
      !
      !   H = MAX(X' (1 - EXP(-k_h dt)), X' - X_t),
      !
      ! taken oldest rings first; rate_h is the rate that removes H over
      ! the step, and the heartwood factor (an output) is rate_h / k_h.
      ! Otherwise heartwood forms at k_h.
      rate_h = nsc_heartwood_rate(n)
      IF (l_nsc_pipe_target) THEN
        nsc_sapwood_target(l,n)   = x_t
        nsc_heartwood_factor(l,n) = 1.0
        IF (nsc_heartwood_model == nsc_heartwood_pipe) THEN
          x_turn = x_s * EXP(-g_turn / r_gamma)
          IF (x_turn > x_min) THEN
            h_step = MAX(x_turn * (1.0 - EXP(-rate_h / r_gamma)),              &
                         x_turn - x_t)
            h_step = MIN(h_step, x_turn * (1.0 - h_keep))
            rate_h = -LOG(1.0 - h_step / x_turn) * r_gamma
            IF (nsc_heartwood_rate(n) > 0.0)                                   &
              nsc_heartwood_factor(l,n) = rate_h / nsc_heartwood_rate(n)
          END IF
        END IF
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
      ! its count and state go on. With the ring cohorts the conductance
      ! ratio does not act: the renewal is the new wood's share of the
      ! sapwood, and nothing is queued.
      IF (x_s > x_min) THEN
        IF (l_nsc_sapwood_rings) THEN
          nsc_xylem_renew_pft(l,n) = g_wood_prod / x_s
        ELSE
          nsc_xylem_renew_pft(l,n) = nsc_xylem_cond_ratio(n) * g_wood_prod    &
                                     / x_s
          ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n)                  &
                                     + nsc_xylem_renew_pft(l,n) / r_gamma
        END IF
      END IF

      ! Ring cohorts: TRIFFID's turnover takes all rings alike, the new
      ! wood joins the youngest ring and heartwood takes the oldest (as the
      ! pool below); the stem's loss follows at the next update. Rings not
      ! yet set start from the stem's loss before the area dilution still
      ! pending.
      IF (l_nsc_sapwood_rings) THEN
        IF (ximpair_ring_frac(l,n,1) < 0.0)                                    &
          CALL ximpair_ring_init(l, n,                                         &
                                 MIN(ximpair_lock(l,n)                         &
                                     / MAX(ximpair_area_w(l,n), 1.0e-6), 1.0))
        CALL ximpair_ring_grow(l, n, x_s * EXP(-g_turn / r_gamma),             &
                               g_wood_prod / r_gamma, rate_h / r_gamma,        &
                               l_ring_end, g_year, w_cw)
        ! The new ring's share of the conductance weight, per 360 days:
        ! the share of the stem's loss it dilutes.
        nsc_xylem_renew_cw(l,n) = w_cw * r_gamma
      END IF

      ! Sapwood pool: the oldest sapwood becomes heartwood, TRIFFID's
      ! wood turnover takes its share and the new ring is added,
      !
      !   dX_s/dt = G - (k_h + g_t) X_s,   k_h = nsc_heartwood_rate,
      !
      ! so in equilibrium X_s = G / (k_h + g_t) (with the pipe-model
      ! target k_h is g k_h).
      IF (l_nsc_sapwood_pool)                                                  &
        nsc_sapwood(l,n) = nsc_sapwood(l,n)                                    &
                           * EXP(-(rate_h + g_turn) / r_gamma)                 &
                           + g_wood_prod / r_gamma

      ! Xylem repair (the Trugman response): the release the store step
      ! claimed for repair replaces damaged sapwood above the period's
      ! peak imposed loss in place, so neither the sapwood nor TRIFFID's
      ! wood changes. The peak then starts again for the next period.
      nsc_xylem_repair_pft(l,n) = 0.0
      IF (nsc_damage_response == nsc_damage_trugman .AND.                      &
          l_nsc_sapwood_rings) THEN
        CALL ximpair_ring_repair(l, n, nsc_sapwood(l,n),                       &
                                 nsc_repair_claim(l,n) / r_gamma, r_rep)
        nsc_xylem_repair_pft(l,n) = r_rep * r_gamma
      END IF
      IF (l_nsc_sapwood_rings) ximpair_ring_peak(l,n) = 0.0

      ! End of the ring year: its net production S settles the count,
      !
      !   D <- MAX(D + MIN(S, 0), 0),   S <- 0,
      !
      ! so D keeps only the over-count (growth counted before a loss that
      ! was not regrown) for next year's growth to repay, and a lasting
      ! loss of wood does not hold back later rings. With the sapwood
      ! pool a net loss of wood also leaves the sapwood, |S| (all rings
      ! alike). With the pipe-model target the ring year's peak LAI
      ! becomes last year's.
      IF (l_ring_end) THEN
        nsc_wood_lost(l,n) = MAX(nsc_wood_lost(l,n)                            &
                                    + MIN(nsc_net_wood_growth(l,n), 0.0), 0.0)
        IF (l_nsc_sapwood_pool)                                                &
          nsc_sapwood(l,n) = MAX(nsc_sapwood(l,n)                              &
                                 + MIN(nsc_net_wood_growth(l,n), 0.0), 0.0)
        nsc_net_wood_growth(l,n) = 0.0
        IF (l_nsc_pipe_target) THEN
          nsc_lai_peak_prev(l,n) = nsc_lai_peak(l,n)
          nsc_lai_peak(l,n)      = 0.0
        END IF
      END IF

      ! Pipe-model target: the sapwood is at most the wood (all rings
      ! alike), and the heartwood is the rest,
      !
      !   H = C_W - X_s.
      IF (l_nsc_pipe_target) THEN
        nsc_sapwood(l,n)   = MIN(nsc_sapwood(l,n), MAX(woodc(l,n), 0.0))
        nsc_heartwood(l,n) = MAX(woodc(l,n), 0.0) - nsc_sapwood(l,n)
      END IF

    END DO
  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_nsc_ring_renewal


END MODULE triffid_nsc_mod
