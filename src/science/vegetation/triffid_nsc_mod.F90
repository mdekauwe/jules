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
! goes to growth. Growth respiration is not charged at release: the JULES
! NPP the store receives is already net of it.
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
          nsc_dc_pft, nsc_starv_pft, nsc_npp_day, nsc_turg_acc,               &
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
                              pft_only )

USE jules_surface_types_mod, ONLY: npft
USE jules_vegetation_mod, ONLY: l_ximpair_rec_nsc, xylem_impairment_memory,  &
                                l_triffid
USE pftparm, ONLY: g_leaf_0, pft_xylem_impairment_model,                     &
                   nsc_flush, nsc_kc, nsc_reserve_frac, nsc_tau_load,        &
                   nsc_wmax_fac, nsc_sapwood_frac, nsc_hold_plc,             &
                   nsc_g_starv, ximpair_psi_growth, nsc_repair_window,       &
                   nsc_plc_background, nsc_la_sa, nsc_wood_density,          &
                   nsc_psi_turg_zero
USE jules_surface_mod, ONLY: cmass
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
  f_turg, c_acc_max
    ! Repair: loss of conductance, conducting xylem (kg C m-2), leaf and fine
    ! root of the intact plant (kg C m-2), turnover TRIFFID charges against
    ! growth (per 360 days), share of W, growth gate, repairable loss.
    ! Turgor-limited sink: turgor factor, accessible store when full (kg C
    ! m-2).
LOGICAL :: l_rep, l_mem, l_turg

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
n2 = npft
IF (PRESENT(pft_only)) THEN
  n1 = pft_only
  n2 = pft_only
END IF

DO n = n1,n2

  l_mem = pft_xylem_impairment_model(n) == xylem_impairment_memory
  l_rep = l_ximpair_rec_nsc .AND. l_mem
  IF (l_rep) CALL ximpair_memory_alloc( land_pts )

  a_mean = MIN(dt / MAX(nsc_tau_load(n), TINY(1.0_real_jlslsm)), 1.0)
  l_turg = nsc_psi_turg_zero(n) > -1.0e29

  DO t = 1,trif_pts
    l = trif_index(t)

    c_max = nsc_flush(n) * (leafc(l,n) + rootc(l,n))
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

    ! Michaelis-Menten loading on the accessible store.
    c_acc = MAX(nsc_pft(l,n) - c_res, 0.0)
    w = 0.0
    IF (l_turg) THEN
      ! Turgor-limited sink, normalised so the full store releases W_max.
      f_turg = 1.0
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

    ! Overflow above the capacity goes to growth.
    IF (nsc_pft(l,n) > c_max) THEN
      w = w + (nsc_pft(l,n) - c_max) / dt
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
    IF (l_rep .AND. w > 0.0) THEN
      IF (nsc_la_sa(n) > 0.0) THEN
        ! Pipe model on the balanced leaf area (sapwood area x height).
        x_s = lai_bal(l,n) / nsc_la_sa(n) * ht(l,n) * nsc_wood_density(n)     &
              * cmass
      ELSE
        x_s = nsc_sapwood_frac(n) * woodc(l,n)
      END IF
      IF (plc_rep > 0.0 .AND. x_s > x_min) THEN
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
        u = x_s / (x_s + l_tot) * MAX(1.0 - turn_l / w, 0.0)
        repair(l,n) = MIN(f_grow * u * w, plc_rep * x_s / dt)
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
                               psi_rz )

USE jules_vegetation_mod, ONLY: l_triffid, triffid_period
USE pftparm, ONLY: nsc_psi_turg_full, nsc_psi_turg_zero

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
REAL(KIND=real_jlslsm) :: period, f

IF (nsc_psi_turg_zero(ft) < -1.0e29) RETURN
CALL triffid_nsc_alloc( land_pts )

period = 86400.0
IF (l_triffid) period = 86400.0 * REAL(triffid_period)

DO j = 1,veg_pts
  l = veg_index(j)
  f = (1.0e-6 * psi_rz(l) - nsc_psi_turg_zero(ft))                         &
      / (nsc_psi_turg_full(ft) - nsc_psi_turg_zero(ft))
  f = MIN(MAX(f, 0.0), 1.0)
  IF (nsc_turg_acc(l,ft) < -1.0e29) nsc_turg_acc(l,ft) = 0.0
  nsc_turg_acc(l,ft) = nsc_turg_acc(l,ft) + f * timestep / period
END DO

END SUBROUTINE triffid_nsc_turgor

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
