! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Leaf shedding: the canopy caps of the leaf turnover options.
!
! Each option caps the phenological state (LAI over the balanced LAI) that
! phenol allows; leaves above the cap are shed at the phenological leaf-drop
! rate and go to litter, and the canopy regrows when the cap rises. The
! options, each off by default:
!
!   l_phen_psi           drought: the predawn (root-zone) water potential
!                        below the turgor loss point for days in a row (Xu
!                        et al. 2016)
!   l_phen_fuse          leaves as hydraulic fuses: the canopy falls with
!                        the plant's daily water status (Wolfe et al. 2016)
!   l_phen_heat          leaf heat damage, the thermal load of Arnold et al.
!                        (2025)
!   nsc_damage_response  1 (Trugman): the canopy follows the share of the
!                        sapwood's conductance that still works (Trugman
!                        et al. 2018)
!
! The regrowth of leaves shed by the drought and heat caps is paid from the
! TRIFFID store (phen_shed_debt holds the share still owed); the regrowth
! under the damage cap is free, as the canopy recovers with the wood.
! leaf_shed_combine sets how the heat cap joins the others. TRIFFID's own
! drought and cold turnover (dgl_dm, fsmc_of, dgl_dt, tleaf_of in leaf_lit)
! acts on the turnover rate instead and stays there.
!
! The daily updates are called from sf_stom every time step; phenol calls
! leaf_shed_caps. The states are in dumps.
!
! How the pieces fit together:
!
!   sf_stom (every step)  -> phen_psi_update     updates phen_fbleaf  (S2)
!                         -> phen_fuse_update    updates shed_f_psi   (S3)
!                         -> heat_damage_update  updates heat_dam_pft (S5)
!   phenol (TRIFFID step) -> leaf_shed_caps      turns the states into a
!                                                cap on phen, sheds above
!                                                it and books the debt
!   triffid_nsc_mod       -> leaf_shed_climate_cap  the climate's cap for
!                                                the pipe-model target
!
! Every cap is a share of the full canopy (0-1, 1 = no cap). Each is kept
! at 0.01 or more, the same floor as phenol's leaf drop, so the canopy
! never falls to zero and can always regrow.
!
! Changes from gce d7805aa that change the answers of existing runs (all
! other runs are unchanged):
!
!   - l_phen_psi: the predawn value is the root-zone psi at the day's last
!     step (was the daily mean), frozen days change neither count, and the
!     store pays only to regrow the leaves this cap shed (was all regrowth
!     while l_phen_psi was on).
!   - nsc_damage_response = 1 with budburst (phen_bb_on): a canopy shed to
!     the 0.01 floor by the damage cap no longer goes dormant; only the
!     cold (or TRIFFID drought) leaf drop does that.
!   - l_ximpair_leaf_loss and ximpair_leaf_sens are retired (the Trugman
!     response replaces them); old namelists still read, but
!     l_ximpair_leaf_loss = .TRUE. stops the run.
!
! References
!
!   Arnold, P. A. et al. (2025) A framework for modelling thermal load
!     sensitivity across life. Global Change Biology 31: e70315.
!     doi:10.1111/gcb.70315
!   Schoolfield, R. M., Sharpe, P. J. H. and Magnuson, C. E. (1981)
!     Non-linear regression of biological temperature-dependent rate
!     models based on absolute reaction-rate theory. Journal of
!     Theoretical Biology 88: 719-731.
!   Trugman, A. T. et al. (2018) Tree carbon allocation explains forest
!     drought-kill and recovery patterns. Ecology Letters 21: 1552-1560.
!   Wolfe, B. T., Sperry, J. S. and Kursar, T. A. (2016) Does leaf shedding
!     protect stems from cavitation during seasonal droughts? A test of
!     the hydraulic fuse hypothesis. New Phytologist 212: 1007-1018.
!   Xu, X., Medvigy, D., Powers, J. S., Becknell, J. M. and Guan, K. (2016)
!     Diversity in plant hydraulic traits explains seasonal and
!     inter-annual variations of vegetation dynamics in seasonally dry
!     tropical forests. New Phytologist 212: 80-95.
!
! MDK, 09/11/26

MODULE leaf_shed_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

! Only the routines phenol, sf_stom and the store need, and the states
! the dumps and the output need, are public; shed_paid_caps, repair_rate
! and root_zone_frozen are helpers used only here.
PRIVATE
PUBLIC :: leaf_shed_alloc, leaf_shed_caps, leaf_shed_climate_cap,             &
          leaf_shed_heat_clear, phen_psi_update, phen_fuse_update,            &
          heat_damage_update
PUBLIC :: phen_nlow, phen_nhigh, phen_fbleaf, shed_f_psi, shed_nup,           &
          shed_psi_pd, shed_psi_md, shed_psi_plant, heat_dam_pft,             &
          phen_shed_debt, phen_max_pft, phen_regrow_c

! All the states are (land points, PFTs), allocated by leaf_shed_alloc on
! first use (phenol calls it every time, the updates and the dump reading
! when needed). With every option off they stay at their starting values.

! Drought (l_phen_psi), in dumps: consecutive dry and wet days, and the
! leaf biomass factor (0-1, the share of the full canopy allowed).
! The day counts are REAL so they go through the dumps like the other
! states.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_nlow(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_nhigh(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_fbleaf(:,:)

! Fuse (l_phen_fuse), in dumps: the share of the full canopy allowed (0-1),
! days in a row on which the day's share was above it, and the day so far's
! highest root-zone psi and lowest leaf psi (Pa; -1e30 and 1e30 = none yet).
! The two running extremes are dumped so a run restarted mid-day carries
! on with the same day.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shed_f_psi(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shed_nup(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shed_psi_pd(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shed_psi_md(:,:)

! Fuse diagnostic: the last day's plant water status (MPa; 0 = not yet).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: shed_psi_plant(:,:)

! Heat (l_phen_heat), in dumps: the canopy's accumulated heat damage
! (fraction of the damage threshold, 0-1).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: heat_dam_pft(:,:)

! Share of the full canopy shed by the drought and heat caps and not yet
! regrown (0-1); its regrowth is paid from the store. In dumps.
! phenol adds to it when it sheds below the free (damage) cap and takes
! from it as the canopy regrows; it never exceeds 1 - phen.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_shed_debt(:,:)

! Diagnostics of the last phenology call: the cap all options together put
! on the phenological state (0-1; 1 = no cap), and the leaf carbon of the
! regrowth charged to the store (kg C m-2 per 360 days).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_max_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_regrow_c(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEAF_SHED_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Allocate the shedding states (full canopy, nothing owed, no damage) if
! they are not yet allocated.
!
! It is safe to call every time step: it returns at once after the first
! call. Values read from a dump overwrite these starting values.
! ---------------------------------------------------------------------
SUBROUTINE leaf_shed_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: n_land_pts

! All the states are allocated together, so one of them tells.
IF (ALLOCATED(phen_fbleaf)) RETURN

ALLOCATE(phen_nlow(n_land_pts, npft), phen_nhigh(n_land_pts, npft),           &
         phen_fbleaf(n_land_pts, npft))
ALLOCATE(shed_f_psi(n_land_pts, npft), shed_nup(n_land_pts, npft),            &
         shed_psi_pd(n_land_pts, npft), shed_psi_md(n_land_pts, npft),        &
         shed_psi_plant(n_land_pts, npft))
ALLOCATE(heat_dam_pft(n_land_pts, npft), phen_shed_debt(n_land_pts, npft))
ALLOCATE(phen_max_pft(n_land_pts, npft), phen_regrow_c(n_land_pts, npft))

! Start with no dry or wet days counted, the full canopy allowed (1), the
! day's psi extremes empty (the sentinels -1e30 and 1e30, so the first
! step of the day replaces them), no heat damage and nothing owed.
phen_nlow(:,:)      = 0.0
phen_nhigh(:,:)     = 0.0
phen_fbleaf(:,:)    = 1.0
shed_f_psi(:,:)     = 1.0
shed_nup(:,:)       = 0.0
shed_psi_pd(:,:)    = -1.0e30
shed_psi_md(:,:)    = 1.0e30
shed_psi_plant(:,:) = 0.0
heat_dam_pft(:,:)   = 0.0
phen_shed_debt(:,:) = 0.0
phen_max_pft(:,:)   = 1.0
phen_regrow_c(:,:)  = 0.0

END SUBROUTINE leaf_shed_alloc

! ---------------------------------------------------------------------
! The caps on the phenological state of PFT n, for phenol.
!
! phen_max_free is the cap whose regrowth is free (the damage cap,
! nsc_damage_response = 1); phen_max adds the caps whose regrowth the store
! pays for (drought and heat, leaf_shed_climate_cap):
!
!   leaf_shed_combine = 0:  phen_max = MIN(free, drought, heat)
!   leaf_shed_combine = 1:  phen_max = MIN(free, drought) heat
!
! Each cap is at least 0.01. lai_sap_cap is the canopy the damage cap
! allows (m2 m-2; < 0 = no cap), for output. No damage cap in spin-up or
! before the sapwood rings are set.
!
! phenol uses the two caps to split what it sheds: leaves shed down to
! phen_max_free were lost to damage and regrow free; leaves shed below it
! were lost to drought or heat and go to phen_shed_debt.
! ---------------------------------------------------------------------
SUBROUTINE leaf_shed_caps( land_pts, veg_pts, veg_index, n, lai_bal,          &
                           phen_max, phen_max_free, lai_sap_cap )

USE jules_vegetation_mod, ONLY: nsc_damage_response, nsc_damage_trugman
USE xylem_impairment_memory_mod, ONLY: ximpair_ring_working
USE model_time_mod, ONLY: is_spinup

INTEGER, INTENT(IN) ::                                                         &
  ! Total number of land points.
  land_pts,                                                                    &
  ! Number of vegetated points.
  veg_pts,                                                                     &
  ! Index of vegetated points on the land grid.
  veg_index(land_pts),                                                         &
  ! Plant functional type.
  n

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Balanced growth LAI.
  lai_bal(land_pts)

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Largest phenological state all the caps allow.
  phen_max(land_pts),                                                          &
  ! The same from the free cap only.
  phen_max_free(land_pts),                                                     &
  ! Canopy the damage cap allows (m2 m-2; < 0 = no cap).
  lai_sap_cap(land_pts)

INTEGER :: j, l

! Working share of the sapwood's conductance (< 0 = not set), and the
! drought and heat caps.
REAL(KIND=real_jlslsm) :: f_work, f_dry, f_heat

! No cap anywhere until an option sets one, including at land points that
! are not vegetated (not in veg_index).
phen_max(:)      = 1.0
phen_max_free(:) = 1.0
lai_sap_cap(:)   = -1.0

DO j = 1,veg_pts
  l = veg_index(j)

  ! The damage cap: the canopy follows the share of the sapwood's
  ! conductance that still works (the Trugman response; Trugman et al.
  ! 2018).
  ! In spin-up the xylem damage is held, so the cap is left off; before
  ! the first ring year ends ximpair_ring_working is < 0 and there is no
  ! cap either.
  IF (nsc_damage_response == nsc_damage_trugman .AND. .NOT. is_spinup) THEN
    f_work = ximpair_ring_working(l, n)
    IF (f_work >= 0.0) THEN
      phen_max_free(l) = MAX(f_work, 0.01)
      lai_sap_cap(l)   = f_work * lai_bal(l)
    END IF
  END IF

  ! The paid caps (drought, heat; 1 where off), joined with the damage
  ! cap.
  CALL shed_paid_caps( l, n, f_dry, f_heat )
  phen_max(l) = join_caps( phen_max_free(l), f_dry, f_heat )

  ! Keep the joint cap for output (allocated by leaf_shed_alloc, which
  ! phenol calls first).
  IF (ALLOCATED(phen_max_pft)) phen_max_pft(l,n) = phen_max(l)
END DO

END SUBROUTINE leaf_shed_caps

! ---------------------------------------------------------------------
! The share of the full canopy the climate allows at point l for PFT n:
! the drought and heat caps together (as in leaf_shed_caps), without the
! damage cap; 1 with no option on. The store uses it for the leaf area the
! sapwood serves, so leaves shed for damage do not lower that target.
!
! (Were the damage cap included, damage would shed leaves, the smaller
! canopy would lower the sapwood target, and the damage would then look
! smaller than it is.)
! ---------------------------------------------------------------------
FUNCTION leaf_shed_climate_cap( l, n ) RESULT( f_cap )

INTEGER, INTENT(IN) :: l, n

REAL(KIND=real_jlslsm) :: f_cap

REAL(KIND=real_jlslsm) :: f_dry, f_heat

! Joined the same way as in leaf_shed_caps, with the damage cap set to 1.
CALL shed_paid_caps( l, n, f_dry, f_heat )
f_cap = join_caps( 1.0_real_jlslsm, f_dry, f_heat )

END FUNCTION leaf_shed_climate_cap

! ---------------------------------------------------------------------
! Clear the heat damage carried by the leaves phenol has just shed (phen
! from phen_old to phen_new) at point l for PFT n. phen_max_free is the
! damage cap of this call.
!
! The heat cap 1 - shed_heat_sens dam asks the canopy to drop the
! damaged leaves. Once they have gone, their damage has gone with them:
! the leaves left (and those regrown later) are undamaged, so the damage
! falls by the share shed for heat over shed_heat_sens,
!
!   dam = MAX(dam - d_heat / shed_heat_sens, 0),
!
! and the cap rises again as the canopy is shed down to it. The share shed
! for heat, d_heat, is what was shed below the cap the other options
! alone would allow (join_caps with no heat cap), so shedding for drought
! or damage clears none of it. Without this, damage repairs so slowly
! once high (repair_rate) that one heatwave would cap the canopy for years.
! ---------------------------------------------------------------------
SUBROUTINE leaf_shed_heat_clear( l, n, phen_max_free, phen_old, phen_new )

USE jules_vegetation_mod, ONLY: l_phen_heat
USE pftparm, ONLY: is_set, ctmax_1h, shed_heat_sens

INTEGER, INTENT(IN) :: l, n

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! The damage cap of this call (0-1).
  phen_max_free,                                                               &
  ! Phenological state before and after the shedding.
  phen_old, phen_new

! Drought and heat caps, the cap without heat, and the share shed for heat.
REAL(KIND=real_jlslsm) :: f_dry, f_heat, f_noheat, d_heat

IF (.NOT. l_phen_heat .OR. .NOT. ALLOCATED(heat_dam_pft)) RETURN
IF (.NOT. is_set(ctmax_1h(n)) .OR. shed_heat_sens(n) <= 0.0) RETURN

CALL shed_paid_caps( l, n, f_dry, f_heat )
f_noheat = join_caps( phen_max_free, f_dry, 1.0_real_jlslsm )
d_heat   = MAX(MIN(phen_old, f_noheat) - phen_new, 0.0)
IF (d_heat > 0.0)                                                              &
  heat_dam_pft(l,n) = MAX(heat_dam_pft(l,n) - d_heat / shed_heat_sens(n), 0.0)

END SUBROUTINE leaf_shed_heat_clear

! ---------------------------------------------------------------------
! Join the damage cap (f_free) with the drought and heat caps, as set by
! leaf_shed_combine:
!
!   0: MIN(f_free, f_dry, f_heat)    the tightest limit wins;
!   1: MIN(f_free, f_dry) f_heat     heat damage is a loss of leaves on
!                                    top of whatever the water supply
!                                    allows, so the caps multiply.
!
! The one place the rule is written, so phenol's cap and the store's
! climate cap (f_free = 1) cannot fall out of step. At least 0.01.
! ---------------------------------------------------------------------
FUNCTION join_caps( f_free, f_dry, f_heat ) RESULT( f_cap )

USE jules_vegetation_mod, ONLY: leaf_shed_combine

REAL(KIND=real_jlslsm), INTENT(IN) :: f_free, f_dry, f_heat

REAL(KIND=real_jlslsm) :: f_cap

IF (leaf_shed_combine == 1) THEN
  f_cap = MAX(MIN(f_free, f_dry) * f_heat, 0.01)
ELSE
  f_cap = MIN(f_free, f_dry, f_heat)
END IF

END FUNCTION join_caps

! ---------------------------------------------------------------------
! The drought cap (the smaller of l_phen_psi and l_phen_fuse) and the heat
! cap (l_phen_heat, MAX(1 - shed_heat_sens damage, shed_heat_fmin)) at point
! l for PFT n, each at least 0.01; 1 where the option is off.
!
! An option counts only when both its switch is on and its key PFT
! parameter is set (phen_psi_tlp, shed_la50, ctmax_1h), so in a run with
! several PFTs it can act on some PFTs and not others.
! ---------------------------------------------------------------------
SUBROUTINE shed_paid_caps( l, n, f_dry, f_heat )

USE jules_vegetation_mod, ONLY: l_phen_psi, l_phen_fuse, l_phen_heat
USE pftparm, ONLY: is_set, phen_psi_tlp, shed_la50, ctmax_1h,                 &
                   shed_heat_sens, shed_heat_fmin

INTEGER, INTENT(IN) :: l, n

REAL(KIND=real_jlslsm), INTENT(OUT) :: f_dry, f_heat

f_dry  = 1.0
f_heat = 1.0
! The states are not there yet (no update has run): no caps.
IF (.NOT. ALLOCATED(phen_fbleaf)) RETURN

! Both drought options can be on at once; the tighter one wins.
IF (l_phen_psi .AND. is_set(phen_psi_tlp(n)))                                 &
  f_dry = MIN(f_dry, MAX(phen_fbleaf(l,n), 0.01))
IF (l_phen_fuse .AND. is_set(shed_la50(n)))                                   &
  f_dry = MIN(f_dry, MAX(shed_f_psi(l,n), 0.01))
! Heat: the canopy loses shed_heat_sens of its leaves per unit damage, but
! keeps at least shed_heat_fmin.
IF (l_phen_heat .AND. is_set(ctmax_1h(n)))                                    &
  f_heat = MAX(1.0 - shed_heat_sens(n) * heat_dam_pft(l,n),                   &
               shed_heat_fmin(n), 0.01)

END SUBROUTINE shed_paid_caps

! ---------------------------------------------------------------------
! Drought (l_phen_psi; Xu et al. 2016, Notes S1 Eqn S7), so a
! droughted plant loses water through less leaf area. Called from sf_stom
! every time step for PFT ft; acts at the day's last time step for PFTs
! with phen_psi_tlp set.
!
! The root-zone water potential (MPa), a soil value and so the predawn
! water potential, is compared with the turgor loss point phen_psi_tlp
! (MPa, negative). A day below phen_psi_tlp is dry, and a day above
! 0.5 phen_psi_tlp is wet; each count of days in a row goes back to 0 when
! the run breaks. After phen_ndays dry days in a row, the leaf biomass
! factor phen_fbleaf falls by phen_dfbleaf a day; after phen_ndays wet days
! in a row it rises by phen_dfbleaf a day. It stays between phen_fbleaf_min
! and 1.
!
! A day whose root zone is frozen at that step (frozen share of its soil
! water above ximpair_frozen_frac, as for the xylem impairment) changes
! neither count: JULES sets the soil water potential from the unfrozen
! water, so frozen soil would read as drought.
!
! Days between 0.5 phen_psi_tlp and phen_psi_tlp are neither dry nor wet:
! they break both runs, so the canopy holds where it is. Once a run has
! reached phen_ndays, every further day of it moves phen_fbleaf again, so
! a long drought keeps shedding down to phen_fbleaf_min.
! ---------------------------------------------------------------------
SUBROUTINE phen_psi_update( land_pts, veg_pts, veg_index, ft, timestep,        &
                            psi_rz )

USE model_time_mod, ONLY: current_time
USE pftparm, ONLY: is_set, phen_psi_tlp, phen_dfbleaf, phen_ndays,            &
                   phen_fbleaf_min

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

! The day's root-zone (predawn) water potential (MPa).
REAL(KIND=real_jlslsm) :: psi_d

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PHEN_PSI_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Only at the day's last time step, for PFTs with phen_psi_tlp set.
! current_time%time is the seconds into the day at the start of this
! step, so the step that ends at midnight passes (half a second spare for
! rounding). This works for any time step that divides the day.
IF (.NOT. is_set(phen_psi_tlp(ft)) .OR.                                        &
    REAL(current_time%time) + timestep < 86400.0 - 0.5) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL leaf_shed_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)
  ! A frozen root zone leaves both counts as they are.
  IF (root_zone_frozen(l)) CYCLE
  ! Pa to MPa, the units of phen_psi_tlp.
  psi_d = 1.0e-6 * psi_rz(l)

  ! Count the dry and wet days in a row; a day that is not dry (or not
  ! wet) starts that count again.
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

  ! Shed or regrow a step once a run is long enough. A day cannot be both
  ! dry and wet, so at most one of the two counts is above 0.
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
! Leaves as hydraulic fuses (l_phen_fuse; Wolfe et al. 2016). Called from
! sf_stom every time step for PFT ft; does nothing unless shed_la50 is set.
!
! Through the day it keeps the highest root-zone psi (the predawn value)
! and the lowest leaf psi (midday). At the day's last time step the plant
! water status is their mean,
!
!   psi_plant = 0.5 (psi_predawn + psi_midday)   (MPa),
!
! as Wolfe et al. took the mean of predawn and midday leaf psi. The share
! of the full canopy for the day falls linearly with it, at 0.5 at
! shed_la50 (their LA50):
!
!   f_day = 0.5 + (psi_plant - shed_la50) / shed_la_width,
!
! held between shed_fmin and 1. The allowed share shed_f_psi drops to f_day
! at once, and rises to it only after phen_ndays days in a row with f_day
! above it, so one wet day does not refill the canopy. Steps with a frozen
! root zone are left out (as in phen_psi_update), and a day with none left
! changes nothing. Leaf psi of 0 (a step without a leaf water potential)
! is not used; the midday value is then the predawn value.
!
! shed_la_width is the psi range (MPa) over which f_day goes from 0 to 1,
! so f_day is 1 at shed_la50 + 0.5 shed_la_width and 0 at
! shed_la50 - 0.5 shed_la_width. Shedding acts fast and regrowth slowly,
! as leaves are lost in days but take weeks to grow.
! ---------------------------------------------------------------------
SUBROUTINE phen_fuse_update( land_pts, veg_pts, veg_index, ft, timestep,       &
                             psi_rz, psi_leaf )

USE model_time_mod, ONLY: current_time
USE pftparm, ONLY: is_set, shed_la50, shed_la_width, shed_fmin, phen_ndays

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
  ! Leaf water potential (Pa).
  psi_leaf(land_pts)

INTEGER :: j, l

! Predawn and midday values and the plant water status (MPa).
REAL(KIND=real_jlslsm) :: psi_pd, psi_md, psi_p

! The day's share of the full canopy (0-1).
REAL(KIND=real_jlslsm) :: f_day

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PHEN_FUSE_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. is_set(shed_la50(ft))) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL leaf_shed_alloc( land_pts )

! Every step: keep the day's wettest root zone (predawn) and driest leaf
! (midday) so far, in Pa. Leaf psi is negative whenever it has been
! solved; 0 means no value this step.
DO j = 1,veg_pts
  l = veg_index(j)
  IF (.NOT. root_zone_frozen(l)) THEN
    shed_psi_pd(l,ft) = MAX(shed_psi_pd(l,ft), psi_rz(l))
    IF (psi_leaf(l) < 0.0) shed_psi_md(l,ft) = MIN(shed_psi_md(l,ft),         &
                                                   psi_leaf(l))
  END IF
END DO

! The day's update at its last time step (as in phen_psi_update).
IF (REAL(current_time%time) + timestep < 86400.0 - 0.5) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO j = 1,veg_pts
  l = veg_index(j)
  ! Still at its sentinel: every step of the day was frozen, so skip it.
  IF (shed_psi_pd(l,ft) > -1.0e29) THEN
    ! Pa to MPa. The midday value cannot be wetter than predawn; with no
    ! leaf psi all day it is taken as the predawn value.
    psi_pd = 1.0e-6 * shed_psi_pd(l,ft)
    psi_md = psi_pd
    IF (shed_psi_md(l,ft) < 1.0e29)                                            &
      psi_md = MIN(1.0e-6 * shed_psi_md(l,ft), psi_pd)
    psi_p  = 0.5 * (psi_pd + psi_md)
    shed_psi_plant(l,ft) = psi_p

    ! The day's canopy share, linear in psi_plant.
    f_day = 0.5 + (psi_p - shed_la50(ft)) / shed_la_width(ft)
    f_day = MIN(MAX(f_day, shed_fmin(ft)), 1.0)

    ! Drier than the allowed share: drop to it now and restart the count.
    ! Wetter: count the day, and rise to it after phen_ndays such days (to
    ! that day's value, not the highest of the run). The same: the run
    ! is broken.
    IF (f_day < shed_f_psi(l,ft)) THEN
      shed_f_psi(l,ft) = f_day
      shed_nup(l,ft)   = 0.0
    ELSE IF (f_day > shed_f_psi(l,ft)) THEN
      shed_nup(l,ft) = shed_nup(l,ft) + 1.0
      IF (shed_nup(l,ft) >= phen_ndays(ft)) shed_f_psi(l,ft) = f_day
    ELSE
      shed_nup(l,ft) = 0.0
    END IF
  END IF
  ! Empty the day's extremes for the next day.
  shed_psi_pd(l,ft) = -1.0e30
  shed_psi_md(l,ft) = 1.0e30
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE phen_fuse_update

! ---------------------------------------------------------------------
! Leaf heat damage (l_phen_heat) of PFT ft over one time step, after the
! thermal load sensitivity framework of Arnold et al. (2025). Called from
! sf_stom every time step; does nothing for a PFT without ctmax_1h set.
!
! Damage accrues above tcrit_heat at the reciprocal of the time to reach
! the damage threshold, from the thermal death time curve
!
!   log10(t_h) = (ctmax_1h - T) / z_tdt   (t_h in hours, T in degC),
!
! so t_h = 1 h at T = ctmax_1h and ten times longer for each z_tdt degrees
! cooler. Repair follows a Sharpe-Schoolfield Arrhenius function of
! temperature (repair_rate), slowed as the damage builds. The damage is
! then diluted by the leaf turnover g_leaf_0 (new leaves are undamaged)
! and, with shed_heat_tau set, decays with that e-folding time (days).
! The temperature is chosen by heat_phen_temp.
!
! The damage is a fraction of the threshold: 1 means the leaves have taken
! a full lethal dose. Each step adds timestep / t_h of it, so a leaf held
! at ctmax_1h reaches 1 in an hour.
! ---------------------------------------------------------------------
SUBROUTINE heat_damage_update( land_pts, veg_pts, veg_index, ft, timestep,     &
                               tstar, t_can, t_leaf_sun )

USE pftparm, ONLY: is_set, g_leaf_0, ctmax_1h, z_tdt, tcrit_heat, krep_ref,    &
                   krep_dec, krep_t_hi, shed_heat_tau
USE jules_vegetation_mod, ONLY: heat_phen_temp, heat_temp_tstar,               &
                                heat_temp_sun
USE conversions_mod, ONLY: zerodegc

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
  ! Surface temperature (K).
  tstar(land_pts),                                                             &
  ! Canopy leaf temperature (K).
  t_can(land_pts),                                                             &
  ! Sunlit leaf temperature (K).
  t_leaf_sun(land_pts)

INTEGER :: j, l

! Temperature driving the damage (K, degC).
REAL(KIND=real_jlslsm) :: t_dam, tc

! Damage added and repaired this step (fraction of the threshold).
REAL(KIND=real_jlslsm) :: d_dam, d_rep

! Fraction of the damage cleared by turnover and decay this step.
REAL(KIND=real_jlslsm) :: f_clear

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='HEAT_DAMAGE_UPDATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. is_set(ctmax_1h(ft))) THEN
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF
CALL leaf_shed_alloc( land_pts )

DO j = 1,veg_pts
  l = veg_index(j)

  ! The leaf temperature that drives the damage: the surface (tstar), the
  ! sunlit leaves (two-leaf canopy) or, by default, the canopy leaves.
  IF (heat_phen_temp == heat_temp_tstar) THEN
    t_dam = tstar(l)
  ELSE IF (heat_phen_temp == heat_temp_sun) THEN
    t_dam = t_leaf_sun(l)
  ELSE
    t_dam = t_can(l)
  END IF
  ! A leaf temperature below 150 K was never solved (a configuration
  ! without leaf temperatures leaves it at 0): use the surface
  ! temperature, as sf_stom does.
  IF (t_dam < 150.0) t_dam = tstar(l)
  tc = t_dam - zerodegc

  ! Damage this step: the step's length (h) over the time to the
  ! threshold at this temperature, t_h = 10**((ctmax_1h - T) / z_tdt).
  ! None below tcrit_heat.
  d_dam = 0.0
  IF (tc > tcrit_heat(ft))                                                     &
    d_dam = (timestep / 3600.0) / 10.0**((ctmax_1h(ft) - tc) / z_tdt(ft))

  ! No damage held and none added (most steps): nothing to repair or
  ! clear, so skip the repair rate's exponentials.
  IF (heat_dam_pft(l,ft) <= 0.0 .AND. d_dam <= 0.0) CYCLE

  ! Repair this step. repair_rate is in % of the threshold per minute:
  ! 0.01 turns % into a fraction, timestep / 60 the step into minutes.
  d_rep = 0.01 * (timestep / 60.0)                                             &
          * repair_rate(t_dam, krep_ref(ft), krep_dec(ft), krep_t_hi(ft),      &
                        heat_dam_pft(l,ft))

  heat_dam_pft(l,ft) = MIN(1.0, MAX(0.0, heat_dam_pft(l,ft) + d_dam - d_rep))

  ! Clear part of the damage: the share of the leaves replaced by the
  ! background turnover this step (g_leaf_0 is per 360 days), plus decay
  ! with e-folding time shed_heat_tau (days) when it is set.
  f_clear = g_leaf_0(ft) * timestep / (360.0 * 86400.0)
  IF (is_set(shed_heat_tau(ft)) .AND. shed_heat_tau(ft) > 0.0)                 &
    f_clear = f_clear + timestep / (shed_heat_tau(ft) * 86400.0)
  heat_dam_pft(l,ft) = heat_dam_pft(l,ft) * (1.0 - MIN(f_clear, 1.0))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE heat_damage_update

! ---------------------------------------------------------------------
! Repair rate (% of the damage threshold per minute) as a function of
! temperature: the Sharpe-Schoolfield Arrhenius model (Schoolfield et al.
! 1981) scaled by the repair rate coefficient k_ref, as in Arnold et al.
! (2025). The coefficient falls once damage builds (their Eqn S6, with the
! undamaged share in percent), and is 0 at full damage. Above t_hi repair
! is deactivated.
!
!   k_rep = k_ref' exp(T_A (1/T_ref - 1/T))
!           / (1 + exp(T_AL (1/T - 1/T_L)) + exp(T_AH (1/T_H - 1/T))),
!
! with k_ref' = k_ref k_ref^(k_dec / alive). The numerator is the
! Arrhenius rise with temperature (k_rep = k_ref' at T_ref with no
! deactivation); the two terms below switch repair off below T_L and above
! T_H. As k_ref < 1, k_ref^(k_dec / alive) falls towards 0 as alive falls.
! ---------------------------------------------------------------------
FUNCTION repair_rate( t_k, k_ref, k_dec, t_hi, damage ) RESULT( k_rep )

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Leaf temperature (K).
  t_k,                                                                         &
  ! Repair rate at t_ref (% of the damage threshold per minute).
  k_ref,                                                                       &
  ! Steepness of the fall in repair with accumulated damage.
  k_dec,                                                                       &
  ! High temperature threshold of repair (K).
  t_hi,                                                                        &
  ! Accumulated damage (fraction of the threshold).
  damage

REAL(KIND=real_jlslsm) :: k_rep

! Coefficient after the damage, and the undamaged share (%).
REAL(KIND=real_jlslsm) :: k_ref_now, alive

! Arrhenius shape (Arnold et al. 2025, Box 1): temperature, low and high
! deactivation temperatures, low threshold and reference temperature (K).
REAL(KIND=real_jlslsm), PARAMETER :: t_arr    = 3516.25
REAL(KIND=real_jlslsm), PARAMETER :: t_arr_lo = 50000.0
REAL(KIND=real_jlslsm), PARAMETER :: t_arr_hi = 83333.0
REAL(KIND=real_jlslsm), PARAMETER :: t_lo     = 283.65
REAL(KIND=real_jlslsm), PARAMETER :: t_ref    = 293.15

! The coefficient after the damage, from Eqn S6 at every damage level so
! it is continuous (an undamaged leaf repairs at k_ref k_ref^(k_dec/100));
! at full damage nothing is left to repair.
alive = 100.0 * (1.0 - damage)
IF (alive <= 0.0) THEN
  k_ref_now = 0.0
ELSE
  k_ref_now = k_ref * k_ref**(k_dec / alive)
END IF

! Arrhenius rise over the low- and high-temperature deactivation.
k_rep = k_ref_now * EXP(t_arr * (1.0 / t_ref - 1.0 / t_k))                     &
        / (1.0 + EXP(t_arr_lo * (1.0 / t_k - 1.0 / t_lo))                      &
               + EXP(t_arr_hi * (1.0 / t_hi - 1.0 / t_k)))

END FUNCTION repair_rate

! ---------------------------------------------------------------------
! Whether the root zone at point l is frozen this step: the frozen share of
! its soil water (worked out in physiol for the PFT whose stomata are being
! solved) above ximpair_frozen_frac.
!
! The same test stops xylem impairment damage in frozen soil (d7805aa), so
! the drought options and the impairment agree on what counts as frozen.
! Before physiol first sets ximpair_root_frozen, nothing is frozen.
! ---------------------------------------------------------------------
LOGICAL FUNCTION root_zone_frozen( l )

USE jules_vegetation_mod, ONLY: ximpair_frozen_frac
USE xylem_impairment_memory_mod, ONLY: ximpair_root_frozen

INTEGER, INTENT(IN) :: l

root_zone_frozen = .FALSE.
IF (ALLOCATED(ximpair_root_frozen))                                            &
  root_zone_frozen = ximpair_root_frozen(l) > ximpair_frozen_frac

END FUNCTION root_zone_frozen

END MODULE leaf_shed_mod
