! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_memory_mod
! Memory model of xylem impairment (pft_xylem_impairment_model = 3): the
! loss of conductance a drought locks in, its spread over new cover, and,
! with ring renewal and the ring cohorts (l_nsc_sapwood_rings), the stem's
! annual sapwood rings, each with its share of the sapwood and its loss.
! The ring functions give the stem's loss, the sapwood-area factor's sums
! and the Trugman response's repair. The state below is carried between
! updates; ximpair_lock, ximpair_repair_pend, ximpair_area_w (with
! l_veg_compete) and the ring arrays are written to and read from dumps.
!
! MDK, 06/11/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_MEMORY_MOD'

PRIVATE

PUBLIC :: ximpair_ring_init, ximpair_ring_loss_mean, ximpair_ring_grow,        &
          ximpair_ring_repair, ximpair_ring_cond_intact,                       &
          ximpair_ring_repairable, ximpair_ring_working,                       &
          ximpair_ring_cond_mean, ximpair_ring_cond_ref, ximpair_ring_next_year
PUBLIC :: ximpair_memory_alloc, ximpair_dilute_area,                           &
          leaf_conductance_impaired_memory_jls,                                &
          psi_leaf_impaired_memory,                                            &
          update_xylem_impairment_memory

PRIVATE :: k_intact, antideriv, psi_at_k, ring_cond, ring_age, ring_steady,   &
           ring_q,                                                             &
           ring_sums, ring_hi, ring_hit, ring_add_intact

! ximpair_memory_alloc allocates the first three arrays below.

! Locked-in loss of conductivity, 1 - k_cap/kmax, of each PFT at each land
! point. Read by hydraulic mortality, TRIFFID and its NSC routines, and
! the initialisation, and written to output. Not allocated until the first
! update (none during spin-up), and then read as zero where not yet set.
!
! MDK, 09/11/26
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_lock(:,:)

! Stem recovery by ring renewal (l_ximpair_rec_nsc): the share of the
! stem's conductance renewed by the wood TRIFFID grew, queued at the
! TRIFFID step by triffid_nsc_ring_renewal and applied (then cleared) at
! the next update, where it dilutes the stem's loss (with
! l_nsc_ring_additive it is taken off the loss).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_repair_pend(:,:)

! Factor (0-1) by which TRIFFID's gain of PFT area has diluted the stem
! damage (ximpair_dilute_area) since the last update, applied to the
! impaired kmax (kmax_impaired_pft, held by the surface code) at the next
! update and then reset to 1. In the dump, so a restart between a TRIFFID
! step and the next update is exact.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_area_w(:,:)

! Ring cohorts of the sapwood (l_nsc_sapwood_rings): for each annual ring
! of the conducting sapwood, youngest first (ring 1 is the ring being
! built), its share of the sapwood (< 0 in ring 1 = not yet set) and its
! loss of conductance (0-1), and the calendar year in which the youngest
! ring grows (< 0 = not yet set), which gives every ring's age in winters.
! Rings older than nsc_n_ring are merged into the last. In dumps;
! allocated only with the switch, with nsc_n_ring rings.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_frac(:,:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_loss(:,:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_year(:,:)

! Ring cohorts: each ring's sapwood in two parts, so wood a drought finds
! intact (grown, or gained as cover, since the ring was damaged) takes
! that drought's loss in full instead of being averaged into the damaged
! wood. The more damaged part holds the share d (0-1) of the ring at loss
! L_hi, the rest has the loss L_lo, and the ring's loss is their mean,
!
!   L = d L_hi + (1 - d) L_lo   (ximpair_ring_loss).
!
! ximpair_ring_dmg holds d (1 = one part, as for a ring set up with one
! loss) and ximpair_ring_loss_lo holds L_lo; L_hi follows from L
! (ring_hi). In dumps; allocated with the rings.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_dmg(:,:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_loss_lo(:,:,:)

! Ring cohorts: the highest loss of conductance (0-1) that the water
! potential has imposed on the rings since the last TRIFFID step, the
! damage the drought would put back at once. Xylem repair (the Trugman
! response) restores only the damage above it, and triffid_nsc_ring_renewal
! resets it to 0 after the repair. In dumps; allocated with the rings.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_ring_peak(:,:)

! Ring cohorts with the pipe heartwood rule (nsc_heartwood_model = 1): the
! stand's mean heartwood rate (per 360 days), a running mean of the rate
! at each TRIFFID step with an e-folding time of one year, set by
! triffid_nsc_ring_renewal. It is above nsc_heartwood_rate where the rule
! holds the sapwood at its target, and ring_steady builds the reference
! stand with it. < 0 = not set (nsc_heartwood_rate is used). In dumps
! with the rule; allocated with the rings.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_hw_rate_mean(:,:)

! With l_ximpair_frozen_hold: the share (0-1) of the root zone's soil
! water that is frozen, each layer weighted by its root fraction, at each
! land point for the PFT whose stomata are being solved. Set by physiol
! before the stomatal scheme on every step, so not in dumps.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE, PUBLIC :: ximpair_root_frozen(:)

CONTAINS

! ---------------------------------------------------------------------
! Xylem impairment, memory model: spread each PFT's drought damage over
! its new area after TRIFFID (the vegetation dynamics scheme) changes how
! much of the grid box the PFT covers (with l_veg_compete).
!
! New cover is new growth with undamaged xylem, so each loss of
! conductance L (0-1) is scaled by the share of the cover that was there
! before,
!
!   w = frac_old (1 - m) / frac_new,
!
! where frac_old and frac_new are the PFT's cover before and after the
! step and m is the fraction of the old cover killed in the step by
! hydraulic mortality (0 without it). Trees that die carry the mean
! damage, and the cover that grows in their place is undamaged, so with
! mortality w is applied wherever the cover is above that of the
! survivors, even if the net change in cover is 0 or less. Without
! mortality, only cover gained is diluted: a PFT that loses cover keeps
! the damage of what is left.
!
! The scaling applies to the stem's loss (ximpair_lock, and the damaged
! kmax through ximpair_area_w at the next update) and to the leaf
! segment's loss. The stem's share renewed by new rings and waiting to be
! applied, rho, multiplies the loss, L (1 - rho), and is the same for old
! and new cover, so it is not scaled. With additive recovery
! (l_nsc_ring_additive), MAX(L - rho, 0), the share renews only the cover
! that grew it, so it is scaled by w as well.
!
! The survivors' cover frac_old (1 - m) is approximate where TRIFFID does
! not remove exactly m frac_old (without l_ht_compete, or at the minimum
! cover frac_seed).
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_dilute_area( land_pts, trif_pts, trif_index, frac, dfrac,   &
                                m_kill )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype
USE jules_vegetation_mod, ONLY: xylem_impairment_memory
USE jules_vegetation_mod, ONLY: l_nsc_ring_additive
USE pftparm, ONLY: pft_xylem_impairment_model
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points, of TRIFFID points, and their land points.
INTEGER, INTENT(IN) :: land_pts, trif_pts, trif_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Areal fractions after the TRIFFID step.
  frac(land_pts,ntype),                                                        &
  ! Change in the areal fractions over the step.
  dfrac(land_pts,npft)
! Fraction of each PFT's old cover killed in the step (hydraulic
! mortality, l_trif_hydmort).
REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL :: m_kill(land_pts,npft)

! PFT, TRIFFID point and land point counters.
INTEGER :: n, t, l
! The old cover's share of the new cover, the share killed and the
! surviving cover.
REAL(KIND=real_jlslsm) :: w, m, frac_surv

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XIMPAIR_DILUTE_AREA'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( ALLOCATED(ximpair_lock) ) THEN
  CALL ximpair_memory_alloc( land_pts )
  ! Natural PFTs only: lotka changes their cover (crop cover is set
  ! elsewhere).
  DO n = 1,nnpft
    IF ( pft_xylem_impairment_model(n) /= xylem_impairment_memory ) CYCLE
    DO t = 1,trif_pts
      l = trif_index(t)
      IF ( frac(l,n) <= 0.0 ) CYCLE
      m = 0.0
      IF ( PRESENT(m_kill) ) m = m_kill(l,n)
      frac_surv = MAX(frac(l,n) - dfrac(l,n), 0.0) * (1.0 - m)
      IF ( m > 0.0 ) THEN
        IF ( frac(l,n) <= frac_surv ) CYCLE
      ELSE
        IF ( dfrac(l,n) <= 0.0 ) CYCLE
      END IF
      w = frac_surv / frac(l,n)
      ximpair_area_w(l,n)      = ximpair_area_w(l,n) * w
      ximpair_lock(l,n)        = ximpair_lock(l,n) * w
      IF ( ALLOCATED(ximpair_kcap_leaf) )                                      &
        ximpair_kcap_leaf(l,n) = 1.0 - (1.0 - ximpair_kcap_leaf(l,n)) * w
      ! Additive recovery: the queued share renews the old cover only.
      IF ( l_nsc_ring_additive )                                               &
        ximpair_repair_pend(l,n) = ximpair_repair_pend(l,n) * w
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ximpair_dilute_area

! ---------------------------------------------------------------------
! Xylem impairment, memory model: create the stored damage state for every
! land point and PFT if it does not exist yet (at the first update, or
! when a dump or initial condition is read). The stem's loss starts at 0,
! nothing is waiting to be renewed and the area dilution factor is 1. The
! sapwood's annual rings are created only with l_nsc_sapwood_rings, and
! are marked as not yet set.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_memory_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_rings, nsc_n_ring

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='XIMPAIR_MEMORY_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_repair_pend)) THEN
  ALLOCATE(ximpair_repair_pend(n_land_pts, npft))
  ximpair_repair_pend(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(ximpair_area_w)) THEN
  ALLOCATE(ximpair_area_w(n_land_pts, npft))
  ximpair_area_w(:,:) = 1.0
END IF

! The ring cohorts, only with them, not yet set.
IF (l_nsc_sapwood_rings .AND. .NOT. ALLOCATED(ximpair_ring_frac)) THEN
  ALLOCATE(ximpair_ring_frac(n_land_pts, npft, nsc_n_ring))
  ALLOCATE(ximpair_ring_loss(n_land_pts, npft, nsc_n_ring))
  ALLOCATE(ximpair_ring_year(n_land_pts, npft))
  ALLOCATE(ximpair_ring_peak(n_land_pts, npft))
  ALLOCATE(ximpair_hw_rate_mean(n_land_pts, npft))
  ALLOCATE(ximpair_ring_dmg(n_land_pts, npft, nsc_n_ring))
  ALLOCATE(ximpair_ring_loss_lo(n_land_pts, npft, nsc_n_ring))
  ximpair_ring_dmg(:,:,:)     = 1.0
  ximpair_ring_loss_lo(:,:,:) = 0.0
  ximpair_ring_frac(:,:,:) = -1.0
  ximpair_ring_loss(:,:,:) = 0.0
  ximpair_ring_year(:,:)   = -1.0
  ximpair_ring_peak(:,:)   = 0.0
  ximpair_hw_rate_mean(:,:) = -1.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ximpair_memory_alloc

! ---------------------------------------------------------------------
! Xylem impairment, memory model: set up the sapwood's annual rings at a
! point the first time they are used (l_nsc_sapwood_rings).
!
! The rings start in the steady state of the ring dynamics (ring_steady),
! the youngest ring grown as far as the season has gone (ring_q), and
! every ring starts with the stem's loss of conductance loss0 (0-1). The
! youngest ring grows in the current calendar year, or in the next one
! if the set-up falls after this year's ring year has ended (it is then
! empty; ximpair_ring_next_year).
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_ring_init( l, pft, loss0 )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! The stem's loss of conductance (0-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: loss0

! Steady-state ring shares (not normalised).
REAL(KIND=real_jlslsm) :: s(nsc_n_ring)

ximpair_ring_year(l,pft) = REAL(ximpair_ring_next_year(pft))

CALL ring_steady(l, pft, ring_q(l, pft), s)
ximpair_ring_frac(l,pft,:) = s(:) / SUM(s)
ximpair_ring_loss(l,pft,:) = loss0
ximpair_ring_dmg(l,pft,:)     = 1.0
ximpair_ring_loss_lo(l,pft,:) = 0.0

END SUBROUTINE ximpair_ring_init

! ---------------------------------------------------------------------
! Work out the ring shares of a stand in steady state under the ring
! dynamics of ximpair_ring_grow, used to set up a point's rings and the
! reference stand that the sapwood-area factor compares with.
!
! New wood G (kg C m-2 per 360 days) joins the youngest ring, TRIFFID's
! wood turnover g_t (g_wood, per 360 days) takes from all rings alike,
! and heartwood forms from the oldest rings at k_h X_s, with k_h the
! stand's mean heartwood rate: nsc_heartwood_rate, or with the pipe
! heartwood rule the point's running mean ximpair_hw_rate_mean, which is
! higher where the rule holds the sapwood at its target. A ring that grew
! m years ago then holds G EXP(-g_t m), and the sapwood holds
!
!   X_s = G / (k_h + g_t),
!
! so, in units of one year's new wood, the full rings j = 2, 3, ... have
! shares
!
!   s_j = EXP(-g_t (j - 1)),
!
! up to a total of 1 / (k_h + g_t), the last of them part full. The
! youngest ring (j = 1) has the share q, the part of this year's ring
! grown so far (0-1). Rings beyond nsc_n_ring are merged into the last,
! as ximpair_ring_grow does, so the last ring holds all the wood older
! than nsc_n_ring - 1 years. Without turnover there are 1 / k_h rings of
! equal carbon.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
SUBROUTINE ring_steady( l, pft, q, s )

USE pftparm, ONLY: nsc_heartwood_rate
USE trif, ONLY: g_wood
USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! The youngest ring's share (0-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: q
! Ring shares, youngest first (not normalised).
REAL(KIND=real_jlslsm), INTENT(OUT) :: s(nsc_n_ring)

! Ring counter.
INTEGER :: j
! Wood turnover rate and the stand's heartwood rate (per 360 days), the
! full rings' carbon still to place and a full ring's share.
REAL(KIND=real_jlslsm) :: g_t, k_h, rest, full

g_t = MAX(g_wood(pft), 0.0)
k_h = nsc_heartwood_rate(pft)
IF (ALLOCATED(ximpair_hw_rate_mean)) THEN
  IF (ximpair_hw_rate_mean(l,pft) >= 0.0) k_h = ximpair_hw_rate_mean(l,pft)
END IF
rest = 1.0 / MAX(k_h + g_t, 1.0e-6)
s(:) = 0.0
s(1) = q
DO j = 2, nsc_n_ring
  IF (rest <= 0.0) EXIT
  full = EXP(-g_t * REAL(j - 1))
  s(j) = MIN(full, rest)
  rest = rest - s(j)
END DO
IF (nsc_n_ring > 1) THEN
  s(nsc_n_ring) = s(nsc_n_ring) + MAX(rest, 0.0)
ELSE
  s(1) = 1.0
END IF

END SUBROUTINE ring_steady

! ---------------------------------------------------------------------
! Return the calendar year in which a ring started today grows: next
! year once this year's ring year has ended (on or after day
! nsc_grow_end_doy or, without it, the last day of the year), otherwise
! this year. At the TRIFFID step that holds the ring year's end this is
! next year, or this year when the step has run on into a new calendar
! year.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_next_year( pft ) RESULT( year )

USE pftparm, ONLY: is_set, nsc_grow_end_doy
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap

! PFT.
INTEGER, INTENT(IN) :: pft
! Calendar year.
INTEGER :: year

! Day of year, days in this year and the ring year's end day.
INTEGER :: doy, n_days, doy_end

n_days  = days_in_year(current_time%year, l_360, l_leap)
doy     = day_of_year(current_time%year, current_time%month,                  &
                      current_time%day, l_360, l_leap)
doy_end = n_days
IF (is_set(nsc_grow_end_doy(pft)))                                             &
  doy_end = MIN(NINT(nsc_grow_end_doy(pft)), n_days)

year = current_time%year
IF (doy >= doy_end) year = year + 1

END FUNCTION ximpair_ring_next_year

! ---------------------------------------------------------------------
! Return the part (0-1) of the youngest sapwood ring that a stand in
! steady state has grown by this day: the fraction of the growing season
! (nsc_season_start, the leaf flush or nsc_grow_start_doy, to
! nsc_grow_end_doy, as the store's release), or of the calendar year if
! those are not set, that has passed; 0 for a ring not grown yet and 1
! for an older youngest ring.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
FUNCTION ring_q( l, pft ) RESULT( q )

USE pftparm, ONLY: is_set, nsc_season_start, nsc_grow_end_doy
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year, days_in_year
USE time_info_mod, ONLY: l_360, l_leap

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! The youngest ring's share (0-1).
REAL(KIND=real_jlslsm) :: q

! The youngest ring's age and the day of year.
INTEGER :: a1, doy

a1  = ring_age(l, pft, 1)
doy = day_of_year(current_time%year, current_time%month, current_time%day,   &
                  l_360, l_leap)
IF (a1 < 0) THEN
  q = 0.0
ELSE IF (a1 > 0) THEN
  q = 1.0
ELSE IF (is_set(nsc_season_start(pft)) .AND.                                   &
         is_set(nsc_grow_end_doy(pft))) THEN
  q = (REAL(doy) - nsc_season_start(pft))                                      &
      / MAX(nsc_grow_end_doy(pft) - nsc_season_start(pft), 1.0)
ELSE
  q = REAL(doy) / REAL(days_in_year(current_time%year, l_360, l_leap))
END IF
q = MAX(MIN(q, 1.0), 0.0)

END FUNCTION ring_q

! ---------------------------------------------------------------------
! Return the age of sapwood ring i (1 = youngest), which sets how well
! the ring conducts (ring_cond), as the number of winters it has passed,
!
!   age = (current year - year the youngest ring grows) + i - 1.
!
! When the ring year ends, a new empty youngest ring starts that grows
! next year, so its age is -1 (not grown yet), and this year's ring, now
! ring 2, keeps age 0 until the calendar year ends. Rings get older over
! winter, whatever the tree's growth calendar.
!
! MDK, 31/10/26
! ---------------------------------------------------------------------
FUNCTION ring_age( l, pft, i ) RESULT( age )

USE model_time_mod, ONLY: current_time

! Land point, PFT and ring.
INTEGER, INTENT(IN) :: l, pft, i
! Age (winters; < 0 for a ring not grown yet).
INTEGER :: age

age = current_time%year - NINT(ximpair_ring_year(l,pft)) + i - 1

END FUNCTION ring_age

! ---------------------------------------------------------------------
! Return how well a sapwood ring of age a (winters) conducts water
! compared with a new ring (0-1), used to weight each ring in the stem's
! conductance. Older rings conduct less,
!
!   c(a) = EXP(-a / e),
!
! where e (nsc_ring_cond_efold, winters) is how fast it falls, so the
! youngest rings carry most of the flow whatever the number of rings
! followed. If e is not set, every ring conducts alike, c = 1. With a
! reference age a0 above 0, c is taken relative to the ring of age a0,
!
!   c(a) = EXP(-(a - a0) / e),
!
! which gives the same ratios between rings without underflow.
!
! c = 0 for a ring not grown yet (a < 0). In a ring-porous tree
! (nsc_ring_func_doy set), this year's ring (a = 0) also has c = 0 before
! day nsc_ring_func_doy of the calendar year, the day its new vessels
! start to conduct.
!
! MDK, 04/11/26
! ---------------------------------------------------------------------
FUNCTION ring_cond( pft, age, a0 ) RESULT( c )

USE pftparm, ONLY: is_set, nsc_ring_cond_efold, nsc_ring_func_doy
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap

! PFT, the ring's age and the reference age (winters).
INTEGER, INTENT(IN) :: pft, age, a0
! Relative conductivity (0-1).
REAL(KIND=real_jlslsm) :: c

! Day of year now.
INTEGER :: doy

c = 0.0
IF (age < 0) RETURN
IF (age == 0 .AND. is_set(nsc_ring_func_doy(pft))) THEN
  doy = day_of_year(current_time%year, current_time%month,                     &
                    current_time%day, l_360, l_leap)
  IF (doy < NINT(nsc_ring_func_doy(pft))) RETURN
END IF

c = 1.0
IF (is_set(nsc_ring_cond_efold(pft))) THEN
  IF (a0 > 0) THEN
    c = EXP(-REAL(age - a0) / nsc_ring_cond_efold(pft))
  ELSE
    c = EXP(-REAL(age) / nsc_ring_cond_efold(pft))
  END IF
END IF

END FUNCTION ring_cond

! ---------------------------------------------------------------------
! Work out three sums over the sapwood rings at point l that the ring
! functions below share, each ring i weighted by its share of the sapwood
! f_i and how well it conducts c_i (ring_cond),
!
!   sw  = SUM_i f_i c_i,
!   swl = SUM_i f_i c_i L_i,
!   swf = SUM_i f_i c_i (1 - L_i),
!
! where L_i is the ring's loss of conductance (0-1).
!
! With l_rel, c_i is taken relative to the youngest ring that conducts.
! The ratios between the sums do not
! change, and that ring keeps the weight f_i, so a very small e in
! ring_cond cannot make every weight 0.
!
! MDK, 31/10/26
! ---------------------------------------------------------------------
SUBROUTINE ring_sums( l, pft, l_rel, sw, swl, swf )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Weights relative to the youngest conducting ring.
LOGICAL, INTENT(IN) :: l_rel
! The sums.
REAL(KIND=real_jlslsm), INTENT(OUT) :: sw, swl, swf

! Ring counter, a ring's age and the reference age.
INTEGER :: i, age, a0
! A ring's weight.
REAL(KIND=real_jlslsm) :: w

a0 = 0
IF (l_rel) THEN
  a0 = -1
  DO i = 1, nsc_n_ring
    IF (ximpair_ring_frac(l,pft,i) <= 0.0) CYCLE
    age = ring_age(l, pft, i)
    IF (ring_cond(pft, age, 0) > 0.0 .OR. age > 0) THEN
      a0 = MAX(age, 0)
      EXIT
    END IF
  END DO
  a0 = MAX(a0, 0)
END IF

sw  = 0.0
swl = 0.0
swf = 0.0
DO i = 1, nsc_n_ring
  IF (ximpair_ring_frac(l,pft,i) <= 0.0) CYCLE
  w = ximpair_ring_frac(l,pft,i) * ring_cond(pft, ring_age(l, pft, i), a0)
  sw  = sw  + w
  swl = swl + w * ximpair_ring_loss(l,pft,i)
  swf = swf + w * (1.0 - ximpair_ring_loss(l,pft,i))
END DO

END SUBROUTINE ring_sums

! ---------------------------------------------------------------------
! Return the loss of conductance (0-1) of the more damaged part of ring
! i at point l (ximpair_ring_dmg),
!
!   L_hi = (L - (1 - d) L_lo) / d,
!
! with L the ring's loss, d the part's share and L_lo the loss of the
! rest; L_lo where the ring has no damaged part (d = 0).
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ring_hi( l, pft, i ) RESULT( l_hi )

! Land point, PFT and ring.
INTEGER, INTENT(IN) :: l, pft, i
! The more damaged part's loss (0-1).
REAL(KIND=real_jlslsm) :: l_hi

! The part's share.
REAL(KIND=real_jlslsm) :: d

d = ximpair_ring_dmg(l,pft,i)
IF (d > 1.0e-6) THEN
  l_hi = (ximpair_ring_loss(l,pft,i)                                           &
          - (1.0 - d) * ximpair_ring_loss_lo(l,pft,i)) / d
  l_hi = MIN(MAX(l_hi, ximpair_ring_loss_lo(l,pft,i)), 1.0)
ELSE
  l_hi = ximpair_ring_loss_lo(l,pft,i)
END IF

END FUNCTION ring_hi

! ---------------------------------------------------------------------
! Damage ring i at point l by a step that imposes the loss of
! conductance l_imp (0-1). Each part of the ring takes at least that
! loss,
!
!   L_hi <- MAX(L_hi, l_imp),   L_lo <- MAX(L_lo, l_imp),
!
! so wood a drought finds intact takes its loss in full. Where the less
! damaged part reaches the other, the ring is one part again (d = 1).
! With d = 1 this is L <- MAX(L, l_imp).
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ring_hit( l, pft, i, l_imp )

! Land point, PFT and ring.
INTEGER, INTENT(IN) :: l, pft, i
! The loss the step imposes (0-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: l_imp

! The parts' losses and the more damaged part's share.
REAL(KIND=real_jlslsm) :: l_hi, l_lo, d

d    = ximpair_ring_dmg(l,pft,i)
l_hi = MAX(ring_hi(l, pft, i), l_imp)
l_lo = MAX(ximpair_ring_loss_lo(l,pft,i), l_imp)
IF (d >= 1.0) THEN
  ximpair_ring_loss(l,pft,i) = MAX(ximpair_ring_loss(l,pft,i), l_imp)
ELSE IF (l_lo >= l_hi) THEN
  ximpair_ring_loss(l,pft,i)    = l_lo
  ximpair_ring_dmg(l,pft,i)     = 1.0
  ximpair_ring_loss_lo(l,pft,i) = 0.0
ELSE
  ximpair_ring_loss(l,pft,i)    = d * l_hi + (1.0 - d) * l_lo
  ximpair_ring_loss_lo(l,pft,i) = l_lo
END IF

END SUBROUTINE ring_hit

! ---------------------------------------------------------------------
! Mix intact wood into ring i at point l, so the ring's old wood keeps
! the share w (0-1) of it: new wood joining the youngest ring, or new
! cover joining the stand. The intact wood joins the less damaged part,
!
!   L <- w L,   d <- w d,
!   L_lo <- L_lo w (1 - d) / (w (1 - d) + 1 - w),
!
! so the more damaged part keeps its loss.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ring_add_intact( l, pft, i, w )

! Land point, PFT and ring.
INTEGER, INTENT(IN) :: l, pft, i
! The old wood's share of the ring after the mixing (0-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: w

! The less damaged part's old wood, and the part after the mixing
! (shares of the ring).
REAL(KIND=real_jlslsm) :: f_lo, f_new

f_lo  = w * (1.0 - ximpair_ring_dmg(l,pft,i))
f_new = f_lo + 1.0 - w
IF (f_new > 0.0) THEN
  ximpair_ring_loss_lo(l,pft,i) = ximpair_ring_loss_lo(l,pft,i) * f_lo / f_new
ELSE
  ximpair_ring_loss_lo(l,pft,i) = 0.0
END IF
ximpair_ring_loss(l,pft,i) = ximpair_ring_loss(l,pft,i) * w
ximpair_ring_dmg(l,pft,i)  = ximpair_ring_dmg(l,pft,i) * w

END SUBROUTINE ring_add_intact

! ---------------------------------------------------------------------
! Return the stem's loss of conductance (0-1) from drought damage to its
! annual rings, which sets the stem's conductance cap in the memory model
! of xylem impairment. Each ring's loss L_i is weighted by its share of
! the sapwood f_i and by how well it conducts c_i (ring_cond),
!
!   L = SUM_i f_i c_i L_i / SUM_i f_i c_i.
!
! If no ring conducts (only a ring-porous tree's new ring, before its
! vessels conduct, holds sapwood), L is the mean of the older rings'
! losses weighted by f_i alone, or 0 if they hold no sapwood.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_loss_mean( l, pft ) RESULT( loss )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Loss of conductance (0-1).
REAL(KIND=real_jlslsm) :: loss

! Ring counter.
INTEGER :: i
! The sums, and the plain sums of shares and share-weighted losses.
REAL(KIND=real_jlslsm) :: sw, swl, swf, sum_f, sum_fl

CALL ring_sums(l, pft, .TRUE., sw, swl, swf)
loss = 0.0
IF (sw > 0.0) THEN
  loss = swl / sw
ELSE
  sum_f  = 0.0
  sum_fl = 0.0
  DO i = 2, nsc_n_ring
    IF (ximpair_ring_frac(l,pft,i) <= 0.0) CYCLE
    sum_f  = sum_f + ximpair_ring_frac(l,pft,i)
    sum_fl = sum_fl + ximpair_ring_frac(l,pft,i) * ximpair_ring_loss(l,pft,i)
  END DO
  IF (sum_f > 0.0) loss = sum_fl / sum_f
END IF

END FUNCTION ximpair_ring_loss_mean

! ---------------------------------------------------------------------
! Return the sapwood's mean ability to conduct water compared with new
! wood (0-1), which scales the plant's conductance with its sapwood
! (nsc_k_area),
!
!   SUM_i f_i c_i,
!
! where f_i is ring i's share of the sapwood and c_i how well it conducts
! on this day (ring_cond). With X_s the sapwood carbon and X_i = f_i X_s
! that of ring i, SUM_i X_i c_i is X_s times this value.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_cond_mean( l, pft ) RESULT( share )

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Weighted share (0-1).
REAL(KIND=real_jlslsm) :: share

! The other sums.
REAL(KIND=real_jlslsm) :: swl, swf

CALL ring_sums(l, pft, .FALSE., share, swl, swf)

END FUNCTION ximpair_ring_cond_mean

! ---------------------------------------------------------------------
! Return the share of the sapwood that still conducts water (0-1), for
! output. Each ring's share f_i counts by the part that is undamaged,
! 1 - L_i, and by how well it conducts compared with a new ring, c_i
! (ring_cond),
!
!   X_f / X_s = SUM_i f_i (1 - L_i) c_i,
!
! where X_f is the conducting sapwood and X_s all the sapwood.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_cond_intact( l, pft ) RESULT( share )

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Functional share of the sapwood (0-1).
REAL(KIND=real_jlslsm) :: share

! The other sums.
REAL(KIND=real_jlslsm) :: sw, swl

CALL ring_sums(l, pft, .FALSE., sw, swl, share)

END FUNCTION ximpair_ring_cond_intact

! ---------------------------------------------------------------------
! Return the mean ability to conduct water (0-1) of the sapwood of a
! stand in steady state on the same day, which nsc_k_area compares with
! the point's own rings to scale the plant's conductance.
!
! The steady stand has the ring shares s_j of ring_steady (the ring
! dynamics with TRIFFID's wood turnover and the stand's mean heartwood
! rate) at the ages of the point's rings (ring_age), its youngest ring
! grown as far as the season has gone (ring_q). With c(a_j) how well a
! ring of its age conducts (ring_cond),
!
!   c_ref = SUM_j s_j c(a_j) / SUM_j s_j,
!
! so a stand at its target sapwood X_t with this ring structure has
! SUM_i X_i c_i = X_t c_ref on any date, whether its heartwood forms at
! nsc_heartwood_rate or faster under the pipe heartwood rule.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_cond_ref( l, pft ) RESULT( c_ref )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Equilibrium mean conductivity (0-1).
REAL(KIND=real_jlslsm) :: c_ref

! Ring counter and the youngest ring's age.
INTEGER :: j, a1
! The steady stand's ring shares and their sum.
REAL(KIND=real_jlslsm) :: s(nsc_n_ring), sum_s

a1 = ring_age(l, pft, 1)
CALL ring_steady(l, pft, ring_q(l, pft), s)
c_ref = 0.0
sum_s = 0.0
DO j = 1, nsc_n_ring
  IF (s(j) <= 0.0) CYCLE
  c_ref = c_ref + s(j) * ring_cond(pft, a1 + j - 1, 0)
  sum_s = sum_s + s(j)
END DO
IF (sum_s > 0.0) c_ref = c_ref / sum_s

END FUNCTION ximpair_ring_cond_ref

! ---------------------------------------------------------------------
! Xylem impairment, memory model: update the sapwood's annual rings over
! one TRIFFID step (l_nsc_sapwood_rings).
!
! 1. The step's new wood dg (kg C m-2) joins the youngest ring
!    undamaged, which dilutes that ring's loss of conductance.
! 2. Sapwood turns to heartwood, oldest rings first. The amount is
!
!      x_turn (1 - EXP(-f_heart)),
!
!    where x_turn (kg C m-2) is the sapwood before the step, after
!    TRIFFID's turnover (taken from all rings alike), and f_heart is the
!    heartwood rate times the step length. The sapwood pool loses the
!    same amount in total.
! 3. At the end of the ring year (l_shift), the two oldest rings merge
!    and a new, empty youngest ring starts, to grow in calendar year
!    g_year. With only one ring, that ring carries on.
! 4. Optionally, the share of the stem's conductance weight that the new
!    wood holds after step 2,
!
!      w_new = (dg / X) c_1 / SUM_i f_i c_i,
!
!    with X the sapwood after the step and c_1 how well the youngest ring
!    conducts on this day: the part of the stem's loss the new wood
!    dilutes (0 for a ring-porous tree's new ring before its vessels
!    conduct).
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_ring_grow( l, pft, x_turn, dg, f_heart, l_shift, g_year,    &
                              w_new )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Sapwood before the step and the new wood (kg C m-2); heartwood rate
! times the step (-).
REAL(KIND=real_jlslsm), INTENT(IN) :: x_turn, dg, f_heart
! The ring year ends with this step.
LOGICAL, INTENT(IN) :: l_shift
! The calendar year in which a ring started at this step grows.
INTEGER, INTENT(IN) :: g_year
! The new wood's share of the conductance weight (0-1).
REAL(KIND=real_jlslsm), INTENT(OUT), OPTIONAL :: w_new

! Ring counter and the number of rings.
INTEGER :: i, n
! The ring sums (ring_sums).
REAL(KIND=real_jlslsm) :: sw, swl, swf
! Sapwood of each ring (kg C m-2), ring 1 with the new wood, heartwood still
! to take, the part taken from a ring, the total, and the merged rings'
! share.
REAL(KIND=real_jlslsm) :: a(nsc_n_ring), x1_new, h, take, tot, f_old
! Merged rings: the shares in the more and less damaged parts, and the
! less damaged part's loss.
REAL(KIND=real_jlslsm) :: f_hi, f_lo, lo_new

n = nsc_n_ring
a(:) = MAX(ximpair_ring_frac(l,pft,:), 0.0) * x_turn

! The new wood is intact: it dilutes the youngest ring's loss, and joins
! its less damaged part.
IF (dg > 0.0) THEN
  x1_new = a(1) + dg
  CALL ring_add_intact(l, pft, 1, a(1) / x1_new)
  a(1) = x1_new
END IF

! Heartwood from the oldest rings first.
h = x_turn * (1.0 - EXP(-f_heart))
DO i = n, 1, -1
  IF (h <= 0.0) EXIT
  take = MIN(a(i), h)
  a(i) = a(i) - take
  h = h - take
END DO
tot = SUM(a)
IF (tot > 0.0) ximpair_ring_frac(l,pft,:) = a(:) / tot

! The new wood's share of the conductance weight.
IF (PRESENT(w_new)) THEN
  w_new = 0.0
  IF (dg > 0.0 .AND. tot > 0.0) THEN
    CALL ring_sums(l, pft, .FALSE., sw, swl, swf)
    IF (sw > 0.0) w_new = MIN(dg / tot, 1.0)                                   &
                          * ring_cond(pft, ring_age(l, pft, 1), 0) / sw
  END IF
END IF

! End of the ring year: the oldest two rings merge and a new ring starts.
IF (l_shift) THEN
  IF (n > 1) THEN
    f_old = ximpair_ring_frac(l,pft,n) + ximpair_ring_frac(l,pft,n-1)
    IF (f_old > 0.0) THEN
      ! The merged ring's parts: the more damaged parts together, and the
      ! rest together, each weighted by its carbon.
      f_hi = ximpair_ring_frac(l,pft,n) * ximpair_ring_dmg(l,pft,n)            &
             + ximpair_ring_frac(l,pft,n-1) * ximpair_ring_dmg(l,pft,n-1)
      f_lo = f_old - f_hi
      lo_new = 0.0
      IF (f_lo > 0.0)                                                          &
        lo_new = (ximpair_ring_frac(l,pft,n)                                   &
                  * (1.0 - ximpair_ring_dmg(l,pft,n))                          &
                  * ximpair_ring_loss_lo(l,pft,n)                              &
                  + ximpair_ring_frac(l,pft,n-1)                               &
                  * (1.0 - ximpair_ring_dmg(l,pft,n-1))                        &
                  * ximpair_ring_loss_lo(l,pft,n-1)) / f_lo
      ximpair_ring_loss(l,pft,n) = (ximpair_ring_frac(l,pft,n)                 &
                                    * ximpair_ring_loss(l,pft,n)               &
                                    + ximpair_ring_frac(l,pft,n-1)             &
                                    * ximpair_ring_loss(l,pft,n-1)) / f_old
      ximpair_ring_dmg(l,pft,n)     = MIN(MAX(f_hi / f_old, 0.0), 1.0)
      ximpair_ring_loss_lo(l,pft,n) = lo_new
    END IF
    ximpair_ring_frac(l,pft,n) = f_old
    DO i = n - 1, 2, -1
      ximpair_ring_frac(l,pft,i)    = ximpair_ring_frac(l,pft,i-1)
      ximpair_ring_loss(l,pft,i)    = ximpair_ring_loss(l,pft,i-1)
      ximpair_ring_dmg(l,pft,i)     = ximpair_ring_dmg(l,pft,i-1)
      ximpair_ring_loss_lo(l,pft,i) = ximpair_ring_loss_lo(l,pft,i-1)
    END DO
    ximpair_ring_frac(l,pft,1)    = 0.0
    ximpair_ring_loss(l,pft,1)    = 0.0
    ximpair_ring_dmg(l,pft,1)     = 1.0
    ximpair_ring_loss_lo(l,pft,1) = 0.0
  END IF
  ximpair_ring_year(l,pft) = REAL(g_year)
END IF

END SUBROUTINE ximpair_ring_grow

! ---------------------------------------------------------------------
! Xylem impairment, memory model: replace damaged sapwood with new xylem in
! the annual rings (Trugman response). Repair is costed by the conductance
! it restores, as new wood: a unit of carbon restores the conductance of a
! unit of new ring (c = 1). With X the sapwood (kg C m-2), X_i = f_i X ring
! i's carbon, c_i how well it conducts (ring_cond), L_i its loss of
! conductance and L_p the highest loss imposed since the last TRIFFID step
! (ximpair_ring_peak), the conductance repair can restore, in kg C of new
! ring, is
!
!   D = SUM_i X_i c_i MAX(L_i - L_p, 0).
!
! Damage at or below L_p is left, as the drought would put it back at the
! next time step. r of new xylem (kg C m-2), at most r_max and D, restores
! each ring's share of D, so its loss becomes
!
!   L_i -> L_i - r MAX(L_i - L_p, 0) / D,
!
! at least L_p (L_p in every ring when r = D). No ring changes size, so the
! sapwood keeps its age structure and a ring conducts at most as well as it
! did undamaged. r is returned; r_max - r is claimed but not used.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ximpair_ring_repair( l, pft, x_s, r_max, r )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Sapwood and the most new xylem the repair may use (kg C m-2).
REAL(KIND=real_jlslsm), INTENT(IN) :: x_s, r_max
! New xylem used to replace damaged wood (kg C m-2).
REAL(KIND=real_jlslsm), INTENT(OUT) :: r

! Ring counter.
INTEGER :: i
! Each ring's carbon, and the conductance repair can restore in it and in
! all rings (kg C m-2 of new ring).
REAL(KIND=real_jlslsm) :: a(nsc_n_ring), d(nsc_n_ring), d_tot
! A ring's loss before the repair, and the factor its loss above L_p is
! cut by.
REAL(KIND=real_jlslsm) :: l_old, phi

r = 0.0
IF (x_s <= 0.0 .OR. r_max <= 0.0) RETURN
IF (ximpair_ring_frac(l,pft,1) < 0.0) RETURN

a(:) = MAX(ximpair_ring_frac(l,pft,:), 0.0) * x_s
DO i = 1, nsc_n_ring
  d(i) = a(i) * ring_cond(pft, ring_age(l, pft, i), 0)                         &
         * MAX(ximpair_ring_loss(l,pft,i) - ximpair_ring_peak(l,pft), 0.0)
END DO
d_tot = SUM(d)
IF (d_tot <= 0.0) RETURN
r = MIN(r_max, d_tot)

! Each ring's repairable loss falls by the share r / D. In a ring in two
! parts, each part's loss above L_p falls by the same factor phi (the
! less damaged part's, if it is above L_p), so the more damaged part is
! repaired most.
DO i = 1, nsc_n_ring
  IF (d(i) <= 0.0) CYCLE
  l_old = ximpair_ring_loss(l,pft,i)
  ximpair_ring_loss(l,pft,i) = MIN(MAX(ximpair_ring_loss(l,pft,i)             &
                                       - r * (ximpair_ring_loss(l,pft,i)       &
                                       - ximpair_ring_peak(l,pft)) / d_tot,    &
                                       ximpair_ring_peak(l,pft)), 1.0)
  IF (ximpair_ring_dmg(l,pft,i) < 1.0) THEN
    phi = (ximpair_ring_loss(l,pft,i) - ximpair_ring_peak(l,pft))              &
          / (l_old - ximpair_ring_peak(l,pft))
    IF (ximpair_ring_loss_lo(l,pft,i) > ximpair_ring_peak(l,pft))              &
      ximpair_ring_loss_lo(l,pft,i) = ximpair_ring_peak(l,pft)                 &
                                      + phi * (ximpair_ring_loss_lo(l,pft,i)   &
                                      - ximpair_ring_peak(l,pft))
    IF (ring_hi(l, pft, i) <= ximpair_ring_loss_lo(l,pft,i)) THEN
      ximpair_ring_dmg(l,pft,i)     = 1.0
      ximpair_ring_loss_lo(l,pft,i) = 0.0
    END IF
  END IF
END DO

END SUBROUTINE ximpair_ring_repair

! ---------------------------------------------------------------------
! Return the conductance repair could restore, per unit of sapwood carbon
! and in units of new ring (0-1),
!
!   SUM_i f_i c_i MAX(L_i - L_p, 0),
!
! with f_i ring i's share of the sapwood, c_i how well it conducts
! (ring_cond), L_i its loss of conductance and L_p the highest loss the
! water potential has imposed since the last TRIFFID step
! (ximpair_ring_peak): damage at or below L_p would come back at once
! while the drought lasts. Times the sapwood X_s it is the carbon that
! repair needs (ximpair_ring_repair). 0 where the rings are not yet set.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_repairable( l, pft ) RESULT( share )

USE jules_vegetation_mod, ONLY: nsc_n_ring

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Repairable conductance per unit of sapwood carbon (0-1).
REAL(KIND=real_jlslsm) :: share

! Ring counter.
INTEGER :: i

share = 0.0
IF (.NOT. ALLOCATED(ximpair_ring_frac)) RETURN
IF (ximpair_ring_frac(l,pft,1) < 0.0) RETURN
DO i = 1, nsc_n_ring
  share = share + MAX(ximpair_ring_frac(l,pft,i), 0.0)                         &
                  * ring_cond(pft, ring_age(l, pft, i), 0)                     &
                  * MAX(ximpair_ring_loss(l,pft,i) - ximpair_ring_peak(l,pft), &
                        0.0)
END DO
share = MIN(MAX(share, 0.0), 1.0)

END FUNCTION ximpair_ring_repairable

! ---------------------------------------------------------------------
! Return the share of the sapwood's conductance that still works (0-1),
!
!   1 - L = SUM_i f_i c_i (1 - L_i) / SUM_i f_i c_i,
!
! with L the stem's loss (ximpair_ring_loss_mean), f_i ring i's share of the
! sapwood, c_i how well it conducts (ring_cond) and L_i its loss of
! conductance; -1 where the rings are not yet set.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
FUNCTION ximpair_ring_working( l, pft ) RESULT( share )

! Land point and PFT.
INTEGER, INTENT(IN) :: l, pft
! Working share of the conductance (0-1; -1 = not set).
REAL(KIND=real_jlslsm) :: share

share = -1.0
IF (.NOT. ALLOCATED(ximpair_ring_frac)) RETURN
IF (ximpair_ring_frac(l,pft,1) < 0.0 .OR. ximpair_ring_year(l,pft) < 0.0)     &
  RETURN
share = 1.0 - ximpair_ring_loss_mean(l, pft)

END FUNCTION ximpair_ring_working

! *********************************************************************
! Xylem impairment, memory model (pft_xylem_impairment_model = 3).
!
! In drought, air enters xylem vessels (embolism) and they stop carrying
! water. The vessels do not refill while the xylem is under tension, so
! the model keeps the conductance lost at the driest water potential met
! as a cap on the vulnerability curve,
!
!   k(psi) = MIN(k_intact(psi), k_cap),   k_cap = kmax (1 - L),
!
! where k_intact is the undamaged curve (mol m-2 s-1 Pa-1), psi the water
! potential (Pa, below 0) and L the loss of conductance (0-1). The dry end
! of the curve is unchanged and only the wet end is lowered, as in Sperry
! et al. (2017) and Venturas et al. (2018). k_cap is held in the damaged
! kmax array (kmax_impaired_pft, leaf area basis). conductance_b and
! conductance_c are not changed.
!
! Each step, update_xylem_impairment_memory first lets the cap recover
! (new wood renewing the stem, new leaves renewing the leaf segment, or
! a reset once a year, each optional) and then lowers it to the undamaged
! curve's conductance at the water potential that drives the damage.
!
! The capped curve is the target curve of Mackay et al. (2015, Water
! Resources Research, eqs 11-12) in the TREES model, Kcav = Ksat
! f(psi_min) with the undamaged curve below psi_min. TREES fits a new
! Weibull curve to it, while this model uses the capped curve as it is.
!
! With the cap, the transpiration E (mol H2O m-2 s-1) carried from the
! root-zone water potential psi_r to the leaf water potential psi_l is
!
!   E = k_cap (psi_r - MAX(psi_l, psi_s))
!       + A(psi_l) - A(MIN(psi_r, psi_s)),
!
! where the second line counts only if psi_l < psi_s. psi_s is the water
! potential at which k_intact = k_cap: wetter than psi_s the cap holds,
! drier than psi_s the undamaged curve does. A(psi_1) - A(psi_2) is the
! integral of k_intact from psi_1 to psi_2.
! *********************************************************************

! ---------------------------------------------------------------------
! Return the conductance of the undamaged xylem (mol m-2 s-1 Pa-1) at n
! water potentials psi (Pa) at one point, which the memory model caps to
! represent drought damage. It uses the cumulative Weibull curve,
! kmax EXP(-(psi / b)^c), or the SOX curve, kmax / (1 + (psi / b)^c).
!
! MDK, 01/11/26
! ---------------------------------------------------------------------
FUNCTION k_intact( pft, n, psi, kmax, b, c ) RESULT( k )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance

INTEGER, INTENT(IN) :: pft, n
REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n), kmax, b, c
REAL(KIND=real_jlslsm) :: k(n)

IF (pft_conductance_model(pft) == CW_conductance) THEN
  k = kmax * EXP( -(psi/b)**c )
ELSE
  k = kmax / (1 + (psi/b)**c)
END IF

END FUNCTION k_intact

! ---------------------------------------------------------------------
! Return A(psi), a function of water potential (Pa) that the memory
! model's solver uses to find the flow the undamaged xylem carries between
! two water potentials,
!
!   A(psi_1) - A(psi_2) = integral of k_intact from psi_1 to psi_2
!                         (mol m-2 s-1),
!
! so dA/dpsi = -k_intact. It uses the incomplete gamma function for the
! cumulative Weibull curve and the 2F1 function for the SOX curve, as the
! Newton solvers in psi_leaf_CW_jls and psi_leaf_SOX_jls do.
!
! MDK, 24/10/26
! ---------------------------------------------------------------------
FUNCTION antideriv( pft, n, psi, kmax, b, c ) RESULT( a )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance
USE xylem_hydraulics_CW_jls_mod, ONLY: incomplete_gamma
USE xylem_hydraulics_SOX_jls_mod, ONLY: SOX_2F1

INTEGER, INTENT(IN) :: pft, n
REAL(KIND=real_jlslsm), INTENT(IN) :: psi(n), kmax, b, c
REAL(KIND=real_jlslsm) :: a(n)

IF (pft_conductance_model(pft) == CW_conductance) THEN
  ! Note: the negative sign is present because b is negative.
  a = kmax * (-b/c) * incomplete_gamma(n, 1/c, (psi/b)**c)
ELSE
  a = - kmax * psi * SOX_2F1(n, 1 + 1/c, -(psi/b)**c)
END IF

END FUNCTION antideriv

! ---------------------------------------------------------------------
! Return the water potential (Pa) at which the undamaged curve falls to
! r kmax (0 < r < 1), the inverse of k_intact, which the memory model's
! solver uses to find where the damaged and undamaged curves meet.
! ---------------------------------------------------------------------
FUNCTION psi_at_k( pft, r, b, c ) RESULT( psi )

USE pftparm, ONLY: pft_conductance_model
USE jules_vegetation_mod, ONLY: CW_conductance

INTEGER, INTENT(IN) :: pft
REAL(KIND=real_jlslsm), INTENT(IN) :: r, b, c
REAL(KIND=real_jlslsm) :: psi

IF (pft_conductance_model(pft) == CW_conductance) THEN
  psi = b * (-LOG(r))**(1/c)
ELSE
  psi = b * (1/r - 1)**(1/c)
END IF

END FUNCTION psi_at_k

! ---------------------------------------------------------------------
! Work out the xylem conductance (mol m-2 s-1 Pa-1) of a plant damaged by
! drought at a given water potential at every land point, under the
! memory model of xylem impairment. It is the capped curve,
!
!   k = MIN(k_intact(psi), MAX(k_cap, kcrit)),
!
! where kcrit is the lowest conductance the model allows.
! ---------------------------------------------------------------------
SUBROUTINE leaf_conductance_impaired_memory_jls( pft,                          &
                                                 land_pnts,                    &
                                                 water_potential,              &
                                                 kmax_ref,                     &
                                                 kcap,                         &
                                                 kcrit,                        &
                                                 conductance_b,                &
                                                 conductance_c,                &
                                               ! INTENT OUT
                                                 leaf_k                        &
  )

USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) ::                                                         &
  pft                                                                          &
                            ! Plant functional type index
, land_pnts
                            ! Number of land points

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  water_potential(land_pnts)                                                   &
                            ! Water potential for each land point (Pa)
, kmax_ref(land_pnts)                                                          &
                            ! Unimpaired maximum xylem conductance
                            ! (mol m-2 s-1 Pa-1)
, kcap(land_pnts)                                                              &
                            ! Conductance cap, kmax * (1 - PLC_mem)
                            ! (mol m-2 s-1 Pa-1)
, kcrit(land_pnts)                                                             &
                            ! Critical xylem conductance (mol m-2 s-1 Pa-1)
, conductance_b(land_pnts)                                                     &
                            ! Conductance parameter b (Pa)
, conductance_c(land_pnts)
                            ! Conductance parameter c

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_k(land_pnts)
                            ! Leaf conductance for each land point
                            ! (mol m-2 s-1 Pa-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER ::                                                 &
  RoutineName='LEAF_CONDUCTANCE_IMPAIRED_MEMORY_JLS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL leaf_conductance_jls( pft, land_pnts, water_potential, kmax_ref, kcrit,   &
                           conductance_b, conductance_c, leaf_k )
leaf_k = MIN(leaf_k, MAX(kcap, kcrit))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE leaf_conductance_impaired_memory_jls


! ---------------------------------------------------------------------
! For each candidate transpiration rate E of the stomatal optimiser, work out
! the leaf water potential (Pa) and the leaf's conductance
! (mol m-2 s-1 Pa-1) of a plant damaged by drought, under the memory
! model of xylem impairment, at points where the stomata are open.
!
! The flow on the capped curve is written out exactly (see the equation
! for E above). From the root-zone water potential psi_r down to psi_s
! the conductance is the cap, so where E is small enough the leaf water
! potential is
!
!   psi_l = psi_r - E / k_cap.
!
! For larger E the rest of the path is on the undamaged curve, and psi_l
! is found by Newton's method, each E on its own. With
! som_psi_solver = 3 the capped curve is inverted from the PFT's supply
! table (a precomputed table of flow against water potential) instead.
!
! The conductance returned is MIN(k_intact(psi_l), k_cap), no lower than
! kcrit.
!
! MDK, 25/10/26
! ---------------------------------------------------------------------
SUBROUTINE psi_leaf_impaired_memory( pft,                                      &
                                     n_e_leaf,                                 &
                                     land_pnts,                                &
                                     open_pnts,                                &
                                     veg_index,                                &
                                     open_index,                               &
                                     e_leaf,                                   &
                                     root_zone_psi,                            &
                                     kmax_ref,                                 &
                                     kcap,                                     &
                                     kcrit,                                    &
                                     conductance_b,                            &
                                     conductance_c,                            &
                                   ! INTENT OUT
                                     leaf_psi,                                 &
                                     leaf_k                                    &
  )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut
USE supply_table_mod, ONLY: psi_leaf_capped_lut_jls

INTEGER, INTENT(IN) ::                                                         &
  pft, n_e_leaf, land_pnts, open_pnts, veg_index(land_pnts),                   &
  open_index(land_pnts)
                            ! See psi_leaf_impaired_jls.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  e_leaf(n_e_leaf, open_pnts)                                                  &
                            ! Transpiration rates for each open point
                            ! (mol H2O m-2 leaf s-1)
, root_zone_psi(land_pnts)                                                     &
                            ! Water potential in the root zone (Pa)
, kmax_ref(land_pnts), kcap(land_pnts), kcrit(land_pnts)                       &
, conductance_b(land_pnts), conductance_c(land_pnts)
                            ! See leaf_conductance_impaired_memory_jls.

REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  leaf_psi(n_e_leaf, open_pnts)                                                &
                            ! Leaf water potential for each open point (Pa)
, leaf_k(n_e_leaf, open_pnts)
                            ! Leaf conductance for each open point
                            ! (mol m-2 s-1 Pa-1)

! Local variables
INTEGER :: i, j, l

REAL(KIND=real_jlslsm) :: kmx, kc, bb, cc, psi_r, psi_s, psi_top, psi_floor
                            ! Per-point curve, cap, root zone psi, psi where
                            ! the intact curve meets the cap, top of the
                            ! intact-curve part of the path, and lower
                            ! bound for the Newton iteration.
REAL(KIND=real_jlslsm) :: e_capped
                            ! Transpiration carried by the capped part of the
                            ! path when psi_l < psi_s.
REAL(KIND=real_jlslsm) :: a_top(1)
                            ! A(psi_top).
REAL(KIND=real_jlslsm) :: psi_vec(n_e_leaf), e_vec(n_e_leaf),                  &
                          k_vec(n_e_leaf), dpsi_vec(n_e_leaf)
LOGICAL :: l_capped(n_e_leaf)
                            ! Samples solved within the capped part.
LOGICAL :: l_nr_active(n_e_leaf)
                            ! Newton: samples still iterating (each stops
                            ! on its own convergence, so the answer for one
                            ! E does not depend on the rest of the call).

INTEGER, PARAMETER :: max_nr_iter = 20
REAL(KIND=real_jlslsm), PARAMETER :: psi_tol = 1.0
                            ! Newton-Raphson convergence tolerance (Pa).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PSI_LEAF_IMPAIRED_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! som_psi_solver = 3: the same capped curve inverted exactly from the PFT's
! supply table (psi_leaf_capped_lut_jls). Every point is on the PFT curve
! (b, c), which the memory model does not change.
IF ( som_psi_solver == psi_solver_lut ) THEN
  CALL psi_leaf_capped_lut_jls( pft, n_e_leaf, land_pnts, open_pnts,           &
                                veg_index, open_index, e_leaf,                 &
                                root_zone_psi, kmax_ref, kcap, kcrit,          &
                                leaf_psi, leaf_k )
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

DO j = 1, open_pnts
  l = veg_index(open_index(j))

  kmx   = kmax_ref(l)
  kc    = MAX(MIN(kcap(l), kmx), kcrit(l))
  bb    = conductance_b(l)
  cc    = conductance_c(l)
  psi_r = root_zone_psi(l)
  psi_floor = psi_r - 5.0 * ABS(bb)

  ! Water potential where the intact curve meets the cap (+HUGE when there
  ! is no cap, so the whole path is on the intact curve).
  IF (kc < kmx * (1.0 - 1.0e-6)) THEN
    psi_s = psi_at_k(pft, kc / kmx, bb, cc)
  ELSE
    psi_s = HUGE(1.0_real_jlslsm)
  END IF

  !-------------------------------------------------------------------------
  ! Capped part: constant conductance kc from psi_r down to psi_s.
  !-------------------------------------------------------------------------
  IF (psi_r > psi_s) THEN
    psi_vec(:) = psi_r - e_leaf(:,j) / kc
    l_capped(:) = psi_vec(:) >= psi_s
    e_capped = kc * (psi_r - psi_s)
  ELSE
    l_capped(:) = .FALSE.
    e_capped = 0.0
  END IF

  !-------------------------------------------------------------------------
  ! Intact part: solve e = e_capped + A(psi_l) - A(psi_top) by Newton-
  ! Raphson (dE/dpsi_l = -k_intact(psi_l)). E is concave and decreasing in
  ! psi_l, so starting from psi_top the iterates approach the root
  ! monotonically from above.
  !-------------------------------------------------------------------------
  IF (.NOT. ALL(l_capped)) THEN
    psi_top = MIN(psi_r, psi_s)
    a_top = antideriv(pft, 1, (/ psi_top /), kmx, bb, cc)

    WHERE (.NOT. l_capped) psi_vec(:) = psi_top

    l_nr_active(:) = .NOT. l_capped(:)
    DO i = 1, max_nr_iter
      e_vec(:) = e_capped + antideriv(pft, n_e_leaf, psi_vec, kmx, bb, cc)     &
                 - a_top(1)
      k_vec(:) = MAX(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc),            &
                     TINY(1.0_real_jlslsm))
      dpsi_vec(:) = (e_vec(:) - e_leaf(:,j)) / k_vec(:)
      WHERE (l_capped) dpsi_vec(:) = 0.0

      ! Only the samples on the intact part of the path are iterated (and
      ! bounded above by psi_top): the capped samples were solved exactly
      ! above and lie above psi_top.
      ! Each sample stops when its own step is below psi_tol.
      WHERE (l_nr_active)                                                      &
        psi_vec(:) = MIN(MAX(psi_vec(:) + dpsi_vec(:), psi_floor), psi_top)
      WHERE (ABS(dpsi_vec(:)) < psi_tol) l_nr_active(:) = .FALSE.
      IF ( .NOT. ANY(l_nr_active(:)) ) EXIT
    END DO
  END IF

  leaf_psi(:,j) = psi_vec(:)
  leaf_k(:,j) = MAX(MIN(k_intact(pft, n_e_leaf, psi_vec, kmx, bb, cc), kc),    &
                    kcrit(l))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE psi_leaf_impaired_memory

! ---------------------------------------------------------------------
! Xylem impairment, memory model: update the conductance cap k_cap (mol
! m-2 s-1 Pa-1, leaf area basis) at the end of a step, recovery first and
! then damage.
!
! 1. Cover gained since the last update is undamaged, so the stem's loss
!    of conductance is diluted by the factor ximpair_dilute_area found.
! 2. Recovery of the stem by new wood (l_ximpair_rec_nsc). TRIFFID
!    renews a share rho of the stem's conductance with new rings and
!    leaves it waiting for this update. New wood replaces old xylem at
!    random, so the stem's loss L (0-1) becomes
!
!      L (1 - rho),
!
!    or with l_nsc_ring_additive the new wood replaces damaged xylem
!    first, as in the MEDFATE model, and L becomes MAX(L - rho, 0). With
!    the annual rings (l_nsc_sapwood_rings) nothing waits: the stem's
!    loss is the rings' weighted mean (ximpair_ring_loss_mean).
! 3. Recovery of the leaf segment (l_ximpair_seg_memory). The leaf's
!    loss falls as leaves are replaced,
!
!      L_leaf -> L_leaf EXP(-dt / tau),
!
!    with dt the step length and tau = ximpair_tau_leaf (years). With
!    l_ximpair_rec_leaf_turn the leaf recovers through new leaves in
!    TRIFFID and phenology instead.
! 4. Damage, only where the stomata are open (the xylem is under
!    tension). The cap is lowered to the undamaged conductance at the
!    water potential psi_x that drives the damage,
!
!      k_cap = MIN(k_cap, k_intact(psi_x)),
!
!    where psi_x is chosen by ximpair_psi_driver: 1 leaf, 2 mean of leaf
!    and root zone, 3 root zone (close to predawn, as Mackay et al. 2015
!    use), 4 stem (the leaf end of the stem segment). With the annual
!    rings, the wood of every ring that conducts takes at least the loss
!    this gives (ring_hit).
!
! With l_ximpair_seg_memory (with l_som_plant_segments), k_cap is the
! stem segment's cap and the leaf segment has its own cap. Each is
! damaged by the water potential at its leaf end on its own curve: the
! stem by psi_stem and the leaf by psi_leaf (ximpair_psi_driver is not
! used). The leaf's damage then fades within a few years, while the stem
! keeps its damage until new rings renew it. With the leaf classes of the
! two-leaf canopy (psi_leaf_cls and w_cls present) the leaf's cap is
! lowered to their leaf-area-weighted mean relative conductance,
!
!   k_leaf = sum over c of w_c f_leaf(psi_leaf_c),
!
! with w_c the class's share of leaf area and f_leaf the leaf segment's
! intact curve, so the sunlit and shaded leaves each lose conductance at
! their own water potential. The stem is damaged by psi_stem as before.
!
! With l_ximpair_frozen_hold (the default), a point whose root zone is
! frozen (ximpair_root_frozen > ximpair_frozen_frac) takes no damage on
! the step:
! its soil water potential then comes from the unfrozen water alone and
! says nothing about tension in the dormant xylem. Recovery (steps 1-4)
! goes on.
!
! k_cap is kept between kcrit and kmax. The stem's loss, 1 - k_cap / kmax,
! is stored for TRIFFID and, with l_trif_hydmort, for hydraulic
! mortality.
!
! MDK, 12/11/26
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_memory( n_land_pts, n_open_pts,             &
                                           open_index, pft, psi_leaf,          &
                                           psi_root, psi_stem, kcap,           &
                                           psi_leaf_cls, w_cls )

USE pftparm, ONLY: kmax_pft, kcrit, conductance_b_pft, conductance_c_pft,      &
                   ximpair_psi_driver,                                         &
                   ximpair_tau_leaf, seg_rel_k
USE jules_vegetation_mod, ONLY: ximpair_driver_leaf, ximpair_driver_mean,      &
                                ximpair_driver_root, ximpair_driver_stem,      &
                                l_ximpair_seg_memory, l_ximpair_rec_nsc,       &
                                l_ximpair_rec_leaf_turn, l_trif_hydmort,       &
                                l_som_seg_whole
USE jules_vegetation_mod, ONLY: l_nsc_ring_additive
USE jules_vegetation_mod, ONLY: l_nsc_sapwood_rings, nsc_n_ring
USE jules_vegetation_mod, ONLY: l_ximpair_frozen_hold, ximpair_frozen_frac
USE xylem_impairment_seg_state_mod, ONLY: ximpair_kcap_leaf,                   &
                                          ximpair_psi_stem,                    &
                                          ximpair_seg_state_alloc
USE triffid_hydmort_mod, ONLY: hydmort_accumulate
USE jules_surface_types_mod, ONLY: npft
USE model_time_mod, ONLY: current_time, timestep_len, is_spinup
USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls
USE logging_mod, ONLY: log_warn

USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points, of points with open stomata, the land point of
! each, and the plant functional type.
INTEGER, INTENT(IN) :: n_land_pts, n_open_pts, open_index(n_open_pts), pft

! Leaf and root-zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_leaf(n_land_pts), psi_root(n_land_pts)
! Stem water potential (Pa; leaf-side end of the stem segment, = psi_root
! without segments).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_stem(n_land_pts)

! Conductance cap, kmax_pft * (1 - PLC_mem) (leaf basis).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: kcap(n_land_pts)

REAL(KIND=real_jlslsm), INTENT(IN), OPTIONAL ::                                &
  ! Leaf water potential of the sunlit (1) and shaded (2) leaves of the
  ! two-leaf canopy (Pa).
  psi_leaf_cls(n_land_pts,2),                                                  &
  ! Their shares of the leaf area (sum 1).
  w_cls(n_land_pts,2)

! Point counters, an error code and a ring counter.
INTEGER :: j, l, errcode, i_r
! The leaf classes are given.
LOGICAL :: l_cls
! Damage is held at points with a mostly frozen root zone.
LOGICAL :: l_frozen
! Relative intact conductance of the leaf segment at this step's leaf
! water potential.
REAL(KIND=real_jlslsm) :: f_leaf
! Queued ring share at one point.
REAL(KIND=real_jlslsm) :: rep_l

! Damage-driving water potential (Pa), the intact conductance there, and
! the PFT's intact curve (kmax, kcrit, b in Pa, c) on each point.
REAL(KIND=real_jlslsm) ::                                                      &
  psi_x(n_land_pts), k_x(n_land_pts), kmax_pts(n_land_pts),                    &
  kcrit_pts(n_land_pts), b_pts(n_land_pts), c_pts(n_land_pts)
! This timestep in years, for the leaf segment's recovery.
REAL(KIND=real_jlslsm) :: dt_years
! kcrit / kmax (floor of the leaf segment's cap).
REAL(KIND=real_jlslsm) :: frac_crit

! Seconds in a year.
REAL(KIND=real_jlslsm), PARAMETER :: sec_per_year = 365.25 * 86400.0
! Fraction of kmax - kcrit by which the whole-plant cap is held above
! kcrit.
REAL(KIND=real_jlslsm), PARAMETER :: kcap_crit_margin = 1.0e-3

! The warning that the rings' year was reset has been given.
LOGICAL, SAVE :: l_ring_year_warned = .FALSE.



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_MEMORY'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Damage-driving water potential.
SELECT CASE (ximpair_psi_driver(pft))
CASE (ximpair_driver_leaf)
  psi_x(:) = psi_leaf(:)
CASE (ximpair_driver_mean)
  psi_x(:) = 0.5 * (psi_leaf(:) + psi_root(:))
CASE (ximpair_driver_root)
  psi_x(:) = psi_root(:)
CASE (ximpair_driver_stem)
  psi_x(:) = psi_stem(:)
CASE DEFAULT
  errcode = 101  !  a hard error
  CALL ereport(RoutineName, errcode,                                           &
               'ximpair_psi_driver should be leaf (1), mean (2), root (3) ' // &
               'or stem (4)')
END SELECT

kmax_pts(:)  = kmax_pft(pft)
kcrit_pts(:) = kcrit(pft)
b_pts(:)     = conductance_b_pft(pft)
c_pts(:)     = conductance_c_pft(pft)

! Ring cohorts: set at the first update, all rings with the stem's loss;
! rings from a dump without their year grow from the current year. The
! youngest ring grows this year or next, so outside spin-up a year after
! that (rings from a dump written later than the run's start) would give
! the rings negative ages: their ages are then counted from today.
IF (l_nsc_sapwood_rings) THEN
  CALL ximpair_memory_alloc( n_land_pts )
  DO l = 1, n_land_pts
    IF (ximpair_ring_frac(l,pft,1) < 0.0)                                      &
      CALL ximpair_ring_init(l, pft, MAX(1.0 - kcap(l) / kmax_pts(l), 0.0))
    IF (ximpair_ring_year(l,pft) < 0.0) THEN
      ximpair_ring_year(l,pft) = REAL(current_time%year)
    ELSE IF (.NOT. is_spinup .AND.                                             &
             NINT(ximpair_ring_year(l,pft)) > current_time%year + 1) THEN
      IF (.NOT. l_ring_year_warned)                                            &
        CALL log_warn(RoutineName,                                             &
                      "The sapwood rings' year is after next year " //         &
                      "(a dump from a later date?): the ring ages are " //     &
                      "counted from today.")
      l_ring_year_warned = .TRUE.
      ximpair_ring_year(l,pft) = REAL(ximpair_ring_next_year(pft))
    END IF
  END DO
END IF

! Area gained by TRIFFID since the last update (ximpair_dilute_area): the
! stem's loss of conductance diluted by the new, intact area.
IF ( ALLOCATED(ximpair_area_w) ) THEN
  IF (l_nsc_sapwood_rings) THEN
    ! Every ring of the old cover keeps its loss, diluted alike; the new
    ! cover's wood is intact and joins each ring's less damaged part.
    DO l = 1, n_land_pts
      IF (ximpair_area_w(l,pft) < 1.0 .AND.                                    &
          ximpair_ring_frac(l,pft,1) >= 0.0) THEN
        DO i_r = 1, nsc_n_ring
          CALL ring_add_intact(l, pft, i_r, ximpair_area_w(l,pft))
        END DO
      END IF
    END DO
  ELSE
    WHERE ( ximpair_area_w(:,pft) < 1.0 )                                      &
      kcap(:) = kmax_pts(:) - (kmax_pts(:) - kcap(:)) * ximpair_area_w(:,pft)
  END IF
  ximpair_area_w(:,pft) = 1.0
END IF

!-----------------------------------------------------------------------------
! Recovery: the stem by ring renewal (l_ximpair_rec_nsc), the leaf segment
! by leaf turnover, and the annual reset.
!-----------------------------------------------------------------------------

IF (l_ximpair_rec_nsc) THEN

  CALL ximpair_memory_alloc( n_land_pts )

  ! Ring renewal: repair_pend is the share rho of the stem's conductance
  ! renewed by the new rings since the last update. New and old xylem are
  ! replaced alike (at random), so the whole stem's loss is diluted,
  !
  !   L -> L (1 - rho),
  !
  ! or with l_nsc_ring_additive the new ring replaces damaged xylem first,
  ! as in MEDFATE, as if the sapwood retired to heartwood were the
  ! embolised part (old rings, vessels sealed by tyloses),
  !
  !   L -> MAX(L - rho, 0).
  !
  ! Points with nothing queued are left untouched.
  !
  ! With the ring cohorts the rings carry the loss instead (TRIFFID adds
  ! the new wood to the youngest ring and retires the oldest), so the cap
  ! follows them and nothing is queued.
  IF (l_nsc_sapwood_rings) THEN
    DO l = 1, n_land_pts
      kcap(l) = kmax_pts(l) * (1.0 - ximpair_ring_loss_mean(l, pft))
    END DO
  ELSE
    DO l = 1, n_land_pts
      rep_l = ximpair_repair_pend(l,pft)
      IF (rep_l <= 0.0) CYCLE
      IF (l_nsc_ring_additive) THEN
        kcap(l) = kmax_pts(l) * (1.0 - MAX(1.0 - kcap(l) / kmax_pts(l)        &
                                           - rep_l, 0.0))
      ELSE
        kcap(l) = kmax_pts(l) * (1.0 - MAX(1.0 - kcap(l) / kmax_pts(l), 0.0)  &
                                       * (1.0 - MIN(rep_l, 1.0)))
      END IF
    END DO
  END IF
  ximpair_repair_pend(:,pft) = 0.0

END IF

! This timestep in years, for the leaf segment's recovery.
dt_years = REAL(timestep_len) / sec_per_year

! With the segments on the PFT curve (l_som_seg_whole) the stem psi is
! kept for the psi_stem_pft output only (the leaf state is not used).
IF (l_som_seg_whole) THEN
  IF (.NOT. ALLOCATED(ximpair_psi_stem)) THEN
    ALLOCATE(ximpair_psi_stem(n_land_pts, npft))
    ximpair_psi_stem(:,:) = 0.0
  END IF
  ximpair_psi_stem(:,pft) = psi_stem(:)
END IF

IF (l_ximpair_seg_memory) THEN
  CALL ximpair_seg_state_alloc( n_land_pts )
  ximpair_psi_stem(:,pft) = psi_stem(:)
  ! With leaf renewal (l_ximpair_rec_leaf_turn) the leaf segment recovers
  ! by new leaves instead (triffid_nsc_turgor, phenol).
  IF (ximpair_tau_leaf(pft) > 0.0 .AND. .NOT. l_ximpair_rec_leaf_turn) THEN
    ximpair_kcap_leaf(:,pft) = 1.0 - (1.0 - ximpair_kcap_leaf(:,pft))          &
                               * EXP(-dt_years / ximpair_tau_leaf(pft))
  END IF
END IF

!-----------------------------------------------------------------------------
! Damage: the cap can not exceed the intact conductance at the damage
! driver. Only points with open stomata (i.e. under tension) are damaged,
! and with l_ximpair_frozen_hold not those whose root zone is mostly
! frozen.
!-----------------------------------------------------------------------------
l_frozen = l_ximpair_frozen_hold .AND. ALLOCATED(ximpair_root_frozen)

IF (l_ximpair_seg_memory) THEN
  ! Per segment: the stem by psi_stem on the stem curve, the leaf by
  ! psi_leaf on the leaf curve (the segments' intact Weibull curves), or
  ! with the leaf classes by the weighted mean of their conductances.
  l_cls = PRESENT(psi_leaf_cls) .AND. PRESENT(w_cls)
  DO j = 1, n_open_pts
    l = open_index(j)

    IF (l_frozen) THEN
      IF (ximpair_root_frozen(l) > ximpair_frozen_frac) CYCLE
    END IF

    IF (l_nsc_sapwood_rings) THEN
      ! Every ring takes at least the imposed loss.
      CALL ring_damage(l, 1.0 - seg_rel_k(pft, 2, psi_stem(l)))
      kcap(l) = kmax_pts(l) * (1.0 - ximpair_ring_loss_mean(l, pft))
    ELSE
      kcap(l) = MIN(kcap(l), kmax_pts(l) * seg_rel_k(pft, 2, psi_stem(l)))
    END IF
    IF (l_cls) THEN
      f_leaf = w_cls(l,1) * seg_rel_k(pft, 3, psi_leaf_cls(l,1))               &
               + w_cls(l,2) * seg_rel_k(pft, 3, psi_leaf_cls(l,2))
    ELSE
      f_leaf = seg_rel_k(pft, 3, psi_leaf(l))
    END IF
    ximpair_kcap_leaf(l,pft) = MIN(ximpair_kcap_leaf(l,pft), f_leaf)
  END DO
  frac_crit = kcrit(pft) / MAX(kmax_pft(pft), TINY(1.0_real_jlslsm))
  ximpair_kcap_leaf(:,pft) = MAX(MIN(ximpair_kcap_leaf(:,pft), 1.0), frac_crit)
ELSE
  CALL leaf_conductance_jls( pft, n_land_pts, psi_x, kmax_pts, kcrit_pts,      &
                             b_pts, c_pts, k_x )
  DO j = 1, n_open_pts
    l = open_index(j)

    IF (l_frozen) THEN
      IF (ximpair_root_frozen(l) > ximpair_frozen_frac) CYCLE
    END IF

    IF (l_nsc_sapwood_rings) THEN
      ! Every ring takes at least the imposed loss.
      CALL ring_damage(l, 1.0 - k_x(l) / kmax_pts(l))
      kcap(l) = kmax_pts(l) * (1.0 - ximpair_ring_loss_mean(l, pft))
    ELSE
      kcap(l) = MIN(kcap(l), k_x(l))
    END IF
  END DO
END IF

kcap(:) = MAX(MIN(kcap(:), kmax_pts(:)), kcrit_pts(:))

! On the whole-plant curve the cap is held just above kcrit: at kcrit
! every transpiration, E = 0 included, gives k = kcrit, which fails the
! optimiser's k > kcrit test, so the plant would stay shut until
! recovery. (With per-segment memory the root keeps the path above kcrit.)
IF (.NOT. l_ximpair_seg_memory)                                                &
  kcap(:) = MAX(kcap(:), kcrit_pts(:) + kcap_crit_margin                       &
                                        * (kmax_pts(:) - kcrit_pts(:)))

IF (.NOT. ALLOCATED(ximpair_lock)) THEN
  ALLOCATE(ximpair_lock(n_land_pts, npft))
  ximpair_lock(:,:) = 0.0
END IF
ximpair_lock(:,pft) = 1.0 - kcap(:) / kmax_pts(:)
! Worst damage of the TRIFFID period, for hydraulic mortality.
IF (l_trif_hydmort) CALL hydmort_accumulate( n_land_pts, pft,                  &
                                             ximpair_lock(:,pft) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS

! ---------------------------------------------------------------------
! Xylem impairment, memory model: apply this step's damage to the annual
! rings at point l. Both parts of each ring (ring_hit) take at least the
! loss
!
!   MIN(loss_imp, 1 - kcrit / kmax),
!
! where loss_imp (0-1) is the loss this step's water potential gives and
! the second term is the loss at the lowest conductance allowed. The
! highest such loss since the last TRIFFID step is kept
! (ximpair_ring_peak).
!
! A youngest ring that does not conduct yet (not grown, or a ring-porous
! tree's new ring before its vessels conduct) is not damaged.
!
! MDK, 05/11/26
! ---------------------------------------------------------------------
SUBROUTINE ring_damage( l, loss_imp )

! Land point.
INTEGER, INTENT(IN) :: l
! The loss of conductance the step imposes (0-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: loss_imp

! The loss at kcrit.
REAL(KIND=real_jlslsm) :: loss_max
! The first ring damaged, and a ring counter.
INTEGER :: i_first, i_r

loss_max = 1.0 - kcrit_pts(l) / kmax_pts(l)
ximpair_ring_peak(l,pft) = MAX(ximpair_ring_peak(l,pft),                       &
                               MIN(loss_imp, loss_max))
i_first = 1
IF (ring_age(l, pft, 1) < 0 .OR.                                               &
    ring_cond(pft, ring_age(l, pft, 1), 0) <= 0.0 .AND.                        &
    ring_age(l, pft, 1) == 0) i_first = 2
DO i_r = i_first, nsc_n_ring
  CALL ring_hit(l, pft, i_r, MIN(loss_imp, loss_max))
END DO

END SUBROUTINE ring_damage

END SUBROUTINE update_xylem_impairment_memory

END MODULE xylem_impairment_memory_mod
