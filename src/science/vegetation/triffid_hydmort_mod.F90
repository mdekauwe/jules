! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE triffid_hydmort_mod
!-----------------------------------------------------------------------------
! Hydraulic (drought) mortality in TRIFFID (l_trif_hydmort).
!
! A grid box (or a site's footprint) is a stand of trees, but the model
! simulates one average tree. Its lasting loss of xylem conductance,
! p = ximpair_lock (memory impairment model: held until new xylem is built),
! is read as the damage of the stand. The trees differ in the damage they
! can survive, and these lethal thresholds are spread logistically about
! hydmort_mu with width hydmort_width, so the fraction of the stand that has
! passed its own threshold at damage p is
!   S(p) = [s((p - mu)/w) - s((b - mu)/w)] / [1 - s((b - mu)/w)],  p > b,
!   S(p) = 0,  p <= b,          s(z) = 1/(1+e^-z),
! with b = hydmort_plc_bg the background loss of a normal year (native
! embolism), which kills no tree: without it the logistic tail of a wide
! spread kills a little every year. This is an expected value; there are no
! random draws. 88 % loss is then not "all dead": with mu = 0.88 half the
! stand is.
!
! Each tree dies once, so mortality is a ratchet. At each TRIFFID step,
! with p_peak the worst damage of the period and p_ref the damage the
! surviving stand has already withstood,
!   m = max(0, S(p_peak) - S(p_ref)) / (1 - S(p_ref))
! is the fraction of the surviving cover that dies in the step, and
!   g = r_gamma m   (per 360 days)
! is the matching area disturbance. p_ref is the worst damage the stand has
! been through, a running maximum: each tree has a fixed threshold, so the
! trees whose thresholds lie below p_ref are already dead, and a later
! drought kills only if it takes the damage past p_ref (it does not fall
! back when the xylem is repaired, which would kill the same threshold band
! twice). A lasting damage that does not rise kills no more. Kills below
! hydmort_m_min are ignored, so the logistic tail at trivial damage does not
! move a run away from one with the switch off.
!
! The killed cover goes to litter through TRIFFID's mass-balance litter:
! - with l_veg_compete, g is added to the disturbance in lotka (alongside
!   g_area and the store's starvation term), and the space is refilled by
!   competition;
! - with fixed fractions (l_veg_compete = F, e.g. monoculture runs), the
!   killed fraction is moved to bare soil here, and is not refilled.
!
! The design notes are local, not in the repository.
! p_ref, the stand killed so far and the period's peak damage are in dumps.
! A dump written without them starts p_ref at the current damage (the stand
! is taken to have withstood it).
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: hydmort_plc_peak, hydmort_pref, hydmort_g_pft, hydmort_m_pft,        &
          hydmort_dead_pft, hydmort_alloc, hydmort_accumulate,                 &
          triffid_hydmort, triffid_hydmort_cover, hydmort_s

! Worst damage (loss of conductance) since the last TRIFFID step
! (< 0 = none seen yet: the damage at the step is used).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: hydmort_plc_peak(:,:)
! Damage the surviving stand has withstood (< 0 = not set: the damage at
! the first TRIFFID step).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: hydmort_pref(:,:)
! Diagnostics of the last TRIFFID step: the area disturbance (per 360 days)
! and the fraction of the PFT's cover killed in the step.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: hydmort_g_pft(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: hydmort_m_pft(:,:)
! Fraction of the stand killed since the start of the run, 1 - prod(1 - m).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: hydmort_dead_pft(:,:)

! Smallest fraction of the cover killed in a step that is applied.
REAL(KIND=real_jlslsm), PARAMETER :: hydmort_m_min = 1.0e-6

CHARACTER(LEN=*), PARAMETER :: ModuleName='TRIFFID_HYDMORT_MOD'

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state and diagnostics (n_land_pts, npft), if not done.
! ---------------------------------------------------------------------
SUBROUTINE hydmort_alloc( n_land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDMORT_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(hydmort_plc_peak)) THEN
  ALLOCATE(hydmort_plc_peak(n_land_pts, npft))
  hydmort_plc_peak(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(hydmort_pref)) THEN
  ALLOCATE(hydmort_pref(n_land_pts, npft))
  hydmort_pref(:,:) = -1.0
END IF
IF (.NOT. ALLOCATED(hydmort_g_pft)) THEN
  ALLOCATE(hydmort_g_pft(n_land_pts, npft))
  hydmort_g_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(hydmort_m_pft)) THEN
  ALLOCATE(hydmort_m_pft(n_land_pts, npft))
  hydmort_m_pft(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(hydmort_dead_pft)) THEN
  ALLOCATE(hydmort_dead_pft(n_land_pts, npft))
  hydmort_dead_pft(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE hydmort_alloc

! ---------------------------------------------------------------------
! Fraction of the stand dead at damage p: the logistic S(p), from 0 at the
! background loss b.
! ---------------------------------------------------------------------
FUNCTION hydmort_s( p, mu, w, b ) RESULT( s )

REAL(KIND=real_jlslsm), INTENT(IN) :: p, mu, w, b
REAL(KIND=real_jlslsm) :: s
REAL(KIND=real_jlslsm) :: s0

s0 = 1.0 / (1.0 + EXP(-(b - mu) / w))
s  = 1.0 / (1.0 + EXP(-(MAX(p, b) - mu) / w))
s  = MAX(0.0, (s - s0) / (1.0 - s0))

END FUNCTION hydmort_s

! ---------------------------------------------------------------------
! Keep the worst damage of the TRIFFID period (called from the memory
! impairment update, every timestep, after ximpair_lock is set).
! ---------------------------------------------------------------------
SUBROUTINE hydmort_accumulate( n_land_pts, pft, lock )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts, pft
                          ! Number of land points; plant functional type.
                          !
                          ! MDK, 24/10/26
REAL(KIND=real_jlslsm), INTENT(IN) :: lock(n_land_pts)
                          ! Lasting loss of conductance (ximpair_lock).
                          !
                          ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='HYDMORT_ACCUMULATE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL hydmort_alloc( n_land_pts )
hydmort_plc_peak(:,pft) = MAX(hydmort_plc_peak(:,pft), lock(:))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE hydmort_accumulate

! ---------------------------------------------------------------------
! At a TRIFFID step: the fraction of each PFT's cover killed (m) and the
! matching disturbance (g = r_gamma m), from the period's peak damage.
! ---------------------------------------------------------------------
SUBROUTINE triffid_hydmort( land_pts, trif_pts, trif_index, r_gamma, lock )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE pftparm, ONLY: hydmort_mu, hydmort_width, hydmort_plc_bg

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: land_pts, trif_pts
                          ! Number of land points, and of points on which
                          ! TRIFFID may operate.
                          !
                          ! MDK, 24/10/26
INTEGER, INTENT(IN) :: trif_index(land_pts)
                          ! Land points on which TRIFFID may operate.
                          !
                          ! MDK, 24/10/26
REAL(KIND=real_jlslsm), INTENT(IN) :: r_gamma
                          ! Inverse TRIFFID timestep (/360 days).
                          !
                          ! MDK, 24/10/26
REAL(KIND=real_jlslsm), INTENT(IN) :: lock(land_pts,npft)
                          ! Lasting loss of conductance (ximpair_lock).

INTEGER :: l, n, t
REAL(KIND=real_jlslsm) :: p_peak, s_peak, s_ref, m
                          ! Peak damage of the period, the fractions of the
                          ! stand dead at it and at p_ref, and the fraction
                          ! of the surviving cover killed in this step.
                          !
                          ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_HYDMORT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL hydmort_alloc( land_pts )

! Natural PFTs only (lotka): crop PFTs (n > nnpft) take no hydraulic
! mortality. MDK, 25/10/26
DO n = 1, nnpft
  DO t = 1, trif_pts
    l = trif_index(t)
    p_peak = MAX(hydmort_plc_peak(l,n), lock(l,n))
    IF (hydmort_pref(l,n) < 0.0) hydmort_pref(l,n) = lock(l,n)
    s_peak = hydmort_s(p_peak, hydmort_mu(n), hydmort_width(n),                &
                       hydmort_plc_bg(n))
    s_ref  = hydmort_s(hydmort_pref(l,n), hydmort_mu(n), hydmort_width(n),     &
                       hydmort_plc_bg(n))
    m = MAX(0.0, s_peak - s_ref) / MAX(1.0 - s_ref, 1.0e-6)
    m = MIN(m, 1.0)
    IF (m < hydmort_m_min) m = 0.0
    hydmort_m_pft(l,n)    = m
    hydmort_g_pft(l,n)    = r_gamma * m
    hydmort_dead_pft(l,n) = 1.0 - (1.0 - hydmort_dead_pft(l,n)) * (1.0 - m)
    ! The survivors have withstood the worst damage so far.
    hydmort_pref(l,n)     = MAX(hydmort_pref(l,n), p_peak)
    hydmort_plc_peak(l,n) = lock(l,n)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_hydmort

! ---------------------------------------------------------------------
! Fixed fractions (l_veg_compete = F): the killed cover becomes bare soil.
! dfrac (< 0) carries the loss, so TRIFFID's mass-balance litter sends the
! killed area's vegetation carbon and N to the soil.
! ---------------------------------------------------------------------
SUBROUTINE triffid_hydmort_cover( land_pts, trif_pts, trif_index, frac,        &
                                  frac_na, frac_nofire, dfrac, dfrac_na,       &
                                  dfrac_nofire )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype, soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: land_pts, trif_pts
                          ! Number of land points, and of points on which
                          ! TRIFFID may operate.
                          !
                          ! MDK, 24/10/26
INTEGER, INTENT(IN) :: trif_index(land_pts)
                          ! Land points on which TRIFFID may operate.
                          !
                          ! MDK, 24/10/26
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac(land_pts,ntype)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac_na(land_pts,ntype)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac_nofire(land_pts,ntype)
                          ! Areal fractions, as in triffid.
                          !
                          ! MDK, 24/10/26
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac_na(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac_nofire(land_pts,npft)
                          ! Their changes over the TRIFFID step.
                          !
                          ! MDK, 24/10/26

INTEGER :: l, n, t
REAL(KIND=real_jlslsm) :: dlost
                          ! Fraction of the grid box killed in this step.
                          !
                          ! MDK, 24/10/26

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_HYDMORT_COVER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1, nnpft   ! natural PFTs only (MDK, 25/10/26)
  DO t = 1, trif_pts
    l = trif_index(t)
    dlost = hydmort_m_pft(l,n) * frac(l,n)
    IF (dlost <= 0.0) CYCLE
    frac(l,n)          = frac(l,n) - dlost
    frac_na(l,n)       = frac_na(l,n) - dlost
    frac_nofire(l,n)   = frac_nofire(l,n) - dlost
    dfrac(l,n)         = dfrac(l,n) - dlost
    dfrac_na(l,n)      = dfrac_na(l,n) - dlost
    dfrac_nofire(l,n)  = dfrac_nofire(l,n) - dlost
    frac(l,soil)        = frac(l,soil) + dlost
    frac_na(l,soil)     = frac_na(l,soil) + dlost
    frac_nofire(l,soil) = frac_nofire(l,soil) + dlost
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE triffid_hydmort_cover

END MODULE triffid_hydmort_mod
