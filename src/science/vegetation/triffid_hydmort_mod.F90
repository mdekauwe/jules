! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE triffid_hydmort_mod
!-----------------------------------------------------------------------------
! Hydraulic (drought) mortality in TRIFFID (l_trif_hydmort): the share of a
! PFT's cover that dies when drought damage to its xylem passes what its
! trees can survive.
!
! A grid box, or a site's footprint, holds a stand of trees, but the model
! simulates one average tree. Its lasting loss of xylem conductance, p
! (0-1; ximpair_lock from the memory impairment model, which holds the
! loss until new xylem is built), is taken as the damage of the whole
! stand. Trees differ in the damage they survive: their lethal thresholds
! follow a logistic spread with centre mu = hydmort_mu and width
! w = hydmort_width (both as fractions of conductance lost). The fraction
! of the stand that has passed its own threshold at damage p is
!
!   S(p) = [s((p - mu)/w) - s((b - mu)/w)] / [1 - s((b - mu)/w)],  p > b,
!   S(p) = 0,                                                      p <= b,
!
! with s(z) = 1 / (1 + EXP(-z)) the logistic function and b =
! hydmort_plc_bg the background loss of a normal year (native embolism),
! which kills no tree. Without b the tail of a wide spread would kill a
! few trees every year. S is an expected fraction; there are no random
! draws. A loss of mu kills a little under half the stand, so a large loss
! does not mean the whole stand is dead.
!
! Each tree dies once, so the dead fraction can only grow. At each TRIFFID
! step, with p_peak the worst damage since the last step and p_ref the
! worst damage the stand had survived before it, the fraction of the
! surviving cover that dies in the step is
!
!   m = MAX(0, S(p_peak) - S(p_ref)) / (1 - S(p_ref)),
!
! and the matching area disturbance is
!
!   g = r_gamma m   (per 360 days),
!
! with r_gamma the inverse TRIFFID step (per 360 days). p_ref is a running
! maximum. The trees whose thresholds lie below it are already dead, so a
! later drought kills only if it takes the damage past p_ref, and p_ref
! does not fall when the xylem is repaired, which would kill the same
! trees twice. A lasting damage that does not rise kills no more trees.
! m below hydmort_m_min is set to 0, so the logistic tail at small damage
! does not make a run differ from one with mortality off.
!
! The dead cover goes to litter through TRIFFID's litter calculation,
! which turns the carbon of a shrinking area into litter:
!   - with competition (l_veg_compete), g is added to the area disturbance
!     in lotka, with the standard disturbance g_area and the store's
!     starvation term, and competition refills the space;
!   - with fixed fractions (l_veg_compete = F), the dead fraction becomes
!     bare soil (triffid_hydmort_cover) and is not refilled.
!
! p_ref, the fraction of the stand dead so far and the worst damage since
! the last TRIFFID step are kept in dumps. A dump without them starts
! p_ref at the current damage, as if the stand has already survived it.
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
! Allocates the mortality state and diagnostics for all land points and
! PFTs, if not already allocated. The worst damage and p_ref start at -1,
! which marks them as not yet set.
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
! Return the fraction of the stand that has died from hydraulic damage p
! (0-1), S(p) in the module header, which triffid_hydmort uses to set
! the mortality: 0 at or below the background loss b, and above it rising
! along the logistic curve with centre mu and width w, towards 1 as p grows.
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
! Keeps the worst damage since the last TRIFFID step for PFT pft: the
! running maximum of the lasting loss of conductance (ximpair_lock, 0-1).
! Called every time step from the memory impairment update, after
! ximpair_lock is set.
!
! MDK, 24/10/26
! ---------------------------------------------------------------------
SUBROUTINE hydmort_accumulate( n_land_pts, pft, lock )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points; plant functional type.
INTEGER, INTENT(IN) :: n_land_pts, pft
! Lasting loss of conductance (ximpair_lock).
REAL(KIND=real_jlslsm), INTENT(IN) :: lock(n_land_pts)

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
! Kill the part of each natural PFT's stand that drought damage to its
! xylem has made unable to survive. At each TRIFFID step it works out the
! fraction m of the surviving cover that dies, and the matching area
! disturbance g = r_gamma m (per 360 days), from the worst damage since
! the last step (equations in the module header). The worst damage is
! the larger of the running peak and the current lasting loss. Where
! p_ref is not yet set it starts at the current loss. After the step
! p_ref rises to the worst damage if that is higher, the running peak
! restarts from the current loss, and the fraction of the stand dead
! since the start of the run becomes 1 - (1 - dead)(1 - m). Crop PFTs
! take no hydraulic mortality.
!
! MDK, 25/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_hydmort( land_pts, trif_pts, trif_index, r_gamma, lock )

USE jules_surface_types_mod, ONLY: npft, nnpft
USE pftparm, ONLY: hydmort_mu, hydmort_width, hydmort_plc_bg

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points, and of points on which TRIFFID may operate.
INTEGER, INTENT(IN) :: land_pts, trif_pts
! Land points on which TRIFFID may operate.
INTEGER, INTENT(IN) :: trif_index(land_pts)
! Inverse TRIFFID timestep (/360 days).
REAL(KIND=real_jlslsm), INTENT(IN) :: r_gamma
! Lasting loss of conductance (ximpair_lock).
REAL(KIND=real_jlslsm), INTENT(IN) :: lock(land_pts,npft)

INTEGER :: l, n, t
! Peak damage of the period, the fractions of the stand dead at it and at p_ref,
! and the fraction of the surviving cover killed in this step.
REAL(KIND=real_jlslsm) :: p_peak, s_peak, s_ref, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_HYDMORT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL hydmort_alloc( land_pts )

! Natural PFTs only (lotka): crop PFTs (n > nnpft) take no hydraulic
! mortality.
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
! Turn the cover killed by drought into bare soil when the vegetation
! fractions are fixed (l_veg_compete = F), where TRIFFID competition
! cannot do it. It moves the cover killed in this TRIFFID step, m times
! the PFT's fraction, from each natural PFT to bare soil, where it is not
! refilled. The loss is also taken off the PFT's
! change in fraction over the step (dfrac, which becomes more negative),
! so TRIFFID's litter calculation sends the dead area's vegetation carbon
! and nitrogen to the soil.
!
! MDK, 25/10/26
! ---------------------------------------------------------------------
SUBROUTINE triffid_hydmort_cover( land_pts, trif_pts, trif_index, frac,        &
                                  frac_na, frac_nofire, dfrac, dfrac_na,       &
                                  dfrac_nofire )

USE jules_surface_types_mod, ONLY: npft, nnpft, ntype, soil

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

! Number of land points, and of points on which TRIFFID may operate.
INTEGER, INTENT(IN) :: land_pts, trif_pts
! Land points on which TRIFFID may operate.
INTEGER, INTENT(IN) :: trif_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac(land_pts,ntype)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac_na(land_pts,ntype)
! Areal fractions, as in triffid.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: frac_nofire(land_pts,ntype)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac(land_pts,npft)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac_na(land_pts,npft)
! Their changes over the TRIFFID step.
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: dfrac_nofire(land_pts,npft)

INTEGER :: l, n, t
! Fraction of the grid box killed in this step.
REAL(KIND=real_jlslsm) :: dlost

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='TRIFFID_HYDMORT_COVER'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO n = 1, nnpft   ! natural PFTs only
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
