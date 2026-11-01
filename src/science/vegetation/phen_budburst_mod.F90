! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE phen_budburst_mod

!-----------------------------------------------------------------------------
! Budburst for deciduous TRIFFID phenology, used for a PFT when phen_bb_doy
! or phen_bb_gdd is set: leaves come out once a year, on a set day or
! after enough spring warmth, and the canopy stays bare through winter.
!
! In the standard phenology (phenol) leaves regrow whenever the daily leaf
! turnover falls below twice its minimum, so a deciduous canopy can regrow
! in a mild winter spell, and spring regrowth uses the same temperature
! threshold as autumn leaf fall (tleaf_of, dgl_dt). With budburst on:
!
!   1. Once the canopy has fallen to its floor (phen = 0.01, the fraction
!      of the balanced LAI) after this year's budburst, the PFT is
!      dormant. phenol grows no leaves while it is dormant, nor before
!      this year's budburst.
!   2. Budburst is on day of year phen_bb_doy, or on the day the growing
!      degree-days reach phen_bb_gdd (degC d). The growing degree-days are
!      summed from 1 January over each time step,
!
!        GDD = SUM MAX(T_air - phen_bb_tbase, 0) dt / 86400,
!
!      with T_air the air temperature (degC), phen_bb_tbase the base
!      temperature (degC) and dt the time step (s).
!   3. Autumn leaf fall is phenol's standard cold leaf drop, held off for
!      phen_bb_min_days after budburst. The drop responds to cold nights,
!      so without the hold a cool spell could stop the spring flush, or
!      shed a summer canopy and start dormancy.
!
! The year starts on 1 January, which suits northern-hemisphere seasons. A
! point that never reaches phen_bb_gdd in a year stays leafless that year,
! as a site too cold for the PFT, so set budburst only for deciduous PFTs.
!
! For each land point and PFT the module keeps, in dumps, the growing
! degree-days since 1 January (degC d), the day of year of this year's
! budburst (-1 until it happens) and whether the PFT is dormant (1) or
! not (0). The budburst day is public so other schemes, such as the
! store's leaf flush, can start at the model's budburst.
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: phen_gdd, phen_bb_day, phen_dormant, phen_bb_alloc,                  &
          phen_bb_update, phen_bb_on, phen_bb_blocked, phen_bb_in_leaf

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_gdd(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_bb_day(:,:)
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: phen_dormant(:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PHEN_BUDBURST_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! True when budburst is on for PFT n, that is when phen_bb_doy or
! phen_bb_gdd is set.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
LOGICAL FUNCTION phen_bb_on( n )

USE pftparm, ONLY: phen_bb_doy, phen_bb_gdd
USE pftparm, ONLY: is_set

IMPLICIT NONE

INTEGER, INTENT(IN) :: n

phen_bb_on = is_set(phen_bb_doy(n)) .OR. is_set(phen_bb_gdd(n))

END FUNCTION phen_bb_on

!-----------------------------------------------------------------------------
! Allocates the budburst state, if not already allocated, at its values
! for a fresh start: no degree-days, no budburst yet (-1), not dormant.
!-----------------------------------------------------------------------------
SUBROUTINE phen_bb_alloc( land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='PHEN_BB_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (.NOT. ALLOCATED(phen_gdd)) THEN
  ALLOCATE(phen_gdd(land_pts, npft))
  phen_gdd(:,:) = 0.0
END IF

IF (.NOT. ALLOCATED(phen_bb_day)) THEN
  ALLOCATE(phen_bb_day(land_pts, npft))
  phen_bb_day(:,:) = -1.0
END IF

IF (.NOT. ALLOCATED(phen_dormant)) THEN
  ALLOCATE(phen_dormant(land_pts, npft))
  phen_dormant(:,:) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE phen_bb_alloc

!-----------------------------------------------------------------------------
! Track spring warming for PFT ft and mark the day its leaves come out,
! so phenol knows when the PFT may grow leaves. Called every time step,
! it adds the step's growing degree-days and marks budburst. It does
! nothing for a PFT without budburst.
!
! On the first time step of 1 January the degree-days go back to 0 and
! the budburst day to -1 (not yet). Budburst happens on the first step on
! which the day of year reaches phen_bb_doy or the degree-days reach
! phen_bb_gdd. That day is stored as this year's budburst day, and the PFT
! is no longer dormant.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE phen_bb_update( land_pts, veg_pts, veg_index, ft, timestep, tair )

USE pftparm, ONLY: phen_bb_doy, phen_bb_gdd, phen_bb_tbase
USE pftparm, ONLY: is_set
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap
USE water_constants_mod, ONLY: tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
  ! Number of land points.
  land_pts,                                                                    &
  ! Number of points with this PFT.
  veg_pts,                                                                     &
  ! Index of those points on the land grid.
  veg_index(land_pts),                                                         &
  ! Plant functional type.
  ft

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Model time step (s).
  timestep,                                                                    &
  ! Air temperature (K).
  tair(land_pts)

INTEGER :: j, l
! Day of the year.
REAL(KIND=real_jlslsm) :: doy
! First time step of the year.
LOGICAL :: l_new_year

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PHEN_BB_UPDATE'

IF (.NOT. phen_bb_on(ft)) RETURN
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL phen_bb_alloc( land_pts )

doy = REAL(day_of_year(current_time%year, current_time%month,                  &
                       current_time%day, l_360, l_leap))

! First time step of the year: restart the degree-days and budburst.
l_new_year = current_time%month == 1 .AND. current_time%day == 1 .AND.         &
             REAL(current_time%time) < timestep

DO j = 1,veg_pts
  l = veg_index(j)
  IF (l_new_year) THEN
    phen_gdd(l,ft)    = 0.0
    phen_bb_day(l,ft) = -1.0
  END IF

  IF (is_set(phen_bb_gdd(ft)))                                                 &
    phen_gdd(l,ft) = phen_gdd(l,ft)                                            &
                     + MAX(tair(l) - tm - phen_bb_tbase(ft), 0.0)              &
                       * timestep / 86400.0

  IF (phen_bb_day(l,ft) < 0.0) THEN
    IF ((is_set(phen_bb_doy(ft)) .AND. doy >= phen_bb_doy(ft)) .OR.            &
        (is_set(phen_bb_gdd(ft)) .AND. phen_gdd(l,ft) >= phen_bb_gdd(ft)))     &
    THEN
      phen_bb_day(l,ft)  = doy
      phen_dormant(l,ft) = 0.0
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE phen_bb_update

!-----------------------------------------------------------------------------
! Used by phenol: true when PFT n at land point l may not grow leaves,
! because it is dormant or this year's budburst has not happened yet.
! False for a PFT without budburst.
!-----------------------------------------------------------------------------
LOGICAL FUNCTION phen_bb_blocked( l, n )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, n

phen_bb_blocked = .FALSE.
IF (.NOT. phen_bb_on(n) .OR. .NOT. ALLOCATED(phen_bb_day)) RETURN
phen_bb_blocked = phen_bb_day(l,n) < 0.0 .OR. phen_dormant(l,n) > 0.5

END FUNCTION phen_bb_blocked

!-----------------------------------------------------------------------------
! Used by phenol: true when PFT n at land point l is less than
! phen_bb_min_days after this year's budburst, so the cold leaf drop is
! held off. False without budburst, before this year's budburst, or when
! phen_bb_min_days <= 0.
!-----------------------------------------------------------------------------
LOGICAL FUNCTION phen_bb_in_leaf( l, n )

USE pftparm, ONLY: phen_bb_min_days
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, n
REAL(KIND=real_jlslsm) :: doy

phen_bb_in_leaf = .FALSE.
IF (.NOT. phen_bb_on(n) .OR. .NOT. ALLOCATED(phen_bb_day)) RETURN
IF (phen_bb_day(l,n) < 0.0 .OR. phen_bb_min_days(n) <= 0.0) RETURN
doy = REAL(day_of_year(current_time%year, current_time%month,                  &
                       current_time%day, l_360, l_leap))
phen_bb_in_leaf = doy < phen_bb_day(l,n) + phen_bb_min_days(n)

END FUNCTION phen_bb_in_leaf

END MODULE phen_budburst_mod
