! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE phen_budburst_mod

!-----------------------------------------------------------------------------
! Budburst for deciduous TRIFFID phenology (phen_bb_doy or phen_bb_gdd set).
!
! Standard phenol lets leaves regrow whenever the day's leaf turnover is
! below twice its minimum, so a deciduous canopy regrows in mild winter
! spells and its spring flush is set by the same temperature threshold
! as leaf fall (tleaf_of, dgl_dt). Here a PFT is dormant from the autumn
! fall of its canopy to the floor (phen 0.01) until budburst, and phenol
! does not regrow leaves while it is dormant or before this year's
! budburst. Budburst is on day of year phen_bb_doy, or when the growing
! degree-days from 1 January, the sum of the daily mean air temperature
! above phen_bb_tbase, reach phen_bb_gdd (e.g. oak about 730 degC d above
! 0 degC at Alice Holt, from MODIS green-up). Leaf fall is the standard
! cold leaf drop, held off for phen_bb_min_days after budburst: the
! half-hourly cold hinge it averages counts cold nights, which otherwise
! stop the spring flush and can drop a summer canopy (and so start
! dormancy) in a cool spell.
! The year starts on 1 January (northern-hemisphere seasons). A point that
! never reaches phen_bb_gdd in a year stays leafless that year, as a site
! too cold for the PFT; set the gate only for deciduous PFTs.
!
! State (pft, land; in dumps with the gate on):
!   phen_gdd      growing degree-days since 1 January (degC d)
!   phen_bb_day   day of year of this year's budburst; -1 = not yet
!   phen_dormant  1 = dormant (canopy fallen, waiting for budburst)
! phen_bb_day is public so other schemes (e.g. the store's leaf flush) can
! start at the model's budburst.
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
! Is the budburst gate on for PFT n?
!-----------------------------------------------------------------------------
LOGICAL FUNCTION phen_bb_on( n )

USE pftparm, ONLY: phen_bb_doy, phen_bb_gdd
USE pftparm, ONLY: is_set   ! MDK, 24/10/26

IMPLICIT NONE

INTEGER, INTENT(IN) :: n

phen_bb_on = is_set(phen_bb_doy(n)) .OR. is_set(phen_bb_gdd(n))

END FUNCTION phen_bb_on

!-----------------------------------------------------------------------------
! Allocate the state at its fresh-start values.
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
! Every time step: degree-days and budburst of PFT ft.
!-----------------------------------------------------------------------------
SUBROUTINE phen_bb_update( land_pts, veg_pts, veg_index, ft, timestep, tair )

USE pftparm, ONLY: phen_bb_doy, phen_bb_gdd, phen_bb_tbase
USE pftparm, ONLY: is_set   ! MDK, 24/10/26
USE model_time_mod, ONLY: current_time
USE datetime_utils_mod, ONLY: day_of_year
USE time_info_mod, ONLY: l_360, l_leap
USE water_constants_mod, ONLY: tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) ::                                                         &
  land_pts,                                                                    &
    ! Number of land points.
  veg_pts,                                                                     &
    ! Number of points with this PFT.
  veg_index(land_pts),                                                         &
    ! Index of those points on the land grid.
  ft
    ! Plant functional type.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
    ! Model time step (s).
  tair(land_pts)
    ! Air temperature (K).

INTEGER :: j, l
REAL(KIND=real_jlslsm) :: doy
    ! Day of the year.
    !
    ! MDK, 24/10/26
LOGICAL :: l_new_year
    ! First time step of the year.
    !
    ! MDK, 24/10/26

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
! In phenol: no leaf regrowth for PFT n at land point l while dormant or
! before this year's budburst.
!-----------------------------------------------------------------------------
LOGICAL FUNCTION phen_bb_blocked( l, n )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, n

phen_bb_blocked = .FALSE.
IF (.NOT. phen_bb_on(n) .OR. .NOT. ALLOCATED(phen_bb_day)) RETURN
phen_bb_blocked = phen_bb_day(l,n) < 0.0 .OR. phen_dormant(l,n) > 0.5

END FUNCTION phen_bb_blocked

!-----------------------------------------------------------------------------
! In phenol: PFT n at land point l is within phen_bb_min_days of this
! year's budburst (no cold leaf drop).
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
