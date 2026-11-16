! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    SUBROUTINE ROOT_FRAC---------------------------------------------

MODULE root_frac_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ROOT_FRAC_MOD'

CONTAINS
! Subroutine Interface:
SUBROUTINE root_frac (ft,nlayer,dz,rootd,f_root)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE pftparm, ONLY: fsmc_mod
USE pftparm, ONLY: rootd_max_ft, is_set   ! MDK, 03/11/26
IMPLICIT NONE

! Description:
!     Calculates the fraction of the total plant roots within each
!     soil layer.

! Documentation : UM Documentation Paper 25


! Subroutine arguments:
!   Scalar arguments with intent(IN) :
INTEGER, INTENT(IN) ::                                                         &
 ft,                                                                           &
                      ! IN Plant functional type.
 nlayer
                      ! IN Number of soil layers.

REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
 dz(nlayer),                                                                   &
                      ! IN Soil layer thicknesses (m).
 rootd                ! IN Rootdepth (m).

!   Array arguments with intent(OUT) :
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
 f_root(nlayer)       ! OUT Fraction of roots in each soil
                      !     layer.
! Local scalars:
INTEGER ::                                                                     &
 n                    ! WORK Loop counters

REAL(KIND=real_jlslsm) ::                                                      &
 ftot,                                                                         &
                      ! WORK Normalisation factor.
 ztot,                                                                         &
                      ! WORK Total depth of soil (m).
 z1,z2,                                                                        &
                      ! WORK Depth of the top and bottom of the
                      !      soil layers (m).
 rootd_max,                                                                    &
                      ! WORK Maximum rooting depth of this PFT (m),
                      !      rootd_max_ft(ft); no roots below it.
 rootd_lim            ! WORK Depth the roots reach with fsmc_mod = 1
                      !      (m): the shallower of rootd and
                      !      rootd_max_ft(ft), or rootd if that is
                      !      unset.
                      !
                      !      MDK, 16/11/26

! Local parameters:
REAL(KIND=real_jlslsm), PARAMETER :: p = 1.0
                      ! WORK Power describing depth dependence
                      !           of the root density profile.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ROOT_FRAC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
z2 = 0.0
ztot = 0.0
ftot = 0.0

IF (fsmc_mod(ft) == 1) THEN
  ! Uniform root density down to rootd_lim. With a maximum rooting depth
  ! (rootd_max_ft) the roots stop at the shallower of rootd and
  ! rootd_max_ft; without one rootd_lim is rootd, as before.
  !
  ! MDK, 16/11/26
  rootd_lim = rootd
  IF ( is_set(rootd_max_ft(ft)) ) rootd_lim = MIN(rootd, rootd_max_ft(ft))

  ! Each layer's roots are its thickness above rootd_lim: all of it for a
  ! layer above rootd_lim, none for a layer below it, and the part above
  ! rootd_lim for the layer it falls in.
  DO n = 1,nlayer
    z1 = z2
    z2 = z2 + dz(n)
    ztot = ztot + dz(n)

    IF (z1 > rootd_lim) THEN
      f_root(n) = 0.0
    ELSE IF (z2 < rootd_lim) THEN
      f_root(n) = dz(n)
    ELSE
      f_root(n) = rootd_lim - z1
    END IF

  END DO

  ! Divide by the rooted depth, so the fractions add to 1. A column
  ! shallower than rootd_lim is rooted to its bottom.
  ftot = MIN(rootd_lim, ztot)
  f_root(:) = f_root(:) / ftot
ELSE IF ( is_set(rootd_max_ft(ft)) ) THEN
  ! Exponential root density with e-folding depth rootd, cut off at the
  ! maximum rooting depth rootd_max (rootd_max_ft). A layer from z1 to z2
  ! holds the roots between z1 and MIN(z2, rootd_max),
  !
  !   f_root(n) = exp(-p * z1 / rootd) - exp(-p * MIN(z2, rootd_max) / rootd),
  !
  ! and a layer wholly below rootd_max holds none.
  !
  ! MDK, 16/11/26
  rootd_max = rootd_max_ft(ft)

  DO n = 1,nlayer
    z1 = z2
    z2 = z2 + dz(n)
    ztot = ztot + dz(n)

    IF (z1 >= rootd_max) THEN
      f_root(n) = 0.0
    ELSE
      f_root(n) = EXP(-p * z1 / rootd)                                         &
                  - EXP(-p * MIN(z2, rootd_max) / rootd)
    END IF

  END DO

  ! Divide by the roots in the whole rooted depth,
  !
  !   ftot = 1 - exp(-p * MIN(ztot, rootd_max) / rootd),
  !
  ! so the fractions add to 1 and the roots cut off below rootd_max go to
  ! the layers above it in proportion. A rootd_max at or below the
  ! bottom of the column gives the profile without the cut-off.
  ftot = 1.0 - EXP(-p * MIN(ztot, rootd_max) / rootd)
  DO n = 1,nlayer
    f_root(n) = f_root(n) / ftot
  END DO
ELSE
  DO n = 1,nlayer
    z1 = z2
    z2 = z2 + dz(n)
    ztot = ztot + dz(n)
    f_root(n) = EXP(-p * z1 / rootd) - EXP(-p * z2 / rootd)

  END DO

  ftot = 1.0 - EXP(-p * ztot / rootd)
  DO n = 1,nlayer
    f_root(n) = f_root(n) / ftot
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE root_frac

END MODULE root_frac_mod
