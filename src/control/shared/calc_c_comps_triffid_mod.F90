! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE calc_c_comps_triffid_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PUBLIC calc_c_comps_triffid, lai_bal_from_ht, ht_from_lai_bal, ht_from_wood,  &
       live_stemwood, sapwood_c

! The TRIFFID allometry in one place: balanced LAI and height, live
! (respiring) stemwood, and the conducting sapwood of the store. Each
! function keeps the arithmetic order of the expressions it stands for, so
! the answers do not change to the bit.
!
! MDK, 20/10/26

!-----------------------------------------------------------------------------
! Description:
!   Calculates carbon contents from vegetation height
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in VEGETATION
!
! Code Description:
!   Language: Fortran 90.
!   https://code.metoffice.gov.uk/trac/jules/wiki/AddingNewSubroutines
!   17 Jan. 2015
!-----------------------------------------------------------------------------
CONTAINS

SUBROUTINE calc_c_comps_triffid(n, ht, lai_bal_pft, leaf, root, wood, c_veg)

USE jules_vegetation_mod, ONLY: l_trait_phys

USE pftparm, ONLY: a_ws, eta_sl, a_wl, b_wl, lma, sigl

USE jules_surface_mod, ONLY: cmass

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Arguments with INTENT(IN).
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                         &
n                   !  PFT number.

REAL(KIND=real_jlslsm), INTENT(IN)    ::                                       &
ht                  ! IN Vegetation height (m).

!-----------------------------------------------------------------------------
! Arguments with INTENT(OUT).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm), INTENT(OUT)   ::                                       &
lai_bal_pft                                                                    &
                      ! OUT Balanced leaf area index
,leaf                                                                          &
                      ! OUT Leaf biomass for balanced LAI (kg C/m2).
,root                                                                          &
                      ! OUT Root biomass (kg C/m2).
,wood                                                                          &
                      ! OUT Woody biomass (kg C/m2).
,c_veg
                      ! OUT Total carbon content of vegetation (kg C/m2).

!-----------------------------------------------------------------------------
!end of header

lai_bal_pft = lai_bal_from_ht(n, ht)   ! MDK, 20/10/26
IF (l_trait_phys) THEN
  leaf = cmass * lma(n) * lai_bal_pft
ELSE
  leaf = sigl(n) * lai_bal_pft
END IF
root = leaf
wood = a_wl(n) * (lai_bal_pft**b_wl(n))
c_veg = leaf + root + wood

END SUBROUTINE calc_c_comps_triffid

!-----------------------------------------------------------------------------
! Balanced LAI from canopy height (m):
!
!   lai_bal = (a_ws eta_sl h / a_wl)^(1 / (b_wl - 1))
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION lai_bal_from_ht(n, ht) RESULT(lai_bal)

USE pftparm, ONLY: a_ws, eta_sl, a_wl, b_wl

IMPLICIT NONE

INTEGER, INTENT(IN) :: n                       ! PFT number.
REAL(KIND=real_jlslsm), INTENT(IN) :: ht       ! Canopy height (m).
REAL(KIND=real_jlslsm) :: lai_bal              ! Balanced LAI (m2 m-2).

lai_bal = (a_ws(n) * eta_sl(n) * ht / a_wl(n))**(1.0 / (b_wl(n) - 1.0))

END FUNCTION lai_bal_from_ht

!-----------------------------------------------------------------------------
! Canopy height (m) from the balanced LAI, the inverse of lai_bal_from_ht:
!
!   h = a_wl / (a_ws eta_sl) lai_bal^(b_wl - 1)
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION ht_from_lai_bal(n, lai_bal) RESULT(ht)

USE pftparm, ONLY: a_ws, eta_sl, a_wl, b_wl

IMPLICIT NONE

INTEGER, INTENT(IN) :: n                       ! PFT number.
REAL(KIND=real_jlslsm), INTENT(IN) :: lai_bal  ! Balanced LAI (m2 m-2).
REAL(KIND=real_jlslsm) :: ht                   ! Canopy height (m).

ht = a_wl(n) / (a_ws(n) * eta_sl(n)) * (lai_bal**(b_wl(n) - 1.0))

END FUNCTION ht_from_lai_bal

!-----------------------------------------------------------------------------
! Canopy height (m) from the wood carbon (kg C m-2), with
! wood = a_wl lai_bal^b_wl:
!
!   h = wood / (a_ws eta_sl) (a_wl / wood)^(1 / b_wl)
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION ht_from_wood(n, wood) RESULT(ht)

USE pftparm, ONLY: a_ws, eta_sl, a_wl, b_wl

IMPLICIT NONE

INTEGER, INTENT(IN) :: n                       ! PFT number.
REAL(KIND=real_jlslsm), INTENT(IN) :: wood     ! Wood carbon (kg C m-2).
REAL(KIND=real_jlslsm) :: ht                   ! Canopy height (m).

ht = wood / (a_ws(n) * eta_sl(n)) * (a_wl(n) / wood)**(1.0 / b_wl(n))

END FUNCTION ht_from_wood

!-----------------------------------------------------------------------------
! Live (respiring) stemwood carbon (kg C m-2), eta_sl h LAI, the standard
! JULES stem (1 / a_ws of the wood at the balanced LAI). Used for stem
! nitrogen and respiration and for the conducting wood that xylem
! recovery renews.
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION live_stemwood(n, ht, lai) RESULT(stemc)

USE pftparm, ONLY: eta_sl

IMPLICIT NONE

INTEGER, INTENT(IN) :: n                       ! PFT number.
REAL(KIND=real_jlslsm), INTENT(IN) :: ht       ! Canopy height (m).
REAL(KIND=real_jlslsm), INTENT(IN) :: lai      ! LAI (m2 m-2).
REAL(KIND=real_jlslsm) :: stemc                ! Live stemwood (kg C m-2).

stemc = eta_sl(n) * ht * lai

END FUNCTION live_stemwood

!-----------------------------------------------------------------------------
! Conducting sapwood carbon (kg C m-2), the xylem X of the TRIFFID store:
! with nsc_la_sa > 0 the pipe model on the balanced leaf area (sapwood area
! LAI_bal / nsc_la_sa times height, wood density and carbon fraction), as in
! FATES; otherwise the fraction nsc_sapwood_frac of the wood carbon
! (Trugman et al. 2018, Table S2).
!
! MDK, 20/10/26
!-----------------------------------------------------------------------------
ELEMENTAL FUNCTION sapwood_c(n, ht, lai_bal, wood) RESULT(x_s)

USE pftparm, ONLY: nsc_la_sa, nsc_sapwood_frac, nsc_wood_density

USE jules_surface_mod, ONLY: cmass

IMPLICIT NONE

INTEGER, INTENT(IN) :: n                       ! PFT number.
REAL(KIND=real_jlslsm), INTENT(IN) :: ht       ! Canopy height (m).
REAL(KIND=real_jlslsm), INTENT(IN) :: lai_bal  ! Balanced LAI (m2 m-2).
REAL(KIND=real_jlslsm), INTENT(IN) :: wood     ! Wood carbon (kg C m-2).
REAL(KIND=real_jlslsm) :: x_s                  ! Sapwood carbon (kg C m-2).

IF (nsc_la_sa(n) > 0.0) THEN
  x_s = lai_bal / nsc_la_sa(n) * ht * nsc_wood_density(n) * cmass
ELSE
  x_s = nsc_sapwood_frac(n) * wood
END IF

END FUNCTION sapwood_c

END MODULE calc_c_comps_triffid_mod
