! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE leaf_resp_mod

! *********************************************************************
! Leaf dark respiration, Rd, before any light inhibition, for the three
! options of leaf_resp_model:
!
!   0 (leaf_resp_jules): Rd = fd * Vcmax(T), so Rd follows the Vcmax
!     temperature curve and falls above its optimum.
!   1 (leaf_resp_heskel): Rd = fd * Vcmax25 * f(T), with the temperature
!     response of Heskel et al. (2016),
!
!       f(T) = exp( b * (T - 25) + c * (T^2 - 25^2) ),   T in degC
!
!   2 (leaf_resp_globresp, C3 leaves): Rd at 25 degC from the GlobResp data
!     of Atkin et al. (2015), with acclimation to the growth temperature,
!     times f(T) (Huntingford et al. 2017, Eqn 1). For the top leaf,
!
!       Rd25 = r0 + r1 * N_a - r2 * T_G   (umol CO2 m-2 s-1)
!
!     with N_a its nitrogen per leaf area (g N m-2) and T_G the growth
!     temperature (degC). A lower leaf has Rd25 in proportion to its
!     nitrogen, as Vcmax. C4 leaves keep option 0.
!
! Options 1 and 2 change the leaves only. JULES scales root and wood
! maintenance respiration with the canopy leaf Rd, so sf_stom scales them
! back to option 0 by rd_default / rd_dark of the top leaf at the canopy
! temperature (leaf_rd gives both).
!
! References:
! Atkin et al., 2015, New Phytol., 206: 614-636.
! Heskel et al., 2016, PNAS, 113: 3832-3837.
! Huntingford et al., 2017, Nat. Commun., 8: 1602.
!
! MDK, 16/11/26
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: rd_heskel_fac, rd25_globresp, leaf_rd

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LEAF_RESP_MOD'

! Grams per kilogram.
REAL(KIND=real_jlslsm), PARAMETER :: g_per_kg = 1000.0

! Moles per micromole.
REAL(KIND=real_jlslsm), PARAMETER :: mol_per_umol = 1.0e-6

CONTAINS

!#############################################################################
!#############################################################################

ELEMENTAL FUNCTION rd_heskel_fac( t_leaf ) RESULT( f_rd )

! The temperature factor of leaf dark respiration of Heskel et al. (2016),
! 1 at 25 degC,
!
!   f(T) = exp( b * (T - 25) + c * (T^2 - 25^2) ),   T in degC
!
! with b = rd_heskel_b and c = rd_heskel_c.
!
! MDK, 16/11/26

USE conversions_mod, ONLY: zerodegc
USE jules_vegetation_mod, ONLY: rd_heskel_b, rd_heskel_c

IMPLICIT NONE

! Leaf temperature (K).
REAL(KIND=real_jlslsm), INTENT(IN) :: t_leaf

REAL(KIND=real_jlslsm) :: f_rd

! Leaf temperature (degC).
REAL(KIND=real_jlslsm) :: t_c

t_c  = t_leaf - zerodegc
f_rd = EXP( rd_heskel_b * (t_c - 25.0) + rd_heskel_c * (t_c**2 - 625.0) )

END FUNCTION rd_heskel_fac

!#############################################################################
!#############################################################################

SUBROUTINE rd25_globresp( ft, land_pts, veg_pts, veg_index, nleaf_top,         &
                          t_growth, rd25_per_nleaf )

! GlobResp leaf dark respiration at 25 degC (leaf_resp_model 2, C3 leaves),
! per unit of leaf nitrogen in the units of nleaf, so that a leaf with
! nitrogen nleaf has
!
!   Rd25 = rd25_per_nleaf * nleaf   (mol CO2 m-2 s-1)
!
! For the top leaf (Huntingford et al. 2017, Eqn 1),
!
!   Rd25_top = r0 + r1 * N_a - r2 * T_G   (umol CO2 m-2 s-1)
!
! no lower than 0, with r0 = rd_r0 (per PFT), r1 = rd_globresp_r1,
! r2 = rd_globresp_r2, T_G the growth temperature (degC) and N_a the top
! leaf's nitrogen per leaf area (g N m-2). With l_trait_phys, nleaf_top is
! already in g N m-2 (nmass * LMA * 1000). Without it, nleaf_top is the leaf
! N:C ratio (kg N kg C-1), so
!
!   N_a = nleaf_top * sigl * 1000
!
! with sigl the leaf carbon per leaf area (kg C m-2). rd25_per_nleaf is 0
! for other options and for C4 PFTs.
!
! MDK, 16/11/26

USE conversions_mod, ONLY: zerodegc
USE jules_vegetation_mod, ONLY: leaf_resp_model, leaf_resp_globresp,           &
                                rd_globresp_r1, rd_globresp_r2, l_trait_phys
USE pftparm, ONLY: c3, rd_r0, sigl

IMPLICIT NONE

! PFT index.
INTEGER, INTENT(IN) :: ft

! Number of land points.
INTEGER, INTENT(IN) :: land_pts

! Number of vegetated points.
INTEGER, INTENT(IN) :: veg_pts

! Land index of each vegetated point.
INTEGER, INTENT(IN) :: veg_index(land_pts)

! Nitrogen of the top leaf (g N m-2 with l_trait_phys, else kg N kg C-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: nleaf_top(land_pts)

! Growth temperature (K).
REAL(KIND=real_jlslsm), INTENT(IN) :: t_growth(land_pts)

! Rd at 25 degC per unit of nleaf (mol CO2 m-2 s-1 per unit of nleaf).
REAL(KIND=real_jlslsm), INTENT(OUT) :: rd25_per_nleaf(land_pts)

INTEGER :: l, m

! Top leaf's nitrogen per leaf area (g N m-2).
REAL(KIND=real_jlslsm) :: narea_top

! Top leaf's Rd at 25 degC (umol CO2 m-2 s-1).
REAL(KIND=real_jlslsm) :: rd25_top

rd25_per_nleaf(:) = 0.0

IF ( leaf_resp_model /= leaf_resp_globresp .OR. c3(ft) /= 1 ) RETURN

DO m = 1,veg_pts
  l = veg_index(m)

  IF ( l_trait_phys ) THEN
    narea_top = nleaf_top(l)
  ELSE
    narea_top = nleaf_top(l) * sigl(ft) * g_per_kg
  END IF

  rd25_top = MAX( rd_r0(ft) + rd_globresp_r1 * narea_top                       &
                  - rd_globresp_r2 * (t_growth(l) - zerodegc), 0.0 )

  IF ( nleaf_top(l) > 0.0 ) THEN
    rd25_per_nleaf(l) = mol_per_umol * rd25_top / nleaf_top(l)
  END IF

END DO

END SUBROUTINE rd25_globresp

!#############################################################################
!#############################################################################

SUBROUTINE leaf_rd( ft, land_pts, veg_pts, veg_index, rd_default, f_vcmax_t,   &
                    rd_temp_fac, nleaf, rd25_per_nleaf, rd_dark )

! Leaf dark respiration before light inhibition for the chosen
! leaf_resp_model, from option 0's value rd_default = fd * Vcmax(T):
!
!   0, and 2 for C4: rd_dark = rd_default
!   1:               rd_dark = rd_default / f_V(T) * f(T)
!                            = fd * Vcmax25 * f(T)
!   2 (C3):          rd_dark = rd25_per_nleaf * nleaf * f(T)
!
! with f_V the Vcmax temperature factor (f_vcmax_t), so Vcmax25 includes any
! acclimation and nitrogen allocation, and f the factor of Heskel et al.
! (2016) at the leaf temperature (rd_temp_fac).
!
! MDK, 16/11/26

USE jules_vegetation_mod, ONLY: leaf_resp_model, leaf_resp_heskel,            &
                                leaf_resp_globresp
USE pftparm, ONLY: c3

IMPLICIT NONE

! PFT index.
INTEGER, INTENT(IN) :: ft

! Number of land points.
INTEGER, INTENT(IN) :: land_pts

! Number of vegetated points.
INTEGER, INTENT(IN) :: veg_pts

! Land index of each vegetated point.
INTEGER, INTENT(IN) :: veg_index(land_pts)

! Rd of option 0, fd * Vcmax at the leaf temperature (mol CO2 m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: rd_default(land_pts)

! Vcmax temperature factor at the leaf temperature (Vcmax / Vcmax25).
REAL(KIND=real_jlslsm), INTENT(IN) :: f_vcmax_t(land_pts)

! Heskel temperature factor of Rd at the leaf temperature.
REAL(KIND=real_jlslsm), INTENT(IN) :: rd_temp_fac(land_pts)

! Leaf nitrogen (g N m-2 with l_trait_phys, else kg N kg C-1).
REAL(KIND=real_jlslsm), INTENT(IN) :: nleaf(land_pts)

! GlobResp Rd at 25 degC per unit of nleaf (rd25_globresp).
REAL(KIND=real_jlslsm), INTENT(IN) :: rd25_per_nleaf(land_pts)

! Leaf dark respiration before light inhibition (mol CO2 m-2 s-1).
REAL(KIND=real_jlslsm), INTENT(OUT) :: rd_dark(land_pts)

INTEGER :: l, m

rd_dark(:) = 0.0

DO m = 1,veg_pts
  l = veg_index(m)

  IF ( leaf_resp_model == leaf_resp_heskel ) THEN
    rd_dark(l) = rd_default(l) / MAX( f_vcmax_t(l), TINY(1.0_real_jlslsm) )   &
                 * rd_temp_fac(l)
  ELSE IF ( leaf_resp_model == leaf_resp_globresp .AND. c3(ft) == 1 ) THEN
    rd_dark(l) = rd25_per_nleaf(l) * nleaf(l) * rd_temp_fac(l)
  ELSE
    rd_dark(l) = rd_default(l)
  END IF

END DO

END SUBROUTINE leaf_rd

END MODULE leaf_resp_mod
