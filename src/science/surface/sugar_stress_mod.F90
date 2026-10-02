! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE sugar_stress_mod

! *********************************************************************
! Feedbacks between the SUGAR NSC pool (sugar_mod) and the stomatal
! optimisation (stom_opt_jls_mod), and water stress on SUGAR growth.
!
! The optimisation maximises gain - cost, gain = A / A_max. Here the value
! of carbon depends on whether the plant can use or store it:
!   som_nsc_feedback = 1 (weight): gain = w A / A_max, with
!     w = 1 - (1 - nsc_w_min) MIN(1, f_nsc / nsc_f_full)**nsc_w_k,
!     so carbon is worth less as the pool fills. This only moves the
!     optimum where the hydraulic cost is not ~0.
!   som_nsc_feedback = 2 (cap): gross A beyond the sink demand S
!     (sugar_sink_demand: SUGAR growth and respiration plus filling the
!     pool to nsc_f_full over nsc_tau_fill) has no value:
!     gain = (A - MAX(0, A + Rd - S)) / A_max. The gain is flat beyond S,
!     so any hydraulic cost puts the optimum where A first meets S.
!   som_nsc_feedback = 3: both.
! som_nsc_cap_curv < 1 softens the cap to a co-limitation (as Collatz),
! theta Ae**2 - (A + S) Ae + A S = 0, so carbon above S keeps some value
! and the stomata still follow the weather when the cap binds. The sink
! demand includes the leaf flush (l_sugar_leaf_flush, leaf_flush_update).
! A_max stays the uncapped, unweighted maximum, so carbon the plant cannot
! use counts as worth less, not renormalised to a full gain. With w = 1 and
! S = HUGE the gain is exactly the original. nsc_f_full = 0 turns the
! feedback off for a PFT (e.g. grasses with no storage).
!
! l_sugar_turgor: SUGAR structural growth is scaled by
!   f_turgor = 1 / (1 + exp(sf_growth (psi_g50 - psi))),
! as growth (cell expansion) stops at milder water potentials than
! photosynthesis (Hsiao 1973; Muller et al. 2011). psi is the root-zone
! water potential, the predawn proxy: stem growth happens mostly at night
! when turgor recovers, and cessation thresholds are given for predawn psi.
! Growth respiration goes with growth; maintenance does not. The value is
! kept for the next step's NSC cap (f_turgor_prev; 1 on the first step and
! after a restart).
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: nsc_gain_weight, turgor_growth_factor, nsc_gain_value,              &
          f_turgor_prev, sugar_stress_init, leaf_flush_update, flush_rate

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: f_turgor_prev(:,:)
    ! Turgor growth factor of the previous time step (land_pts, npft).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: lai_prev(:,:)
    ! LAI of the previous time step (land_pts, npft); < 0 until first set.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: flush_rate(:,:)
    ! Leaf-flush construction, new structural leaf carbon smoothed over
    ! tau_flush (kg C m-2 s-1; land_pts, npft). l_sugar_leaf_flush.
REAL(KIND=real_jlslsm), PARAMETER :: tau_flush = 5.0 * 86400.0
    ! Smoothing time of the leaf flush (s): prescribed LAI changes in daily
    ! steps; the exponential average keeps the total carbon.

CONTAINS

!-----------------------------------------------------------------------------
! Allocate f_turgor_prev on first use.
!-----------------------------------------------------------------------------
SUBROUTINE sugar_stress_init( land_pts, npft )
INTEGER, INTENT(IN) :: land_pts, npft
IF ( .NOT. ALLOCATED(f_turgor_prev) ) THEN
  ALLOCATE( f_turgor_prev(land_pts, npft) )
  ALLOCATE( lai_prev(land_pts, npft) )
  ALLOCATE( flush_rate(land_pts, npft) )
  f_turgor_prev(:,:) = 1.0
  lai_prev(:,:)      = -1.0
  flush_rate(:,:)    = 0.0
END IF
END SUBROUTINE sugar_stress_init

!-----------------------------------------------------------------------------
! Update the leaf flush of PFT ft at land point l from this step's LAI
! (l_sugar_leaf_flush): new leaf carbon lma cmass max(0, dLAI), as a rate,
! exponentially averaged over tau_flush. Call once per time step.
!-----------------------------------------------------------------------------
SUBROUTINE leaf_flush_update( l, ft, lai, leafc_per_lai, timestep )
INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: lai, leafc_per_lai, timestep
    ! LAI, leaf carbon per unit LAI (kg C m-2), time step (s).
REAL(KIND=real_jlslsm) :: rate, a
IF ( lai_prev(l,ft) < 0.0 ) lai_prev(l,ft) = lai
rate = leafc_per_lai * MAX( lai - lai_prev(l,ft), 0.0 ) / timestep
a    = MIN( 1.0, timestep / tau_flush )
flush_rate(l,ft) = flush_rate(l,ft) + a * ( rate - flush_rate(l,ft) )
lai_prev(l,ft)   = lai
END SUBROUTINE leaf_flush_update

!-----------------------------------------------------------------------------
! Gain weight from the NSC pool (som_nsc_feedback = 1 or 3).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION nsc_gain_weight( ft, f_nsc ) RESULT( w )
USE pftparm, ONLY: nsc_f_full, nsc_w_min, nsc_w_k
INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: f_nsc
IF ( nsc_f_full(ft) > 0.0 ) THEN
  w = 1.0 - ( 1.0 - nsc_w_min(ft) )                                            &
            * MIN( 1.0, MAX( f_nsc, 0.0 ) / nsc_f_full(ft) )**nsc_w_k(ft)
ELSE
  w = 1.0
END IF
END FUNCTION nsc_gain_weight

!-----------------------------------------------------------------------------
! Turgor limit on SUGAR growth (l_sugar_turgor). psi in Pa.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION turgor_growth_factor( ft, psi ) RESULT( f )
USE pftparm, ONLY: psi_g50, sf_growth
INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi
REAL(KIND=real_jlslsm) :: x
! Argument in MPa; bounded so EXP cannot overflow.
x = MIN( sf_growth(ft) * ( psi_g50(ft) - psi ) * 1.0e-6, 50.0 )
f = 1.0 / ( 1.0 + EXP(x) )
END FUNCTION turgor_growth_factor

!-----------------------------------------------------------------------------
! Carbon gain of a leaf state before normalisation: al + g_off (net or gross
! A, see l_som_gain_gross), less the gross A that the sink demand cap does not
! count, times w. With som_nsc_cap_curv = 1 the counted gross A is
! min(A, cap); with theta < 1 it is the smaller root of
! theta Ae**2 - (A + cap) Ae + A cap = 0, a smooth co-limitation. Exactly
! al + g_off when w = 1 and cap = HUGE.
!-----------------------------------------------------------------------------
ELEMENTAL REAL(KIND=real_jlslsm) FUNCTION nsc_gain_value( al, g_off, rd, w,    &
                                                          cap ) RESULT( g )
USE jules_vegetation_mod, ONLY: som_nsc_cap_curv
REAL(KIND=real_jlslsm), INTENT(IN) :: al, g_off, rd, w, cap
REAL(KIND=real_jlslsm) :: ag, ae, th
th = som_nsc_cap_curv
IF ( th >= 1.0 .OR. cap >= HUGE(cap) ) THEN
  g = w * ( al + g_off - MAX( 0.0_real_jlslsm, al + rd - cap ) )
ELSE
  ag = MAX( al + rd, 0.0_real_jlslsm )
  ae = ( ( ag + cap ) - SQRT( MAX( ( ag + cap )**2 - 4.0 * th * ag * cap,      &
                                   0.0_real_jlslsm ) ) ) / ( 2.0 * th )
  g = w * ( al + g_off - ( ag - ae ) )
END IF
END FUNCTION nsc_gain_value

END MODULE sugar_stress_mod
