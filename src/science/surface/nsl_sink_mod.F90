! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE nsl_sink_mod

!-----------------------------------------------------------------------------
! Sink-limited nonstomatal limitation of photosynthesis (l_som_nsl_sink).
!
! A source-sink imbalance pool S (kg C m-2) fills with gross photosynthesis
! A and empties with the sink demand U:
!   dS/dt = A - U,
!   U = U_max S/(S + K) [ m Q10^((T - 25)/10) + (1 - m) g(psi) h(T) ],
! where m is the maintenance fraction of U_max (Q10 = 2), g the turgor limit
! on growth from the root-zone (predawn) water potential,
!   g = 1 / (1 + exp(sf (psi_g50 - psi))), psi in MPa,
! and h a cold limit on growth (logistic, half at 5 C, Lempereur et al.
! 2015, New Phytol. 207: 579: Q. ilex stem growth stops below about 5 C and
! below a predawn psi of about -1.1 MPa while GPP stays positive).
! Photosynthetic capacity (Vcmax, Jmax, so Rd) is scaled by
!   f = 1 - S/S0, S0 = U_max tau,
! the sugar-regulated NSL of Dewar et al. (2022, New Phytol. 233: 639,
! Eqn 3a), applied as a state outside the stomatal optimisation (which
! keeps its hydraulic cost). When growth stops (drought or cold) the pool
! fills and capacity falls until A is drawn down to the sink demand; after
! rain the pool has to drain before capacity recovers. K = 0.02 S0.
!
! The pool is held here (not in the dump): it starts empty at the start of
! a run. Profit max (stomata_model = 4) with Farquhar only.
!-----------------------------------------------------------------------------

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC nsl_sink_factor, nsl_sink_update

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: s_pool(:,:)
    ! Imbalance pool (kg C m-2), (land_pts, npft).

REAL(KIND=real_jlslsm), PARAMETER ::                                           &
  q10_maint = 2.0,                                                             &
    ! Q10 of the maintenance sink.
  t_cold_half = 5.0,                                                           &
    ! Temperature (C) at which the cold limit halves growth.
  t_cold_width = 1.0,                                                          &
    ! Width (K) of the cold limit.
  k_frac = 0.02,                                                               &
    ! Half-saturation of the sink, as a fraction of S0.
  gcm2d_to_kgm2s = 1.0e-3 / 86400.0
    ! g C m-2 d-1 to kg C m-2 s-1.

CONTAINS

!#############################################################################

SUBROUTINE nsl_sink_alloc( land_pts )

USE jules_surface_types_mod, ONLY: npft

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts

IF ( .NOT. ALLOCATED(s_pool) ) THEN
  ALLOCATE( s_pool(land_pts, npft) )
  s_pool(:,:) = 0.0
END IF

END SUBROUTINE nsl_sink_alloc

!#############################################################################

SUBROUTINE nsl_sink_factor( ft, land_pts, veg_pts, veg_index, f_vc )

! Scale the capacity factor f_vc by 1 - S/S0 (veg points only).

USE pftparm, ONLY: nsl_sink_umax, nsl_sink_tau

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: f_vc(land_pts)

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: s0

CALL nsl_sink_alloc( land_pts )
s0 = nsl_sink_umax(ft) * gcm2d_to_kgm2s * nsl_sink_tau(ft) * 86400.0

DO m = 1,veg_pts
  l = veg_index(m)
  f_vc(l) = f_vc(l) * MAX( 1.0 - s_pool(l,ft) / s0, 0.0 )
END DO

END SUBROUTINE nsl_sink_factor

!#############################################################################

SUBROUTINE nsl_sink_update( ft, land_pts, veg_pts, veg_index, gpp, tair,       &
                            psi_rz )

! Advance the pool one timestep with this timestep's GPP.

USE pftparm, ONLY: nsl_sink_umax, nsl_sink_tau, nsl_sink_maint,               &
                   nsl_sink_psi50, nsl_sink_sf
USE timestep_mod, ONLY: timestep
USE conversions_mod, ONLY: zerodegc

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  gpp(land_pts),                                                               &
    ! Gross primary productivity (kg C m-2 s-1).
  tair(land_pts),                                                              &
    ! Air temperature (K).
  psi_rz(land_pts)
    ! Root-zone water potential (Pa).

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: s0, umax, tc, g_turgor, h_cold, f_maint, u_sink,     &
                          psi_mpa

CALL nsl_sink_alloc( land_pts )
umax = nsl_sink_umax(ft) * gcm2d_to_kgm2s
s0   = umax * nsl_sink_tau(ft) * 86400.0

DO m = 1,veg_pts
  l = veg_index(m)
  tc       = tair(l) - zerodegc
  psi_mpa  = MIN(psi_rz(l), 0.0) * 1.0e-6
  g_turgor = 1.0 / ( 1.0 + EXP( nsl_sink_sf(ft)                                &
                                * ( nsl_sink_psi50(ft) * 1.0e-6 - psi_mpa ) ) )
  h_cold   = 1.0 / ( 1.0 + EXP( (t_cold_half - tc) / t_cold_width ) )
  f_maint  = q10_maint ** ( 0.1 * (tc - 25.0) )
  u_sink   = umax * s_pool(l,ft) / ( s_pool(l,ft) + k_frac * s0 )              &
             * ( nsl_sink_maint(ft) * f_maint                                  &
                 + (1.0 - nsl_sink_maint(ft)) * g_turgor * h_cold )
  s_pool(l,ft) = s_pool(l,ft) + timestep * ( MAX(gpp(l), 0.0) - u_sink )
  s_pool(l,ft) = MIN( MAX(s_pool(l,ft), 0.0), s0 )
END DO

END SUBROUTINE nsl_sink_update

END MODULE nsl_sink_mod
