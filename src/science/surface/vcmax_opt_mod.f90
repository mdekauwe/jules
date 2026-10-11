! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE vcmax_opt_mod

! Optimal photosynthetic capacity (vcmax_opt_model = 1, Farquhar PFTs,
! can_rad_mod = 7), after Sabot et al. (2022, PCE 45: 2744) and Medlyn
! (1996): the N in photosynthetic capacity, N = V25 / n_alloc_vcmax +
! J25 / n_alloc_jmax (Mercado et al. 2018, Eq. 5), is fixed at its PFT
! value, and J25:V25 is chosen every n_day_vcmax_update days to maximise the
! antecedent canopy gross photosynthesis. No N cycle is needed.
!
! Each timestep, the gross photosynthesis of the sunlit and shaded leaves
! is evaluated at n_vopt_jv candidate values of J25:V25 (jv_lo to jv_hi,
! Sabot's JVmin and JVmax), each with its own Vcmax and Jmax on the N
! budget and the leaves' own absorbed light, temperature and ci of the
! timestep; their sum is added to an exponential moving average per
! candidate (time constant n_day_vcmax_opt days). At an update, J25:V25 is
! the maximum of these averages (a parabola through the best candidate
! and its neighbours). So the optimum is that of the actual half-hours,
! rather than of their mean conditions. The rate of each leaf is JULES's,
! MIN(Ac, Aj), at the realised ci of the timestep (not re-optimised for
! each candidate), so that the objective at the current J25:V25 is the
! realised gross photosynthesis. Sabot's An <= 0.5 Vcmax (a TPU proxy) is
! not used: JULES has no TPU limit, and in cool weather, where Vcmax(T) is
! small and Km low, the realised rate often exceeds 0.5 Vcmax, so the cap
! moved N to Rubisco that the model then could not use. The objective is
! gross photosynthesis, so that a lower Vcmax is not favoured for its
! lower Rd.
!
! The capacity is applied as J25:V25 and, with photo_jv_model = jv_scale,
! a factor on Vcmax at the reference temperature (vopt_fac, set by sf_stom
! for the current PFT and read by calc_photo_parameters); with
! photo_jv_model = jv_ntotal, calc_photo_parameters moves V25 with J25:V25
! on the same budget itself. Leaf N, and the respiration scaled from it,
! are unchanged; Rd follows Vcmax through fd. Until the first update the
! PFT parameters are used. The state (vopt_jv, vopt_asum, vopt_tsince) is
! written to dumps.

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

PUBLIC :: n_vopt_jv, vopt_fac, vopt_jv, vopt_asum, vopt_tsince, vopt_v25d,     &
          vopt_jvd, vopt_vleaf, vopt_alloc, vopt_budget_fac, vopt_accumulate,  &
          vopt_update, vopt_vref

! Number, bounds of the candidate J25:V25, mol e- mol-1 CO2 (bounds as
! JVmin, JVmax in Sabot et al. 2022, TractLSM default_params).
INTEGER, PARAMETER :: n_vopt_jv = 21
REAL(KIND=real_jlslsm), PARAMETER :: jv_lo = 1.0, jv_hi = 3.0

! Factor on Vcmax at the reference temperature for the PFT being
! processed in sf_stom (1 = PFT parameters).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: vopt_fac(:)

! State, (land_pts, npft): the optimised J25:V25 (0 until the first
! update), the time since the last update (s), and, (land_pts, npft,
! n_vopt_jv), the moving average of canopy gross photosynthesis (mol CO2
! m-2 s-1) at each candidate J25:V25.
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: vopt_jv(:,:), vopt_tsince(:,:),   &
                                             vopt_asum(:,:,:)

! Diagnostics, (land_pts, npft), for all Farquhar PFTs whatever
! vcmax_opt_model: the top-leaf Vcmax25 and J25:V25 in use, and the
! top-leaf Vcmax at the sunlit leaf temperature (mol CO2 m-2 s-1, with any
! soil-water factor; can_rad_mod = 7).
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE :: vopt_v25d(:,:), vopt_jvd(:,:),    &
                                             vopt_vleaf(:,:)

CONTAINS

! ---------------------------------------------------------------------
! Allocate the state and diagnostics, if not done.
! ---------------------------------------------------------------------
SUBROUTINE vopt_alloc( land_pts )

USE jules_surface_types_mod, ONLY: npft

INTEGER, INTENT(IN) :: land_pts

IF ( ALLOCATED(vopt_fac) ) RETURN
ALLOCATE( vopt_fac(land_pts) )
ALLOCATE( vopt_jv(land_pts,npft), vopt_tsince(land_pts,npft) )
ALLOCATE( vopt_asum(land_pts,npft,n_vopt_jv) )
ALLOCATE( vopt_v25d(land_pts,npft), vopt_jvd(land_pts,npft),                   &
          vopt_vleaf(land_pts,npft) )
vopt_fac(:)       = 1.0
vopt_jv(:,:)      = 0.0
vopt_tsince(:,:)  = 0.0
vopt_asum(:,:,:)  = 0.0
vopt_v25d(:,:)    = 0.0
vopt_jvd(:,:)     = 0.0
vopt_vleaf(:,:)   = 0.0

END SUBROUTINE vopt_alloc

! ---------------------------------------------------------------------
! Candidate k of J25:V25.
! ---------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION jv_cand( k )

INTEGER, INTENT(IN) :: k

jv_cand = jv_lo + REAL(k - 1) * ( jv_hi - jv_lo ) / REAL(n_vopt_jv - 1)

END FUNCTION jv_cand

! ---------------------------------------------------------------------
! Ratio of V25 at J25:V25 = jv_to to V25 at J25:V25 = jv_from on the N
! budget N = V25 / n_alloc_vcmax + J25 / n_alloc_jmax.
! ---------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION vopt_budget_fac( jv_from, jv_to )

USE jules_vegetation_mod, ONLY: n_alloc_jmax, n_alloc_vcmax

REAL(KIND=real_jlslsm), INTENT(IN) :: jv_from, jv_to

vopt_budget_fac = ( 1.0 / n_alloc_vcmax + jv_from / n_alloc_jmax )             &
                  / ( 1.0 / n_alloc_vcmax + jv_to / n_alloc_jmax )

END FUNCTION vopt_budget_fac

! ---------------------------------------------------------------------
! Add this timestep's canopy gross photosynthesis at each candidate J25:V25
! to the moving averages of PFT ft at land point l. For each of the ncls
! leaf classes: vc, jm, the class Vcmax and Jmax at its temperature (mol
! CO2 m-2 s-1, at the current J25:V25, jv_now); i2, the class radiation to
! PSII (mol e- m-2 s-1); ci, ccp, km (Pa). theta: light curvature.
! ---------------------------------------------------------------------
SUBROUTINE vopt_accumulate( l, ft, dt, ncls, vc, jm, i2, ci, ccp, km, jv_now,  &
                            theta )

USE conversions_mod, ONLY: rsec_per_day
USE jules_vegetation_mod, ONLY: n_day_vcmax_opt

INTEGER, INTENT(IN) :: l, ft, ncls
REAL(KIND=real_jlslsm), INTENT(IN) :: dt, vc(ncls), jm(ncls), i2(ncls),        &
                                      ci(ncls), ccp(ncls), km(ncls), jv_now,   &
                                      theta

INTEGER :: c, k
REAL(KIND=real_jlslsm) :: a, ag, fv, fj, vx, jx, je, ac, aj, jvk

a = 1.0 - EXP( -dt / ( n_day_vcmax_opt * rsec_per_day ) )

DO k = 1,n_vopt_jv
  jvk = jv_cand(k)
  fv  = vopt_budget_fac( jv_now, jvk )
  fj  = fv * jvk / jv_now
  ag  = 0.0
  DO c = 1,ncls
    IF ( i2(c) <= 0.0 .OR. vc(c) <= 0.0 .OR. ci(c) <= ccp(c) ) CYCLE
    vx = fv * vc(c)
    jx = fj * jm(c)
    je = ( i2(c) + jx - SQRT( ( i2(c) + jx )**2 - 4.0 * theta * i2(c) * jx ) ) &
         / ( 2.0 * theta )
    ac = vx * ( ci(c) - ccp(c) ) / ( ci(c) + km(c) )
    aj = 0.25 * je * ( ci(c) - ccp(c) ) / ( ci(c) + 2.0 * ccp(c) )
    ag = ag + MIN( ac, aj )
  END DO
  vopt_asum(l,ft,k) = vopt_asum(l,ft,k) + a * ( ag - vopt_asum(l,ft,k) )
END DO

END SUBROUTINE vopt_accumulate

! ---------------------------------------------------------------------
! Advance the clock of PFT ft at land point l by dt and, where an update is
! due, set vopt_jv to the J25:V25 of the largest moving average (kept as
! it is where there has been no photosynthesis).
! ---------------------------------------------------------------------
SUBROUTINE vopt_update( l, ft, dt )

USE conversions_mod, ONLY: rsec_per_day
USE jules_vegetation_mod, ONLY: n_day_vcmax_update

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: dt

INTEGER :: k, kb
REAL(KIND=real_jlslsm) :: y0, y1, y2, den, step

vopt_tsince(l,ft) = vopt_tsince(l,ft) + dt
IF ( vopt_tsince(l,ft) < n_day_vcmax_update * rsec_per_day ) RETURN
vopt_tsince(l,ft) = 0.0

kb = 1
DO k = 2,n_vopt_jv
  IF ( vopt_asum(l,ft,k) > vopt_asum(l,ft,kb) ) kb = k
END DO
IF ( vopt_asum(l,ft,kb) <= 0.0 ) RETURN

vopt_jv(l,ft) = jv_cand(kb)
IF ( kb > 1 .AND. kb < n_vopt_jv ) THEN
  ! Vertex of the parabola through the best candidate and its neighbours.
  y0   = vopt_asum(l,ft,kb-1)
  y1   = vopt_asum(l,ft,kb)
  y2   = vopt_asum(l,ft,kb+1)
  den  = y0 - 2.0 * y1 + y2
  step = ( jv_hi - jv_lo ) / REAL(n_vopt_jv - 1)
  IF ( den < 0.0 ) vopt_jv(l,ft) = vopt_jv(l,ft)                               &
                                   + 0.5 * step * ( y0 - y2 ) / den
END IF

END SUBROUTINE vopt_update

! ---------------------------------------------------------------------
! Vcmax at the reference temperature (mol CO2 m-2 s-1) of PFT ft for leaf
! N nleaf, as calc_photo_parameters (vcmax_ref, before any factor).
! ---------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION vopt_vref( ft, nleaf )

USE jules_vegetation_mod, ONLY: l_trait_phys
USE pftparm, ONLY: neff, vint, vsl

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: nleaf

IF ( l_trait_phys ) THEN
  vopt_vref = ( vsl(ft) * nleaf + vint(ft) ) * 1.0e-6
ELSE
  vopt_vref = neff(ft) * nleaf
END IF

END FUNCTION vopt_vref

END MODULE vcmax_opt_mod
