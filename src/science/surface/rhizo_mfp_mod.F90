! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

! *****************************************************************************
! Nonlinear soil-to-root (rhizosphere) link, l_som_rhizo_mfp.
!
! With l_som_rhizo_series the soil link is linear in E: each layer conducts
! k_i = B_i K_i(psi_soil,i), with the soil conductivity K at the bulk layer
! psi, and psi_root = psi_src - E / sum k_i. Steady radial flow to a root
! through soil whose K falls towards the root surface is instead, per layer
! (Gardner 1960; Schroeder et al. 2008; Carminati & Javaux 2020),
!   E_i = B_i ( Phi_i(h_soil,i) - Phi_i(h_root) ),
!   Phi_i(h) = integral of K_i(h') dh' from h to the driest state,
! the matric flux potential (Kirchhoff transform), with h the suction (m)
! and B_i the geometric factor of soil_to_root_conductance (Bonan et al.
! 2014 eq. A23, the steady-state cylinder 2 pi L / ln(r_s/r_r)), so that
! soil_to_root_k = B_i K_i. All layers feed one root node, so for a given E
! the root suction h_root solves
!   sum_i B_i ( Phi_i(h_soil,i) - Phi_i(h_root) ) = E,
! with layers drier than the root giving nothing unless l_som_rhizo_hr. The
! left side increases with h_root, so the root is unique; the soil can
! supply at most sum_i B_i Phi_i(h_soil,i). The marginal soil conductance
! dE/dpsi_root = sum_i B_i K_i(h_root) / (rho_water g) is what the hydraulic
! cost and kcrit see. At small E or in wet soil this is the linear link.
!
! Phi has no closed form for van Genuchten-Mualem, so Phi / (K_sat sathh)
! is tabulated as a function of x = h / sathh, once for each distinct
! exponent b (it depends on b only), with JULES's own K(Se) (hyd_con_vg or
! hyd_con_ch, including their clamps) and Se(h) (hh_from_sthu). Values
! between nodes use cubic Hermite interpolation in ln x with the exact
! slope dPhi/d ln x = -K x, so Phi and K stay consistent.
!
! Units (per m2 ground): h, sathh (m); K, K_sat (kg m-2 s-1 per unit head
! gradient); B (m-1); E (mol m-2 s-1, as el_pft and ksr); conductances
! mol m-2 s-1 Pa-1 (as ksr). Per layer c_i = B_i K_sat,i sathh_i / m_h2o
! (mol m-2 s-1), so E_i = c_i (G(x_soil,i) - G(x_root,i)) with G the
! dimensionless table.
!
! The state is per PFT (physiol sets it after smc_ext with rhizo_set_state;
! rz_ft is the PFT of the current solve). A leaf path that carries a share
! s of the plant (set_ksr_path) sees E / s against the whole root system
! and has conductance s dE/dpsi, as the linear ksr_path = s ksr.
! *****************************************************************************
MODULE rhizo_mfp_mod

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE

INTEGER, PARAMETER :: rk = SELECTED_REAL_KIND(15, 300)
      ! Working precision of this module: double. The model reals are single
      ! (real_jlslsm), but E_i = c_i (G_s,i - G_i(h_root)) is a small
      ! difference of two numbers 1e3-1e6 x E in wet soil, so the tables,
      ! state and solver are double and only the interface is real_jlslsm.
PUBLIC :: rhizo_set_state, rhizo_drop, rhizo_path_drop, rhizo_k_zero,          &
          rhizo_uptake, rhizo_psi0, rhizo_active, rz_share, rz_ft,             &
          rhizo_geometry, rhizo_active_ft

CHARACTER(LEN=*), PARAMETER :: ModuleName = 'RHIZO_MFP_MOD'

INTEGER :: rz_ft = 0
      ! PFT of the current stomatal solve (set by physiol before sf_stom).
REAL(KIND=rk), ALLOCATABLE :: rz_share(:)
      ! Share of the plant carried by the current leaf path (set_ksr_path).

LOGICAL, ALLOCATABLE :: rz_on(:,:)
      ! (land_pts, npft): state set for this point and PFT.
REAL(KIND=rk), ALLOCATABLE ::                                         &
  rz_c(:,:,:),                                                                 &
      ! (land_pts, nl, npft): B K_sat sathh / m_h2o (mol m-2 s-1).
  rz_sathh(:,:,:),                                                             &
      ! Saturated suction, 1/alpha for van Genuchten (m).
  rz_hs(:,:,:),                                                                &
      ! Bulk soil suction of the layer (m).
  rz_gs(:,:,:),                                                                &
      ! G at the bulk soil suction.
  rz_h0(:,:)
      ! (land_pts, npft): root suction at zero flow (m).
INTEGER, ALLOCATABLE :: rz_itab(:,:,:)
      ! Table of the layer's exponent b.
INTEGER :: rz_nl = 0

! Tables of G(x) = Phi / (K_sat sathh), one per distinct b.
INTEGER, PARAMETER :: n_tab_max = 32, n_g = 4001
REAL(KIND=rk), PARAMETER :: x_lo = 1.0e-6_rk
      ! Bottom of the table in x (near saturation, K ~ K_sat below it).
REAL(KIND=rk), PARAMETER :: x_cap = 1.0e12_rk
      ! Cap on the driest x (overflow guard).
REAL(KIND=rk), PARAMETER :: sthu_min_rz = 0.01_rk
      ! Driest state, as smc_ext's sthu_min for psi_from_sthu.
INTEGER :: n_tab = 0
REAL(KIND=rk) :: tab_b(n_tab_max), tab_u_lo(n_tab_max),              &
                          tab_du(n_tab_max), tab_x_max(n_tab_max)
REAL(KIND=rk) :: tab_g(n_g, n_tab_max), tab_dg(n_g, n_tab_max)
      ! G and dG/du (u = ln x) at the nodes.

REAL(KIND=rk), PARAMETER :: m_h2o_rz = 0.018015_rk
      ! kg mol-1, as physiol's m_h2o_rs.

INTEGER, PARAMETER :: max_iter = 60

! Solver statistics, printed every n_report solves.
INTEGER(KIND=8) :: n_solve = 0, n_iter_sum = 0
INTEGER :: n_iter_max = 0
INTEGER(KIND=8), PARAMETER :: n_report = 2000000

CONTAINS

! *****************************************************************************
! Geometric factor B_i (m-1) of soil_to_root_conductance (smc_ext), i.e.
! soil_to_root_k / soil_k: the Bonan et al. (2014) eq. A23 cylinder per unit
! root length, 4 pi / ln(rho_r dz / (f_root m_root)), times the root length
! per m2 ground, f_root m_root / (rho_r pi r_r^2).
! *****************************************************************************
FUNCTION rhizo_geometry( ft, nl, f_root, root_mass ) RESULT( b_geom )

USE pftparm, ONLY: root_radi_pft, rootc_density_pft
USE conversions_mod, ONLY: pi
USE jules_soil_mod, ONLY: dzsoil

IMPLICIT NONE

INTEGER, INTENT(IN) :: ft, nl
REAL(KIND=real_jlslsm), INTENT(IN) :: f_root(nl), root_mass
REAL(KIND=real_jlslsm) :: b_geom(nl)
INTEGER :: n
REAL(KIND=real_jlslsm) :: arg

b_geom(:) = 0.0
DO n = 1, nl
  IF ( f_root(n) <= 0.0 .OR. root_mass <= 0.0 ) CYCLE
  arg = rootc_density_pft(ft) * dzsoil(n) / (f_root(n) * root_mass)
  IF ( arg <= 1.0 ) CYCLE
  b_geom(n) = 4.0 * pi / LOG(arg) * (f_root(n) * root_mass)                   &
              / (rootc_density_pft(ft) * pi * root_radi_pft(ft)**2)
END DO

END FUNCTION rhizo_geometry

! *****************************************************************************
! Set the rhizosphere state of PFT ft at the points of index (physiol, after
! smc_ext): per layer the geometry b_geom (m-1), the soil psi (Pa, as
! soil_wp, i.e. after any bounds), and the soil hydraulic parameters.
! *****************************************************************************
SUBROUTINE rhizo_set_state( land_pts, nl, npft, ft, npts, index, b_geom,       &
                            psi_soil, sathh, bexp, k_sat )

USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: g
USE jules_vegetation_mod, ONLY: l_som_rhizo_hr

IMPLICIT NONE

INTEGER, INTENT(IN) :: land_pts, nl, npft, ft, npts, index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  b_geom(land_pts,nl), psi_soil(land_pts,nl), sathh(land_pts,nl),              &
  bexp(land_pts,nl), k_sat(land_pts,nl)

INTEGER :: j, l, n, it
REAL(KIND=rk) :: x, h_lo, h_hi, k_g
LOGICAL :: l_ok, l_rooted

IF ( .NOT. ALLOCATED(rz_on) ) THEN
  rz_nl = nl
  ALLOCATE( rz_on(land_pts,npft), rz_c(land_pts,nl,npft),                      &
            rz_sathh(land_pts,nl,npft), rz_hs(land_pts,nl,npft),               &
            rz_gs(land_pts,nl,npft), rz_h0(land_pts,npft),                     &
            rz_itab(land_pts,nl,npft) )
  rz_on(:,:) = .FALSE.
  rz_c(:,:,:) = 0.0
  rz_sathh(:,:,:) = 1.0
  rz_hs(:,:,:) = 0.0
  rz_gs(:,:,:) = 0.0
  rz_h0(:,:) = 0.0
  rz_itab(:,:,:) = 1
END IF
IF ( .NOT. ALLOCATED(rz_share) ) THEN
  ALLOCATE( rz_share(land_pts) )
  rz_share(:) = 1.0
END IF

rz_on(:,ft) = .FALSE.
DO j = 1, npts
  l = index(j)
  DO n = 1, nl
    rz_itab(l,n,ft)  = table_index(REAL(bexp(l,n), rk))
    rz_sathh(l,n,ft) = MAX(REAL(sathh(l,n), rk), TINY(1.0_rk))
    rz_c(l,n,ft)     = MAX(REAL(b_geom(l,n), rk), 0.0_rk)                      &
                       * MAX(REAL(k_sat(l,n), rk), 0.0_rk)                     &
                       * rz_sathh(l,n,ft) / m_h2o_rz
    rz_hs(l,n,ft)    = MAX(-REAL(psi_soil(l,n), rk), 0.0_rk)                   &
                       / (REAL(rho_water, rk) * REAL(g, rk))
    x = MIN(rz_hs(l,n,ft) / rz_sathh(l,n,ft), tab_x_max(rz_itab(l,n,ft)))
    rz_hs(l,n,ft)    = x * rz_sathh(l,n,ft)
    rz_gs(l,n,ft)    = g_of_x(rz_itab(l,n,ft), x)
  END DO
  rz_on(l,ft) = .TRUE.

  ! Root suction at zero flow: the wettest rooted layer, or with
  ! redistribution the balance sum c_i (G_s,i - G_i(h)) = 0, which lies
  ! between the wettest and driest rooted layers.
  h_lo = HUGE(1.0_rk)
  h_hi = 0.0
  l_rooted = .FALSE.
  DO n = 1, nl
    IF ( rz_c(l,n,ft) > 0.0 ) THEN
      h_lo = MIN(h_lo, rz_hs(l,n,ft))
      h_hi = MAX(h_hi, rz_hs(l,n,ft))
      l_rooted = .TRUE.
    END IF
  END DO
  IF ( .NOT. l_rooted ) THEN
    rz_h0(l,ft) = 0.0         ! no rooted layer: no supply (see rhizo_drop)
  ELSE
    rz_h0(l,ft) = h_lo
    IF ( l_som_rhizo_hr .AND. h_hi > h_lo ) THEN
      CALL solve_root( l, ft, 0.0_rk, h_lo, h_hi, rz_h0(l,ft),                 &
                       k_g, it, l_ok )
    END IF
  END IF
END DO

END SUBROUTINE rhizo_set_state

! *****************************************************************************
! .TRUE. if the nonlinear link applies at land point l for the current PFT.
! *****************************************************************************
LOGICAL FUNCTION rhizo_active( l )

USE jules_vegetation_mod, ONLY: l_som_rhizo_mfp

IMPLICIT NONE

INTEGER, INTENT(IN) :: l

rhizo_active = .FALSE.
IF ( .NOT. l_som_rhizo_mfp ) RETURN
IF ( .NOT. ALLOCATED(rz_on) ) RETURN
IF ( rz_ft < 1 ) RETURN
rhizo_active = rz_on(l,rz_ft)

END FUNCTION rhizo_active

! *****************************************************************************
! As rhizo_active for PFT ft (DESICA's commit runs after the PFT loop).
! *****************************************************************************
LOGICAL FUNCTION rhizo_active_ft( l, ft )

USE jules_vegetation_mod, ONLY: l_som_rhizo_mfp

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft

rhizo_active_ft = .FALSE.
IF ( .NOT. l_som_rhizo_mfp ) RETURN
IF ( .NOT. ALLOCATED(rz_on) ) RETURN
rhizo_active_ft = rz_on(l,ft)

END FUNCTION rhizo_active_ft

! *****************************************************************************
! Root water potential at zero flow (Pa) of point l, PFT ft.
! *****************************************************************************
REAL(KIND=real_jlslsm) FUNCTION rhizo_psi0( l, ft )

USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: g

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft

rhizo_psi0 = REAL(-REAL(rho_water, rk) * REAL(g, rk) * rz_h0(l,ft), real_jlslsm)

END FUNCTION rhizo_psi0

! *****************************************************************************
! Whole-plant soil link at point l, PFT ft, for ground-area transpiration
! e_g (mol m-2 s-1): the drop dpsi = psi_root(0) - psi_root(e_g) >= 0 (Pa)
! and the marginal conductance k_g = dE/dpsi_root (mol m-2 s-1 Pa-1).
! l_ok = .FALSE. if the soil cannot supply e_g (then dpsi = HUGE, k_g = 0).
! *****************************************************************************
SUBROUTINE rhizo_drop( l, ft, e_g, dpsi, k_g, l_ok )

USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: g
USE jules_print_mgr, ONLY: jules_message, jules_print

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: e_g
REAL(KIND=real_jlslsm), INTENT(OUT) :: dpsi, k_g
LOGICAL, INTENT(OUT) :: l_ok

INTEGER :: n, it
REAL(KIND=rk) :: h_hi, h_r, k_r

IF ( e_g <= 0.0 ) THEN
  dpsi = 0.0
  k_g  = REAL(k_soil( l, ft, rz_h0(l,ft) ), real_jlslsm)
  l_ok = .TRUE.
  RETURN
END IF

h_hi = 0.0
DO n = 1, rz_nl
  IF ( rz_c(l,n,ft) > 0.0 )                                                    &
    h_hi = MAX(h_hi, rz_sathh(l,n,ft) * tab_x_max(rz_itab(l,n,ft)))
END DO

CALL solve_root( l, ft, REAL(e_g, rk), rz_h0(l,ft), h_hi, h_r, k_r, it, l_ok )

n_solve    = n_solve + 1
n_iter_sum = n_iter_sum + it
n_iter_max = MAX(n_iter_max, it)
IF ( MOD(n_solve, n_report) == 0 ) THEN
  WRITE(jules_message,'(A,I0,A,F6.2,A,I0)') 'rhizo_mfp: ', n_solve,           &
    ' root solves, mean iterations ',                                          &
    REAL(n_iter_sum) / REAL(n_solve), ', max ', n_iter_max
  CALL jules_print('rhizo_mfp_mod', jules_message)
END IF

IF ( l_ok ) THEN
  dpsi = REAL(REAL(rho_water, rk) * REAL(g, rk) * (h_r - rz_h0(l,ft)),         &
              real_jlslsm)
  k_g  = REAL(k_r, real_jlslsm)
ELSE
  dpsi = HUGE(1.0_real_jlslsm)
  k_g  = 0.0
END IF

END SUBROUTINE rhizo_drop

! *****************************************************************************
! As rhizo_drop for a leaf path at point l (current PFT rz_ft) that carries
! the share rz_share(l) of the plant: e_path and k_path in the units of the
! path (as ksr_path).
! *****************************************************************************
SUBROUTINE rhizo_path_drop( l, e_path, dpsi, k_path, l_ok )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l
REAL(KIND=real_jlslsm), INTENT(IN) :: e_path
REAL(KIND=real_jlslsm), INTENT(OUT) :: dpsi, k_path
LOGICAL, INTENT(OUT) :: l_ok

REAL(KIND=real_jlslsm) :: s, k_g

s = REAL(rz_share(l), real_jlslsm)
IF ( s <= 0.0 ) THEN
  dpsi = 0.0
  k_path = 0.0
  l_ok = e_path <= 0.0
  RETURN
END IF
CALL rhizo_drop( l, rz_ft, e_path / s, dpsi, k_g, l_ok )
k_path = s * k_g

END SUBROUTINE rhizo_path_drop

! *****************************************************************************
! Marginal soil conductance of the current leaf path at zero flow (in the
! path's units, as ksr_path).
! *****************************************************************************
REAL(KIND=real_jlslsm) FUNCTION rhizo_k_zero( l )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l

rhizo_k_zero = REAL(rz_share(l) * k_soil( l, rz_ft, rz_h0(l,rz_ft) ),        &
                    real_jlslsm)

END FUNCTION rhizo_k_zero

! *****************************************************************************
! Per-layer uptake (mol m-2 s-1, >= 0) at point l, PFT ft, for ground-area
! transpiration e_g: c_i (G_s,i - G_i(h_root))^+. Zero everywhere if the
! soil cannot supply e_g or e_g <= 0.
! *****************************************************************************
SUBROUTINE rhizo_uptake( l, ft, e_g, q )

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=real_jlslsm), INTENT(IN) :: e_g
REAL(KIND=real_jlslsm), INTENT(OUT) :: q(:)

INTEGER :: n, it
REAL(KIND=rk) :: h_hi, h_r, k_g
LOGICAL :: l_ok

q(:) = 0.0
IF ( e_g <= 0.0 ) RETURN
h_hi = 0.0
DO n = 1, rz_nl
  IF ( rz_c(l,n,ft) > 0.0 )                                                    &
    h_hi = MAX(h_hi, rz_sathh(l,n,ft) * tab_x_max(rz_itab(l,n,ft)))
END DO
CALL solve_root( l, ft, REAL(e_g, rk), rz_h0(l,ft), h_hi, h_r, k_g, it, l_ok )
IF ( .NOT. l_ok ) RETURN
DO n = 1, rz_nl
  IF ( rz_c(l,n,ft) <= 0.0 ) CYCLE
  q(n) = REAL(MAX(rz_c(l,n,ft) * ( rz_gs(l,n,ft)                              &
             - g_of_x(rz_itab(l,n,ft), h_r / rz_sathh(l,n,ft)) ), 0.0_rk),     &
             real_jlslsm)
END DO

END SUBROUTINE rhizo_uptake

! =============================================================================
! Internal routines
! =============================================================================

! Root suction h_r (m) with F(h_r) = sum_i c_i d_i(h_r) = target, where
! d_i = G_s,i - G_i(h_r / sathh_i) (clipped at 0 without redistribution),
! searched in [h_lo, h_hi] with F(h_lo) <= target. F increases with h_r,
! and is concave where the set of supplying layers does not change, so
! Newton from the left converges monotonically; a bracket and log-space
! bisection guard the kinks and the flat dry tail. k_g = F'(h_r) / (rho g).
SUBROUTINE solve_root( l, ft, target, h_lo_in, h_hi_in, h_r, k_g, it, l_ok )

USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: g

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=rk), INTENT(IN) :: target, h_lo_in, h_hi_in
REAL(KIND=rk), INTENT(OUT) :: h_r, k_g
INTEGER, INTENT(OUT) :: it
LOGICAL, INTENT(OUT) :: l_ok

REAL(KIND=rk) :: lo, hi, h, hn, f, fp, f_hi, tol

it = 0
l_ok = .TRUE.
lo = MAX(h_lo_in, 0.0_rk)
hi = MAX(h_hi_in, lo)

CALL eval_f( l, ft, hi, f_hi, fp )
IF ( f_hi < target ) THEN
  ! Beyond what the soil can supply.
  l_ok = .FALSE.
  h_r = hi
  k_g = 0.0
  RETURN
END IF

! Converged when E is met to 1e-6 (relative), or, for the zero-flow
! balance (target 0), to 1e-9 of the supply; or when the Newton step is
! below 1e-9 of h (the table's interpolation noise, relative to the supply,
! can otherwise stall the residual test in wet soil, where the supply is
! 1e3-1e5 x E).
IF ( target > 0.0 ) THEN
  tol = 1.0e-6_rk * target
ELSE
  tol = 1.0e-9_rk * MAX(f_hi, TINY(1.0_rk))
END IF
h = lo
CALL eval_f( l, ft, h, f, fp )
f = f - target
DO it = 1, max_iter
  IF ( ABS(f) <= tol ) EXIT
  IF ( fp > 0.0 ) THEN
    hn = h - f / fp
  ELSE
    hn = hi
  END IF
  IF ( hn <= lo .OR. hn >= hi ) THEN
    IF ( lo > 0.0 .AND. hi > 4.0 * lo ) THEN
      hn = SQRT(lo * hi)
    ELSE
      hn = 0.5 * (lo + hi)
    END IF
  ELSE IF ( ABS(hn - h) <= 1.0e-9_rk * MAX(h, 1.0e-3_rk) ) THEN
    h = hn
    CALL eval_f( l, ft, h, f, fp )
    EXIT
  END IF
  h = hn
  CALL eval_f( l, ft, h, f, fp )
  f = f - target
  IF ( f > 0.0 ) THEN
    hi = h
  ELSE
    lo = h
  END IF
  IF ( hi - lo <= 1.0e-12_rk * hi ) EXIT
END DO
it = MIN(it, max_iter)
h_r = h
k_g = fp / (REAL(rho_water, rk) * REAL(g, rk))

END SUBROUTINE solve_root

! F(h) = sum_i c_i d_i and dF/dh = sum_i c_i K_i/K_sat,i / sathh_i over the
! supplying layers (all rooted layers with redistribution). Without
! redistribution a layer joins once the root is at least as dry as it (the
! right-hand derivative, so Newton can leave the zero-flow root).
SUBROUTINE eval_f( l, ft, h, f, fp )

USE jules_vegetation_mod, ONLY: l_som_rhizo_hr

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=rk), INTENT(IN) :: h
REAL(KIND=rk), INTENT(OUT) :: f, fp

INTEGER :: n, itab
REAL(KIND=rk) :: x, gx, dgdu, d

f  = 0.0
fp = 0.0
DO n = 1, rz_nl
  IF ( rz_c(l,n,ft) <= 0.0 ) CYCLE
  IF ( .NOT. l_som_rhizo_hr .AND. h < rz_hs(l,n,ft) ) CYCLE
  itab = rz_itab(l,n,ft)
  x = h / rz_sathh(l,n,ft)
  CALL g_and_kr( itab, x, gx, dgdu )
  d = rz_gs(l,n,ft) - gx
  IF ( .NOT. l_som_rhizo_hr ) d = MAX(d, 0.0_rk)
  f = f + rz_c(l,n,ft) * d
  ! dG/dh = -kr / sathh, so dd/dh = kr / sathh = -dgdu / (x sathh), with
  ! dgdu the slope of the same interpolant (exactly consistent with G).
  IF ( x > x_lo ) THEN
    fp = fp - rz_c(l,n,ft) * dgdu / (x * rz_sathh(l,n,ft))
  ELSE
    fp = fp - rz_c(l,n,ft) * tab_dg(1,itab) / (x_lo * rz_sathh(l,n,ft))
  END IF
END DO

END SUBROUTINE eval_f

! Marginal soil conductance (mol m-2 s-1 Pa-1, ground) at root suction h:
! with no redistribution only the layers at least as wet as the root.
REAL(KIND=rk) FUNCTION k_soil( l, ft, h )

USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: g
USE jules_vegetation_mod, ONLY: l_som_rhizo_hr

IMPLICIT NONE

INTEGER, INTENT(IN) :: l, ft
REAL(KIND=rk), INTENT(IN) :: h

INTEGER :: n
REAL(KIND=rk) :: x, gx, dgdu

k_soil = 0.0
DO n = 1, rz_nl
  IF ( rz_c(l,n,ft) <= 0.0 ) CYCLE
  IF ( .NOT. l_som_rhizo_hr .AND. h < rz_hs(l,n,ft) ) CYCLE
  x = h / rz_sathh(l,n,ft)
  IF ( x <= x_lo ) THEN
    k_soil = k_soil - rz_c(l,n,ft) * tab_dg(1,rz_itab(l,n,ft))                 &
             / (x_lo * rz_sathh(l,n,ft))
  ELSE
    CALL g_and_kr( rz_itab(l,n,ft), x, gx, dgdu )
    k_soil = k_soil - rz_c(l,n,ft) * dgdu / (x * rz_sathh(l,n,ft))
  END IF
END DO
k_soil = k_soil / (REAL(rho_water, rk) * REAL(g, rk))

END FUNCTION k_soil

! Relative conductivity K/K_sat at x = h / sathh, as JULES: van Genuchten
! (hh_from_sthu inverted, hyd_con_vg with its Se clamps) or Brooks-Corey
! (hyd_con_ch).
REAL(KIND=rk) FUNCTION kr_of_x( b, x )

USE jules_soil_mod, ONLY: l_vg_soil

IMPLICIT NONE

REAL(KIND=rk), INTENT(IN) :: b, x

REAL(KIND=rk), PARAMETER :: l_wag = 0.5_rk, se_min = 0.05_rk,        &
                                     se_max = 0.95_rk
REAL(KIND=rk) :: se, sd, k

IF ( x <= 0.0 ) THEN
  kr_of_x = 1.0
  RETURN
END IF
IF ( l_vg_soil ) THEN
  se = ( 1.0 + x**((b + 1.0) / b) )**(-1.0 / (b + 1.0))
  sd = MIN(MAX(se, se_min), se_max)
  k  = sd**l_wag * ( 1.0 - (1.0 - sd**(b + 1.0))**(1.0 / (b + 1.0)) )**2
  IF ( se < se_min ) THEN
    k = k / se_min * se
  ELSE IF ( se > se_max ) THEN
    k = k + (1.0 - k) / (1.0 - se_max) * (se - se_max)
  END IF
  kr_of_x = k
ELSE
  IF ( x <= 1.0 ) THEN
    kr_of_x = 1.0
  ELSE
    se = x**(-1.0 / b)
    kr_of_x = se**(2.0 * b + 3.0)
  END IF
END IF

END FUNCTION kr_of_x

! x = h / sathh at the driest state (Se = sthu_min_rz).
REAL(KIND=rk) FUNCTION x_driest( b )

USE jules_soil_mod, ONLY: l_vg_soil

IMPLICIT NONE

REAL(KIND=rk), INTENT(IN) :: b

IF ( l_vg_soil ) THEN
  x_driest = ( sthu_min_rz**(-b - 1.0) - 1.0 )**(b / (b + 1.0))
ELSE
  x_driest = sthu_min_rz**(-b)
END IF
x_driest = MIN(x_driest, x_cap)

END FUNCTION x_driest

! Index of the G table for exponent b, building it on first use:
! G(x) = integral of kr dx' from x to x_driest(b), on a uniform grid in
! u = ln x, integrated from the dry end with Simpson's rule per interval.
INTEGER FUNCTION table_index( b )

USE ereport_mod, ONLY: ereport

IMPLICIT NONE

REAL(KIND=rk), INTENT(IN) :: b

INTEGER :: i, k, errcode
REAL(KIND=rk) :: u_hi, du, u0, u1, um, x0, x1, xm

DO i = 1, n_tab
  IF ( tab_b(i) == b ) THEN
    table_index = i
    RETURN
  END IF
END DO

!$OMP CRITICAL (rhizo_mfp_table)
IF ( n_tab >= n_tab_max ) THEN
  errcode = 101
  CALL ereport('rhizo_mfp_mod:table_index', errcode,                           &
               'too many distinct soil b values for the rhizosphere tables')
END IF
n_tab = n_tab + 1
i = n_tab
tab_b(i)     = b
tab_x_max(i) = x_driest(b)
tab_u_lo(i)  = LOG(x_lo)
u_hi         = LOG(MAX(tab_x_max(i), 10.0_rk * x_lo))
du           = (u_hi - tab_u_lo(i)) / REAL(n_g - 1)
tab_du(i)    = du
tab_g(n_g,i) = 0.0
x1 = EXP(u_hi)
tab_dg(n_g,i) = -kr_of_x(b, x1) * x1
DO k = n_g - 1, 1, -1
  u0 = tab_u_lo(i) + REAL(k - 1) * du
  u1 = u0 + du
  um = 0.5 * (u0 + u1)
  x0 = EXP(u0)
  x1 = EXP(u1)
  xm = EXP(um)
  ! dG/du = -kr(x) x
  tab_dg(k,i) = -kr_of_x(b, x0) * x0
  tab_g(k,i)  = tab_g(k+1,i) + du / 6.0                                        &
                * ( kr_of_x(b, x0) * x0 + 4.0 * kr_of_x(b, xm) * xm            &
                    + kr_of_x(b, x1) * x1 )
END DO
!$OMP END CRITICAL (rhizo_mfp_table)
table_index = i

END FUNCTION table_index

! G(x) from table itab.
REAL(KIND=rk) FUNCTION g_of_x( itab, x )

IMPLICIT NONE

INTEGER, INTENT(IN) :: itab
REAL(KIND=rk), INTENT(IN) :: x
REAL(KIND=rk) :: dgdu

CALL g_and_kr( itab, x, g_of_x, dgdu )

END FUNCTION g_of_x

! G(x) and dG/du by cubic Hermite interpolation in u = ln x (dG/du is the
! interpolant's own slope, so K = -dG/du / x is consistent with G). Wetter
! than the table: G grows linearly with K = K(x_lo); drier: G = 0.
SUBROUTINE g_and_kr( itab, x, gx, dgdu )

IMPLICIT NONE

INTEGER, INTENT(IN) :: itab
REAL(KIND=rk), INTENT(IN) :: x
REAL(KIND=rk), INTENT(OUT) :: gx, dgdu

INTEGER :: k
REAL(KIND=rk) :: u, t, t2, t3, du, h00, h10, h01, h11

IF ( x <= x_lo ) THEN
  gx   = tab_g(1,itab) - tab_dg(1,itab) / x_lo * (x_lo - MAX(x, 0.0_rk))
  dgdu = tab_dg(1,itab) / x_lo * MAX(x, 0.0_rk)
  RETURN
END IF
IF ( x >= tab_x_max(itab) ) THEN
  gx   = 0.0
  dgdu = 0.0
  RETURN
END IF
du = tab_du(itab)
u  = LOG(x)
t  = (u - tab_u_lo(itab)) / du
k  = MIN(MAX(INT(t) + 1, 1), n_g - 1)
t  = t - REAL(k - 1)
t2 = t * t
t3 = t2 * t
h00 = 2.0 * t3 - 3.0 * t2 + 1.0
h10 = t3 - 2.0 * t2 + t
h01 = -2.0 * t3 + 3.0 * t2
h11 = t3 - t2
gx = h00 * tab_g(k,itab) + h10 * du * tab_dg(k,itab)                           &
     + h01 * tab_g(k+1,itab) + h11 * du * tab_dg(k+1,itab)
! Slope of the interpolant, dG/du = -kr x (no powers: K from the table).
dgdu = ( (6.0_rk * t2 - 6.0_rk * t) * (tab_g(k,itab) - tab_g(k+1,itab)) / du  &
         + (3.0_rk * t2 - 4.0_rk * t + 1.0_rk) * tab_dg(k,itab)               &
         + (3.0_rk * t2 - 2.0_rk * t) * tab_dg(k+1,itab) )

END SUBROUTINE g_and_kr

END MODULE rhizo_mfp_mod
