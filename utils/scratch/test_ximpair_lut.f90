! Unit check of the capped lookup-table solvers for the memory xylem
! impairment model (branch ximpair_lut) against the Newton-Raphson solvers,
! at fixed caps fr = 1, 0.7, 0.3, 0.05, over a grid of root-zone psi and E.
!   (1) segments: leaf_psi_segments_lut_jls vs leaf_psi_segments_jls
!       (kcap_frac = fr; leaf_k, leaf_k_intact, stem_psi compared)
!   (2) no segments: leaf_psi_impaired_memory with som_psi_solver = 3
!       (leaf_psi_capped_lut_jls) vs 2 (Newton-Raphson)
! FR-Pue curve parameters (as printed by best_seg rec3's run.log).
! Build: see utils/scratch/build_test_ximpair_lut.sh. Prints max |dpsi| (MPa)
! and max |dk|/kmax over feasible samples (k > kcrit in both) and over all.
! MDK, 18/10/26
!   (3) whole-plant table down to the failure point k = sl_kfail_frac * kmax
!       (kcrit if 0), Weibull (pft 1) and SOX (pft 2, same P50/P88): table
!       and Newton psi_leaf against an exact forward integral, the failure
!       flag (k > k_fail) and E_crit found by bisection on k(E) as
!       stom_opt_supply_loss does.
! Usage: test_ximpair_lut [sl_kfail_frac]  (default 0; the tables are built
! once, so run once per value).
! MDK, 20/10/26
PROGRAM test_ximpair_lut

USE um_types, ONLY: real_jlslsm
USE pftparm, ONLY: seg_kfac, conductance_b_seg, conductance_c_seg,             &
                   conductance_b_pft, conductance_c_pft, pft_conductance_model, &
                   kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut,               &
                                psi_solver_newton, sl_kfail_frac
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_segments_jls,                  &
                                       leaf_psi_segments_lut_jls, ksr_path,    &
                                       supply_lut_psi, supply_lut_f
USE xylem_hydraulics_jls_mod, ONLY: leaf_psi_jls
USE xylem_impairment_memory_mod, ONLY: leaf_psi_impaired_memory

IMPLICIT NONE

INTEGER, PARAMETER :: rk = real_jlslsm
INTEGER, PARAMETER :: n_e = 200, n_r = 13, n_fr = 4
REAL(KIND=rk), PARAMETER :: fr_list(n_fr) = [1.0_rk, 0.7_rk, 0.3_rk, 0.05_rk]
REAL(KIND=rk), PARAMETER :: kmax0 = 2.4e-9_rk, kcrit0 = 0.05_rk * kmax0
                            ! Big-leaf kmax ~ kmax_pft * LAI (2) at FR-Pue.
REAL(KIND=rk), PARAMETER :: e_max = 1.2e-2_rk
                            ! mol m-2 s-1: well past what the path supplies.
INTEGER :: veg_index(1), open_index(1), i, ir, ifr, nfeas
REAL(KIND=rk) :: e(n_e,1), psi_r(1), kmax(1), kcrit(1), kfr(1), kcap(1),      &
                 b(1), c(1)
REAL(KIND=rk) :: psi_n(n_e,1), k_n(n_e,1), ki_n(n_e,1), st_n(n_e,1)
REAL(KIND=rk) :: psi_t(n_e,1), k_t(n_e,1), ki_t(n_e,1), st_t(n_e,1)
REAL(KIND=rk) :: dpsi_f, dk_f, dki_f, dst_f, dpsi_a, dk_a, dmax_rel
! Section (3). MDK, 20/10/26
INTEGER, PARAMETER :: n_pl = 400, n_simp = 20000, n_r3 = 4
REAL(KIND=rk), PARAMETER :: psi_r3(n_r3) = [-0.1e6_rk, -1.0e6_rk, -3.0e6_rk, &
                                            -5.0e6_rk]
INTEGER :: ip, k, nflag, nfeas3, it
CHARACTER(LEN=32) :: arg
REAL(KIND=rk) :: f_fail, psi_crit, psi_l, e3(1,1), e_ex, psi_tt, k_tt,       &
                 psi_nn(1,1), k_nn(1,1), dpsi_t3, dpsi_n3, ecrit_ex, e_lo,   &
                 e_hi, e_mid, ecrit_err, psi_bot, kmx3(1), kcr3(1), b3(1), &
                 c3(1), pr3(1)

sl_kfail_frac = 0.0_rk
IF ( COMMAND_ARGUMENT_COUNT() >= 1 ) THEN
  CALL GET_COMMAND_ARGUMENT(1, arg)
  READ(arg,*) sl_kfail_frac
END IF
WRITE(*,'(A,ES10.3)') 'sl_kfail_frac = ', sl_kfail_frac

ALLOCATE( seg_kfac(2,3), conductance_b_seg(2,3), conductance_c_seg(2,3),     &
          conductance_b_pft(2), conductance_c_pft(2), pft_conductance_model(2),&
          kcrit_fractional_loss(2) )
seg_kfac(1,:)          = [2.0_rk, 4.0_rk, 4.0_rk]
conductance_b_seg(1,:) = [-2924382.75_rk, -6999276.0_rk, -6999276.0_rk]
conductance_c_seg(1,:) = [1.58152008_rk, 7.37631893_rk, 7.37631893_rk]
conductance_b_pft(1)   = -6999276.0_rk
conductance_c_pft(1)   = 7.37631893_rk
pft_conductance_model(1) = 1
kcrit_fractional_loss(1) = 0.95_rk
! pft 2: SOX with FR-Pue's P50 = -6.66, P88 = -7.75 MPa (pftparm_io_mod:
! b = P50, c = ln(1/0.12 - 1) / ln(P88/P50)). Segments as pft 1 (unused).
! MDK, 20/10/26
seg_kfac(2,:) = seg_kfac(1,:)
conductance_b_seg(2,:) = conductance_b_seg(1,:)
conductance_c_seg(2,:) = conductance_c_seg(1,:)
conductance_b_pft(2)   = -6.66e6_rk
conductance_c_pft(2)   = LOG(1.0_rk / 0.12_rk - 1.0_rk) / LOG(7.75_rk / 6.66_rk)
pft_conductance_model(2) = 2
kcrit_fractional_loss(2) = 0.95_rk
ALLOCATE( ksr_path(1) )
ksr_path(:) = -1.0_rk   ! no soil link (l_som_rhizo_series is off)

veg_index(1) = 1
open_index(1) = 1
kmax(1)  = kmax0
kcrit(1) = kcrit0
b(1) = conductance_b_pft(1)
c(1) = conductance_c_pft(1)
DO i = 1, n_e
  e(i,1) = e_max * REAL(i - 1, rk) / REAL(n_e - 1, rk)
END DO

WRITE(*,'(A)') '(1) segments: table vs Newton (4 it.), psi_root -0.1..-6 MPa'
WRITE(*,'(A)') '  fr    n_feas  max|dpsi_leaf| max|dk|/kmax max|dk_int|/kmax' // &
               ' max|dpsi_stem|  [feasible] | max|dpsi| max|dk|/kmax [all]'
DO ifr = 1, n_fr
  dpsi_f = 0.0; dk_f = 0.0; dki_f = 0.0; dst_f = 0.0; dpsi_a = 0.0; dk_a = 0.0
  nfeas = 0
  kfr(1) = fr_list(ifr)
  DO ir = 1, n_r
    psi_r(1) = -0.1e6_rk - 5.9e6_rk * REAL(ir - 1, rk) / REAL(n_r - 1, rk)
    CALL leaf_psi_segments_jls( 1, n_e, 1, 1, veg_index, open_index, e,      &
                                psi_r, kmax, kcrit, psi_n, k_n,               &
                                kcap_frac = kfr, leaf_k_intact = ki_n,        &
                                stem_psi = st_n )
    CALL leaf_psi_segments_lut_jls( 1, n_e, 1, 1, veg_index, open_index, e,  &
                                    psi_r, kmax, psi_t, k_t,                  &
                                    kcap_frac = kfr, leaf_k_intact = ki_t,    &
                                    stem_psi = st_t )
    DO i = 1, n_e
      dpsi_a = MAX(dpsi_a, ABS(psi_t(i,1) - psi_n(i,1)))
      dk_a   = MAX(dk_a, ABS(k_t(i,1) - k_n(i,1)))
      IF ( k_n(i,1) > kcrit0 .AND. k_t(i,1) > kcrit0 ) THEN
        nfeas = nfeas + 1
        dpsi_f = MAX(dpsi_f, ABS(psi_t(i,1) - psi_n(i,1)))
        dk_f   = MAX(dk_f, ABS(k_t(i,1) - k_n(i,1)))
        dki_f  = MAX(dki_f, ABS(ki_t(i,1) - ki_n(i,1)))
        dst_f  = MAX(dst_f, ABS(st_t(i,1) - st_n(i,1)))
      END IF
    END DO
  END DO
  WRITE(*,'(F6.2,I8,4ES14.3,A,2ES12.3)') fr_list(ifr), nfeas, dpsi_f * 1.0e-6,&
        dk_f / kmax0, dki_f / kmax0, dst_f * 1.0e-6, '   |', dpsi_a * 1.0e-6,  &
        dk_a / kmax0
END DO

WRITE(*,'(A)') '(2) no segments: table (solver 3) vs Newton (solver 2, 20 it.)'
WRITE(*,'(A)') '  fr    n_feas  max|dpsi_leaf| max|dk|/kmax  [feasible] |' // &
               ' max|dpsi| max|dk|/kmax [all]'
DO ifr = 1, n_fr
  dpsi_f = 0.0; dk_f = 0.0; dpsi_a = 0.0; dk_a = 0.0
  nfeas = 0
  kcap(1) = fr_list(ifr) * kmax0
  DO ir = 1, n_r
    psi_r(1) = -0.1e6_rk - 5.9e6_rk * REAL(ir - 1, rk) / REAL(n_r - 1, rk)
    som_psi_solver = psi_solver_newton
    CALL leaf_psi_impaired_memory( 1, n_e, 1, 1, veg_index, open_index, e,   &
                                   psi_r, kmax, kcap, kcrit, b, c, psi_n, k_n )
    som_psi_solver = psi_solver_lut
    CALL leaf_psi_impaired_memory( 1, n_e, 1, 1, veg_index, open_index, e,   &
                                   psi_r, kmax, kcap, kcrit, b, c, psi_t, k_t )
    DO i = 1, n_e
      dpsi_a = MAX(dpsi_a, ABS(psi_t(i,1) - psi_n(i,1)))
      dk_a   = MAX(dk_a, ABS(k_t(i,1) - k_n(i,1)))
      IF ( k_n(i,1) > kcrit0 .AND. k_t(i,1) > kcrit0 ) THEN
        nfeas = nfeas + 1
        dpsi_f = MAX(dpsi_f, ABS(psi_t(i,1) - psi_n(i,1)))
        dk_f   = MAX(dk_f, ABS(k_t(i,1) - k_n(i,1)))
      END IF
    END DO
  END DO
  WRITE(*,'(F6.2,I8,2ES14.3,A,2ES12.3)') fr_list(ifr), nfeas,                &
        dpsi_f * 1.0e-6, dk_f / kmax0, '   |', dpsi_a * 1.0e-6, dk_a / kmax0
END DO

! Largest relative psi error where it matters most: report the worst sample
! of the segment case at fr = 0.3 (diagnostic of where the differences are).
dmax_rel = 0.0
kfr(1) = 0.3_rk
psi_r(1) = -3.0e6_rk
CALL leaf_psi_segments_jls( 1, n_e, 1, 1, veg_index, open_index, e, psi_r,    &
                            kmax, kcrit, psi_n, k_n, kcap_frac = kfr )
CALL leaf_psi_segments_lut_jls( 1, n_e, 1, 1, veg_index, open_index, e,      &
                                psi_r, kmax, psi_t, k_t, kcap_frac = kfr )
WRITE(*,'(A)') 'fr = 0.3, psi_root = -3 MPa: E (mmol m-2 s-1), psi Newton,' // &
               ' psi table (MPa), k Newton, k table (/kmax)'
DO i = 1, n_e, 10
  WRITE(*,'(F8.3,2F10.4,2F9.4)') e(i,1) * 1.0e3, psi_n(i,1) * 1.0e-6,        &
        psi_t(i,1) * 1.0e-6, k_n(i,1) / kmax0, k_t(i,1) / kmax0
END DO

! (3) Whole-plant table down to the failure point. For each psi_root, psi_l
! runs from psi_root down to 1.3 x (psi_crit - psi_root) below psi_root,
! psi_crit where f = f_fail. E_exact(psi_l) = kmax * integral(f, psi_l,
! psi_root) by composite Simpson (n_simp cells). Table (solver 3) and Newton
! (solver 2) psi_leaf at E_exact are compared with psi_l where psi_l >=
! psi_crit; the failure flag k > k_fail with the exact one everywhere.
! E_crit: bisection on k(E) > k_fail with the table (as
! stom_opt_supply_loss) against E_exact(psi_crit).
! MDK, 20/10/26
WRITE(*,'(A)') '(3) whole plant, failure at f_fail = sl_kfail_frac (kcrit if 0)'
WRITE(*,'(A)') ' model  psi_r  psi_crit  n_feas max|dpsi|tab max|dpsi|Newt' // &
               ' flag_miss  E_crit rel.err  psi_tab(E_crit_exact)'
f_fail = 1.0_rk - kcrit_fractional_loss(1)
IF ( sl_kfail_frac > 0.0 ) f_fail = sl_kfail_frac
DO ip = 1, 2
  IF ( pft_conductance_model(ip) == 2 ) THEN
    psi_crit = conductance_b_pft(ip) * (1.0_rk / f_fail - 1.0_rk)            &
               **(1.0_rk / conductance_c_pft(ip))
  ELSE
    psi_crit = conductance_b_pft(ip) * (-LOG(f_fail))                        &
               **(1.0_rk / conductance_c_pft(ip))
  END IF
  kmx3(1) = kmax0
  kcr3(1) = f_fail * kmax0
  b3(1) = conductance_b_pft(ip)
  c3(1) = conductance_c_pft(ip)
  DO ir = 1, n_r3
    pr3(1) = psi_r3(ir)
    dpsi_t3 = 0.0; dpsi_n3 = 0.0; nflag = 0; nfeas3 = 0
    DO k = 1, n_pl
      psi_l = pr3(1) + 1.3_rk * (psi_crit - pr3(1)) * REAL(k, rk) / REAL(n_pl, rk)
      e_ex = kmax0 * s_exact(ip, psi_l, pr3(1))
      psi_tt = supply_lut_psi(ip, pr3(1), e_ex / kmax0)
      k_tt = kmax0 * supply_lut_f(ip, psi_tt)
      IF ( (k_tt > f_fail * kmax0) .NEQV. (psi_l > psi_crit) ) nflag = nflag + 1
      IF ( psi_l >= psi_crit ) THEN
        nfeas3 = nfeas3 + 1
        dpsi_t3 = MAX(dpsi_t3, ABS(psi_tt - psi_l))
        e3(1,1) = e_ex
        som_psi_solver = psi_solver_newton
        CALL leaf_psi_jls( ip, 1, 1, 1, veg_index, open_index, e3, pr3, kmx3, &
                           kcr3, b3, c3, psi_nn, k_nn )
        dpsi_n3 = MAX(dpsi_n3, ABS(psi_nn(1,1) - psi_l))
      END IF
    END DO
    ! E_crit by bisection on the table's k(E), as stom_opt_supply_loss.
    ecrit_ex = kmax0 * s_exact(ip, psi_crit, pr3(1))
    e_lo = 0.0; e_hi = 3.0_rk * ecrit_ex
    DO it = 1, 80
      e_mid = 0.5_rk * (e_lo + e_hi)
      psi_tt = supply_lut_psi(ip, pr3(1), e_mid / kmax0)
      IF ( kmax0 * supply_lut_f(ip, psi_tt) > f_fail * kmax0 ) THEN
        e_lo = e_mid
      ELSE
        e_hi = e_mid
      END IF
    END DO
    ecrit_err = e_lo / ecrit_ex - 1.0_rk
    psi_bot = supply_lut_psi(ip, pr3(1), ecrit_ex / kmax0)
    WRITE(*,'(I4,2F9.3,I7,2ES14.3,I8,ES14.3,F12.4)') pft_conductance_model(ip),&
          pr3(1) * 1.0e-6, psi_crit * 1.0e-6, nfeas3, dpsi_t3 * 1.0e-6,       &
          dpsi_n3 * 1.0e-6, nflag, ecrit_err, psi_bot * 1.0e-6
  END DO
END DO

CONTAINS

  ! integral(f, psi_lo, psi_hi) for pft ip, composite Simpson, n_simp cells.
  ! MDK, 20/10/26
  REAL(KIND=rk) FUNCTION s_exact(ip, psi_lo, psi_hi)
  INTEGER, INTENT(IN) :: ip
  REAL(KIND=rk), INTENT(IN) :: psi_lo, psi_hi
  INTEGER :: m
  REAL(KIND=rk) :: hh, acc
  hh = (psi_hi - psi_lo) / REAL(n_simp, rk)
  acc = f_ex(ip, psi_lo) + f_ex(ip, psi_hi)
  DO m = 1, n_simp - 1
    IF ( MOD(m, 2) == 1 ) THEN
      acc = acc + 4.0_rk * f_ex(ip, psi_lo + m * hh)
    ELSE
      acc = acc + 2.0_rk * f_ex(ip, psi_lo + m * hh)
    END IF
  END DO
  s_exact = acc * hh / 3.0_rk
  END FUNCTION s_exact

  REAL(KIND=rk) FUNCTION f_ex(ip, psi)
  INTEGER, INTENT(IN) :: ip
  REAL(KIND=rk), INTENT(IN) :: psi
  IF ( pft_conductance_model(ip) == 2 ) THEN
    f_ex = 1.0_rk / (1.0_rk + (ABS(psi / conductance_b_pft(ip)))             &
                              **conductance_c_pft(ip))
  ELSE
    f_ex = EXP(-(ABS(psi / conductance_b_pft(ip)))**conductance_c_pft(ip))
  END IF
  END FUNCTION f_ex

END PROGRAM test_ximpair_lut
