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
PROGRAM test_ximpair_lut

USE um_types, ONLY: real_jlslsm
USE pftparm, ONLY: seg_kfac, conductance_b_seg, conductance_c_seg,             &
                   conductance_b_pft, conductance_c_pft, pft_conductance_model, &
                   kcrit_fractional_loss
USE jules_vegetation_mod, ONLY: som_psi_solver, psi_solver_lut,               &
                                psi_solver_newton
USE xylem_hydraulics_CW_jls_mod, ONLY: leaf_psi_segments_jls,                  &
                                       leaf_psi_segments_lut_jls, ksr_path
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

ALLOCATE( seg_kfac(1,3), conductance_b_seg(1,3), conductance_c_seg(1,3),     &
          conductance_b_pft(1), conductance_c_pft(1), pft_conductance_model(1),&
          kcrit_fractional_loss(1) )
seg_kfac(1,:)          = [2.0_rk, 4.0_rk, 4.0_rk]
conductance_b_seg(1,:) = [-2924382.75_rk, -6999276.0_rk, -6999276.0_rk]
conductance_c_seg(1,:) = [1.58152008_rk, 7.37631893_rk, 7.37631893_rk]
conductance_b_pft(1)   = -6999276.0_rk
conductance_c_pft(1)   = 7.37631893_rk
pft_conductance_model(1) = 1
kcrit_fractional_loss(1) = 0.95_rk
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

END PROGRAM test_ximpair_lut
