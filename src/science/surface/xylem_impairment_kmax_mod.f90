! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE xylem_impairment_kmax_mod

! kmax-impaired xylem impairment model (pft_xylem_impairment_model = 1).
!
! The plant's kmax falls with the conductance the tree has at its leaves
! and in its root zone, and the vulnerability curve is refitted to the lowered
! kmax, so every point has its own kmax, b and c. Damage is not recovered.
!
! Each timestep (once, after the humidity iteration) the target kmax is
!
!   k_new = (w_l * k_leaf + (1 - w_l) * k_root) / (kmax' - k_thr) * kmax'
!
! with k_leaf and k_root the conductances of the impaired curve at the
! leaf and root zone water potentials, w_l = ximpair_leaf_weight, kmax' the
! current impaired kmax and k_thr = ximpair_threshold. The new kmax is
!
!   kmax'' = w_n * k_new + (1 - w_n) * kmax',  w_n = ximpair_new_kmax_weight,
!
! limited to [kcrit, kmax'].
!
! The curve is refitted to the impaired kmax so that it passes through the
! intact curve at b'. For the cumulative Weibull curve k = kmax * EXP(-(psi
! / b)^c),
!
!   b' = b * (1 - LN(kmax' / kmax))^(1 / c),
!   c' = c * (b' / b)^c,
!
! so the impaired curve meets the intact one at b' with the same slope and
! never exceeds it. For the SOX curve k = kmax / (1 + (psi / b)^c),
!
!   b' = b * (2 * kmax / kmax' - 1)^(1 / c),
!   c' = c * (kmax / kmax') * (b' / b)^c.
!
! The refit depends only on kmax' and the PFT's curve, so kmax' is the only
! state of the model (kmax_impaired_pft). It acts with the big-leaf and
! two-leaf canopies, the Newton leaf water potential solver and no plant
! segments.
!
! MDK, 21/11/26

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='XYLEM_IMPAIRMENT_KMAX_MOD'

PRIVATE

PUBLIC :: update_xylem_impairment_kmax,                                        &
          xylem_kmax_refit_curve

CONTAINS

! ---------------------------------------------------------------------
! The vulnerability curve (b in Pa, c) refitted to an impaired kmax, at
! each of n points.
!
! MDK, 21/11/26
! ---------------------------------------------------------------------
SUBROUTINE xylem_kmax_refit_curve( pft, n, impaired_k_max,                     &
                                   impaired_conductance_b,                     &
                                   impaired_conductance_c )

USE pftparm, ONLY: pft_conductance_model, kmax_pft, conductance_b_pft,        &
                   conductance_c_pft
USE jules_vegetation_mod, ONLY: CW_conductance, SOX_conductance

USE ereport_mod, ONLY: ereport

INTEGER, INTENT(IN) :: pft
INTEGER, INTENT(IN) :: n

! The impaired kmax of each point (mol m-2 s-1 Pa-1, leaf basis).
REAL(KIND=real_jlslsm), INTENT(IN) :: impaired_k_max(n)

! The refitted curve of each point.
REAL(KIND=real_jlslsm), INTENT(OUT) :: impaired_conductance_b(n)
REAL(KIND=real_jlslsm), INTENT(OUT) :: impaired_conductance_c(n)

! The impaired share of kmax, held in (0, 1].
REAL(KIND=real_jlslsm) :: ratio(n)

INTEGER :: errcode

CHARACTER(LEN=*), PARAMETER :: RoutineName='XYLEM_KMAX_REFIT_CURVE'

ratio(:) = MIN( MAX( impaired_k_max(:) / kmax_pft(pft),                       &
                     TINY(1.0_real_jlslsm) ), 1.0_real_jlslsm )

SELECT CASE ( pft_conductance_model(pft) )

CASE ( CW_conductance )
  impaired_conductance_b(:) = conductance_b_pft(pft)                          &
                              * (1.0 - LOG(ratio(:)))                         &
                                **(1.0 / conductance_c_pft(pft))
  impaired_conductance_c(:) = conductance_c_pft(pft)                          &
                              * (impaired_conductance_b(:)                    &
                                 / conductance_b_pft(pft))                    &
                                **conductance_c_pft(pft)

CASE ( SOX_conductance )
  impaired_conductance_b(:) = conductance_b_pft(pft)                          &
                              * (2.0 / ratio(:) - 1.0)                        &
                                **(1.0 / conductance_c_pft(pft))
  impaired_conductance_c(:) = conductance_c_pft(pft) / ratio(:)               &
                              * (impaired_conductance_b(:)                    &
                                 / conductance_b_pft(pft))                    &
                                **conductance_c_pft(pft)

CASE DEFAULT
  errcode = 101
  CALL ereport(RoutineName, errcode,                                           &
    'pft_conductance_model should be CW_conductance (1) or SOX_conductance (2)')

END SELECT

END SUBROUTINE xylem_kmax_refit_curve

! ---------------------------------------------------------------------
! Lower the impaired kmax at the points with open stomata (see the module
! header). The conductances at the leaf and in the root zone are those of
! the impaired curve at psi_leaf and psi_root (Pa).
!
! MDK, 21/11/26
! ---------------------------------------------------------------------
SUBROUTINE update_xylem_impairment_kmax( n_land_pts, n_open_pts, open_index,   &
                                         pft, psi_leaf, psi_root,              &
                                         impaired_k_max )

USE pftparm, ONLY: ximpair_leaf_weight, ximpair_new_kmax_weight,              &
                   ximpair_threshold, kcrit
USE xylem_hydraulics_jls_mod, ONLY: leaf_conductance_jls

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: n_land_pts
INTEGER, INTENT(IN) :: n_open_pts
INTEGER, INTENT(IN) :: open_index(n_open_pts)
INTEGER, INTENT(IN) :: pft

! Water potentials of the leaf and the root zone (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_leaf(n_land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root(n_land_pts)

! The impaired kmax (mol m-2 s-1 Pa-1, leaf basis).
REAL(KIND=real_jlslsm), INTENT(IN OUT) :: impaired_k_max(n_land_pts)

REAL(KIND=real_jlslsm) :: b_imp(n_land_pts)
REAL(KIND=real_jlslsm) :: c_imp(n_land_pts)
REAL(KIND=real_jlslsm) :: kcrit_pts(n_land_pts)
REAL(KIND=real_jlslsm) :: k_leaf(n_land_pts)
REAL(KIND=real_jlslsm) :: k_root(n_land_pts)
REAL(KIND=real_jlslsm) :: k_new

INTEGER :: j, l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_XYLEM_IMPAIRMENT_KMAX'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( n_open_pts > 0 ) THEN

  ! The current impaired curve.
  CALL xylem_kmax_refit_curve( pft, n_land_pts, impaired_k_max, b_imp, c_imp )
  kcrit_pts(:) = kcrit(pft)

  CALL leaf_conductance_jls( pft, n_land_pts, psi_leaf, impaired_k_max,        &
                             kcrit_pts, b_imp, c_imp, k_leaf )
  CALL leaf_conductance_jls( pft, n_land_pts, psi_root, impaired_k_max,        &
                             kcrit_pts, b_imp, c_imp, k_root )

  DO j = 1, n_open_pts
    l = open_index(j)

    ! The target kmax, a weighted mean of the leaf and root zone conductance.
    k_new = ximpair_leaf_weight(pft) * k_leaf(l)                               &
            + (1.0 - ximpair_leaf_weight(pft)) * k_root(l)

    ! The conductance loss threshold at which damage starts.
    k_new = k_new / (impaired_k_max(l) - ximpair_threshold(pft))               &
            * impaired_k_max(l)

    ! The weighted mean of the new and the previous kmax.
    k_new = ximpair_new_kmax_weight(pft) * k_new                               &
            + (1.0 - ximpair_new_kmax_weight(pft)) * impaired_k_max(l)

    ! Damage is not recovered, and the plant keeps its critical conductance.
    k_new = MIN( impaired_k_max(l), k_new )
    impaired_k_max(l) = MAX( k_new, kcrit(pft) )

  END DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_xylem_impairment_kmax

END MODULE xylem_impairment_kmax_mod
