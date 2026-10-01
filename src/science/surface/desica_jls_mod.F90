! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE desica_jls_mod

! *********************************************************************
! DESICA plant hydraulics (stomata_model = stomata_desica).
!
! Stomata follow the Tuzet et al. (2003) closure on leaf water potential,
!   fw = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi_leaf))),
!   gs = g1 fw An / ca   (set in leaf_limits as ci = ca (1 - 1.6/(g1 fw))).
! DESICA and CABLE take psi_leaf from the previous timestep. At 30 min that
! oscillates (the leaf pool equilibrates in ~minutes, so a step at high gs
! drives psi_leaf down, the next step shuts, and so on), so sf_stom instead
! solves fw = Tuzet(psi_leaf at the end of the step, for the E that fw
! gives) by bisection, in the big-leaf block of sf_stom: the right-hand
! side falls as fw rises, so the root is unique. This is the limit of
! DESICA's (Python) repeated solves within a step.
!
! psi_leaf and psi_stem are prognostic and follow Xu et al. (2016) Notes S1:
!   C_stem dpsi_stem/dt = Q - J                     (Eqn S4a)
!   C_leaf dpsi_leaf/dt = J - E                     (Eqn S4b)
!   J = k_stem2leaf (psi_stem - psi_leaf - psi_h)   (Eqn S2a)
!   Q = k_root2stem (psi_soil - psi_stem), Q = 0 if psi_soil < psi_stem (S1a)
! Each equation is solved analytically with the other pool held fixed
! (Eqn S4c): psi_leaf first with psi_stem fixed, J from mass conservation,
! then psi_stem with psi_soil fixed. Xu et al. use 600 s steps, so the
! JULES timestep is sub-stepped to <= 600 s, with E held fixed.
!
! Differences from Xu et al. (as in DESICA): the gs model is Tuzet, not
! the Cowan-Farquhar optimisation (Eqn S3c, S6), and there is no
! non-stomatal Vcmax/Jmax down-regulation (Eqn S5). Shared with the JULES
! stomatal optimisation (profit max): the plant conductance is
! kmax_pft * LAI, on the PFT's P50/P88 vulnerability curve evaluated at
! psi_stem (Eqn S2b uses psi_stem too), and the soil water potential is
! psi_root_zone, the uptake-weighted root-zone value; there is no soil-root
! resistance (as in the profit max). DESICA's placement of the stem store
! halfway along the plant is som_leaf_resist_frac (0.5 = DESICA): the
! leaf-side conductance is k_plant / som_leaf_resist_frac and the root-side
! one k_plant / (1 - som_leaf_resist_frac).
!
! Soil water extraction stays E (as in CABLE-DESICA); Q is a diagnostic.
! The state is not written to the dump: on a (re)start psi_leaf and
! psi_stem start at psi_root_zone.
!
! References:
! Tuzet et al. (2003) Plant Cell Environ. 26: 1097-1116.
! Xu et al. (2016) New Phytol. 212: 80-95, doi:10.1111/nph.14009.
! De Kauwe et al. (2020) Biogeosciences 17: 3589-3612 (CABLE-DESICA).
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: desica_fw, desica_hydraulics, tuzet_fw,                              &
          psi_leaf_desica, psi_stem_desica, flux_root_desica, flux_sap_desica

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  psi_leaf_desica(:,:),                                                        &
                            ! Leaf water potential (Pa), (land_pts, npft).
  psi_stem_desica(:,:),                                                        &
                            ! Stem water potential (Pa).
  flux_root_desica(:,:),                                                       &
                            ! Root water uptake Q, timestep mean
                            ! (kg m-2 ground s-1).
  flux_sap_desica(:,:)
                            ! Sap flow stem -> leaf J, timestep mean
                            ! (kg m-2 ground s-1).

LOGICAL, ALLOCATABLE, SAVE :: l_desica_init(:,:)
                            ! State has been initialised at this point/PFT.

REAL(KIND=real_jlslsm), PARAMETER :: dt_max = 600.0
                            ! Longest sub-step (s), Xu et al. (2016).
REAL(KIND=real_jlslsm), PARAMETER :: psi_leaf_min = -20.0e6
                            ! Lower bound on psi_leaf (Pa), as CABLE-DESICA.
REAL(KIND=real_jlslsm), PARAMETER :: lai_min = 1.0e-3
                            ! Below this LAI there is no leaf pool.
REAL(KIND=real_jlslsm), PARAMETER :: mol_h2o = 0.018
                            ! kg H2O per mol.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DESICA_JLS_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! Tuzet factor from the psi_leaf of the previous timestep. Initialises the
! state to psi_root_zone at points not seen before.
!-----------------------------------------------------------------------------
SUBROUTINE desica_fw( ft, land_pts, veg_pts, veg_index, psi_root_zone, fw )

USE jules_surface_types_mod, ONLY: npft
INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root_zone(land_pts)
                            ! Root zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(OUT) :: fw(land_pts)
                            ! Tuzet factor (0-1).

INTEGER :: l, m

IF ( .NOT. ALLOCATED(psi_leaf_desica) ) THEN
  ALLOCATE( psi_leaf_desica(land_pts,npft), psi_stem_desica(land_pts,npft),   &
            flux_root_desica(land_pts,npft), flux_sap_desica(land_pts,npft),  &
            l_desica_init(land_pts,npft) )
  psi_leaf_desica(:,:)  = 0.0
  psi_stem_desica(:,:)  = 0.0
  flux_root_desica(:,:) = 0.0
  flux_sap_desica(:,:)  = 0.0
  l_desica_init(:,:)    = .FALSE.
END IF

fw(:) = 1.0
DO m = 1,veg_pts
  l = veg_index(m)
  IF ( .NOT. l_desica_init(l,ft) ) THEN
    psi_leaf_desica(l,ft) = MIN(psi_root_zone(l), 0.0)
    psi_stem_desica(l,ft) = MIN(psi_root_zone(l), 0.0)
    l_desica_init(l,ft)   = .TRUE.
  END IF
  fw(l) = tuzet_fw(ft, psi_leaf_desica(l,ft))
END DO

END SUBROUTINE desica_fw

!-----------------------------------------------------------------------------
! Tuzet et al. (2003) factor for a leaf water potential psi (Pa).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION tuzet_fw( ft, psi ) RESULT( fw )

USE pftparm, ONLY: sf_tuzet, psi_f_tuzet

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: psi

REAL(KIND=real_jlslsm) :: psi_f_mpa

psi_f_mpa = psi_f_tuzet(ft) * 1.0e-6
fw = ( 1.0 + EXP(sf_tuzet(ft) * psi_f_mpa) )                                   &
     / ( 1.0 + EXP(sf_tuzet(ft) * (psi_f_mpa - psi * 1.0e-6)) )
fw = MAX(0.0, MIN(1.0, fw))

END FUNCTION tuzet_fw

!-----------------------------------------------------------------------------
! Advance psi_leaf and psi_stem over the timestep for the canopy
! transpiration el (Xu et al. 2016, Eqns S1, S2, S4).
!-----------------------------------------------------------------------------
SUBROUTINE desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,     &
                              lai, canht, psi_root_zone, el, l_commit,         &
                              psi_leaf, leaf_k )

USE pftparm, ONLY: kmax_pft, conductance_b, conductance_c, p50,               &
                   cap_leaf, cap_stem
USE jules_vegetation_mod, ONLY: som_leaf_resist_frac
USE planet_constants_mod, ONLY: g
USE water_constants_mod, ONLY: rho_water

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep,                                                                    &
                            ! Model timestep (s).
  lai(land_pts),                                                               &
                            ! Leaf area index.
  canht(land_pts),                                                             &
                            ! Canopy height (m), for the gravity term.
  psi_root_zone(land_pts),                                                     &
                            ! Root zone (uptake-weighted) soil water
                            ! potential (Pa).
  el(land_pts)
                            ! Canopy transpiration (mol H2O m-2 s-1).
LOGICAL, INTENT(IN) :: l_commit
                            ! Store the new state (else only project it).
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  psi_leaf(land_pts),                                                          &
                            ! Leaf water potential at the end of the step (Pa).
  leaf_k(land_pts)
                            ! Plant conductance at psi_stem
                            ! (mol m-2 ground s-1 Pa-1).

INTEGER :: l, m, n, n_sub
REAL(KIND=real_jlslsm) ::                                                      &
  dt, kmax_c, k_plant, k_leaf, k_root, c_leaf, c_stem, psi_h, psi_stem_min,    &
  pl, ps, pl_new, ps_new, ap, bp, j_sap, q_root, q_sum, j_sum, e

n_sub = MAX(1, CEILING(timestep / dt_max))
dt    = timestep / REAL(n_sub)
psi_stem_min = 2.0 * p50(ft)        ! As CABLE-DESICA.

psi_leaf(:) = 0.0
leaf_k(:)   = 0.0

DO m = 1,veg_pts
  l = veg_index(m)

  pl    = psi_leaf_desica(l,ft)
  ps    = psi_stem_desica(l,ft)
  psi_h = rho_water * g * MAX(canht(l), 0.0)

  IF ( lai(l) < lai_min ) THEN
    ! No leaves (and so no store, which scales with LAI): the stem is at
    ! equilibrium with the soil and psi_leaf follows the stem minus the
    ! gravity drop (Xu et al. 2016).
    ps = MIN(psi_root_zone(l), 0.0)
    psi_leaf(l) = ps - psi_h
    IF ( l_commit ) THEN
      psi_stem_desica(l,ft)  = ps
      psi_leaf_desica(l,ft)  = psi_leaf(l)
      flux_root_desica(l,ft) = 0.0
      flux_sap_desica(l,ft)  = 0.0
    END IF
    CYCLE
  END IF

  e      = MAX(el(l), 0.0)
  kmax_c = kmax_pft(ft) * lai(l)
  c_leaf = cap_leaf(ft) * lai(l)
  c_stem = cap_stem(ft) * lai(l)
  q_sum  = 0.0
  j_sum  = 0.0

  DO n = 1,n_sub
    ! Plant conductance on the PFT vulnerability curve at psi_stem
    ! (Eqn S2b), split either side of the stem store.
    k_plant = kmax_c * EXP(-(ABS(ps / conductance_b(ft)))**conductance_c(ft))
    k_plant = MAX(k_plant, 1.0e-6 * kmax_c)
    k_leaf  = k_plant / som_leaf_resist_frac
    k_root  = k_plant / (1.0 - som_leaf_resist_frac)

    ! Leaf (Eqn S4c), psi_stem fixed: dpl/dt = ap pl + bp.
    ap     = -k_leaf / c_leaf
    bp     = (k_leaf * (ps - psi_h) - e) / c_leaf
    pl_new = ((ap * pl + bp) * EXP(ap * dt) - bp) / ap
    pl_new = MAX(pl_new, psi_leaf_min)
    ! Sap flow from mass conservation of the leaf pool.
    j_sap  = (pl_new - pl) * c_leaf / dt + e

    ! Stem (Eqn S4a), psi_soil fixed. Uptake only while the soil is wetter
    ! than the stem, no leak to the soil (Eqn S1a).
    q_root = 0.0
    IF ( psi_root_zone(l) > ps ) THEN
      ap     = -k_root / c_stem
      bp     = (k_root * psi_root_zone(l) - j_sap) / c_stem
      ps_new = ((ap * ps + bp) * EXP(ap * dt) - bp) / ap
      q_root = (ps_new - ps) * c_stem / dt + j_sap
    END IF
    IF ( q_root <= 0.0 ) THEN
      q_root = 0.0
      ps_new = ps - j_sap * dt / c_stem
    END IF
    ps_new = MAX(ps_new, psi_stem_min)

    pl    = pl_new
    ps    = ps_new
    q_sum = q_sum + q_root
    j_sum = j_sum + j_sap
  END DO

  IF ( l_commit ) THEN
    psi_leaf_desica(l,ft)  = pl
    psi_stem_desica(l,ft)  = ps
    flux_root_desica(l,ft) = q_sum / REAL(n_sub) * mol_h2o
    flux_sap_desica(l,ft)  = j_sum / REAL(n_sub) * mol_h2o
  END IF

  psi_leaf(l) = pl
  leaf_k(l)   = k_plant
END DO

END SUBROUTINE desica_hydraulics

END MODULE desica_jls_mod
