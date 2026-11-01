! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in SURFACE

MODULE desica_jls_mod

! *********************************************************************
! DESICA plant hydraulics (stomata_model = stomata_desica): the plant has
! a water store in its leaves and one in its stem, and the stomata close
! as the leaf dries.
!
! Stomata. The stomatal conductance gs follows Tuzet et al. (2003), with a
! factor fw (0 to 1) that falls as the leaf water potential psi_leaf
! falls:
!
!   fw = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi_leaf)))
!   gs = g1 fw An / ca
!
! where sf (MPa-1) sets how steeply the stomata close, psi_f (MPa) is the
! leaf water potential at which fw is about a half, g1 is a slope, An the
! net photosynthesis and ca the CO2 at the leaf. In the code this enters
! leaf_limits as ci = ca (1 - 1.6 / (g1 fw)).
!
! DESICA and CABLE take psi_leaf from the previous timestep. At JULES
! timesteps this swings from step to step: the leaf store fills or empties
! within minutes, so a step with open stomata drives psi_leaf down, the
! stomata shut in the next step, and so on. sf_stom instead finds the fw
! that agrees with the psi_leaf it produces by the end of the step,
!
!   fw = Tuzet factor at the psi_leaf reached with the transpiration E
!        that fw gives,
!
! by bisection in its big-leaf code. The right-hand side falls as fw
! rises, so there is one answer. This is what DESICA's repeated solves
! within a step tend to.
!
! Plant water. psi_leaf and the stem water potential psi_stem (Pa) are
! carried from step to step and follow Xu et al. (2016, Notes S1):
!
!   C_stem dpsi_stem/dt = Q - J                        (Eqn S4a)
!   C_leaf dpsi_leaf/dt = J - E                        (Eqn S4b)
!   J = k_leaf (psi_stem - psi_leaf - psi_h)           (Eqn S2a)
!   Q = k_root (psi_soil - psi_stem), and Q = 0 if
!       the soil is drier than the stem                (Eqn S1a)
!
! where C_stem and C_leaf are the stores' capacitances (mol m-2 Pa-1), Q
! the water taken up by the roots, J the sap flow from stem to leaf and E
! the transpiration (all mol m-2 s-1), k_leaf and k_root the conductances
! on the leaf and root side of the stem store (mol m-2 s-1 Pa-1), psi_h
! the fall in water potential from lifting water to the canopy height,
! and psi_soil the soil water potential (Pa).
!
! Each equation is solved exactly with the other store held fixed (Eqn
! S4c): psi_leaf first with psi_stem fixed, then J from the leaf store's
! water balance, then psi_stem with psi_soil fixed. Xu et al. use 600 s
! steps, so the JULES timestep is split into steps of 600 s or less, with
! E held fixed. psi_stem is kept above 2 P50 and psi_leaf above -20 MPa,
! as in CABLE-DESICA. Where a limit would be passed, J and E are cut so
! each store keeps its water balance. (DESICA and CABLE instead reset psi
! to the limit after the solve, which creates water.) The routine returns
! the transpiration the plant can deliver, and sf_stom lowers the canopy
! conductance gc to give it.
!
! Differences from Xu et al. (2016), as in DESICA: the stomata follow
! Tuzet, with no Cowan-Farquhar optimisation (their Eqns S3c and S6), and
! there is no drop in Vcmax and Jmax as the leaf dries (their Eqn S5).
!
! Shared with the stomatal optimisation (profit max): the plant's
! conductance is kmax_pft times LAI on the PFT's vulnerability curve
! (P50, P88), here at psi_stem as in Xu et al. Eqn S2b. psi_soil is
! psi_root_zone, the root-zone water potential weighted by conductance.
! The conductance from the soil to the root surface, ksr
! (l_som_rhizo_series), is in series before the root. The stem store sits
! part way along the plant, at som_leaf_resist_frac (0.5 in DESICA): the
! leaf-side conductance is k_xylem / som_leaf_resist_frac and the
! root-side one k_xylem / (1 - som_leaf_resist_frac), in series with ksr.
! The vulnerability curve acts on the xylem only. Conductances are held
! fixed over each short step.
!
! Water budget. sf_stom finds gs with an expected E and stores the step's
! inputs (desica_store_inputs). Once the surface fluxes are final, sf_evap
! moves psi_leaf and psi_stem on with the actual transpiration T
! (desica_commit) and takes the root uptake Q from the soil in place of
! T. The plant's stores take up the difference:
!
!   dW/dt = Q - T      (kg m-2 s-1)
!
! where W is the water in the plant. Where the plant cannot deliver T
! (a psi limit), the shortfall is added to Q. Where the soil cannot supply
! Q, Q is cut and the stem store gives the water (desica_cut_uptake). So
! the water in soil, plant and air adds up.
!
! psi_leaf and psi_stem are written to the dump ('psi_leaf_desica',
! 'psi_stem_desica', Pa). A value of 0 or above, such as an initial value
! of 0, means "not set", and the stores then start at psi_root_zone.
!
! References:
! Tuzet et al. (2003) Plant Cell Environ. 26: 1097-1116.
! Xu et al. (2016) New Phytol. 212: 80-95, doi:10.1111/nph.14009.
! De Kauwe et al. (2020) Biogeosciences 17: 3589-3612 (CABLE-DESICA).
!
! MDK, 24/10/26
! *********************************************************************

USE um_types, ONLY: real_jlslsm
USE water_constants_mod, ONLY: Water_kg_per_mol

IMPLICIT NONE

PRIVATE
PUBLIC :: desica_alloc, desica_fw, desica_hydraulics, tuzet_fw,                &
          desica_store_inputs, desica_commit, desica_cut_uptake,               &
          psi_leaf_desica, psi_stem_desica, flux_root_desica, flux_sap_desica, &
          dw_plant_desica

REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  ! Leaf water potential (Pa), (land_pts, npft).
  psi_leaf_desica(:,:),                                                        &
  ! Stem water potential (Pa).
  psi_stem_desica(:,:),                                                        &
  ! Root water uptake Q, timestep mean (kg m-2 ground s-1).
  flux_root_desica(:,:),                                                       &
  ! Sap flow stem -> leaf J, timestep mean (kg m-2 ground s-1).
  flux_sap_desica(:,:),                                                        &
  ! Change of plant water store, Q - T (kg m-2 ground s-1).
  dw_plant_desica(:,:),                                                        &
  ! This step's LAI, canopy height (m), root zone water potential (Pa) and
  ! soil-to-root conductance (mol m-2 s-1 Pa-1, < 0: none), stored by sf_stom
  ! for desica_commit.
  lai_desica(:,:), ht_desica(:,:), psi_soil_desica(:,:), ksr_desica(:,:)

! Longest sub-step (s), Xu et al. (2016).
REAL(KIND=real_jlslsm), PARAMETER :: dt_max = 600.0
! Lower bound on psi_leaf (Pa), as CABLE-DESICA.
REAL(KIND=real_jlslsm), PARAMETER :: psi_leaf_min = -20.0e6
! Below this LAI there is no leaf pool.
REAL(KIND=real_jlslsm), PARAMETER :: lai_min = 1.0e-3
! kg H2O per mol (JULES's value, so mol/kg round trips with the rest of the code
! close).
REAL(KIND=real_jlslsm), PARAMETER :: mol_h2o = Water_kg_per_mol

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DESICA_JLS_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! Allocate the module's arrays (land_pts, npft) unless they are already
! allocated. The water potentials start at 0, meaning "not set".
!-----------------------------------------------------------------------------
SUBROUTINE desica_alloc( land_pts )

USE jules_surface_types_mod, ONLY: npft

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: land_pts

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='DESICA_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( .NOT. ALLOCATED(psi_leaf_desica) ) THEN
  ALLOCATE( psi_leaf_desica(land_pts,npft), psi_stem_desica(land_pts,npft),    &
            flux_root_desica(land_pts,npft), flux_sap_desica(land_pts,npft),   &
            dw_plant_desica(land_pts,npft), lai_desica(land_pts,npft),         &
            ht_desica(land_pts,npft), psi_soil_desica(land_pts,npft),          &
            ksr_desica(land_pts,npft) )
  psi_leaf_desica(:,:)  = 0.0
  psi_stem_desica(:,:)  = 0.0
  flux_root_desica(:,:) = 0.0
  flux_sap_desica(:,:)  = 0.0
  dw_plant_desica(:,:)  = 0.0
  lai_desica(:,:)       = 0.0
  ht_desica(:,:)        = 0.0
  psi_soil_desica(:,:)  = 0.0
  ksr_desica(:,:)       = -1.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE desica_alloc

!-----------------------------------------------------------------------------
! Work out how much the stomata close because the leaf is dry, as the
! Tuzet factor fw (0 to 1) at the leaf water potential of the end of the
! previous timestep, for the stomatal conductance of this timestep (DESICA
! stomata). Where the leaf and stem water potentials are not yet set, they
! first start at the root-zone water potential (no wetter than -1 Pa).
!-----------------------------------------------------------------------------
SUBROUTINE desica_fw( ft, land_pts, veg_pts, veg_index, psi_root_zone, fw )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
! Root zone water potential (Pa).
REAL(KIND=real_jlslsm), INTENT(IN) :: psi_root_zone(land_pts)
! Tuzet factor (0-1).
REAL(KIND=real_jlslsm), INTENT(OUT) :: fw(land_pts)

INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='DESICA_FW'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL desica_alloc( land_pts )

fw(:) = 1.0
DO m = 1,veg_pts
  l = veg_index(m)
  ! Not set (no dump/initial value): start at the soil water potential.
  IF ( psi_stem_desica(l,ft) >= 0.0 .OR. psi_leaf_desica(l,ft) >= 0.0 ) THEN
    psi_leaf_desica(l,ft) = MIN(psi_root_zone(l), -1.0)
    psi_stem_desica(l,ft) = MIN(psi_root_zone(l), -1.0)
  END IF
  fw(l) = tuzet_fw(ft, psi_leaf_desica(l,ft))
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE desica_fw

!-----------------------------------------------------------------------------
! Return the Tuzet et al. (2003) factor fw (0 to 1), which scales the
! stomatal conductance down as the leaf dries, for a leaf water potential
! psi (Pa):
!
!   fw = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi)))
!
! with sf = sf_tuzet (MPa-1) and psi_f = psi_f_tuzet, both potentials in
! MPa.
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
! Move the leaf and stem water potentials on over the timestep for the
! canopy transpiration el (mol m-2 s-1), following Xu et al. (2016, Eqns
! S1, S2 and S4); see the module header for the equations.
!
! The timestep is split into steps of 600 s or less. In each, the xylem
! conductance is taken from the vulnerability curve at the stem's water
! potential, then the leaf store is solved with the stem held fixed, and
! the stem store with the soil held fixed. Where the leaf would pass
! -20 MPa or the stem 2 P50, the flows are cut to keep it at the limit.
! With almost no leaves (LAI < lai_min) the stem takes the soil's water
! potential and the leaf the stem's less the lift to the canopy top.
!
! It returns the leaf water potential at the end of the step, the
! plant's conductance from soil to leaf (kplant, mol m-2 s-1 Pa-1) and
! the transpiration the plant can deliver (el_hyd, mol m-2 s-1). With
! l_commit the new water potentials and the step's mean root uptake and
! sap flow are stored; otherwise the call only looks ahead and stores
! nothing.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE desica_hydraulics( ft, land_pts, veg_pts, veg_index, timestep,      &
                              lai, canht, psi_root_zone, ksr, el, l_commit,    &
                              psi_leaf, leaf_k, el_hyd )

USE pftparm, ONLY: kmax_pft, conductance_b_pft, conductance_c_pft, p50,        &
                   cap_leaf, cap_stem
USE jules_vegetation_mod, ONLY: som_leaf_resist_frac
USE planet_constants_mod, ONLY: g
USE water_constants_mod, ONLY: rho_water

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  ! Model timestep (s).
  timestep,                                                                    &
  ! Leaf area index.
  lai(land_pts),                                                               &
  ! Canopy height (m), for the gravity term.
  canht(land_pts),                                                             &
  ! Root zone (uptake-weighted) soil water potential (Pa).
  psi_root_zone(land_pts),                                                     &
  ! Soil-to-root conductance, per m2 ground (mol m-2 s-1 Pa-1); < 0: no soil
  ! link.
  ksr(land_pts),                                                               &
  ! Canopy transpiration (mol H2O m-2 s-1).
  el(land_pts)
! Store the new state (else only project it).
LOGICAL, INTENT(IN) :: l_commit
REAL(KIND=real_jlslsm), INTENT(OUT) ::                                         &
  ! Leaf water potential at the end of the step (Pa).
  psi_leaf(land_pts),                                                          &
  ! Kplant: soil-to-leaf conductance, ksr in series with the xylem at psi_stem
  ! (mol m-2 ground s-1 Pa-1).
  leaf_k(land_pts),                                                            &
  ! Transpiration the plant can deliver over the step (mol m-2 s-1): el, or less
  ! where psi_stem or psi_leaf reach their lower bounds.
  el_hyd(land_pts)

INTEGER :: l, m, n, n_sub
REAL(KIND=real_jlslsm) ::                                                      &
  ! Sub-step (s); canopy kmax, xylem, leaf-side, root-side and soil-to-leaf
  ! conductances (mol m-2 s-1 Pa-1); leaf and stem capacitances (mol m-2 Pa-1);
  ! gravity drop and stem bound (Pa); leaf and stem psi at the start and end of
  ! a sub-step (Pa); the pool solution's a, b and EXP(a dt); sap flow, its cap,
  ! root uptake and their step sums, and transpiration (mol m-2 s-1); and ksr
  ! here.
  dt, kmax_c, k_xylem, k_leaf, k_root, k_plant, c_leaf, c_stem, psi_h,         &
  psi_stem_min,                                                                &
  pl, ps, pl_new, ps_new, ap, bp, ex, j_sap, j_cap, q_root, q_sum, j_sum,      &
  e, e_sub, e_sum, ksr_l

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='DESICA_HYDRAULICS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n_sub = MAX(1, CEILING(timestep / dt_max))
dt    = timestep / REAL(n_sub)
psi_stem_min = 2.0 * p50(ft)        ! As CABLE-DESICA.

psi_leaf(:) = 0.0
leaf_k(:)   = 0.0
el_hyd(:)   = MAX(el(:), 0.0)

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
  ksr_l = ksr(l)
  c_leaf = cap_leaf(ft) * lai(l)
  c_stem = cap_stem(ft) * lai(l)
  q_sum  = 0.0
  j_sum  = 0.0
  e_sum  = 0.0

  DO n = 1,n_sub
    ! Xylem conductance on the PFT vulnerability curve at psi_stem
    ! (Eqn S2b), split either side of the stem store; the soil-to-root ksr
    ! is in series on the root side (soil -> root -> stem store).
    k_xylem = kmax_c * EXP(-(ABS(ps / conductance_b_pft(ft)))                  &
                           **conductance_c_pft(ft))
    k_xylem = MAX(k_xylem, 1.0e-6 * kmax_c)
    k_leaf  = k_xylem / som_leaf_resist_frac
    k_root  = k_xylem / (1.0 - som_leaf_resist_frac)
    k_plant = k_xylem
    IF ( ksr(l) >= 0.0 ) THEN
      k_root  = 1.0 / ( 1.0 / k_root + 1.0 / MAX(ksr_l, 1.0e-6 * kmax_c) )
      k_plant = 1.0 / ( 1.0 / k_xylem + 1.0 / MAX(ksr_l, 1.0e-6 * kmax_c) )
    END IF

    ! The lower bounds on psi_leaf and psi_stem are met by capping the
    ! fluxes, not the potentials, so that the stores stay in balance
    ! (C dpsi = in - out): clamping psi after the solve would let water
    ! leave a store that no longer falls. For a pool dx/dt = a x + b over dt,
    ! x ends at x_end for
    !
    !   b = a (x_end - x e^(a dt)) / (e^(a dt) - 1).

    ! Leaf (Eqn S4c), psi_stem fixed:
    !
    !   dpl/dt = ap pl + bp.
    !
    ! If psi_leaf would fall below its bound, the pathway cannot deliver E:
    ! cut the transpiration to what keeps psi_leaf at the bound.
    e_sub  = e
    ap     = -k_leaf / c_leaf
    ex     = EXP(ap * dt)
    bp     = (k_leaf * (ps - psi_h) - e_sub) / c_leaf
    pl_new = ((ap * pl + bp) * ex - bp) / ap
    IF ( pl_new < psi_leaf_min ) THEN
      bp     = ap * (psi_leaf_min - pl * ex) / (ex - 1.0)
      e_sub  = MAX(0.0, k_leaf * (ps - psi_h) - bp * c_leaf)
      bp     = (k_leaf * (ps - psi_h) - e_sub) / c_leaf
      pl_new = MAX(((ap * pl + bp) * ex - bp) / ap, psi_leaf_min)
    END IF
    ! Sap flow from mass conservation of the leaf pool.
    j_sap  = (pl_new - pl) * c_leaf / dt + e_sub

    ! Stem (Eqn S4a), psi_soil fixed. Uptake only while the soil is wetter
    ! than the stem, no leak to the soil (Eqn S1a).
    q_root = 0.0
    j_cap  = j_sap
    IF ( psi_root_zone(l) > ps ) THEN
      ap     = -k_root / c_stem
      ex     = EXP(ap * dt)
      bp     = (k_root * psi_root_zone(l) - j_sap) / c_stem
      ps_new = ((ap * ps + bp) * ex - bp) / ap
      IF ( ps_new < psi_stem_min ) THEN
        ! The store would empty past its bound: the sap flow is what keeps
        ! psi_stem at the bound.
        bp     = ap * (psi_stem_min - ps * ex) / (ex - 1.0)
        j_cap  = k_root * psi_root_zone(l) - bp * c_stem
        ps_new = psi_stem_min
      END IF
      q_root = (ps_new - ps) * c_stem / dt + MIN(j_sap, j_cap)
    END IF
    IF ( q_root <= 0.0 ) THEN
      q_root = 0.0
      j_cap  = MIN(j_sap, (ps - psi_stem_min) * c_stem / dt)
      ps_new = ps - MIN(j_sap, j_cap) * dt / c_stem
    END IF
    IF ( j_cap < j_sap ) THEN
      ! Less sap reaches the leaf: re-balance the leaf pool for that flow,
      ! and cut E further if psi_leaf would pass its bound.
      j_sap  = MAX(j_cap, 0.0)
      pl_new = pl + (j_sap - e_sub) * dt / c_leaf
      IF ( pl_new < psi_leaf_min ) THEN
        e_sub  = j_sap + (pl - psi_leaf_min) * c_leaf / dt
        pl_new = psi_leaf_min
      END IF
    END IF

    pl    = pl_new
    ps    = ps_new
    q_sum = q_sum + q_root
    j_sum = j_sum + j_sap
    e_sum = e_sum + e_sub
  END DO
  el_hyd(l) = e_sum / REAL(n_sub)

  IF ( l_commit ) THEN
    psi_leaf_desica(l,ft)  = pl
    psi_stem_desica(l,ft)  = ps
    flux_root_desica(l,ft) = q_sum / REAL(n_sub) * mol_h2o
    flux_sap_desica(l,ft)  = j_sum / REAL(n_sub) * mol_h2o
  END IF

  psi_leaf(l) = pl
  leaf_k(l)   = k_plant
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE desica_hydraulics

!-----------------------------------------------------------------------------
! Keep this step's leaf area, canopy height, root-zone water potential and
! soil-to-root conductance, so desica_commit can use them later in the
! step. Called from sf_stom.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE desica_store_inputs( ft, land_pts, veg_pts, veg_index, lai, canht,  &
                                psi_root_zone, ksr )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: lai(land_pts), canht(land_pts),          &
                                      psi_root_zone(land_pts), ksr(land_pts)
                            ! As desica_hydraulics.
INTEGER :: l, m

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='DESICA_STORE_INPUTS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL desica_alloc( land_pts )

DO m = 1,veg_pts
  l = veg_index(m)
  lai_desica(l,ft)      = lai(l)
  ht_desica(l,ft)       = canht(l)
  psi_soil_desica(l,ft) = psi_root_zone(l)
  ksr_desica(l,ft)      = ksr(l)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE desica_store_inputs

!-----------------------------------------------------------------------------
! Move the leaf and stem water potentials on with the actual
! transpiration of PFT ft, t_stom (kg m-2 s-1), and return q_soil, the
! water the roots take from the soil (kg m-2 s-1). Transpiration the plant
! could not deliver, because a water potential reached its limit, is
! added to q_soil, so it comes straight from the soil. Called from
! sf_evap once the surface fluxes are final.
!
! MDK, 24/10/26
!-----------------------------------------------------------------------------
SUBROUTINE desica_commit( ft, land_pts, npts, pts_index, timestep, t_stom,     &
                          q_soil )

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

INTEGER, INTENT(IN) :: ft, land_pts, npts, pts_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, t_stom(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: q_soil(land_pts)

INTEGER :: l, m
! Transpiration (mol m-2 s-1), and the desica_hydraulics outputs: psi_leaf (Pa),
! leaf_k and the deliverable transpiration.
REAL(KIND=real_jlslsm) :: el(land_pts), psi_l(land_pts), k_l(land_pts),        &
                          el_hyd(land_pts)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER :: RoutineName='DESICA_COMMIT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL desica_alloc( land_pts )

q_soil(:) = 0.0
el(:)     = 0.0
DO m = 1,npts
  l = pts_index(m)
  el(l) = MAX(t_stom(l), 0.0) / mol_h2o
END DO

CALL desica_hydraulics( ft, land_pts, npts, pts_index, timestep,               &
                        lai_desica(:,ft), ht_desica(:,ft),                     &
                        psi_soil_desica(:,ft), ksr_desica(:,ft), el, .TRUE.,   &
                        psi_l, k_l, el_hyd )

DO m = 1,npts
  l = pts_index(m)
  ! Transpiration the plant could not deliver (psi at a bound) comes
  ! straight from the soil, so no water is created.
  q_soil(l) = flux_root_desica(l,ft) + MAX(el(l) - el_hyd(l), 0.0) * mol_h2o
  flux_root_desica(l,ft) = q_soil(l)
  dw_plant_desica(l,ft)  = q_soil(l) - MAX(t_stom(l), 0.0)
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE desica_commit

!-----------------------------------------------------------------------------
! Take water from the stem store when the soil cannot supply all of the
! root uptake, so the plant still transpires what the stomata allowed.
! sf_evap calls it at point l of PFT ft. The uptake is cut by dq
! (kg m-2 s-1) and the stem store gives that water instead, so its water
! potential falls by
!
!   dq dt / (M_w C_stem)      (Pa)
!
! where dt is the timestep (s), M_w the molar mass of water (kg mol-1)
! and C_stem the stem store's capacitance (mol m-2 Pa-1).
!-----------------------------------------------------------------------------
SUBROUTINE desica_cut_uptake( ft, l, timestep, dq )

USE pftparm, ONLY: cap_stem

INTEGER, INTENT(IN) :: ft, l
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, dq

REAL(KIND=real_jlslsm) :: c_stem

c_stem = cap_stem(ft) * lai_desica(l,ft)
flux_root_desica(l,ft) = flux_root_desica(l,ft) - dq
dw_plant_desica(l,ft)  = dw_plant_desica(l,ft) - dq
IF ( c_stem > 0.0 ) THEN
  psi_stem_desica(l,ft) = psi_stem_desica(l,ft) - dq / mol_h2o * timestep      &
                                                  / c_stem
END IF

END SUBROUTINE desica_cut_uptake

END MODULE desica_jls_mod
