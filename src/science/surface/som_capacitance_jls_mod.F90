! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE som_capacitance_jls_mod

! *********************************************************************
! Stem water store for the stomatal optimisation (profit max,
! l_som_plant_capacitance), the "psi_leaf target" form of
! NOTES_capacitance_sketch.md (after Martinetti et al. 2026, New Phytol.,
! doi:10.1111/nph.71143):
!
!   soil --(root side)--> [stem store, C] --(leaf side)--> leaf
!
! With l_som_plant_segments the root side is the root segment and the leaf
! side the stem and leaf segments in series, each on its own curve.
! Without, both sides are on the PFT's whole-plant curve, split by
! som_leaf_resist_frac (f) as in DESICA: kmax / (1 - f) below the store and
! kmax / f above it. Works for the big leaf (one class) and the two-leaf
! model (sunlit and shaded classes, sharing the root side).
!
! 1. The profit-max optimisation is unchanged: steady state, from the root
!    zone through the whole plant (l_som_plant_segments). It sets the leaf
!    water potential each leaf class accepts, psi_leaf*.
! 2. The store sets the transpiration that holds each class at psi_leaf*:
!    E_c = E_down(psi_s -> psi_leaf*), through the stem and leaf segments
!    of the class, capped at the class's gl_max. The store is advanced
!    implicitly, C (psi_s' - psi_s) / dt = Q(psi_s') - sum_c E_c(psi_s'),
!    with the root uptake Q through the shared root side (no hydraulic
!    redistribution, Q >= 0). The store is full in the morning, so E is
!    higher than the steady state; as it empties, E falls below it; it
!    refills at night. With C -> 0 this is the steady-state model (for
!    two leaf classes, with the root side shared by the classes).
! 3. Water accounting (as DESICA, sharing its state arrays): sf_stom stores
!    the step's inputs; sf_evap advances psi_s with the actual
!    transpiration T (som_cap_commit) and takes Q = T + C dpsi_s / dt from
!    the soil, so the plant store holds the difference.
!
! The state is psi_stem_desica (Pa, dumped); C is cap_stem * LAI (mol m-2
! ground Pa-1); conductances are mol m-2 ground s-1 Pa-1 and E is
! mol m-2 ground s-1, as in the profit max.
!
! Two optional refinements from the Puechabon Q. ilex measurements:
! - Storage resistance (cap_stem_k > 0; Salomon et al. 2017, R_S): the
!   store hangs off the xylem node between the root and leaf sides through
!   a conductance k_s * LAI. Over a step (implicit) it then acts as a source
!   at psi_s with conductance k_e = 1 / (1/k_s + dt/C) onto the node, and
!   releases F = k_e (psi_s - psi_x); with no resistance k_e = C/dt and the
!   node is the store.
! - Two-phase capacitance (cap_stem_dry_frac < 1; Salomon et al. 2020):
!   C falls to cap_stem_dry_frac of cap_stem below the store potential
!   cap_stem_psi_brk (smooth step, 0.1 MPa wide), from psi_s at the start of
!   the step.
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: som_cap_store, som_cap_psi_leaf, som_cap_commit,                    &
          som_cap_supply_set, som_cap_supply_off, som_cap_supply,             &
          som_cap_supply_on

! Form 2 (som_cap_form = 2): the optimiser's supply curve comes from the
! store. sf_stom sets these for the leaf class being optimised
! (som_cap_supply_set) and leaf_psi_jls then calls som_cap_supply.
LOGICAL, SAVE :: som_cap_supply_on = .FALSE.
INTEGER, SAVE :: sup_ft = 0
REAL(KIND=real_jlslsm), SAVE :: sup_dt = 0.0
REAL(KIND=real_jlslsm), ALLOCATABLE, SAVE ::                                   &
  sup_ke(:), sup_ps_n(:), sup_psi_soil(:)
                            ! Class share of the store's effective
                            ! conductance k_e (mol m-2 s-1 Pa-1), the store
                            ! psi at the start of the step and the root-zone
                            ! psi (Pa).

REAL(KIND=real_jlslsm), PARAMETER :: lai_min = 1.0e-3
                            ! Below this LAI there is no store.
REAL(KIND=real_jlslsm), PARAMETER :: mol_h2o = 0.018
                            ! kg H2O per mol.
REAL(KIND=real_jlslsm), PARAMETER :: x_lim = 10.0
                            ! Deepest point of a curve used, x = (psi/b)^c
                            ! (k/kmax = e^-10): below it a segment carries
                            ! no more flow, and the incomplete gamma series
                            ! (40 terms) stays accurate.
INTEGER, PARAMETER :: n_bisect = 50
                            ! Bisection iterations (psi to ~1e-15 of the
                            ! bracket).

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SOM_CAPACITANCE_JLS_MOD'

CONTAINS

!-----------------------------------------------------------------------------
! Weibull b (Pa) and c of segment iseg (1 root, 2 stem, 3 leaf): the
! segment's own curve with l_som_plant_segments, else the PFT's.
!-----------------------------------------------------------------------------
SUBROUTINE seg_curve( ft, iseg, bs, cs )

USE pftparm, ONLY: conductance_b, conductance_c, conductance_b_seg,          &
                   conductance_c_seg
USE jules_vegetation_mod, ONLY: l_som_plant_segments

INTEGER, INTENT(IN) :: ft, iseg
REAL(KIND=real_jlslsm), INTENT(OUT) :: bs, cs

IF ( l_som_plant_segments ) THEN
  bs = conductance_b_seg(ft,iseg)
  cs = conductance_c_seg(ft,iseg)
ELSE
  bs = conductance_b(ft)
  cs = conductance_c(ft)
END IF

END SUBROUTINE seg_curve

!-----------------------------------------------------------------------------
! Maximum conductance of segment iseg as a multiple of the whole-plant
! kmax. Without segments, 1 is the root side, kmax / (1 - f), and 3 the
! whole leaf side, kmax / f (2, the stem, is then not used).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION seg_fac( ft, iseg )

USE pftparm, ONLY: seg_kfac
USE jules_vegetation_mod, ONLY: l_som_plant_segments, som_leaf_resist_frac

INTEGER, INTENT(IN) :: ft, iseg

IF ( l_som_plant_segments ) THEN
  seg_fac = seg_kfac(ft,iseg)
ELSE IF ( iseg == 1 ) THEN
  seg_fac = 1.0 / (1.0 - som_leaf_resist_frac)
ELSE
  seg_fac = 1.0 / som_leaf_resist_frac
END IF

END FUNCTION seg_fac

!-----------------------------------------------------------------------------
! Lowest psi (Pa) of segment iseg's curve that is used, b x_lim^(1/c).
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION psi_floor( ft, iseg )

INTEGER, INTENT(IN) :: ft, iseg

REAL(KIND=real_jlslsm) :: bs, cs

CALL seg_curve( ft, iseg, bs, cs )
psi_floor = bs * x_lim**(1.0 / cs)

END FUNCTION psi_floor

!-----------------------------------------------------------------------------
! Flow (mol m-2 s-1) through segment iseg of PFT ft, maximum conductance
! kmx, from psi_in to psi_out: the integral of the cumulative Weibull,
! as leaf_psi_segments_jls, with both ends limited to psi_floor. Zero for
! psi_out >= psi_in.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION seg_flow( ft, iseg, kmx, psi_in, psi_out )

USE xylem_hydraulics_CW_jls_mod, ONLY: incomplete_gamma

INTEGER, INTENT(IN) :: ft, iseg
REAL(KIND=real_jlslsm), INTENT(IN) :: kmx, psi_in, psi_out

REAL(KIND=real_jlslsm) :: bs, cs, g(2)

seg_flow = 0.0
IF ( psi_out >= psi_in ) RETURN
CALL seg_curve( ft, iseg, bs, cs )
g  = incomplete_gamma(2, 1.0/cs,                                               &
                      [ (MIN(MAX(psi_out, psi_floor(ft,iseg)), 0.0)/bs)**cs,   &
                        (MIN(MAX(psi_in,  psi_floor(ft,iseg)), 0.0)/bs)**cs ])
seg_flow = MAX(kmx * (-bs/cs) * (g(1) - g(2)), 0.0)

END FUNCTION seg_flow

!-----------------------------------------------------------------------------
! Transpiration of a leaf class (whole-plant kmax kmax_c) from the store at
! psi_s to the leaf at psi_tgt, through its leaf side. With segments, the
! stem and leaf segments in series: bisection on the stem-leaf junction
! potential.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION e_down( ft, kmax_c, psi_s, psi_tgt )

USE jules_vegetation_mod, ONLY: l_som_plant_segments

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax_c, psi_s, psi_tgt

INTEGER :: it
REAL(KIND=real_jlslsm) :: lo, hi, mid, e1, e2

e_down = 0.0
IF ( psi_tgt >= psi_s .OR. kmax_c <= 0.0 ) RETURN
IF ( .NOT. l_som_plant_segments ) THEN
  e_down = seg_flow( ft, 3, kmax_c * seg_fac(ft,3), psi_s, psi_tgt )
  RETURN
END IF
! g(mid) = E_stem(psi_s -> mid) - E_leaf(mid -> psi_tgt) falls with mid:
! > 0 at mid = psi_tgt, < 0 at mid = psi_s.
lo = psi_tgt
hi = psi_s
DO it = 1,n_bisect
  mid = 0.5 * (lo + hi)
  e1  = seg_flow( ft, 2, kmax_c * seg_fac(ft,2), psi_s, mid )
  e2  = seg_flow( ft, 3, kmax_c * seg_fac(ft,3), mid, psi_tgt )
  IF ( e1 > e2 ) THEN
    lo = mid
  ELSE
    hi = mid
  END IF
END DO
e_down = 0.5 * (e1 + e2)

END FUNCTION e_down

!-----------------------------------------------------------------------------
! Root uptake (mol m-2 s-1) from the root zone at psi_soil to the store at
! psi_s, canopy conductance kmax_c: no hydraulic redistribution.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION q_root( ft, kmax_c, psi_soil, psi_s )

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax_c, psi_soil, psi_s

q_root = seg_flow( ft, 1, kmax_c * seg_fac(ft,1), psi_soil, psi_s )

END FUNCTION q_root

!-----------------------------------------------------------------------------
! Store state at a point: psi_s at the start of the step (initialised to
! psi_soil where not set, >= 0) and the capacitance.
!-----------------------------------------------------------------------------
SUBROUTINE store_state( ft, l, lai, psi_soil, ps_n, c_store )

USE pftparm, ONLY: cap_stem, cap_stem_dry_frac, cap_stem_psi_brk
USE desica_jls_mod, ONLY: psi_stem_desica

INTEGER, INTENT(IN) :: ft, l
REAL(KIND=real_jlslsm), INTENT(IN) :: lai, psi_soil
REAL(KIND=real_jlslsm), INTENT(OUT) :: ps_n, c_store

REAL(KIND=real_jlslsm), PARAMETER :: w_brk = 0.1e6
                            ! Width of the two-phase step (Pa).
REAL(KIND=real_jlslsm) :: x

ps_n = psi_stem_desica(l,ft)
IF ( ps_n >= 0.0 ) ps_n = MIN(psi_soil, 0.0)
c_store = cap_stem(ft) * lai
IF ( cap_stem_dry_frac(ft) < 1.0 ) THEN
  x = MAX(MIN((ps_n - cap_stem_psi_brk(ft)) / w_brk, 50.0), -50.0)
  c_store = c_store * ( cap_stem_dry_frac(ft)                                 &
                        + (1.0 - cap_stem_dry_frac(ft)) / (1.0 + EXP(-x)) )
END IF

END SUBROUTINE store_state

!-----------------------------------------------------------------------------
! Effective conductance (mol m-2 s-1 Pa-1) of the store onto the xylem node
! over a step dt: C/dt without storage resistance, else 1/(1/k_s + dt/C)
! with k_s = cap_stem_k * LAI.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION store_ke( ft, c_store, lai, dt )

USE pftparm, ONLY: cap_stem_k

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: c_store, lai, dt

store_ke = 0.0
IF ( c_store <= 0.0 ) RETURN
IF ( cap_stem_k(ft) > 0.0 .AND. lai > 0.0 ) THEN
  store_ke = 1.0 / ( 1.0 / (cap_stem_k(ft) * lai) + dt / c_store )
ELSE
  store_ke = c_store / dt
END IF

END FUNCTION store_ke

!-----------------------------------------------------------------------------
! sf_stom: the transpiration e_cls (mol m-2 ground s-1) of each of n_cls
! leaf classes that holds it at its accepted psi_leaf psi_tgt, capped at
! e_cap, with the store advanced implicitly (not committed). kmax_cls is the
! whole-plant kmax of each class (they share the root segment, whose
! conductance is the sum). Points with no store keep e_star.
!-----------------------------------------------------------------------------
SUBROUTINE som_cap_store( ft, land_pts, veg_pts, veg_index, n_cls, timestep, &
                          lai, psi_root_zone, kmax_cls, psi_tgt, e_star,     &
                          e_cap, e_cls, psi_s_new )

USE desica_jls_mod, ONLY: desica_alloc
USE pftparm, ONLY: pft_conductance_model
USE ereport_mod, ONLY: ereport

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts), n_cls
REAL(KIND=real_jlslsm), INTENT(IN) ::                                          &
  timestep, lai(land_pts), psi_root_zone(land_pts),                            &
  kmax_cls(land_pts,n_cls), psi_tgt(land_pts,n_cls), e_star(land_pts,n_cls),   &
  e_cap(land_pts,n_cls)
REAL(KIND=real_jlslsm), INTENT(OUT) :: e_cls(land_pts,n_cls),                  &
                                       psi_s_new(land_pts)

INTEGER :: l, m, ic, it, errcode
REAL(KIND=real_jlslsm) :: ps_n, c_store, psi_soil, kmax_c, lo, hi, ps, f,      &
                          e_sum, ke

CALL desica_alloc( land_pts )
IF ( pft_conductance_model(ft) /= 1 ) THEN
  errcode = 101
  CALL ereport('som_cap_store', errcode, 'l_som_plant_capacitance is ' //    &
               'coded for the cumulative Weibull (pft_conductance_model = 1)')
END IF
e_cls(:,:)   = e_star(:,:)
psi_s_new(:) = 0.0

DO m = 1,veg_pts
  l = veg_index(m)
  psi_soil = MIN(psi_root_zone(l), 0.0)
  CALL store_state( ft, l, lai(l), psi_soil, ps_n, c_store )
  psi_s_new(l) = ps_n
  IF ( lai(l) < lai_min .OR. c_store <= 0.0 ) CYCLE
  kmax_c = SUM(kmax_cls(l,:))
  ke = store_ke( ft, c_store, lai(l), timestep )

  ! ps is the xylem node (the store itself without storage resistance).
  ! f(ps) = k_e (ps - ps_n) - Q(ps) + sum_c E_c(ps) rises with ps.
  ! hi: Q = 0 and the store term >= 0, so f >= 0. lo: below every
  ! psi_tgt E_c = 0 and the store term <= 0, so f <= 0.
  hi = MAX(psi_soil, ps_n)
  lo = MIN(ps_n, psi_soil, MINVAL(psi_tgt(l,:)))
  DO it = 1,n_bisect
    ps = 0.5 * (lo + hi)
    e_sum = 0.0
    DO ic = 1,n_cls
      e_sum = e_sum + MIN(e_down( ft, kmax_cls(l,ic), ps, psi_tgt(l,ic) ),    &
                          e_cap(l,ic))
    END DO
    f = ke * (ps - ps_n) - q_root( ft, kmax_c, psi_soil, ps ) + e_sum
    IF ( f > 0.0 ) THEN
      hi = ps
    ELSE
      lo = ps
    END IF
  END DO
  ps = 0.5 * (lo + hi)
  psi_s_new(l) = ps
  DO ic = 1,n_cls
    e_cls(l,ic) = MIN(e_down( ft, kmax_cls(l,ic), ps, psi_tgt(l,ic) ),       &
                      e_cap(l,ic))
  END DO
END DO

END SUBROUTINE som_cap_store

!-----------------------------------------------------------------------------
! Leaf water potential of a class (whole-plant kmax kmax_c) transpiring e
! from the store at psi_s, through its leaf side: each segment's outlet by
! bisection, down to the segment's psi_floor.
!-----------------------------------------------------------------------------
REAL(KIND=real_jlslsm) FUNCTION som_cap_psi_leaf( ft, kmax_c, psi_s, e )

USE jules_vegetation_mod, ONLY: l_som_plant_segments

INTEGER, INTENT(IN) :: ft
REAL(KIND=real_jlslsm), INTENT(IN) :: kmax_c, psi_s, e

INTEGER :: iseg, iseg1, it
REAL(KIND=real_jlslsm) :: psi_in, lo, hi, mid

iseg1 = 2
IF ( .NOT. l_som_plant_segments ) iseg1 = 3
psi_in = psi_s
DO iseg = iseg1,3
  ! No flow: the leaf is at the store's potential.
  IF ( e <= 0.0 .OR. kmax_c <= 0.0 ) EXIT
  lo = MIN(psi_floor(ft,iseg), psi_in)
  hi = psi_in
  DO it = 1,n_bisect
    mid = 0.5 * (lo + hi)
    IF ( seg_flow( ft, iseg, kmax_c * seg_fac(ft,iseg), psi_in, mid ) > e )   &
    THEN
      lo = mid
    ELSE
      hi = mid
    END IF
  END DO
  psi_in = 0.5 * (lo + hi)
END DO
som_cap_psi_leaf = psi_in

END FUNCTION som_cap_psi_leaf

!-----------------------------------------------------------------------------
! sf_evap: advance the store with the actual transpiration t_stom of PFT ft
! (kg m-2 s-1) and return the root uptake to take from the soil, q_soil
! (kg m-2 s-1), Q = T + C dpsi_s / dt (so water is conserved). If the root
! cannot supply T with the store at its floor (psi_floor of the root
! side), the shortfall comes straight from the soil, as in
! DESICA. Uses the inputs desica_store_inputs saved in sf_stom.
!-----------------------------------------------------------------------------
SUBROUTINE som_cap_commit( ft, land_pts, npts, pts_index, timestep, t_stom,   &
                           q_soil )

USE pftparm, ONLY: kmax_pft
USE desica_jls_mod, ONLY: desica_alloc, psi_stem_desica, flux_root_desica,  &
                          flux_sap_desica, dw_plant_desica, lai_desica,      &
                          psi_soil_desica

INTEGER, INTENT(IN) :: ft, land_pts, npts, pts_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, t_stom(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: q_soil(land_pts)

INTEGER :: l, m, it
REAL(KIND=real_jlslsm) :: t, ps_n, c_store, psi_soil, kmax_c, lo, hi, ps, f,   &
                          ke

CALL desica_alloc( land_pts )
q_soil(:) = 0.0

DO m = 1,npts
  l = pts_index(m)
  t = MAX(t_stom(l), 0.0) / mol_h2o
  psi_soil = MIN(psi_soil_desica(l,ft), 0.0)
  CALL store_state( ft, l, lai_desica(l,ft), psi_soil, ps_n, c_store )

  IF ( lai_desica(l,ft) < lai_min .OR. c_store <= 0.0 ) THEN
    ps = psi_soil
    q_soil(l) = t
  ELSE
    kmax_c = kmax_pft(ft) * lai_desica(l,ft)
    ke = store_ke( ft, c_store, lai_desica(l,ft), timestep )
    ! ps is the xylem node; the store releases F = k_e (ps_n - ps).
    ! f(ps) = k_e (ps - ps_n) - Q(ps) + T rises with ps; f(hi) >= 0.
    hi = MAX(psi_soil, ps_n)
    lo = MAX(ps_n - t / ke, MIN(psi_floor(ft,1), psi_soil))
    lo = MIN(lo, hi)
    f  = ke * (lo - ps_n) - q_root( ft, kmax_c, psi_soil, lo ) + t
    IF ( f > 0.0 ) THEN
      ps = lo
    ELSE
      DO it = 1,n_bisect
        ps = 0.5 * (lo + hi)
        f  = ke * (ps - ps_n) - q_root( ft, kmax_c, psi_soil, ps ) + t
        IF ( f > 0.0 ) THEN
          hi = ps
        ELSE
          lo = ps
        END IF
      END DO
      ps = 0.5 * (lo + hi)
    END IF
    q_soil(l) = MAX(t + ke * (ps - ps_n), 0.0)
    ! The store loses F dt: psi_s' = ps_n - F dt / C (= ps without storage
    ! resistance).
    ps = ps_n + ke * (ps - ps_n) * timestep / c_store
  END IF

  q_soil(l) = q_soil(l) * mol_h2o
  psi_stem_desica(l,ft)  = ps
  flux_root_desica(l,ft) = q_soil(l)
  flux_sap_desica(l,ft)  = MAX(t_stom(l), 0.0)
  dw_plant_desica(l,ft)  = q_soil(l) - MAX(t_stom(l), 0.0)
END DO

END SUBROUTINE som_cap_commit


!-----------------------------------------------------------------------------
! Form 2: set the store seen by the optimiser for one leaf class (share of
! the store's capacitance; the class's own kmax carries its share of the
! root side), or switch it off.
!-----------------------------------------------------------------------------
SUBROUTINE som_cap_supply_set( ft, land_pts, veg_pts, veg_index, timestep,    &
                               lai, psi_root_zone, share )

USE desica_jls_mod, ONLY: desica_alloc

INTEGER, INTENT(IN) :: ft, land_pts, veg_pts, veg_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: timestep, lai(land_pts),                &
                                      psi_root_zone(land_pts), share(land_pts)

INTEGER :: l, m
REAL(KIND=real_jlslsm) :: ps_n, c_store

CALL desica_alloc( land_pts )
IF ( .NOT. ALLOCATED(sup_ke) ) THEN
  ALLOCATE( sup_ke(land_pts), sup_ps_n(land_pts), sup_psi_soil(land_pts) )
END IF
sup_ke(:)       = 0.0
sup_ps_n(:)     = 0.0
sup_psi_soil(:) = 0.0
DO m = 1,veg_pts
  l = veg_index(m)
  sup_psi_soil(l) = MIN(psi_root_zone(l), 0.0)
  CALL store_state( ft, l, lai(l), sup_psi_soil(l), ps_n, c_store )
  sup_ps_n(l) = ps_n
  IF ( lai(l) >= lai_min )                                                    &
    sup_ke(l) = store_ke( ft, c_store, lai(l), timestep ) * share(l)
END DO
sup_ft = ft
sup_dt = timestep
som_cap_supply_on = .TRUE.

END SUBROUTINE som_cap_supply_set

SUBROUTINE som_cap_supply_off()
som_cap_supply_on = .FALSE.
END SUBROUTINE som_cap_supply_off

!-----------------------------------------------------------------------------
! Form 2 supply curve (called by leaf_psi_jls in place of the steady-state
! solvers): for each sampled transpiration E of a class with whole-plant
! kmax, the xylem node psi_x at the store from
!   k_e (psi_x - psi_s) = Q_root(psi_x) - E               (implicit, Newton)
! (k_e = C/dt without storage resistance, when psi_x is the store's
! end-of-step psi), then the leaf through the leaf side at E (Newton per
! segment), and the conductance the hydraulic cost uses, k = -dE/dpsi_leaf,
! carried down as in leaf_psi_segments_jls from
! dpsi_x/dE = -1 / (k_e + k_root(psi_x)). C -> 0 gives the steady-state
! series model (dpsi_x/dE = -1 / k_root).
!-----------------------------------------------------------------------------
SUBROUTINE som_cap_supply( pft, n_e_leaf, land_pts, open_pnts, veg_index,     &
                           open_index, e_leaf, kmax, leaf_psi, leaf_k )

USE jules_vegetation_mod, ONLY: l_som_plant_segments

INTEGER, INTENT(IN) :: pft, n_e_leaf, land_pts, open_pnts
INTEGER, INTENT(IN) :: veg_index(land_pts), open_index(land_pts)
REAL(KIND=real_jlslsm), INTENT(IN) :: e_leaf(n_e_leaf, open_pnts),            &
                                      kmax(land_pts)
REAL(KIND=real_jlslsm), INTENT(OUT) :: leaf_psi(n_e_leaf, open_pnts),         &
                                       leaf_k(n_e_leaf, open_pnts)

INTEGER, PARAMETER :: max_it = 40
REAL(KIND=real_jlslsm), PARAMETER :: k_floor = 1.0e-12, tol = 1.0e-7

INTEGER :: i, j, l, it, iseg, iseg1, ft
REAL(KIND=real_jlslsm) :: e, ps, hi, lo, g, kr, kmx_r, c_dt, d, psi_in,      &
                          psi_out, kmx, k_in, k_out, e_cur, bs, cs, pfl,      &
                          a, b, step, g_lo

ft = sup_ft
iseg1 = 2
IF ( .NOT. l_som_plant_segments ) iseg1 = 3

DO j = 1,open_pnts
  l = veg_index(open_index(j))
  c_dt  = sup_ke(l)
  kmx_r = kmax(l) * seg_fac(ft,1)
  CALL seg_curve( ft, 1, bs, cs )
  hi = MAX(sup_psi_soil(l), sup_ps_n(l))
  pfl = MIN(psi_floor(ft,1), sup_psi_soil(l))
  ps = sup_ps_n(l)
  IF ( c_dt <= 0.0 ) ps = sup_psi_soil(l)
  DO i = 1,n_e_leaf
    e = MAX(e_leaf(i,j), 0.0)
    ! Node: g(ps) = k_e (ps - ps_n) - Q(ps) + E rises with ps. Safeguarded
    ! Newton on the bracket [pfl, hi] (g(hi) >= 0), starting from the
    ! previous sample's ps; bisection whenever a step leaves the bracket.
    ! If g(pfl) > 0 the root and store cannot supply E: ps = pfl.
    a = pfl
    b = hi
    g_lo = c_dt * (a - sup_ps_n(l)) - q_root( ft, kmax(l), sup_psi_soil(l), a ) &
           + e
    IF ( g_lo > 0.0 ) THEN
      ps = pfl
    ELSE
      ps = MIN(MAX(ps, a), b)
      DO it = 1,max_it
        IF ( ps < sup_psi_soil(l) ) THEN
          kr = kmx_r * EXP( -(MAX(ps, pfl)/bs)**cs )
        ELSE
          kr = 0.0
        END IF
        g = c_dt * (ps - sup_ps_n(l)) - q_root( ft, kmax(l), sup_psi_soil(l), ps ) &
            + e
        IF ( g > 0.0 ) THEN
          b = ps
        ELSE
          a = ps
        END IF
        IF ( ABS(g) <= tol * MAX(e, kmx_r * 1.0e3) .OR. b - a < 1.0 ) EXIT
        IF ( c_dt + kr > 0.0 ) THEN
          step = ps - g / (c_dt + kr)
        ELSE
          step = a - 1.0
        END IF
        IF ( step <= a .OR. step >= b ) step = 0.5 * (a + b)
        ps = step
      END DO
    END IF
    IF ( ps < sup_psi_soil(l) ) THEN
      kr = kmx_r * EXP( -(MAX(ps, pfl)/bs)**cs )
    ELSE
      kr = 0.0
    END IF
    d = -1.0 / MAX(c_dt + kr, k_floor * kmx_r)

    ! Leaf side at E, from the store.
    psi_in = ps
    DO iseg = iseg1,3
      kmx = kmax(l) * seg_fac(ft,iseg)
      CALL seg_curve( ft, iseg, bs, cs )
      k_in = MAX(kmx * EXP( -(MAX(psi_in, psi_floor(ft,iseg))/bs)**cs ),      &
                 k_floor * kmx)
      ! Outlet: h(po) = seg_flow(psi_in, po) - E falls with po; bracket
      ! [floor, psi_in] (h(psi_in) = -E <= 0), safeguarded Newton. If the
      ! segment cannot carry E even to its floor, the outlet is the floor.
      a = MIN(psi_floor(ft,iseg), psi_in)
      b = psi_in
      IF ( e <= 0.0 ) THEN
        psi_out = psi_in
      ELSE IF ( seg_flow( ft, iseg, kmx, psi_in, a ) <= e ) THEN
        psi_out = a
      ELSE
        psi_out = MIN(MAX(psi_in - e / k_in, a), b)
        DO it = 1,max_it
          e_cur = seg_flow( ft, iseg, kmx, psi_in, psi_out )
          IF ( e_cur > e ) THEN
            a = psi_out
          ELSE
            b = psi_out
          END IF
          IF ( ABS(e - e_cur) <= tol * e .OR. b - a < 1.0 ) EXIT
          k_out = MAX(kmx * EXP( -(MAX(psi_out, psi_floor(ft,iseg))/bs)**cs ),&
                      k_floor * kmx)
          step = psi_out - (e - e_cur) / k_out
          IF ( step <= a .OR. step >= b ) step = 0.5 * (a + b)
          psi_out = step
        END DO
      END IF
      k_out = MAX(kmx * EXP( -(MAX(psi_out, psi_floor(ft,iseg))/bs)**cs ),    &
                  k_floor * kmx)
      d = (k_in * d - 1.0) / k_out
      psi_in = psi_out
    END DO
    CALL seg_curve( ft, 1, bs, cs )

    leaf_psi(i,j) = psi_in
    IF ( d < 0.0 .AND. ABS(d) < HUGE(1.0_real_jlslsm) ) THEN
      leaf_k(i,j) = -1.0 / d
    ELSE
      leaf_k(i,j) = 0.0
    END IF
  END DO
END DO

END SUBROUTINE som_cap_supply

END MODULE som_capacitance_jls_mod
