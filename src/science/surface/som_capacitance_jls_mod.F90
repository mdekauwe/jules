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
! *********************************************************************

USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

PRIVATE
PUBLIC :: som_cap_store, som_cap_psi_leaf, som_cap_commit

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

USE pftparm, ONLY: cap_stem
USE desica_jls_mod, ONLY: psi_stem_desica

INTEGER, INTENT(IN) :: ft, l
REAL(KIND=real_jlslsm), INTENT(IN) :: lai, psi_soil
REAL(KIND=real_jlslsm), INTENT(OUT) :: ps_n, c_store

ps_n = psi_stem_desica(l,ft)
IF ( ps_n >= 0.0 ) ps_n = MIN(psi_soil, 0.0)
c_store = cap_stem(ft) * lai

END SUBROUTINE store_state

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
REAL(KIND=real_jlslsm) :: ps_n, c_store, psi_soil, kmax_c, lo, hi, ps, f, e_sum

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

  ! f(ps) = C (ps - ps_n) / dt - Q(ps) + sum_c E_c(ps) rises with ps.
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
    f = c_store * (ps - ps_n) / timestep - q_root( ft, kmax_c, psi_soil, ps ) &
        + e_sum
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
REAL(KIND=real_jlslsm) :: t, ps_n, c_store, psi_soil, kmax_c, lo, hi, ps, f

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
    ! f(ps) = C (ps - ps_n) / dt - Q(ps) + T rises with ps; f(hi) >= 0.
    hi = MAX(psi_soil, ps_n)
    lo = MAX(ps_n - t * timestep / c_store, MIN(psi_floor(ft,1), psi_soil))
    lo = MIN(lo, hi)
    f  = c_store * (lo - ps_n) / timestep - q_root( ft, kmax_c, psi_soil, lo ) &
         + t
    IF ( f > 0.0 ) THEN
      ps = lo
    ELSE
      DO it = 1,n_bisect
        ps = 0.5 * (lo + hi)
        f  = c_store * (ps - ps_n) / timestep                                  &
             - q_root( ft, kmax_c, psi_soil, ps ) + t
        IF ( f > 0.0 ) THEN
          hi = ps
        ELSE
          lo = ps
        END IF
      END DO
      ps = 0.5 * (lo + hi)
    END IF
    q_soil(l) = MAX(t + c_store * (ps - ps_n) / timestep, 0.0)
  END IF

  q_soil(l) = q_soil(l) * mol_h2o
  psi_stem_desica(l,ft)  = ps
  flux_root_desica(l,ft) = q_soil(l)
  flux_sap_desica(l,ft)  = MAX(t_stom(l), 0.0)
  dw_plant_desica(l,ft)  = q_soil(l) - MAX(t_stom(l), 0.0)
END DO

END SUBROUTINE som_cap_commit

END MODULE som_capacitance_jls_mod
