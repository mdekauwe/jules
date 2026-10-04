! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! This module contains variables used for reading in pftparm data
! and initialisations

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.


MODULE pftparm_io

USE max_dimensions, ONLY:                                                      &
  npft_max
USE missing_data_mod, ONLY: imdi, rmdi
USE um_types, ONLY: real_jlslsm

IMPLICIT NONE

!---------------------------------------------------------------------
! Set up variables to use in IO (a fixed size version of each array
! in pftparm that we want to initialise).
!---------------------------------------------------------------------
#if !defined(UM_JULES)
INTEGER ::                                                                     &
  fsmc_mod_io(npft_max) = imdi

LOGICAL ::                                                                     &
  calc_rz_psi_io(npft_max) = .FALSE.

REAL(KIND=real_jlslsm) ::                                                      &
  canht_ft_io(npft_max) = rmdi,                                                &
  lai_io(npft_max) = rmdi,                                                     &
  min_gl_pft_io(npft_max) = rmdi,                                              & ! JBaguley
  min_rootc_pft_io(npft_max) = rmdi,                                           & ! JBaguley
  psi_close_io(npft_max) = rmdi,                                               &
  psi_open_io(npft_max) = rmdi,                                                &
  root_psi_crit_io(npft_max) = rmdi,                                           & ! JBaguley
  root_radi_pft_io(npft_max) = rmdi,                                           & ! JBaguley
  rootc_density_pft_io(npft_max) = rmdi                                          ! JBaguley
#endif

INTEGER ::                                                                     &
  c3_io(npft_max) = imdi,                                                      &
  irrig_pft_io(npft_max) = imdi,                                               &
  orient_io(npft_max) = imdi,                                                  &
  pft_conductance_model_io(npft_max) = imdi,                                     & ! JBaguley
  seg_root_vc_io(npft_max) = 0
      ! Root segment vulnerability curve (l_som_plant_segments):
      !   0: from p50_root_io / p88_root_io as given
      !   1: from p50_root_io alone, following Christoffersen et al.
      !      (2016, Geosci. Model Dev. 9: 4227-4255; TFS v.1-Hydro), who use
      !      the inverse polynomial of Manzoni et al. (2013a) for the
      !      fraction of maximum xylem conductivity (their Eqn 4),
      !        FMC_x(psi_x) = ( 1 + (psi_x / P50_x)^a_x )^-1,
      !      with the slope of the PLC curve at P50 from their tropical
      !      synthesis (Table 1),
      !        S = 54.4 (-P50 [MPa])^-1.17   (% MPa-1).
      !      For this FMC the PLC slope at P50 is 100 a / (4 |P50|), so
      !        a = 4 |P50| S / 100 = 2.176 |P50|^-0.17,
      !      and FMC = 0.12 at P88 = P50 (1/0.12 - 1)^(1/a). The segment's
      !      cumulative Weibull then passes through this P50 and P88 (exact
      !      there, approximate in the tails). p88_root_io must be unset.
      !      NOTE: S is fitted to bench-dehydration data for tropical
      !      upland trees, extrapolated here to roots; for resistant P50
      !      it gives much shallower curves than measured roots (see 2).
      !   2: from p50_root_io alone, with P88 from a fit to root
      !      vulnerability curves in the Xylem Functional Traits database
      !      (Choat et al. 2012 and updates; download 2026-09-28; root
      !      organ, air-injection methods excluded, n = 93):
      !        |P88| = 2.61 |P50|^0.70   (MPa),
      !      e.g. P50 -3 -> P88 -5.6, P50 -5 -> -8.1, P50 -8.5 -> -11.7 MPa
      !      (cf. the one Q. ilex root curve, Madrid, P50 -4.98 / P88 -9.61).
      !      The segment's cumulative Weibull passes through this P50/P88.
      !      p88_root_io must be unset.

REAL(KIND=real_jlslsm) ::                                                      &
  a_wl_io(npft_max) = rmdi,                                                    &
  a_ws_io(npft_max) = rmdi,                                                    &
  act_jmax_io(npft_max) = rmdi,                                                &
  act_vcmax_io(npft_max) = rmdi,                                               &
  aef_io(npft_max) = rmdi,                                                     &
  albsnc_max_io(npft_max) = rmdi,                                              &
  albsnc_min_io(npft_max) = rmdi,                                              &
  albsnf_max_io(npft_max) = rmdi,                                              &
  albsnf_maxl_io(npft_max) = rmdi,                                             &
  albsnf_maxu_io(npft_max) = rmdi,                                             &
  alpha_io(npft_max) = rmdi,                                                   &
  alpha_elec_io(npft_max) = rmdi,                                              &
  alnir_io(npft_max) = rmdi,                                                   &
  alnirl_io(npft_max) = rmdi,                                                  &
  alniru_io(npft_max) = rmdi,                                                  &
  alpar_io(npft_max) = rmdi,                                                   &
  alparl_io(npft_max) = rmdi,                                                  &
  alparu_io(npft_max) = rmdi,                                                  &
  avg_ba_io(npft_max) = rmdi,                                                  &
  b_wl_io(npft_max) = rmdi,                                                    &
  can_struct_a_io(npft_max) = rmdi,                                            &
  catch0_io(npft_max) = rmdi,                                                  &
  ccleaf_min_io(npft_max) = rmdi,                                              &
  ccleaf_max_io(npft_max) = rmdi,                                              &
  ccwood_min_io(npft_max) = rmdi,                                              &
  ccwood_max_io(npft_max) = rmdi,                                              &
  ci_st_io(npft_max) = rmdi,                                                   &
  dcatch_dlai_io(npft_max) = rmdi,                                             &
  deact_jmax_io(npft_max) = rmdi,                                              &
  deact_vcmax_io(npft_max) = rmdi,                                             &
  dfp_dcuo_io(npft_max) = rmdi,                                                &
  dgl_dm_io(npft_max) = rmdi,                                                  &
  dgl_dt_io(npft_max) = rmdi,                                                  &
  dqcrit_io(npft_max) = rmdi,                                                  &
  ds_jmax_io(npft_max) = rmdi,                                                 &
  ds_vcmax_io(npft_max) = rmdi,                                                &
  dust_veg_scj_io(npft_max) = rmdi,                                            &
  dz0v_dh_io(npft_max) = rmdi,                                                 &
  emis_pft_io(npft_max) = rmdi,                                                &
  eta_sl_io(npft_max) = rmdi,                                                  &
  f0_io(npft_max) = rmdi,                                                      &
  fef_bc_io(npft_max) = rmdi,                                                  &
  fef_c2h4_io(npft_max) = rmdi,                                                &
  fef_c2h6_io(npft_max) = rmdi,                                                &
  fef_c3h8_io(npft_max) = rmdi,                                                &
  fef_ch4_io(npft_max) = rmdi,                                                 &
  fef_co_io(npft_max) = rmdi,                                                  &
  fef_co2_io(npft_max) = rmdi,                                                 &
  fef_dms_io(npft_max) = rmdi,                                                 &
  fef_hcho_io(npft_max) = rmdi,                                                &
  fef_mecho_io(npft_max) = rmdi,                                               &
  fef_nh3_io(npft_max) = rmdi,                                                 &
  fef_nox_io(npft_max) = rmdi,                                                 &
  fef_oc_io(npft_max) = rmdi,                                                  &
  fef_so2_io(npft_max) = rmdi,                                                 &
  fd_io(npft_max) = rmdi,                                                      &
  fire_mort_io(npft_max) = rmdi,                                               &
  fl_o3_ct_io(npft_max) = rmdi,                                                &
  fsmc_of_io(npft_max) = rmdi,                                                 &
  fsmc_p0_io(npft_max) = rmdi,                                                 &
  sug_g0_io(npft_max) = rmdi,                                                  &
  g1_stomata_io(npft_max) = rmdi,                                              &
  g_leaf_0_io(npft_max) = rmdi,                                                &
  glmin_io(npft_max) = rmdi,                                                   &
  gpp_st_io(npft_max) = rmdi,                                                  &
  sug_grec_io(npft_max) = rmdi,                                                &
  gsoil_f_io(npft_max) = rmdi,                                                 &
  hw_sw_io(npft_max) = rmdi,                                                   &
  ief_io(npft_max) = rmdi,                                                     &
  infil_f_io(npft_max) = rmdi,                                                 &
  jv25_ratio_io(npft_max) = rmdi,                                              &
  kcrit_fractional_loss_io(npft_max) = rmdi,                                   & ! JBaguley
  kcrit_io(npft_max) = rmdi,                                                   & ! JBaguley
  kext_io(npft_max) = rmdi,                                                    &
  kmax_pft_io(npft_max) = rmdi,                                                & ! JBaguley
      ! Whole-plant hydraulic conductance per unit leaf area, in
      ! mmol m-2 s-1 MPa-1 (e.g. ~0.9 at FR-Pue from E / (LAI * (psi_pd -
      ! psi_md))). Converted to kmax_pft (mol m-2 s-1 Pa-1) below.
  kn_io(npft_max) = rmdi,                                                      &
  knl_io(npft_max) = rmdi,                                                     &
  kpar_io(npft_max) = rmdi,                                                    &
  lai_alb_lim_io(npft_max) = rmdi,                                             &
  lma_io(npft_max) = rmdi,                                                     &
  mef_io(npft_max) = rmdi,                                                     &
  neff_io(npft_max) = rmdi,                                                    &
  nl0_io(npft_max) = rmdi,                                                     &
  nmass_io(npft_max) = rmdi,                                                   &
  nr_nl_io(npft_max) = rmdi,                                                   &
  ns_nl_io(npft_max) = rmdi,                                                   &
  nsw_io(npft_max) = rmdi,                                                     &
  nr_io(npft_max) = rmdi,                                                      &
  omega_io(npft_max) = rmdi,                                                   &
  omegal_io(npft_max) = rmdi,                                                  &
  omegau_io(npft_max) = rmdi,                                                  &
  omnir_io(npft_max) = rmdi,                                                   &
  omnirl_io(npft_max) = rmdi,                                                  &
  omniru_io(npft_max) = rmdi,                                                  &
  p50_io(npft_max) = rmdi,                                                     & ! JBaguley
  p88_io(npft_max) = rmdi,                                                     & ! JBaguley
  ! Root / stem / leaf segments (l_som_plant_segments): resistance shares
  ! (defaults from Wang et al. 2019, root:stem:leaf kmax 1000:2000:2000) and
  ! optional segment P50/P88 (Pa; missing = the PFT's p50_io/p88_io).
  seg_frac_root_io(npft_max) = 0.5,                                            &
  seg_frac_stem_io(npft_max) = 0.25,                                           &
  seg_frac_leaf_io(npft_max) = 0.25,                                           &
  p50_root_io(npft_max) = rmdi,                                                &
  p50_stem_io(npft_max) = rmdi,                                                &
  p50_leaf_io(npft_max) = rmdi,                                                &
  p88_root_io(npft_max) = rmdi,                                                &
  p88_stem_io(npft_max) = rmdi,                                                &
  p88_leaf_io(npft_max) = rmdi,                                                &
  ! Cuticular leaf conductance (mmol H2O m-2 leaf s-1), the floor used when
  ! l_som_cuticular_floor (default 3, SurEau-Ecos Q. ilex, Ruffault 2022).
  gcut_io(npft_max) = 3.0,                                                     &
  ! Nonstomatal limitation (l_som_nsl): onset and zero point of the
  ! leaf-psi ramp on photosynthesis (Pa; missing = pftparm defaults 0 and
  ! -3 MPa, i.e. Dewar et al. 2022 Eqn 3(b)). Onset at the turgor loss
  ! point is an optional variant.
  psi_nsl_onset_io(npft_max) = rmdi,                                           &
  psi_nsl0_io(npft_max) = rmdi,                                                &
  ! Curvature exponent of the soil moisture stress factor (missing = 1,
  ! linear).
  fsmc_q_io(npft_max) = rmdi,                                                  &
  ! Soil-water down-regulation of Vcmax/Jmax (l_som_vcmax_psi; missing =
  ! pftparm defaults -2 MPa and 2 MPa-1).
  psi_vcmax_f_io(npft_max) = rmdi,                                             &
  sf_vcmax_io(npft_max) = rmdi,                                                &
  psi_vcmax_fmin_io(npft_max) = rmdi,                                          &
  nsl_sink_umax_io(npft_max) = rmdi,                                          &
  nsl_sink_tau_io(npft_max) = rmdi,                                           &
  nsl_sink_maint_io(npft_max) = rmdi,                                         &
  nsl_sink_psi50_io(npft_max) = rmdi,                                         &
  nsl_sink_sf_io(npft_max) = rmdi,                                            &
  ! DESICA (stomata_model = 5). Tuzet et al. (2003) closure
  ! fw = (1 + exp(sf psi_f)) / (1 + exp(sf (psi_f - psi_leaf))), with
  ! gs = g1 fw An / ca; defaults are the CABLE-DESICA evergreen broadleaf
  ! values (De Kauwe et al. 2020). Capacitances per unit leaf area (mmol H2O
  ! m-2 leaf MPa-1); leaf default from Xu et al. (2016) Table S3
  ! (1.5e-3 kg m-2 MPa-1), stem default ~Q. ilex (SurEau, Ruffault 2022).
  g1_tuzet_io(npft_max) = 4.19,                                                &
  sf_tuzet_io(npft_max) = 2.0,                                                 &
  psi_f_tuzet_io(npft_max) = -2.05e6,                                          &
  cap_leaf_io(npft_max) = 83.3,                                                &
  cap_stem_io(npft_max) = 3000.0,                                              &
  q10_leaf_io(npft_max) = rmdi,                                                &
  r_grow_io(npft_max) = rmdi,                                                  &
  rmass_io(npft_max) = rmdi,                                                   & ! JBaguley
  rootd_ft_io(npft_max) = rmdi,                                                &
  sigl_io(npft_max) = rmdi,                                                    &
  tef_io(npft_max) = rmdi,                                                     &
  tleaf_of_io(npft_max) = rmdi,                                                &
  tlow_io(npft_max) = rmdi,                                                    &
  tupp_io(npft_max) = rmdi,                                                    &
  vint_io(npft_max) = rmdi,                                                    &
  vsl_io(npft_max) = rmdi,                                                     &
  sug_yg_io(npft_max) = rmdi,                                                  &
  leaf_crit_io(npft_max) = rmdi,                                               & ! JBaguley
  z0hm_pft_io(npft_max) = rmdi,                                                &
  z0hm_classic_pft_io(npft_max) = rmdi,                                        &
  z0v_io(npft_max) = rmdi,                                                     &
  sox_a_io(npft_max) = rmdi,                                                   &
  sox_p50_io(npft_max) = rmdi,                                                 &
  sox_rp_min_io(npft_max) = rmdi
!---------------------------------------------------------------------
! Set up a namelist for reading and writing these arrays
!---------------------------------------------------------------------
NAMELIST  / jules_pftparm/                                                     &
#if !defined(UM_JULES)
  calc_rz_psi_io,  canht_ft_io,      lai_io,                                   & ! JBaguley
  fsmc_mod_io,     psi_close_io,     psi_open_io,                              &
  min_gl_pft_io,   min_rootc_pft_io, root_psi_crit_io,                         & ! JBaguley
  root_radi_pft_io,rootc_density_pft_io,                                       & ! JBaguley
#endif
  a_wl_io,         a_ws_io,          aef_io,                                   &
  act_jmax_io,     act_vcmax_io,     albsnc_max_io,                            &
  albsnc_min_io,   albsnf_max_io,    albsnf_maxl_io,                           &
  albsnf_maxu_io,  alpha_io,         alpha_elec_io,                            &
  alnir_io,        alnirl_io,        alniru_io,                                &
  alpar_io,        alparl_io,        alparu_io,                                &
  avg_ba_io,       b_wl_io,          c3_io,                                    &
  can_struct_a_io, catch0_io,        ccleaf_max_io,                            &
  ccleaf_min_io,   ccwood_max_io,    ccwood_min_io,                            &
  ci_st_io,        pft_conductance_model_io,                                   & ! JBaguley
  seg_root_vc_io,                                                              &
  dcatch_dlai_io,  deact_jmax_io,    deact_vcmax_io,                           &
  dfp_dcuo_io,     dgl_dm_io,        dgl_dt_io,                                &
  dqcrit_io,       ds_jmax_io,       ds_vcmax_io,                              &
  dust_veg_scj_io, dz0v_dh_io,       emis_pft_io,                              &
  eta_sl_io,       f0_io,            fef_bc_io,                                &
  fef_ch4_io,      fef_co_io,        fef_co2_io,                               &
  fef_nox_io,      fef_oc_io,        fef_so2_io,                               &
  fef_c2h4_io,     fef_c2h6_io,      fef_c3h8_io,                              &
  fef_hcho_io,     fef_mecho_io,                                               &
  fef_nh3_io,      fef_dms_io,                                                 &
  fd_io,           fire_mort_io,     fl_o3_ct_io,                              &
  fsmc_of_io,      fsmc_p0_io,       sug_g0_io,                                &
  g1_stomata_io,   g_leaf_0_io,      glmin_io,                                 &
  gpp_st_io,       sug_grec_io,      gsoil_f_io,                               &
  hw_sw_io,        ief_io,           infil_f_io,                               &
  irrig_pft_io,                                                                &
  jv25_ratio_io,   kcrit_fractional_loss_io,                                   & ! JBaguley
  kext_io,         kmax_pft_io,                                                & ! JBaguley
  kn_io,                                    &
  knl_io,          kpar_io,          lai_alb_lim_io,                           &
  lma_io,          mef_io,           neff_io,                                  &
  nl0_io,          nmass_io,         nr_io,                                    &
  nr_nl_io,        ns_nl_io,         nsw_io,                                   &
  omega_io,        omegal_io,        omegau_io,                                &
  omnir_io,        omnirl_io,        omniru_io,                                &
  orient_io,       p50_io,           p88_io,                                   & ! JBaguley
  seg_frac_root_io, seg_frac_stem_io, seg_frac_leaf_io,                        &
  p50_root_io,     p50_stem_io,      p50_leaf_io,                              &
  p88_root_io,     p88_stem_io,      p88_leaf_io,                              &
  gcut_io,         psi_nsl_onset_io, psi_nsl0_io,                              &
  fsmc_q_io,       psi_vcmax_f_io,   sf_vcmax_io,      psi_vcmax_fmin_io,      &
  nsl_sink_umax_io, nsl_sink_tau_io, nsl_sink_maint_io,                    &
  nsl_sink_psi50_io, nsl_sink_sf_io,                                  &
  g1_tuzet_io,     sf_tuzet_io,      psi_f_tuzet_io,                           &
  cap_leaf_io,     cap_stem_io,                                                &
  q10_leaf_io,      r_grow_io,                                &
  rmass_io,        rootd_ft_io,      sigl_io,                                  & !JBaguley
  tef_io,          tleaf_of_io,      tlow_io,                                  &
  tupp_io,         vint_io,          vsl_io,                                   &
  sug_yg_io,       leaf_crit_io,     z0hm_pft_io,                              & !JBaguley
  z0hm_classic_pft_io,               z0v_io,                                   &
  sox_a_io,        sox_p50_io,       sox_rp_min_io

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PFTPARM_IO'

CONTAINS

#if defined(UM_JULES)
SUBROUTINE read_nml_jules_pftparm (unitnumber)

! Description:
!  Read the JULES_PFTPARM namelist

USE setup_namelist, ONLY: setup_nml_type
USE check_iostat_mod, ONLY:  check_iostat
USE UM_parcore,       ONLY:  mype
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: unitnumber
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_JULES_PFTPARM'
INTEGER(KIND=jpim), PARAMETER          :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER          :: zhook_out = 1
CHARACTER(LEN=errormessagelength) :: iomessage

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 2
INTEGER, PARAMETER :: n_int = 5 * npft_max ! = the INTEGER arrays in my_namelist
INTEGER, PARAMETER :: n_real = 140 * npft_max ! = the REAL arrays in my_namelist

TYPE :: my_namelist
  SEQUENCE
  INTEGER :: c3_io(npft_max)
  INTEGER :: irrig_pft_io(npft_max)
  INTEGER :: orient_io(npft_max)
  INTEGER :: pft_conductance_model_io(npft_max) ! JBaguley
  INTEGER :: seg_root_vc_io(npft_max)
  REAL(KIND=real_jlslsm) :: a_wl_io(npft_max)
  REAL(KIND=real_jlslsm) :: a_ws_io(npft_max)
  REAL(KIND=real_jlslsm) :: act_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: act_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: aef_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnc_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnc_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_maxl_io(npft_max)
  REAL(KIND=real_jlslsm) :: albsnf_maxu_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpha_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpha_elec_io(npft_max)
  REAL(KIND=real_jlslsm) :: alnir_io(npft_max)
  REAL(KIND=real_jlslsm) :: alnirl_io(npft_max)
  REAL(KIND=real_jlslsm) :: alniru_io(npft_max)
  REAL(KIND=real_jlslsm) :: alpar_io(npft_max)
  REAL(KIND=real_jlslsm) :: alparl_io(npft_max)
  REAL(KIND=real_jlslsm) :: alparu_io(npft_max)
  REAL(KIND=real_jlslsm) :: avg_ba_io(npft_max)
  REAL(KIND=real_jlslsm) :: b_wl_io(npft_max)
  REAL(KIND=real_jlslsm) :: can_struct_a_io(npft_max)
  REAL(KIND=real_jlslsm) :: catch0_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccleaf_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccleaf_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccwood_min_io(npft_max)
  REAL(KIND=real_jlslsm) :: ccwood_max_io(npft_max)
  REAL(KIND=real_jlslsm) :: ci_st_io(npft_max)
  REAL(KIND=real_jlslsm) :: dcatch_dlai_io(npft_max)
  REAL(KIND=real_jlslsm) :: deact_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: deact_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: dfp_dcuo_io(npft_max)
  REAL(KIND=real_jlslsm) :: dgl_dm_io(npft_max)
  REAL(KIND=real_jlslsm) :: dgl_dt_io(npft_max)
  REAL(KIND=real_jlslsm) :: dqcrit_io(npft_max)
  REAL(KIND=real_jlslsm) :: ds_jmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: ds_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: dust_veg_scj_io(npft_max)
  REAL(KIND=real_jlslsm) :: dz0v_dh_io(npft_max)
  REAL(KIND=real_jlslsm) :: emis_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: eta_sl_io(npft_max)
  REAL(KIND=real_jlslsm) :: f0_io(npft_max)
  REAL(KIND=real_jlslsm) :: fd_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_bc_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_ch4_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_co_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_co2_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_nox_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_oc_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_so2_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_c2h4_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_c2h6_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_c3h8_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_hcho_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_mecho_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_nh3_io(npft_max)
  REAL(KIND=real_jlslsm) :: fef_dms_io(npft_max)
  REAL(KIND=real_jlslsm) :: fire_mort_io(npft_max)
  REAL(KIND=real_jlslsm) :: fl_o3_ct_io(npft_max)
  REAL(KIND=real_jlslsm) :: fsmc_of_io(npft_max)
  REAL(KIND=real_jlslsm) :: fsmc_p0_io(npft_max)
  REAL(KIND=real_jlslsm) :: sug_g0_io(npft_max)
  REAL(KIND=real_jlslsm) :: g1_stomata_io(npft_max)
  REAL(KIND=real_jlslsm) :: g_leaf_0_io(npft_max)
  REAL(KIND=real_jlslsm) :: glmin_io(npft_max)
  REAL(KIND=real_jlslsm) :: gpp_st_io(npft_max)
  REAL(KIND=real_jlslsm) :: sug_grec_io(npft_max)
  REAL(KIND=real_jlslsm) :: gsoil_f_io(npft_max)
  REAL(KIND=real_jlslsm) :: hw_sw_io(npft_max)
  REAL(KIND=real_jlslsm) :: ief_io(npft_max)
  REAL(KIND=real_jlslsm) :: infil_f_io(npft_max)
  REAL(KIND=real_jlslsm) :: jv25_ratio_io(npft_max)
  REAL(KIND=real_jlslsm) :: kcrit_fractional_loss_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: kext_io(npft_max)
  REAL(KIND=real_jlslsm) :: kmax_pft_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: kn_io(npft_max)
  REAL(KIND=real_jlslsm) :: knl_io(npft_max)
  REAL(KIND=real_jlslsm) :: kpar_io(npft_max)
  REAL(KIND=real_jlslsm) :: lai_alb_lim_io(npft_max)
  REAL(KIND=real_jlslsm) :: lma_io(npft_max)
  REAL(KIND=real_jlslsm) :: mef_io(npft_max)
  REAL(KIND=real_jlslsm) :: neff_io(npft_max)
  REAL(KIND=real_jlslsm) :: nl0_io(npft_max)
  REAL(KIND=real_jlslsm) :: nmass_io(npft_max)
  REAL(KIND=real_jlslsm) :: nr_io(npft_max)
  REAL(KIND=real_jlslsm) :: nr_nl_io(npft_max)
  REAL(KIND=real_jlslsm) :: ns_nl_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsw_io(npft_max)
  REAL(KIND=real_jlslsm) :: omega_io(npft_max)
  REAL(KIND=real_jlslsm) :: omegal_io(npft_max)
  REAL(KIND=real_jlslsm) :: omegau_io(npft_max)
  REAL(KIND=real_jlslsm) :: omnir_io(npft_max)
  REAL(KIND=real_jlslsm) :: omnirl_io(npft_max)
  REAL(KIND=real_jlslsm) :: omniru_io(npft_max)
  REAL(KIND=real_jlslsm) :: p50_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: p88_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: seg_frac_root_io(npft_max)
  REAL(KIND=real_jlslsm) :: gcut_io(npft_max)
  REAL(KIND=real_jlslsm) :: psi_nsl_onset_io(npft_max)
  REAL(KIND=real_jlslsm) :: psi_nsl0_io(npft_max)
  REAL(KIND=real_jlslsm) :: fsmc_q_io(npft_max)
  REAL(KIND=real_jlslsm) :: psi_vcmax_f_io(npft_max)
  REAL(KIND=real_jlslsm) :: sf_vcmax_io(npft_max)
  REAL(KIND=real_jlslsm) :: psi_vcmax_fmin_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsl_sink_umax_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsl_sink_tau_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsl_sink_maint_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsl_sink_psi50_io(npft_max)
  REAL(KIND=real_jlslsm) :: nsl_sink_sf_io(npft_max)
  REAL(KIND=real_jlslsm) :: g1_tuzet_io(npft_max)
  REAL(KIND=real_jlslsm) :: sf_tuzet_io(npft_max)
  REAL(KIND=real_jlslsm) :: psi_f_tuzet_io(npft_max)
  REAL(KIND=real_jlslsm) :: cap_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: cap_stem_io(npft_max)
  REAL(KIND=real_jlslsm) :: seg_frac_stem_io(npft_max)
  REAL(KIND=real_jlslsm) :: seg_frac_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: p50_root_io(npft_max)
  REAL(KIND=real_jlslsm) :: p50_stem_io(npft_max)
  REAL(KIND=real_jlslsm) :: p50_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: p88_root_io(npft_max)
  REAL(KIND=real_jlslsm) :: p88_stem_io(npft_max)
  REAL(KIND=real_jlslsm) :: p88_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: q10_leaf_io(npft_max)
  REAL(KIND=real_jlslsm) :: r_grow_io(npft_max)
  REAL(KIND=real_jlslsm) :: rmass_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: rootd_ft_io(npft_max)
  REAL(KIND=real_jlslsm) :: sigl_io(npft_max)
  REAL(KIND=real_jlslsm) :: tef_io(npft_max)
  REAL(KIND=real_jlslsm) :: tleaf_of_io(npft_max)
  REAL(KIND=real_jlslsm) :: tlow_io(npft_max)
  REAL(KIND=real_jlslsm) :: tupp_io(npft_max)
  REAL(KIND=real_jlslsm) :: vint_io(npft_max)
  REAL(KIND=real_jlslsm) :: vsl_io(npft_max)
  REAL(KIND=real_jlslsm) :: sug_yg_io(npft_max)
  REAL(KIND=real_jlslsm) :: leaf_crit_io(npft_max) ! JBaguley
  REAL(KIND=real_jlslsm) :: z0hm_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: z0hm_classic_pft_io(npft_max)
  REAL(KIND=real_jlslsm) :: z0v_io(npft_max)
  REAL(KIND=real_jlslsm) :: sox_a_io(npft_max)
  REAL(KIND=real_jlslsm) :: sox_p50_io(npft_max)
  REAL(KIND=real_jlslsm) :: sox_rp_min_io(npft_max)
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in = n_int,               &
                    n_real_in = n_real)

IF (mype == 0) THEN

  READ (UNIT = unitnumber, NML = jules_pftparm, IOSTAT = errorstatus,          &
        IOMSG = iomessage)
  CALL check_iostat(errorstatus, "namelist jules_pftparm",                     &
           iomessage)

  my_nml % a_wl_io        = a_wl_io
  my_nml % a_ws_io        = a_ws_io
  my_nml % act_jmax_io    = act_jmax_io
  my_nml % act_vcmax_io   = act_vcmax_io
  my_nml % aef_io         = aef_io
  my_nml % albsnc_max_io  = albsnc_max_io
  my_nml % albsnc_min_io  = albsnc_min_io
  my_nml % albsnf_max_io  = albsnf_max_io
  my_nml % albsnf_maxl_io = albsnf_maxl_io
  my_nml % albsnf_maxu_io = albsnf_maxu_io
  my_nml % alpha_io       = alpha_io
  my_nml % alpha_elec_io  = alpha_elec_io
  my_nml % alnir_io       = alnir_io
  my_nml % alnirl_io      = alnirl_io
  my_nml % alniru_io      = alniru_io
  my_nml % alpar_io       = alpar_io
  my_nml % alparl_io      = alparl_io
  my_nml % alparu_io      = alparu_io
  my_nml % avg_ba_io      = avg_ba_io
  my_nml % b_wl_io        = b_wl_io
  my_nml % c3_io          = c3_io
  my_nml % can_struct_a_io = can_struct_a_io
  my_nml % catch0_io      = catch0_io
  my_nml % ccleaf_min_io  = ccleaf_min_io
  my_nml % ccleaf_max_io  = ccleaf_max_io
  my_nml % ccwood_min_io  = ccwood_min_io
  my_nml % ccwood_max_io  = ccwood_max_io
  my_nml % ci_st_io       = ci_st_io
  my_nml % pft_conductance_model_io = pft_conductance_model_io ! JBaguley
  my_nml % seg_root_vc_io = seg_root_vc_io
  my_nml % dcatch_dlai_io = dcatch_dlai_io
  my_nml % deact_jmax_io  = deact_jmax_io
  my_nml % deact_vcmax_io = deact_vcmax_io
  my_nml % dfp_dcuo_io    = dfp_dcuo_io
  my_nml % dgl_dm_io      = dgl_dm_io
  my_nml % dgl_dt_io      = dgl_dt_io
  my_nml % dqcrit_io      = dqcrit_io
  my_nml % ds_jmax_io     = ds_jmax_io
  my_nml % ds_vcmax_io    = ds_vcmax_io
  my_nml % dust_veg_scj_io = dust_veg_scj_io
  my_nml % dz0v_dh_io     = dz0v_dh_io
  my_nml % emis_pft_io    = emis_pft_io
  my_nml % eta_sl_io      = eta_sl_io
  my_nml % f0_io          = f0_io
  my_nml % fd_io          = fd_io
  my_nml % fef_bc_io      = fef_bc_io
  my_nml % fef_ch4_io     = fef_ch4_io
  my_nml % fef_co_io      = fef_co_io
  my_nml % fef_co2_io     = fef_co2_io
  my_nml % fef_nox_io     = fef_nox_io
  my_nml % fef_oc_io      = fef_oc_io
  my_nml % fef_so2_io     = fef_so2_io
  my_nml % fef_c2h4_io    = fef_c2h4_io
  my_nml % fef_c2h6_io    = fef_c2h6_io
  my_nml % fef_c3h8_io    = fef_c3h8_io
  my_nml % fef_hcho_io    = fef_hcho_io
  my_nml % fef_mecho_io   = fef_mecho_io
  my_nml % fef_nh3_io     = fef_nh3_io
  my_nml % fef_dms_io     = fef_dms_io
  my_nml % fire_mort_io   = fire_mort_io
  my_nml % fl_o3_ct_io    = fl_o3_ct_io
  my_nml % fsmc_of_io     = fsmc_of_io
  my_nml % fsmc_p0_io     = fsmc_p0_io
  my_nml % sug_g0_io      = sug_g0_io
  my_nml % g1_stomata_io  = g1_stomata_io
  my_nml % g_leaf_0_io    = g_leaf_0_io
  my_nml % glmin_io       = glmin_io
  my_nml % gpp_st_io      = gpp_st_io
  my_nml % sug_grec_io    = sug_grec_io
  my_nml % gsoil_f_io     = gsoil_f_io
  my_nml % hw_sw_io       = hw_sw_io
  my_nml % ief_io         = ief_io
  my_nml % infil_f_io     = infil_f_io
  my_nml % irrig_pft_io   = irrig_pft_io
  my_nml % jv25_ratio_io  = jv25_ratio_io
  my_nml % kcrit_fractional_loss_io = kcrit_fractional_loss_io ! JBaguley
  my_nml % kext_io        = kext_io
  my_nml % kmax_pft_io    = kmax_pft_io ! JBaguley
  my_nml % kn_io          = kn_io
  my_nml % knl_io         = knl_io
  my_nml % kpar_io        = kpar_io
  my_nml % lai_alb_lim_io = lai_alb_lim_io
  my_nml % lma_io         = lma_io
  my_nml % mef_io         = mef_io
  my_nml % neff_io        = neff_io
  my_nml % nl0_io         = nl0_io
  my_nml % nmass_io       = nmass_io
  my_nml % nr_io          = nr_io
  my_nml % nr_nl_io       = nr_nl_io
  my_nml % ns_nl_io       = ns_nl_io
  my_nml % nsw_io         = nsw_io
  my_nml % omega_io       = omega_io
  my_nml % omegal_io      = omegal_io
  my_nml % omegau_io      = omegau_io
  my_nml % omnir_io       = omnir_io
  my_nml % omnirl_io      = omnirl_io
  my_nml % omniru_io      = omniru_io
  my_nml % orient_io      = orient_io
  my_nml % p50_io         = p50_io ! JBaguley
  my_nml % p88_io         = p88_io ! JBaguley
  my_nml % seg_frac_root_io = seg_frac_root_io
  my_nml % gcut_io        = gcut_io
  my_nml % psi_nsl_onset_io = psi_nsl_onset_io
  my_nml % psi_nsl0_io    = psi_nsl0_io
  my_nml % fsmc_q_io      = fsmc_q_io
  my_nml % psi_vcmax_f_io = psi_vcmax_f_io
  my_nml % sf_vcmax_io    = sf_vcmax_io
  my_nml % psi_vcmax_fmin_io = psi_vcmax_fmin_io
  my_nml % nsl_sink_umax_io = nsl_sink_umax_io
  my_nml % nsl_sink_tau_io = nsl_sink_tau_io
  my_nml % nsl_sink_maint_io = nsl_sink_maint_io
  my_nml % nsl_sink_psi50_io = nsl_sink_psi50_io
  my_nml % nsl_sink_sf_io = nsl_sink_sf_io
  my_nml % g1_tuzet_io    = g1_tuzet_io
  my_nml % sf_tuzet_io    = sf_tuzet_io
  my_nml % psi_f_tuzet_io = psi_f_tuzet_io
  my_nml % cap_leaf_io    = cap_leaf_io
  my_nml % cap_stem_io    = cap_stem_io
  my_nml % seg_frac_stem_io = seg_frac_stem_io
  my_nml % seg_frac_leaf_io = seg_frac_leaf_io
  my_nml % p50_root_io    = p50_root_io
  my_nml % p50_stem_io    = p50_stem_io
  my_nml % p50_leaf_io    = p50_leaf_io
  my_nml % p88_root_io    = p88_root_io
  my_nml % p88_stem_io    = p88_stem_io
  my_nml % p88_leaf_io    = p88_leaf_io
  my_nml % q10_leaf_io    = q10_leaf_io
  my_nml % r_grow_io      = r_grow_io
  my_nml % rmass_io       = rmass_io ! JBaguley
  my_nml % rootd_ft_io    = rootd_ft_io
  my_nml % sigl_io        = sigl_io
  my_nml % tef_io         = tef_io
  my_nml % tleaf_of_io    = tleaf_of_io
  my_nml % tlow_io        = tlow_io
  my_nml % tupp_io        = tupp_io
  my_nml % vint_io        = vint_io
  my_nml % vsl_io         = vsl_io
  my_nml % sug_yg_io      = sug_yg_io
  my_nml % leaf_crit_io   = leaf_crit_io
  my_nml % z0hm_pft_io    = z0hm_pft_io
  my_nml % z0hm_classic_pft_io = z0hm_classic_pft_io
  my_nml % z0v_io         = z0v_io
  my_nml % sox_a_io       = sox_a_io
  my_nml % sox_p50_io     = sox_p50_io
  my_nml % sox_rp_min_io  = sox_rp_min_io
END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  a_wl_io         = my_nml % a_wl_io
  a_ws_io         = my_nml % a_ws_io
  act_jmax_io     = my_nml % act_jmax_io
  act_vcmax_io    = my_nml % act_vcmax_io
  aef_io          = my_nml % aef_io
  albsnc_max_io   = my_nml % albsnc_max_io
  albsnc_min_io   = my_nml % albsnc_min_io
  albsnf_max_io   = my_nml % albsnf_max_io
  albsnf_maxl_io  = my_nml % albsnf_maxl_io
  albsnf_maxu_io  = my_nml % albsnf_maxu_io
  alpha_io        = my_nml % alpha_io
  alpha_elec_io   = my_nml % alpha_elec_io
  alnir_io        = my_nml % alnir_io
  alnirl_io       = my_nml % alnirl_io
  alniru_io       = my_nml % alniru_io
  alpar_io        = my_nml % alpar_io
  alparl_io       = my_nml % alparl_io
  alparu_io       = my_nml % alparu_io
  avg_ba_io       = my_nml % avg_ba_io
  b_wl_io         = my_nml % b_wl_io
  c3_io           = my_nml % c3_io
  can_struct_a_io = my_nml % can_struct_a_io
  catch0_io       = my_nml % catch0_io
  ccleaf_min_io   = my_nml % ccleaf_min_io
  ccleaf_max_io   = my_nml % ccleaf_max_io
  ccwood_min_io   = my_nml % ccwood_min_io
  ccwood_max_io   = my_nml % ccwood_max_io
  ci_st_io        = my_nml % ci_st_io
  pft_conductance_model_io = my_nml % pft_conductance_model_io ! JBaguley
  seg_root_vc_io = my_nml % seg_root_vc_io
  dcatch_dlai_io  = my_nml % dcatch_dlai_io
  deact_jmax_io   = my_nml % deact_jmax_io
  deact_vcmax_io  = my_nml % deact_vcmax_io
  dfp_dcuo_io     = my_nml % dfp_dcuo_io
  dgl_dm_io       = my_nml % dgl_dm_io
  dgl_dt_io       = my_nml % dgl_dt_io
  dqcrit_io       = my_nml % dqcrit_io
  ds_jmax_io      = my_nml % ds_jmax_io
  ds_vcmax_io     = my_nml % ds_vcmax_io
  dust_veg_scj_io = my_nml % dust_veg_scj_io
  dz0v_dh_io      = my_nml % dz0v_dh_io
  emis_pft_io     = my_nml % emis_pft_io
  eta_sl_io       = my_nml % eta_sl_io
  f0_io           = my_nml % f0_io
  fd_io           = my_nml % fd_io
  fef_bc_io       = my_nml % fef_bc_io
  fef_ch4_io      = my_nml % fef_ch4_io
  fef_co_io       = my_nml % fef_co_io
  fef_co2_io      = my_nml % fef_co2_io
  fef_nox_io      = my_nml % fef_nox_io
  fef_oc_io       = my_nml % fef_oc_io
  fef_so2_io      = my_nml % fef_so2_io
  fef_c2h4_io     = my_nml % fef_c2h4_io
  fef_c2h6_io     = my_nml % fef_c2h6_io
  fef_c3h8_io     = my_nml % fef_c3h8_io
  fef_hcho_io     = my_nml % fef_hcho_io
  fef_mecho_io    = my_nml % fef_mecho_io
  fef_nh3_io      = my_nml % fef_nh3_io
  fef_dms_io      = my_nml % fef_dms_io
  fire_mort_io    = my_nml % fire_mort_io
  fl_o3_ct_io     = my_nml % fl_o3_ct_io
  fsmc_of_io      = my_nml % fsmc_of_io
  fsmc_p0_io      = my_nml % fsmc_p0_io
  g1_stomata_io   = my_nml % g1_stomata_io
  sug_g0_io       = my_nml % sug_g0_io
  g_leaf_0_io     = my_nml % g_leaf_0_io
  glmin_io        = my_nml % glmin_io
  gpp_st_io       = my_nml % gpp_st_io
  sug_grec_io     = my_nml % sug_grec_io
  gsoil_f_io      = my_nml % gsoil_f_io
  hw_sw_io        = my_nml % hw_sw_io
  ief_io          = my_nml % ief_io
  infil_f_io      = my_nml % infil_f_io
  irrig_pft_io    = my_nml % irrig_pft_io
  jv25_ratio_io   = my_nml % jv25_ratio_io
  kcrit_fractional_loss_io = my_nml % kcrit_fractional_loss_io ! JBaguley
  kext_io         = my_nml % kext_io
  kmax_pft_io     = my_nml % kmax_pft_io ! JBaguley
  kn_io           = my_nml % kn_io
  knl_io          = my_nml % knl_io
  kpar_io         = my_nml % kpar_io
  lai_alb_lim_io  = my_nml % lai_alb_lim_io
  lma_io          = my_nml % lma_io
  mef_io          = my_nml % mef_io
  neff_io         = my_nml % neff_io
  nl0_io          = my_nml % nl0_io
  nmass_io        = my_nml % nmass_io
  nr_io           = my_nml % nr_io
  nr_nl_io        = my_nml % nr_nl_io
  ns_nl_io        = my_nml % ns_nl_io
  nsw_io          = my_nml % nsw_io
  omega_io        = my_nml % omega_io
  omegal_io       = my_nml % omegal_io
  omegau_io       = my_nml % omegau_io
  omnir_io        = my_nml % omnir_io
  omnirl_io       = my_nml % omnirl_io
  omniru_io       = my_nml % omniru_io
  orient_io       = my_nml % orient_io
  p50_io          = my_nml % p50_io ! JBaguley
  p88_io          = my_nml % p88_io ! JBaguley
  seg_frac_root_io = my_nml % seg_frac_root_io
  gcut_io         = my_nml % gcut_io
  psi_nsl_onset_io = my_nml % psi_nsl_onset_io
  psi_nsl0_io     = my_nml % psi_nsl0_io
  fsmc_q_io       = my_nml % fsmc_q_io
  psi_vcmax_f_io  = my_nml % psi_vcmax_f_io
  sf_vcmax_io     = my_nml % sf_vcmax_io
  psi_vcmax_fmin_io = my_nml % psi_vcmax_fmin_io
  nsl_sink_umax_io = my_nml % nsl_sink_umax_io
  nsl_sink_tau_io = my_nml % nsl_sink_tau_io
  nsl_sink_maint_io = my_nml % nsl_sink_maint_io
  nsl_sink_psi50_io = my_nml % nsl_sink_psi50_io
  nsl_sink_sf_io = my_nml % nsl_sink_sf_io
  g1_tuzet_io     = my_nml % g1_tuzet_io
  sf_tuzet_io     = my_nml % sf_tuzet_io
  psi_f_tuzet_io  = my_nml % psi_f_tuzet_io
  cap_leaf_io     = my_nml % cap_leaf_io
  cap_stem_io     = my_nml % cap_stem_io
  seg_frac_stem_io = my_nml % seg_frac_stem_io
  seg_frac_leaf_io = my_nml % seg_frac_leaf_io
  p50_root_io     = my_nml % p50_root_io
  p50_stem_io     = my_nml % p50_stem_io
  p50_leaf_io     = my_nml % p50_leaf_io
  p88_root_io     = my_nml % p88_root_io
  p88_stem_io     = my_nml % p88_stem_io
  p88_leaf_io     = my_nml % p88_leaf_io
  q10_leaf_io     = my_nml % q10_leaf_io
  r_grow_io       = my_nml % r_grow_io
  rmass_io        = my_nml % rmass_io ! JBaguley
  rootd_ft_io     = my_nml % rootd_ft_io
  sigl_io         = my_nml % sigl_io
  tef_io          = my_nml % tef_io
  tleaf_of_io     = my_nml % tleaf_of_io
  tlow_io         = my_nml % tlow_io
  tupp_io         = my_nml % tupp_io
  vint_io         = my_nml % vint_io
  vsl_io          = my_nml % vsl_io
  sug_yg_io       = my_nml % sug_yg_io
  leaf_crit_io    = my_nml % leaf_crit_io
  z0hm_pft_io     = my_nml % z0hm_pft_io
  z0hm_classic_pft_io = my_nml % z0hm_classic_pft_io
  z0v_io          = my_nml % z0v_io
  sox_a_io        = my_nml % sox_a_io
  sox_p50_io      = my_nml % sox_p50_io
  sox_rp_min_io   = my_nml % sox_rp_min_io
END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_nml_jules_pftparm
#endif


SUBROUTINE init_pftparm_allocated()

!No USE statements other than Dr Hook
USE ereport_mod, ONLY: ereport
USE parkind1,    ONLY: jprb, jpim
USE yomhook,     ONLY: lhook, dr_hook

USE pftparm, ONLY:                                                             &
! namelist variables:
#if !defined(UM_JULES)
  calc_rz_psi,     fsmc_mod,         psi_close,                                & ! JBaguley
  psi_open,        min_rootc_pft,    root_psi_crit,                            & ! JBaguley
  root_radi_pft,   rootc_density_pft, min_gl_pft,                              & ! JBaguley
  gcut,            psi_nsl_onset,    psi_nsl0,         fsmc_q,                 &
  psi_vcmax_f,     sf_vcmax,         psi_vcmax_fmin,                           &
  nsl_sink_umax, nsl_sink_tau, nsl_sink_maint,                       &
  nsl_sink_psi50, nsl_sink_sf,                                     &
  g1_tuzet,        sf_tuzet,                                                   &
  psi_f_tuzet,     cap_leaf,         cap_stem,                                 &
#endif
  a_wl,            a_ws,             aef,                                      &
  act_jmax,        act_vcmax,        albsnc_max,                               &
  albsnc_min,      albsnf_max,       albsnf_maxl,                              &
  albsnf_maxu,     alpha,            alpha_elec,                               &
  alnir,           alnirl,           alniru,                                   &
  alpar,           alparl,           alparu,                                   &
  avg_ba,          b_wl,             c3,                                       &
  can_struct_a,    catch0,           ccleaf_max,                               &
  ccleaf_min,      ccwood_max,       ccwood_min,                               &
  ci_st,           pft_conductance_model,                                      & ! JBaguley
  conductance_b,   conductance_c,                                              & ! JBaguley
  seg_kfac,        conductance_b_seg, conductance_c_seg,                       &
  dcatch_dlai,                                                                 &
  deact_jmax,      deact_vcmax,      dfp_dcuo,                                 &
  dgl_dm,          dgl_dt,           dqcrit,                                   &
  ds_jmax,         ds_vcmax,         dust_veg_scj,                             &
  dz0v_dh,         emis_pft,         eta_sl,                                   &
  f0,              fd,               fef_bc,                                   &
  fef_ch4,         fef_co,           fef_co2,                                  &
  fef_nox,         fef_oc,           fef_so2,                                  &
  fef_c2h4,        fef_c2h6,         fef_c3h8,                                 &
  fef_hcho,        fef_mecho,        fef_nh3,                                  &
  fef_dms,                                                                     &
  fire_mort,       fl_o3_ct,         fsmc_of,                                  &
  fsmc_p0,         sug_g0,           g1_stomata,                               &
  g_leaf_0,        glmin,            gpp_st,                                   &
  sug_grec,        gsoil_f,          hw_sw,                                    &
  ief,             infil_f,          jv25_ratio,                               &
  kcrit_fractional_loss,             kcrit,                                    & ! JBaguley
  kext,            kmax_pft,                                                   & ! JBaguley
  kn,               knl,                                                       &
  kpar,            lai_alb_lim,      lma,                                      &
  mef,             neff,             nl0,                                      &
  nmass,           nr,               nr_nl,                                    &
  ns_nl,           nsw,              omega,                                    &
  omegal,          omegau,           omnir,                                    &
  omnirl,          omniru,           orient,                                   &
  P50,             P88,                                                        & ! JBaguley
  q10_leaf,        r_grow,           rmass,                                    &
  rootd_ft,        sigl,             tef,                                      & !JBaguley
  tleaf_of,        tlow,             tupp,                                     &
  vint,            vsl,              sug_yg,                                   &
  leaf_crit,       z0v,                                                        & !JBaguley
  sox_a,           sox_p50,          sox_rp_min



USE c_irrigation_mod, ONLY: irrig_tile
USE c_z0h_z0m,    ONLY: z0h_z0m,  z0h_z0m_classic

USE jules_surface_types_mod, ONLY: npft

USE jules_vegetation_mod, ONLY: l_som_plant_segments, l_som_nsl,              &
                                l_som_root_supply, l_som_vcmax_psi,            &
                                l_som_nsl_sink
USE jules_soil_mod, ONLY: l_bound_soil_wp

IMPLICIT NONE

INTEGER(KIND=jpim) :: i = 0
INTEGER :: errcode
INTEGER :: iseg
REAL(KIND=real_jlslsm) :: seg_frac(3), p50_seg(3), p88_seg(3)
REAL(KIND=real_jlslsm) :: a_root   ! Christoffersen root curve shape (seg_root_vc_io = 1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_PFTPARM_ALLOCATED'

!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Radiation and albedo parameters.
albsnc_max(:)   = albsnc_max_io(1:npft)
albsnc_min(:)   = albsnc_min_io(1:npft)
albsnf_max(:)   = albsnf_max_io(1:npft)
albsnf_maxl(:)  = albsnf_maxl_io(1:npft)
albsnf_maxu(:)  = albsnf_maxu_io(1:npft)
alnir(:)        = alnir_io(1:npft)
alnirl(:)       = alnirl_io(1:npft)
alniru(:)       = alniru_io(1:npft)
alpar(:)        = alpar_io(1:npft)
alparl(:)       = alparl_io(1:npft)
alparu(:)       = alparu_io(1:npft)
kext(:)         = kext_io(1:npft)
kpar(:)         = kpar_io(1:npft)
lai_alb_lim(:)  = lai_alb_lim_io(1:npft)
omega(:)        = omega_io(1:npft)
omegal(:)       = omegal_io(1:npft)
omegau(:)       = omegau_io(1:npft)
omnir(:)        = omnir_io(1:npft)
omnirl(:)       = omnirl_io(1:npft)
omniru(:)       = omniru_io(1:npft)
orient(:)       = orient_io(1:npft)

! Photosynthesis and respiration parameters.
act_jmax(:)   = act_jmax_io(1:npft)
act_vcmax(:)  = act_vcmax_io(1:npft)
alpha(:)        = alpha_io(1:npft)
alpha_elec(:) = alpha_elec_io(1:npft)
c3(:)           = c3_io(1:npft)
can_struct_a(:) = can_struct_a_io(1:npft)
deact_jmax(:) = deact_jmax_io(1:npft)
deact_vcmax(:)= deact_vcmax_io(1:npft)
dqcrit(:)       = dqcrit_io(1:npft)
ds_jmax(:)    = ds_jmax_io(1:npft)
ds_vcmax(:)   = ds_vcmax_io(1:npft)
f0(:)           = f0_io(1:npft)
fd(:)           = fd_io(1:npft)
g1_stomata(:) = g1_stomata_io(1:npft)
jv25_ratio(:) = jv25_ratio_io(1:npft)
kn(:)           = kn_io(1:npft)
knl(:)          = knl_io(1:npft)
neff(:)         = neff_io(1:npft)
nl0(:)          = nl0_io(1:npft)
nr_nl(:)        = nr_nl_io(1:npft)
ns_nl(:)        = ns_nl_io(1:npft)
r_grow(:)       = r_grow_io(1:npft)
tlow(:)         = tlow_io(1:npft)
tupp(:)         = tupp_io(1:npft)

! Trait physiology parameters
hw_sw(:)        = hw_sw_io(1:npft)
lma(:)          = lma_io(1:npft)
nmass(:)        = nmass_io(1:npft)
nr(:)           = nr_io(1:npft)
nsw(:)          = nsw_io(1:npft)
q10_leaf(:)     = q10_leaf_io(1:npft)
rmass(:)        = rmass_io(1:npft) ! JBaguley
vint(:)         = vint_io(1:npft)
vsl(:)          = vsl_io(1:npft)

! Allometric and other parameters.
a_wl(:)         = a_wl_io(1:npft)
a_ws(:)         = a_ws_io(1:npft)
b_wl(:)         = b_wl_io(1:npft)
eta_sl(:)       = eta_sl_io(1:npft)
sigl(:)         = sigl_io(1:npft)

! Phenology parameters.
dgl_dm(:)       = dgl_dm_io(1:npft)
dgl_dt(:)       = dgl_dt_io(1:npft)
fsmc_of(:)      = fsmc_of_io(1:npft)
g_leaf_0(:)     = g_leaf_0_io(1:npft)
tleaf_of(:)     = tleaf_of_io(1:npft)

! Hydrological, thermal and other "physical" characteristics.
#if !defined(UM_JULES)
calc_rz_psi(:)      = calc_rz_psi_io(1:npft) ! JBaguley
fsmc_mod(:)         = fsmc_mod_io(1:npft)
min_gl_pft(:)       = min_gl_pft_io(1:npft) ! JBaguley
gcut(:)             = gcut_io(1:npft)
! Nonstomatal limitation ramp: keep the pftparm_alloc defaults where unset.
WHERE (ABS(psi_nsl_onset_io(1:npft) - rmdi) > EPSILON(1.0))                  &
  psi_nsl_onset(:) = psi_nsl_onset_io(1:npft)
WHERE (ABS(psi_nsl0_io(1:npft) - rmdi) > EPSILON(1.0))                       &
  psi_nsl0(:) = psi_nsl0_io(1:npft)
WHERE (ABS(fsmc_q_io(1:npft) - rmdi) > EPSILON(1.0)) fsmc_q(:) = fsmc_q_io(1:npft)
IF ( ANY(fsmc_q(:) <= 0.0) ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode, 'fsmc_q_io must be > 0.')
END IF
IF ( l_som_nsl .AND. ( ANY(psi_nsl_onset(:) > 0.0) .OR.                       &
                       ANY(psi_nsl0(:) >= psi_nsl_onset(:)) ) ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                         &
               'l_som_nsl needs psi_nsl0_io < psi_nsl_onset_io <= 0.')
END IF
! Soil-water down-regulation of capacity: keep the defaults where unset.
WHERE (ABS(psi_vcmax_f_io(1:npft) - rmdi) > EPSILON(1.0))                    &
  psi_vcmax_f(:) = psi_vcmax_f_io(1:npft)
WHERE (ABS(sf_vcmax_io(1:npft) - rmdi) > EPSILON(1.0))                       &
  sf_vcmax(:) = sf_vcmax_io(1:npft)
WHERE (ABS(psi_vcmax_fmin_io(1:npft) - rmdi) > EPSILON(1.0))                 &
  psi_vcmax_fmin(:) = psi_vcmax_fmin_io(1:npft)
IF ( l_som_vcmax_psi .AND. ( ANY(psi_vcmax_fmin(:) < 0.0) .OR.               &
                             ANY(psi_vcmax_fmin(:) >= 1.0) ) ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                         &
               'l_som_vcmax_psi needs 0 <= psi_vcmax_fmin_io < 1.')
END IF
WHERE (ABS(nsl_sink_umax_io(1:npft) - rmdi) > EPSILON(1.0))                   &
  nsl_sink_umax(:) = nsl_sink_umax_io(1:npft)
WHERE (ABS(nsl_sink_tau_io(1:npft) - rmdi) > EPSILON(1.0))                    &
  nsl_sink_tau(:) = nsl_sink_tau_io(1:npft)
WHERE (ABS(nsl_sink_maint_io(1:npft) - rmdi) > EPSILON(1.0))                  &
  nsl_sink_maint(:) = nsl_sink_maint_io(1:npft)
WHERE (ABS(nsl_sink_psi50_io(1:npft) - rmdi) > EPSILON(1.0))                  &
  nsl_sink_psi50(:) = nsl_sink_psi50_io(1:npft)
WHERE (ABS(nsl_sink_sf_io(1:npft) - rmdi) > EPSILON(1.0))                     &
  nsl_sink_sf(:) = nsl_sink_sf_io(1:npft)
IF ( l_som_nsl_sink .AND. ( ANY(nsl_sink_umax(:) <= 0.0) .OR.             &
     ANY(nsl_sink_tau(:) <= 0.0) .OR. ANY(nsl_sink_maint(:) < 0.0) .OR.       &
     ANY(nsl_sink_maint(:) > 1.0) .OR. ANY(nsl_sink_sf(:) <= 0.0) ) ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                         &
               'l_som_nsl_sink needs nsl_sink_umax, _tau, _sf > 0 and 0 <= _maint <= 1.')
END IF
IF ( l_som_vcmax_psi .AND. ( ANY(psi_vcmax_f(:) >= 0.0) .OR.                  &
                             ANY(sf_vcmax(:) <= 0.0) ) ) THEN
  errcode = 101
  CALL ereport(RoutineName, errcode,                                         &
               'l_som_vcmax_psi needs psi_vcmax_f_io < 0 and sf_vcmax_io > 0.')
END IF
g1_tuzet(:)         = g1_tuzet_io(1:npft)
sf_tuzet(:)         = sf_tuzet_io(1:npft)
psi_f_tuzet(:)      = psi_f_tuzet_io(1:npft)
! mmol m-2 leaf MPa-1 -> mol m-2 leaf Pa-1, as for kmax_pft.
cap_leaf(:)         = cap_leaf_io(1:npft) * 1.0e-9
cap_stem(:)         = cap_stem_io(1:npft) * 1.0e-9
min_rootc_pft(:)    = min_rootc_pft_io(1:npft) ! JBaguley
psi_close(:)        = psi_close_io(1:npft)
psi_open(:)         = psi_open_io(1:npft)
root_psi_crit(:)    = root_psi_crit_io(1:npft) ! JBaguley
root_radi_pft(:)    = root_radi_pft_io(1:npft) ! JBaguley
rootc_density_pft(:)= rootc_density_pft_io(1:npft) ! JBaguley
! The root supply limit uses soil_to_root_k, which smc_ext computes only for
! fsmc_mod = 2; with fsmc_mod 0/1 it is zero and the stomata shut silently.
IF ( l_som_root_supply ) THEN
  IF ( ANY(fsmc_mod(:) /= 2) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,                                         &
                 'l_som_root_supply needs fsmc_mod_io = 2 for every PFT.')
  END IF
  IF ( ANY(min_rootc_pft(:) <= 0.0) .OR. ANY(root_radi_pft(:) <= 0.0) .OR.    &
       ANY(rootc_density_pft(:) <= 0.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,                                         &
                 'l_som_root_supply needs min_rootc_pft_io, root_radi_pft_io ' &
                 // 'and rootc_density_pft_io > 0.')
  END IF
  ! With the soil psi bound on (l_bound_soil_wp), a hard clamp at psi_close
  ! means no layer is seen drier than psi_close. If psi_close were above
  ! root_psi_crit, root supply (psi - root_psi_crit > 0) would never run
  ! out. Keep psi_close <= root_psi_crit; never tie them.
  IF ( l_bound_soil_wp ) THEN
    IF ( ANY(psi_close(:) > root_psi_crit(:)) ) THEN
      errcode = 101
      CALL ereport(RoutineName, errcode,                                       &
                   'l_som_root_supply with l_bound_soil_wp needs '            //&
                   'psi_close_io <= root_psi_crit_io for every PFT.')
    END IF
  END IF
END IF
#endif
catch0(:)       = catch0_io(1:npft)
dcatch_dlai(:)  = dcatch_dlai_io(1:npft)
dust_veg_scj(:) = dust_veg_scj_io(1:npft)
dz0v_dh(:)      = dz0v_dh_io(1:npft)
emis_pft(:)     = emis_pft_io(1:npft)
fsmc_p0(:)      = fsmc_p0_io(1:npft)
glmin(:)        = glmin_io(1:npft)
gsoil_f(:)      = gsoil_f_io(1:npft)
infil_f(:)      = infil_f_io(1:npft)
irrig_tile(1:npft)      = irrig_pft_io(1:npft)
rootd_ft(:)     = rootd_ft_io(1:npft)
z0v(:)          = z0v_io(1:npft)
z0h_z0m(1:npft) = z0hm_pft_io(1:npft)
z0h_z0m_classic(1:npft) = z0hm_classic_pft_io(1:npft)


! Ozone damage parameters.
dfp_dcuo(:)     = dfp_dcuo_io(1:npft)
fl_o3_ct(:)     = fl_o3_ct_io(1:npft)

! BVOC emission parameters.
aef(:)          = aef_io(1:npft)
ci_st(:)        = ci_st_io(1:npft)
gpp_st(:)       = gpp_st_io(1:npft)
ief(:)          = ief_io(1:npft)
mef(:)          = mef_io(1:npft)
tef(:)          = tef_io(1:npft)

! INFERNO combustion parameters
avg_ba(:)       = avg_ba_io(1:npft)
ccleaf_max(:)   = ccleaf_max_io(1:npft)
ccleaf_min(:)   = ccleaf_min_io(1:npft)
ccwood_max(:)   = ccwood_max_io(1:npft)
ccwood_min(:)   = ccwood_min_io(1:npft)
fire_mort(:)  = fire_mort_io(1:npft)

! INFERNO emission parameters
fef_bc(:)       = fef_bc_io(1:npft)
fef_ch4(:)      = fef_ch4_io(1:npft)
fef_co(:)       = fef_co_io(1:npft)
fef_co2(:)      = fef_co2_io(1:npft)
fef_nox(:)      = fef_nox_io(1:npft)
fef_oc(:)       = fef_oc_io(1:npft)
fef_so2(:)      = fef_so2_io(1:npft)
fef_c2h4(:)     = fef_c2h4_io(1:npft)
fef_c2h6(:)     = fef_c2h6_io(1:npft)
fef_c3h8(:)     = fef_c3h8_io(1:npft)
fef_hcho(:)     = fef_hcho_io(1:npft)
fef_mecho(:)    = fef_mecho_io(1:npft)
fef_nh3(:)      = fef_nh3_io(1:npft)
fef_dms(:)      = fef_dms_io(1:npft)

! SUGAR parameters
sug_g0(:)       = sug_g0_io(1:npft)
sug_grec(:)     = sug_grec_io(1:npft)
sug_yg(:)       = sug_yg_io(1:npft)

! stomatal optimisation model JBaguley
leaf_crit(:)    = leaf_crit_io(1:npft)
pft_conductance_model(:) = pft_conductance_model_io(1:npft)
kcrit_fractional_loss(:) = kcrit_fractional_loss_io(1:npft)
! kmax_pft_io is in mmol m-2 s-1 MPa-1; the model works in
! mol m-2 s-1 Pa-1: x 1e-3 (mmol -> mol) x 1e-6 (MPa-1 -> Pa-1)
kmax_pft(:)     = kmax_pft_io(1:npft) * 1.0e-9
P50(:)          = p50_io(1:npft)
P88(:)          = p88_io(1:npft)
! SOX parameters
sox_a(:)        = sox_a_io(1:npft)
sox_p50(:)      = sox_p50_io(1:npft)
sox_rp_min(:)   = sox_rp_min_io(1:npft)

! ---------------------------------------------------------------------
! The conductance model parameters and critical conductance while not
!  input by the user are constant through out the simulation. To save
!  repeated calculations they are calculated here.
! ---------------------------------------------------------------------
DO i = 1, npft

  ! Calculate critical conductance
  kcrit(i) = (1-kcrit_fractional_loss(i))*kmax_pft(i)

  ! First check that the input values of P50 and P88 are valid if a
  ! conductance model has been chosen (pft_conductance_model = 1 or 2).
  ! The values of P50 and P88 must satisfy the following condition:
  !   0.0 > P50 > P88
  SELECT CASE (pft_conductance_model(i))
  CASE (1,2)
    ! kmax_pft_io is in mmol m-2 s-1 MPa-1 (plant values ~0.1-10). A value
    ! below 1e-4 is almost certainly still in the old mol m-2 s-1 Pa-1
    ! units (e.g. 0.4e-9), which would silently mean ~no conductance.
    IF (kmax_pft_io(i) <= 0.0 .OR. kmax_pft_io(i) < 1.0e-4) THEN
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
               'kmax_pft_io must be > 1e-4, in mmol m-2 s-1 MPa-1 (e.g. ' //   &
               '0.9), not mol m-2 s-1 Pa-1 (e.g. 0.9e-9).')
    END IF
    ! Raise error if P50 is not less than zero.
    IF (0.0 <= P50(i)) THEN
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
               'P50 must be less than zero.')
    ! Raise error if P88 is not less than P50.
    ELSE IF (P50(i) <= P88(i)) THEN
      errcode = 101  !  a hard error
      CALL ereport(RoutineName, errcode,                                       &
               'P88 must be less than P50 (and zero).')
    END IF
  CASE DEFAULT
    ! No conductance model selected so no need to test parameters.
  END SELECT

  SELECT CASE (pft_conductance_model(i))
  CASE (1) ! Cumulative Weibull xylem conductance model
    ! Model equation:
    !   k(psi) = kmax * exp(-(psi/b)^c)

    ! Equation for the shape parameter c:
    !      ln(ln(1-0.5)/ln(1-0.88))
    ! c = --------------------------
    !          (ln(P50/P88))
    ! NOTE: this used to be hardcoded as the literal -0.78135, which does not
    ! equal ln(ln(1-0.5)/ln(1-0.88)) (that evaluates to -1.11805): the
    ! literal was actually ln(ln(1-0.5)/ln(1-0.78)), i.e. it located the
    ! P88 parameter at 78% loss of conductivity (22% conductance remaining)
    ! rather than the intended 88% loss (12% remaining). This flattened the
    ! vulnerability curve, most severely beyond P88 (e.g. ~23x too much
    ! conductance retained at 9 MPa below P50 for a P50/P88 spacing of
    ! ~1 MPa), understating hydraulic cost/risk under severe water stress.
    conductance_c(i) = LOG( LOG(1.0 - 0.5)/LOG(1.0 - 0.88) )/LOG( P50(i)/P88(i) )

    ! Equation for the sensitivity parameter b:
    !             P50                 P50
    ! b = ------------------ = -----------------
    !      (-ln(0.5))^(1/c)     (0.69315)^(1/c)
    conductance_b(i) = P50(i)/( 0.69315**(1/conductance_c(i)) )

  CASE (2) ! SOX xylem conductance model
    ! Model equation:
    !   k(psi) = kmax / (1 + (psi/P50)^c)

    ! By the definition of the conductance model used by SOX, the
    !  sensitivity parameter b is equal to the P50 parameter.
    conductance_b(i) = P50(i)

    ! Equation for the shape parameter c:
    !      ln(1/0.12 - 1)
    ! c = ----------------
    !       ln(P88/P50)
    ! NOTE: this used to be hardcoded as the literal 1.2657 = ln(1/0.22 - 1),
    ! i.e. it located P88 at 78% loss of conductivity (22% remaining) rather
    ! than the intended 88% loss (12% remaining) - see the matching note in
    ! the CW branch above, which has the same error.
    conductance_c(i) = LOG( 1.0/0.12 - 1.0 )/LOG( P88(i)/P50(i) )


  CASE DEFAULT
    ! No conductance model selected so leave the parameters set to 0.0.
  END SELECT
ENDDO

! ---------------------------------------------------------------------
! Root / stem / leaf segments (l_som_plant_segments), as in GDAY
! (gs_opt.c setup_plant). The segments share the whole-plant resistance in
! proportion to seg_frac, so segment s has maximum conductance
! kmax * sum(seg_frac) / seg_frac(s) and the three in series give kmax when
! well watered. Each segment has its own cumulative Weibull curve from its
! P50/P88, defaulting to the PFT's (one set of traits for the whole plant).
! Only used with pft_conductance_model = 1 (checked in jules_vegetation).
! ---------------------------------------------------------------------
DO i = 1, npft
  seg_frac(:) = [ seg_frac_root_io(i), seg_frac_stem_io(i),                   &
                  seg_frac_leaf_io(i) ]
  p50_seg(:)  = [ p50_root_io(i), p50_stem_io(i), p50_leaf_io(i) ]
  p88_seg(:)  = [ p88_root_io(i), p88_stem_io(i), p88_leaf_io(i) ]
  IF ( l_som_plant_segments .AND. pft_conductance_model(i) /= 1 ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,                                         &
                 'l_som_plant_segments needs pft_conductance_model = 1.')
  END IF
  IF ( ANY(seg_frac(:) <= 0.0) ) THEN
    errcode = 101
    CALL ereport(RoutineName, errcode,                                         &
                 'seg_frac_root_io/stem/leaf must all be > 0.')
  END IF
  DO iseg = 1,3
    seg_kfac(i,iseg) = SUM(seg_frac(:)) / seg_frac(iseg)
    ! Root curve from P50 alone (seg_root_vc_io = 1; see its declaration).
    IF ( iseg == 1 .AND. ( seg_root_vc_io(i) == 1 .OR.                         &
                           seg_root_vc_io(i) == 2 ) ) THEN
      IF ( ABS(p50_seg(1) - rmdi) < EPSILON(1.0) .OR.                          &
           ABS(p88_seg(1) - rmdi) >= EPSILON(1.0) ) THEN
        errcode = 101
        CALL ereport(RoutineName, errcode,                                     &
                     'seg_root_vc_io = 1 or 2 needs p50_root_io set and '   // &
                     'p88_root_io unset.')
      END IF
      IF ( seg_root_vc_io(i) == 1 ) THEN
        ! Christoffersen et al. (2016)
        a_root = 2.176 * ( ABS(p50_seg(1)) * 1.0e-6 )**(-0.17)
        p88_seg(1) = p50_seg(1) * ( 1.0 / 0.12 - 1.0 )**( 1.0 / a_root )
      ELSE
        ! XFT root fit, |P88| = 2.61 |P50|^0.70 (MPa)
        p88_seg(1) = -2.61e6 * ( ABS(p50_seg(1)) * 1.0e-6 )**0.70
      END IF
    ELSE IF ( iseg == 1 .AND. seg_root_vc_io(i) /= 0 ) THEN
      errcode = 101
      CALL ereport(RoutineName, errcode, 'seg_root_vc_io should be 0, 1 or 2.')
    END IF
    IF ( ABS(p50_seg(iseg) - rmdi) < EPSILON(1.0) ) p50_seg(iseg) = P50(i)
    IF ( ABS(p88_seg(iseg) - rmdi) < EPSILON(1.0) ) p88_seg(iseg) = P88(i)
    IF ( pft_conductance_model(i) == 1 ) THEN
      IF ( 0.0 <= p50_seg(iseg) .OR. p50_seg(iseg) <= p88_seg(iseg) ) THEN
        errcode = 101
        CALL ereport(RoutineName, errcode,                                     &
                     'segment P50/P88 must satisfy 0 > P50 > P88.')
      END IF
      conductance_c_seg(i,iseg) = LOG( LOG(1.0 - 0.5)/LOG(1.0 - 0.88) )      &
                                  / LOG( p50_seg(iseg)/p88_seg(iseg) )
      conductance_b_seg(i,iseg) = p50_seg(iseg)                               &
                                  / ( 0.69315**(1/conductance_c_seg(i,iseg)) )
    ELSE
      conductance_b_seg(i,iseg) = conductance_b(i)
      conductance_c_seg(i,iseg) = conductance_c(i)
    END IF
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_pftparm_allocated

END MODULE pftparm_io
