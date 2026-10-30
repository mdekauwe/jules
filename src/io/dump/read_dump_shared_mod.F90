! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

MODULE read_dump_shared_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE read_dump_shared(file_name, identifiers, nvars, FILE, var_ids)

!Others
USE io_constants, ONLY:                                                        &
  mode_read, max_dim_var, max_sdf_name_len

USE parallel_mod, ONLY:                                                        &
  is_master_task

USE string_utils_mod, ONLY:                                                    &
  to_string

USE file_mod, ONLY:                                                            &
  file_handle, file_open, file_introspect,                                     &
  file_inquire_dim, file_inquire_var

USE logging_mod,      ONLY: log_fatal, log_warn
USE get_dim_info_mod, ONLY: get_dim_info

IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Check that the given file is a dump compatible with the current
!   run, and read the given identifiers from it. This is code shared by both
!   standalone JULES & Rivers.
!   Note that the reading of the JULES dump is done by the master task and the
!   results scattered to other tasks. This means that dumps written with
!   different amounts of tasks should be interchangable.
!
! Code Owner: Please refer to ModuleLeaders.txt
! This file belongs in TECHNICAL
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
CHARACTER(LEN=*) :: file_name  ! The dump file
CHARACTER(LEN=*) :: identifiers(:)  ! The model identifiers for the
                                    ! variables to define
INTEGER :: nvars
  ! The number of variables we are processing

TYPE(file_handle) :: FILE
  ! The opened file

INTEGER :: var_ids(SIZE(identifiers))
  ! The ids of the variables in the dump file

! Local parameters.
LOGICAL, PARAMETER :: l_reading_true = .TRUE.
  ! A value of .TRUE. that is passed to argument l_reading of subroutine
  ! get_dim_info to show that it is being called in connection with reading
  ! (rather than writing) a dump.

! Work variables
INTEGER :: dim_size_file
  ! The size of the dimension currently being
  ! processed in the file

! Used when defining dimensions and variables
INTEGER :: ndims
  ! The number of dimensions the current variable has
CHARACTER(LEN=max_sdf_name_len) :: dim_names(max_dim_var)
  ! The dimension names the current variable should use
INTEGER :: dim_sizes(max_dim_var)
  ! The dimension sizes for the current variable
INTEGER :: dim_ids(max_dim_var)
  ! The dimension ids for the current variable as
  ! calculated from file_inquire_dim
LOGICAL :: is_record_dim
  ! Detects if the current dimension is a record dim
INTEGER :: ndims_file
  ! The number of dimensions the variable has in file
  ! Compared to ndims above for each variable
INTEGER :: dim_ids_file(max_dim_var)
  ! The ids of the dimensions the variable has in the file
  ! Compared to dim_ids above to verify the variable has the
  ! correct dimensions

LOGICAL :: is_record_var
  ! Indicates if a variable uses the record dimension

LOGICAL :: l_read_from_dump
  ! Used to bypass checking the dimensions of ancil variables that are
  ! not to be read from the dump.

INTEGER :: i, j  ! Loop counters

nvars = SIZE(identifiers)

!-----------------------------------------------------------------------------
! In the master task only, we open the file and check that the correct
! dimensions exist and are of a size compatible with this run
!-----------------------------------------------------------------------------
IF ( is_master_task() ) THEN
  !-----------------------------------------------------------------------------
  ! We use the lowest level file API here, as we don't want to impose the input
  ! grid
  !-----------------------------------------------------------------------------
  FILE=file_open(file_name, mode_read)

  ! We want to auto-detect the dimensions and variables in the file
  CALL file_introspect(FILE)

  DO i = 1,nvars

    !-----------------------------------------------------------------------
    ! Get information about the dimensions used by the variable.
    ! The argument l_reading_true shows that we are reading (not writing) a
    ! dump.
    !-----------------------------------------------------------------------
    CALL get_dim_info( l_reading_true, identifiers(i), ndims,  dim_sizes,      &
                       dim_names, l_read_from_dump )

    !-------------------------------------------------------------------------
    ! Check the dimensions exist and have the correct size
    !-------------------------------------------------------------------------
    IF ( l_read_from_dump ) THEN
      DO j = 1,ndims
        ! Retrive information about the dimension from the file we store the id
        ! for use outside this loop
        CALL file_inquire_dim(                                                 &
          FILE, dim_names(j), dim_ids(j), dim_size_file, is_record_dim         &
        )

        ! Check that we found a dimension
        IF ( dim_ids(j) < 0 )                                                  &
          CALL log_fatal("read_dump_shared",                                   &
                         "Could not find expected dimension '" //              &
                         TRIM(dim_names(j)) // "' in dump file")

        ! Check that the dimension is not a record dimension (there shouldn't
        ! be one in dump files).
        IF ( is_record_dim )                                                   &
          CALL log_fatal("read_dump_shared",                                   &
                         "Dimension '" // TRIM(dim_names(j)) // "' is a " //   &
                         "record dimension - should not exist in dump file")

        ! Check that the dimension has the correct size
        IF ( dim_size_file /= dim_sizes(j) )                                   &
          CALL log_fatal("read_dump_shared",                                   &
                         "Dimension '" // TRIM(dim_names(j)) // "' has " //    &
                         "size incompatible with current run (required: " //   &
                         TRIM(to_string(dim_sizes(j))) // ", found: " //       &
                         TRIM(to_string(dim_size_file)) // ")")
      END DO  ! dims

      !-----------------------------------------------------------------------
      ! Check that the variable exists and has the correct dimensions
      !-----------------------------------------------------------------------
      ! Retrieve information about the variable from the file
      ! (missing_ok: the variables that may be missing are handled below.)
      ! MDK, 18/10/26
      CALL file_inquire_var(                                                   &
        FILE, identifiers(i), var_ids(i), ndims_file, dim_ids_file,            &
        is_record_var, missing_ok = .TRUE.                                     &
      )

      ! Check that we found a variable
      IF ( var_ids(i) < 1 ) THEN
        SELECT CASE ( identifiers(i) )
        CASE ( 'rivers_outflow_rp' )
          CALL log_warn("read_dump_shared",                                    &
             "rivers_outflow_rp is not in the dump so will be initialised " // &
             "to zero.")
        CASE ( 'smcl' )
          ! Older dumps: smcl is then rebuilt from sthuf as before.
          ! MDK, 17/10/26
          CALL log_warn("read_dump_shared",                                    &
             "smcl is not in the dump so will be calculated from sthuf " //    &
             "(restart exact only to rounding).")
        CASE ( 'phen_asteps' )
          ! Dumps written without it: phenology counts from the restart,
          ! as before. MDK, 24/10/26
          CALL log_warn("read_dump_shared",                                    &
             "phen_asteps is not in the dump: phenology steps count " //       &
             "from the restart (restart not exact).")
        CASE ( 'trif_asteps', 'npp_acc', 'resp_w_acc', 'g_leaf_acc',         &
               'g_leaf_phen_acc', 'gpp_acc', 'resp_s_acc' )
          ! TRIFFID's timer and accumulators: dumps written without them. The
          ! TRIFFID steps then count from the restart, as before.
          ! MDK, 22/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump: TRIFFID steps " //  &
             "count from the restart (restart not exact).")
        CASE ( 'fvcmax_psi_state', 'cpl_gl_sun', 'cpl_gl_shd',                &
               'cpl_gbw_sun', 'cpl_gbw_shd', 'cpl_ra_ca',                     &
               't_leaf_sun_pft', 't_leaf_shd_pft' )
          ! (t_leaf_*: the leaf temperatures' first guess, l_leaf_temp; 0 =
          ! not yet set, so the first step starts from tstar. MDK, 24/10/26)
          ! sf_stom state: dumps written without it start from 'not yet set'
          ! (the restart is then not exact). MDK, 22/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts " //       &
             "'not yet set' (restart not exact).")
        CASE ( 'nsc_store', 'nsc_npp_mean', 'nsc_npp_wt', 'nsc_npp_day',      &
               'nsc_turg_acc', 'nsc_t_run', 'nsc_ft_mean',                     &
               'ximpair_repair_pend',                                          &
               'phen_psi_acc', 'phen_nlow', 'phen_nhigh', 'phen_fbleaf',       &
               'nsc_regrow_acc', 'hydmort_pref', 'hydmort_dead',               &
               'hydmort_plc_peak' )
          ! TRIFFID store state: dumps written without the store.
          ! Hydraulic mortality state: dumps written without it (the stand
          ! is taken to have withstood its current damage). MDK, 21/10/26
          ! MDK, 18/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts from " //  &
             "its 'not yet set' value (store filled at the first step).")
        CASE ( 'nsc_wood_lost', 'nsc_net_wood_growth' )
          ! Ring renewal: dumps written without them start with no wood
          ! lost and no net wood growth in the ring year (a restart is
          ! then not exact within a ring year).
          !
          ! MDK, 30/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts at 0 " //  &
             "(restart not exact within a ring year).")
        CASE ( 'nsc_sapwood' )
          ! Sapwood pool: dumps written without it start from the pipe
          ! model at the first ring step.
          !
          ! MDK, 30/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts from " //  &
             "the pipe model at the first ring step.")
        CASE ( 'ximpair_psi_stem', 'ximpair_area_w' )
          ! Dumps written without them: the stem psi of the PLC output is
          ! unset until the first update, and the area dilution of the
          ! impaired kmax since the last update is lost (restart not exact
          ! with l_veg_compete). MDK, 25/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts " //       &
             "'not yet set' (restart not exact).")
        CASE ( 'cut_w_store', 'cut_psi_stem', 'cut_c_store' )
          ! Cuticular-loss store: dumps written without it start it full
          ! (the restart is then not exact).
          !
          ! MDK, 22/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so the store " //    &
             "starts full (restart not exact).")
        CASE ( 'phen_gdd', 'phen_bb_day', 'phen_dormant' )
          ! Budburst state: dumps written without the budburst gate.
          ! MDK, 20/10/26
          CALL log_warn("read_dump_shared",                                    &
             TRIM(identifiers(i)) // " is not in the dump so starts from " //  &
             "its 'not yet set' value (no budburst yet this year).")
        CASE DEFAULT
          CALL log_fatal("read_dump_shared",                                   &
                         "Failed to find requested variable '" //              &
                         TRIM(identifiers(i)) // "' in dump file")
        END SELECT
      ELSE
          ! Check that the number of dimensions match
        IF ( ndims_file /= ndims )                                             &
           CALL log_fatal("read_dump_shared",                                  &
                          "Variable '" // TRIM(identifiers(i)) // "' has " //  &
                          "incorrect number of dimensions in dump file (" //   &
                          "expected: " // TRIM(to_string(ndims)) // ", " //    &
                          "found: " // TRIM(to_string(ndims_file)) // ")")

        ! Check that the dimension ids match
        IF ( .NOT. ALL(dim_ids(1:ndims) == dim_ids_file(1:ndims)) )            &
           CALL log_fatal("read_dump_shared",                                  &
                          "Variable '" // TRIM(identifiers(i)) // "' has " //  &
                          "incorrect dimensions in dump file")
      END IF ! ( var_ids(i) < 1 )
    END IF  !  l_read_from_dump

  END DO  ! vars

END IF  ! MASTER TASK

RETURN

END SUBROUTINE read_dump_shared
END MODULE read_dump_shared_mod
