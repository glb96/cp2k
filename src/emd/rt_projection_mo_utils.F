!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Function related to MO projection in RTP calculations
!> \author Guillaume Le Breton 04.2023
! **************************************************************************************************
MODULE rt_projection_mo_utils
   USE cp_control_types,                ONLY: dft_control_type,&
                                              rtp_control_type,&
                                              proj_mo_type
   USE cp_files,                        ONLY: close_file,&
                                              open_file                                           
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_release
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_type,&
                                              cp_fm_struct_release
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_iter_string,&
                                              cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_generate_filename,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_dbcsr_operations,             ONLY: cp_dbcsr_sm_fm_multiply
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_constants,                 ONLY: proj_mo_ref_scf,&
                                              proj_mo_ref_xas_tdp
   USE input_section_types,             ONLY: &
        section_get_ival, section_get_keyword, section_release, section_type, section_vals_get, &
        section_vals_get_subs_vals, section_vals_type, section_vals_val_get, section_vals_val_set
   USE kinds,                           ONLY: default_string_length,&
                                              dp
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: deallocate_mo_set,&
                                              duplicate_mo_set,&
                                              get_mo_set,&
                                              mo_set_type
   USE rt_propagation_types,            ONLY: rt_prop_type
   USE particle_types,                  ONLY: particle_type
   USE message_passing,                 ONLY: mp_para_env_type                                           
   USE qs_mo_io,                        ONLY: read_mos_restart_low

#include "./../base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_projection_rtp'

   PUBLIC :: init_mo_projection, compute_and_write_proj_mo

CONTAINS

! **************************************************************************************************
!> \brief Initialize the mo projection objects for time dependent run
!> \param qs_env ...
!> \param rtp_control ...
!> \author Guillaume Le Breton (04.2023)
! **************************************************************************************************
   SUBROUTINE init_mo_projection(qs_env, rtp_control)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      
      TYPE(section_vals_type), POINTER                   :: input, proj_mo_section, print_key
      TYPE(cp_logger_type), POINTER                      :: logger 
      INTEGER                                            :: i, j, nrep, n_rep_val, reftype, nbr_mo_td_max
      TYPE(proj_mo_type), POINTER                        :: proj_mo
      INTEGER, DIMENSION(:), POINTER                     :: tmp_ints
      TYPE(mo_set_type), DIMENSION(:),  POINTER          :: mos


      NULLIFY (rtp_control%proj_mo_list, tmp_ints, proj_mo, logger, &
               input, proj_mo_section, print_key, mos)
      
      IF (.NOT. rtp_control%fixed_ions) & 
           CALL cp_abort(__LOCATION__, &
                "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO not implemented for EMD yet.")

      CALL get_qs_env(qs_env, input=input, mos=mos)
      
      proj_mo_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO")

      ! Read the input section and load the reference MOs
      CALL section_vals_get(proj_mo_section, n_repetition=nrep)
      ALLOCATE (rtp_control%proj_mo_list(nrep))

      DO i = 1, nrep
        NULLIFY (rtp_control%proj_mo_list(i)%proj_mo)
        ALLOCATE (rtp_control%proj_mo_list(i)%proj_mo)
        proj_mo => rtp_control%proj_mo_list(i)%proj_mo
        
        CALL section_vals_val_get(proj_mo_section, "REFERENCE_TYPE", i_rep_section=i, &
                                   i_val=reftype)

        CALL section_vals_val_get(proj_mo_section, "REF_MO_FILE_NAME", i_rep_section=i, &
                                   c_val=proj_mo%ref_mo_file_name)
         
        IF (reftype == proj_mo_ref_scf) THEN
          ! If no reference .wfn is provided, using the restart SCF file:
          IF (proj_mo%ref_mo_file_name == "DEFAULT") THEN
             CALL section_vals_val_get(input, "DFT%WFN_RESTART_FILE_NAME", n_rep_val=n_rep_val)
             IF (n_rep_val > 0) THEN
               CALL section_vals_val_get(input, "DFT%WFN_RESTART_FILE_NAME", c_val=proj_mo%ref_mo_file_name)
             ELSE
                 !try to read from the filename that is generated automatically from the printkey
                  print_key => section_vals_get_subs_vals(input, "DFT%SCF%PRINT%RESTART")
                  logger => cp_get_default_logger()
                  proj_mo%ref_mo_file_name = cp_print_key_generate_filename(logger, print_key, &
                                                      extension=".wfn", my_local=.FALSE.)
             END IF
          END IF
          
          CALL section_vals_val_get(proj_mo_section, "REF_MO_INDEX", i_rep_section=i, &
                                      i_vals=tmp_ints)
          ALLOCATE (proj_mo%ref_mo_index, SOURCE=tmp_ints(:))
          CALL section_vals_val_get(proj_mo_section, "REF_MO_SPIN", i_rep_section=i, &
                                      i_val=proj_mo%ref_mo_spin)
          
          ! Read the SCF mos and store the one required
          CALL read_reference_mo_from_wfn(qs_env, proj_mo) 

        ELSE IF (reftype == proj_mo_ref_xas_tdp) THEN
          IF (proj_mo%ref_mo_file_name == "DEFAULT") THEN 
             CALL cp_abort(__LOCATION__, &
                "Input error in DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO. "// &
                "For REFERENCE_TYPE XAS_TDP one must define the name "// &
                "of the .wfn file to read the reference MO from. Please define REF_MO_FILE_NAME.")
          END IF
          ! Read XAS TDP mos
        END IF
        
        WRITE(*, *) 'default value for REF_MO_FILE_NAME = ', TRIM(proj_mo%ref_mo_file_name)
         
        
        ! Initialize the other parameters related to the TD mos.
        CALL section_vals_val_get(proj_mo_section, "SUM_ON_ALL_REF", i_rep_section=i, &
                                   l_val=proj_mo%sum_on_all_ref)

        CALL section_vals_val_get(proj_mo_section, "TD_MO_SPIN", i_rep_section=i, &
                                   i_val=proj_mo%td_mo_spin)
        IF (proj_mo%td_mo_spin > SIZE(mos)) &
           CALL cp_abort(__LOCATION__, &
               "You asked to project the time dependent BETA spin  while the "// &
               "real time DFT run has only one spin defined. "// &
               "Please set TD_MO_SPIN to 1 or use UKS.")

        CALL section_vals_val_get(proj_mo_section, "TD_MO_INDEX", i_rep_section=i, &
                                   i_vals=tmp_ints)

        ALLOCATE (proj_mo%td_mo_index, SOURCE=tmp_ints(:))
        nbr_mo_td_max = mos(proj_mo%td_mo_spin)%mo_coeff%matrix_struct%ncol_global
        IF (proj_mo%td_mo_index(1) == -1) THEN
          DEALLOCATE(proj_mo%td_mo_index)
          ALLOCATE(proj_mo%td_mo_index(nbr_mo_td_max))
          DO j = 1, nbr_mo_td_max
            proj_mo%td_mo_index(i) = i
          END DO
        ELSE
          DO j = 1, SIZE(proj_mo%td_mo_index)
            IF (proj_mo%td_mo_index(i) > nbr_mo_td_max) &
              CALL cp_abort(__LOCATION__, &
                     "The MO number available in the Time Dependent run "// &
                     "is smaller than the MO number you have required in TD_MO_INDEX.")
          END DO
        END IF

      END DO 

   END SUBROUTINE init_mo_projection

! **************************************************************************************************
!> \brief Read the MO from .wfn file and store the required mos for TD projections
!> \param  qs_env...
!> \param  proj_mo...
!> \author Guillaume Le Breton (04.2023)
! **************************************************************************************************
   SUBROUTINE read_reference_mo_from_wfn(qs_env, proj_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(proj_mo_type), POINTER                        :: proj_mo
            
      LOGICAL                                            :: is_file
      INTEGER                                            :: restart_unit, natom, ispin, nspins, &
                                                            nbr_ref_mo, real_mo_index, mo_index, &
                                                            nbr_mo_max, i
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_qs, mo_ref_temp
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_fm_struct_type), POINTER                   :: mo_ref_fmstruct
      TYPE(cp_fm_type)                                   :: mo_coeff_temp

      
      NULLIFY (mo_qs, mo_ref_temp, qs_kind_set, particle_set, para_env, dft_control, &
               mo_ref_fmstruct, matrix_s)

      CALL get_qs_env(qs_env, &
                  qs_kind_set=qs_kind_set, &
                  particle_set=particle_set, &
                  dft_control=dft_control, &
                  matrix_s_kp=matrix_s, &
                  mos=mo_qs, &
                  para_env=para_env)

      natom = SIZE(particle_set, 1)
      nspins = SIZE(mo_qs)

      ALLOCATE (mo_ref_temp(nspins))
      DO ispin = 1, nspins
         CALL duplicate_mo_set(mo_ref_temp(ispin), mo_qs(ispin))
      END DO

      IF (para_env%is_source()) THEN 
         INQUIRE (FILE=TRIM(proj_mo%ref_mo_file_name), exist=is_file)
         IF (.NOT. is_file) & 
            CALL cp_abort(__LOCATION__, &
            "Reference file not found! Name of the file CP2K looked for: " // TRIM(proj_mo%ref_mo_file_name))

         CALL open_file(file_name=proj_mo%ref_mo_file_name, &
                     file_action="READ", &
                     file_form="UNFORMATTED", &
                     file_status="OLD", &
                     unit_number=restart_unit)
      END IF

      CALL read_mos_restart_low(mo_ref_temp, para_env=para_env, qs_kind_set=qs_kind_set, &
                             particle_set=particle_set, natom=natom, &
                             rst_unit=restart_unit, multiplicity=dft_control%multiplicity)
      
      IF (para_env%is_source()) CALL close_file(unit_number=restart_unit)
      
      IF (proj_mo%ref_mo_spin > SIZE(mo_ref_temp)) &
         CALL cp_abort(__LOCATION__, &
            "You asked as reference spin the BETA one while the "// &
            "reference .wfn file has only one spin. Use a reference .wfn "// &
            "with 2 spins separated or set REF_MO_SPIN to 1")

      ! Store only the mos required
      nbr_mo_max = mo_ref_temp(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%ncol_global
      IF (proj_mo%ref_mo_index(1) == -1) THEN
         DEALLOCATE(proj_mo%ref_mo_index)
         ALLOCATE(proj_mo%ref_mo_index(nbr_mo_max))
         DO i = 1, nbr_mo_max
            proj_mo%ref_mo_index(i) = i
         END DO
      ELSE
         DO i = 1, SIZE(proj_mo%ref_mo_index)
              IF (proj_mo%ref_mo_index(i) > nbr_mo_max) &
                 CALL cp_abort(__LOCATION__, &
                     "The MO number available in the Reference SCF "// &
                     "is smaller than the MO number you have required in REF_MO_INDEX.")
         END DO
      END IF
      nbr_ref_mo = SIZE(proj_mo%ref_mo_index)
      WRITE(*, *) 'proj_mo%ref_mo_index:', proj_mo%ref_mo_index

      IF (nbr_ref_mo > nbr_mo_max) &
         CALL cp_abort(__LOCATION__, &
         "The number of reference mo is larger then the total number of available one in the .wfn file.")
      
      ! Store
      ALLOCATE(proj_mo%mo_ref(nbr_ref_mo))
      CALL cp_fm_struct_create(mo_ref_fmstruct, &
         context=mo_ref_temp(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%context, &
         nrow_global=mo_ref_temp(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%nrow_global, &
         ncol_global=1)
      CALL cp_fm_create(mo_coeff_temp, mo_ref_fmstruct, 'mo_ref')

      DO mo_index = 1, nbr_ref_mo
        real_mo_index = proj_mo%ref_mo_index(mo_index)
        IF (real_mo_index >  nbr_mo_max) &
           CALL cp_abort(__LOCATION__, &
              "One of reference mo index is larger then the total number of available mo in the .wfn file.")

        ! fill with the reference mo values
        CALL cp_fm_create(proj_mo%mo_ref(mo_index), mo_ref_fmstruct, 'mo_ref')
        CALL cp_fm_to_fm(mo_ref_temp(proj_mo%ref_mo_spin)%mo_coeff, mo_coeff_temp, &
                                ncol=1, &
                                source_start=real_mo_index, &
                                target_start=1)

        ! multiply with overlap matrix to save time later on
        CALL cp_dbcsr_sm_fm_multiply(matrix_s(1, 1)%matrix, mo_coeff_temp,  proj_mo%mo_ref(mo_index), ncol=1)
        !CALL cp_fm_to_fm(mo_coeff_temp, proj_mo%mo_ref(mo_index))
      END DO

      ! Clean temporary variables 
      DO ispin = 1, nspins
        CALL deallocate_mo_set(mo_ref_temp(ispin))
      END DO
      DEALLOCATE(mo_ref_temp)

      CALL cp_fm_struct_release(mo_ref_fmstruct)
      CALL cp_fm_release(mo_coeff_temp)


   END SUBROUTINE read_reference_mo_from_wfn

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!>        and write the results.
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_and_write_proj_mo(qs_env, mos_new, proj_mo, n_proj)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(proj_mo_type)                                 :: proj_mo
      INTEGER                                            :: n_proj
      
      INTEGER                                            :: i, nbr_ref_mo, nbr_ref_td
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: popu, phase
      REAL(KIND=dp)                                      :: popu_tot
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: input, proj_mo_section
      TYPE(cp_logger_type), POINTER                      :: logger

      NULLIFY( dft_control, input, proj_mo_section, logger)
      
      logger => cp_get_default_logger()

      CALL get_qs_env(qs_env, &
               dft_control=dft_control, &
               input=input)
    
      proj_mo_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO", i_rep_section=n_proj)

      ! STOP if not the required time step
      IF (.NOT. BTEST(cp_print_key_should_output(logger%iter_info, &
                                             proj_mo_section, "PRINT"), &
                                             cp_p_file)) &
         RETURN

      IF (.NOT. dft_control%rtp_control%fixed_ions) &
           CALL cp_abort(__LOCATION__, &
                "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO not implemented for EMD yet.")

      nbr_ref_mo = SIZE(proj_mo%ref_mo_index)
      nbr_ref_td = SIZE(proj_mo%td_mo_index)

      ALLOCATE(popu(nbr_ref_td))
      ALLOCATE(phase(nbr_ref_td))
      DO i = 1, nbr_ref_mo
         CALL compute_proj_mo(popu, phase, mos_new, proj_mo, i)

         IF (proj_mo%sum_on_all_ref) THEN
            popu_tot = SUM(popu)
            CALL write_proj_mo(qs_env, proj_mo_section, proj_mo, i, popu_tot=popu_tot)
         ELSE

            CALL write_proj_mo(qs_env, proj_mo_section, proj_mo, i, popu=popu, phase=phase)
         END IF
      END DO
      
      DEALLOCATE(popu)
      DEALLOCATE(phase)
      

   END SUBROUTINE compute_and_write_proj_mo

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!> \param  ... 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_proj_mo(popu, phase, mos_new, proj_mo, i)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: popu, phase
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(proj_mo_type)                                 :: proj_mo
      INTEGER                                            :: i, handle

      CHARACTER(len=*), PARAMETER                        :: routineN = 'rt_compute_proj_mo'
      INTEGER                                            :: j, nbr_ref_td, spin_td
      REAL(KIND=dp)                                      :: real_proj, imag_proj
      TYPE(cp_fm_struct_type), POINTER                   :: mo_ref_fmstruct
      TYPE(cp_fm_type)                                   :: mo_coeff_temp


      CALL timeset(routineN, handle)

      nbr_ref_td = SIZE(popu)
      spin_td = proj_mo%td_mo_spin

      CALL cp_fm_struct_create(mo_ref_fmstruct, &
         context=mos_new(1)%matrix_struct%context, &
         nrow_global=mos_new(1)%matrix_struct%nrow_global, &
         ncol_global=1)
      CALL cp_fm_create(mo_coeff_temp, mo_ref_fmstruct, 'mo_temp')


      DO j = 1, nbr_ref_td
         ! Real part of the projection:
         real_proj = 0.0_dp
         CALL cp_fm_to_fm(mos_new(2*spin_td-1), mo_coeff_temp, &
                                ncol=1, &
                                source_start=proj_mo%td_mo_index(j), &
                                target_start=1)
         
         CALL cp_fm_trace(mo_coeff_temp, proj_mo%mo_ref(i), real_proj) 

         ! Imaginary part of the projection
         imag_proj = 0.0_dp
         CALL cp_fm_to_fm(mos_new(2*spin_td), mo_coeff_temp, &
                                ncol=1, &
                                source_start=proj_mo%td_mo_index(j), &
                                target_start=1)

         CALL cp_fm_trace(mo_coeff_temp, proj_mo%mo_ref(i), imag_proj)
         
         ! Store the result
         phase(j) = ATAN2(imag_proj, real_proj) ! in radians
         popu(j) = real_proj**2 + imag_proj**2

      END DO
      
      CALL cp_fm_struct_release(mo_ref_fmstruct)
      CALL cp_fm_release(mo_coeff_temp)
      
      CALL timestop(handle)

   END SUBROUTINE compute_proj_mo

! **************************************************************************************************
!> \brief Write the projection of the current MO coefficients on reference ones
!> \param  ... 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE write_proj_mo(qs_env, proj_mo_section, proj_mo, i, popu, phase, popu_tot)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(section_vals_type), POINTER                   :: proj_mo_section
      TYPE(proj_mo_type)                                 :: proj_mo
      INTEGER                                            :: i 

      REAL(KIND=dp), DIMENSION(:), OPTIONAL              :: popu, phase
      REAL(KIND=dp), OPTIONAL                            :: popu_tot

      CHARACTER(LEN=default_string_length)               :: filename, ext
      INTEGER                                            :: output_unit, print_unit, j
      TYPE(cp_logger_type), POINTER                      :: logger


      NULLIFY( logger)

      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      IF (output_unit > 0) THEN
         ! Filename is update wrt the reference MO number
         ext = "-REF-"// &
               TRIM(ACHAR(IACHAR('0') + proj_mo%ref_mo_index(i)))// &
               ".dat"
         print_unit = cp_print_key_unit_nr(logger, proj_mo_section, "PRINT", &
                                           extension=TRIM(ext))

         IF (print_unit /= output_unit) THEN
            INQUIRE (UNIT=print_unit, NAME=filename)
            WRITE (UNIT=output_unit, FMT="(/,T2,A,2(/,T3,A),/)") &
               "PROJECTION MO", "The projection of the TD MOs is done in the file:", &
               TRIM(filename)
            WRITE (UNIT=print_unit, FMT="(/,(T2,A,T40,I6))") &
               "Real time propagation step:", qs_env%sim_step
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,T2,A)") "PROJECTION MO"
         END IF

         WRITE (print_unit, "(T3,A)") & 
            "Projection on the MO number "// &
            TRIM(ACHAR(IACHAR('0') + proj_mo%ref_mo_index(i)))// &
            " from the reference file "// &
            TRIM(proj_mo%ref_mo_file_name)

         IF (proj_mo%sum_on_all_ref) THEN
            WRITE (print_unit, "(T3, A, E16.8)") & 
               "The sum over all the TD MOs population:", popu_tot
         ELSE
            WRITE (print_unit, "(T3,A)") &
               "For each TD MOs required is printed: Population & Phase [rad] "
            DO j = 1, SIZE(popu)
               WRITE (print_unit, "(T5,2(E16.8, E16.8, 1X))")  popu(j), phase(j)
            END DO 
         END IF

         CALL cp_print_key_finished_output(print_unit, logger, proj_mo_section, "PRINT")

      END IF 

   END SUBROUTINE write_proj_mo

END MODULE rt_projection_mo_utils

