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
                                              cp_fm_to_fm
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_type,&
                                              cp_fm_struct_release
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_constants,                 ONLY: proj_mo_ref_scf,&
                                              proj_mo_ref_xas_tdp
   USE input_section_types,             ONLY: &
        section_get_ival, section_get_keyword, section_release, section_type, section_vals_get, &
        section_vals_get_subs_vals, section_vals_type, section_vals_val_get, section_vals_val_set
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
      INTEGER                                            :: i, nrep, n_rep_val, reftype
      TYPE(proj_mo_type), POINTER                        :: proj_mo
      INTEGER, DIMENSION(:), POINTER                     :: tmp_ints
      TYPE(mo_set_type), DIMENSION(:), OPTIONAL, POINTER :: mos


      NULLIFY (rtp_control%proj_mo_list, tmp_ints, proj_mo, logger, &
               input, proj_mo_section, print_key, mos)

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
        CALL section_vals_val_get(proj_mo_section, "TD_MO_INDEX", i_rep_section=i, &
                                   i_vals=tmp_ints)
        ALLOCATE (proj_mo%td_mo_index, SOURCE=tmp_ints(:))
        CALL section_vals_val_get(proj_mo_section, "TD_MO_SPIN", i_rep_section=i, &
                                   i_val=proj_mo%td_mo_spin)
        IF (proj_mo%td_mo_spin > SIZE(mos)) &
           CALL cp_abort(__LOCATION__, &
               "You asked to project the time dependent BETA spin  while the "// &
               "real time DFT run has only one spin defined. "// &
               "Please set TD_MO_SPIN to 1 or use UKS.")
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
            
      INTEGER                                            :: restart_unit, natom, ispin, nspins, &
                                                            nbr_ref_mo, real_mo_index, mo_index, &
                                                            nbr_mo_max, i
      TYPE(mo_set_type), DIMENSION(:), POINTER           :: mo_qs, mo_ref
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(mp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(cp_fm_struct_type), POINTER                   :: mo_ref_fmstruct

      
      NULLIFY (mo_qs, mo_ref, qs_kind_set, particle_set, para_env, dft_control, &
               mo_ref_fmstruct)

      CALL get_qs_env(qs_env, &
                  qs_kind_set=qs_kind_set, &
                  particle_set=particle_set, &
                  dft_control=dft_control, &
                  mos=mo_qs, &
                  para_env=para_env)
      natom = SIZE(particle_set, 1)
      nspins = SIZE(mo_qs)

      ALLOCATE (mo_ref(nspins))
      DO ispin = 1, nspins
         CALL duplicate_mo_set(mo_ref(ispin), mo_qs(ispin))
      END DO

      IF (para_env%is_source()) &
         CALL open_file(file_name=proj_mo%ref_mo_file_name, &
                     file_action="READ", &
                     file_form="UNFORMATTED", &
                     file_status="OLD", &
                     unit_number=restart_unit)
      
      CALL read_mos_restart_low(mo_ref, para_env=para_env, qs_kind_set=qs_kind_set, &
                             particle_set=particle_set, natom=natom, &
                             rst_unit=restart_unit, multiplicity=dft_control%multiplicity)
      
      IF (para_env%is_source()) CALL close_file(unit_number=restart_unit)
      
      IF (proj_mo%ref_mo_spin > SIZE(mo_ref)) &
         CALL cp_abort(__LOCATION__, &
            "You asked as reference spin the BETA one while the "// &
            "reference .wfn file has only one spin. Use a reference .wfn "// &
            "with 2 spins separated or set REF_MO_SPIN to 1")

      ! Store only the mos required
      nbr_mo_max = mo_ref(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%ncol_global
      IF (proj_mo%ref_mo_index(1) == -1) THEN
         DEALLOCATE(proj_mo%ref_mo_index)
         ALLOCATE(proj_mo%ref_mo_index(nbr_mo_max))
         DO i = 1, nbr_mo_max
            proj_mo%ref_mo_index(i) = i
         END DO
      END IF
      nbr_ref_mo = SIZE(proj_mo%ref_mo_index)
      WRITE(*, *) 'proj_mo%ref_mo_index:', proj_mo%ref_mo_index

      IF (nbr_ref_mo > nbr_mo_max) &
         CALL cp_abort(__LOCATION__, &
         "The number of reference mo is larger then the total number of available one in the .wfn file.")

      ! initialize the ref_mo which is saved
      CALL cp_fm_struct_create(mo_ref_fmstruct, &
         context=mo_ref(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%context, &
         nrow_global=mo_ref(proj_mo%ref_mo_spin)%mo_coeff%matrix_struct%nrow_global, &
         ncol_global=nbr_ref_mo)
      CALL cp_fm_create(proj_mo%mo_ref, mo_ref_fmstruct, 'mo_ref')

      ! fill with the reference mo values (column per column)
      DO mo_index = 1, nbr_ref_mo
        real_mo_index = proj_mo%ref_mo_index(mo_index)
        IF (real_mo_index >  nbr_mo_max) &
           CALL cp_abort(__LOCATION__, &
              "One of reference mo index is larger then the total number of available mo in the .wfn file.")

        CALL cp_fm_to_fm(mo_ref(proj_mo%ref_mo_spin)%mo_coeff, proj_mo%mo_ref, &
                                ncol=1, &
                                source_start=real_mo_index, &
                                target_start=mo_index)
      END DO

      ! Clean temporary variables 
      DO ispin = 1, nspins
        CALL deallocate_mo_set(mo_ref(ispin))
      END DO
      CALL cp_fm_struct_release(mo_ref_fmstruct)

   END SUBROUTINE read_reference_mo_from_wfn

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!>        and write the results.
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_and_write_proj_mo(qs_env, mos_new, proj_mo)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(proj_mo_type), POINTER                        :: proj_mo
      
      INTEGER                                            :: nbr_ref_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: popu, phase

      NULLIFY(matrix_s)

      nbr_ref_mo = SIZE(proj_mo%ref_mo_index)
      ALLOCATE(popu(SIZE(nbr_ref_mo))
      ALLOCATE(phase(SIZE(nbr_ref_mo))

      CALL get_qs_env(qs_env, matrix_s_kp=matrix_s)
     

      CALL compute_proj_mo(popu, phase, mo_td, proj_mo, matrix_s)

      
      DEALLOCATE(popu)
      DEALLOCATE(phase)
      

   END SUBROUTINE compute_and_write_proj_mo

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!> \param  ... 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_proj_mo(popu, phase, mos_new, proj_mo, matrix_s)
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE           :: popu, phase
      TYPE(cp_fm_type), DIMENSION(:), POINTER            :: mos_new
      TYPE(proj_mo_type), POINTER                        :: proj_mo
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_s

      INTEGER                                            :: i, nbr_ref_mo
      REAL(KIND=dp)                                      :: real_proj, imag_proj

      nbr_ref_mo = SIZE(proj_mo%ref_mo_index)
      spin_td = proj_mo%td_mo_spin

      DO i = 1, nbr_ref_mo
         ! Real part of the projection:
         real_proj = 0.0_dp
         mos_new(2*spin_td-1)
         ! Imaginary part of the projection
         imag_proj = 0.0_dp
         
         mos_new(2*spin_td) 

      END DO


      WRITE(*, *) 'for compilation'
   END SUBROUTINE compute_proj_mo

! **************************************************************************************************
!> \brief Write the projection of the current MO coefficients on reference ones
!> \param  ... 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE write_proj_mo()
      WRITE(*, *) 'for compilation'
   END SUBROUTINE write_proj_mo

END MODULE rt_projection_mo_utils

