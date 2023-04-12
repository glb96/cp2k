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
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_print_key_generate_filename
   USE input_constants,                 ONLY: proj_mo_ref_scf,&
                                              proj_mo_ref_xas_tdp
   USE input_section_types,             ONLY: &
        section_get_ival, section_get_keyword, section_release, section_type, section_vals_get, &
        section_vals_get_subs_vals, section_vals_type, section_vals_val_get, section_vals_val_set
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE rt_propagation_types,            ONLY: rt_prop_type

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
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE init_mo_projection(qs_env, rtp_control)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(rtp_control_type), POINTER                    :: rtp_control
      
      TYPE(section_vals_type), POINTER                   :: input, proj_mo_section, print_key
      TYPE(cp_logger_type), POINTER                      :: logger
      INTEGER                                            :: i, nrep, n_rep_val, reftype
      TYPE(proj_mo_type), POINTER                        :: proj_mo
      INTEGER, DIMENSION(:), POINTER                     :: tmp_ints

      NULLIFY (rtp_control%proj_mo_list, tmp_ints)

      CALL get_qs_env(qs_env, input=input)
      
      proj_mo_section => section_vals_get_subs_vals(input, "DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO")

      ! Read the input section and load the reference MOs
      CALL section_vals_get(proj_mo_section, n_repetition=nrep)

      ALLOCATE (rtp_control%proj_mo_list(nrep))
      DO i = 1, nrep
        NULLIFY (rtp_control%proj_mo_list(i)%proj_mo)
        ALLOCATE (rtp_control%proj_mo_list(i)%proj_mo)
        proj_mo => rtp_control%proj_mo_list(i)%proj_mo
        NULLIFY (proj_mo%mo_ref)
        
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
          ! Read SCF mos
        

          CALL section_vals_val_get(proj_mo_section, "REF_MO_INDEX", i_rep_section=i, &
                                      i_vals=tmp_ints)
          ALLOCATE (proj_mo%ref_mo_index, SOURCE=tmp_ints(:))
          CALL section_vals_val_get(proj_mo_section, "REF_MO_SPIN", i_rep_section=i, &
                                      i_val=proj_mo%ref_mo_spin)

        ELSE IF (reftype == proj_mo_ref_xas_tdp) THEN
          IF (proj_mo%ref_mo_file_name == "DEFAULT") THEN 
             CALL cp_abort(__LOCATION__, &
                "Input error in DFT%REAL_TIME_PROPAGATION%PRINT%PROJECTION_MO. "// &
                "For REFERENCE_TYPE XAS_TDP one must define the name "// &
                "of the .wfn file to read the reference MO from. Please define REF_MO_FILE_NAME.")
          END IF
          ! Read XAS TDP mos
        END IF
        
        WRITE(*, *) 'default value for REF_MO_FILE_NAME = ', TRIM(proj_mo%ref_mo_file_name), ' ok?' 
         
        

        CALL section_vals_val_get(proj_mo_section, "SUM_ON_ALL_REF", i_rep_section=i, &
                                   l_val=proj_mo%sum_on_all_ref)
        CALL section_vals_val_get(proj_mo_section, "TD_MO_INDEX", i_rep_section=i, &
                                   i_vals=tmp_ints)
        ALLOCATE (proj_mo%td_mo_index, SOURCE=tmp_ints(:))
        CALL section_vals_val_get(proj_mo_section, "TD_MO_SPIN", i_rep_section=i, &
                                   i_val=proj_mo%td_mo_spin)
      END DO 


   END SUBROUTINE init_mo_projection

! **************************************************************************************************
!> \brief Read the MO from .wfn file
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE read_reference_mo_from_wfn()
      WRITE(*, *) 'for compilation'
   END SUBROUTINE read_reference_mo_from_wfn

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!>        and write the results.
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_and_write_proj_mo()
      WRITE(*, *) 'for compilation'
   END SUBROUTINE compute_and_write_proj_mo

! **************************************************************************************************
!> \brief Compute the projection of the current MO coefficients on reference ones
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE compute_proj_mo()
      WRITE(*, *) 'for compilation'
   END SUBROUTINE compute_proj_mo

! **************************************************************************************************
!> \brief Write the projection of the current MO coefficients on reference ones
!> \param  ... UNMATCHED_PROCEDURE_ARGUMENT: please check 
!> \author Guillaume Le Breton
! **************************************************************************************************
   SUBROUTINE write_proj_mo()
      WRITE(*, *) 'for compilation'
   END SUBROUTINE write_proj_mo

END MODULE rt_projection_mo_utils

