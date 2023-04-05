!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2023 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the field section of the input
!> \par History
!>      02.2017 moved out of input_cp2k_dft [JHU]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_projection_rtp

   USE input_constants,                 ONLY: proj_mo_ref_scf,&
                                              proj_mo_ref_xas_tdp


   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_projection_rtp'

   PUBLIC :: create_projection_rtp_section

! **************************************************************************************************
!> \brief creates the section for time dependent nonperiodic fields
!> \param section ...
!> \author Florian Schiffmann
! **************************************************************************************************
   SUBROUTINE create_projection_rtp_section(section)
      TYPE(section_type), POINTER                        :: section

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, __LOCATION__, name="PROJECTION_MO" // &
                          description="Projects the time dependent MO " // &
                       "coefficient on reference one. You can define " //  &
                       "several section like this to project the TD-MOs " // &
                       "on different reference MO.", &
                        repeats=.TRUE.)

      NULLIFY (keyword, subsection)

      CALL keyword_create(keyword, __LOCATION__, name="REFERENCE_TYPE", &
                           description="Type of the reference MO. ", &
                           enum_c_vals=s2a("SCF", "XAS_TDP"), &
                           usage="ENVELOP SCF", &
                           default_i_val=proj_mo_ref_scf, &
                           enum_desc=s2a("The reference MO is from an SCF calculation. " // &
                        "Defining MO_REF_NUMBER is recommanded to define the reference MO " // &
                        "among the  available one. ", &
                        "The reference MO is from an XAS_TDP analysis. Only one MO is " // & 
                        "available and thus will be used as reference.", &
                           enum_i_vals=(/proj_mo_ref_scf, proj_mo_ref_xas_tdp/))

      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)
      
      


   END SUBROUTINE create_projection_rtp_section


END MODULE input_cp2k_projection_rtp
