!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Initialize a QM/MM calculation with Force-Mixing
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_create
  USE cp_subsys_methods,               ONLY: cp_subsys_create
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_release
  USE bibliography,                    ONLY: Golze2013,&
                                             Laino2005,&
                                             Laino2006,&
                                             cite_reference
  USE cell_methods,                    ONLY: write_cell
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_subsys_methods,               ONLY: create_small_subsys
  USE cp_subsys_types,                 ONLY: cp_subsys_type, cp_subsys_release
  USE fist_environment,                ONLY: fist_init
  USE fist_environment_types,          ONLY: fist_env_create,&
                                             fist_env_get,&
                                             fist_env_set,fist_env_release,&
                                             fist_environment_type
  USE fist_nonbond_env_types,          ONLY: fist_nonbond_env_type
  USE force_env_methods,               ONLY: force_env_create
  USE force_env_types,                 ONLY: force_env_get,&
                                             force_env_p_type,&
                                             force_env_release,&
                                             force_env_type
  USE global_types,                    ONLY: global_environment_type
  USE header,                          ONLY: qmmm_header
  USE input_constants,                 ONLY: &
       do_fist, do_multipole_section_off, do_multipole_section_on, &
       do_qmmm_center_every_step, do_qmmm_center_never, &
       do_qmmm_center_pbc_aware, do_qmmm_center_setup_only, do_qmmm_none, &
       do_qs, do_qmmm
  USE input_section_types,             ONLY: section_vals_get,&
                                             section_vals_get_subs_vals,&
                                             section_vals_release,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE pw_env_types,                    ONLY: pw_env_type
  USE qmmmx_util,                      ONLY: setup_force_mixing_qmmm_sections,&
                                             update_force_mixing_labels
  USE qmmm_init,                       ONLY: &
       assign_mm_charges_and_radius, move_or_add_atoms, &
       print_image_charge_info, print_qmmm_charges, print_qmmm_links, &
       qmmm_init_gaussian_type, qmmm_init_periodic_potential, &
       qmmm_init_potential, setup_origin_mm_cell, setup_qmmm_links, &
       setup_qmmm_vars_mm, setup_qmmm_vars_qm
  USE qmmm_links_methods,              ONLY: qmmm_link_Imomm_coord
  USE qmmm_pw_grid,                    ONLY: qmmm_pw_grid_init
  USE qmmmx_types,                     ONLY: qmmmx_env_type, qmmmx_env_release
  USE qmmm_create,                     ONLY: qmmm_env_create
  USE qmmm_types,                      ONLY: qmmm_env_type, qmmm_env_get, qmmm_env_release
  USE qmmm_types_low,                  ONLY: &
       add_set_release, add_set_type, &
       qmmm_env_mm_create, qmmm_env_mm_type, &
       qmmm_env_qm_create, qmmm_env_qm_type, qmmm_links_type,&
       qmmm_env_qm_release, qmmm_env_mm_release
  USE qs_environment,                  ONLY: qs_init
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_create,&
                                             qs_environment_type,&
                                             set_qs_env,qs_env_release
  USE termination,                     ONLY: stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_create'

  PUBLIC :: qmmmx_env_create

CONTAINS

! *****************************************************************************
!> \brief ...
!> \param force_env ...
!> \param root_section ...
!> \param para_env ...
!> \param globenv ...
!> \param force_env_section ...
!> \param subsys_section ...
!> \param use_motion_section ...
!> \param dummy_force_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      02.2012 created [noam]
!> \author Noam Bernstein
! *****************************************************************************
  SUBROUTINE  qmmmx_env_create(qmmmx_env, root_section, para_env, globenv,&
       force_env_section, subsys_section, use_motion_section, error)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(section_vals_type), POINTER         :: root_section
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(global_environment_type), POINTER   :: globenv
    TYPE(section_vals_type), POINTER         :: force_env_section, &
                                                subsys_section
    LOGICAL, INTENT(IN)                      :: use_motion_section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_create', &
      routineP = moduleN//':'//routineN

    TYPE(cp_subsys_type), POINTER            :: subsys
    TYPE(section_vals_type), POINTER         :: qmmm_core_section, &
                                                qmmm_extended_section, &
                                                qmmm_section
    TYPE(qmmm_env_type), POINTER             :: dummy_qmmm_env

    NULLIFY(dummy_qmmm_env)

    qmmm_section => section_vals_get_subs_vals(force_env_section,"QMMM",error=error)

    CALL qmmm_env_create(dummy_qmmm_env, root_section, para_env, globenv,&
                         force_env_section, qmmm_section, subsys_section, use_motion_section, &
                         ignore_outside_box = .TRUE., error=error)
    CALL qmmm_env_get(dummy_qmmm_env, subsys=subsys, error=error)

    CALL update_force_mixing_labels(subsys, qmmm_section, error=error)

    ! using CUR_INDICES and CUR_LABELS, create appropriate QM_KIND sections for two QM/MM calculations
    CALL setup_force_mixing_qmmm_sections(subsys, qmmm_section, qmmm_core_section, qmmm_extended_section, error=error)


    ALLOCATE(qmmmx_env)
    CALL qmmm_env_create(qmmmx_env%core, root_section, para_env, globenv,&
         force_env_section, qmmm_core_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE., error=error)

    CALL qmmm_env_create(qmmmx_env%ext, root_section, para_env, globenv,&
         force_env_section, qmmm_extended_section, subsys_section, use_motion_section, &
         ignore_outside_box = .TRUE., error=error)

    CALL section_vals_release(qmmm_core_section, error=error)
    CALL section_vals_release(qmmm_extended_section, error=error)
    CALL qmmm_env_release(dummy_qmmm_env,error=error)

  END SUBROUTINE qmmmx_env_create

END MODULE qmmmx_create
