!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Calculates QM/MM energy and forces with Force-Mixing
!> \par History
!>      2015 Factored out of force_env_methods.F
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_force
  USE atprop_types,                    ONLY: atprop_init,&
                                             atprop_type
  USE cell_types,                      ONLY: cell_clone,&
                                             cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             init_cell,&
                                             real_to_scaled,&
                                             scaled_to_real
  USE constraint_fxd,                  ONLY: fix_atom_control
  USE constraint_vsite,                ONLY: vsite_force_control
  USE cp_iter_types,                   ONLY: cp_iteration_info_copy_iter
  USE cp_output_handling,              ONLY: cp_print_key_finished_output,&
                                             cp_print_key_unit_nr,&
                                             low_print_level
  USE cp_para_env,                     ONLY: cp_para_env_retain
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             cp_results_mp_bcast,&
                                             get_results,&
                                             test_for_result
  USE cp_result_types,                 ONLY: cp_result_copy,&
                                             cp_result_create,&
                                             cp_result_p_type,&
                                             cp_result_release,&
                                             cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_p_type,&
                                             cp_subsys_set,&
                                             cp_subsys_type
  USE eip_environment_types,           ONLY: eip_env_retain,&
                                             eip_environment_type
  USE eip_silicon,                     ONLY: eip_bazant,&
                                             eip_lenosky
  USE external_potential_methods,      ONLY: add_external_potential
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_env_retain,&
                                             fist_environment_type
  USE fist_force,                      ONLY: fist_calc_energy_force
  USE force_env_types,                 ONLY: &
       force_env_get, force_env_get_natom, force_env_p_type, force_env_set, &
       force_env_type, use_eip_force, use_fist_force, use_mixed_force, &
       use_prog_name, use_qmmm, use_qmmmx, use_qs_force
  USE force_env_utils,                 ONLY: rescale_forces,&
                                             write_forces,&
                                             write_stress_tensor
  USE force_fields_util,               ONLY: get_generic_info
  USE fp_methods,                      ONLY: fp_eval
  USE fparser,                         ONLY: EvalErrType,&
                                             evalf,&
                                             evalfd,&
                                             finalizef,&
                                             initf,&
                                             parsef
  USE global_types,                    ONLY: global_environment_type,&
                                             globenv_retain
  USE input_constants,                 ONLY: &
       debug_run, do_fm_mom_conserv_QM, do_fm_mom_conserv_buffer, &
       do_fm_mom_conserv_core, do_fm_mom_conserv_equal_a, &
       do_fm_mom_conserv_equal_f, do_fm_mom_conserv_none, mix_coupled, &
       mix_generic, mix_linear_combination, mix_minimum, mix_restrained, &
       use_bazant_eip, use_lenosky_eip
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_retain,&
                                             section_vals_type,&
                                             section_vals_val_get,&
                                             section_vals_val_set
  USE kahan_sum,                       ONLY: accurate_sum
  USE kinds,                           ONLY: default_path_length,&
                                             default_string_length,&
                                             dp
  USE message_passing,                 ONLY: mp_sum,&
                                             mp_sync
  USE metadynamics_types,              ONLY: meta_env_retain,&
                                             meta_env_type
  USE mixed_energy_types,              ONLY: mixed_energy_type,&
                                             mixed_force_type
  USE mixed_environment_types,         ONLY: get_mixed_env,&
                                             mixed_env_retain,&
                                             mixed_environment_type
  USE mixed_environment_utils,         ONLY: get_subsys_map_index,&
                                             mixed_map_forces
  USE particle_list_types,             ONLY: particle_list_p_type,&
                                             particle_list_type
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qmmm_force,                      ONLY: qmmm_calc_energy_force
  USE qmmm_types,                      ONLY: qmmm_env_retain,&
                                             qmmm_env_type, qmmm_env_get
  USE qmmmx_types,                     ONLY: qmmmx_env_retain,&
                                             qmmmx_env_type
  USE qmmm_types_low,                  ONLY: force_mixing_label_QM_core,&
                                             force_mixing_label_QM_dynamics,&
                                             force_mixing_label_buffer
  USE qmmmx_util,                      ONLY: apply_qmmmx_translate
  USE qmmm_util,                       ONLY: apply_qmmm_translate,&
                                             apply_qmmm_unwrap,&
                                             apply_qmmm_wrap,&
                                             qmmm_force_mixing_active
  USE qs_energy,                       ONLY: qs_energies
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_retain,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_force,                        ONLY: qs_calc_energy_force
  USE restraint,                       ONLY: restraint_control
  USE string_utilities,                ONLY: compress
  USE termination,                     ONLY: print_warning
  USE virial_types,                    ONLY: cp_virial,&
                                             sym_virial,&
                                             virial_create,&
                                             virial_p_type,&
                                             virial_release,&
                                             virial_type,&
                                             zero_virial

  USE cell_types,                      ONLY: cell_type
  USE cp_output_handling,              ONLY: cp_iter_string,&
                                             cp_p_file,&
                                             cp_print_key_finished_output,&
                                             cp_print_key_should_output,&
                                             cp_print_key_unit_nr
  USE cp_result_methods,               ONLY: cp_results_erase,&
                                             get_results,&
                                             put_results
  USE cp_result_types,                 ONLY: cp_result_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get
  USE fist_force,                      ONLY: fist_calc_energy_force
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: default_string_length,&
                                             dp
  USE particle_types,                  ONLY: particle_type
  USE physcon,                         ONLY: debye
  USE qmmm_gpw_energy,                 ONLY: qmmm_el_coupling
  USE qmmm_gpw_forces,                 ONLY: qmmm_forces
  USE qmmm_links_methods,              ONLY: qmmm_added_chrg_coord,&
                                             qmmm_added_chrg_forces,&
                                             qmmm_link_Imomm_coord,&
                                             qmmm_link_Imomm_forces
  USE qmmm_types,                      ONLY: qmmm_env_type
  USE qmmmx_types,                      ONLY: qmmmx_env_type
  USE qmmm_types_low,                  ONLY: qmmm_links_type
  USE qmmm_util,                       ONLY: apply_qmmm_walls, apply_qmmm_translate
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env
  USE qs_force,                        ONLY: qs_calc_energy_force
  USE qs_ks_qmmm_methods,              ONLY: ks_qmmm_env_rebuild

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_force'

  PUBLIC :: qmmmx_calc_energy_force

CONTAINS

! *****************************************************************************
!> \brief calculates the qm/mm energy and forces
!> \param force_env ...
!> \param calc_force if also the forces should be calculated
!> \param require_consistent_energy_force ...
!> \param linres ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \par History
!>      05.2004 created [fawzi]
!> \author Fawzi Mohamed
! *****************************************************************************
  SUBROUTINE qmmmx_calc_energy_force(qmmmx_env,force_env_section,calc_force,&
         consistent_energies,linres,error)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    LOGICAL, INTENT(IN)                      :: calc_force, consistent_energies, linres
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_calc_energy_force', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: ip, &
                                                mom_conserv_min_label, &
                                                mom_conserv_n, &
                                                mom_conserv_region, &
                                                mom_conserv_type
    INTEGER, POINTER                         :: cur_indices(:), cur_labels(:)
    LOGICAL                                  ::  failure
    REAL(dp)                                 :: delta_a(3), delta_f(3), &
                                                mom_conserv_mass, total_f(3)
    TYPE(cp_subsys_type), POINTER            :: subsys_primary, &
                                                subsys_qmmm_core, &
                                                subsys_qmmm_extended
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particles_primary, &
                                                particles_qmmm_core, &
                                                particles_qmmm_extended

    IF (consistent_energies) THEN
      CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
          routineP,"qmmm_energy_and_forces got energy_consistency=.TRUE. but force mixing is active. "//&
          CPSourceFileRef,&
          error,failure)
    ENDIF

    ! Possibly translate the system
    CALL apply_qmmmx_translate(qmmmx_env, error)

    ! actual energy force calculation
    CALL qmmmx_calc_energy_force_low(qmmmx_env%ext, force_env_section, calc_force, consistent_energies, linres, error=error)
    CALL qmmmx_calc_energy_force_low(qmmmx_env%core, force_env_section, calc_force, consistent_energies, linres, error=error)

    ! get forces from subsys of each sub force env
    CALL qmmm_env_get(qmmmx_env%core,subsys=subsys_qmmm_core,error=error)
    CALL qmmm_env_get(qmmmx_env%ext ,subsys=subsys_qmmm_extended,error=error)

    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%RESTART_INFO%INDICES",i_vals=cur_indices,error=error)
    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%RESTART_INFO%LABELS",i_vals=cur_labels,error=error)

    particles_qmmm_extended => subsys_qmmm_extended%particles%els
    particles_qmmm_core => subsys_qmmm_core%particles%els
    DO ip=1,SIZE(cur_indices)
       IF (cur_labels(ip) >= force_mixing_label_QM_dynamics) THEN ! this is a QM atom
         ! copy (QM) force from extended calculation
         particles_qmmm_core(cur_indices(ip))%f=particles_qmmm_extended(cur_indices(ip))%f
       ENDIF
    END DO

    ! zero momentum
    CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%MOMENTUM_CONSERVATION_TYPE",&
                              i_val=mom_conserv_type,error=error)
    IF (mom_conserv_type /= do_fm_mom_conserv_none) THEN
       CALL section_vals_val_get(force_env_section,"QMMM%FORCE_MIXING%MOMENTUM_CONSERVATION_REGION",&
                                 i_val=mom_conserv_region,error=error)

       IF (mom_conserv_region == do_fm_mom_conserv_core) THEN
          mom_conserv_min_label = force_mixing_label_QM_core
       ELSEIF (mom_conserv_region == do_fm_mom_conserv_QM) THEN
          mom_conserv_min_label = force_mixing_label_QM_dynamics
       ELSEIF (mom_conserv_region == do_fm_mom_conserv_buffer) THEN
          mom_conserv_min_label = force_mixing_label_buffer
       ELSE
          CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
               routineP,"Got unknown MOMENTUM_CONSERVATION_REGION (not CORE, QM, or BUFFER) !"//&
             CPSourceFileRef,&
            error,failure)
       ENDIF

       total_f = 0.0_dp
       DO ip=1, SIZE(particles_qmmm_core)
         total_f(1:3) = total_f(1:3) + particles_qmmm_core(ip)%f(1:3)
       END DO
       IF (mom_conserv_type == do_fm_mom_conserv_equal_f) THEN
          mom_conserv_n = COUNT(cur_labels >= mom_conserv_min_label)
          delta_f = total_f/mom_conserv_n
          DO ip=1, SIZE(cur_indices)
            IF (cur_labels(ip) >= mom_conserv_min_label) THEN
              particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - delta_f
            ENDIF
          END DO
       ELSE IF (mom_conserv_type == do_fm_mom_conserv_equal_a) THEN
          mom_conserv_mass = 0.0_dp
          DO ip=1, SIZE(cur_indices)
           IF (cur_labels(ip) >= mom_conserv_min_label) &
             mom_conserv_mass = mom_conserv_mass + particles_qmmm_core(cur_indices(ip))%atomic_kind%mass
          END DO
          delta_a = total_f/mom_conserv_mass
          DO ip=1, SIZE(cur_indices)
            IF (cur_labels(ip) >= mom_conserv_min_label) THEN
              particles_qmmm_core(cur_indices(ip))%f = particles_qmmm_core(cur_indices(ip))%f - &
                 particles_qmmm_core(cur_indices(ip))%atomic_kind%mass * delta_a
            ENDIF
          END DO
       ENDIF
    ENDIF

    CALL qmmm_env_get(qmmmx_env%ext, subsys=subsys_primary, error=error)
    particles_primary => subsys_primary%particles%els
    DO ip=1,SIZE(particles_qmmm_core)
       particles_primary(ip)%f=particles_qmmm_core(ip)%f
    END DO

  END SUBROUTINE qmmmx_calc_energy_force

  ! *****************************************************************************
  ! *****************************************************************************
 SUBROUTINE qmmmx_calc_energy_force_low(qmmm_env,force_env_section,calc_force,consistent_energies,linres,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(section_vals_type), POINTER         :: force_env_section
    LOGICAL, INTENT(IN)                      :: calc_force, &
                                                consistent_energies, linres
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_calc_energy_force_low', &
      routineP = moduleN//':'//routineN

    TYPE(cell_type), POINTER                 :: mm_cell
    TYPE(cp_subsys_type), POINTER            :: subsys_qm, subsys_mm
    INTEGER, DIMENSION(:), POINTER           :: qm_atom_index
    REAL(dp), DIMENSION(:,:),ALLOCATABLE      :: saved_pos
    !CHARACTER(default_string_length)         :: restart_filename, &
    !                                            restart_filename_suffix, &
    !                                            restart_history_filename

    NULLIFY(mm_cell,subsys_qm,subsys_mm,qm_atom_index)

    !TODO: rewrite global restart file name
    !! maybe use this CALL get_qs_env(qmmmx_env%core%qs_env, input=force_env_section, error=error)
    !
    !CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
    !  c_val=restart_filename, error=error)
    !CALL section_vals_val_get(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
    !  c_val=restart_history_filename, error=error)
    !WRITE (unit=restart_filename_suffix, fmt=*) isubf
    !restart_filename_suffix= ADJUSTL(restart_filename_suffix)
    !CALL force_env_get(force_env%sub_force_env(isubf)%force_env, force_env_section=force_env_section,error=error)
    !CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
    !     c_val=TRIM(restart_filename)//"-SubForceEnv-"//TRIM(restart_filename_suffix), error=error)
    !CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
    !     c_val=TRIM(restart_history_filename)//"-SubForceEnv-"//TRIM(restart_filename_suffix), error=error)

    ! wrap positions before QM/MM calculation.
    ! Required if diffusion causes atoms outside of periodic box get added to QM
    CALL fist_env_get(qmmm_env%fist_env, cell=mm_cell,subsys=subsys_mm,error=error)
    CALL get_qs_env(qmmm_env%qs_env, cp_subsys=subsys_qm,error=error)
    qm_atom_index => qmmm_env%qm%qm_atom_index
    CALL apply_qmmm_wrap(subsys_mm, mm_cell, subsys_qm, qm_atom_index, saved_pos, error)

    ! actual energy force calculation
    CALL qmmm_calc_energy_force(qmmm_env, calc_force, consistent_energies, linres, translate=.FALSE., error=error)

    ! restore unwrapped positions
    CALL apply_qmmm_unwrap(subsys_mm, subsys_qm, qm_atom_index, saved_pos, error)

    !TODO: restore global restart file name
    !CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART%FILENAME", &
    !     c_val=TRIM(restart_filename), error=error)
    !CALL section_vals_val_set(force_env_section, "DFT%SCF%PRINT%RESTART_HISTORY%FILENAME", &
    !     c_val=TRIM(restart_history_filename), error=error)

  END SUBROUTINE qmmmx_calc_energy_force_low

END MODULE qmmmx_force
