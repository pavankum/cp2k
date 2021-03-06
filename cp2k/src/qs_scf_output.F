!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

MODULE qs_scf_output
   USE atomic_kind_types,               ONLY: atomic_kind_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: cp_dbcsr_p_type
   USE cp_dbcsr_output,                 ONLY: cp_dbcsr_write_sparse_matrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_units,                        ONLY: cp_unit_from_cp2k
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kahan_sum,                       ONLY: accurate_sum
   USE kinds,                           ONLY: dp
   USE machine,                         ONLY: m_flush
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt,&
                                              kcalmol
   USE ps_implicit_types,               ONLY: MIXED_BC,&
                                              MIXED_PERIODIC_BC,&
                                              NEUMANN_BC,&
                                              PERIODIC_BC
   USE pw_env_types,                    ONLY: pw_env_type
   USE pw_poisson_types,                ONLY: pw_poisson_implicit
   USE qmmm_image_charge,               ONLY: print_image_coefficients
   USE qs_charges_types,                ONLY: qs_charges_type
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_mo_io,                        ONLY: write_mo_set
   USE qs_mo_methods,                   ONLY: calculate_magnitude,&
                                              calculate_orthonormality
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_scf_types,                    ONLY: qs_scf_env_type,&
                                              special_diag_method_nr
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_scf_output'

   PUBLIC :: qs_scf_loop_info, &
             qs_scf_print_summary, &
             qs_scf_loop_print, &
             qs_scf_outer_loop_info, &
             qs_scf_initial_info, &
             qs_scf_write_mos

CONTAINS

! **************************************************************************************************
!> \brief writes a summary of information after scf
!> \param output_unit ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE qs_scf_print_summary(output_unit, qs_env)
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_print_summary', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: nelectron_total
      LOGICAL                                            :: gapw, gapw_xc, qmmm
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_scf_env_type), POINTER                     :: scf_env

      NULLIFY (rho, energy, dft_control, scf_env, qs_charges)
      CALL get_qs_env(qs_env=qs_env, rho=rho, energy=energy, dft_control=dft_control, &
                      scf_env=scf_env, qs_charges=qs_charges)

      gapw = dft_control%qs_control%gapw
      gapw_xc = dft_control%qs_control%gapw_xc
      qmmm = qs_env%qmmm
      nelectron_total = scf_env%nelectron

      CALL qs_scf_print_scf_summary(output_unit, rho, qs_charges, energy, nelectron_total, &
                                    dft_control, qmmm, qs_env, gapw, gapw_xc)

   END SUBROUTINE qs_scf_print_summary

! **************************************************************************************************
!> \brief writes basic information at the beginning of an scf run
!> \param output_unit ...
!> \param mos ...
!> \param dft_control ...
! **************************************************************************************************
   SUBROUTINE qs_scf_initial_info(output_unit, mos, dft_control)
      INTEGER                                            :: output_unit
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(dft_control_type), POINTER                    :: dft_control

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_initial_info', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: homo, ispin, nao, nelectron_spin, nmo

! print occupation numbers

      IF (output_unit > 0) THEN
         DO ispin = 1, dft_control%nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, &
                            homo=homo, &
                            nelectron=nelectron_spin, &
                            nao=nao, &
                            nmo=nmo)
            IF (dft_control%nspins > 1) THEN
               WRITE (UNIT=output_unit, FMT="(/,T2,A,I2)") "Spin", ispin
            END IF
            WRITE (UNIT=output_unit, FMT="(/,(T2,A,T71,I10))") &
               "Number of electrons:", nelectron_spin, &
               "Number of occupied orbitals:", homo, &
               "Number of molecular orbitals:", nmo
         END DO
         WRITE (UNIT=output_unit, FMT="(/,T2,A,T71,I10)") &
            "Number of orbital functions:", nao
      END IF

   END SUBROUTINE qs_scf_initial_info

! **************************************************************************************************
!> \brief writes out the mos in the scf loop if needed
!> \param mos ...
!> \param atomic_kind_set ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param dft_section ...
! **************************************************************************************************
   SUBROUTINE qs_scf_write_mos(mos, atomic_kind_set, qs_kind_set, particle_set, dft_section)
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(section_vals_type), POINTER                   :: dft_section

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_write_mos', &
         routineP = moduleN//':'//routineN

      IF (SIZE(mos) > 1) THEN
         CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                           dft_section, spin="ALPHA", last=.FALSE.)
         CALL write_mo_set(mos(2)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                           dft_section, spin="BETA", last=.FALSE.)
      ELSE
         CALL write_mo_set(mos(1)%mo_set, atomic_kind_set, qs_kind_set, particle_set, 4, &
                           dft_section, last=.FALSE.)
      END IF

   END SUBROUTINE qs_scf_write_mos
! **************************************************************************************************
!> \brief writes basic information obtained in a scf outer loop step
!> \param output_unit ...
!> \param scf_control ...
!> \param scf_env ...
!> \param energy ...
!> \param total_steps ...
!> \param should_stop ...
!> \param outer_loop_converged ...
! **************************************************************************************************

   SUBROUTINE qs_scf_outer_loop_info(output_unit, scf_control, scf_env, &
                                     energy, total_steps, should_stop, outer_loop_converged)
      INTEGER                                            :: output_unit
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER                                            :: total_steps
      LOGICAL, INTENT(IN)                                :: should_stop, outer_loop_converged

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_outer_loop_info', &
         routineP = moduleN//':'//routineN

      REAL(KIND=dp)                                      :: outer_loop_eps

      outer_loop_eps = SQRT(SUM(scf_env%outer_scf%gradient(:, scf_env%outer_scf%iter_count)**2))/ &
                       SIZE(scf_env%outer_scf%gradient, 1)
      IF (output_unit > 0) WRITE (output_unit, '(/,T3,A,I4,A,E10.2,A,F22.10)') &
         "outer SCF iter = ", scf_env%outer_scf%iter_count, &
         " RMS gradient = ", outer_loop_eps, " energy =", energy%total

      IF (outer_loop_converged) THEN
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
            "outer SCF loop converged in", scf_env%outer_scf%iter_count, &
            " iterations or ", total_steps, " steps"
      END IF

      IF (scf_env%outer_scf%iter_count > scf_control%outer_scf%max_scf &
          .OR. should_stop) THEN
         IF (output_unit > 0) WRITE (output_unit, '(T3,A,I4,A,I4,A,/)') &
            "outer SCF loop FAILED to converge after ", &
            scf_env%outer_scf%iter_count, " iterations or ", total_steps, " steps"
      END IF

   END SUBROUTINE qs_scf_outer_loop_info

! **************************************************************************************************
!> \brief writes basic information obtained in a scf step
!> \param scf_env ...
!> \param output_unit ...
!> \param just_energy ...
!> \param t1 ...
!> \param t2 ...
!> \param energy ...
! **************************************************************************************************

   SUBROUTINE qs_scf_loop_info(scf_env, output_unit, just_energy, t1, t2, energy)

      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      INTEGER                                            :: output_unit
      LOGICAL                                            :: just_energy
      REAL(KIND=dp)                                      :: t1, t2
      TYPE(qs_energy_type), POINTER                      :: energy

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_info', &
         routineP = moduleN//':'//routineN

      IF ((output_unit > 0) .AND. scf_env%print_iter_line) THEN
         IF (just_energy) THEN
            WRITE (UNIT=output_unit, &
                   FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,16X,F20.10)") &
               scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
               t2-t1, energy%total
         ELSE
            IF ((ABS(scf_env%iter_delta) < 1.0E-8_dp) .OR. &
                (ABS(scf_env%iter_delta) >= 1.0E5_dp)) THEN
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,1X,ES14.4,1X,F20.10,1X,ES9.2)") &
                  scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
                  t2-t1, scf_env%iter_delta, energy%total, energy%total-energy%tot_old
            ELSE
               WRITE (UNIT=output_unit, &
                      FMT="(T2,I5,1X,A,T20,E8.2,1X,F6.1,1X,F14.8,1X,F20.10,1X,ES9.2)") &
                  scf_env%iter_count, TRIM(scf_env%iter_method), scf_env%iter_param, &
                  t2-t1, scf_env%iter_delta, energy%total, energy%total-energy%tot_old
            END IF
         END IF
      END IF

   END SUBROUTINE qs_scf_loop_info

! **************************************************************************************************
!> \brief writes rather detailed summary of densities and energies
!>      after the SCF
!> \param output_unit ...
!> \param rho ...
!> \param qs_charges ...
!> \param energy ...
!> \param nelectron_total ...
!> \param dft_control ...
!> \param qmmm ...
!> \param qs_env ...
!> \param gapw ...
!> \param gapw_xc ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE qs_scf_print_scf_summary(output_unit, rho, qs_charges, energy, nelectron_total, &
                                       dft_control, qmmm, qs_env, gapw, gapw_xc)
      INTEGER, INTENT(IN)                                :: output_unit
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(qs_charges_type), POINTER                     :: qs_charges
      TYPE(qs_energy_type), POINTER                      :: energy
      INTEGER, INTENT(IN)                                :: nelectron_total
      TYPE(dft_control_type), POINTER                    :: dft_control
      LOGICAL, INTENT(IN)                                :: qmmm
      TYPE(qs_environment_type), POINTER                 :: qs_env
      LOGICAL, INTENT(IN)                                :: gapw, gapw_xc

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_print_scf_summary', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: bc, handle, ispin, psolver
      REAL(kind=dp)                                      :: implicit_ps_ehartree, tot1_h, tot1_s
      REAL(KIND=dp), DIMENSION(:), POINTER               :: tot_rho_r
      TYPE(pw_env_type), POINTER                         :: pw_env

      NULLIFY (tot_rho_r, pw_env)
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      psolver = pw_env%poisson_env%parameters%solver

      IF (output_unit > 0) THEN
         CALL qs_rho_get(rho, tot_rho_r=tot_rho_r)
         IF (.NOT. (dft_control%qs_control%semi_empirical .OR. dft_control%qs_control%dftb .OR. &
                    dft_control%qs_control%scptb)) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T41,2F20.10))") &
               "Electronic density on regular grids: ", &
               accurate_sum(tot_rho_r), &
               accurate_sum(tot_rho_r)+nelectron_total, &
               "Core density on regular grids:", &
               qs_charges%total_rho_core_rspace, &
               qs_charges%total_rho_core_rspace-REAL(nelectron_total+dft_control%charge, dp)
            IF (gapw) THEN
               tot1_h = qs_charges%total_rho1_hard(1)
               tot1_s = qs_charges%total_rho1_soft(1)
               DO ispin = 2, dft_control%nspins
                  tot1_h = tot1_h+qs_charges%total_rho1_hard(ispin)
                  tot1_s = tot1_s+qs_charges%total_rho1_soft(ispin)
               END DO
               WRITE (UNIT=output_unit, FMT="((T3,A,T41,2F20.10))") &
                  "Hard and soft densities (Lebedev):", &
                  tot1_h, tot1_s
               WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
                  "Total Rho_soft + Rho1_hard - Rho1_soft (r-space): ", &
                  accurate_sum(tot_rho_r)+tot1_h-tot1_s, &
                  "Total charge density (r-space):      ", &
                  accurate_sum(tot_rho_r)+tot1_h-tot1_s &
                  +qs_charges%total_rho_core_rspace, &
                  "Total Rho_soft + Rho0_soft (g-space):", &
                  qs_charges%total_rho_gspace
            ELSE
               WRITE (UNIT=output_unit, FMT="(T3,A,T41,F20.10)") &
                  "Total charge density on r-space grids:     ", &
                  accurate_sum(tot_rho_r)+ &
                  qs_charges%total_rho_core_rspace, &
                  "Total charge density g-space grids:     ", &
                  qs_charges%total_rho_gspace
            END IF
         END IF
         IF (dft_control%qs_control%semi_empirical) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Core-core repulsion energy [eV]:               ", energy%core_overlap*evolt, &
               "Core Hamiltonian energy [eV]:                  ", energy%core*evolt, &
               "Two-electron integral energy [eV]:             ", energy%hartree*evolt, &
               "Electronic energy [eV]:                        ", &
               (energy%core+0.5_dp*energy%hartree)*evolt
            IF (energy%dispersion /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Dispersion energy [eV]:                     ", energy%dispersion*evolt
         ELSEIF (dft_control%qs_control%dftb) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Core Hamiltonian energy:                       ", energy%core, &
               "Repulsive potential energy:                    ", energy%repulsive, &
               "Electronic energy:                             ", energy%hartree, &
               "Dispersion energy:                             ", energy%dispersion
            IF (energy%dftb3 /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "DFTB3 3rd order energy:                     ", energy%dftb3
         ELSEIF (dft_control%qs_control%scptb) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Repulsive pair potential energy:               ", energy%repulsive, &
               "Zeroth order Hamiltonian energy:               ", energy%core, &
               "Kinetic energy:                                ", energy%kinetic, &
               "Charge fluctuation energy:                     ", energy%hartree, &
               "London dispersion energy:                      ", energy%dispersion
         ELSE
            IF (dft_control%do_admm) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "Overlap energy of the core charge distribution:", energy%core_overlap, &
                  "Self energy of the core charge distribution:   ", energy%core_self, &
                  "Core Hamiltonian energy:                       ", energy%core, &
                  "Hartree energy:                                ", energy%hartree, &
                  "Exchange-correlation energy:                   ", energy%exc+energy%exc_aux_fit
            ELSE
               IF (psolver .EQ. pw_poisson_implicit) THEN
                  implicit_ps_ehartree = pw_env%poisson_env%implicit_env%ehartree
                  bc = pw_env%poisson_env%parameters%ps_implicit_params%boundary_condition
                  SELECT CASE (bc)
                  CASE (MIXED_PERIODIC_BC, MIXED_BC)
                     WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                        "Overlap energy of the core charge distribution:", energy%core_overlap, &
                        "Self energy of the core charge distribution:   ", energy%core_self, &
                        "Core Hamiltonian energy:                       ", energy%core, &
                        "Hartree energy:                                ", implicit_ps_ehartree, &
                        "Electric enthalpy:                             ", energy%hartree, &
                        "Exchange-correlation energy:                   ", energy%exc
                  CASE (PERIODIC_BC, NEUMANN_BC)
                     WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                        "Overlap energy of the core charge distribution:", energy%core_overlap, &
                        "Self energy of the core charge distribution:   ", energy%core_self, &
                        "Core Hamiltonian energy:                       ", energy%core, &
                        "Hartree energy:                                ", energy%hartree, &
                        "Exchange-correlation energy:                   ", energy%exc
                  END SELECT
               ELSE
                  WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                     "Overlap energy of the core charge distribution:", energy%core_overlap, &
                     "Self energy of the core charge distribution:   ", energy%core_self, &
                     "Core Hamiltonian energy:                       ", energy%core, &
                     "Hartree energy:                                ", energy%hartree, &
                     "Exchange-correlation energy:                   ", energy%exc
               END IF
            END IF
            IF (energy%e_hartree /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,/,T3,A,T56,F25.14)") &
               "Coulomb Electron-Electron Interaction Energy ", &
               "- Already included in the total Hartree term ", energy%e_hartree
            IF (energy%ex /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Hartree-Fock Exchange energy:                  ", energy%ex
            IF (energy%dispersion /= 0.0_dp) &
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Dispersion energy:                             ", energy%dispersion
            IF (gapw) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "GAPW| Exc from hard and soft atomic rho1:      ", energy%exc1, &
                  "GAPW| local Eh = 1 center integrals:           ", energy%hartree_1c
            END IF
            IF (gapw_xc) THEN
               WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
                  "GAPW_XC| Exc from hard and soft atomic rho1:      ", energy%exc1
            END IF
         END IF
         IF (dft_control%smear) THEN
            WRITE (UNIT=output_unit, FMT="((T3,A,T56,F25.14))") &
               "Electronic entropic energy:", energy%kTS
            WRITE (UNIT=output_unit, FMT="((T3,A,T56,F25.14))") &
               "Fermi energy:", energy%efermi
         END IF
         IF (dft_control%dft_plus_u) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "DFT+U energy:", energy%dft_plus_u
         END IF
         IF (dft_control%do_sccs) THEN
            WRITE (UNIT=output_unit, FMT="(/,T3,A,T56,F25.14)") &
               "SCCS| Hartree energy of solute and solvent [Hartree]", energy%sccs_hartree
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14,/,T3,A,T61,F20.3)") &
               "SCCS| Polarisation energy                  [Hartree]", energy%sccs_pol, &
               "SCCS|                                      [kcal/mol]", &
               cp_unit_from_cp2k(energy%sccs_pol, "kcalmol")
         END IF
         IF (qmmm) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "QM/MM Electrostatic energy:                    ", energy%qmmm_el
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
                  "QM/MM image charge energy:                ", energy%image_charge
            ENDIF
         END IF
         IF (dft_control%qs_control%mulliken_restraint) THEN
            WRITE (UNIT=output_unit, FMT="(T3,A,T56,F25.14)") &
               "Mulliken restraint energy: ", energy%mulliken
         END IF
         IF (dft_control%qs_control%semi_empirical) THEN
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Total energy [eV]:                             ", energy%total*evolt
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Atomic reference energy [eV]:                  ", energy%core_self*evolt, &
               "Heat of formation [kcal/mol]:                  ", &
               (energy%total+energy%core_self)*kcalmol
         ELSE
            WRITE (UNIT=output_unit, FMT="(/,(T3,A,T56,F25.14))") &
               "Total energy:                                  ", energy%total
         END IF
         IF (qmmm) THEN
            IF (qs_env%qmmm_env_qm%image_charge) THEN
               CALL print_image_coefficients(qs_env%image_coeff, qs_env)
            ENDIF
         ENDIF
         CALL m_flush(output_unit)
      END IF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_print_scf_summary

! **************************************************************************************************
!> \brief collects the 'heavy duty' printing tasks out of the SCF loop
!> \param qs_env ...
!> \param scf_env ...
!> \param para_env ...
!> \par History
!>      03.2006 created [Joost VandeVondele]
! **************************************************************************************************
   SUBROUTINE qs_scf_loop_print(qs_env, scf_env, para_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(cp_para_env_type), POINTER                    :: para_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'qs_scf_loop_print', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: after, handle, ic, ispin, iw
      LOGICAL                                            :: do_kpoints, omit_headers
      REAL(KIND=dp)                                      :: mo_mag_max, mo_mag_min, orthonormality
      TYPE(cp_dbcsr_p_type), DIMENSION(:, :), POINTER    :: matrix_ks, matrix_p, matrix_s
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: dft_section, input, scf_section

      logger => cp_get_default_logger()
      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env=qs_env, input=input, dft_control=dft_control, &
                      do_kpoints=do_kpoints)

      dft_section => section_vals_get_subs_vals(input, "DFT")
      scf_section => section_vals_get_subs_vals(dft_section, "SCF")

      CALL section_vals_val_get(input, "DFT%PRINT%AO_MATRICES%OMIT_HEADERS", l_val=omit_headers)
      DO ispin = 1, dft_control%nspins

         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "PRINT%AO_MATRICES/DENSITY"), cp_p_file)) THEN
            CALL get_qs_env(qs_env, rho=rho)
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%AO_MATRICES/DENSITY", &
                                      extension=".Log")
            CALL section_vals_val_get(dft_section, "PRINT%AO_MATRICES%NDIGITS", i_val=after)
            after = MIN(MAX(after, 1), 16)
            DO ic = 1, SIZE(matrix_p, 2)
               CALL cp_dbcsr_write_sparse_matrix(matrix_p(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                 output_unit=iw, omit_headers=omit_headers)
            END DO
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              "PRINT%AO_MATRICES/DENSITY")
         END IF

         IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                              dft_section, "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX"), cp_p_file)) THEN
            iw = cp_print_key_unit_nr(logger, dft_section, "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX", &
                                      extension=".Log")
            CALL section_vals_val_get(dft_section, "PRINT%AO_MATRICES%NDIGITS", i_val=after)
            after = MIN(MAX(after, 1), 16)
            CALL get_qs_env(qs_env=qs_env, matrix_ks_kp=matrix_ks)
            DO ic = 1, SIZE(matrix_ks, 2)
               IF (dft_control%qs_control%semi_empirical) THEN
                  CALL cp_dbcsr_write_sparse_matrix(matrix_ks(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                    scale=evolt, output_unit=iw, omit_headers=omit_headers)
               ELSE
                  CALL cp_dbcsr_write_sparse_matrix(matrix_ks(ispin, ic)%matrix, 4, after, qs_env, para_env, &
                                                    output_unit=iw, omit_headers=omit_headers)
               END IF
            END DO
            CALL cp_print_key_finished_output(iw, logger, dft_section, &
                                              "PRINT%AO_MATRICES/KOHN_SHAM_MATRIX")
         END IF

      ENDDO

      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           scf_section, "PRINT%MO_ORTHONORMALITY"), cp_p_file)) THEN
         IF (do_kpoints) THEN
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_ORTHONORMALITY", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A)') &
                  " K-points: Maximum deviation from MO S-orthonormality not determined"
            ENDIF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_ORTHONORMALITY")
         ELSE
            CALL get_qs_env(qs_env, mos=mos)
            IF (scf_env%method == special_diag_method_nr) THEN
               CALL calculate_orthonormality(orthonormality, mos)
            ELSE
               CALL get_qs_env(qs_env=qs_env, matrix_s_kp=matrix_s)
               CALL calculate_orthonormality(orthonormality, mos, matrix_s(1, 1)%matrix)
            END IF
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_ORTHONORMALITY", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A,T61,E20.4)') &
                  " Maximum deviation from MO S-orthonormality", orthonormality
            ENDIF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_ORTHONORMALITY")
         END IF
      ENDIF
      IF (BTEST(cp_print_key_should_output(logger%iter_info, &
                                           scf_section, "PRINT%MO_MAGNITUDE"), cp_p_file)) THEN
         IF (do_kpoints) THEN
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_MAGNITUDE", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A)') &
                  " K-points: Minimum/Maximum MO magnitude not determined"
            ENDIF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_MAGNITUDE")
         ELSE
            CALL get_qs_env(qs_env, mos=mos)
            CALL calculate_magnitude(mos, mo_mag_min, mo_mag_max)
            iw = cp_print_key_unit_nr(logger, scf_section, "PRINT%MO_MAGNITUDE", &
                                      extension=".scfLog")
            IF (iw > 0) THEN
               WRITE (iw, '(T8,A,T41,2E20.4)') &
                  " Minimum/Maximum MO magnitude ", mo_mag_min, mo_mag_max
            ENDIF
            CALL cp_print_key_finished_output(iw, logger, scf_section, &
                                              "PRINT%MO_MAGNITUDE")
         END IF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE qs_scf_loop_print

END MODULE qs_scf_output
