!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_types
 USE fist_environment_types,          ONLY: fist_env_release,&
                                             fist_env_retain,&
                                             fist_environment_type
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_retain,&
                                             qmmm_env_qm_type
  USE qs_environment_types,            ONLY: qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE fist_energy_types,               ONLY: fist_energy_type
  USE fist_environment_types,          ONLY: fist_env_get,&
                                             fist_env_release,&
                                             fist_environment_type
  USE cp_subsys_types,                 ONLY: cp_subsys_get,&
                                             cp_subsys_type,&
                                             pack_subsys_particles
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_env_release,&
                                             qs_environment_type
  USE qmmm_types,                      ONLY: qmmm_env_release,&
                                             qmmm_env_retain,&
                                             qmmm_env_type,&
                                             qmmm_env_get
  USE kinds,                           ONLY: dp

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_types'

  PUBLIC :: qmmmx_env_type, qmmmx_env_get, qmmmx_env_retain, qmmmx_env_release

  ! *****************************************************************************
  TYPE qmmmx_env_type
     INTEGER                                                 :: ref_count = 1
     TYPE(qmmm_env_type), POINTER                            :: core => Null()
     TYPE(qmmm_env_type), POINTER                            :: ext => Null()
  END TYPE qmmmx_env_type

CONTAINS

! *****************************************************************************
! *****************************************************************************
  SUBROUTINE qmmmx_env_get(qmmmx_env,subsys,potential_energy,kinetic_energy,error)
    TYPE(qmmmx_env_type), POINTER             :: qmmmx_env
    TYPE(cp_subsys_type), OPTIONAL, POINTER  :: subsys
    REAL(KIND=dp), INTENT(OUT), OPTIONAL     :: potential_energy, &
                                                kinetic_energy
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_get', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(qs_energy_type), POINTER            :: qs_energy
    TYPE(fist_energy_type), POINTER          :: thermo
    failure=.FALSE.
    NULLIFY(qs_energy, thermo)

    CPPrecondition(ASSOCIATED(qmmmx_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(qmmmx_env%ref_count>0,cp_failure_level,routineP,error,failure)

    ! get the underlying energies from primary subsys.  This is the only subsys
    ! for conventional QM/MM, and force-mixing knows to put relevant energy there.
    CALL qmmm_env_get(qmmmx_env%ext,&
                      kinetic_energy=kinetic_energy,&
                      potential_energy=potential_energy,&
                      subsys=subsys,&
                      error=error)

  END SUBROUTINE qmmmx_env_get


! *****************************************************************************
!> \brief ...
!> \param qmmmx_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_retain(qmmmx_env,error)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qmmmx_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(qmmmx_env%ref_count>0,cp_failure_level,routineP,error)
       qmmmx_env%ref_count = qmmmx_env%ref_count+1
    END IF
  END SUBROUTINE qmmmx_env_retain

! *****************************************************************************
!> \brief releases the given qmmmx_env (see doc/ReferenceCounting.html)
!> \param qmmm_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmmx_env_release(qmmmx_env,error)
    TYPE(qmmmx_env_type), POINTER             :: qmmmx_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(qmmmx_env)) THEN
       CPPreconditionNoFail(qmmmx_env%ref_count>0,cp_failure_level,routineP,error)
       qmmmx_env%ref_count = qmmmx_env%ref_count-1
       IF (qmmmx_env%ref_count==0) THEN
          CALL qmmm_env_release(qmmmx_env%core, error)
          CALL qmmm_env_release(qmmmx_env%ext, error)
          DEALLOCATE(qmmmx_env)
       END IF
    END IF
    NULLIFY(qmmmx_env)
  END SUBROUTINE qmmmx_env_release

END MODULE qmmmx_types
