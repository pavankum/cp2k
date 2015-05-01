!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmm_types
  USE fist_environment_types,          ONLY: fist_env_release,&
                                             fist_env_retain,&
                                             fist_environment_type
  USE qmmm_types_low,                  ONLY: qmmm_env_qm_release,&
                                             qmmm_env_qm_retain,&
                                             qmmm_env_qm_type
  USE qs_environment_types,            ONLY: qs_env_release,&
                                             qs_env_retain,&
                                             qs_environment_type
#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmm_types'

  PUBLIC :: qmmm_env_type, qmmm_env_retain, qmmm_env_release, qmmm_env_create

 ! *****************************************************************************
  TYPE qmmm_env_type
     INTEGER                                                 :: ref_count = -1
     TYPE(qs_environment_type), POINTER                      :: qs_env => Null()
     TYPE(fist_environment_type), POINTER                    :: fist_env => Null()
     TYPE(qmmm_env_qm_type), POINTER                         :: qm => Null()
  END TYPE qmmm_env_type

CONTAINS

! *****************************************************************************
!> \brief TODO
!> \param qmmm_env ...
!> \param qs_env ...
!> \param fist_env ...
!> \param qmmm_env_qm ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_create(qmmm_env, qs_env, fist_env, qmmm_env_qm, error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(fist_environment_type), POINTER     :: fist_env
    TYPE(qmmm_env_qm_type), POINTER          :: qmmm_env_qm
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_create', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(qmmm_env)) STOP "qmmm_env_create: qmmm_env already associated"
    IF(.NOT. ASSOCIATED(qs_env)) STOP "qmmm_env_create: qs_env not associated"
    IF(.NOT. ASSOCIATED(fist_env)) STOP "qmmm_env_create: fist_env not associated"
    IF(.NOT. ASSOCIATED(qmmm_env_qm)) STOP "qmmm_env_create: qmmm_env_qm not associated"

    ALLOCATE(qmmm_env)
    qmmm_env%ref_count = 1
    qmmm_env%qs_env => qs_env
    qmmm_env%fist_env => fist_env
    qmmm_env%qm => qmmm_env_qm

    CALL qs_env_retain(qs_env, error)
    CALL fist_env_retain(fist_env, error)
    CALL qmmm_env_qm_retain(qmmm_env_qm, error)

  END SUBROUTINE qmmm_env_create

! *****************************************************************************
!> \brief ...
!> \param qmmm_env ...
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_retain(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_retain', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    CPPrecondition(ASSOCIATED(qmmm_env),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
       qmmm_env%ref_count=qmmm_env%ref_count+1
    END IF
  END SUBROUTINE qmmm_env_retain

! *****************************************************************************
!> \brief releases the given qmmm_env (see doc/ReferenceCounting.html)
!> \param qmmm_env the object to release
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
!> \author Ole Schuett
! *****************************************************************************
  SUBROUTINE qmmm_env_release(qmmm_env,error)
    TYPE(qmmm_env_type), POINTER             :: qmmm_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmm_env_release', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure

    failure=.FALSE.

    IF (ASSOCIATED(qmmm_env)) THEN
       CPPreconditionNoFail(qmmm_env%ref_count>0,cp_failure_level,routineP,error)
       qmmm_env%ref_count=qmmm_env%ref_count-1
       IF (qmmm_env%ref_count==0) THEN
          CALL qs_env_release(qmmm_env%qs_env, error)
          CALL fist_env_release(qmmm_env%fist_env, error)
          CALL qmmm_env_qm_release(qmmm_env%qm, error)
          DEALLOCATE(qmmm_env)
       END IF
    END IF
    NULLIFY(qmmm_env)
  END SUBROUTINE qmmm_env_release

END MODULE qmmm_types
