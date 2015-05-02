!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \author Ole Schuett
! *****************************************************************************
MODULE qmmmx_types
  USE qmmm_types,                      ONLY: qmmm_env_release,&
                                             qmmm_env_retain,&
                                             qmmm_env_type

#include "./common/cp_common_uses.f90"

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qmmmx_types'

  PUBLIC :: qmmmx_env_type, qmmmx_env_retain, qmmmx_env_release, qmmmx_env_create

  ! *****************************************************************************
  TYPE qmmmx_env_type
     INTEGER                                                 :: ref_count = -1
     TYPE(qmmm_env_type), POINTER                            :: core => Null()
     TYPE(qmmm_env_type), POINTER                            :: ext => Null()
  END TYPE qmmmx_env_type

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
  SUBROUTINE qmmmx_env_create(qmmmx_env,  error)
    TYPE(qmmmx_env_type), POINTER            :: qmmmx_env
    TYPE(qmmm_env_type), POINTER             :: qmmm_core
    TYPE(qmmm_env_type), POINTER             :: qmmm_ext
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'qmmmx_env_create', &
      routineP = moduleN//':'//routineN

    IF(ASSOCIATED(qmmmx_env)) STOP "qmmmx_env_create: qmmmx_env already associated"
    IF(.NOT. ASSOCIATED(qmmm_core)) STOP "qmmmx_env_create: qmmm_core not associated"
    IF(.NOT. ASSOCIATED(qmmm_ext)) STOP "qmmmx_env_create: qmmm_ext not associated"

    ALLOCATE(qmmmx_env)
    qmmmx_env%ref_count = 1
    qmmmx_env%core => qmmm_core
    qmmmx_env%ext => qmmm_ext

    CALL qmmm_env_retain(qmmm_core, error)
    CALL qmmm_env_retain(qmmm_ext, error)

  END SUBROUTINE qmmmx_env_create

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
