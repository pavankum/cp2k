!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2015  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Helper routines shared by multiple PAO modules
!> \author Ole Schuett
! *****************************************************************************
MODULE pao_utils
  USE cp_blacs_env,                    ONLY: cp_blacs_env_type
  USE cp_dbcsr_diag,                   ONLY: cp_dbcsr_syevd
  USE cp_dbcsr_interface,              ONLY: &
       cp_dbcsr_create, cp_dbcsr_init, cp_dbcsr_iterator, &
       cp_dbcsr_iterator_blocks_left, cp_dbcsr_iterator_next_block, &
       cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, cp_dbcsr_multiply, &
       cp_dbcsr_nfullcols_total, cp_dbcsr_nfullrows_total, cp_dbcsr_release, &
       cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dm_ls_scf_types,                 ONLY: ls_scf_env_type
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: create_rng_stream,&
                                             delete_rng_stream,&
                                             next_random_number,&
                                             rng_stream_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
#include "./base/base_uses.f90"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pao_utils'

  PUBLIC :: fill_randomly

CONTAINS

! *****************************************************************************
!> \brief TODO
!> \param qs_env ...
!> \param ls_scf_env ...
! *****************************************************************************
  SUBROUTINE print_eigenvalues_PS(qs_env, ls_scf_env)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(ls_scf_env_type)                    :: ls_scf_env

    TYPE(cp_dbcsr_type)                      :: matrix_x

!TYPE(cp_dbcsr_type)                      :: matrix_PS
!CALL cp_dbcsr_init(matrix_PS)
!CALL cp_dbcsr_create(matrix_PS, template=ls_scf_env%matrix_s, matrix_type="N")
!CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_p(1), ls_scf_env%matrix_s,&
!                       0.0_dp, matrix_PS)
!CALL print_eigenvalues(qs_env, matrix_PS)
!CALL cp_dbcsr_release(matrix_PS)

    CALL cp_dbcsr_init(matrix_x)
    CALL cp_dbcsr_create(matrix_x, template=ls_scf_env%matrix_p(1), matrix_type="N")
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, ls_scf_env%matrix_s_sqrt, ls_scf_env%matrix_p(1),&
                           0.0_dp, matrix_x)
    CALL cp_dbcsr_multiply("N", "N", 1.0_dp, matrix_x, ls_scf_env%matrix_s_sqrt, &
                          0.0_dp, matrix_x)
    CALL print_eigenvalues(qs_env, matrix_x)
    CALL cp_dbcsr_release(matrix_x)
  END SUBROUTINE print_eigenvalues_PS


! *****************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param matrix ...
! *****************************************************************************
  SUBROUTINE print_eigenvalues(qs_env, matrix)
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_dbcsr_type)                      :: matrix

    INTEGER                                  :: k
    REAL(dp), DIMENSION(:), POINTER          :: eigenvalues
    TYPE(cp_blacs_env_type), POINTER         :: blacs_env
    TYPE(cp_dbcsr_type)                      :: eigenvectors
    TYPE(cp_para_env_type), POINTER          :: para_env

    CPASSERT(cp_dbcsr_nfullrows_total(matrix)==cp_dbcsr_nfullcols_total(matrix))
    k = cp_dbcsr_nfullrows_total(matrix)
    CALL cp_dbcsr_init(eigenvectors)
    CALL cp_dbcsr_create(eigenvectors, template=matrix, matrix_type="N")
    ALLOCATE(eigenvalues(k))
    CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)
    CALL cp_dbcsr_syevd(matrix, eigenvectors, eigenvalues, para_env, blacs_env)
    WRITE (*,*) "eigenvalues(PS)", eigenvalues
    WRITE (*,*) "SUM(eigenvalues(PS))", SUM(eigenvalues)
    DEALLOCATE(eigenvalues)
    CALL cp_dbcsr_release(eigenvectors)
  END SUBROUTINE print_eigenvalues


! *****************************************************************************
!> \brief Fills a DBCSR matrix with random noise.
!> \param matrix ...
!> \param seed ...
! *****************************************************************************
  SUBROUTINE fill_randomly(matrix, seed)
    TYPE(cp_dbcsr_type)                      :: matrix
    INTEGER                                  :: seed

    CHARACTER(len=*), PARAMETER :: routineN = 'fill_randomly', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: acol, arow, handle, i, j
    REAL(dp), DIMENSION(:, :), POINTER       :: block
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed
    TYPE(cp_dbcsr_iterator)                  :: iter
    TYPE(rng_stream_type), POINTER           :: rng_stream

    NULLIFY(rng_stream)
    CALL timeset(routineN,handle)

    CALL cp_dbcsr_iterator_start(iter, matrix)
    DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
       CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block)
       ! a "unique" reproducible seed for each block
       initial_seed(1,:) = REAL(seed*arow,KIND=dp)
       initial_seed(2,:) = REAL(seed*acol,KIND=dp)
       initial_seed(3,:) = REAL(seed*SIZE(block),KIND=dp)
       CALL create_rng_stream(rng_stream, seed=initial_seed, name="pao")
       DO i=1, SIZE(block, 1)
          DO j=1, SIZE(block, 2)
             block(i,j) = next_random_number(rng_stream)
          ENDDO
       ENDDO
       CALL delete_rng_stream(rng_stream)
    END DO
    CALL cp_dbcsr_iterator_stop(iter)

    CALL timestop(handle)
  END SUBROUTINE fill_randomly

END MODULE pao_utils
