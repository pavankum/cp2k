!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2016  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Common framework for a linear parametrization of the potential.
!> \author Ole Schuett
! **************************************************************************************************
MODULE pao_param_linpot
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_interface,              ONLY: &
        cp_dbcsr_col_block_sizes, cp_dbcsr_create, cp_dbcsr_distribution, cp_dbcsr_frobenius_norm, &
        cp_dbcsr_get_block_p, cp_dbcsr_init, cp_dbcsr_iterator, cp_dbcsr_iterator_blocks_left, &
        cp_dbcsr_iterator_next_block, cp_dbcsr_iterator_start, cp_dbcsr_iterator_stop, &
        cp_dbcsr_p_type, cp_dbcsr_release, cp_dbcsr_reserve_diag_blocks, cp_dbcsr_type, &
        dbcsr_distribution_mp, dbcsr_mp_group
   USE kinds,                           ONLY: dp
   USE mathlib,                         ONLY: diamat_all
   USE message_passing,                 ONLY: mp_min,&
                                              mp_sum
   USE pao_input,                       ONLY: pao_fock_param,&
                                              pao_rotinv_param
   USE pao_linpot_full,                 ONLY: linpot_full_calc_term,&
                                              linpot_full_count_terms
   USE pao_linpot_rotinv,               ONLY: linpot_rotinv_calc_term,&
                                              linpot_rotinv_count_terms
   USE pao_param_fock,                  ONLY: pao_calc_U_block_fock,&
                                              pao_param_finalize_fock,&
                                              pao_param_init_fock
   USE pao_types,                       ONLY: pao_env_type
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: pao_param_init_linpot, pao_param_finalize_linpot, pao_calc_U_linpot
   PUBLIC :: pao_param_count_linpot, pao_param_initguess_linpot

CONTAINS

! **************************************************************************************************
!> \brief Initialize the linear potential parametrization
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_init_linpot(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_init_linpot'

      INTEGER                                            :: acol, arow, handle, iatom, ikind, kterm, &
                                                            n, natoms, nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_pri, col_blk_size, row_blk_size
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_V_term, block_V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)
      CALL pao_param_init_fock(pao, qs_env)

      CALL get_qs_env(qs_env, &
                      matrix_s=matrix_s, &
                      dft_control=dft_control, &
                      particle_set=particle_set, &
                      natom=natoms)

      IF (dft_control%nspins /= 1) CPABORT("open shell not yet implemented")

      ! figure out number of potential terms
      ALLOCATE (row_blk_size(natoms), col_blk_size(natoms))
      DO iatom = 1, natoms
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL pao_param_count_linpot(pao, qs_env, ikind, nterms)
         col_blk_size(iatom) = nterms
      ENDDO

      ! allocate matrix_V_terms
      blk_sizes_pri => cp_dbcsr_col_block_sizes(matrix_s(1)%matrix)
      row_blk_size = blk_sizes_pri**2
      CALL cp_dbcsr_init(pao%matrix_V_terms)
      CALL cp_dbcsr_create(pao%matrix_V_terms, &
                           template=matrix_s(1)%matrix, &
                           matrix_type="N", &
                           row_blk_size=row_blk_size, &
                           col_blk_size=col_blk_size)
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_V_terms)
      DEALLOCATE (row_blk_size, col_blk_size)

      ! calculate, normalize, and store potential terms as rows of block_V_terms
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_V_terms)
         iatom = arow
         n = blk_sizes_pri(iatom)
         CPASSERT(n*n == SIZE(block_V_terms, 1))
         nterms = SIZE(block_V_terms, 2)
         DO kterm = 1, nterms
            block_V_term(1:n, 1:n) => block_V_terms(:, kterm) ! map column into matrix
            CALL linpot_calc_term(pao, qs_env, iatom, kterm, block_V_term)
         ENDDO
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      CALL pao_param_linpot_regularizer(pao)

      IF (pao%precondition) &
         CALL pao_param_linpot_preconditioner(pao)

      CALL timestop(handle)
   END SUBROUTINE pao_param_init_linpot

! **************************************************************************************************
!> \brief Builds the regularization metric matrix_R
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_linpot_regularizer(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_linpot_regularizer'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, &
                                                            nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_nterms
      LOGICAL                                            :: found
      REAL(dp)                                           :: v, w
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_R, V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Building linpot regularizer"

      blk_sizes_nterms => cp_dbcsr_col_block_sizes(pao%matrix_V_terms)

      ! build regularization metric
      CALL cp_dbcsr_init(pao%matrix_R)
      CALL cp_dbcsr_create(pao%matrix_R, &
                           template=pao%matrix_V_terms, &
                           matrix_type="N", &
                           row_blk_size=blk_sizes_nterms, &
                           col_blk_size=blk_sizes_nterms, &
                           name="PAO matrix_R")
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_R)

      ! fill matrix_R
      CALL cp_dbcsr_iterator_start(iter, pao%matrix_R)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_R)
         iatom = arow; CPASSERT(arow == acol)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=V_terms, found=found)
         CPASSERT(ASSOCIATED(V_terms))
         nterms = SIZE(V_terms, 2)

         ! build overlap matrix
         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(V_terms), V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         block_R = 0.0_dp
         DO k = 1, nterms
            v = pao%linpot_regu_delta/S_evals(k)
            w = pao%linpot_regu_strength*MIN(1.0_dp, ABS(v))
            DO i = 1, nterms
            DO j = 1, nterms
               block_R(i, j) = block_R(i, j)+w*S_evecs(i, k)*S_evecs(j, k)
            ENDDO
            ENDDO
         ENDDO

         ! clean up
         DEALLOCATE (S, S_evals, S_evecs)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE pao_param_linpot_regularizer

! **************************************************************************************************
!> \brief Builds the preconditioner matrix_precon and matrix_precon_inv
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_linpot_preconditioner(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_linpot_preconditioner'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, &
                                                            nterms
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes_nterms
      LOGICAL                                            :: found
      REAL(dp)                                           :: eval_capped
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_precon, block_precon_inv, &
                                                            block_V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Building linpot preconditioner"

      blk_sizes_nterms => cp_dbcsr_col_block_sizes(pao%matrix_V_terms)

      CALL cp_dbcsr_init(pao%matrix_precon)
      CALL cp_dbcsr_create(pao%matrix_precon, &
                           template=pao%matrix_V_terms, &
                           matrix_type="N", &
                           row_blk_size=blk_sizes_nterms, &
                           col_blk_size=blk_sizes_nterms, &
                           name="PAO matrix_precon")
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_precon)

      CALL cp_dbcsr_init(pao%matrix_precon_inv)
      CALL cp_dbcsr_create(pao%matrix_precon_inv, template=pao%matrix_precon, name="PAO matrix_precon_inv")
      CALL cp_dbcsr_reserve_diag_blocks(pao%matrix_precon_inv)

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_V_terms)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_V_terms)
         iatom = arow; CPASSERT(arow == acol)

         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_precon, row=iatom, col=iatom, block=block_precon, found=found)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_precon_inv, row=iatom, col=iatom, block=block_precon_inv, found=found)
         CPASSERT(ASSOCIATED(block_precon))
         CPASSERT(ASSOCIATED(block_precon_inv))

         nterms = SIZE(block_V_terms, 2)
         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(block_V_terms), block_V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         ! construct 1/Sqrt(S) and Sqrt(S)
         block_precon = 0.0_dp
         block_precon_inv = 0.0_dp
         DO k = 1, nterms
            eval_capped = MAX(pao%linpot_precon_delta, S_evals(k)) ! too small eigenvalues are hurtful
            DO i = 1, nterms
            DO j = 1, nterms
               block_precon(i, j) = block_precon(i, j)+S_evecs(i, k)*S_evecs(j, k)/SQRT(eval_capped)
               block_precon_inv(i, j) = block_precon_inv(i, j)+S_evecs(i, k)*S_evecs(j, k)*SQRT(eval_capped)
            ENDDO
            ENDDO
         ENDDO

         DEALLOCATE (S, S_evecs, S_evals)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      CALL timestop(handle)
   END SUBROUTINE pao_param_linpot_preconditioner

! **************************************************************************************************
!> \brief Finalize the linear potential parametrization
!> \param pao ...
! **************************************************************************************************
   SUBROUTINE pao_param_finalize_linpot(pao)
      TYPE(pao_env_type), POINTER                        :: pao

      CALL pao_param_finalize_fock(pao)
      CALL cp_dbcsr_release(pao%matrix_V_terms)
      CALL cp_dbcsr_release(pao%matrix_R)

      IF (pao%precondition) THEN
         CALL cp_dbcsr_release(pao%matrix_precon)
         CALL cp_dbcsr_release(pao%matrix_precon_inv)
      ENDIF

   END SUBROUTINE pao_param_finalize_linpot

! **************************************************************************************************
!> \brief Returns the number of potential terms for given atomic kind
!> \param pao ...
!> \param qs_env ...
!> \param ikind ...
!> \param nparams ...
! **************************************************************************************************
   SUBROUTINE pao_param_count_linpot(pao, qs_env, ikind, nparams)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: ikind
      INTEGER, INTENT(OUT)                               :: nparams

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         CALL linpot_full_count_terms(qs_env, ikind, nterms=nparams)
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_count_terms(qs_env, ikind, nterms=nparams)
      CASE DEFAULT
         CPABORT("unkown parameterization")
      END SELECT

   END SUBROUTINE pao_param_count_linpot

! **************************************************************************************************
!> \brief Calculate new matrix U and optinally its gradient G
!> \param pao ...
!> \param penalty ...
!> \param matrix_M ...
!> \param matrix_G ...
! **************************************************************************************************
   SUBROUTINE pao_calc_U_linpot(pao, penalty, matrix_M, matrix_G)
      TYPE(pao_env_type), POINTER                        :: pao
      REAL(dp), INTENT(INOUT)                            :: penalty
      TYPE(cp_dbcsr_type), OPTIONAL                      :: matrix_M, matrix_G

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pao_calc_U_linpot'

      INTEGER                                            :: acol, arow, group, handle, iatom, n
      LOGICAL                                            :: found
      REAL(dp), DIMENSION(:), POINTER                    :: vec_M2, vec_V
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_G, block_M1, block_M2, block_R, &
                                                            block_U, block_V, block_V_terms, &
                                                            block_X
      REAL(KIND=dp)                                      :: min_gap, regu_energy
      TYPE(cp_dbcsr_iterator)                            :: iter

      CALL timeset(routineN, handle)

      min_gap = HUGE(dp)
      regu_energy = 0.0_dp
      group = dbcsr_mp_group(dbcsr_distribution_mp(cp_dbcsr_distribution(pao%matrix_U)))

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_R, row=iatom, col=iatom, block=block_R, found=found)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_U, row=iatom, col=iatom, block=block_U, found=found)
         CPASSERT(ASSOCIATED(block_R) .AND. ASSOCIATED(block_U))
         n = SIZE(block_U, 1)

         ! calculate potential V
         ALLOCATE (vec_V(n*n))
         vec_V(:) = 0.0_dp
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=block_V_terms, found=found)
         IF (SIZE(block_V_terms) > 0) & ! protect against corner-case of zero pao parameters
            vec_V = MATMUL(block_V_terms, block_X(:, 1))
         block_V(1:n, 1:n) => vec_V(:) ! map vector into matrix

         ! symmetrize
         IF (MAXVAL(ABS(block_V-TRANSPOSE(block_V))/MAX(1.0_dp, MAXVAL(ABS(block_V)))) > 1e-12) &
            CPABORT("block_V not symmetric")
         block_V = 0.5_dp*(block_V+TRANSPOSE(block_V)) ! symmetrize exactly

         ! regularization energy
         regu_energy = regu_energy+DOT_PRODUCT(block_X(:, 1), MATMUL(block_R, block_X(:, 1)))

         IF (.NOT. PRESENT(matrix_G)) THEN
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U)

         ELSE ! TURNING POINT (if calc grad) ------------------------------------
            CPASSERT(PRESENT(matrix_M))
            CALL cp_dbcsr_get_block_p(matrix=pao%matrix_G, row=iatom, col=iatom, block=block_G, found=found)
            CPASSERT(ASSOCIATED(block_G))
            CALL cp_dbcsr_get_block_p(matrix=matrix_M, row=iatom, col=iatom, block=block_M1, found=found)
            ! don't check ASSOCIATED(block_M1), it might have been filtered out.

            ALLOCATE (vec_M2(n*n))
            block_M2(1:n, 1:n) => vec_M2(:) ! map vector into matrix
            CALL pao_calc_U_block_fock(pao, iatom, penalty, min_gap, block_V, block_U, block_M1, block_M2)
            IF (MAXVAL(ABS(block_M2-TRANSPOSE(block_M2))) > 1e-14_dp) CPABORT("Expect symmetric matrix")
            IF (SIZE(block_V_terms) > 0) THEN ! protect against corner-case of zero pao parameters
               block_G(:, 1) = MATMUL(vec_M2, block_V_terms)
               block_G = block_G+2.0_dp*MATMUL(block_R, block_X) ! regularization gradient
            ENDIF
            DEALLOCATE (vec_M2)
         ENDIF

         DEALLOCATE (vec_V)

      END DO
      CALL cp_dbcsr_iterator_stop(iter)

      ! sum penalty energies across ranks
      CALL mp_sum(penalty, group)
      CALL mp_sum(regu_energy, group)

      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| linpot regularization energy:", regu_energy
      penalty = penalty+regu_energy

      ! print smallest homo-lumo gap encountered by fock-layer
      CALL mp_min(min_gap, group)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| min_gap:", min_gap

      CALL timestop(handle)

   END SUBROUTINE pao_calc_U_linpot

! **************************************************************************************************
!> \brief Internal routine, calculates terms in potential parametrization
!> \param pao ...
!> \param qs_env ...
!> \param iatom ...
!> \param kterm ...
!> \param block_V ...
! **************************************************************************************************
   SUBROUTINE linpot_calc_term(pao, qs_env, iatom, kterm, block_V)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(IN)                                :: iatom, kterm
      REAL(dp), DIMENSION(:, :), INTENT(OUT)             :: block_V

      SELECT CASE (pao%parameterization)
      CASE (pao_fock_param)
         CALL linpot_full_calc_term(kterm, block_V)
      CASE (pao_rotinv_param)
         CALL linpot_rotinv_calc_term(qs_env, iatom, kterm, block_V)
      CASE DEFAULT
         CPABORT("unkown parameterization")
      END SELECT

   END SUBROUTINE linpot_calc_term

! **************************************************************************************************
!> \brief Calculate intial guess for matrix_X
!> \param pao ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE pao_param_initguess_linpot(pao, qs_env)
      TYPE(pao_env_type), POINTER                        :: pao
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'pao_param_initguess_linpot'

      INTEGER                                            :: acol, arow, handle, i, iatom, j, k, n, &
                                                            nterms
      LOGICAL                                            :: found
      REAL(dp)                                           :: norm, w
      REAL(dp), ALLOCATABLE, DIMENSION(:)                :: S_evals
      REAL(dp), ALLOCATABLE, DIMENSION(:, :)             :: S, S_evecs, S_inv
      REAL(dp), DIMENSION(:), POINTER                    :: V_guess_vec
      REAL(dp), DIMENSION(:, :), POINTER                 :: block_H0, block_ks, block_X, V_guess, &
                                                            V_terms
      TYPE(cp_dbcsr_iterator)                            :: iter
      TYPE(cp_dbcsr_p_type), DIMENSION(:), POINTER       :: matrix_ks

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)

      CALL cp_dbcsr_iterator_start(iter, pao%matrix_X)
      DO WHILE (cp_dbcsr_iterator_blocks_left(iter))
         CALL cp_dbcsr_iterator_next_block(iter, arow, acol, block_X)
         iatom = arow; CPASSERT(arow == acol)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_V_terms, row=iatom, col=iatom, block=V_terms, found=found)
         CPASSERT(ASSOCIATED(V_terms))
         nterms = SIZE(V_terms, 2)
         IF (nterms == 0) CYCLE ! protect against corner-case of zero pao parameters

         ! retrieve block_ks and block_H0
         CALL cp_dbcsr_get_block_p(matrix=matrix_ks(1)%matrix, row=iatom, col=iatom, block=block_ks, found=found)
         CALL cp_dbcsr_get_block_p(matrix=pao%matrix_H0, row=iatom, col=iatom, block=block_H0, found=found)
         CPASSERT(ASSOCIATED(block_H0) .AND. ASSOCIATED(block_ks))

         ! use diagonal block of Kohn-Sham matrix as initial guess
         n = SIZE(block_ks, 1)
         ALLOCATE (V_guess_vec(n*n))
         V_guess(1:n, 1:n) => V_guess_vec
         V_guess = block_ks-block_H0

         ! build overlap matrix
         ALLOCATE (S(nterms, nterms))
         S(:, :) = MATMUL(TRANSPOSE(V_terms), V_terms)

         ! diagonalize S
         ALLOCATE (S_evals(nterms), S_evecs(nterms, nterms))
         S_evecs(:, :) = S
         CALL diamat_all(S_evecs, S_evals)

         ! calculate Tikhonov regularized inverse
         ALLOCATE (S_inv(nterms, nterms))
         S_inv(:, :) = 0.0_dp
         DO k = 1, nterms
            w = S_evals(k)/(S_evals(k)**2+pao%linpot_init_delta)
            DO i = 1, nterms
            DO j = 1, nterms
               S_inv(i, j) = S_inv(i, j)+w*S_evecs(i, k)*S_evecs(j, k)
            ENDDO
            ENDDO
         ENDDO

         ! perform fit
         block_X(:, 1) = MATMUL(MATMUL(S_inv, TRANSPOSE(V_terms)), V_guess_vec)

         ! clean up
         DEALLOCATE (V_guess_vec, S, S_evecs, S_evals, S_inv)
      ENDDO
      CALL cp_dbcsr_iterator_stop(iter)

      norm = cp_dbcsr_frobenius_norm(pao%matrix_X)
      IF (pao%iw > 0) WRITE (pao%iw, *) "PAO| Calculated initial guess, norm(X):", norm

      CALL timestop(handle)
   END SUBROUTINE pao_param_initguess_linpot

END MODULE pao_param_linpot
