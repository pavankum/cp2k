/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"

#include <libxstream_begin.h>
#include <cstdio>
#if defined(_OPENMP)
# include <omp.h>
#endif
#include <libxstream_end.h>

#define LIBMICSMM_MAX_MATRIX_SIZE LIBXSTREAM_MAX(LIBXSTREAM_MAX( \
  LIBMICSMM_MAX_M * LIBMICSMM_MAX_K,  \
  LIBMICSMM_MAX_N * LIBMICSMM_MAX_K), \
  LIBMICSMM_MAX_M * LIBMICSMM_MAX_N)

#define LIBMICSMM_MAX_MNK LIBXSTREAM_MAX(LIBXSTREAM_MAX( \
  LIBMICSMM_MAX_M, LIBMICSMM_MAX_N), \
  LIBMICSMM_MAX_K)

#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
LIBXSTREAM_EXTERN_C void MKL_Simatcopy(const char, const char, size_t, size_t, float, float*, size_t, size_t);
LIBXSTREAM_EXTERN_C void MKL_Dimatcopy(const char, const char, size_t, size_t, double, double*, size_t, size_t);
#endif


namespace libmicsmm_transpose_private {

LIBXSTREAM_TARGET(mic) void mkl_imatcopy(size_t m, size_t n, float* matrix)
{
#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
  MKL_Simatcopy('R', 'T', m, n, 1.f, matrix, n, m);
#endif
}


LIBXSTREAM_TARGET(mic) void mkl_imatcopy(size_t m, size_t n, double* matrix)
{
#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
  MKL_Dimatcopy('R', 'T', m, n, 1.0, matrix, n, m);
#endif
}


template<typename T, typename U>
LIBXSTREAM_TARGET(mic) void kernel(const U *LIBXSTREAM_RESTRICT stack, U m, U n, T *LIBXSTREAM_RESTRICT matrix)
{
#if defined(LIBXSTREAM_DEBUG) && defined(_OPENMP)
  const double start = omp_get_wtime();
#endif
  const libxstream_argument* arg = 0;
  size_t stacksize = 0;
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_get_argument(stack, &arg));
  LIBXSTREAM_CHECK_CALL_RETURN(libxstream_get_shape(arg, &stacksize));

#if defined(_OPENMP)
# pragma omp parallel for
#endif
  for (U s = 0; s < stacksize; ++s) {
    T *const mat = matrix + stack[s];

#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
    mkl_imatcopy(static_cast<size_t>(m), static_cast<size_t>(n), mat);
#else
    LIBXSTREAM_ALIGNED(T tmp[LIBMICSMM_MAX_MATRIX_SIZE], LIBXSTREAM_MAX_SIMD);

# if defined(LIBMICSMM_USE_LOOPHINTS) && defined(__INTEL_COMPILER)
    // Best coice would be LIBMICSMM_MAX_MNK. However, the
    // trip count must be a compile-time constant with no further
    // macro invocation (actual trip can be larger anyways).
#   pragma loop_count min(1), max(LIBMICSMM_MAX_M), avg(23)
# endif
    for (U i = 0; i < m; ++i) {
# if defined(__INTEL_COMPILER)
#   if defined(LIBMICSMM_USE_LOOPHINTS)
#     pragma loop_count min(1), max(LIBMICSMM_MAX_N), avg(23)
#   endif
#     pragma simd
# elif defined(_OPENMP)
#     pragma omp simd
# endif
      for (U j = 0; j < n; ++j) {
        tmp[j*m+i] = mat[i*n+j];
      }
    }

# if defined(LIBMICSMM_USE_LOOPHINTS) && defined(__INTEL_COMPILER)
#   pragma loop_count min(1), max(LIBMICSMM_MAX_M), avg(23)
# endif
    for (U i = 0; i < m; ++i) {
# if defined(__INTEL_COMPILER)
#   if defined(LIBMICSMM_USE_LOOPHINTS)
#     pragma loop_count min(1), max(LIBMICSMM_MAX_N), avg(23)
#   endif
#     pragma simd aligned(tmp:1)
# elif defined(_OPENMP)
#     pragma omp simd aligned(tmp:1)
# endif
      for (U j = 0; j < n; ++j) {
        mat[i*n+j] = tmp[i*n+j];
      }
    }
#endif
  }

#if defined(LIBXSTREAM_DEBUG) && defined(_OPENMP)
  static double duration = 0, problemsize = 0;
  const double stop = omp_get_wtime();
  if (start < stop) {
#   pragma omp atomic
    duration += stop - start;
#   pragma omp atomic
    problemsize += 2ul * m * n * sizeof(T) * stacksize;
    LIBXSTREAM_PRINT_INFO("libsmm_acc_transpose: %.f GB/s", problemsize / (1E9 * duration));
  }
#endif
}


template<typename T, bool Complex, typename U>
int transpose(const U* stack, U offset, U nblocks, U m, U n, void* data, void* stream)
{
  LIBXSTREAM_PRINT_INFOCTX("type=%s offset=%i size=%i m=%i n=%i buffer=0x%lx stream=0x%lx",
    dbcsr_elem<T,COmplex>::name(), offset, nblks, m, n,
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(buffer)),
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(stream)));
  LIBXSTREAM_CHECK_CONDITION(
    stack && 0 <= offset && 0 <= nblocks
    && LIBMICSMM_MAX_MATRIX_SIZE >= (m * n)
    && 0 <= m && 0 <= n
    && data && stream);

  if (1 < m || 1 < n) {
    /*const*/ libxstream_function function = reinterpret_cast<libxstream_function>(kernel<T,U>);
    const size_t stacksize = nblocks;
    libxstream_argument* signature = 0;
    LIBXSTREAM_CHECK_CALL(libxstream_fn_create_signature(&signature, 4));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_input(signature, 0, stack + offset, libxstream_type2value<U>::value, 1, &stacksize));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_input(signature, 1, &m, libxstream_type2value<U>::value, 0, 0));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_input(signature, 2, &n, libxstream_type2value<U>::value, 0, 0));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_inout(signature, 3, data, libxstream_type2value<T>::value, 1, 0/*unknown*/));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_call(function, signature, static_cast<libxstream_stream*>(stream), LIBXSTREAM_CALL_DEFAULT));
    LIBXSTREAM_CHECK_CALL(libxstream_fn_destroy_signature(signature));
  }

  return LIBXSTREAM_ERROR_NONE;
}

} // namespace libmicsmm_transpose_private


LIBXSTREAM_EXTERN_C int libsmm_acc_transpose(void* trs_stack, int offset, int nblks, void* buffer, int datatype, int m, int n, void* stream)
{
  int result = LIBXSTREAM_ERROR_NONE;

#if defined(LIBMICSMM_USE_PRETRANSPOSE)
  const int *const stack = static_cast<const int*>(trs_stack);
  switch(static_cast<dbcsr_elem_type>(datatype)) {
    case DBCSR_ELEM_F32: {
      result = libmicsmm_transpose_private::transpose<float,false>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case DBCSR_ELEM_F64: {
      result = libmicsmm_transpose_private::transpose<double,false>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case DBCSR_ELEM_C32: {
      LIBXSTREAM_ASSERT(false/*TODO: not implemented yet*/);
      result = LIBXSTREAM_ERROR_CONDITION;
    } break;
    case DBCSR_ELEM_C64: {
      LIBXSTREAM_ASSERT(false/*TODO: not implemented yet*/);
      result = LIBXSTREAM_ERROR_CONDITION;
    } break;
    default:
      result = LIBXSTREAM_ERROR_CONDITION;
  }
#endif

  return result;
}

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
