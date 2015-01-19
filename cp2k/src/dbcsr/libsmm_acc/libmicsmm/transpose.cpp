/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"
#include <libxstream.hpp>

#if defined(LIBXSTREAM_OFFLOAD)
# pragma offload_attribute(push,target(mic))
#endif

#if defined(_OPENMP)
# include <omp.h>
#endif

#if defined(LIBXSTREAM_OFFLOAD)
# pragma offload_attribute(pop)
#endif

#define LIBMICSMM_MAX_MATRIX_SIZE LIBXSTREAM_MAX(LIBXSTREAM_MAX( \
  LIBMICSMM_MAX_M * LIBMICSMM_MAX_K,  \
  LIBMICSMM_MAX_N * LIBMICSMM_MAX_K), \
  LIBMICSMM_MAX_M * LIBMICSMM_MAX_N)

#define LIBMICSMM_MAX_MNK LIBXSTREAM_MAX(LIBXSTREAM_MAX( \
  LIBMICSMM_MAX_M, LIBMICSMM_MAX_N), \
  LIBMICSMM_MAX_K)

#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
extern "C" LIBXSTREAM_EXPORT void MKL_Simatcopy(const char, const char, size_t, size_t, float, float*, size_t, size_t);
extern "C" LIBXSTREAM_EXPORT void MKL_Dimatcopy(const char, const char, size_t, size_t, double, double*, size_t, size_t);
#endif


namespace libmicsmm_transpose_private {

LIBXSTREAM_EXPORT void mkl_imatcopy(size_t m, size_t n, float* matrix)
{
#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
  MKL_Simatcopy('R', 'T', m, n, 1.f, matrix, n, m);
#endif
}


LIBXSTREAM_EXPORT void mkl_imatcopy(size_t m, size_t n, double* matrix)
{
#if defined(LIBMICSMM_USE_MKLTRANS) && defined(__MKL)
  MKL_Dimatcopy('R', 'T', m, n, 1.0, matrix, n, m);
#endif
}


template<typename T, typename U>
LIBXSTREAM_EXPORT void kernel(const U *LIBXSTREAM_RESTRICT stack, U offset, U nblocks, U m, U n, T *LIBXSTREAM_RESTRICT matrix)
{
#if defined(LIBXSTREAM_TEST) && (0 != (2*LIBXSTREAM_TEST+1)/2) && defined(_OPENMP)
  const double start = omp_get_wtime();
#endif
  const U *const offsets = stack + offset;

#if defined(_OPENMP)
# pragma omp parallel for
#endif
  for (U s = 0; s < nblocks; ++s) {
    T *const mat = matrix + offsets[s];

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

#if defined(LIBXSTREAM_TEST) && (0 != (2*LIBXSTREAM_TEST+1)/2) && defined(_OPENMP)
  static double duration = 0, problemsize = 0;
  const double stop = omp_get_wtime();
  if (start < stop) {
#   pragma omp atomic
    duration += stop - start;
#   pragma omp atomic
    problemsize += 2ul * m * n * sizeof(T) * nblocks;
    fprintf(stderr, "PRF libsmm_acc_transpose: %.f GB/s\n", problemsize / (1E9 * duration));
  }
#endif
}


template<typename T, typename U>
int transpose(const U* stack, U offset, U nblocks, U m, U n, void* data, void* stream)
{
  LIBXSTREAM_CHECK_CONDITION(
    stack && 0 <= offset && 0 <= nblocks
    && LIBMICSMM_MAX_MATRIX_SIZE >= (m * n)
    && 0 <= m && 0 <= n
    && data && stream);
  const int result = static_cast<libxstream_stream*>(stream)->reset();

  if (1 < m || 1 < n) {
    LIBXSTREAM_OFFLOAD_BEGIN(stream, offset, nblocks, m, n, stack, data)
    {
      const U offset = val<const U,0>();
      const U nblocks = val<const U,1>();
      const U m = val<const U,2>();
      const U n = val<const U,3>();
      const U *const stack = ptr<const U,4>();
      T* buffer = ptr<T,5>();

#if defined(LIBXSTREAM_OFFLOAD)
      if (0 <= LIBXSTREAM_OFFLOAD_DEVICE) {
        if (LIBXSTREAM_OFFLOAD_READY) {
#         pragma offload LIBXSTREAM_OFFLOAD_TARGET_SIGNAL \
            in(offset, nblocks, m, n) \
            in(stack: length(0) alloc_if(false) free_if(false)) \
            inout(buffer: length(0) alloc_if(false) free_if(false))
          kernel(stack, offset, nblocks, m, n, buffer);
        }
        else {
#         pragma offload LIBXSTREAM_OFFLOAD_TARGET_WAIT \
            in(offset, nblocks, m, n) \
            in(stack: length(0) alloc_if(false) free_if(false)) \
            inout(buffer: length(0) alloc_if(false) free_if(false))
          kernel(stack, offset, nblocks, m, n, buffer);
        }
      }
      else
#endif
      {
        kernel(stack, offset, nblocks, m, n, buffer);
      }
    }
    LIBXSTREAM_OFFLOAD_END(false)
  }

  return result;
}

} // namespace libmicsmm_transpose_private


extern "C" int libsmm_acc_transpose(void* trs_stack, int offset, int nblks, void* buffer, int datatype, int m, int n, void* stream)
{
#if defined(LIBXSTREAM_DEBUG)
  fprintf(stderr, "DBG libsmm_acc_transpose: offset=%i nblocks=%i m=%i n=%i buffer=0x%lx stream=0x%lx\n", offset, nblks, m, n,
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(buffer)),
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(stream)));
#endif
  const int *const stack = static_cast<const int*>(trs_stack);
  int result = LIBXSTREAM_ERROR_NONE;

  switch(static_cast<dbcsr_elem_type>(datatype)) {
    case DBCSR_ELEM_F32: {
      result = libmicsmm_transpose_private::transpose<float>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case DBCSR_ELEM_F64: {
      result = libmicsmm_transpose_private::transpose<double>(stack, offset, nblks, m, n, buffer, stream);
    } break;
    case DBCSR_ELEM_C32: {
      result = LIBXSTREAM_ERROR_CONDITION;
    } break;
    case DBCSR_ELEM_C64: {
      result = LIBXSTREAM_ERROR_CONDITION;
    } break;
    default:
      result = LIBXSTREAM_ERROR_CONDITION;
  }

  return result;
}

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
