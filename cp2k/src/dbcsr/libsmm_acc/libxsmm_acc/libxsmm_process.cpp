/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#if defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#include "libxsmm_acc.hpp"

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_begin.h>
#endif
#include <algorithm>
#include <cassert>
#include <cstdlib>
#if defined(LIBXSMM_ACC_OPENMP)
# include <omp.h>
#endif
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_end.h>
#endif

#if !defined(__LIBXSMM)
LIBXSMM_ACC_EXTERN_C LIBXSMM_ACC_TARGET(mic) void LIBXSMM_ACC_FSYMBOL(dgemm)(
  const char*, const char*, const int*, const int*, const int*,
  const double*, const double*, const int*, const double*, const int*,
  const double*, double*, const int*);
LIBXSMM_ACC_EXTERN_C LIBXSMM_ACC_TARGET(mic) void LIBXSMM_ACC_FSYMBOL(sgemm)(
  const char*, const char*, const int*, const int*, const int*,
  const float*, const float*, const int*, const float*, const int*,
  const float*, float*, const int*);
#endif

#define LIBXSMM_ACC_MAX_RESULT_SIZE (LIBXSMM_ACC_MAX_M * LIBXSMM_ACC_MAX_N)


namespace libxsmm_process_private {

#if defined(LIBXSMM_ACC_OPENMP) && defined(LIBXSMM_ACC_SYNCHRONIZATION) && (1 < (LIBXSMM_ACC_SYNCHRONIZATION))
LIBXSMM_ACC_TARGET(mic) class LIBXSMM_ACC_TARGET(mic) lock_type {
public:
  lock_type() {
    for (int i = 0; i < (LIBXSMM_ACC_SYNCHRONIZATION); ++i) omp_init_lock(m_lock + i);
  }
  
  ~lock_type() {
    for (int i = 0; i < (LIBXSMM_ACC_SYNCHRONIZATION); ++i) omp_destroy_lock(m_lock + i);
  }

public:
  void acquire(const void* address) {
    omp_set_lock(m_lock + LIBXSMM_ACC_HASH2(address, LIBXSMM_ACC_ALIGNED_MAX, LIBXSMM_ACC_SYNCHRONIZATION));
  }

  void release(const void* address) {
    omp_unset_lock(m_lock + LIBXSMM_ACC_HASH2(address, LIBXSMM_ACC_ALIGNED_MAX, LIBXSMM_ACC_SYNCHRONIZATION));
  }

private:
  omp_lock_t m_lock[LIBXSMM_ACC_SYNCHRONIZATION];
} lock;
#endif


template<typename T, typename U>
class LIBXSMM_ACC_TARGET(mic) smm_type {
public:
  typedef void (*xmm_function_type)(const T*, const T*, T*);
  typedef void (*smm_function_type)(U, U, U, U, const T*, const T*, T*, xmm_function_type);

public:
  smm_type(U def_mnk, U m, U n, U k)
#if defined(__LIBXSMM)
    : m_xmm_function(0 != def_mnk ? static_cast<xmm_function_type>(libxsmm_mm_dispatch<T>(m, n, k)) : static_cast<xmm_function_type>(0))
    , m_smm_function(0 != def_mnk
      ? ((LIBXSMM_MAX_MNK) >= (m * n * k) ? (0 != m_xmm_function ? smm_type::xmm : smm_type::imm) : smm_type::bmm)
      : ((LIBXSMM_MAX_MNK) >= (m * n * k) ? smm_type::amm : smm_type::bmm))
#else
    : m_xmm_function(0), m_smm_function(smm_type::bmm)
#endif
  {}

public:
  void zero_c(T *LIBXSMM_ACC_RESTRICT c, U size) const {
    LIBXSMM_ACC_ASSUME_ALIGNED(c, LIBXSMM_ACC_ALIGNED_STORES);
    LIBXSMM_ACC_PRAGMA_LOOP_COUNT(1, LIBXSMM_ACC_LOOP_MAX_M*LIBXSMM_ACC_LOOP_MAX_N, LIBXSMM_ACC_LOOP_AVG_M*LIBXSMM_ACC_LOOP_AVG_N)
    for (U i = 0; i < size; ++i) c[i] = 0;
  }

  void copy_c(const T *LIBXSMM_ACC_RESTRICT c, T *LIBXSMM_ACC_RESTRICT out, U m, U n, U ldc) const {
#if defined(LIBXSMM_ACC_OPENMP) && defined(LIBXSMM_ACC_SYNCHRONIZATION) && (0 < (LIBXSMM_ACC_SYNCHRONIZATION))
# if (1 == (LIBXSMM_ACC_SYNCHRONIZATION))
#   pragma omp critical(libxsmm_process)
# else
    lock.acquire(out);
# endif
#endif
    {
      LIBXSMM_ACC_ASSUME_ALIGNED(c, LIBXSMM_ACC_ALIGNED_STORES);
      for (U j = 0; j < n; ++j) {
        LIBXSMM_ACC_PRAGMA_LOOP_COUNT(1, LIBXSMM_ACC_LOOP_MAX_M, LIBXSMM_ACC_LOOP_AVG_M)
        for (U i = 0; i < m; ++i) {
          const T value = c[j*ldc+i];
#if defined(LIBXSMM_ACC_OPENMP) && (!defined(LIBXSMM_ACC_SYNCHRONIZATION) || (0 == (LIBXSMM_ACC_SYNCHRONIZATION)))
#         pragma omp atomic
#endif
          out[j*m+i] += value;
        }
      }
    }
#if defined(LIBXSMM_ACC_OPENMP) && defined(LIBXSMM_ACC_SYNCHRONIZATION) && (1 < (LIBXSMM_ACC_SYNCHRONIZATION))
    lock.release(out);
#endif
  }

  void operator()(U m, U n, U k, U ldc, const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c) const {
    LIBXSMM_ACC_ASSERT(m_smm_function);
    m_smm_function(m, n, k, ldc, a, b, c, m_xmm_function);
  }

private:
#if defined(__LIBXSMM)
  LIBXSMM_ACC_TARGET(mic) static void xmm(U, U, U, U, const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c, xmm_function_type xmm_function) {
    LIBXSMM_ACC_ASSERT(xmm_function);
    xmm_function(a, b, c);
  }

  LIBXSMM_ACC_TARGET(mic) static void imm(U m, U n, U k, U, const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c, xmm_function_type) {
    LIBXSMM_ACC_ASSERT((LIBXSMM_MAX_MNK) >= (m * n * k));
    libxsmm_imm(m, n, k, a, b, c);
  }

  LIBXSMM_ACC_TARGET(mic) static void amm(U m, U n, U k, U, const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c, xmm_function_type) {
    LIBXSMM_ACC_ASSERT((LIBXSMM_MAX_MNK) >= (m * n * k));
    const xmm_function_type xmm_function = libxsmm_mm_dispatch<T>(m, n, k);

    if (xmm_function) {
      xmm_function(a, b, c);
    }
    else {
      libxsmm_imm(m, n, k, a, b, c);
    }
  }
#endif

  LIBXSMM_ACC_TARGET(mic) static void bmm(U m, U n, U k, U ldc, const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c, xmm_function_type) {
    blasmm(m, n, k, ldc, a, b, c);
  }

  static void blasmm(U m, U n, U k, U ldc, const float* a, const float* b, float* c) {
    static float alpha = 1.f, beta = 1.f;
    static char trans = 'N';
    int im = static_cast<int>(m), in = static_cast<int>(n), ik = static_cast<int>(k), ildc = static_cast<int>(ldc);
    LIBXSMM_ACC_FSYMBOL(sgemm)(&trans, &trans, &im, &in, &ik, &alpha, const_cast<float*>(a), &im, const_cast<float*>(b), &ik, &beta, c, &ildc);
  }

  static void blasmm(U m, U n, U k, U ldc, const double* a, const double* b, double* c) {
    static double alpha = 1.0, beta = 1.0;
    static char trans = 'N';
    int im = static_cast<int>(m), in = static_cast<int>(n), ik = static_cast<int>(k), ildc = static_cast<int>(ldc);
    LIBXSMM_ACC_FSYMBOL(dgemm)(&trans, &trans, &im, &in, &ik, &alpha, const_cast<double*>(a), &im, const_cast<double*>(b), &ik, &beta, c, &ildc);
  }

private:
  mutable/*offload attribute*/ xmm_function_type m_xmm_function;
  smm_function_type m_smm_function;
};



template<size_t N, typename T, typename U>
LIBXSMM_ACC_TARGET(mic) void work_basic(const U *LIBXSMM_ACC_RESTRICT stack, size_t stacksize, const smm_type<T,U>& smm,
  const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c)
{
  const int nstacksize = static_cast<int>(stacksize * N);

#if defined(LIBXSMM_ACC_OPENMP)
# pragma omp parallel for LIBXSMM_ACC_SCHEDULE
#endif
  for (int s = 0; s < nstacksize; s += ((LIBXSMM_ACC_NLOCAL) * N)) {
    LIBXSMM_ACC_ALIGNED(T tmp[LIBXSMM_ACC_MAX_RESULT_SIZE], LIBXSMM_ACC_ALIGNED_MAX);
    const int end = s + std::min(static_cast<int>((LIBXSMM_ACC_NLOCAL) * N), nstacksize - s);
    U i = s, kc = stack[s+5], next = kc;
 
    do {
      const U m = stack[i+0], n = stack[i+1], ldc = LIBXSMM_ACC_ALIGN_VALUE(m, sizeof(T), LIBXSMM_ACC_ALIGNED_STORES);
      smm.zero_c(tmp, ldc * n);

      for (;;) {
        const U k = stack[i+2], ka = stack[i+3], kb = stack[i+4];
        smm(m, n, k, ldc, a + ka - 1, b + kb - 1, tmp);
        i += N;

        if (i < nstacksize) {
          next = stack[i+5];
          if (next != kc || end <= i) {
            break;
          }
        }
        else {
          break;
        }
      }

      smm.copy_c(tmp, c + kc - 1, m, n, ldc);
      kc = next;
    }
    while(i < end);
  }
}


template<size_t N, typename T, typename U>
LIBXSMM_ACC_TARGET(mic) void work_planned(const U *LIBXSMM_ACC_RESTRICT stack, U stacksize, const smm_type<T,U>& smm,
  const T *LIBXSMM_ACC_RESTRICT a, const T *LIBXSMM_ACC_RESTRICT b, T *LIBXSMM_ACC_RESTRICT c)
{
  const U nstacksize = static_cast<U>(N) * stacksize;
  const U plansize = 32768;
  U colspan[plansize];

  for (U s = 0; s < nstacksize;) {
    int size = 0;

    colspan[0] = s;
    for (; size < (plansize - 1) && s < nstacksize; s += N) {
      for (U kc0 = stack[s+5], kc1 = (((s + N) < nstacksize) ? stack[s+N+5] : (kc0 + 1)); kc0 == kc1; s += N, kc0 = kc1, kc1 = stack[s+N+5]);
      colspan[++size] = s + N;
    }

#if defined(LIBXSMM_ACC_OPENMP)
#   pragma omp parallel for LIBXSMM_ACC_SCHEDULE
#endif
    for (int i = 0; i < size; ++i) {
      const U j0 = colspan[i], j1 = colspan[i+1];
      const U kc = stack[j0+5];
      LIBXSMM_ACC_ASSERT(j1 <= nstacksize);
      LIBXSMM_ACC_ALIGNED(T tmp[LIBXSMM_ACC_MAX_RESULT_SIZE], LIBXSMM_ACC_ALIGNED_MAX);
      const U m = stack[i+0], n = stack[i+1], ldc = LIBXSMM_ACC_ALIGN_VALUE(m, sizeof(T), LIBXSMM_ACC_ALIGNED_STORES);
      smm.zero_c(tmp, ldc * n);

      for (U j = j0; j < j1; j += N) {
        LIBXSMM_ACC_ASSERT(j < nstacksize);
        LIBXSMM_ACC_ASSERT(kc == stack[j+5]);
        const U k = stack[j+2], ka = stack[j+3], kb = stack[j+4];
        smm(m, n, k, ldc, a + ka - 1, b + kb - 1, tmp);
      }

      smm.copy_c(tmp, c + kc - 1, m, n, ldc);
    }
  }
}


template<size_t N, typename T, typename U>
LIBXSMM_ACC_TARGET(mic) void context(const U* stack, const U* stacksize, const U* def_mnk, const U* max_m, const U* max_n, const U* max_k, const T* a, const T* b, T* c)
{
  const smm_type<T,U> smm(*def_mnk, *max_m, *max_n, *max_k);
#if defined(LIBXSMM_ACC_NLOCAL) && (0 < (LIBXSMM_ACC_NLOCAL))
  work_basic<LIBXSMM_ACC_NPARAMS,T,U>(stack, *stacksize, smm, a, b, c);
#elif defined(LIBXSMM_ACC_NLOCAL) && (0 == (LIBXSMM_ACC_NLOCAL))
  work_planned<LIBXSMM_ACC_NPARAMS,T,U>(stack, *stacksize, smm, a, b, c);
#else
  LIBXSMM_ACC_ASSERT(false/*TODO: not yet implemented.*/);
#endif
}


template<typename T, libxsmm_acc_bool_type Complex, typename U>
int process(const U* stack, U stacksize, U nparams, U def_mnk, U max_m, U max_n, U max_k, const void* a_data, const void* b_data, void* c_data, void* stream)
{
  LIBXSMM_ACC_CHECK_CONDITION(
    0 != stack && 0 <= stacksize && LIBXSMM_ACC_NPARAMS == nparams &&
    0 <= max_m && 0 <= max_n && 0 <= max_k &&
    0 != a_data && 0 != b_data && 0 != c_data);

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
  if (stream) {
    const size_t shape = stacksize;
    libxstream_argument* signature = 0;
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_signature(&signature));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 0,      stack, libxstream_map_to<U>::type(), 1, &shape));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 1, &stacksize, libxstream_map_to<U>::type(), 0, 0));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 2,   &def_mnk, libxstream_map_to<U>::type(), 0, 0));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 3,     &max_m, libxstream_map_to<U>::type(), 0, 0));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 4,     &max_n, libxstream_map_to<U>::type(), 0, 0));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 5,     &max_k, libxstream_map_to<U>::type(), 0, 0));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 6,     a_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 7,     b_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_inout(signature, 8,     c_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
    const libxstream_function libxsmm_process_function = reinterpret_cast<libxstream_function>(context<LIBXSMM_ACC_NPARAMS,T,U>);
    LIBXSMM_ACC_CHECK_CALL_ASSERT(libxstream_fn_call(libxsmm_process_function, signature, static_cast<libxstream_stream*>(stream), LIBXSTREAM_CALL_DEFAULT));
  }
  else
#endif
  {
    context<LIBXSMM_ACC_NPARAMS>(stack, &stacksize, &def_mnk, &max_m, &max_n, &max_k,
      static_cast<const T*>(a_data), static_cast<const T*>(b_data), static_cast<T*>(c_data));
  }

  return LIBXSMM_ACC_ERROR_NONE;
}

} // namespace libxsmm_process_private

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
// workaround for issue "cannot find address of function" (use unoptimized build or apply mic attribute globally)
const libxstream_function libxsmm_process_function = reinterpret_cast<libxstream_function>(libxsmm_process_private::context<LIBXSMM_ACC_NPARAMS,double,int>);
#endif


extern "C" int libsmm_acc_process(void* param_stack, int stacksize, int nparams, int datatype, void* a_data, void* b_data, void* c_data, int max_m, int max_n, int max_k, int def_mnk, void* stream)
{
#if defined(LIBXSMM_ACC_PRETRANSPOSE)
  LIBXSMM_ACC_ASSERT(false/*TODO: implement C = A * B which is assuming that B is pre-transposed (B^T).*/);
#endif
#if defined(__RECONFIGURE) && defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
  const int mflops = 0 != stream ? static_cast<int>(2E-6 * stacksize * max_m * max_n * max_k + 0.5) : LIBXSMM_ACC_ACCDRV_MIN_MFLOPS_PERSTACK;
  int result = (LIBXSMM_ACC_ACCDRV_MIN_MFLOPS_PERSTACK) <= mflops ? LIBXSMM_ACC_ERROR_NONE : LIBXSMM_ACC_NOT_SUPPORTED;
#else
  int result = LIBXSMM_ACC_ERROR_NONE;
#endif

  if (LIBXSMM_ACC_ERROR_NONE == result) {
    const int *const stack = static_cast<const int*>(param_stack);

    switch(static_cast<libxsmm_acc_elem_type>(datatype)) {
      case LIBXSMM_ACC_ELEM_F32: {
#if 0
        result = libxsmm_process_private::process<float,false>(stack, stacksize, nparams, def_mnk, max_m, max_n, max_k, a_data, b_data, c_data, stream);
#else
        result = LIBXSMM_ACC_NOT_SUPPORTED;
#endif
        result = LIBXSMM_ACC_ERROR_CONDITION;
      } break;
      case LIBXSMM_ACC_ELEM_F64: {
        result = libxsmm_process_private::process<double,false>(stack, stacksize, nparams, def_mnk, max_m, max_n, max_k, a_data, b_data, c_data, stream);
      } break;
      case LIBXSMM_ACC_ELEM_C32: {
        result = LIBXSMM_ACC_ERROR_CONDITION;
      } break;
      case LIBXSMM_ACC_ELEM_C64: {
        result = LIBXSMM_ACC_ERROR_CONDITION;
      } break;
      default:
        result = LIBXSMM_ACC_ERROR_CONDITION;
    }
  }

  return result;
}

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
