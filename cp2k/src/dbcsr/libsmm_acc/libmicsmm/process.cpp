/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"
#include <vector>

#include <libxstream_begin.h>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#if defined(_OPENMP)
# include <omp.h>
#endif
#include <libxstream_end.h>

#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM)
# include <libxsmm.h>
# if (0 != LIBXSMM_ROW_MAJOR)
#   error Please compile LIBXSMM using "make ROW_MAJOR=0 ..."!
# endif
#else
LIBXSTREAM_EXTERN_C LIBXSTREAM_TARGET(mic) void LIBXSTREAM_FSYMBOL(dgemm)(
  const char*, const char*, const int*, const int*, const int*,
  const double*, const double*, const int*, const double*, const int*,
  const double*, double*, const int*);
LIBXSTREAM_EXTERN_C LIBXSTREAM_TARGET(mic) void LIBXSTREAM_FSYMBOL(sgemm)(
  const char*, const char*, const int*, const int*, const int*,
  const float*, const float*, const int*, const float*, const int*,
  const float*, float*, const int*);
#endif

#if defined(LIBMICSMM_STACKSIZE) && (0 < (LIBMICSMM_STACKSIZE))
LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int*, void*);
LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(const int*, void*);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers);
#endif

#define LIBMICSMM_MAX_RESULT_SIZE (LIBMICSMM_MAX_M * LIBMICSMM_MAX_N)
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM) && (0 < (LIBXSMM_ALIGNED_STORES))
# define LIBMICSMM_ALIGNMENT LIBXSMM_ALIGNED_STORES
#else
# define LIBMICSMM_ALIGNMENT LIBXSTREAM_MAX_SIMD
#endif


LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int* driver, void* error)
{
  LIBXSTREAM_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(driver, error);
#if defined(LIBMICSMM_STACKSIZE) && (0 < (LIBMICSMM_STACKSIZE))
  const int stacksize = LIBMICSMM_STACKSIZE;
  LIBXSTREAM_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(&stacksize, error);
  LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams) = LIBMICSMM_POSTERIOR_STREAMS;
  LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers) = LIBMICSMM_POSTERIOR_BUFFERS;
  LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams) = LIBMICSMM_PRIORITY_STREAMS;
  LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers) = LIBMICSMM_PRIORITY_BUFFERS;
#endif
}


namespace libmicsmm_process_private {

#if defined(_OPENMP) && defined(LIBMICSMM_SYNCHRONIZATION) && (1 < (LIBMICSMM_SYNCHRONIZATION))
LIBXSTREAM_TARGET(mic) class LIBXSTREAM_TARGET(mic) lock_type {
public:
  lock_type() {
    for (size_t i = 0; i < (LIBMICSMM_SYNCHRONIZATION); ++i) omp_init_lock(m_lock + i);
  }
  
  ~lock_type() {
    for (size_t i = 0; i < (LIBMICSMM_SYNCHRONIZATION); ++i) omp_destroy_lock(m_lock + i);
  }

public:
  void acquire(const void* id) {
    omp_set_lock(m_lock + reinterpret_cast<uintptr_t>(id) % LIBMICSMM_SYNCHRONIZATION);
  }

  void release(const void* id) {
    omp_unset_lock(m_lock + reinterpret_cast<uintptr_t>(id) % LIBMICSMM_SYNCHRONIZATION);
  }

private:
  omp_lock_t m_lock[LIBMICSMM_SYNCHRONIZATION];
} lock;
#endif // LIBXSTREAM_TARGET(mic)


template<typename T, typename U>
class LIBXSTREAM_TARGET(mic) smm_type {
public:
  typedef void (*xmm_function_type)(const T*, const T*, T*);
  typedef void (*smm_function_type)(U, U, U, U, const T*, const T*, T*, xmm_function_type);

public:
  smm_type(U m, U n, U k)
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM)
    : m_xmm_function(libxsmm_mm_dispatch<T>(m, n, k))
    , m_smm_function((LIBXSMM_MAX_MNK) >= (m * n * k)
      ? (0 != m_xmm_function ? smm_type::xmm : smm_type::imm)
      : smm_type::bmm)
#else
    : m_xmm_function(0), m_smm_function(smm_type::bmm)
#endif
    , m_m(m), m_n(n), m_k(k)
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM) && (0 < (LIBXSMM_ALIGNED_STORES))
    , m_ldc(LIBXSTREAM_ALIGN_VALUE(U, T, m, LIBXSMM_ALIGNED_STORES))
#else
    , m_ldc(m)
#endif
  {}

public:
  void zero_c(T *LIBXSTREAM_RESTRICT c) const {
    const U size = m_n * m_ldc;
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM) && (0 < (LIBXSMM_ALIGNED_STORES))
    LIBXSMM_ASSUME_ALIGNED(c, LIBMICSMM_ALIGNMENT);
#endif
    LIBXSTREAM_PRAGMA_LOOP_COUNT(1, LIBMICSMM_MAX_RESULT_SIZE, 23*23)
    for (U i = 0; i < size; ++i) c[i] = 0;
  }

  void copy_c(const T *LIBXSTREAM_RESTRICT c, T *LIBXSTREAM_RESTRICT out) const {
#if defined(_OPENMP) && defined(LIBMICSMM_SYNCHRONIZATION) && (0 < (LIBMICSMM_SYNCHRONIZATION))
# if (1 == (LIBMICSMM_SYNCHRONIZATION))
#   pragma omp critical(libmicsmm_process)
# else
    lock.acquire(c);
# endif
#endif
    {
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM) && (0 < (LIBXSMM_ALIGNED_STORES))
      LIBXSMM_ASSUME_ALIGNED(c, LIBMICSMM_ALIGNMENT);
#endif
      for (U j = 0; j < m_n; ++j) {
        LIBXSTREAM_PRAGMA_LOOP_COUNT(1, LIBMICSMM_MAX_M, 23)
        for (U i = 0; i < m_m; ++i) {
          const T value = c[j*m_ldc+i];
#if defined(_OPENMP) && (!defined(LIBMICSMM_SYNCHRONIZATION) || (0 == (LIBMICSMM_SYNCHRONIZATION)))
#         pragma omp atomic
#endif
          out[j*m_m+i] += value;
        }
      }
    }
#if defined(_OPENMP) && defined(LIBMICSMM_SYNCHRONIZATION) && (1 < (LIBMICSMM_SYNCHRONIZATION))
    lock.release(c);
#endif
  }

  void operator()(const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c) const {
    LIBXSTREAM_ASSERT(m_smm_function);
    m_smm_function(m_m, m_n, m_k, m_ldc, a, b, c, m_xmm_function);
  }

private:
  LIBXSTREAM_TARGET(mic) static void bmm(U m, U n, U k, U ldc, const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c, xmm_function_type) {
#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM)
    LIBXSTREAM_ASSERT((LIBXSMM_MAX_MNK) < (m * n * k));
    libxsmm_blasmm(m, n, k, a, b, c);
    libxstream_use_sink(&ldc);
#else
    blasmm(m, n, k, ldc, a, b, c);
#endif
  }

#if defined(LIBMICSMM_LIBXSMM) && defined(__LIBXSMM)
  LIBXSTREAM_TARGET(mic) static void xmm(U, U, U, U, const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c, xmm_function_type xmm_function) {
    LIBXSTREAM_ASSERT(xmm_function);
    xmm_function(a, b, c);
  }

  LIBXSTREAM_TARGET(mic) static void imm(U m, U n, U k, U, const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c, xmm_function_type) {
    LIBXSTREAM_ASSERT((LIBXSMM_MAX_MNK) >= (m * n * k));
    libxsmm_imm(m, n, k, a, b, c);
  }

#else /*no LIBXSMM*/

  static void blasmm(U m, U n, U k, U ldc, const float* a, const float* b, float* c) {
    static float alpha = 1.f, beta = 1.f;
    static char trans = 'N';
    int im = static_cast<int>(m), in = static_cast<int>(n), ik = static_cast<int>(k), ildc = static_cast<int>(ldc);
    LIBXSTREAM_FSYMBOL(sgemm)(&trans, &trans, &im, &in, &ik, &alpha, const_cast<float*>(a), &im, const_cast<float*>(b), &ik, &beta, c, &ildc);
  }

  static void blasmm(U m, U n, U k, U ldc, const double* a, const double* b, double* c) {
    static double alpha = 1.0, beta = 1.0;
    static char trans = 'N';
    int im = static_cast<int>(m), in = static_cast<int>(n), ik = static_cast<int>(k), ildc = static_cast<int>(ldc);
    LIBXSTREAM_FSYMBOL(dgemm)(&trans, &trans, &im, &in, &ik, &alpha, const_cast<double*>(a), &im, const_cast<double*>(b), &ik, &beta, c, &ildc);
  }
#endif

private:
  mutable/*offload attribute*/ xmm_function_type m_xmm_function;
  smm_function_type m_smm_function;
  U m_m, m_n, m_k, m_ldc;
};


template<size_t N, typename T, typename U>
LIBXSTREAM_TARGET(mic) void kernel(const U *LIBXSTREAM_RESTRICT stack, LIBXSTREAM_INVAL(U) max_m, LIBXSTREAM_INVAL(U) max_n, LIBXSTREAM_INVAL(U) max_k,
  const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c)
{
  size_t stacksize = 0;
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_get_shape(0/*current context*/, 0/*stack*/, &stacksize));
  LIBXSTREAM_PRINT_INFO("libsmm_acc_process (" LIBXSTREAM_DEVICE_NAME "): stacksize=%lu max_m=%i max_n=%i max_k=%i", static_cast<unsigned long>(stacksize),
    LIBXSTREAM_GETVAL(max_m), LIBXSTREAM_GETVAL(max_n), LIBXSTREAM_GETVAL(max_k));

#if defined(LIBXSTREAM_PRINT) && defined(_OPENMP)
  const double start = omp_get_wtime();
#endif
#if defined(LIBMICSMM_THREADPRIVATE) && defined(_OPENMP)
# if 1 == (LIBMICSMM_THREADPRIVATE) // native OpenMP TLS
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
# pragma omp threadprivate(tmp)
#else
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static LIBXSTREAM_TLS T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
# endif
#else // without OpenMP nothing needs to be thread-local due to a single-threaded program
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
#endif
  const smm_type<T,U> smm(LIBXSTREAM_GETVAL(max_m), LIBXSTREAM_GETVAL(max_n), LIBXSTREAM_GETVAL(max_k));
  const int nstacksize = static_cast<int>(stacksize * N);

#if defined(_OPENMP)
# pragma omp parallel for schedule(LIBMICSMM_SCHEDULE)
#endif
  for (int n = 0; n < nstacksize; n += ((LIBMICSMM_NLOCAL) * N)) {
#if !defined(LIBMICSMM_THREADPRIVATE)
    LIBXSTREAM_ALIGNED(T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
#endif
    const int end = n + std::min(static_cast<int>((LIBMICSMM_NLOCAL) * N), nstacksize - n);
    U i = n, kc = stack[n+5], next = kc;
 
    do {
      smm.zero_c(tmp);

      for (;;) {
        const U ka = stack[i+3], kb = stack[i+4];
        smm(a + ka - 1, b + kb - 1, tmp);
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

      smm.copy_c(tmp, c + kc - 1);
      kc = next;
    }
    while(i < end);
  }

#if defined(LIBXSTREAM_PRINT) && defined(_OPENMP)
  static double duration = 0, flops = 0;
  const double stop = omp_get_wtime();
  if (start < stop) {
#   pragma omp atomic
    duration += stop - start;
#   pragma omp atomic
    flops += static_cast<double>(2ul * LIBXSTREAM_GETVAL(max_m) * LIBXSTREAM_GETVAL(max_n) * LIBXSTREAM_GETVAL(max_k) * stacksize);
    LIBXSTREAM_PRINT_INFO("libsmm_acc_process (" LIBXSTREAM_DEVICE_NAME "): %.f GFLOP/s", flops / (1E9 * duration));
  }
#endif
}


#if defined(LIBMICSMM_PLANSIZE) && (0 < (LIBMICSMM_PLANSIZE))
template<size_t N, typename T, typename U>
LIBXSTREAM_TARGET(mic) void kernel_plan(const U *LIBXSTREAM_RESTRICT stack, LIBXSTREAM_INVAL(U) max_m, LIBXSTREAM_INVAL(U) max_n, LIBXSTREAM_INVAL(U) max_k,
  const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c)
{
  size_t stacksize = 0;
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_get_shape(0/*current context*/, 0/*stack*/, &stacksize));
  LIBXSTREAM_PRINT_INFO("libsmm_acc_process (" LIBXSTREAM_DEVICE_NAME "): stacksize=%lu max_m=%i max_n=%i max_k=%i", static_cast<unsigned long>(stacksize),
    LIBXSTREAM_GETVAL(max_m), LIBXSTREAM_GETVAL(max_n), LIBXSTREAM_GETVAL(max_k));

#if defined(LIBXSTREAM_PRINT) && defined(_OPENMP)
  const double start = omp_get_wtime();
#endif
#if defined(LIBMICSMM_THREADPRIVATE) && defined(_OPENMP)
# if 1 == (LIBMICSMM_THREADPRIVATE) // native OpenMP TLS
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
# pragma omp threadprivate(tmp)
#else
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static LIBXSTREAM_TLS T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
# endif
#else // without OpenMP nothing needs to be thread-local due to a single-threaded program
  LIBXSTREAM_TARGET(mic) LIBXSTREAM_ALIGNED(static T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
#endif
  const smm_type<T,U> smm(LIBXSTREAM_GETVAL(max_m), LIBXSTREAM_GETVAL(max_n), LIBXSTREAM_GETVAL(max_k));
  const U nstacksize = static_cast<U>(stacksize * N);
  U colspan[LIBMICSMM_PLANSIZE];

  for (U n = 0; n < nstacksize;) {
    int size = 0;

    colspan[0] = n;
    for (; size < (LIBMICSMM_PLANSIZE - 1) && n < nstacksize; n += N) {
      for (U kc0 = stack[n+5], kc1 = (((n + N) < nstacksize) ? stack[n+N+5] : (kc0 + 1)); kc0 == kc1; n += N, kc0 = kc1, kc1 = stack[n+N+5]);
      colspan[++size] = n + N;
    }
    LIBXSTREAM_PRINT_INFO("libsmm_acc_process (" LIBXSTREAM_DEVICE_NAME "): parallel=%lu", static_cast<unsigned long>(size));

#if defined(_OPENMP)
#   pragma omp parallel for schedule(LIBMICSMM_SCHEDULE)
#endif
    for (int i = 0; i < size; ++i) {
      const U j0 = colspan[i], j1 = colspan[i+1];
      const U kc = stack[j0+5] - 1;
      LIBXSTREAM_ASSERT(j1 <= nstacksize);
#if !defined(LIBMICSMM_THREADPRIVATE)
      LIBXSTREAM_ALIGNED(T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBMICSMM_ALIGNMENT);
#endif
      smm.zero_c(tmp);

      for (U j = j0; j < j1; j += N) {
        LIBXSTREAM_ASSERT(j < nstacksize);
        LIBXSTREAM_ASSERT(kc == stack[j+5] - 1);
        LIBXSTREAM_ASSERT(LIBXSTREAM_GETVAL(max_m) == stack[j+0]);
        LIBXSTREAM_ASSERT(LIBXSTREAM_GETVAL(max_n) == stack[j+1]);
        LIBXSTREAM_ASSERT(LIBXSTREAM_GETVAL(max_k) == stack[j+2]);
        const U ka = stack[j+3] - 1, kb = stack[j+4] - 1;
        smm(a + ka, b + kb, tmp);
      }

      smm.copy_c(tmp, c + kc);
    }
  }

#if defined(LIBXSTREAM_PRINT) && defined(_OPENMP)
  static double duration = 0, flops = 0;
  const double stop = omp_get_wtime();
  if (start < stop) {
#   pragma omp atomic
    duration += stop - start;
#   pragma omp atomic
    flops += static_cast<double>(2ul * LIBXSTREAM_GETVAL(max_m) * LIBXSTREAM_GETVAL(max_n) * LIBXSTREAM_GETVAL(max_k) * stacksize);
    LIBXSTREAM_PRINT_INFO("libsmm_acc_process (" LIBXSTREAM_DEVICE_NAME "): %.f GFLOP/s", flops / (1E9 * duration));
  }
#endif
}
#endif // LIBMICSMM_PLANSIZE


template<typename T, bool Complex, typename U>
int process(const U* stack, U stacksize, U nparams, U max_m, U max_n, U max_k, const void* a_data, const void* b_data, void* c_data,
  U def_mnk, void* stream)
{
  LIBXSTREAM_PRINT_INFO("libsmm_acc_process (host): type=%s homogeneous=%s stack=0x%llx a=0x%llx b=0x%llx c=0x%llx stream=0x%llx",
    dbcsr_elem<T,Complex>::name(), 1 == def_mnk ? "true" : "false", reinterpret_cast<unsigned long long>(stack),
    reinterpret_cast<unsigned long long>(a_data), reinterpret_cast<unsigned long long>(b_data), reinterpret_cast<unsigned long long>(c_data),
    reinterpret_cast<unsigned long long>(stream));
  LIBXSTREAM_CHECK_CONDITION(
    stack && 0 <= stacksize && 0 <= nparams
    && 0 <= max_m && 0 <= max_n && 0 <= max_k
    && a_data && b_data && c_data && stream
    && LIBMICSMM_NPARAMS == nparams
    && 1 == def_mnk);

  const size_t shape = stacksize;
  libxstream_argument* signature = 0;
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_signature(&signature));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 0,  stack, libxstream_map_to<U>::type(), 1, &shape));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 1, &max_m, libxstream_map_to<U>::type(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 2, &max_n, libxstream_map_to<U>::type(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 3, &max_k, libxstream_map_to<U>::type(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 4, a_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 5, b_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_inout(signature, 6, c_data, libxstream_map_to<T>::type(), 1, 0/*unknown*/));
#if !defined(LIBMICSMM_PLANSIZE) || (0 >= (LIBMICSMM_PLANSIZE))
  const libxstream_function libmicsmm_process_function = reinterpret_cast<libxstream_function>(kernel<LIBMICSMM_NPARAMS,T,U>);
#else
  const libxstream_function libmicsmm_process_function = reinterpret_cast<libxstream_function>(kernel_plan<LIBMICSMM_NPARAMS,T,U>);
#endif
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_call(libmicsmm_process_function, signature, static_cast<libxstream_stream*>(stream), LIBXSTREAM_CALL_DEFAULT));

  return LIBXSTREAM_ERROR_NONE;
}

} // namespace libmicsmm_process_private

// workaround for issue "cannot find address of function" (use unoptimized build or apply mic attribute globally)
#if !defined(LIBMICSMM_PLANSIZE) || (0 >= (LIBMICSMM_PLANSIZE))
const libxstream_function libmicsmm_process_function = reinterpret_cast<libxstream_function>(libmicsmm_process_private::kernel<LIBMICSMM_NPARAMS,double,int>);
#else
const libxstream_function libmicsmm_process_function = reinterpret_cast<libxstream_function>(libmicsmm_process_private::kernel_plan<LIBMICSMM_NPARAMS,double,int>);
#endif

extern "C" int libsmm_acc_process(void* param_stack, int stacksize, int nparams, int datatype, void* a_data, void* b_data, void* c_data, int max_m, int max_n, int max_k, int def_mnk, void* stream)
{
#if defined(LIBMICSMM_PRETRANSPOSE)
  LIBXSTREAM_ASSERT(false/*TODO: implement C = A * B which is assuming that B is pre-transposed (B^T).*/);
#endif
  const int *const stack = static_cast<const int*>(param_stack);
  int result = LIBXSTREAM_ERROR_NONE;

  switch(static_cast<dbcsr_elem_type>(datatype)) {
    case DBCSR_ELEM_F32: {
#if 0
      result = libmicsmm_process_private::process<float,false>(stack, stacksize, nparams, max_m, max_n, max_k, a_data, b_data, c_data, def_mnk, stream);
#endif
      result = LIBXSTREAM_ERROR_CONDITION;
    } break;
    case DBCSR_ELEM_F64: {
      result = libmicsmm_process_private::process<double,false>(stack, stacksize, nparams, max_m, max_n, max_k, a_data, b_data, c_data, def_mnk, stream);
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
