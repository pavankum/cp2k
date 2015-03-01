/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"
#include <cstdlib>
#include <vector>

#include <libxstream_begin.h>
#include <cstdio>
#if defined(LIBMICSMM_USE_LIBXSMM) && (2 != (LIBMICSMM_USE_LIBXSMM+1)) && defined(__LIBXSMM) && defined(__MIC__)
# include <libxsmm.h>
#endif
#if defined(LIBMICSMM_USE_MKLSMM) && defined(__MKL)
# if !defined(MKL_DIRECT_CALL_SEQ) && !defined(MKL_DIRECT_CALL)
#   define MKL_DIRECT_CALL_SEQ
# endif
# include <mkl.h>
#endif
#if defined(_OPENMP)
# include <omp.h>
#endif
#include <libxstream_end.h>

#define LIBMICSMM_MAX_RESULT_SIZE (LIBMICSMM_MAX_M * LIBMICSMM_MAX_N)


namespace libmicsmm_process_private {

template<typename T, typename U>
class LIBXSTREAM_TARGET(mic) smm_type {
  U M, N, K, LDC;
#if defined(LIBMICSMM_USE_LIBXSMM) && (2 != (LIBMICSMM_USE_LIBXSMM+1)) && defined(__LIBXSMM) && defined(__MIC__) && (0 != LIBXSMM_COL_MAJOR)
  libxsmm_mm_dispatch<T> xmm;
#endif

public:
  smm_type(U M_, U N_, U K_, U LDC_ = 0)
    : M(M_), N(N_), K(K_), LDC(0 == LDC_ ? M_ : LDC_)
#if defined(LIBMICSMM_USE_LIBXSMM) && (2 != (LIBMICSMM_USE_LIBXSMM+1)) && defined(__LIBXSMM) && defined(__MIC__) && (0 != LIBXSMM_COL_MAJOR)
    , xmm(LDC == M ? libxsmm_mm_dispatch<T>(M_, N_, K_) : libxsmm_mm_dispatch<T>())
#endif
  {}

  void zero_c(T *LIBXSTREAM_RESTRICT c) const {
    LIBXSTREAM_ASSUME_ALIGNED(c, LIBXSTREAM_MAX_SIMD);
#if defined(__INTEL_COMPILER)
# if defined(LIBMICSMM_USE_LOOPHINTS)
#   pragma loop_count min(1), max(LIBMICSMM_MAX_RESULT_SIZE), avg(23*23)
# endif
#   pragma simd aligned(c:1)
#elif defined(_OPENMP)
#   pragma omp simd aligned(c:1)
#endif
    for (U i = 0; i < (N * LDC); ++i) {
      c[i] = 0;
    }
  }

  void copy_c(const T *LIBXSTREAM_RESTRICT c, T *LIBXSTREAM_RESTRICT out) const {
    LIBXSTREAM_ASSUME_ALIGNED(c, LIBXSTREAM_MAX_SIMD);
#if defined(__INTEL_COMPILER)
# if defined(LIBMICSMM_USE_LOOPHINTS)
#   pragma loop_count min(1), max(LIBMICSMM_MAX_N), avg(23)
# endif
#   pragma vector nontemporal(out)
#   pragma simd collapse(2)
#elif defined(_OPENMP)
#   pragma omp simd collapse(2)
#endif
    for (U j = 0; j < N; ++j) {
#if defined(__INTEL_COMPILER)
#     pragma unroll(16)
#endif
      for (U i = 0; i < M; ++i) {
#if defined(_OPENMP)
#       pragma omp atomic
#endif
        out[j*M+i] += c[j*LDC+i];
      }
    }
  }

#if defined(MKL_DIRECT_CALL_SEQ)
  void blasmm(const float* a, const float* b, float* c) const {
    static float alpha = 1.f, beta = 1.f;
    static char trans = 'N';
    int m = M, n = N, k = K, ldc = LDC;
    sgemm(&trans, &trans, &m, &n, &k, &alpha, const_cast<float*>(a), &m, const_cast<float*>(b), &k, &beta, c, &ldc);
  }

  void blasmm(const double* a, const double* b, double* c) const {
    static double alpha = 1.0, beta = 1.0;
    static char trans = 'N';
    int m = M, n = N, k = K, ldc = LDC;
    dgemm(&trans, &trans, &m, &n, &k, &alpha, const_cast<double*>(a), &m, const_cast<double*>(b), &k, &beta, c, &ldc);
  }
#endif

  void operator()(const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c) const {
#if defined(LIBMICSMM_USE_LIBXSMM) && (2 != (LIBMICSMM_USE_LIBXSMM+1)) && defined(__LIBXSMM) && defined(__MIC__) && (0 != LIBXSMM_COL_MAJOR)
    if (0 != xmm) {
      (*xmm)(a, b, c);
    }
    else if (LIBXSMM_MAX_MNK >= (M * N * K)) {
      libxsmm_imm(M, N, K, a, b, c);
    }
    else {
      libxsmm_blasmm(M, N, K, a, b, c);
    }
#elif defined(MKL_DIRECT_CALL_SEQ)
    blasmm(a, b, c);
#else
# if defined(__INTEL_COMPILER)
#   if defined(LIBMICSMM_USE_LOOPHINTS)
#   pragma loop_count min(1), max(LIBMICSMM_MAX_M), avg(23)
#   endif
#   pragma vector nontemporal(c)
#   pragma simd collapse(2)
# elif defined(_OPENMP)
#   pragma omp simd collapse(2)
# endif
    for (U j = 0; j < M; ++j) {
      for (U i = 0; i < N; ++i) {
        const U index = i * LDC + j;
        T r = c[index];
# if defined(__INTEL_COMPILER)
#       pragma unroll(16)
#       pragma simd reduction(+:r)
# elif defined(_OPENMP)
#       pragma omp simd reduction(+:r)
# endif
        for (U k = 0; k < K; ++k) {
          const T aj = a[k*M+j];
          const T bk = b[i*K+k];
          r += aj * bk;
        }
        c[index] = r;
      }
    }
#endif
  }
};


template<size_t N, typename T, typename U>
LIBXSTREAM_TARGET(mic) void kernel(const U *LIBXSTREAM_RESTRICT stack, U max_m, U max_n, U max_k,
  const T *LIBXSTREAM_RESTRICT a, const T *LIBXSTREAM_RESTRICT b, T *LIBXSTREAM_RESTRICT c)
{
#if defined(LIBMICSMM_USE_PRETRANSPOSE)
  LIBXSTREAM_ASSERT(false/*TODO: implement C = A * B which is assuming that B is pre-transposed (B^T).*/);
#endif
#if defined(LIBXSTREAM_DEBUG) && defined(_OPENMP)
  const double start = omp_get_wtime();
#endif

  const smm_type<T,U> smm(max_m, max_n, max_k/*, LIBMICSMM_MAX_M*/);
  U colspan[LIBMICSMM_MAX_BURST];

  size_t stacksize = 0;
  libxstream_get_shape(0, 0/*stack*/, &stacksize);
  const U n = static_cast<U>(stacksize * N);

  for (U s = N; s <= n;) {
    int size = 0;

    colspan[0] = s - N;
    for (; size < (LIBMICSMM_MAX_BURST - 1) && s <= n; s += N) {
      for (U kc1 = stack[s+5-N], kc2 = stack[s+5]; kc1 == kc2; s += N, kc1 = kc2, kc2 = stack[s+5]);
      colspan[++size] = s < n ? s : n;
    }

#if defined(_OPENMP)
#   pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < size; ++i) {
      LIBXSTREAM_ASSERT(colspan[i] < n && max_m == stack[colspan[i]+0] && max_n == stack[colspan[i]+1] && max_k == stack[colspan[i]+2]);
      const U j0 = colspan[i], j1 = colspan[i+1], kc = stack[j0+5] - 1;

      LIBXSTREAM_ALIGNED(T tmp[LIBMICSMM_MAX_RESULT_SIZE], LIBXSTREAM_MAX_SIMD);
      LIBXSTREAM_ASSERT(0 == (reinterpret_cast<uintptr_t>(tmp) % LIBXSTREAM_MAX_SIMD));
      smm.zero_c(tmp);

      for (U j = j0; j < j1; j += N) {
        LIBXSTREAM_ASSERT(j < n);
        const U ka = stack[j+3] - 1, kb = stack[j+4] - 1;
        smm(a + ka, b + kb, tmp);
      }

      smm.copy_c(tmp, c + kc);
    }
  }

#if defined(LIBXSTREAM_DEBUG) && defined(_OPENMP)
  static double duration = 0, flops = 0;
  const double stop = omp_get_wtime();
  if (start < stop) {
#   pragma omp atomic
    duration += stop - start;
#   pragma omp atomic
    flops += static_cast<double>(2ul * max_m * max_n * max_k * stacksize);
    LIBXSTREAM_PRINT_INFO("libsmm_acc_process: %.f GFLOP/s", flops / (1E9 * duration));
  }
#endif
}


template<typename T, bool Complex, typename U>
int process(const U* stack, U stacksize, U nparams, U max_m, U max_n, U max_k, const void* a_data, const void* b_data, void* c_data,
  U def_mnk, void* stream)
{
  LIBXSTREAM_PRINT_INFOCTX("type=%s size=%i homogeneous=%s max_m=%i max_n=%i max_k=%i a=0x%lx b=0x%lx c=0x%lx stream=0x%lx",
    dbcsr_elem<T,Complex>::name(), stacksize, 1 == def_mnk ? "true" : "false", max_m, max_n, max_k,
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(a_data)),
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(b_data)),
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(c_data)),
    static_cast<unsigned long>(reinterpret_cast<uintptr_t>(stream)));
  LIBXSTREAM_CHECK_CONDITION(
    stack && 0 <= stacksize && 0 <= nparams
    && 0 <= max_m && 0 <= max_n && 0 <= max_k
    && a_data && b_data && c_data && stream
    && LIBMICSMM_NPARAMS == nparams
    && 1 == def_mnk);
#if defined(LIBMICSMM_USE_DUMP)
  static size_t id = 0;
  static const char *const groupname = getenv("LIBMICSMM_DUMP");
  std::vector<U> stack_buffer(stacksize * LIBMICSMM_NPARAMS);
  LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(stack, &stack_buffer[0], stack_buffer.size() * sizeof(U), stream));
  LIBXSTREAM_CHECK_CALL(acc_stream_sync(stream));
  LIBXSTREAM_CHECK_CALL(libsmm_acc_file_save(groupname, "stack", id, &stack_buffer[0], stack_buffer.size() * sizeof(U), &def_mnk, sizeof(def_mnk)));
  U size_a = 0, size_b = 0, size_c = 0;
  for (U i = 0; i < stack_buffer.size(); i += LIBMICSMM_NPARAMS) {
    const U ka1 = stack_buffer[i+3], kb1 = stack_buffer[i+4], kc1 = stack_buffer[i+5];
    size_a = ka1 <= size_a ? size_a : ka1; // size_a = max(ka1, size_a)
    size_b = kb1 <= size_b ? size_b : kb1; // size_b = max(kb1, size_b)
    size_c = kc1 <= size_c ? size_c : kc1; // size_c = max(kc1, size_c)
  }
  size_a += max_m * max_k - 1;
  size_b += max_k * max_n - 1;
  size_c += max_m * max_n - 1;
  std::vector<T> buffer(size_a);
  LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(a_data, &buffer[0], size_a * sizeof(T), stream));
  LIBXSTREAM_CHECK_CALL(acc_stream_sync(stream));
  LIBXSTREAM_CHECK_CALL(libsmm_acc_file_save(groupname, "adata", id, &buffer[0], size_a * sizeof(T), &max_m, sizeof(max_m)));
  buffer.resize(size_b);
  LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(b_data, &buffer[0], size_b * sizeof(T), stream));
  LIBXSTREAM_CHECK_CALL(acc_stream_sync(stream));
  LIBXSTREAM_CHECK_CALL(libsmm_acc_file_save(groupname, "bdata", id, &buffer[0], size_b * sizeof(T), &max_k, sizeof(max_k)));
  buffer.resize(size_c);
  LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(c_data, &buffer[0], size_c * sizeof(T), stream));
  LIBXSTREAM_CHECK_CALL(acc_stream_sync(stream));
  LIBXSTREAM_CHECK_CALL(libsmm_acc_file_save(groupname, "cdata", id, &buffer[0], size_c * sizeof(T), &max_n, sizeof(max_n)));
#endif

  const libxstream_function function = reinterpret_cast<libxstream_function>(kernel<LIBMICSMM_NPARAMS,T,U>);
  const size_t shape = stacksize;
  libxstream_argument* signature = 0;
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_signature(&signature));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 0,  stack, libxstream_type2value<U>::value(), 1, &shape));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 1, &max_m, libxstream_type2value<U>::value(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 2, &max_n, libxstream_type2value<U>::value(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 3, &max_k, libxstream_type2value<U>::value(), 0, 0));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 4, a_data, libxstream_type2value<T>::value(), 1, 0/*unknown*/));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_input(signature, 5, b_data, libxstream_type2value<T>::value(), 1, 0/*unknown*/));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_inout(signature, 6, c_data, libxstream_type2value<T>::value(), 1, 0/*unknown*/));
  LIBXSTREAM_CHECK_CALL_ASSERT(libxstream_fn_call(function, signature, static_cast<libxstream_stream*>(stream), LIBXSTREAM_CALL_DEFAULT));

#if defined(LIBMICSMM_USE_DUMP)
  LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(c_data, &buffer[0], size_c * sizeof(T), stream));
  LIBXSTREAM_CHECK_CALL(acc_stream_sync(stream));
  LIBXSTREAM_CHECK_CALL(libsmm_acc_file_save(groupname, "cgold", id, &buffer[0], size_c * sizeof(T), 0, 0));
  ++id;
#endif // LIBMICSMM_USE_DUMP

  return LIBXSTREAM_ERROR_NONE;
}

} // namespace libmicsmm_process_private


extern "C" int libsmm_acc_process(void* param_stack, int stacksize, int nparams, int datatype, void* a_data, void* b_data, void* c_data, int max_m, int max_n, int max_k, int def_mnk, void* stream)
{
  const int *const stack = static_cast<const int*>(param_stack);
  int result = LIBXSTREAM_ERROR_NONE;

  switch(static_cast<dbcsr_elem_type>(datatype)) {
    case DBCSR_ELEM_F32: {
      //result = libmicsmm_process_private::process<float,false>(stack, stacksize, nparams, max_m, max_n, max_k, a_data, b_data, c_data, def_mnk, stream);
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
