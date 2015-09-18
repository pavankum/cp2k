/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#ifndef LIBXSMM_ACC_HPP
#define LIBXSMM_ACC_HPP

#if defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#include "../include/libsmm_acc.h"

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include "../../../acc/mic/libmicacc.h"
#endif

#if defined(__LIBXSMM)
# include <libxsmm.h>
# if (0 != LIBXSMM_ROW_MAJOR)
#   error Please compile LIBXSMM using "make ROW_MAJOR=0 ..."!
# endif
# define LIBXSMM_ACC_ASSERT(A) assert(A)
# define LIBXSMM_ACC_NOT_SUPPORTED 1
# define LIBXSMM_ACC_ERROR_NONE 0
# define LIBXSMM_ACC_ERROR_CONDITION -2
# define LIBXSMM_ACC_CHECK_CONDITION(CONDITION) if (!(CONDITION)) return LIBXSMM_ACC_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT(EXPRESSION) EXPRESSION
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSMM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSMM_ALIGN_VALUE
# define LIBXSMM_ACC_HASH2 LIBXSMM_HASH2
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSMM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGNED LIBXSMM_ALIGNED
# define LIBXSMM_ACC_ALIGNED_STORES LIBXSMM_ALIGNED_STORES
# define LIBXSMM_ACC_ALIGNED_MAX LIBXSMM_ALIGNED_MAX
# define LIBXSMM_ACC_PREFETCH_DECL LIBXSMM_PREFETCH_DECL
# define LIBXSMM_ACC_PREFETCH_ARGA LIBXSMM_PREFETCH_ARGA
# define LIBXSMM_ACC_PREFETCH_ARGB LIBXSMM_PREFETCH_ARGB
# define LIBXSMM_ACC_PREFETCH_ARGC LIBXSMM_PREFETCH_ARGC
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M LIBXSMM_AVG_M
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N LIBXSMM_AVG_N
# define LIBXSMM_ACC_RESTRICT LIBXSMM_RESTRICT
# define LIBXSMM_ACC_EXTERN_C LIBXSMM_EXTERN_C
# define LIBXSMM_ACC_FSYMBOL LIBXSMM_FSYMBOL
# define LIBXSMM_ACC_MOD2 LIBXSMM_MOD2
# define LIBXSMM_ACC_MAX LIBXSMM_MAX
# define LIBXSMM_ACC_RETARGETABLE LIBXSMM_RETARGETABLE
# if defined(LIBXSMM_OFFLOAD_BUILD)
#   define LIBXSMM_ACC_OFFLOAD_BUILD LIBXSMM_OFFLOAD_BUILD
# endif
#else // defined(__LIBXSTREAM)
# define LIBXSMM_ACC_ASSERT LIBXSTREAM_ASSERT
# define LIBXSMM_ACC_NOT_SUPPORTED LIBXSTREAM_NOT_SUPPORTED
# define LIBXSMM_ACC_ERROR_NONE LIBXSTREAM_ERROR_NONE
# define LIBXSMM_ACC_ERROR_CONDITION LIBXSTREAM_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CONDITION LIBXSTREAM_CHECK_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT LIBXSTREAM_CHECK_CALL_ASSERT
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSTREAM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSTREAM_ALIGN_VALUE
# define LIBXSMM_ACC_HASH2 LIBXSTREAM_HASH2
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSTREAM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGNED LIBXSTREAM_ALIGNED
# define LIBXSMM_ACC_ALIGNED_STORES LIBXSTREAM_MAX_SIMD
# define LIBXSMM_ACC_ALIGNED_MAX LIBXSTREAM_MAX_SIMD
# define LIBXSMM_ACC_PREFETCH_DECL(TYPE, ARG)
# define LIBXSMM_ACC_PREFETCH_ARGA(ARG)
# define LIBXSMM_ACC_PREFETCH_ARGB(ARG)
# define LIBXSMM_ACC_PREFETCH_ARGC(ARG)
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_ACC_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M 23
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_ACC_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N 23
# define LIBXSMM_ACC_RESTRICT LIBXSTREAM_RESTRICT
# define LIBXSMM_ACC_EXTERN_C LIBXSTREAM_EXTERN_C
# define LIBXSMM_ACC_FSYMBOL LIBXSTREAM_FSYMBOL
# define LIBXSMM_ACC_MOD2 LIBXSTREAM_MOD2
# define LIBXSMM_ACC_MAX LIBXSTREAM_MAX
# define LIBXSMM_ACC_RETARGETABLE LIBXSTREAM_RETARGETABLE
# if defined(LIBXSTREAM_OFFLOAD_BUILD)
#   define LIBXSMM_ACC_OFFLOAD_BUILD LIBXSTREAM_OFFLOAD_BUILD
# endif
#endif

/** Upper limits for the supported matrix sizes. */
#define LIBXSMM_ACC_MAX_M 368
#define LIBXSMM_ACC_MAX_N 368
#define LIBXSMM_ACC_MAX_K 368

/** Number of parameters per stack entry. */
#define LIBXSMM_ACC_NPARAMS 7

/** Nested parallelism. */
#if defined(_OPENMP) && defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# define LIBXSMM_ACC_OPENMP
#endif

/**
 * Ensures an amortized synchronization overhead.
 * >0: maximum number of locally processed MM
 * =0: process MM using a plan
 * <0: adaptive MM-processing
 */
#if defined(LIBXSMM_ACC_OPENMP)
# if defined(__MIC__)
#   define LIBXSMM_ACC_NLOCAL 128
# else
#   define LIBXSMM_ACC_NLOCAL 16
# endif
#else
# define LIBXSMM_ACC_NLOCAL 128
#endif

/** OpenMP scheduling policy (and chunk size) */
#if defined(__MIC__)
# define LIBXSMM_ACC_SCHEDULE schedule(dynamic)
#else
# define LIBXSMM_ACC_SCHEDULE
#endif

/**
 * Synchronization mechanism.
 * >1: number of locks (POT)
 * =1: omp critical
 * =0: atomic
 */
#if defined(__MIC__)
# define LIBXSMM_ACC_SYNCHRONIZATION 1
#else
# define LIBXSMM_ACC_SYNCHRONIZATION 16
#endif

#if defined(__RECONFIGURE)
/** Reconfigure property if non-zero. */
# define LIBXSMM_ACC_STACKSIZE 1000000
/** Reconfigure property if non-zero. */
# define LIBXSMM_ACC_MULTREC_LIMIT 64
/** Reconfigure property if non-zero. */
# define LIBXSMM_ACC_COMM_THREAD_LOAD 99
# if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_POSTERIOR_STREAMS 1
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_POSTERIOR_BUFFERS 1
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_PRIORITY_STREAMS 1
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_PRIORITY_BUFFERS 1
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_MIN_MFLOPS_PERSTACK 120
/** Reconfigure property if non-zero. */
#   define LIBXSMM_ACC_ACCDRV_MIN_NFLOPS_PERMM 0
# endif
#endif // defined(__RECONFIGURE)

/*#define LIBXSMM_ACC_PRETRANSPOSE*/
/*#define LIBXSMM_ACC_MKLTRANS*/

typedef bool libxsmm_acc_bool_type;

/** Must match stack_descriptor_type. */
struct libxsmm_acc_stackdesc_type {
  int m, n, k, max_m, max_n, max_k;
  libxsmm_acc_bool_type defined_mnk;
};

enum libxsmm_acc_elem_type {
  LIBXSMM_ACC_ELEM_UNKNOWN = 0,
  LIBXSMM_ACC_ELEM_F32 = 1, LIBXSMM_ACC_ELEM_F64 = 3,
  LIBXSMM_ACC_ELEM_C32 = 5, LIBXSMM_ACC_ELEM_C64 = 7
};

template<typename T, libxsmm_acc_bool_type Complex> struct libxsmm_acc_elem {
                                                    static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_UNKNOWN;
                                                    static const char* name() { return "unknown"; } };
template<> struct libxsmm_acc_elem<float,false>   { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_F32;
                                                    static const char* name() { return "f32"; } };
template<> struct libxsmm_acc_elem<double,false>  { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_F64;
                                                    static const char* name() { return "f64"; } };
template<> struct libxsmm_acc_elem<float,true>    { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_C32;
                                                    static const char* name() { return "c32"; } };
template<> struct libxsmm_acc_elem<double,true>   { static const libxsmm_acc_elem_type type = LIBXSMM_ACC_ELEM_C64;
                                                    static const char* name() { return "c64"; } };

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#endif // LIBXSMM_ACC_HPP
