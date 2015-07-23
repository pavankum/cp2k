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

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include "../include/libsmm_acc.h"
# include "../../../acc/mic/libmicacc.h"
#endif

#if defined(__LIBXSMM)
# include <libxsmm.h>
# if (0 != LIBXSMM_ROW_MAJOR)
#   error Please compile LIBXSMM using "make ROW_MAJOR=0 ..."!
# endif
# define LIBXSMM_ACC_ASSERT assert
# define LIBXSMM_ACC_NOT_SUPPORTED 1
# define LIBXSMM_ACC_ERROR_NONE 0
# define LIBXSMM_ACC_ERROR_CONDITION -2
# define LIBXSMM_ACC_CHECK_CONDITION(CONDITION) if (!(CONDITION)) return LIBXSMM_ACC_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT(EXPRESSION) LIBXSMM_ACC_ASSERT(LIBXSMM_ACC_ERROR_NONE == (EXPRESSION))
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSMM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSMM_ALIGN_VALUE
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSMM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGNED LIBXSMM_ALIGNED
# define LIBXSMM_ACC_ALIGNED_STORES LIBXSMM_ALIGNED_STORES
# define LIBXSMM_ACC_ALIGNED_MAX LIBXSMM_ALIGNED_MAX
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M LIBXSMM_AVG_M
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N LIBXSMM_AVG_N
# define LIBXSMM_ACC_RESTRICT LIBXSMM_RESTRICT
# define LIBXSMM_ACC_EXTERN_C LIBXSMM_EXTERN_C
# define LIBXSMM_ACC_FSYMBOL LIBXSMM_FSYMBOL
# define LIBXSMM_ACC_TARGET LIBXSMM_TARGET
# define LIBXSMM_ACC_MOD LIBXSMM_MOD
# define LIBXSMM_ACC_MAX LIBXSMM_MAX
#else // defined(__LIBXSTREAM)
# define LIBXSMM_ACC_ASSERT LIBXSTREAM_ASSERT
# define LIBXSMM_ACC_NOT_SUPPORTED LIBXSTREAM_NOT_SUPPORTED
# define LIBXSMM_ACC_ERROR_NONE LIBXSTREAM_ERROR_NONE
# define LIBXSMM_ACC_ERROR_CONDITION LIBXSTREAM_ERROR_CONDITION
# define LIBXSMM_ACC_CHECK_CONDITION LIBXSTREAM_CHECK_CONDITION
# define LIBXSMM_ACC_CHECK_CALL_ASSERT LIBXSTREAM_CHECK_CALL_ASSERT
# define LIBXSMM_ACC_PRAGMA_LOOP_COUNT LIBXSTREAM_PRAGMA_LOOP_COUNT
# define LIBXSMM_ACC_ALIGN_VALUE LIBXSTREAM_ALIGN_VALUE
# define LIBXSMM_ACC_ASSUME_ALIGNED LIBXSTREAM_ASSUME_ALIGNED
# define LIBXSMM_ACC_ALIGNED LIBXSTREAM_ALIGNED
# define LIBXSMM_ACC_ALIGNED_STORES LIBXSTREAM_MAX_SIMD
# define LIBXSMM_ACC_ALIGNED_MAX LIBXSTREAM_MAX_SIMD
# define LIBXSMM_ACC_LOOP_MAX_M LIBXSMM_ACC_MAX_M
# define LIBXSMM_ACC_LOOP_AVG_M 23
# define LIBXSMM_ACC_LOOP_MAX_N LIBXSMM_ACC_MAX_N
# define LIBXSMM_ACC_LOOP_AVG_N 23
# define LIBXSMM_ACC_RESTRICT LIBXSTREAM_RESTRICT
# define LIBXSMM_ACC_EXTERN_C LIBXSTREAM_EXTERN_C
# define LIBXSMM_ACC_FSYMBOL LIBXSTREAM_FSYMBOL
# define LIBXSMM_ACC_TARGET LIBXSTREAM_TARGET
# define LIBXSMM_ACC_MOD LIBXSTREAM_MOD
# define LIBXSMM_ACC_MAX LIBXSTREAM_MAX
#endif

/** Upper limits for the supported matrix sizes. */
#define LIBXSMM_ACC_MAX_M 368
#define LIBXSMM_ACC_MAX_N 368
#define LIBXSMM_ACC_MAX_K 368

/** Number of parameters per stack entry. */
#define LIBXSMM_ACC_NPARAMS 7

/**
 * Ensures an amortized synchronization overhead.
 * >0: maximum number of locally processed MM
 * =0: process MM using a plan
 * <0: adaptive MM-processing
 */
#define LIBXSMM_ACC_NLOCAL 128

/** OpenMP scheduling policy (and chunk size) */
#define LIBXSMM_ACC_SCHEDULE dynamic

/**
 * Synchronization mechanism.
 * >1: number of locks (POT)
 * =1: omp critical
 * =0: atomic
 */
#define LIBXSMM_ACC_SYNCHRONIZATION 1

/** Determines if CP2K/ACC is reconfigured. */
#define LIBXSMM_ACC_STACKSIZE 1000000
#define LIBXSMM_ACC_POSTERIOR_STREAMS 1
#define LIBXSMM_ACC_POSTERIOR_BUFFERS 1
#define LIBXSMM_ACC_PRIORITY_STREAMS 1
#define LIBXSMM_ACC_PRIORITY_BUFFERS 1
#define LIBXSMM_ACC_MIN_MFLOPS_PERSTACK 120
#define LIBXSMM_ACC_MIN_NFLOPS_PERMM 0

/*#define LIBXSMM_ACC_PRETRANSPOSE*/
/*#define LIBXSMM_ACC_MKLTRANS*/


typedef enum dbcsr_elem_type {
  DBCSR_ELEM_UNKNOWN = 0,
  DBCSR_ELEM_F32 = 1, DBCSR_ELEM_F64 = 3,
  DBCSR_ELEM_C32 = 5, DBCSR_ELEM_C64 = 7
} dbcsr_elem_type;

template<typename T, bool Complex> struct dbcsr_elem  { static const dbcsr_elem_type type = DBCSR_ELEM_UNKNOWN;
                                                        static const char* name() { return "unknown"; } };
template<> struct dbcsr_elem<float,false>             { static const dbcsr_elem_type type = DBCSR_ELEM_F32;
                                                        static const char* name() { return "f32"; } };
template<> struct dbcsr_elem<double,false>            { static const dbcsr_elem_type type = DBCSR_ELEM_F64;
                                                        static const char* name() { return "f64"; } };
template<> struct dbcsr_elem<float,true>              { static const dbcsr_elem_type type = DBCSR_ELEM_C32;
                                                        static const char* name() { return "c32"; } };
template<> struct dbcsr_elem<double,true>             { static const dbcsr_elem_type type = DBCSR_ELEM_C64;
                                                        static const char* name() { return "c64"; } };

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
#endif // LIBXSMM_ACC_HPP
