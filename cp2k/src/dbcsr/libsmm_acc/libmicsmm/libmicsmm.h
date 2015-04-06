/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#ifndef LIBMICSMM_H
#define LIBMICSMM_H

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "../include/libsmm_acc.h"
#include "../../../acc/mic/libmicacc.h"

/** Upper limits for the supported matrix sizes. */
#define LIBMICSMM_MAX_M 368
#define LIBMICSMM_MAX_N 368
#define LIBMICSMM_MAX_K 368

/** Number of parameters per stack entry. */
#define LIBMICSMM_NPARAMS 7

/**
 * Ensures an amortized synchronization overhead.
 * >0: maximum number of locally processed MM
 * =0: process MM using a plan
 * <0: adaptive MM-processing
 */
#define LIBMICSMM_NLOCAL 128

/** OpenMP scheduling policy (and chunk size) */
#define LIBMICSMM_SCHEDULE dynamic

/**
 * Synchronization mechanism.
 * >1: number of locks
 * =1: omp critical
 * =0: atomic
 */
#define LIBMICSMM_SYNCHRONIZATION 1

/** The kind of thread-private data. */
#define LIBMICSMM_THREADPRIVATE 1

/** Determines if CP2K/ACC is reconfigured. */
#define LIBMICSMM_STACKSIZE 1000000
#define LIBMICSMM_POSTERIOR_STREAMS 1
#define LIBMICSMM_POSTERIOR_BUFFERS 1
#define LIBMICSMM_PRIORITY_STREAMS 1
#define LIBMICSMM_PRIORITY_BUFFERS 1

/** Determines if LIBXSMM is used. */
#define LIBMICSMM_LIBXSMM

/*#define LIBMICSMM_PRETRANSPOSE*/
/*#define LIBMICSMM_MKLTRANS*/


typedef enum dbcsr_elem_type {
  DBCSR_ELEM_UNKNOWN = 0,
  DBCSR_ELEM_F32 = 1, DBCSR_ELEM_F64 = 3,
  DBCSR_ELEM_C32 = 5, DBCSR_ELEM_C64 = 7
} dbcsr_elem_type;

#endif /*defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)*/
#endif /*LIBMICSMM_H*/
