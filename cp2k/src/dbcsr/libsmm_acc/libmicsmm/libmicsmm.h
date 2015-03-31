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

/** Maximum number of matrices potentially processed in parallel. */
#define LIBMICSMM_MAX_BURST 65536

/** Determines if LIBXSMM is used. */
#define LIBMICSMM_USE_LIBXSMM

/** OpenMP scheduling policy (and chunk size) */
#define LIBMICSMM_SCHEDULE dynamic

/** The kind of thread-private data. */
#define LIBMICSMM_THREADPRIVATE 1

/*#define LIBMICSMM_USE_PRETRANSPOSE*/
/*#define LIBMICSMM_USE_MKLTRANS*/


typedef enum dbcsr_elem_type {
  DBCSR_ELEM_UNKNOWN = 0,
  DBCSR_ELEM_F32 = 1, DBCSR_ELEM_F64 = 3,
  DBCSR_ELEM_C32 = 5, DBCSR_ELEM_C64 = 7
} dbcsr_elem_type;

#endif /*defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)*/
#endif /*LIBMICSMM_H*/
