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
#define LIBMICSMM_MAX_M 256
#define LIBMICSMM_MAX_N 256
#define LIBMICSMM_MAX_K 256

/** Number of parameters per stack entry. */
#define LIBMICSMM_NPARAMS 7

/** Maximum number of matrices potentially processed in parallel. */
#define LIBMICSMM_MAX_BURST 16384

#define LIBMICSMM_USE_LOOPHINTS
#define LIBMICSMM_USE_LIBXSMM
//#define LIBMICSMM_USE_MKLTRANS
//#define LIBMICSMM_USE_MKLSMM
//#define LIBMICSMM_USE_DUMP


typedef enum dbcsr_elem_type {
  DBCSR_ELEM_UNKNOWN = 0,
  DBCSR_ELEM_F32 = 1, DBCSR_ELEM_F64 = 3,
  DBCSR_ELEM_C32 = 5, DBCSR_ELEM_C64 = 7
} dbcsr_elem_type;


LIBXSTREAM_EXTERN_C int libsmm_acc_file_save(const char groupname[], const char name[], size_t id, const void* data, size_t data_size, const void* header, size_t header_size);
LIBXSTREAM_EXTERN_C int libsmm_acc_file_load(const char groupname[], const char name[], size_t id, void* data, size_t* data_size, void* header);
LIBXSTREAM_EXTERN_C int libsmm_acc_file_diff(const char groupname[], const char name[], size_t id, const void* data, dbcsr_elem_type elem_type, double* max_diff);

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
#endif // LIBMICSMM_H
