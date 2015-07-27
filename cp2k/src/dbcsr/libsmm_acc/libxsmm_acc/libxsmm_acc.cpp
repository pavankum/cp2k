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
#include <cstdlib>
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM)
# include <libxstream_end.h>
#endif


#if defined(__LIBXSMM) && !(defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_smm_process_mm_stack_s)(const void*, const int*, const int*, const float*, const float*, float*, void*);
LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_smm_process_mm_stack_d)(const void*, const int*, const int*, const double*, const double*, double*, void*);


LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(xsmm_process_mm_stack_s)(const void* stack_descr, const int* params, const int* stack_size, const float* a_data, const float* b_data, float* c_data, void* error)
{
#if 0 // TODO
  const int result = libsmm_acc_process(void *param_stack, stack_size, LIBXSMM_ACC_NPARAMS, DBCSR_ELEM_F32, a_data, b_data, c_data,
    int m_max, int n_max, int k_max, int def_mnk, 0);
#else
  LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_smm_process_mm_stack_s)(stack_descr, params, stack_size, a_data, b_data, c_data, error);
#endif
}


LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(xsmm_process_mm_stack_d)(const void* stack_descr, const int* params, const int* stack_size, const double* a_data, const double* b_data, double* c_data, void* error)
{
#if 0
    // TODO
#else
  LIBXSMM_ACC_FSYMBOL(dbcsr_mm_hostdrv_mp_smm_process_mm_stack_d)(stack_descr, params, stack_size, a_data, b_data, c_data, error);
#endif
}
#endif


#if defined(__RECONFIGURE)
LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int*, void*);
LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(const int*, void*);
extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams);
extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers);
extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams);
extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers);
extern int LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process);


LIBXSMM_ACC_EXTERN_C void LIBXSMM_ACC_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int* driver, void* error)
{
  // make sure to reconfigure *after* the original configuration procedure ran
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(driver, error);

  static const char *const env = getenv("LIBXSMM_ACC_RECONFIGURE");
  static const bool reconfigure = (env && *env) ? (0 != atoi(env)) : true/*default*/;

  if (reconfigure) {
    const int stacksize = LIBXSMM_ACC_STACKSIZE;
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(&stacksize, error);
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams) = LIBXSMM_ACC_POSTERIOR_STREAMS;
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers) = LIBXSMM_ACC_POSTERIOR_BUFFERS;
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams) = LIBXSMM_ACC_PRIORITY_STREAMS;
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers) = LIBXSMM_ACC_PRIORITY_BUFFERS;
    LIBXSMM_ACC_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process) = LIBXSMM_ACC_MIN_NFLOPS_PERMM;
  }
}
#endif

#endif // defined(__LIBXSMM) || (defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC) && defined(__LIBXSTREAM))
