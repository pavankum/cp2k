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
  // make sure to follow the original configuration procedure
  LIBXSMM_ACC_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(driver, error);

  static const char *const env = getenv("LIBXSMM_ACC_RECONFIGURE");
  static const bool reconfigure = (env && *env && 0 != atoi(env))

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
