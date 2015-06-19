/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"

#include <libxstream_begin.h>
#include <cstdlib>
#include <libxstream_end.h>


#if defined(__RECONFIGURE)
LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int*, void*);
LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(const int*, void*);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers);
extern int LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process);


LIBXSTREAM_EXTERN_C void LIBXSTREAM_FSYMBOL(__wrap_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(const int* driver, void* error)
{
  // make sure to follow the original configuration procedure
  LIBXSTREAM_FSYMBOL(__real_dbcsr_config_mp_dbcsr_set_conf_mm_driver)(driver, error);

  static const char *const env = getenv("LIBMICSMM_RECONFIGURE");
  static const bool reconfigure = (env && *env)
    ? 0 != atoi(env)
# if defined(LIBMICSMM_RECONFIGURE)
    : true;
# else
    : false;
# endif

  if (reconfigure) {
    const int stacksize = LIBMICSMM_STACKSIZE;
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_dbcsr_set_conf_mm_stacksize)(&stacksize, error);
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_streams) = LIBMICSMM_POSTERIOR_STREAMS;
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_posterior_buffers) = LIBMICSMM_POSTERIOR_BUFFERS;
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_streams) = LIBMICSMM_PRIORITY_STREAMS;
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_priority_buffers) = LIBMICSMM_PRIORITY_BUFFERS;
    LIBXSTREAM_FSYMBOL(dbcsr_config_mp_accdrv_min_flop_process) = LIBMICSMM_MIN_NFLOPS_PERMM;
  }
}
#endif

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
