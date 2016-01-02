/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2016  CP2K developers group                         *
 *****************************************************************************/

//! **************************************************************************
//!> \author Hans Pabst (Intel Corp.)
//! **************************************************************************

#if defined(__ACC) && defined(__ACC_MIC)

#include "libmicacc.h"
#include <libxstream.h>
#include <stdlib.h>


int acc_get_ndevices(int* n_devices)
{
  LIBXSTREAM_CHECK_CONDITION(0 != n_devices);
  size_t ndevices;
  const int result = libxstream_get_ndevices(&ndevices);
  *n_devices = ndevices;
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_set_active_device(int device_id)
{
  const int result = libxstream_set_active_device(device_id);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_stream_priority_range(int* least, int* greatest)
{
  const int result = libxstream_stream_priority_range(least, greatest);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_stream_create(void** stream_p, const char* name, int priority)
{
  int device = -1, result = libxstream_get_active_device(&device);
  LIBXSTREAM_CHECK_ERROR(result);
  result = libxstream_stream_create((libxstream_stream**)stream_p, device, priority, name);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_stream_destroy(void* stream)
{
  const int result = libxstream_stream_destroy((libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_stream_sync(void* stream)
{
  const int result = libxstream_stream_wait((libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_stream_wait_event(void* stream, void* event)
{
  const int result = libxstream_stream_wait_event((libxstream_stream*)stream, (libxstream_event*)event);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_event_create(void** event_p)
{
  const int result = libxstream_event_create((libxstream_event**)event_p);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_event_destroy(void* event)
{
  const int result = libxstream_event_destroy((libxstream_event*)event);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_event_record(void* event, void* stream)
{
  const int result = libxstream_event_record((libxstream_event*)event, (libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_event_query(void* event, int* has_occured)
{
  const int result = libxstream_event_query((libxstream_event*)event, has_occured);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_event_synchronize(void* event)
{
  const int result = libxstream_event_wait((libxstream_event*)event);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_dev_mem_allocate(void** dev_mem, size_t n)
{
  int device = -1, result = libxstream_get_active_device(&device);
  LIBXSTREAM_CHECK_ERROR(result);
  result = libxstream_mem_allocate(device, dev_mem, n, 0/*automatic*/);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_dev_mem_deallocate(void* dev_mem)
{
  int device = -1, result = libxstream_get_active_device(&device);
  LIBXSTREAM_CHECK_ERROR(result);
  result = libxstream_mem_deallocate(device, dev_mem);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_host_mem_allocate(void** host_mem, size_t n, void* stream)
{
  const int result = libxstream_mem_allocate(-1, host_mem, n, 0);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_host_mem_deallocate(void* host_mem, void* stream)
{
  const int result = libxstream_mem_deallocate(-1, host_mem);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_memcpy_h2d(const void* host_mem, void* dev_mem, size_t count, void* stream)
{
  const int result = libxstream_memcpy_h2d(host_mem, dev_mem, count, (libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_memcpy_d2h(const void* dev_mem, void* host_mem, size_t count, void* stream)
{
  const int result = libxstream_memcpy_d2h(dev_mem, host_mem, count, (libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_memcpy_d2d(const void* devmem_src, void* devmem_dst, size_t count, void* stream)
{
  const int result = libxstream_memcpy_d2d(devmem_src, devmem_dst, count, (libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_memset_zero(void* dev_mem, size_t offset, size_t length, void* stream)
{
  const int result = libxstream_memset_zero((char*)dev_mem + offset, length, (libxstream_stream*)stream);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}


int acc_dev_mem_info(size_t* free, size_t* avail)
{
  int device = -1, result = libxstream_get_active_device(&device);
  LIBXSTREAM_CHECK_ERROR(result);
  result = libxstream_get_meminfo(device, free, avail);
  LIBXSTREAM_ASSERT(LIBXSTREAM_ERROR_NONE == result);
  return result;
}

#endif // defined(__ACC) && defined(__ACC_MIC)
