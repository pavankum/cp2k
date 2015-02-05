/*****************************************************************************
 *  CP2K: A general program to perform molecular dynamics simulations        *
 *  Copyright (C) 2000 - 2015 the CP2K developers group                      *
 *****************************************************************************/
/* Hans Pabst (Intel Corp.)
******************************************************************************/
#if defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)

#include "libmicsmm.hpp"

#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <fstream>
#include <vector>
#include <limits>
#include <cstdio>

#if defined(_OPENMP)
# include <omp.h>
#endif


namespace libmicsmm_private {

void libsmm_acc_file_diff(double& a, double b, int position)
{
#if defined(LIBXSTREAM_DEBUG)
  if (a < b) {
    LIBXSTREAM_PRINT_INFOCTX("%f @ %i", b, position);
    a = b;
  }
#else
  a = std::max(a, b);
#endif
}

} // namespace libmicsmm_private


LIBXSTREAM_EXTERN_C int libsmm_acc_file_save(const char groupname[], const char name[], size_t id, const void* data, size_t data_size, const void* header, size_t header_size)
{
  LIBXSTREAM_PRINT_INFOCTX("name=\"%s\" id=%lu", (name && *name) ? name : "<invalid>", static_cast<unsigned long>(id));
  LIBXSTREAM_CHECK_CONDITION(name && *name && (0 == header_size || header));
  char buffer[1024];
  LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%s%s%s-%d.bin",
    (groupname && *groupname) ? groupname : "",
    (groupname && *groupname) ? "-" : "", name, id);
  std::ofstream file(buffer, std::ios_base::binary);

  if (file) file << data_size << ' ' << header_size;
  if (file) file.write(static_cast<const char*>(header), header_size);
  if (file) file.write(static_cast<const char*>(data), data_size);
  return file ? LIBXSTREAM_ERROR_NONE : LIBXSTREAM_ERROR_RUNTIME;
}


LIBXSTREAM_EXTERN_C int libsmm_acc_file_load(const char groupname[], const char name[], size_t id, void* data, size_t* data_size, void* header)
{
  LIBXSTREAM_PRINT_INFOCTX("name=\"%s\" id=%lu", (name && *name) ? name : "<invalid>", static_cast<unsigned long>(id));
  LIBXSTREAM_CHECK_CONDITION(name && *name);
  char buffer[1024];
  LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%s%s%s-%d.bin",
    (groupname && *groupname) ? groupname : "",
    (groupname && *groupname) ? "-" : "", name, id);
  std::ifstream file(buffer, std::ios_base::binary);
  size_t size = 0, header_size = 0;

  if (!data_size) data_size = &size;
  if (file) file >> *data_size >> header_size;
  if (file) {
    char *const header_data = header ? static_cast<char*>(header) : (sizeof(buffer) >= header_size ? buffer : 0);
    LIBXSTREAM_CHECK_CONDITION(0 != header_data);
    file.read(header_data, header_size);
  }
  if (file && data) file.read(static_cast<char*>(data), *data_size);
  return file ? LIBXSTREAM_ERROR_NONE : LIBXSTREAM_ERROR_RUNTIME;
}


LIBXSTREAM_EXTERN_C int libsmm_acc_file_diff(const char groupname[], const char name[], size_t id, const void* data, dbcsr_elem_type elem_type, double* max_diff)
{
  LIBXSTREAM_PRINT_INFOCTX("name=\"%s\" id=%lu", (name && *name) ? name : "<invalid>", static_cast<unsigned long>(id));
  LIBXSTREAM_CHECK_CONDITION(name && *name && data);
  size_t data_size = 0, read_size = 0;

  if (max_diff) {
    char buffer[1024];
    LIBXSTREAM_SNPRINTF(buffer, sizeof(buffer), "%s%s%s-%d.bin",
      (groupname && *groupname) ? groupname : "",
      (groupname && *groupname) ? "-" : "", name, id);
    std::ifstream file(buffer, std::ios_base::binary);
    size_t header_size = 0;

    if (file) file >> data_size >> header_size;
    if (file) {
      LIBXSTREAM_CHECK_CONDITION(sizeof(buffer) >= header_size);
      file.read(buffer, header_size);
    }

    switch (elem_type) {
      case DBCSR_ELEM_F32: {
        const float* value = static_cast<const float*>(data);
        float expect = 0;
        while (read_size < data_size && file) {
          file.read(reinterpret_cast<char*>(&expect), sizeof(float));
          libmicsmm_private::libsmm_acc_file_diff(*max_diff, std::abs(static_cast<double>(*value) - static_cast<double>(expect)), value - static_cast<const float*>(data));
          read_size += sizeof(float);
          ++value;
        }
      } break;
      case DBCSR_ELEM_F64: {
        const double* value = static_cast<const double*>(data);
        double expect = 0;
        while (read_size < data_size && file) {
          file.read(reinterpret_cast<char*>(&expect), sizeof(double));
          libmicsmm_private::libsmm_acc_file_diff(*max_diff, std::abs(*value - expect), value - static_cast<const double*>(data));
          read_size += sizeof(double);
          ++value;
        }
      } break;
      case DBCSR_ELEM_C32: {
        LIBXSTREAM_CHECK_CONDITION(false/*TODO*/);
      } break;
      case DBCSR_ELEM_C64: {
        LIBXSTREAM_CHECK_CONDITION(false/*TODO*/);
      } break;
      default: { // binary comparison
        const char* value = static_cast<const char*>(data);
        char expect = 0;
        while (read_size < data_size && file) {
          file.read(&expect, sizeof(char));
          libmicsmm_private::libsmm_acc_file_diff(*max_diff, std::abs(*value - expect), value - static_cast<const char*>(data));
          read_size += sizeof(char);
          ++value;
        }
      } break;
    }
  }

  return read_size == data_size ? LIBXSTREAM_ERROR_NONE : LIBXSTREAM_ERROR_RUNTIME;
}


#if defined(LIBMICSMM_USE_STANDALONE) && !defined(LIBMICSMM_USE_DUMP)
struct params_type {
  int param[LIBMICSMM_NPARAMS];
};


template<typename T>
class stream_t {
public:
  stream_t()
    : m_stream(0)
    , m_size_stack(0), m_size_a(0), m_size_b(0), m_size_c(0)
    , m_dev_params(0), m_dev_a(0), m_dev_b(0), m_dev_c(0)
  {}

  explicit stream_t(const char* name)
    : m_stream(0)
    , m_size_stack(0), m_size_a(0), m_size_b(0), m_size_c(0)
    , m_dev_params(0), m_dev_a(0), m_dev_b(0), m_dev_c(0)
  {
    LIBXSTREAM_CHECK_CALL_THROW(acc_stream_create(&m_stream, name, -1));
  }
  
  ~stream_t() {
    if (m_stream) {
      deinit();
      acc_stream_destroy(m_stream);
    }
  }

public:
  operator const void*() const { return m_stream; }

  void swap(stream_t& other) throw() {
    std::swap(m_stream, other.m_stream);
    std::swap(m_size_stack, other.m_size_stack);
    std::swap(m_size_a, other.m_size_a);
    std::swap(m_size_b, other.m_size_b);
    std::swap(m_size_c, other.m_size_c);
    std::swap(m_dev_params, other.m_dev_params);
    std::swap(m_dev_a, other.m_dev_a);
    std::swap(m_dev_b, other.m_dev_b);
    std::swap(m_dev_c, other.m_dev_c);
    m_hst_params.swap(other.m_hst_params);
    m_hst_a.swap(other.m_hst_a);
    m_hst_b.swap(other.m_hst_b);
    m_hst_c.swap(other.m_hst_c);
  }

  int operator()(const char* groupname, size_t id, double& max_diff) {
    size_t size_stack = 0, size_a = 0, size_b = 0, size_c = 0;
    int homogeneous = 0, max_m = 0, max_n = 0, max_k = 0;
    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "stack", id, 0, &size_stack, &homogeneous));
    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "adata", id, 0, &size_a, &max_m));
    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "bdata", id, 0, &size_b, &max_k));
    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "cdata", id, 0, &size_c, &max_n));

    if (this->not_ready(size_stack, size_a, size_b, size_c)) {
#if defined(_OPENMP)
#     pragma omp critical
#endif
      if (this->not_ready(size_stack, size_a, size_b, size_c)) {
        LIBXSTREAM_CHECK_CALL(this->deinit());
        LIBXSTREAM_CHECK_CALL(this->init(size_stack, size_a, size_b, size_c));
      }
    }

    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "stack", id, &m_hst_params[0], 0, 0));
    LIBXSTREAM_CHECK_CALL(acc_memcpy_h2d(&m_hst_params[0], m_dev_params, size_stack, m_stream));

    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "adata", id, &m_hst_a[0], 0, 0));
    LIBXSTREAM_CHECK_CALL(acc_memcpy_h2d(&m_hst_a[0], m_dev_a, size_a, m_stream));

    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "bdata", id, &m_hst_b[0], 0, 0));
    LIBXSTREAM_CHECK_CALL(acc_memcpy_h2d(&m_hst_b[0], m_dev_b, size_b, m_stream));

    LIBXSTREAM_CHECK_CALL(libsmm_acc_file_load(groupname, "cdata", id, &m_hst_c[0], 0, 0));
    LIBXSTREAM_CHECK_CALL(acc_memcpy_h2d(&m_hst_c[0], m_dev_c, size_c, m_stream));

    LIBXSTREAM_CHECK_CALL(libsmm_acc_process(m_dev_params, m_hst_params.size(), LIBMICSMM_NPARAMS,
      dbcsr_elem<T>::type, m_dev_a, m_dev_b, m_dev_c, max_m, max_n, max_k, homogeneous, m_stream));
    LIBXSTREAM_CHECK_CALL(acc_memcpy_d2h(m_dev_c, &m_hst_c[0], size_c, m_stream));

    // synchronize stream in order to materialize the result
    LIBXSTREAM_CHECK_CALL(acc_stream_sync(m_stream));

    return libsmm_acc_file_diff(groupname, "cgold", id, &m_hst_c[0], dbcsr_elem<T>::type, &max_diff);
  }

private:
  stream_t(const stream_t&);
  stream_t& operator=(const stream_t&);

  bool not_ready(size_t size_stack, size_t size_a, size_t size_b, size_t size_c) const {
    return m_size_stack < size_stack || m_size_a < size_a || m_size_b < size_b || m_size_c < size_c;
  }

  int init(size_t size_stack, size_t size_a, size_t size_b, size_t size_c) {
    LIBXSTREAM_ASSERT(sizeof(params_type) == LIBMICSMM_NPARAMS * sizeof(int));
    m_hst_params.resize(size_stack / (LIBMICSMM_NPARAMS * sizeof(int)));
    LIBXSTREAM_ASSERT(m_hst_params.size() * LIBMICSMM_NPARAMS * sizeof(int) == size_stack);
    m_hst_a.resize(size_a / sizeof(T));
    LIBXSTREAM_ASSERT(m_hst_a.size() * sizeof(T) == size_a);
    m_hst_b.resize(size_b / sizeof(T));
    LIBXSTREAM_ASSERT(m_hst_b.size() * sizeof(T) == size_b);
    m_hst_c.resize(size_c / sizeof(T));
    LIBXSTREAM_ASSERT(m_hst_c.size() * sizeof(T) == size_c);
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_allocate(&m_dev_params, m_hst_params.size() * sizeof(params_type)));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_allocate(&m_dev_a, size_a));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_allocate(&m_dev_b, size_b));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_allocate(&m_dev_c, size_c));
    m_size_stack = size_stack;
    m_size_a = size_a;
    m_size_b = size_b;
    m_size_c = size_c;
    return LIBXSTREAM_ERROR_NONE;
  }

  int deinit() {
    if (0 != m_dev_params || 0 != m_dev_a || 0 != m_dev_b || 0 != m_dev_c) {
      LIBXSTREAM_CHECK_CALL(acc_stream_sync(m_stream));
    }
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_deallocate(m_dev_params));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_deallocate(m_dev_c));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_deallocate(m_dev_b));
    LIBXSTREAM_CHECK_CALL(acc_dev_mem_deallocate(m_dev_a));
    return LIBXSTREAM_ERROR_NONE;
  }

private:
  void* m_stream;
  size_t m_size_stack, m_size_a, m_size_b, m_size_c;
  void *m_dev_params, *m_dev_a, *m_dev_b, *m_dev_c;
  std::vector<params_type> m_hst_params;
  std::vector<T> m_hst_a, m_hst_b, m_hst_c;
};


int main(int argc, char* argv[])
{
  typedef double T;
  typedef stream_t<T> stream_type;

  try {
    const char *const groupname = 1 < argc ? argv[1] : getenv("LIBMICSMM_DUMP");
    const int end = std::max(2 < argc ? std::atoi(argv[2]) : std::numeric_limits<int>::max(), 1);

    int ndevices = 0;
    LIBXSTREAM_CHECK_ERROR(acc_get_ndevices(&ndevices));
    LIBXSTREAM_CHECK_CONDITION(0 < ndevices);

    stream_type stream[LIBXSTREAM_MAX_NSTREAMS];
    double max_diff = 0;
    size_t size = 0;

#if defined(_OPENMP)
#   pragma omp parallel
#endif
    for (size_t i = 0; i < end && LIBXSTREAM_ERROR_NONE == libsmm_acc_file_load(groupname, "stack", i, 0, &size, 0); ++i) {
#if defined(_OPENMP)
#     pragma omp single nowait
#endif
      {
#if defined(_OPENMP) && (200203 < _OPENMP)
#       pragma omp task
        {
          const size_t sid = omp_get_thread_num() % LIBXSTREAM_MAX_NSTREAMS;
#else
          const size_t sid = 0;
#endif
          stream_type& current_stream = stream[sid];

          if (0 == current_stream) {
#if defined(_OPENMP)
#           pragma omp critical
#endif
            if (0 == current_stream) {
              LIBXSTREAM_CHECK_CALL_THROW(acc_set_active_device(sid % ndevices));
              char name[128];
              LIBXSTREAM_SNPRINTF(name, sizeof(name), "Stream %d", sid + 1);
              stream_type(name).swap(current_stream);
            }
          }

          LIBXSTREAM_CHECK_CALL_THROW(current_stream(groupname, i, max_diff));

#if defined(_OPENMP) && (200203 < _OPENMP)
        }
#endif
      }
    }

    if (0 == size) throw std::runtime_error("Failed to load Gold data set!");
    fprintf(stdout, "diff = %f\n", max_diff);
  }
  catch(const std::exception& e) {
    fprintf(stderr, "Error: %s\n", e.what());
    return EXIT_FAILURE;
  }
  catch(...) {
    fprintf(stderr, "Error: unknown exception caught!\n");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
#endif // defined(LIBMICSMM_USE_STANDALONE) && !defined(LIBMICSMM_USE_DUMP)

#endif // defined(__ACC) && defined(__ACC_MIC) && defined(__DBCSR_ACC)
