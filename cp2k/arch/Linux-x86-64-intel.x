# Arch file targeting Linux 64-bit using the Intel tool chain
#
PERL     = perl
CPP      = #cpp
AR       = xiar -r

# LIBINT: configure, build, and install
# Version 1.1.5 (tested)
#
# $ env \
#   AR=xiar CC=icc CXX=icpc \
#   ./configure --prefix=$HOME/libint \
#     --with-cc-optflags="-O1 -xHost" \
#     --with-cxx-optflags="-O1 -xHost" \
#     --with-libint-max-am=5 \
#     --with-libderiv-max-am1=4
# $ make
# $ make install
# $ make realclean
#
#LIBINTROOT = $(HOME)/libint

# LIBXC: configure, build, and install
# Version 2.2.2 (tested)
#
# $ autoreconf --force --install
# $ env \
#   AR=xiar \
#   FC=ifort F77=ifort F90=ifort FCFLAGS="-O2 -xHost" \
#   CC=icc CFLAGS="-O2 -xHost" \
#   ./configure --prefix=$HOME/libxc
# $ make
# $ make install
# $ make clean
#
#LIBXCROOT = $(HOME)/libxc

# LIBXSMM: cp2k/tools/build_libsmm/libxsmm or https://github.com/hfp/libxsmm
# Please note that CP2K redistributes a tested version of LIBXSMM
# which is built automatically (if LIBXSMMROOT is set).
#
LIBXSMMROOT = $(TOOLSRC)/build_libsmm/libxsmm

# LIBXSTREAM: cp2k/tools/mic/libxstream or https://github.com/hfp/libxstream
# Please note that CP2K redistributes a tested version of LIBXSTREAM
# which is built automatically (if LIBXSTREAMROOT is set).
#
LIBXSTREAMROOT = $(TOOLSRC)/mic/libxstream

# Diagnostic message to be turned off
DIAG_DISABLE = 8290,8291,10010,10212,11060

# DEFAULTS
#
BEEP ?= 1
JIT ?= 0
SSE ?= 0
AVX ?= 0
MPI ?= 1
OMP ?= 1
ACC ?= 0
OCL ?= 0
SYM ?= 0
DBG ?= 0
IPO ?= 0
MKL ?= 1
MKL_DIRECT ?= 0
MKL_STATIC ?= 1
RECONFIGURE ?= 1
TBBMALLOC ?= 1
MEMKIND ?= 1
NESTED ?= 0

OFFLOAD ?= $(ACC)
ifeq ($(OFFLOAD),0)
  MIC ?= 0
else
  JIT = 0
endif
ifneq ($(MIC),0)
  JIT = 0
endif

ifeq (1,$(shell echo $$((2 > $(DBG)))))
  ifeq (1,$(AVX))
    TARGET = -xAVX
  else ifeq (2,$(AVX))
    TARGET = -xCORE-AVX2
  else ifeq (3,$(AVX))
    ifeq (0,$(MIC))
      TARGET = -xCOMMON-AVX512
    else
      TARGET = -xMIC-AVX512
    endif
  else ifeq (1,$(shell echo $$((2 <= $(SSE)))))
    TARGET = -xSSE$(SSE)
  else ifeq (1,$(SSE))
    TARGET = -xSSE3
  else
    TARGET = -xHost
  endif
endif

# initial build flags
CPPFLAGS  = $(NULL)
CXXFLAGS  = -std=c++0x
CFLAGS    = #
FCFLAGS   = -free -fpp -heap-arrays
LDFLAGS   = #
OPTFLAGS  = $(TARGET)

# workaround for issue "cannot find address of function"
#ATTRIBUTE = mic
#DIAG_DISABLE := $(DIAG_DISABLE),2571,3218

ifeq (0,$(DBG))
  OPTFLAGS  += -O2
  DFLAGS    += -DNDEBUG

  CXXFLAGS  += -fno-alias -ansi-alias #-fp-model fast=2 #precise
  CFLAGS    += -fno-alias -ansi-alias #-fp-model fast=2 #precise
  FCFLAGS   += #-fp-model fast=2 #source
  LDFLAGS   += #
else
  OPTFLAGS  += -O0
  ifeq (2,$(DBG))
    FCFLAGS   += -fpe0 # debugging NaNs
  endif
  SYM = $(DBG)
endif

ifneq (0,$(IPO))
  OPTFLAGS += -ipo-separate
else ifeq (0,$(IPO))
  LDFLAGS += -no-ipo
endif

ifneq (0,$(SYM))
  DFLAGS += -D__USE_CP2K_TRACE
  OPTFLAGS  += -traceback
  ifneq (1,$(SYM))
    CXXFLAGS := -g3 -gdwarf-2 -debug $(CXXFLAGS)
    CFLAGS := -g3 -gdwarf-2 -debug $(CFLAGS)
    FCFLAGS := -g -debug $(FCFLAGS)
  else
    CXXFLAGS := -g -debug $(CXXFLAGS)
    CFLAGS := -g -debug $(CFLAGS)
    FCFLAGS := -g -debug $(FCFLAGS)
  endif
endif

ifneq (0,$(MPI))
  CXX = mpiicpc
  CC  = mpiicc
  FC  = mpiifort
  LD  = mpiifort
  DFLAGS += -D__parallel -D__BLACS -D__SCALAPACK
  #DFLAGS += -D__SCALAPACK2
  ifneq (1,$(MPI))
    DFLAGS += -D__MPI_VERSION=$(MPI)
  endif
else
  CXX = icpc
  CC  = icc
  FC  = ifort
  LD  = ifort
endif

ifneq (0,$(OMP))
  FCFLAGS   += -threads
  LDFLAGS   += -threads
  OPTFLAGS  += -openmp
  ifneq (0,$(NESTED))
    DFLAGS += -D__NESTED_OPENMP
  endif
endif

ifneq (,$(LIBINTROOT))
  DFLAGS  += -D__LIBINT -D__LIBINT_MAX_AM=6 -D__LIBDERIV_MAX_AM1=5
  IFLAGS  += -I$(LIBINTROOT)/include
  LIBS    += $(LIBINTROOT)/lib/libderiv.a $(LIBINTROOT)/lib/libint.a
endif

ifneq (,$(LIBXCROOT))
  DFLAGS  += -D__LIBXC2
  IFLAGS  += -I$(LIBXCROOT)/include
  LIBS    += $(LIBXCROOT)/lib/libxcf90.a $(LIBXCROOT)/lib/libxc.a
endif

ifneq (,$(ELPAROOT))
  ifneq (1,$(ELPA)) # default is ELPA2
    DFLAGS  += -D__ELPA2
    IFLAGS  += -I$(ELPAROOT)/include/elpa/modules
    LIBS    += $(ELPAROOT)/lib/libelpa.a
  else ifeq (1,$(ELPA))
    DFLAGS  += -D__ELPA
    IFLAGS  += -I$(ELPAROOT)/include/elpa/modules
    LIBS    += $(ELPAROOT)/lib/libelpa.a
  endif
endif

ifneq (0,$(TBBMALLOC))
  ifneq (,$(TBBROOT))
    GCC = $(notdir $(shell which gcc 2> /dev/null))
    ifneq (,$(GCC))
      GCC_VERSION = $(shell $(GCC) --version | grep "gcc (GCC)" | sed "s/gcc (GCC) \([0-9]\+\.[0-9]\+\.[0-9]\+\).*$$/\1/")
      GCC_VERSION_MAJOR = $(shell echo "$(VERSION)" | cut -d"." -f1)
      GCC_VERSION_MINOR = $(shell echo "$(VERSION)" | cut -d"." -f2)
      GCC_VERSION_PATCH = $(shell echo "$(VERSION)" | cut -d"." -f3)
      TBBMALLOCLIB = $(wildcard $(TBBROOT)/lib/intel64/gcc$(GCC_VERSION_MAJOR).$(GCC_VERSION_MINOR)/libtbbmalloc_proxy.so)
    endif
    ifeq (,$(TBBMALLOCLIB))
      TBBGCCDIR = $(shell ls -1 "$(TBBROOT)/lib/intel64" | tr "\n" " " | rev | cut -d" " -f2 | rev)
      TBBMALLOCLIB = $(wildcard $(TBBROOT)/lib/intel64/$(TBBGCCDIR)/libtbbmalloc_proxy.so)
    endif
    ifneq (,$(TBBMALLOCLIB))
      LIBS += -L$(dir $(TBBMALLOCLIB)) -ltbbmalloc_proxy
      #FCFLAGS += -heap-arrays
    endif
  endif
else ifneq (,$(TCMALLOCROOT))
  # configured using ./configure --enable-minimal --prefix=<TCMALLOCROOT>
  LIBS += $(TCMALLOCROOT)/lib/libtcmalloc_minimal.a
endif

ifneq (0,$(MEMKIND))
  ifneq (,$(MEMKINDROOT))
    #LIBS += -L$(MEMKINDROOT)/lib -lmemkind
    LIBS += $(MEMKINDROOT)/lib/libmemkind.a
  endif
endif

ifneq (,$(LIBXSMMROOT))
  ifneq (-1,$(JIT))
    LIBXSMM_MNK := " \
      23, \
      6, \
      14 16 29, \
      14 32 29, \
      5 32 13 24 26, \
      9 32 22, \
      64, \
      78, \
      16 29 55, \
      32 29 55, \
      12, \
      4 5 7 9 13 25 26 28 32 45"
  endif
  LIBXSMM_ALIGNED_STORES := 0
  ifneq (0,$(OMP))
    ifneq (0,$(NESTED))
      LIBXSMM_ALIGNED_STORES := 1
    endif
  endif
  ifneq (0,$(ACC))
    ifneq (0,$(OFFLOAD))
      LIBXSMM_ALIGNED_STORES := 1
    endif
  endif
  ifneq (,$(filter %MIC-AVX512,$(TARGET)))
    LIBXSMM_PREFETCH := 1
  else
    LIBXSMM_PREFETCH := 0
  endif
  LIBXSMM_BUILD := $(shell $(MAKE) -f $(LIBXSMMROOT)/Makefile \
    INCDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/include \
    BLDDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/build \
    BINDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/bin \
    OUTDIR=$(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib \
    SYM=$(SYM) DBG=$(DBG) IPO=$(IPO) OFFLOAD=$(OFFLOAD) MIC=$(MIC) \
    ALIGNED_STORES=$(LIBXSMM_ALIGNED_STORES) MNK=$(LIBXSMM_MNK) \
    PREFETCH=$(LIBXSMM_PREFETCH) JIT=$(JIT) \
  >&2)

  DFLAGS  += -D__LIBXSMM
  IFLAGS  += -I$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/include
  LIBS    += $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxsmm/lib/intel64/libxsmm.a
endif

ifneq (0,$(ACC))
  DFLAGS += -D__ACC -D__DBCSR_ACC

  ifeq (0,$(OCL))
    DFLAGS += -D__ACC_MIC
    ifeq (0,$(OFFLOAD))
      OPTFLAGS += -no-offload
      LDFLAGS += -offload-option,mic,ld,"--unresolved-symbols=ignore-all"
    else # also true if OFFLOAD is undefined
      #OPTFLAGS += -offload=mandatory
      # enable OpenMP for OFFLOAD regardless of wether OMP is enabled or not
      MIC_CXFLAGS += -openmp -no-openmp -offload-option,mic,compiler,"-openmp"
      MIC_CCFLAGS += -openmp -no-openmp -offload-option,mic,compiler,"-openmp"
      MIC_FCFLAGS += -openmp -no-openmp -offload-option,mic,compiler,"-openmp"
      MIC_LDFLAGS += -offload-option,mic,ld,"--no-undefined"
      ifneq (,$(ATTRIBUTE))
        MIC_CXFLAGS += -offload-attribute-target=$(ATTRIBUTE)
        MIC_CCFLAGS += -offload-attribute-target=$(ATTRIBUTE)
        #MIC_FCFLAGS += -offload-attribute-target=$(ATTRIBUTE)
      endif
    endif
  else
    DFLAGS  += -D__OPENCL -D__USE_INTEL_CL
    LIBS    += -L/usr/lib64 -lOpenCL -lrt
  endif
else
  ifeq (0,$(OFFLOAD))
    OPTFLAGS += -no-offload
    LDFLAGS += -offload-option,mic,ld,"--unresolved-symbols=ignore-all"
  endif

  # save some build time
  LIBXSTREAMROOT = $(NULL)
endif

ifneq (,$(LIBXSTREAMROOT))
  LIBXSTREAM_BUILD := $(shell $(MAKE) -f $(LIBXSTREAMROOT)/Makefile \
    BLDDIR=$(MAINOBJDIR)/$(ARCH)/$(ONEVERSION)/libxstream \
    OUTDIR=$(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxstream \
    SYM=$(SYM) DBG=$(DBG) IPO=$(IPO) \
    OFFLOAD=$(OFFLOAD) MIC=$(MIC) \
  >&2)

  DFLAGS  += -D__LIBXSTREAM
  IFLAGS  += -I$(LIBXSTREAMROOT)/include
  LIBS    += $(MAINLIBDIR)/$(ARCH)/$(ONEVERSION)/libxstream/libxstream.a
endif

ifeq (1,$(MKL_DIRECT))
  MKL_STATIC = 1
  DFLAGS += -DMKL_DIRECT_CALL_SEQ
endif

ifneq (1,$(MKL))
  ifneq (0,$(MKL)) # smp
    DFLAGS  += -D__MKL -D__FFTSG -D__FFTW3
    IFLAGS  +=-I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
    ifeq (0,$(MKL_STATIC))
      LIBS += -L$(MKLROOT)/lib/intel64
      ifneq (0,$(MPI))
        MIC_LDFLAGS += -offload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64"
        LIBS += -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64
      else
        MIC_LDFLAGS += -offload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread"
        LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread
      endif
    else # static
      ifneq (0,$(MPI))
        MIC_LDFLAGS += -offload-option,mic,ld," \
          --start-group \
            $(MKLROOT)/lib/mic/libmkl_scalapack_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_core.a \
            $(MKLROOT)/lib/mic/libmkl_intel_thread.a \
            $(MKLROOT)/lib/mic/libmkl_blacs_intelmpi_lp64.a \
          --end-group"
        LIBS += \
          -Wl,--start-group \
            $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
            $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a \
            $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
            $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
          -Wl,--end-group
      else
        MIC_LDFLAGS += -offload-option,mic,ld," \
          --start-group \
            $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/mic/libmkl_core.a \
            $(MKLROOT)/lib/mic/libmkl_intel_thread.a \
          --end-group"
        LIBS += \
          -Wl,--start-group \
            $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
            $(MKLROOT)/lib/intel64/libmkl_core.a \
            $(MKLROOT)/lib/intel64/libmkl_intel_thread.a \
          -Wl,--end-group
      endif
    endif
    ifeq (0,$(OMP))
      MIC_LDFLAGS += -offload-option,mic,ld,"-liomp5"
      LIBS += -liomp5
    endif
    MIC_LDFLAGS += -offload-option,mic,ld,"-lpthread -lm"
    LIBS += -lpthread -lm
  endif
else # sequential
  DFLAGS  += -D__MKL -D__FFTSG -D__FFTW3
  IFLAGS  +=-I$(MKLROOT)/include -I$(MKLROOT)/include/fftw
  ifeq (0,$(MKL_STATIC))
    LIBS += -L$(MKLROOT)/lib/intel64
    ifneq (0,$(MPI))
      MIC_LDFLAGS += -offload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64"
      LIBS += -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_lp64
    else
      MIC_LDFLAGS += -offload-option,mic,ld,"-L$(MKLROOT)/lib/mic -lmkl_intel_lp64 -lmkl_core -lmkl_sequential"
      LIBS += -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  else # static
    ifneq (0,$(MPI))
      MIC_LDFLAGS += -offload-option,mic,ld," \
        --start-group \
          $(MKLROOT)/lib/mic/libmkl_scalapack_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_core.a \
          $(MKLROOT)/lib/mic/libmkl_sequential.a \
          $(MKLROOT)/lib/mic/libmkl_blacs_intelmpi_lp64.a \
        --end-group"
      LIBS += \
        -Wl,--start-group \
          $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
          $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/intel64/libmkl_core.a \
          $(MKLROOT)/lib/intel64/libmkl_sequential.a \
          $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
        -Wl,--end-group
    else
      MIC_LDFLAGS += -offload-option,mic,ld," \
        --start-group \
          $(MKLROOT)/lib/mic/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/mic/libmkl_core.a \
          $(MKLROOT)/lib/mic/libmkl_sequential.a \
        --end-group"
      MKL_LIBS = \
        -Wl,--start-group \
          $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
          $(MKLROOT)/lib/intel64/libmkl_core.a \
          $(MKLROOT)/lib/intel64/libmkl_sequential.a \
        -Wl,--end-group
    endif
  endif
  MIC_LDFLAGS += -offload-option,mic,ld,"-lpthread -lm"
  LIBS += -lpthread -lm
endif

ifeq (,$(LIBXSMMROOT))
  ifeq (,$(LIBXSTREAMROOT))
    RECONFIGURE = 0
  else ifeq (0,$(ACC))
    RECONFIGURE = 0
  endif
endif
ifneq (0,$(RECONFIGURE))
  DFLAGS  += -D__RECONFIGURE
  LDFLAGS += -Wl,--wrap=dbcsr_config_mp_dbcsr_set_conf_mm_driver_
  DIAG_DISABLE := $(DIAG_DISABLE),11021
endif

DFLAGS  += -D__INTEL -D__HAS_ISO_C_BINDING
IFLAGS  += # general include paths

# Define __INTEL_COMPILER in case of external preprocessing because some source (pw/fft/fftw3_lib.F)
# toggles code using this symbol, but of course the cpp preprocessor is not defining this symbol.
CPPFLAGS  += #-C $(IFLAGS) $(DFLAGS) -D__INTEL_COMPILER -P -traditional

CXXFLAGS  += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
CFLAGS    += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
FCFLAGS   += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE) $(DFLAGS) $(IFLAGS)
LDFLAGS   += $(OPTFLAGS) -diag-disable $(DIAG_DISABLE)

LIBS += -lstdc++
ifneq (0,$(ACC))
  ifneq (0,$(OFFLOAD))
    LIBS      += $(MIC_LDFLAGS)
    CXXFLAGS  += $(MIC_CXFLAGS)
    CFLAGS    += $(MIC_CCFLAGS)
    FCFLAGS   += $(MIC_FCFLAGS)
    #LDFLAGS   += $(MIC_LDFLAGS)
  endif
endif

ifeq (1,$(shell echo $$((1 <= $(BEEP)))))
mp2_optimize_ri_basis.o: mp2_optimize_ri_basis.F
	$(FC) -c $(FCFLAGS) -O0 $<
qs_vxc_atom.o: qs_vxc_atom.F
	$(FC) -c $(FCFLAGS) -O1 $<
cp_fm_types.o: cp_fm_types.F
	$(FC) -c $(FCFLAGS) -O1 $<
cube_utils.o: cube_utils.F
	$(FC) -c $(FCFLAGS) -O1 $<
ifneq (0,$(OMP))
xc_tpss.o: xc_tpss.F
	$(FC) -c $(filter-out -openmp,$(FCFLAGS)) $<
realspace_grid_types.o: realspace_grid_types.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
matrix_exp.o: matrix_exp.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
cp_dbcsr_operations.o: cp_dbcsr_operations.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
dbcsr_util.o: dbcsr_util.F
	$(FC) -c $(filter-out -heap-arrays,$(FCFLAGS)) $<
endif
endif

ifeq (1,$(shell echo $$((2 <= $(BEEP)))))
ifneq (0,$(OMP))
qs_integrate_potential_product.o: qs_integrate_potential_product.F
	$(FC) -c $(FCFLAGS) -no-openmp $<
dbcsr_work_operations.o: dbcsr_work_operations.F
	$(FC) -c $(FCFLAGS) -no-openmp $<
endif
endif

