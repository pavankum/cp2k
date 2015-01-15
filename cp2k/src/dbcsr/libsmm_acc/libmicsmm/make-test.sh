#!/bin/bash

CXX=$(which icpc 2> /dev/null)

ICCOPT="-O2 -xHost -ansi-alias -DNDEBUG"
GCCOPT="-O2 -march=native -DNDEBUG"
OUT="test"

if [ "" = "$CXX" ] ; then
  OPT=$GCCOPT
  CXX="g++"
else
  OPT=$ICCOPT
fi

if [ "-g" = "$1" ] ; then
  OPT="-O0 -g"
  shift
fi

if [ "$LIBXSMMROOT" != "" ] ; then
  OPT+=" -D__LIBXSMM -I$LIBXSMMROOT/include"
  OUT+=" $LIBXSMMROOT/lib/intel64/libxsmm.a -mkl=sequential"
fi

$CXX -std=c++0x $OPT $* -D__ACC -D__ACC_MIC -D__DBCSR_ACC -DLIBMICSMM_USE_STANDALONE \
  -I../../../acc/mic/libxstream/include -I. \
  ../../../acc/mic/libxstream/src/*.cpp \
  ../../../acc/mic/*.c \
  *.cpp  \
  -o $OUT

