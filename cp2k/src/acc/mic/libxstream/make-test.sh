#!/bin/bash

CXX=$(which icpc 2> /dev/null)

ICCOPT="-O2 -xHost -ansi-alias -DNDEBUG"
GCCOPT="-O2 -march=native -DNDEBUG"

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

$CXX -std=c++0x $OPT $* -pthread \
  -D__ACC -D__ACC_MIC -DLIBXSTREAM_TEST -DLIBXSTREAM_TEST_STANDALONE \
  -Iinclude src/*.cpp \
  -o test
