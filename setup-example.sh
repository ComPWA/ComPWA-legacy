#!/bin/bash

export COMPWA_DIR=$PWD
export PWA_LIBS=$COMPWA_DIR/lib
export EXTERNAL_DIR=/home/pflueger/ComPWA_externals

export QFT=$EXTERNAL_DIR/qft++
export GENEVA=$EXTERNAL_DIR/geneva_1_4_0
export MINUIT2=$EXTERNAL_DIR/Minuit2-5.34.14
#BOOST installation root directory
#export BOOST_ROOT=/home/weidenka/external/boost_1_54_0/build
export BOOST_ROOT=$EXTERNAL_DIR/boost_1_55_0
#BOOST.build root directory (!)
export BOOST_BUILD_PATH=$BOOST_ROOT/share/boost-build
export ROOTSYS=$EXTERNAL_DIR/root_v5.34.23

export PATH=$BOOST_ROOT/bin:$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$GENEVA/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_ROOT/lib:$PWA_LIBS:$QFT/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/weidenka/work/rootAnalysis/DKsKK/lib
#export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENEVA/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_ROOT/lib:$PWA_LIBS:$QFT/lib
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/weidenka/work/rootAnalysis/DKsKK/lib

#echo "check libs:"
#echo $LD_LIBRARY_PATH

#prevents minuit vom using multiple threads
#needed as long as data-reader is not thread-safe
export OMP_NUM_THREADS=1
