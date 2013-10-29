#!/bin/bash

export COMPWA_DIR=/Users/weidenka/work/rootAnalysis/ComPWA
export PWA_LIBS=$COMPWA_DIR/lib
export EXTERNAL_DIR=/Users/weidenka/work/external

export QFT=$EXTERNAL_DIR/qft++
export GENEVA=
export MINUIT2=$EXTERNAL_DIR/Minuit2-5.28.00/build
#BOOST installation root directory
export BOOST_ROOT=$EXTERNAL_DIR/boost_1_54_0/build
#BOOST.build root directory (!)
export BOOST_BUILD_PATH=$BOOST_ROOT/../tools/build/v2
export ROOTSYS=$EXTERNAL_DIR/root-v5.34.09

export PATH=$PATH:$BOOST_ROOT/bin:$ROOTSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENEVA/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_ROOT/lib:$PWA_LIBS:$QFT/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Users/weidenka/work/rootAnalysis/DKsKK/lib
export DYLD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENEVA/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_ROOT/lib:$PWA_LIBS:$QFT/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/Users/weidenka/work/rootAnalysis/DKsKK/lib

#echo "check libs:"
#echo $LD_LIBRARY_PATH

#prevents minuit vom using multiple threads
#needed as long as data-reader is not thread-safe
export OMP_NUM_THREADS=1
