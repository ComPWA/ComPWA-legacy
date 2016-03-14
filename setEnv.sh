#!/bin/bash

export EXTERNAL_DIR=/home/steve/compwa_externals
export COMPWA_DIR=/home/steve/ComPWA

export QFT=$EXTERNAL_DIR/qft++
#export GENEVA=$EXTERNAL_DIR/geneva_1_4_0
export MINUIT2=$EXTERNAL_DIR/Minuit2-5.34.14
#BOOST installation root directory
export BOOST_ROOT=$EXTERNAL_DIR/boost_1_55_0
#BOOST.build root directory (!)
export BOOST_BUILD_PATH=$BOOST_ROOT/share/boost-build
export CLIPS_ROOT=$EXTERNAL_DIR/clips_630

if [ ! $ROOTSYS ]; then
  export ROOTSYS=/cvmfs/fairroot.gsi.de/fairsoft/apr13
fi

export PATH=$BOOST_ROOT/bin:$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$CLIPS_ROOT/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_ROOT/lib:$QFT/lib:$COMPWA_DIR/Tools/RootNeatPlotting/build:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$COMPWA_DIR/lib:$LD_LIBRARY_PATH

#prevents minuit vom using multiple threads
#needed as long as data-reader is not thread-safe
export OMP_NUM_THREADS=1
