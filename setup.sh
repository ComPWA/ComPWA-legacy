
#module load root/5.32.00
export PATH=$PATH:$BOOST_ROOT

export COMPWA_DIR=/data/work/kpha1/weidenka/ComPWA
export PWA_LIBS=$COMPWA_DIR/lib
export EXTERNAL_DIR=$COMPWA_DIR/external

export QFT_LIB=$EXTERNAL_DIR/qft++/lib

export GENEVA=

export MINUIT2=$EXTERNAL_DIR/Minuit2-5.28.00/build

export BOOST_ROOT=$EXTERNAL_DIR/boost_1_48_0/
export BOOST_BUILD_PATH=$BOOST_ROOT/build

export ROOTSYS=/home/weidenka/opt/
export ROOTINCLUDE=$ROOTSYS/include/root

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GENEVA/lib:$ROOTSYS/lib/root:$MINUIT2/lib:$BOOST_BUILD_PATH/lib:$PWA_LIBS:$QFT_LIB

#echo "check libs:"
#echo $LD_LIBRARY_PATH

#prevents minuit vom using multiple threads
#needed as long as data-reader is not thread-safe
export OMP_NUM_THREADS=1
