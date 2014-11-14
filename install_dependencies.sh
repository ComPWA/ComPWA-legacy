#!/bin/bash

#install boost --------------------------------------------------------------------------
function install_boost {
cd ${temp_workdir}

wget http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download -O ${temp_workdir}/boost_1_55_0.tar.gz

tar xf boost_1_55_0.tar.gz

cd boost_1_55_0
boost_root=${install_path}/boost_1_55_0
mkdir ${boost_root}

./bootstrap.sh --prefix=${boost_root}
./b2 -j ${numproc} cxxflags=-std=c++0x install

#install boost_build
cd tools/build/v2
./bootstrap.sh
./b2 install --prefix=${boost_root}
}

#install geneva --------------------------------------------------------------------------
function install_geneva {
cd ${temp_workdir}

wget https://launchpad.net/geneva/trunk/1.4/+download/geneva-v1.4.0.tgz -O ${temp_workdir}/geneva-v1.4.0.tgz

tar xf geneva-v1.4.0.tgz

cd geneva-v1.4.0
geneva_root=${install_path}/geneva_1_4_0

cmake . -DCMAKE_CXX_COMPILER=${compiler_path} -DCMAKE_INSTALL_PREFIX="${geneva_root}" -DBOOST_ROOT="${boost_root}" -DBOOST_INCLUDEDIR="${boost_root}/include/boost" -DBOOST_LIBRARYDIR="${boost_root}/lib"
make -j$numproc
make install
}

#install root --------------------------------------------------------------------------
function install_root {
cd ${temp_workdir}

wget ftp://root.cern.ch/root/root_v5.34.23.source.tar.gz -O ${temp_workdir}/root_v5.34.23.source.tar.gz

tar xf root_v5.34.23.source.tar.gz

root_root=${install_path}/root_v5.34.23
cd root
 
./configure --prefix=${root_root} --etcdir=${root_root}/etc --enable-cxx11 --enable-fftw3 --enable-mathmore --enable-minuit2 --enable-python --enable-roofit --enable-tmva

make -j{$numproc}
make install
}

#install qft++ --------------------------------------------------------------------------
function install_qftpp {
cd ${temp_workdir}

wget http://www-meg.phys.cmu.edu/williams/qft++/qft++.v1.tar.gz -O ${temp_workdir}/qft++.v1.tar.gz

tar xf qft++.v1.tar.gz

cd qft++/src
wget --no-check-certificate https://raw.githubusercontent.com/ComPWA/ComPWA/gh-pages/patches/qft++.patch -O ${temp_workdir}/qft++/src/qft++.patch
patch -p0 < qft++.patch
rm qft++.patch
make

#install
cd ${temp_workdir}
mv qft++ ${install_path}/.
}

#install minuit2 --------------------------------------------------------------------------
function install_minuit2 {
cd ${temp_workdir}

wget http://www.cern.ch/mathlibs/sw/5_34_14/Minuit2/Minuit2-5.34.14.tar.gz -O ${temp_workdir}/Minuit2-5.34.14.tar.gz

tar xf Minuit2-5.34.14.tar.gz

cd Minuit2-5.34.14

./configure --prefix=${install_path}/Minuit2-5.34.14
make
make install
}
