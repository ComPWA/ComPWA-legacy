[![image](https://travis-ci.com/ComPWA/ComPWA.svg?branch=master)](https://travis-ci.com/ComPWA/ComPWA)
[![image](https://codecov.io/gh/ComPWA/ComPWA/branch/master/graph/badge.svg)](https://codecov.io/gh/ComPWA/ComPWA)
[![image](https://www.codefactor.io/repository/github/compwa/compwa/badge)](https://www.codefactor.io/repository/github/compwa/compwa)
[![image](https://scan.coverity.com/projects/13697/badge.svg)](https://scan.coverity.com/projects/compwa-compwa)

![image](https://github.com/ComPWA/ComPWA/blob/master/doc/images/logo.png)

*For ComPWA's **Python interface** (`pycompwa`), see [our website]([compwa.github.io/ComPWA/](https://compwa.github.io/)).*

# About

ComPWA is a project to provide a flexible and modular Partial Wave Analysis
framework for various use-cases. It was initiated and is developed by the Panda
Collaboration (antiProton ANnihilation at DArmstadt) at Fair (Facility for
Antiproton and Ion Research). But ComPWA will not just be used for the Panda
physics program, but also for various other experiments to provide a commonly
used tool which is stable, efficient and provides comparable results. At the
moment there are many PWA-tools on the market, most used just for specific
experiments and specific physics cases, some experiments even have multiple
tools. But why write the same software again and again? E.g. the model
describing physical processes should stay the same independent where and how
there was a measurement of the process. Using the actual same implementation of
the model does not only save a lot of time, it also ensures that two experiments
are comparing the same thing. The same argument holds for optimization-routines
and estimation-functions. It might even allow combined fitting of different
experiments instead of taking the average of the results\!

The natural modularization, following the considerations above, would be to
separate into experiment specific information, physics (models, formalisms),
estimation how good the model fits the data and optimization of free parameters.
The first considerations on this where discussed with experts from different
experiments and different technologies where discussed and tested. The result of
this process is the first requirement document of the new PWA-Framework. This
sketch illustrates the modular concept:

![image](https://github.com/ComPWA/ComPWA/blob/master/doc/images/compwa_modules.png)

# Available Features:

- **Computational Backends:**
  - [x] FunctionTree
  - [ ] Interface to Tensorflow
- **Physics Models:**
  - [x] Helicity Formalism
  - [x] Canonical Formalism
  - [ ] K-Matrix Formalism
- **Data:**
  - [x] Input/Output from/to ROOT files
  - [ ] Input/Output from/to Ascii files
  - [x] Direct IO via python
        (see [pycompwa](https://github.com/ComPWA/pycompwa))
  - [x] Data Generation
    - [x] ROOT phase space generator
    - [x] EvtGen phase space generator
    - [x] Model-based hit&miss Monte Carlo data generation
    - [x] Model-based importance sampled data generation
- **Estimators:**
  - [x] Unbinned Maximum Likelihood
- **Optimizers:**
  - [x] Interface to Minuit2
  - [x] Interface to Geneva
- **User Interfaces/Steering:**
  - [x] C++
  - [x] Python, see [pycompwa](https://github.com/ComPWA/pycompwa)

**ComPWA offers the python interface
[pycompwa](https://github.com/ComPWA/pycompwa) which gives access to most
functionality of ComPWA and provides an expert-system to generate amplitude
models.**

# Prerequisites

ComPWA is supposed to run on most modern unix systems (including MacOS).
The following packages are mandatory:

- git (optional, for easier updates and if you want to contribute)
- cmake ( \> 3.3 )
- gcc (\> 5.1) or clang
- [Boost](http://www.boost.org/users/download/), version \>= 1.54

For a more feature rich installation, the following packages are highly
recommended:

- python3 (requirement of
  [pycompwa](https://github.com/ComPWA/pycompwa))
- [ROOT](http://root.cern.ch/drupal/content/downloading-root),
  version 5.34, 6.08
- [Minuit2](http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/),
  version 5.34, 6.08
- [Geneva](https://launchpad.net/geneva/+download), version 1.8.0

In case that some dependencies are not met on your system use your package
manager to install them or use the manual procedures described below. For MacOS
you can use e.g. [MacPorts](https://www.macports.org) as package manager. You
can also try to use different (newer) versions than the ones stated above,
however those are not tested.

In order to install all dependencies it is probably also useful to have a look
on the
[installation instructions file](https://github.com/ComPWA/ComPWA/blob/master/.travis.yml)
for TravisCI.

# Quick Installation

A detailed guide can be found below. The installation basically boils down to:

```shell
git clone --recurse-submodules https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>
cd <COMPWA_SOURCE_PATH>
mkdir build && cd build
cmake ..
cmake --build .
```

# Examples

The repository contains a couple of
[examples](https://github.com/ComPWA/ComPWA/tree/master/Examples).To learn about
more detailed features of ComPWA you also might want to have a look on the
examples of the
[pycompwa package](https://github.com/ComPWA/pycompwa/tree/master/examples/jupyter).

# Documentation

Source code documentation via Doxygen is provided
[here](https://compwa.github.io/ComPWA/). The master branch is automatically
built using TravisCI. Probably, it is interesting to check out the
[log file](https://travis-ci.com/ComPWA/ComPWA) and the projects TravisCI
configuration file
[travisCI.yml](https://github.com/ComPWA/ComPWA/blob/master/.travis.yml).

# ComPWA installation

## Manual installation of dependencies

- **Boost**:
  to install Boost follow
  [these](http://www.boost.org/doc/libs/1_54_0/more/getting_started/unix-variants.html#easy-build-and-install)
  instructions. The `--prefix=path/to/installation/prefix` option is useful as
  you can use the same path to point ComPWA to this specific Boost installation.

- **ROOT**:
  to install Root follow
  [these](https://root.cern.ch/building-root) instructions.

- **Minuit2**:
  is included in most ROOT installations. In case you want to install a
  stand-alone version follow
  [these](http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/gettingStarted/autoconf.html)
  instructions. In addition you should use
  `./configure --prefix=/your/minuit2/installation/path` followed by
  `make install` to specify an installation directory which ComPWA later needs
  to find the library and to install all needed files to this location.

- **Geneva**:
  to install Geneva follow
  [these](http://www.gemfony.eu/index.php?id=genevainstallation)
  instructions. The `DCMAKE_INSTALL_PREFIX="/where/the/library/goes"` option is
  useful as you can use the same path to point ComPWA to this specific Geneva
  installation:

  ```shell
  cd GENEVA_SOURCE
  mkdir -p build/install
  cd build
  cmake ../ -DCMAKE_INSTALL_PREFIX=./install
  cmake --build .
  make install
  cp install/CMakeModules/FindGeneva.cmake YOUR_COMPWA_PATH/cmake/Modules/
  export GENEVA_ROOT=YOUR_GENEVA_PATH/build/install
  ```

  - Note for Fedora 25:
    The Geneva tests are build by default but might have trouble finding the
    boost test libraries of the Fedora boost package. A workaround is to disable
    them within `YOUR_GENEVA_PATH/CMakeModules/CommonGenevaBuild.cmake, line 55`
    (replace the line with `SET( GENEVA_BUILD_TESTS FALSE )`.
  - Alternatively you can follow the instructions from the Geneva
    [manual](http://www.gemfony.eu/fileadmin/documentation/geneva-manual.pdf).

## Getting ComPWA

Get the most recent version:

```shell
git clone --recurse-submodules git@github.com:ComPWA/ComPWA <COMPWA_SOURCE_PATH>
```

This will clone the repository and its sub-modules to the sub-folder
`<COMPWA_SOURCE_PATH>` within the current directory. For multi-threading ComPWA
uses the parallel stl algorithms of c++17. Unfortunately the current compilers
do not have any implementations for this. Here ComPWA currently relies on
[TBB](https://github.com/01org/tbb) and
[parallelstl](https://github.com/intel/parallelstl), which are included in
ComPWA as git submodules.

## Building ComPWA

ComPWA uses `cmake` as build system. The usual steps to build all libraries and
the test executable are the following:

- Create and enter a build folder (preferably not the ComPWA source folder)

  ```shell
  mkdir build
  cd build
  ```

- Set your compiler if you do not use the system default compiler

  ```shell
  export CC=<path_to_your_compiler>
  export CXX=<path_to_your_compiler>
  ```

- Build the project. You may leave the `DCMAKE_INSTALL_PREFIX` empty, but then
  ComPWA will be installed system wide.

  ```shell
  cmake .. -DCMAKE_INSTALL_PREFIX=<COMPWA_INSTALL_PATH>
  cmake --build .
  make install      # optional
  ctest -C debug    # optional: run test suite
  ```

- You might want to create a pre-configured project for an IDE (e.g. vscode,
  [eclipse](https://www.eclipse.org), Xcode) via:

  ```shell
  cmake -G"Eclipse CDT4 - Unix Makefiles" ../<COMPWA_SOURCE_PATH>
  ```

## Installation via Docker

A [Dockerfile](https://github.com/ComPWA/ComPWA/blob/master/Dockerfile) for
ComPWA is provided. You can use it to build an [docker](https://www.docker.com)
image to run ComPWA. Using such an image ComPWA should run on
[all systems that are supported by docker](https://docs.docker.com/engine/installation/)
including several (commercial) cloud computing services. If you are new to
docker you can have a look on [this](https://prakhar.me/docker-curriculum/)
tutorial.

## System specific notes

### HimsterII / Mogon II

[Mogon2](https://hpc.uni-mainz.de/) is the supercomputer of the Mainz
University. If you work on it, you can fulfill the ComPWA
[installation requirements](#requirements) by loading a series of modules:

```shell
module load devel/CMake/3.9.5
module load toolchain/foss/2017a
module load devel/Boost/1.65.1-foss-2017a
module load ROOT/v6.12-foss-2017a-python3
export CC=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/gcc
export CXX=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/g++
```

Now follow: [Building ComPWA](#building-compwa).
