.. image:: https://travis-ci.org/ComPWA/ComPWA.svg?branch=master
    :target: https://travis-ci.org/ComPWA/ComPWA

.. image:: https://codecov.io/gh/ComPWA/ComPWA/branch/master/graph/badge.svg 
    :target: https://codecov.io/gh/ComPWA/ComPWA

.. image:: https://www.codefactor.io/repository/github/compwa/compwa/badge 
    :target: https://www.codefactor.io/repository/github/compwa/compwa

.. image:: https://scan.coverity.com/projects/13697/badge.svg
    :target: https://scan.coverity.com/projects/compwa-compwa

.. image:: https://github.com/ComPWA/ComPWA/blob/master/doc/images/logo.png

About
=====

ComPWA is a project to provide a flexible and modular Partial Wave Analysis framework for various use-cases. It was initiated and is developed by the Panda Collaboration (antiProton ANnihilation at DArmstadt) at Fair (Facility for Antiproton and Ion Research). But ComPWA will not just be used for the Panda physics program, but also for various other experiments to provide a commonly used tool which is stable, efficient and provides comparable results. At the moment there are many PWA-tools on the market, most used just for specific experiments and specific physics cases, some experiments even have multiple tools. But why write the same software again and again? E.g. the model describing physical processes should stay the same independent where and how there was a measurement of the process. Using the actual same implementation of the model does not only save a lot of time, it also ensures that two experiments are comparing the same thing. The same argument holds for optimization-routines and estimation-functions. It might even allow combined fitting of different experiments instead of taking the average of the results!

The natural modularization, following the considerations above, would be to separate into experiment specific information, physics (models, formalisms), estimation how good the model fits the data and optimization of free parameters. The first considerations on this where discussed with experts from different experiments and different technologies where discussed and tested. The result of this process is the first requirement document of the new PWA-Framework.
This sketch illustrates the modular concept: 

.. image:: https://github.com/ComPWA/ComPWA/blob/master/doc/images/compwa_modules.png

Available Features:

=========================  ===================================================
Physic Models              Helicity Formalism, Canonical Formalism
Data Formats               ROOT
Estimators                 Chi Square, Unbinned-LogLikelihood
Optimizers                 Miniuit2, Geneva
User Interfaces/Steering   C++, Python
=========================  ===================================================

Prerequisites
=============

ComPWA is supposed to run on most modern unix systems (including MacOS). The following packages are mandatory:

* git (optional, for easier updates and if you want to contribute)
* cmake ( > 3.3 )
* gcc (> 5.1) or clang
* `Boost <http://www.boost.org/users/download/>`__\ , version >= 1.54

For a more feature rich installation, the following packages are highly recommended:

* python3 (for ComPWA expert system and python interface as well as a python plotting module)
* `GSL <https://www.gnu.org/software/gsl/>`__
* `ROOT <http://root.cern.ch/drupal/content/downloading-root>`__\ , version 5.34, 6.08
* `Minuit2 <http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/>`__\ , version 5.34, 6.08
* `Geneva <https://launchpad.net/geneva/+download>`__\ , version 1.8.0

In case that some dependencies are not met on your system use your package manager to install them or use the manual procedures described below. For MacOS you can use e.g. `MacPorts <https://www.macports.org>`_ as package manager.
You can also try to use different (newer) versions than the ones stated above, however those are not tested.

In order to install all dependencies it is probably also useful to have a look
on the `installation instructions file <https://github.com/ComPWA/ComPWA/blob/master/.travis.yml>`__ for TravisCI.

Quick Installation
==================
A detailed guide can be found below. The installation basically boils down to:

.. code-block:: shell

   git clone https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>
   cd <COMPWA_SOURCE_PATH> && git submodule init && git submodule update
   mkdir build && cd build
   cmake ../<COMPWA_SOURCE_PATH>
   make
   pip install ./pycompwa --user (optional)

Examples
========
The repository contains a couple of `examples <https://github.com/ComPWA/ComPWA/tree/master/Examples>`_. To learn about more detailed features of ComPWA you also might want to have a look on the examples of the `pycompwa package <https://github.com/ComPWA/ComPWA/tree/master/Examples/jupyter>`_.

Documentation
=============
Source code documentation via Doxygen is provided at ...
The master branch is automatically built using TravisCI. Probably it is interesting to check out the `log file <https://travis-ci.org/ComPWA/ComPWA>`_ and the projects TravisCI configuration file `travisCI.yml <https://github.com/ComPWA/ComPWA/blob/master/.travis.yml>`_.


ComPWA installation
===================
Manual installation of dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Boost**: to install Boost follow 
  `these <http://www.boost.org/doc/libs/1_54_0/more/getting_started/unix-variants.html#easy-build-and-install>`__ 
  instructions. The ``--prefix=path/to/installation/prefix`` option is useful
  as you can use the same path to point ComPWA to this specific Boost
  installation.

* **ROOT**: to install Root follow
  `these <http://root.cern.ch/drupal/content/installing-root-source>`_
  instructions.

* **Minuit2** is included in most ROOT installations. In case you want to
  install a stand-alone version follow
  `these <http://seal.web.cern.ch/seal/snapshot/work-packages/mathlibs/minuit/gettingStarted/autoconf.html>`__
  instructions. In addition you should use
  ``./configure --prefix=/your/minuit2/installation/path`` followed by
  ``make install`` to specify an installation directory which ComPWA later
  needs to find the library and to install all needed files to this location.

* **Geneva**: to install Geneva follow 
  `these <http://www.gemfony.eu/index.php?id=genevainstallation>`__ 
  instructions. The ``DCMAKE_INSTALL_PREFIX="/where/the/library/goes"`` option
  is useful as you can use the same path to point ComPWA to this specific 
  Geneva installation:

  .. code-block:: shell

        cd GENEVA_SOURCE
        mkdir -p build/install
        cd build
        cmake ../ -DCMAKE_INSTALL_PREFIX=./install
        make
        make install
        cp install/CMakeModules/FindGeneva.cmake YOUR_COMPWA_PATH/cmake/Modules/
        export GENEVA_ROOT=YOUR_GENEVA_PATH/build/install

  * Note for Fedora 25: The Geneva tests are build by default but might have trouble finding the boost test libraries of the Fedora boost package. A workaround is to disable them within ``YOUR_GENEVA_PATH/CMakeModules/CommonGenevaBuild.cmake, line 55`` (replace the line with ``SET( GENEVA_BUILD_TESTS FALSE )``.
  * Alternatively you can follow the instructions from the Geneva `manual <http://www.gemfony.eu/fileadmin/documentation/geneva-manual.pdf>`__.


Getting ComPWA
^^^^^^^^^^^^^^

To get the most recent version of the ComPWA framework clone its GitHub repository:

.. code-block:: shell

   git clone --recursive git@github.com:ComPWA/ComPWA <COMPWA_SOURCE_PATH>

this will clone the repository to the subfolder ``<COMPWA_SOURCE_PATH>`` within the current directory.
For multithreading ComPWA uses the parallel stl algorithms of c++17. Unfortunately the current compilers do not have any implementations for this. Here ComPWA currently relies on `TBB <https://github.com/01org/tbb>`_ and `parallelstl <https://github.com/intel/parallelstl>`_\ , which are included in ComPWA as git submodules. 


Building ComPWA
^^^^^^^^^^^^^^^

ComPWA uses ``cmake`` as build system. The usual steps to build all libraries and the test executable are the following:

* Create and enter a build folder (preferably not the ComPWA source folder)

  .. code-block:: shell

     mkdir build
     cd build

* Set your compiler if you do not use the system default compiler

  .. code-block:: shell

     export CC=<path_to_your_compiler> 
     export CXX=<path_to_your_compiler>

* Build the project. You can add ``-DCMAKE_INSTALL_PREFIX=<COMPWA_INSTALL_PATH>`` to specify an install location.

  .. code-block:: shell

     cmake ../<COMPWA_SOURCE_PATH> 
     make
     make install (optional)

Testing the ComPWA installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can run the test suite via:

.. code-block:: shell
   
   make test

or

.. code-block:: shell
   
   ctest

Other
^^^^^

* You can also use cmake to create a preconfigured project for an IDE (e.g.
  `eclipse <https://www.eclipse.org>`__ ):

  .. code-block:: shell
  
     cmake -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_CXX_COMPILER_ARG1=-std=c++14 ../<COMPWA_SOURCE_PATH>

Installation via Docker
^^^^^^^^^^^^^^^^^^^^^^^

A `Dockerfile <https://github.com/ComPWA/ComPWA/blob/master/Dockerfile>`__ for
ComPWA is provided. You can use it to build an 
`docker <https://www.docker.com>`__ image to run ComPWA. Using such an image
ComPWA should run on 
`all systems that are supported by docker <https://docs.docker.com/engine/installation/>`__
including several (commercial) cloud computing services. If you are new to
docker you can have a look on `this <https://prakhar.me/docker-curriculum/>`__
tutorial.

System specific notes
^^^^^^^^^^^^^^^^^^^^^

HimsterII / Mogon II
^^^^^^^^^^^^^^^^^^^^

`Mogon2 <https://hpc.uni-mainz.de/>`__ is the supercomputer of the Mainz
University. If you work on it you can fulfill the ComPWA 
`installation requirements <#requirements>`_ by loading a series of modules:

.. code-block:: shell

   module load devel/CMake/3.9.5
   module load toolchain/foss/2017a
   module load devel/Boost/1.65.1-foss-2017a
   module load numlib/GSL/2.4-foss-2017a
   module load ROOT/v6.12-foss-2017a-python3
   export CC=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/gcc
   export CXX=/cluster/easybuild/broadwell/software/compiler/GCCcore/6.3.0/bin/g++

Now follow: `Building ComPWA`_.
