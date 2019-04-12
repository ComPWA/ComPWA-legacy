Installation
============

.. note::
   A detailed guide can be found below. The installation basically boils down to:

   * `git clone https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>`
   * `cd <COMPWA_SOURCE_PATH> && git submodule init && git submodule update`
   * `mkdir build && cd build`
   * `cmake ../<COMPWA_SOURCE_PATH>`
   * `make`
   * `pip install ./pycompwa --user` (optional)


Prerequisites
-------------

ComPWA is supposed to run on most modern unix systems (including MacOS). The following packages are mandatory:

* git (optional, for easier updates and if you want to contribute)
* cmake ( > 3.3 )
* gcc (> 5.1) or clang
* `Boost <http://www.boost.org/users/download/>`__\ , version >= 1.54

.. tip::
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


Manual installation of dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Boost**: to install Boost follow 
  `these <http://www.boost.org/doc/libs/1_54_0/more/getting_started/unix-variants.html#easy-build-and-install>`__ 
  instructions. The ``--prefix=path/to/installation/prefix`` option is useful
  as you can use the same path to point ComPWA to this specific Boost
  installation.

* **ROOT**: to install Root follow
  `these <http://root.cern.ch/drupal/content/installing-root-source>`__
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

  * Extract Geneva in a folder of your choice, lets call it "YOUR_GENEVA_PATH"
  * Go to ``YOUR_GENEVA_PATH/build``
  * Execute ``mkdir install``
  * Execute ``cmake ../ -DCMAKE_INSTALL_PREFIX="YOUR_GENEVA_PATH/build/install"``
  * Execute ``make``
  * Execute ``make install``
  * Execute ``cp install/CMakeModules/FindGeneva.cmake YOUR_COMPWA_PATH/cmake/Modules/``
  * Before compiling ComPWA, execute ``export GENEVA_ROOT=YOUR_GENEVA_PATH/build/install``
  * Note for Fedora 25: The Geneva tests are build by default but might have trouble finding the boost test libraries of the Fedora boost package. A workaround is to disable them within ``YOUR_GENEVA_PATH/CMakeModules/CommonGenevaBuild.cmake, line 55`` (replace the line with ``SET( GENEVA_BUILD_TESTS FALSE )``.
  * Alternatively you can follow the instructions from the Geneva `manual <http://www.gemfony.eu/fileadmin/documentation/geneva-manual.pdf>`__\ :

    * Go to ``/your/geneva/source/path/build``
    * ``cp ../scripts/genevaConfig.gcfg``
    * Modify ``./genevaConfig.gcfg`` to fit your needs.
    * ``../scripts/prepareBuild.sh ./genevaConfig.gcfg``

ComPWA installation
-------------------

Getting ComPWA
^^^^^^^^^^^^^^

To get the most recent version of the ComPWA framework clone its GitHub repository:

.. code-block:: shell

   git clone --recursive git@github.com:ComPWA/ComPWA <COMPWA_SOURCE_PATH>

this will clone the repository to the subfolder ``<COMPWA_SOURCE_PATH>`` within the current directory.
For multithreading ComPWA uses the parallel stl algorithms of c++17. Unfortunately the current compilers do not have any implementations for this. Here ComPWA currently relies on `TBB <https://github.com/01org/tbb>`_ and `parallelstl <https://github.com/intel/parallelstl>`_\ , which are included in ComPWA as git submodules. 


.. _build-compwa-label:

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

.. _setup-venv-label:

Installing python module
^^^^^^^^^^^^^^^^^^^^^^^^

During the build process an installable python module is created. The installation depends on your system. The most straightforward way would be:

.. code-block:: shell

   pip install ./pycompwa --user

The python module is also copied to the install location after ``make install``. You can find it at ``$CMAKE_INSTALL_PREFIX/share/ComPWA/pycompwa``. If you would like to user a virtual environment you could do something like:

.. code-block:: shell

   pipenv --python 3.xx
   pipenv install ./pycompwa

Here we have used `pipenv <https://github.com/pypa/pipenv>`_. Steps with the python3 default module `venv <https://docs.python.org/3/tutorial/venv.html>`_ are similar. If you would like to use `jupyter <https://jupyter.org/>`_ to perform your analysis you could create a custom jupyter kernel of your virtual environment:

.. code-block:: shell

   pipenv install ipykernel
   pipenv shell
   python -m ipykernel install --user --name=my-pycompwa-kernel
    

Testing the ComPWA installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can run the test suite via:

.. code-block:: shell
   
   make test

or

.. code-block:: shell
   
   ctest

The tests of the python module can be run via:

.. code-block:: shell
   
   cd pycompwa
   python setup.py pytest

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

Now follow :ref:`the build instructions <build-compwa-label>`.

Troubleshooting
---------------

Add content here
