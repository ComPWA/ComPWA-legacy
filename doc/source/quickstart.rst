Quickstart
==========

If you want to use the python interface to ComPWA and/or use the python
modules of ComPWA, setting up a virtual environment (venv) is highly 
recommended. Below are the setup instructions.

Installation & Setup
--------------------

.. note::
   
   Make sure the following software is installed on your system:
   
   * git (optional, for easier updates and if you want to contribute)
   * cmake ( > 3.3 )
   * gcc (> 5.1) or clang
   * `Boost <http://www.boost.org/users/download/>`_\ , version >= 1.54

   Not required, but recommended:

   * python3 + virtualenv (for ComPWA expert system and python interface as well as a python plotting module)
   * `ROOT <http://root.cern.ch/drupal/content/downloading-root>`_\ , version 5.34, 6.08

To install, simply run

.. code-block:: shell

   * `git clone https://github.com/ComPWA/ComPWA.git <COMPWA_SOURCE_PATH>`
   * `cd <COMPWA_SOURCE_PATH> && git submodule init && git submodule update`
   * `mkdir build && cd build`
   * `cmake ../<COMPWA_SOURCE_PATH>`
   * `make`

.. note::

   Here `<COMPWA_SOURCE_PATH>` points to the ComPWA source directory.

**Setup a python virtual environment**

.. code-block:: shell

   virtualenv -p python3 <PATH_OF_YOUR_VENV>
   source <PATH_OF_YOUR_VENV>/bin/activate
   pip install virtualenvwrapper
   source virtualenvwrapper.sh
   add2virtualenv <COMPWA_SOURCE_PATH>/Physics/ExpertSystem
   add2virtualenv <COMPWA_SOURCE_PATH>/Tools
   add2virtualenv <COMPWA_BUILD_DIR>/Tools/PythonInterface

.. note::
   Replace **$PATH_OF_YOUR_VENV** with the path where the venv should be installed.

   `<COMPWA_BUILD_DIR>` points to the ComPWA build directory.
   

**Install requirements for modules**
  
Each python module of ComPWA contains a requirements.txt file. If you want to
use this module simply install the requirements by executing:
  
.. code-block:: shell

   pip install -r <PATH_TO_COMPWA_PYTHON_MODULE>/requirements.txt
    
For example: ``pip install -r Physics/ExpertSystem/requirements.txt``
(assuming you are in the `<COMPWA_SOURCE_PATH>` directory)


Running/Usage
-------------

On how to use ComPWA please refer to the :ref:`Quickstart Example <example_quickstart>`.

This should get you started. You can check some of the other examples to learn
about more detailed features of ComPWA.

We would be happy to recieve some feedback or contributions ;)!