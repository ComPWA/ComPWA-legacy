C++ Modules
===========

All of the c++ code belongs to one of the 4 main module groups in ComPWA,
which are:

- Data
- Physics/Models
- Esimators
- Optimizers

There is another supplementary module group `Tools`, which contains for example
integration algorithms.

Data Modules
------------

The Data modules is responsible for:

- transformation of different data structures
- data input/output
- data generation (generators should be moved from tools to here)


Physics Modules
---------------

The Physics modules

- define Intensities and Amplitudes of physics models/theories
  (i.e. helicty formalism, ...)
- define Kinematics classes, which are responsible for transforming Event based
  data to the kinematic variables that are needed by the underlying theory.

Estimator Modules
-----------------

An estimator is responsible for determining the "closeness" of the model to the data.
Each estimator module has a specific way to estimate the closeness.

Optimizer Modules
-----------------

The Optimizer modules are responsible for finding the "optimal" set of
parameters of a model, by minimizing the "closeness" determined by a specific
Estimator.
Each Optimizer module implements a specific algorithm for finding the optimum.