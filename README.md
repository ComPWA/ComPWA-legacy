[![Build Status](https://travis-ci.org/ComPWA/ComPWA.svg?branch=master)](https://travis-ci.org/ComPWA/ComPWA)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/13697/badge.svg)](https://scan.coverity.com/projects/compwa-compwa)

![ComPWA Logo](https://github.com/ComPWA/ComPWA/blob/master/doc/images/logo.png)

## About
ComPWA is a project to provide a flexible and modular Partial Wave Analysis framework for various use-cases. It was initiated and is developed by the Panda Collaboration (antiProton ANnihilation at DArmstadt) at Fair (Facility for Antiproton and Ion Research). But ComPWA will not just be used for the Panda physics program, but also for various other experiments to provide a commonly used tool which is stable, efficient and provides comparable results. At the moment there are many PWA-tools on the market, most used just for specific experiments and specific physics cases, some experiments even have multiple tools. But why write the same software again and again? E.g. the model describing physical processes should stay the same independent where and how there was a measurement of the process. Using the actual same implementation of the model does not only save a lot of time, it also ensures that two experiments are comparing the same thing. The same argument holds for optimization-routines and estimation-functions. It might even allow combined fitting of different experiments instead of taking the average of the results!

The natural modularization, following the considerations above, would be to separate into experiment specific information, physics (models, formalisms), estimation how good the model fits the data and optimization of free parameters. The first considerations on this where discussed with experts from different experiments and different technologies where discussed and tested. The result of this process is the first requirement document of the new PWA-Framework.
This sketch illustrates the modular concept: 
![ComPWA Modules](https://github.com/ComPWA/ComPWA/wiki/fw.png)

## Available Features :sparkles:
#### Physic Models:
- Helicity Formalism
- Canonical Formalism
#### Data Formats: 
- ROOT
#### Estimators:
- Chi Square
- Unbinned-LogLikelihood
#### Optimizers:
- Miniuit2
- Geneva
#### User Interfaces/Steering:
- C++
- Python

## Example
coming soon

## Installation
A detailed guide can be found [here](https://github.com/ComPWA/ComPWA/wiki/Installation). The installation basically boils down to:
```bash
git clone https://github.com/ComPWA/ComPWA.git
cd ComPWA
git submodule init
git submodule update
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=<COMPWA_INSTALL_PATH> .. # -DCMAKE_INSTALL_PREFIX is optional
make -j2
make install # this step is optional
```

## Documentation
The [Doxygen documentation](http://ComPWA.github.io/ComPWA) of the current master branch is located on the github pages of ComPWA.
The master branch is automatically build using TravisCI. Probably it is interesting to check out the [log file](https://travis-ci.org/ComPWA/ComPWA) and the projects TravisCI configuration file [<code>.travisCI.yml</code>](https://github.com/ComPWA/ComPWA/blob/master/.travis.yml).
