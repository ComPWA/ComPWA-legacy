/*! \mainpage ComPWA Documentation
 *
 * \section about About ComPWA
 *
 * ComPWA aims to provide a framework and toolkit for amplitude analysis in
 * a very general approach with clearly defined interfaces to extend its
 * features. Therefore, it can easily be extended for different experiments
 * and various models while maintaining a common structure, producing
 * comparable output, using the same fit-procedures and quality assurance
 * tools.
 *
 * \subsection structure Overview
 *
 * The basic modularization of the software can be seen in the following
 * picture. The four basic modules (Experiment, Physics, Estimation and
 * Optimization) each provide interfaces which need to be implemented
 * depending on your problem. They are represented in the framework by the
 * folder structure as well as in the build libraries. When attacking a new
 * problem using ComPWA, you might find that you can use existing module
 * implementations and only need to rewrite or extend part of it.
 *
 * \image html ComPWA_Modules.jpg
 * \image latex ComPWA_Modules.eps "ComPWA Modules" width=10cm
 *
 * \section install Installation
 *
 * \subsection step1 Build Tool: Boost.Build
 *
 * The ComPWA libraries and example-executables are build using the build tool
 * provided by the boost library. First, one needs to provide some paths to
 * external tools via environmental libraries, as can be seen in the setEnv
 * blank file. Secondly, you call the configure.pl which looks for available
 * modules and asks which you want to build. Now you can call bjam to compile
 * the code. All libraries will be located in the lib folder, all binaries in
 * the bin folder. Boost.Build uses Jamfiles for the build configuration. The
 * main file used when bjam is called is Jamroot in the root directory of
 * ComPWA. In the subdirectories of the interfaces and the actual
 * implementations you can find Jamfiles to manage the build process.
 *
 * \subsection step2 Used External Libraries
 *
 * The core Framework only needs a boost installation and a compiler
 * supporting c++11. But for the examples and most likely also for starting
 * some fits and plotting some output, it is recommended to have a ROOT
 * installation (for plots and maybe DataReader::Data storage/reading) and a
 * Minuit2
 * (for optimization) installation ready. Besides Minuit2 also Geneva can be
 * used as optimizer if a compatible installation is available.
 *
 * \section starting Getting Started
 * A good point to start learning how ComPWA works is by looking in the
 * Examples folder. Especially Examples/DalitzFit is a very simple fitting
 * routine.
 */

#error Documentation only.

/*!
 * @namespace ComPWA
 * General namespace for the ComPWA modules and classes.
 */

/*!
 * @namespace ComPWA::DataReader
 * Namespace for all ComPWA Data modules.
 */


/*!
 * @namespace ComPWA::Estimator
 * Namespace for all ComPWA Estimator modules.
 */


/*!
 * @namespace ComPWA::Optimizer
 * Namespace for all ComPWA Optimizer modules.
 */

/*!
 * @namespace ComPWA::Optimizer::Geneva
 * Wrapper to the Geneva library optimization tools.
 */

/*!
 * @namespace ComPWA::Optimizer::Minuit2
 * Wrapper to the Minuit2 library optimization tools.
 */

/*!
 * @namespace ComPWA::Physics
 * Namespace for all ComPWA Physics modules.
 */

/*!
 * @namespace ComPWA::Physics::HelicityFormalism
 * Module to calculate helicity amplitudes for general final states.
 */

/*!
 * @namespace ComPWA::Physics::QFT
 * qft++ library written by Mike Williams
 * (http://www-meg.phys.cmu.edu/williams/qft++/index.php/main_page)
 */

/*!
 * @namespace ComPWA::Physics::DecayDynamics
 * Dynamical decay functions (e.g. Breit-Wigner, Flatte ...)
 */
