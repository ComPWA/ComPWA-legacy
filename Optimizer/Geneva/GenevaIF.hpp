//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Wrapper of the Geneva Optimizer library.
/*! \class GenevaIF
 * @file GenevaIF.hpp
 * This class provides a wrapper around the Geneva library. It fulfills the
 * Optimizer interface to be easily adapted to other modules. Parameters for the
 * optimization have to be provided in a config-file, the data needs to be
 * provided with the ControlParameter interface.
*/

#ifndef _GENEVAIF_HPP
#define _GENEVAIF_HPP


// Boost header files go here

// Geneva header files go here

#include <vector>
#include <memory>
#include <iostream>
//#include <boost/shared_ptr.hpp>
#include "Core/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"
#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Optimizer/Geneva/GStartIndividual.hpp"
//#include "Optimizer/Geneva/GArgumentParser.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// Geneva header files go here
/*#include <courtier/GAsioHelperFunctions.hpp>
#include <courtier/GAsioTCPClientT.hpp>
#include <courtier/GAsioTCPConsumerT.hpp>
#include <geneva/GBrokerEA.hpp>
#include <geneva/GEvolutionaryAlgorithm.hpp>
#include <geneva/GIndividual.hpp>
#include <geneva/GMultiThreadedEA.hpp>
#include <common/GCommonEnums.hpp>
#include <common/GSerializationHelperFunctionsT.hpp>
#include <geneva/GOptimizationEnums.hpp>*/
#include "geneva/Go2.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

class GenevaIF : public Optimizer {

public:
  /// Default Constructor (0x0)
  GenevaIF(std::shared_ptr<IEstimator> theData, std::string inConfigFileDir="test/config/");
  virtual std::shared_ptr<FitResult> exec(ParameterList& par) ;

  /** Destructor */
  virtual ~GenevaIF();

  virtual void setServerMode();
  virtual void setClientMode(std::string serverip="localhost", unsigned int serverport=10000);

 protected:

private:
  std::shared_ptr<IEstimator> _myData;
  std::string configFileDir;
  //Gem::Geneva::parMode parallelizationMode;
  Gem::Geneva::execMode parallelizationMode;
  Gem::Common::serializationMode serMode;
  bool clientMode;
  std::string ip;
  unsigned int port;
 // vector<string> paramNames;

  /*
  boost::uint16_t parallelizationMode;
  bool serverMode;
  std::string ip;
  unsigned short port;
  Gem::Common::serializationMode serMode;
  boost::uint16_t nProducerThreads;
  boost::uint16_t nEvaluationThreads;
  std::size_t populationSize;
  std::size_t nParents;
  boost::uint32_t maxIterations;
  long maxMinutes;
  boost::uint32_t reportIteration;
  Gem::Geneva::recoScheme rScheme;
  std::size_t arraySize;
  Gem::Geneva::sortingMode smode;
  boost::uint32_t processingCycles;
  bool returnRegardless;
  boost::uint32_t waitFactor;*/

};

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif /* _GENEVAIF_HPP */
