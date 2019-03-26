// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

/*
 * Copyright (C) Gemfony scientific UG (haftungsbeschraenkt)
 *
 * See the AUTHORS file in the top-level directory for a list of authors.
 *
 * Contact: contact [at] gemfony (dot) com
 *
 * This file is part of the Geneva library collection.
 *
 * Geneva was developed with kind support from Karlsruhe Institute of
 * Technology (KIT) and Steinbuch Centre for Computing (SCC). Further
 * information about KIT and SCC can be found at http://www.kit.edu/english
 * and http://scc.kit.edu .
 *
 * Geneva is free software: you can redistribute and/or modify it under
 * the terms of version 3 of the GNU Affero General Public License
 * as published by the Free Software Foundation.
 *
 * Geneva is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with the Geneva library. If not, see <http://www.gnu.org/licenses/>.
 *
 * For further information on Gemfony scientific and Geneva, visit
 * http://www.gemfony.com .
 */

#include "Optimizer/Geneva/GenevaIF.hpp"
#include "Core/Logging.hpp"
#include "Optimizer/Geneva/GStartIndividual.hpp"
#include "Optimizer/Geneva/GenevaResult.hpp"

#include "geneva/Go2.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

GenevaIF::GenevaIF(std::shared_ptr<ComPWA::Estimator::Estimator> theData,
                   std::string inConfigFileDir)
    : Estimator(theData), ConfigFileDir(inConfigFileDir) {
  LOG(INFO) << "GenevaIF::GenevaIF() | "
               "Starting Geneva interface: config dir="
            << ConfigFileDir;
}

std::shared_ptr<ComPWA::FitResult> GenevaIF::exec(ParameterList &par) {
  std::shared_ptr<GenevaResult> result(new GenevaResult());
  ParameterList initialParList(par);

  using namespace Gem::Geneva;

  // IMPORTANT: for some reason the other constructor (with just a config file)
  // does not work correctly and the program is waiting endlessly! I did not
  // find any way to get it running except by handing the constructor below some
  // dummy arguments...
  int argc(1);
  char temp[] = {'a'};
  char *argv[] = {temp};
  Go2 go(argc, argv, ConfigFileDir + "Go2.json");

  // Initialize a client, if requested
  if (go.clientMode()) {
    std::cout << "Geneva Client waiting for action!" << std::endl;
    go.clientRun();
    return result;
  }

  std::shared_ptr<GStartIndividual> p(new GStartIndividual(Estimator, par));
  go.push_back(p);

  // TODO: here the minimization algorithms should be chosen by the user
  GEvolutionaryAlgorithmFactory ea(ConfigFileDir +
                                   "GEvolutionaryAlgorithm.json");

  GGradientDescentFactory f(ConfigFileDir + "GGradientDescentAlgorithm.json");

  go &ea() & f();

  // Perform the actual optimization
  std::shared_ptr<GStartIndividual> bestIndividual_ptr =
      go.optimize<GStartIndividual>();

  bestIndividual_ptr->getPar(par);
  result->setResult(bestIndividual_ptr);

  // write Parameters back
  ParameterList resultPar;
  bestIndividual_ptr->getPar(resultPar);
  for (unsigned int i = 0; i < par.doubleParameters().size(); i++) {
    if (!par.doubleParameter(i)->isFixed())
      par.doubleParameter(i)->setValue(resultPar.doubleParameter(i)->value());
  }

  return result;
}

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */
