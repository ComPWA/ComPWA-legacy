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
#include "Optimizer/Geneva/GFMinIndividual.hpp"
#include "Optimizer/Geneva/GenevaResult.hpp"

#include "geneva/Go2.hpp"

#include <chrono>

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

GenevaIF::GenevaIF(
    std::vector<ComPWA::Optimizer::Geneva::AlgorithmTypes> AlgorithmOrder_,
    std::string ConfigFileDir_)
    : AlgorithmOrder(AlgorithmOrder_), ConfigFileDir(ConfigFileDir_) {
  LOG(INFO) << "GenevaIF::GenevaIF(): Starting Geneva interface (config dir="
            << ConfigFileDir << ")";
}

GenevaResult GenevaIF::optimize(Estimator::Estimator<double> &Estimator,
                                FitParameterList FitParameters) {
  if (!isValid(FitParameters, Estimator.getParameters()))
    return GenevaResult();

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
    LOG(INFO) << "Geneva Client waiting for action!";
    go.clientRun();
    return GenevaResult();
  }

  std::vector<double> initialpars;
  for (auto const &x : FitParameters)
    initialpars.push_back(x.Value);
  Estimator.updateParametersFrom(initialpars);
  double InitialEstimatorValue(Estimator.evaluate());

  std::shared_ptr<GFMinIndividual> p(
      new GFMinIndividual(Estimator, FitParameters));
  go.push_back(p);

  for (auto AlgorithmType : AlgorithmOrder) {
    switch (AlgorithmType) {
    case AlgorithmTypes::EVOLUTIONARY: {
      GEvolutionaryAlgorithmFactory eaf(ConfigFileDir + "GEvolutionary.json");
      go &eaf();
      break;
    }
    case AlgorithmTypes::PARTICLE_SWARM: {
      GSwarmAlgorithmFactory sf(ConfigFileDir + "GGradientDescent.json");
      go &sf();
      break;
    }
    case AlgorithmTypes::GRADIENT_DECENT: {
      GGradientDescentFactory gdf(ConfigFileDir + "GGradientDescent.json");
      go &gdf();
      break;
    }
    }
  }

  std::chrono::steady_clock::time_point StartTime =
      std::chrono::steady_clock::now();
  // Perform the actual optimization
  std::shared_ptr<GFMinIndividual> bestIndividual_ptr =
      go.optimize<GFMinIndividual>();
  std::chrono::steady_clock::time_point EndTime =
      std::chrono::steady_clock::now();
  auto finalFitPars = getFinalParameters(FitParameters, bestIndividual_ptr);
  return GenevaResult(ComPWA::FitResult{
      FitParameters, finalFitPars, finalFitPars.size(), InitialEstimatorValue,
      std::get<1>(bestIndividual_ptr->getBestKnownPrimaryFitness()),
      std::chrono::duration_cast<std::chrono::seconds>(EndTime - StartTime)});
}

ComPWA::FitParameterList GenevaIF::getFinalParameters(
    const FitParameterList &ParList,
    std::shared_ptr<Gem::Geneva::GFMinIndividual> BestIndividual) const {

  std::vector<double> finalpars;
  BestIndividual->streamline(finalpars);

  FitParameterList FinalParameters;
  size_t counter(0);
  for (auto const &p : ParList) {
    if (p.IsFixed)
      continue;
    // TODO: error estimation
    FinalParameters.push_back(
        ComPWA::FitParameter<double>(p.Name, finalpars[counter], false));
    ++counter;
  }
  return FinalParameters;
}

} // namespace Geneva
} // namespace Optimizer
} // namespace ComPWA
