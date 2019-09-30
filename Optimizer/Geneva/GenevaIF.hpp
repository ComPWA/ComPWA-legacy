// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//! Wrapper of the Geneva Optimizer library.
/*! \class GenevaIF
 * @file GenevaIF.hpp
 * This class provides a wrapper around the Geneva library. It fulfills the
 * Optimizer interface to be easily adapted to other modules. Parameters for the
 * optimization have to be provided in a config-file, the data needs to be
 * provided with the ControlParameter interface.
 */

#ifndef COMPWA_OPTIMIZER_GENEVAIF_HPP_
#define COMPWA_OPTIMIZER_GENEVAIF_HPP_

#include "Estimator/Estimator.hpp"
#include "Optimizer/Geneva/GenevaResult.hpp"
#include "Optimizer/Optimizer.hpp"

#include <memory>

namespace Gem {
namespace Geneva {
class GFMinIndividual;
}
} // namespace Gem

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

enum class AlgorithmTypes { EVOLUTIONARY, PARTICLE_SWARM, GRADIENT_DECENT };

class GenevaIF : public ComPWA::Optimizer::Optimizer<GenevaResult> {

public:
  GenevaIF(std::vector<ComPWA::Optimizer::Geneva::AlgorithmTypes>
               AlgorithmOrder_ = {AlgorithmTypes::EVOLUTIONARY,
                                  AlgorithmTypes::GRADIENT_DECENT},
           std::string ConfigFileDir_ = "./");
  GenevaResult optimize(Estimator::Estimator<double> &Estimator,
                        FitParameterList FitParameters) final;

  virtual ~GenevaIF() = default;

private:
  ComPWA::FitParameterList
  getFinalParameters(const FitParameterList &ParList,
                     std::shared_ptr<Gem::Geneva::GFMinIndividual> min) const;

  std::vector<ComPWA::Optimizer::Geneva::AlgorithmTypes> AlgorithmOrder;
  std::string ConfigFileDir;
};

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif
