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

#include <memory>

#include "Core/Parameter.hpp"
#include "Core/ParameterList.hpp"
#include "Estimator/Estimator.hpp"
#include "Optimizer/Optimizer.hpp"

namespace ComPWA {
namespace Optimizer {
namespace Geneva {

class GenevaIF : public Optimizer {

public:
  GenevaIF(std::shared_ptr<ComPWA::Estimator::Estimator> theData,
           std::string inConfigFileDir = "./");
  std::shared_ptr<ComPWA::FitResult> exec(ParameterList &par) final;

  virtual ~GenevaIF() = default;

private:
  std::shared_ptr<ComPWA::Estimator::Estimator> Estimator;
  std::string ConfigFileDir;
};

} /* namespace Geneva */
} /* namespace Optimizer */
} /* namespace ComPWA */

#endif
