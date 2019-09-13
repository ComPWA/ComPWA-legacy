// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "ChiOneD.hpp"

#include "Core/Function.hpp"

namespace ComPWA {

namespace Estimator {
namespace ChiOneD {

ChiOneD::ChiOneD(std::shared_ptr<ComPWA::Intensity> Intensity_,
                 const Data::DataSet &DataSample_)
    : Intensity(Intensity_), DataSample(DataSample_) {}

double ChiOneD::evaluate() noexcept {
  throw std::runtime_error("ChiOneD::evaluate: currently not implemented!");
}

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */
