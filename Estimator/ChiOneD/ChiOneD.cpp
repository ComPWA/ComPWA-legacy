// Copyright (c) 2013 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include "ChiOneD.hpp"

#include "Core/Intensity.hpp"

namespace ComPWA {

namespace Estimator {
namespace ChiOneD {

ChiOneD::ChiOneD(std::shared_ptr<ComPWA::Intensity> intensity,
                 const std::vector<DataPoint> &points)
    : Intensity(intensity), DataPoints(points) {}

double ChiOneD::evaluate() const {
  throw std::runtime_error("ChiOneD::evaluate: currently not implemented!");
}

} /* namespace ChiOneD */
} /* namespace Estimator */
} /* namespace ComPWA */
