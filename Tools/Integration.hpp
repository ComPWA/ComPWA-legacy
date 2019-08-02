// Copyright (c) 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_TOOLS_INTEGRATION_HPP_
#define COMPWA_TOOLS_INTEGRATION_HPP_

#include "Core/Function.hpp"

namespace ComPWA {

namespace Data {
struct DataSet;
}

namespace Tools {

double integrate(ComPWA::Intensity &intensity,
                 const ComPWA::Data::DataSet &phspsample,
                 double phspVolume = 1.0);

double maximum(ComPWA::Intensity &intensity,
               const ComPWA::Data::DataSet &sample);

} // namespace Tools
} // namespace ComPWA

#endif
