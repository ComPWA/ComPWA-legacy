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

/// Calculate integral and its error
///
/// We follow https://en.wikipedia.org/wiki/Monte_Carlo_integration.
/// The average intensity is given by:
/// \f[
///   \langle f \rangle =\frac{1}{N} \sum_{i=1}^N f(\overline{\mathbf{x}}_i)
/// \f]
/// and the integral estimate Q_N using a sample size of N events is given by:
/// \f[
///   I \approx Q_N = V*\langle f \rangle.
/// \f]
/// The variance is given by
/// \f[
/// \mathrm{Var}(Q_N) = \frac{V^2}{N} \frac{1}{(N-1)} \sum_{i=1}^N \left
/// (f(\overline{\mathbf{x}}_i) - \langle f \rangle \right )^2
/// \f]
/// We use a Kahan summation to improve numerical stability. We return
/// \f$ \sqrt{\mathrm{Var}(Q_N)} \f$ as error.
std::pair<double, double>
integrateWithError(ComPWA::Intensity &intensity,
                   const ComPWA::Data::DataSet &phspsample,
                   double phspVolume = 1.0);

double integrate(ComPWA::Intensity &intensity,
                 const ComPWA::Data::DataSet &phspsample,
                 double phspVolume = 1.0);

double maximum(ComPWA::Intensity &intensity,
               const ComPWA::Data::DataSet &sample);

} // namespace Tools
} // namespace ComPWA

#endif
