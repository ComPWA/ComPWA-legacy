// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_ESTIMATOR_ESTIMATOR_HPP_
#define COMPWA_ESTIMATOR_ESTIMATOR_HPP_

#include "Core/Function.hpp"

namespace ComPWA {
namespace Estimator {

///
/// This class template provides the interface to implementations, which
/// estimate the "closeness" of a Function to a data set, with respect to the
/// parameters of the Function.
///
/// The Estimator is defined as a Function with a return value, but without
/// input arguments.
///
/// Optimizer implementations use the Estimator to find the parameter set, that
/// model the data set optimally.
template <typename OutputType> class Estimator : public Function<OutputType> {};

} // namespace Estimator
} // namespace ComPWA

#endif
