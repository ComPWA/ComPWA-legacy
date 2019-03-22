// Copyright (c) 2013, 2018 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// IIntensity interface class.
///

#ifndef IINTENSITY_HPP_
#define IINTENSITY_HPP_

#include <memory>
#include <vector>

#include "Core/DataPoint.hpp"

namespace ComPWA {

///
/// \class IIntensity
/// Interface class for Intensity.
///
class IIntensity {
public:
  /// Evaluate intensity of model at \p point in phase-space
  virtual double intensity(const DataPoint &point) const = 0;

};

} // namespace ComPWA
#endif

