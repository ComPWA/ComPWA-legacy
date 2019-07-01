// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_INTENSITY_HPP_
#define COMPWA_INTENSITY_HPP_

#include <memory>

#include "Core/FunctionTree.hpp"
#include "Core/ParameterList.hpp"

namespace ComPWA {

struct DataPoint;

///
/// \class Intensity
/// Pure interface class, resembling a real valued function. It can be evaluated
/// using DataPoints, which are elements of the domain of the function. For
/// example it can resemble a transition probability from an initial state to a
/// final state. However note that an Intensity can be, but does not have to be
/// normalized.
///
class OldIntensity : public Optimizable, public FunctionTreeInterface {
public:
  virtual ~OldIntensity() = default;

  /// evaluate intensity of model at \p point in phase-space
  virtual double evaluate(const DataPoint &point) const = 0;
};

} // namespace ComPWA
#endif
