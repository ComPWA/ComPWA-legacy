// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// Efficiency base class
///

#ifndef EFFICIENCY_HPP_
#define EFFICIENCY_HPP_

#include "Core/Logging.hpp"

namespace ComPWA {

struct DataPoint;

///
/// \class Efficiency
/// Base class for efficiency description over the phase space.
///
class Efficiency {
private:
public:
  Efficiency();

  virtual ~Efficiency();

  virtual double evaluate(const DataPoint &point) const = 0;
};

///
/// \class UnitEfficiency
/// Efficiency object with unit efficiency all over the phase space.
///
class UnitEfficiency : public Efficiency {
private:
public:
  UnitEfficiency() { LOG(DEBUG) << "Efficiency: creating UnitEfficiency!"; };

  ~UnitEfficiency(){};

  virtual double evaluate(const DataPoint &point) const { return 1; };
};

} // namespace ComPWA

#endif
