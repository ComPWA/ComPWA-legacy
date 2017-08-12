// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef EFFICIENCY_HPP_
#define EFFICIENCY_HPP_

#include <iostream>
#include <vector>
#include "Core/Logging.hpp"

namespace ComPWA {

class dataPoint;

/**
 *  \class DalitzEfficiency
 *  \brief Virtual efficiency class
 */
class Efficiency {
private:
public:
  Efficiency();

  virtual ~Efficiency();

  virtual double Evaluate(const dataPoint &point) const = 0;
};

/**
 *  \class UnitEfficiency
 *  \brief implementation of virtual class efficiency. Efficiency ist constant
 * one allover the PHSP
 */
class UnitEfficiency : public Efficiency {
private:
public:
  UnitEfficiency() { LOG(debug) << "Efficiency: creating UnitEfficiency!"; };
  ~UnitEfficiency(){};
  virtual double Evaluate(const dataPoint &point) const { return 1; };
};

} /* namespace ComPWA */

#endif /* EFFICIENCY_HPP_ */
