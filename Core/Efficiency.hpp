//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//		Peter Weidenkaff -
//-------------------------------------------------------------------------------
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

  virtual double evaluate(std::vector<double> x) = 0;

  virtual double evaluate(const dataPoint &point) = 0;
};

/**
 *  \class UnitEfficiency
 *  \brief implementation of virtual class efficiency. Efficiency ist constant
 * one allover the PHSP
 */
class UnitEfficiency : public Efficiency {
private:
public:
  UnitEfficiency() { LOG(info) << "Efficiency: creating UnitEfficiency!"; };
  ~UnitEfficiency(){};
  virtual double evaluate(std::vector<double> x) { return 1; };
  virtual double evaluate(const dataPoint &point) { return 1; };
};

} /* namespace ComPWA */

#endif /* EFFICIENCY_HPP_ */
