//-------------------------------------------------------------------------------
// Copyright (c) 2013 Peter Weidenkaff.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Peter Weidenkaff - initial API
//-------------------------------------------------------------------------------

//! dataPoint stores kinematic information of a dalitz plot
/*!
 * @file DataPoint.hpp
 *\class dataPoint
 *      dataPoint is a singleton class which provides a
 *      certain phase space point to all classes of the framework. The point can
 *be set anywhere in
 *      the framework and can be read by any amplitude class.
 */

#ifndef DPPOINT2_HPP_
#define DPPOINT2_HPP_

#include <cstdlib>
#include <math.h>
#include <vector>
#include <map>
#include "Core/Kinematics.hpp"
#include "Core/Efficiency.hpp"
#include "Core/Event.hpp"

namespace ComPWA {

class dataPoint {

public:
  dataPoint();

  /// Construct dataPoint from Event
  dataPoint(const Event &ev);
  
  /// Construct dataPoint from vector of invariant masses
  dataPoint(std::vector<double> vec);
  
  ~dataPoint(){};

  void Reset(unsigned int size);

  std::size_t Size() const { return var.size(); }

  void SetValue(unsigned int pos, double val);
  
  double GetValue(unsigned int num) const;

  std::vector<double>& GetPoint() { return var; };

  void SetWeight(double w) { weight = w; };
  
  double getWeight() const { return weight; };
  
  void SetEfficiency(double e) { eff = e; };
  
  double GetEfficiency() const { return eff; };

protected:
  std::vector<double> var;
  double weight;
  double eff;
  friend std::ostream &operator<<(std::ostream &os, const dataPoint &p);
};

} /* namespace ComPWA */
#endif /*DPPOINT2_HPP_*/
