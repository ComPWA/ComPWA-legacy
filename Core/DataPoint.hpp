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
private:
public:
  //! Default constructor
  dataPoint();

  //! Construct dataPoint from Event
  dataPoint(const Event &ev);
  
  //! Construct dataPoint from vector of invariant masses
  dataPoint(std::vector<double> vec);
  
  ~dataPoint(){};

  void Reset(unsigned int size);

  unsigned int Size() const { return var.size(); }

  //! Set value of coordinate num
  void SetValue(unsigned int pos, double val);
  
  //! Get value of coordinate num
  double GetValue(unsigned int num) const;

  //! Set coordinates by vector
  void SetPoint(std::vector<double> values);
  
  //! Get vector of coordinates
  std::vector<double>& GetPoint() { return var; };

  //! Set weight
  void SetWeight(double w) { weight = w; };
  //! Get weight
  double getWeight() const { return weight; };
  //! Set efficiency
  void SetEfficiency(double e) { eff = e; };
  //! Get efficiency
  double GetEfficiency() const { return eff; };

  static std::vector<double> GetRow(int n, std::vector<dataPoint> v) {
    std::vector<double> ret;
    if (!v.size())
      return ret;
    if (n >= Kinematics::instance()->GetNVars())
      throw std::runtime_error("dataPoint::getRow() | out of range!");
    for (int i = 0; i < v.size(); i++)
      ret.push_back(v.at(i).GetValue(n));
    return ret;
  }

protected:
  std::vector<double> var;
  double weight;
  double eff;
  friend std::ostream &operator<<(std::ostream &os, const dataPoint &p);
};

} /* namespace ComPWA */
#endif /*DPPOINT2_HPP_*/
