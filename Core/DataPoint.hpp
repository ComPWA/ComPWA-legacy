// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// DataPoint class.
///

#ifndef DATAPOINT_HPP_
#define DATAPOINT_HPP_

#include <vector>
#include <ostream>

namespace ComPWA {

///
/// \class DataPoint
/// DataPoint stores the values that are needed to evaluate an AmpIntensity.
/// In case of the HelicityFormalism is would be a triple (m,theta,phi) for
/// each relevant SubSysyem.
/// The interface class IIntensity depends upon it. Dependencies should be 
/// therefore be as minimal as possible.
///
class DataPoint {

public:
  DataPoint();

  void reset(unsigned int size);

  std::size_t size() const { return Values.size(); }

  void setValue(unsigned int pos, double val);
  
  double value(unsigned int num) const;

  void setWeight(double w) { Weight = w; };
  
  double weight() const { return Weight; };
  
  void setEfficiency(double e) { Eff = e; };
  
  double efficiency() const { return Eff; };
  
  /// Reference to interval vector of values
  std::vector<double>& values() { return Values; };

  std::vector<double>::iterator first() { return Values.begin(); }

  std::vector<double>::iterator last() { return Values.end(); }

protected:
  std::vector<double> Values;
  double Weight;
  double Eff;
  friend std::ostream &operator<<(std::ostream &os, const DataPoint &p);
};

} // ns::ComPWA
#endif
