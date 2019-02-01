// Copyright (c) 2015, 2019 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//
// \file
// Kinematics interface class
//

#ifndef IKINEMATICS_HPP_
#define IKINEMATICS_HPP_

#include <string>
#include <vector>

namespace ComPWA {

class DataPoint;
class Event;

class IKinematics {

public:
  /// Convert Event to DataPoint
  virtual void convert(const ComPWA::Event &ev, DataPoint &point) const = 0;

  /// Check if DataPoint is within phase space boundaries
  virtual bool isWithinPhsp(const DataPoint &point) const = 0;

  virtual double phspVolume() const = 0;

};

} // namespace ComPWA
#endif

