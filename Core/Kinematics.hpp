// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_KINEMATICS_HPP_
#define COMPWA_KINEMATICS_HPP_

namespace ComPWA {

struct DataPoint;
struct Event;

/// The Kinematics interface is responsible for converting an Event into a
/// DataPoint.
class Kinematics {
public:
  virtual ~Kinematics() = default;

  virtual DataPoint convert(const ComPWA::Event &event) const = 0;

  virtual std::vector<std::string> getKinematicVariableNames() const = 0;

  /// checks if DataPoint is within phase space boundaries
  virtual bool isWithinPhaseSpace(const DataPoint &point) const = 0;

  virtual double phspVolume() const = 0;
};

} // namespace ComPWA
#endif
