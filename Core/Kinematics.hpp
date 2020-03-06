// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef COMPWA_KINEMATICS_HPP_
#define COMPWA_KINEMATICS_HPP_

#include <string>
#include <vector>

namespace ComPWA {

struct Event;

namespace Data {
struct DataSet;
}

/// The Kinematics interface defines the conversion of Events to a DataSet.
class Kinematics {
public:
  virtual ~Kinematics() = default;

  virtual ComPWA::Data::DataSet
  convert(const std::vector<ComPWA::Event> &Events) const = 0;

  /// checks if DataPoint is within phase space boundaries
  virtual std::vector<ComPWA::Event>
  reduceToPhaseSpace(const std::vector<ComPWA::Event> &Events) const = 0;

  virtual double phspVolume() const = 0;

  /// Get a vector of PIDs of the final state.
  /// This interface allows the user to use the info in this object to interpret
  /// momentum tuples in a date file.
  virtual const std::vector<int> &getFinalStatePIDs() const = 0;
};

} // namespace ComPWA
#endif
