// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//
// \file
// Kinematics inteface class
//

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <complex>
#include <string>
#include <vector>

#include "Core/DataPoint.hpp"
#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/Properties.hpp"
#include "Core/Spin.hpp"
#include "Core/SubSystem.hpp"

namespace ComPWA {

class Kinematics {
public:
  //! Constructor
  Kinematics(std::vector<pid> initial = std::vector<pid>(),
             std::vector<pid> finalS = std::vector<pid>(),
             ComPWA::FourMomentum cmsP4 = ComPWA::FourMomentum(0, 0, 0, 0))
      : InitialState(initial), FinalState(finalS), InitialStateP4(cmsP4),
        HasPhspVolume(false), PhspVolume(1.0){};

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping for which
  /// SubSystems variables needs to be calculated. That instance can then be
  /// passed as (smart) pointer. Note: Not sure if we also should delete the
  /// move constructor.
  Kinematics(const Kinematics &that) = delete;

  /// Convert Event to dataPoint
  virtual void convert(const ComPWA::Event &ev, DataPoint &point) const = 0;

  /// Check if dataPoint is within phase space boundaries
  virtual bool isWithinPhsp(const DataPoint &point) const = 0;

  virtual double phspVolume() const;

  virtual void setPhspVolume(double phsp);

  virtual std::size_t numVariables() const { return VariableNames.size(); }

  virtual std::vector<pid> finalState() const { return FinalState; }

  virtual std::vector<pid> initialState() const { return InitialState; }

  virtual ComPWA::FourMomentum initialStateFourMomentum() const {
    return InitialStateP4;
  }

  virtual unsigned int getDataID(const ComPWA::SubSystem &sys) const = 0;

  virtual std::vector<std::string> variableNames() const {
    return VariableNames;
  }

  virtual std::vector<std::string> variableTitles() const {
    return VariableTitles;
  }

protected:
  std::vector<pid> InitialState;

  std::vector<pid> FinalState;

  // we use a vector instead of a map here, due to cache optimizations
  std::vector<unsigned int> FinalStateEventPositionMapping;

  /// Four momentum of the initial particle reaction
  ComPWA::FourMomentum InitialStateP4;

  /// Names of variabes
  std::vector<std::string> VariableNames;

  /// (Latex) Titles of variables
  std::vector<std::string> VariableTitles;

  virtual double calculatePhspVolume() const = 0;
  bool HasPhspVolume;
  double PhspVolume;
};

} // namespace ComPWA
#endif
