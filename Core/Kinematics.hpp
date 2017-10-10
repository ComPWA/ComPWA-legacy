// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

//
// \file
// Kinematics inteface class
//

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <vector>
#include <string>
#include <complex>

#include "Core/Event.hpp"
#include "Core/Particle.hpp"
#include "Core/SubSystem.hpp"
#include "Core/Properties.hpp"
#include "Core/Spin.hpp"
#include "Core/DataPoint.hpp"

namespace ComPWA {

class Kinematics {
public:
  //! Constructor
  Kinematics(std::vector<pid> initial = std::vector<pid>(),
             std::vector<pid> finalS = std::vector<pid>())
      : _initialState(initial), _finalState(finalS),
        is_PS_area_calculated_(false), PS_area_(0.0){};

  /// Delete copy constructor. For each Kinematics in the analysis only
  /// one instance should exist since Kinematics does the bookkeeping for which
  /// SubSystems variables needs to be calculated. That instance can then be
  /// passed as (smart) pointer. Note: Not sure if we also should delete the
  /// move constructor.
  Kinematics(const Kinematics &that) = delete;

  /// Convert Event to dataPoint
  virtual void EventToDataPoint(const ComPWA::Event &ev,
                                dataPoint &point) const = 0;

  /// Check if dataPoint is within phase space boundaries
  virtual bool IsWithinPhsp(const dataPoint &point) const = 0;

  virtual double GetPhspVolume() const;

  /// Specify a phase space volume
  virtual void SetPhspVolume(double phsp);

  virtual std::size_t GetNVars() const { return _varNames.size(); }

  virtual std::vector<pid> GetFinalState() const { return _finalState; }

  virtual std::vector<pid> GetInitialState() const { return _initialState; }

  virtual ComPWA::FourMomentum GetInitialFourMomentum() const {
    return _initialP4;
  }

  virtual int GetDataID(const ComPWA::SubSystem &sys) = 0;

  virtual std::vector<std::string> GetVarNames() const { return _varNames; }
  
  virtual std::vector<std::string> GetVarTitles() const { return _varTitles; }

protected:
  std::vector<pid> _initialState;

  std::vector<pid> _finalState;

  /// Four momentum of the initial particle reaction
  ComPWA::FourMomentum _initialP4;

  /// Names of variabes
  std::vector<std::string> _varNames;

  /// (Latex) Titles of variables
  std::vector<std::string> _varTitles;

  virtual double calculatePSArea() const = 0;
  bool is_PS_area_calculated_;
  double PS_area_;
};

} // namespace ComPWA
#endif
