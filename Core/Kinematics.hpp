// Copyright (c) 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <vector>
#include <string>
#include <complex>

#include "Core/Event.hpp"
#include "Core/SubSystem.hpp"
#include "Core/PhysConst.hpp"
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
  Kinematics(const Kinematics& that) = delete;
  
  //! converts Event to dataPoint
  virtual void EventToDataPoint(const ComPWA::Event &ev,
                                dataPoint &point) const = 0;

  //! Checks of data point is within phase space boundaries
  virtual bool IsWithinPhsp(const dataPoint &point) const = 0;

  //! calculated the PHSP volume of the current decay by MC integration
  virtual double GetPhspVolume();

  //! calculated the PHSP volume of the current decay by MC integration
  virtual void SetPhspVolume( double phsp );
  
  //! Get number of variables
  virtual std::size_t GetNVars() const { return _varNames.size(); }

  //! Get final state
  virtual std::vector<pid> GetFinalState() { return _finalState; }

  //! Get inital state
  virtual std::vector<pid> GetInitialState() { return _initialState; }

  virtual int GetDataID(const ComPWA::SubSystem s) = 0;
  
protected:
  std::vector<pid> _initialState;
  std::vector<pid> _finalState;

  //! Internal names of variabes
  std::vector<std::string> _varNames;
  //! Latex titles for variables
  std::vector<std::string> _varTitles;

  virtual double calculatePSArea() = 0;
  bool is_PS_area_calculated_;
  double PS_area_;
};

} /* namespace ComPWA */
#endif /* KINEMATICS_HPP_ */
