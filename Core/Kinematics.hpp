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

#ifndef KINEMATICS_HPP_
#define KINEMATICS_HPP_

#include <vector>
#include <string>
#include <complex>

#include "Core/Event.hpp"
#include "Core/PhysConst.hpp"
#include "Core/Spin.hpp"

namespace ComPWA {

class dataPoint;

class Kinematics {
public:
  //! Constructor
  Kinematics(std::vector<pid> initial = std::vector<pid>(),
             std::vector<pid> finalS = std::vector<pid>())
      : _initialState(initial), _finalState(finalS),
        is_PS_area_calculated_(false), PS_area_(0.0){};

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
