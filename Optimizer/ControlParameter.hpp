//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
//
// This file is part of ComPWA, check license.txt for details
//
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
//! Optimizer Data Base-Class.
/*! \class ControlParameter
 * @file ControlParameter.hpp
 * This class provides the interface for the optimizer module to access the data.
 * If the access to the data is provided fulfilling this interface, then one can
 * use the same data and parameter-set with different optimizers.
*/

#ifndef CONTROLPARAMETER_HPP_
#define CONTROLPARAMETER_HPP_

#include <vector>
#include <memory>

#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Amplitude.hpp"

class ControlParameter{

public:
  static std::shared_ptr<ControlParameter> Instance();

  virtual double controlParameter(ParameterList& minPar) =0;

  virtual void resetInstance() { instance_ = std::shared_ptr<ControlParameter>(); };

  virtual std::shared_ptr<FunctionTree> getTree() { return std::shared_ptr<FunctionTree>(); }

  virtual std::shared_ptr<Amplitude> getAmplitude() { return std::shared_ptr<Amplitude>(); }

protected:
  static std::shared_ptr<ControlParameter> instance_;
std::shared_ptr<FunctionTree> f;
  ControlParameter(){
  }

  virtual ~ControlParameter(){ /* nothing */	}
 
};
#endif
