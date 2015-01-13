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
//! Estimator Interface Base-Class.
/*! \class Estimator
 * @file Estimator.hpp
 * This class provides the interface to classes which estimate the "closeness" of
 * the modeled intensity to the data. As it is pure virtual, one needs at least
 * one implementation to provide an estimator for the analysis. If a new estimator
 * is derived from and fulfills this base-class, no change in other modules are
 * necessary to work with the new estimation function. As it is derived from
 * ControlParameter, it can be used in the optimizer modules.
*/

#ifndef _ESTIMATOR_HPP_
#define _ESTIMATOR_HPP_

#include <vector>
#include <memory>

#include "Optimizer/ControlParameter.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"
#include "Core/Amplitude.hpp"

class Estimator : public ControlParameter
{

public:
  //static std::shared_ptr<Estimator> Instance();

  virtual double controlParameter(ParameterList& minPar) = 0;

  virtual bool hasTree() {return 0;}

  virtual std::shared_ptr<FunctionTree> getTree() { return std::shared_ptr<FunctionTree>(); }

  virtual std::shared_ptr<Amplitude> getAmplitude() { return std::shared_ptr<Amplitude>(); }

protected:
    //static std::shared_ptr<Estimator> instance_;

    Estimator(){
    }

    virtual ~Estimator(){
    /* nothing */
    }


};

//std::shared_ptr<Estimator> Estimator::Instance() {
 //   return Estimator::instance_;
//}

//std::shared_ptr<Estimator> Estimator::instance_ = 0;

//Why can i do this here but not inside ControlParameter.hpp?
/*std::shared_ptr<ControlParameter> ControlParameter::Instance() {
    return ControlParameter::instance_;
}

std::shared_ptr<ControlParameter> ControlParameter::instance_ = 0;*/

#endif
