// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
/// SumMinLogLH class.
///

#ifndef _SUMMINLOGLH_HPP
#define _SUMMINLOGLH_HPP

#include <vector>
#include <memory>
#include <string>

#include "Core/Estimator.hpp"
#include "Core/AmpIntensity.hpp"
#include "DataReader/Data.hpp"
#include "Core/Event.hpp"
#include "Core/ParameterList.hpp"
#include "Core/FunctionTree.hpp"

#include "Estimator/MinLogLH/MinLogLH.hpp"

namespace ComPWA {
namespace Estimator {

///
/// \class SumMinLogLH
/// Calculates the combined likelihood of multiple MinLogLH.
///
class SumMinLogLH : public ComPWA::IEstimator {

public:
  SumMinLogLH();

  /// Value of minimum log likelhood function.
  virtual double ControlParameter(ComPWA::ParameterList &par);

  virtual void AddLogLh(std::shared_ptr<ComPWA::Estimator::MinLogLH> logLh) {
    _minLogLh.push_back(logLh);
  }
  
  virtual void UseFunctionTree(bool onoff);

protected:
  std::vector<std::shared_ptr<ComPWA::Estimator::MinLogLH>> _minLogLh;

  std::shared_ptr<ComPWA::FunctionTree> _tree;
};

} // ns::Estimator
} // ns::ComPWA
#endif
