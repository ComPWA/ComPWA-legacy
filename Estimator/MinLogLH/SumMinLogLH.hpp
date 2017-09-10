// Copyright (c) 2013, 2015, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

///
/// \file
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
///
class SumMinLogLH : public ComPWA::IEstimator {

public:
  SumMinLogLH();

protected:
  std::vector<std::shared_ptr<MinLogLH>> _minLogLh;
};

} // ns::Estimator
} // ns::ComPWA
#endif /* SumMinLogLH_hpp */
