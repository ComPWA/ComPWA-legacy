// Copyright (c) 2013, 2017 The ComPWA Team.
// This file is part of the ComPWA framework, check
// https://github.com/ComPWA/ComPWA/license.txt for details.

#include <cassert>
#include <memory>
#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "Core/ParameterList.hpp"
#include "Core/FitParameter.hpp"
#include "Core/Logging.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ComPWA::IEstimator> myData,
                     ComPWA::ParameterList &parList)
    : _myDataPtr(myData), _parList(parList) {
  if (0 == _myDataPtr)
    throw std::runtime_error("MinuitFcn::MinuitFcn() | Data pointer is 0!");
}

MinuitFcn::~MinuitFcn() {}

double MinuitFcn::operator()(const std::vector<double> &x) const {
  std::ostringstream paramOut;

  size_t pos = 0;
  for (auto p : _parList.doubleParameters()) {
    if (p->isFixed()) {
      pos++;
      continue;
    }
    paramOut << x.at(pos) << " "; // print only free parameters
    p->setValue(x.at(pos));
    pos++;
  }
  assert(x.size() == pos && "MinuitFcn::operator() | Number is (internal) "
                            "Minuit parameters and number of ComPWA "
                            "parameters does not match!");

  // Start timing
  clock_t begin = clock();
  double result = _myDataPtr->controlParameter(_parList);
  double sec = double(clock() - begin) / CLOCKS_PER_SEC;

  LOG(info) << "MinuitFcn: -log(L) = " << std::setprecision(10) << result
            << std::setprecision(4) << " Time: " << sec << "s"
            << " nCalls: " << _myDataPtr->status();
  LOG(debug) << "Parameters: " << paramOut.str();

  return result;
}

double MinuitFcn::Up() const {
  return 0.5; // TODO: Setter, LH 0.5, Chi2 1.
}

