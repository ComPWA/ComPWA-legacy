//-------------------------------------------------------------------------------
// Copyright (c) 2013 Mathias Michel.
// All rights reserved. This program and the accompanying materials
// are made available under the terms of the GNU Public License v3.0
// which accompanies this distribution, and is available at
// http://www.gnu.org/licenses/gpl.html
//
// Contributors:
//     Mathias Michel - initial API and implementation
//-------------------------------------------------------------------------------
#include <cassert>
#include <memory>
#include <iostream>
#include <cmath>
#include <iomanip>

#include <boost/chrono.hpp>
namespace bc = boost::chrono;

#include "Core/ParameterList.hpp"
#include "Core/Parameter.hpp"
#include "Core/Logging.hpp"
#include "Optimizer/Minuit2/MinuitFcn.hpp"
#include "Core/ControlParameter.hpp"

using namespace ROOT::Minuit2;

MinuitFcn::MinuitFcn(std::shared_ptr<ComPWA::ControlParameter> myData,
                     ComPWA::ParameterList &parList)
    : _myDataPtr(myData), _parList(parList) {
  if (0 == _myDataPtr)
    throw std::runtime_error("MinuitFcn::MinuitFcn() | Data pointer is 0!");
}

MinuitFcn::~MinuitFcn() {}

double MinuitFcn::operator()(const std::vector<double> &x) const {
  // ParameterList par;
  std::ostringstream paramOut;
  for (unsigned int i = 0; i < x.size(); i++) {
    std::shared_ptr<ComPWA::DoubleParameter> actPat =
        _parList.GetDoubleParameter(i);
    // std::cout<<i<<" "<<actPat->GetName()<<" "<<actPat->GetValue()
    //<<" "<<x[i]<<" "<<actPat->IsFixed()<<std::endl;
    if (!actPat->IsFixed())
      if (x[i] == x[i]) {
        actPat->SetValue(x[i]);
        paramOut << x[i] << " "; // print only free parameters
      }
  }
  bc::system_clock::time_point start = bc::system_clock::now();
  double result = _myDataPtr->controlParameter(_parList);
  bc::duration<double> sec = bc::system_clock::now() - start;

  LOG(info) << "MinuitFcn: -log(L) = " << std::setprecision(10) << result
            << std::setprecision(4) << " Time: " << sec.count() << "s"
            << " nCalls: " << _myDataPtr->nCalls();
  LOG(debug) << "Parameters: " << paramOut.str();

  return result;
}

double MinuitFcn::Up() const {
  return 0.5; // TODO: Setter, LH 0.5, Chi2 1.
}

